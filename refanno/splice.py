import gzip
import os
import subprocess

from .audit import file_sha256, read_json, write_json
from .adjudication import chain_groups


def pick_consensus_models(loci, min_sets=2):
    train = []
    for locus in loci:
        for group in chain_groups(locus):
            if len({model.set_name for model in group}) >= min_sets:
                train.append(max(group, key=lambda m: (m.cds_len, m.key)))
    return train


def bed12_line(model):
    blocks = sorted(model.cds)
    start0 = blocks[0][0] - 1
    sizes = ",".join(str(e - s + 1) for s, e in blocks)
    starts = ",".join(str(s - 1 - start0) for s, e in blocks)
    thick0 = blocks[0][0] - 1
    thick1 = blocks[-1][1]
    return "\t".join([
        model.chrom, str(start0), str(blocks[-1][1]), model.tx_id, "0", model.strand,
        str(thick0), str(thick1), "0", str(len(blocks)), sizes, starts,
    ])


def pick_training_models(models, loci, trusted, min_training=500):
    consensus = pick_consensus_models(loci)
    if len(consensus) >= min_training:
        return consensus, "consensus"
    best = {}
    for model in models:
        key = (model.set_name, model.gene_id)
        current = best.get(key)
        if current is None or model.cds_len > current.cds_len:
            best[key] = model
    trusted_models = [model for (set_name, _), model in best.items() if set_name.lower() in trusted]
    if len(trusted_models) >= min_training:
        return trusted_models, "trusted"
    return list(best.values()), "all"


def _run(command, stdout_path, missing_hint):
    temporary = stdout_path + ".tmp"
    try:
        with open(temporary, "w") as handle:
            subprocess.run(command, stdout=handle, check=True)
        os.replace(temporary, stdout_path)
    except FileNotFoundError as exc:
        raise RuntimeError(missing_hint) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"minisplice 运行失败（{' '.join(command[:2])}），退出码 {exc.returncode}") from exc
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def minisplice_version(executable):
    try:
        result = subprocess.run([executable, "version"], text=True, capture_output=True, timeout=15)
        text = (result.stdout or result.stderr).strip()
        return text.splitlines()[0] if text else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def run_minisplice(executable, genome, threads, outdir, train_models, annot_models, model_path=None):
    os.makedirs(outdir, exist_ok=True)
    output = os.path.join(outdir, "splice_scores.tsv")
    cache_path = os.path.join(outdir, "minisplice.cache.json")
    train_bed = os.path.join(outdir, "minisplice.train.bed12")
    annot_bed = os.path.join(outdir, "minisplice.annot.bed12")
    with open(train_bed, "w") as handle:
        for model in sorted(train_models, key=lambda m: m.key):
            handle.write(bed12_line(model) + "\n")
    with open(annot_bed, "w") as handle:
        for model in sorted(annot_models, key=lambda m: m.key):
            handle.write(bed12_line(model) + "\n")
    fingerprint = {
        "genome_sha256": file_sha256(genome),
        "minisplice_version": minisplice_version(executable),
        "threads": threads,
        "model_sha256": file_sha256(model_path) if model_path else "",
        "train_bed_sha256": "" if model_path else file_sha256(train_bed),
        "annot_bed_sha256": file_sha256(annot_bed),
    }
    if os.path.isfile(output) and os.path.getsize(output) > 0 and os.path.isfile(cache_path):
        try:
            if read_json(cache_path) == fingerprint:
                print(f"minisplice: reuse {output}", flush=True)
                return output
        except (OSError, ValueError):
            pass
    hint = "找不到 minisplice；它不在 conda 仓库中，请从 https://github.com/lh3/minisplice 编译后用 --minisplice 指定可执行文件路径"
    kan = model_path or os.path.join(outdir, "minisplice.kan")
    cali = kan + ".cali"
    if not model_path:
        train_txt = os.path.join(outdir, "minisplice.train.txt")
        _run([executable, "gentrain", train_bed, genome], train_txt, hint)
        command = [executable, "train", "-t", str(threads), "-o", kan + ".tmp", train_txt]
        print("minisplice: " + " ".join(command), flush=True)
        try:
            subprocess.run(command, check=True)
            os.replace(kan + ".tmp", kan)
        except FileNotFoundError as exc:
            raise RuntimeError(hint) from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"minisplice train 运行失败，退出码 {exc.returncode}") from exc
    if not os.path.isfile(cali):
        command = [executable, "predict", "-t", str(threads), "-b", annot_bed, kan, genome]
        print("minisplice calibrate: " + " ".join(command), flush=True)
        _run(command, cali, hint)
    command = [executable, "predict", "-t", str(threads), "-c", cali, kan, genome]
    print("minisplice: " + " ".join(command), flush=True)
    _run(command, output, hint)
    write_json(cache_path, fingerprint)
    return output


def _open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


class SpliceIndex:
    def __init__(self):
        self.sites = {}

    @staticmethod
    def donor_key(chrom, intron, strand):
        start, end = intron
        return (chrom, start - 1, "+", "D") if strand == "+" else (chrom, end, "-", "D")

    @staticmethod
    def acceptor_key(chrom, intron, strand):
        start, end = intron
        return (chrom, end, "+", "A") if strand == "+" else (chrom, start - 1, "-", "A")

    def load(self, path, models):
        targets = set()
        for model in models:
            for intron in model.cds_introns:
                targets.add(self.donor_key(model.chrom, intron, model.strand))
                targets.add(self.acceptor_key(model.chrom, intron, model.strand))
        with _open_text(path) as handle:
            for line in handle:
                fields = line.split("\t")
                if len(fields) < 5:
                    continue
                key = (fields[0], int(fields[1]), fields[2], fields[3])
                if key in targets:
                    score = float(fields[4])
                    if score > self.sites.get(key, float("-inf")):
                        self.sites[key] = score
        return len(self.sites)

    def __len__(self):
        return len(self.sites)

    def intron_scores(self, model, intron):
        donor = self.sites.get(self.donor_key(model.chrom, intron, model.strand))
        acceptor = self.sites.get(self.acceptor_key(model.chrom, intron, model.strand))
        return donor, acceptor
