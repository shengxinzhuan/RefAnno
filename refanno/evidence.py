import bisect
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

from .audit import executable_version, file_sha256, read_json, write_json
from .models import (
    TX_TYPES,
    IsoChain,
    Model,
    ProteinHit,
    attrs,
    blocks_coverage,
    reverse_complement,
)


class Genome:
    def __init__(self, path):
        self.seq = {}
        name = None
        with open(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    self.seq[name] = []
                elif name:
                    self.seq[name].append(line.strip().upper())
        self.seq = {name: "".join(parts) for name, parts in self.seq.items()}

    def cds_sequence(self, model):
        if model.chrom not in self.seq:
            raise ValueError(f"基因组缺少 {model.chrom}")
        seq = "".join(self.seq[model.chrom][s - 1:e] for s, e in model.cds)
        return reverse_complement(seq) if model.strand == "-" else seq

    def splice_strand(self, chrom, intron_list):
        if chrom not in self.seq:
            return ""
        votes = set()
        seq = self.seq[chrom]
        for start, end in intron_list:
            if start < 1 or end > len(seq) or end < start + 3:
                return ""
            motif = seq[start - 1:start + 1] + seq[end - 2:end]
            if motif in {"GTAG", "GCAG", "ATAC"}:
                votes.add("+")
            elif motif in {"CTAC", "CTGC", "GTAT"}:
                votes.add("-")
            else:
                votes.add("?")
        return votes.pop() if len(votes) == 1 else ""


def parse_gff(path, set_name):
    gene_info = {}
    tx_parent = {}
    children = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, kind, start, end, _, strand, _, attr = fields
            a = attrs(attr)
            if kind in {"gene", "pseudogene"} and a.get("ID"):
                gene_info[a["ID"]] = (chrom, strand, a.get("gene_biotype", a.get("biotype", kind)))
            elif kind in TX_TYPES and a.get("ID"):
                tx_parent[a["ID"]] = a.get("Parent", a["ID"]).split(",")[0]
            elif kind in {"exon", "CDS"}:
                parents = a.get("Parent", "").split(",")
                for parent in filter(None, parents):
                    children.append((parent, kind, chrom, int(start), int(end), strand))

    models = {}
    gene_for_tx = dict(tx_parent)
    for parent, kind, chrom, start, end, strand in children:
        tx_id = parent
        gene_id = gene_for_tx.get(parent)
        if gene_id is None:
            if parent in gene_info:
                gene_id = parent
                tx_id = parent + ".t1"
            else:
                gene_id = re.sub(r"\.t\d+$", "", parent)
            gene_for_tx[tx_id] = gene_id
        key = (gene_id, tx_id)
        if key not in models:
            biotype = gene_info.get(gene_id, (chrom, strand, "protein_coding"))[2]
            models[key] = Model(set_name, gene_id, tx_id, chrom, strand, biotype=biotype)
        target = models[key].exons if kind == "exon" else models[key].cds
        target.append((start, end))

    out = []
    for model in models.values():
        model.finish()
        if model.cds:
            out.append(model)
    return out


class SJIndex:
    def __init__(self):
        self.data = defaultdict(lambda: defaultdict(lambda: [0, 0]))

    def load(self, path):
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) < 8:
                    continue
                chrom, start, end, code = fields[:4]
                strand = {"0": ".", "1": "+", "2": "-"}.get(code, code)
                dest = self.data[chrom][(int(start), int(end), strand)]
                dest[0] += int(fields[6])
                dest[1] += int(fields[7])

    def add(self, chrom, intron, strand, unique=0, multi=0):
        dest = self.data[chrom][(intron[0], intron[1], strand)]
        dest[0] += unique
        dest[1] += multi

    def get(self, model, intron):
        found = [0, 0]
        for strand in (model.strand, "."):
            value = self.data[model.chrom].get((intron[0], intron[1], strand))
            if value:
                found[0] += value[0]
                found[1] += value[1]
        return found


def load_bam_junctions(paths, genome, sj, depth=None, min_mapq=10, min_anchor=8):
    stats = []
    if not paths:
        return stats
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("读取 BAM 需要 pysam：conda install -c bioconda pysam") from exc
    for path in paths:
        counts = {"path": os.path.abspath(path), "reads": 0, "mapq_pass": 0, "spliced_reads": 0}
        stats.append(counts)
        with pysam.AlignmentFile(path, "rb") as bam:
            if depth is not None and not set(bam.references) & set(depth.by_chrom):
                raise RuntimeError(f"{path} 与注释没有共同的序列名，请检查 seqid 是否一致")
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                counts["reads"] += 1
                if depth is not None:
                    depth.consume(read)
                if read.mapping_quality < min_mapq:
                    continue
                counts["mapq_pass"] += 1
                cigar = read.cigartuples
                if not any(op == 3 for op, _ in (cigar or ())):
                    continue
                counts["spliced_reads"] += 1
                xs = read.get_tag("XS") if read.has_tag("XS") else "."
                multi = read.get_tag("NH") > 1 if read.has_tag("NH") else False
                ref = read.reference_start + 1
                anchor = 0
                pending = None
                for op, length in cigar:
                    if op in {0, 7, 8}:
                        ref += length
                        anchor += length
                        if pending and anchor >= min_anchor:
                            intron, left_anchor = pending
                            strand = genome.splice_strand(read.reference_name, [intron]) or xs
                            if strand not in {"+", "-"}:
                                strand = "."
                            sj.add(read.reference_name, intron, strand, multi=int(multi), unique=int(not multi))
                            pending = None
                    elif op == 2:
                        ref += length
                    elif op == 3:
                        pending = ((ref, ref + length - 1), anchor) if anchor >= min_anchor and 5 <= length <= 200000 else None
                        ref += length
                        anchor = 0
        if counts["reads"] >= 10000 and counts["spliced_reads"] < 0.01 * counts["mapq_pass"]:
            print(f"警告: {path} 的 {counts['mapq_pass']} 条比对中仅 {counts['spliced_reads']} 条含剪接"
                  f"（{100 * counts['spliced_reads'] / max(counts['mapq_pass'], 1):.2f}%），可能不是 RNA-seq"
                  "（Hi-C/ATAC/基因组测序？）；其覆盖度会按 RNA 深度参与裁决，请确认是否应作为 --bam 输入",
                  file=sys.stderr, flush=True)
    if depth is not None:
        depth.finalize()
    return stats


class ProteinIndex:
    def __init__(self, path="", protein_lengths=None):
        self.by_chrom = defaultdict(list)
        self.starts = {}
        self.max_span = 0
        self.protein_lengths = protein_lengths or {}
        if path:
            self.load(path)

    def __len__(self):
        return sum(len(values) for values in self.by_chrom.values())

    def load(self, path):
        hits = {}
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                a = attrs(fields[8])
                if fields[2] == "mRNA" and a.get("ID"):
                    target = a.get("Target", "").split()
                    if len(target) < 3:
                        continue
                    hit = ProteinHit(a["ID"], fields[0], fields[6], int(fields[3]), int(fields[4]),
                                     float(a.get("Identity", 0)), target[0], int(target[1]), int(target[2]))
                    hits[hit.hit_id] = hit
                    self.by_chrom[hit.chrom].append(hit)
                    self.max_span = max(self.max_span, hit.end - hit.start + 1)
                elif fields[2] == "CDS":
                    for parent in a.get("Parent", "").split(","):
                        if parent in hits:
                            hits[parent].cds.append((int(fields[3]), int(fields[4])))
        for chrom, values in self.by_chrom.items():
            values.sort(key=lambda h: (h.start, h.end, h.hit_id))
            self.starts[chrom] = [h.start for h in values]

    def best(self, model, min_identity=0.0):
        values = self.by_chrom.get(model.chrom, [])
        starts = self.starts.get(model.chrom, [])
        lo = bisect.bisect_left(starts, model.cds[0][0] - self.max_span)
        hi = bisect.bisect_right(starts, model.cds[-1][1])
        best = None
        for hit in values[lo:hi]:
            if hit.end < model.cds[0][0] or hit.strand != model.strand or not hit.cds or hit.identity < min_identity:
                continue
            coverage = min(1.0, blocks_coverage(hit.cds, model.cds) / max(model.cds_len, 1))
            query_coverage = ((hit.query_end - hit.query_start + 1) / self.protein_lengths[hit.target]
                              if hit.target in self.protein_lengths and self.protein_lengths[hit.target] else 0.0)
            candidate = (coverage, query_coverage, hit.identity, hit)
            if best is None or candidate[:3] > best[:3]:
                best = candidate
        return best


def load_protein_lengths(path):
    lengths = {}
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
    return lengths


def run_miniprot(genome, proteins, executable, threads, outdir):
    os.makedirs(outdir, exist_ok=True)
    output = os.path.join(outdir, "miniprot.gff3")
    cache_path = os.path.join(outdir, "miniprot.cache.json")
    fingerprint = {
        "genome_sha256": file_sha256(genome), "proteins_sha256": file_sha256(proteins),
        "miniprot_version": executable_version(executable), "threads": threads,
        "options": ["--gff", "--outn=1"],
    }
    if os.path.isfile(output) and os.path.getsize(output) > 0 and os.path.isfile(cache_path):
        try:
            if read_json(cache_path) == fingerprint:
                print(f"miniprot: reuse {output}", flush=True)
                return output
        except (OSError, ValueError):
            pass
    command = [executable, "-t", str(threads), "--gff", "--outn=1", genome, proteins]
    print("miniprot: " + " ".join(command), flush=True)
    temporary = output + ".tmp"
    try:
        with open(temporary, "w") as handle:
            subprocess.run(command, stdout=handle, check=True)
        os.replace(temporary, output)
        write_json(cache_path, fingerprint)
    except FileNotFoundError as exc:
        raise RuntimeError("找不到 miniprot；请安装后用 --miniprot 指定可执行文件路径") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"miniprot 运行失败，退出码 {exc.returncode}") from exc
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)
    return output


class IsoIndex:
    def __init__(self, chains):
        self.by_chrom = defaultdict(list)
        self.junctions = Counter()
        self.max_span = 0
        seen = set()
        for chain in chains:
            key = (chain.chrom, chain.read_id)
            if not chain.exons or key in seen:
                continue
            seen.add(key)
            self.by_chrom[chain.chrom].append(chain)
            self.max_span = max(self.max_span, chain.end - chain.start + 1)
            for intron in chain.introns:
                self.junctions[(chain.chrom, chain.strand, intron)] += 1
        self.starts = {}
        for chrom, values in self.by_chrom.items():
            values.sort(key=lambda c: (c.start, c.end, c.read_id))
            self.starts[chrom] = [c.start for c in values]

    def __len__(self):
        return sum(len(values) for values in self.by_chrom.values())

    def overlapping(self, chrom, start, end):
        values = self.by_chrom.get(chrom, [])
        starts = self.starts.get(chrom, [])
        lo = bisect.bisect_left(starts, start - self.max_span)
        hi = bisect.bisect_right(starts, end)
        return [c for c in values[lo:hi] if c.end >= start]

    def exact_support(self, model, min_coverage):
        if not model.cds:
            return 0
        count = 0
        target = model.cds_introns
        cds_span = (model.cds[0][0], model.cds[-1][1])
        for chain in self.overlapping(model.chrom, *cds_span):
            if chain.strand != model.strand or chain.introns != target:
                continue
            if blocks_coverage(chain.exons, model.cds) >= min_coverage * model.cds_len:
                count += 1
        return count

    def junction_support(self, model, intron):
        return self.junctions.get((model.chrom, model.strand, intron), 0)

    def utr_compatible(self, model, min_coverage):
        if not model.cds:
            return []
        out = []
        cds_span = (model.cds[0][0], model.cds[-1][1])
        model_introns = set(model.cds_introns)
        for chain in self.overlapping(model.chrom, *cds_span):
            if chain.strand != model.strand or not model_introns.issubset(chain.introns):
                continue
            if any(s <= ce and e >= cs for cs, ce in model.cds for s, e in chain.introns):
                continue
            if blocks_coverage(chain.exons, model.cds) >= min_coverage * model.cds_len:
                out.append(chain)
        return out


def load_bed12(path):
    chains = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip().split("\t")
            if len(p) < 12 or p[5] not in {"+", "-"}:
                continue
            chrom_start = int(p[1])
            sizes = [int(x) for x in p[10].rstrip(",").split(",")]
            offsets = [int(x) for x in p[11].rstrip(",").split(",")]
            exons = [(chrom_start + off + 1, chrom_start + off + size)
                     for size, off in zip(sizes, offsets) if size > 0]
            if exons:
                chains.append(IsoChain(p[3], p[0], p[5], sorted(exons), int(p[4])))
    return chains


def load_iso_bam(path, genome, min_mapq=10, min_anchor=20):
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("读取 BAM 需要 pysam：conda install -c bioconda pysam") from exc
    chains = []
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                continue
            exons, read_introns, valid = [], [], True
            ref = read.reference_start + 1
            start = None
            anchor = 0
            pending = None
            anchors = []
            for op, length in read.cigartuples or []:
                if op in {0, 7, 8}:
                    if start is None:
                        start = ref
                    if pending is not None:
                        anchors[pending][1] = length
                        pending = None
                    ref += length
                    anchor += length
                elif op == 2:
                    ref += length
                elif op == 3:
                    if start is None or not 5 <= length <= 200000:
                        valid = False
                        break
                    exons.append((start, ref - 1))
                    read_introns.append((ref, ref + length - 1))
                    anchors.append([anchor, -1])
                    pending = len(anchors) - 1
                    ref += length
                    start = None
                    anchor = 0
            if start is not None:
                exons.append((start, ref - 1))
            if not valid or not exons or any(left < min_anchor or right < min_anchor for left, right in anchors):
                continue
            strand = genome.splice_strand(read.reference_name, read_introns)
            if not read_introns:
                try:
                    strand = read.get_tag("ts")
                except KeyError:
                    try:
                        strand = read.get_tag("XS")
                    except KeyError:
                        strand = ""
            if strand in {"+", "-"}:
                chains.append(IsoChain(read.query_name, read.reference_name, strand, exons, read.mapping_quality))
    return chains


class DepthIndex:
    """BAM 顺序扫描时按 CDS 块累计深度；唯一/多重比对分开统计，扫描结束后一次求和。"""

    def __init__(self, paths, models=()):
        self.bams = list(paths)
        self.blocks = {}
        self.by_chrom = {}
        self.index = None
        self.results = {}
        self.cache = {}
        for model in models:
            for start, end in model.cds:
                self.add_block(model.chrom, start, end)

    def add_block(self, chrom, start, end):
        key = (chrom, start, end)
        if key not in self.blocks:
            size = end - start + 1
            self.blocks[key] = ([0] * (size + 1), [0] * (size + 1))
            self.by_chrom.setdefault(chrom, []).append(key)

    def build_index(self):
        self.index = {}
        for chrom, keys in self.by_chrom.items():
            keys.sort()
            starts = [key[1] - 1 for key in keys]
            longest = max(key[2] - key[1] + 1 for key in keys)
            self.index[chrom] = (starts, keys, longest)

    def consume(self, read):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return
        if self.index is None:
            self.build_index()
        entry = self.index.get(read.reference_name)
        if entry is None:
            return
        blocks = read.get_blocks()
        if not blocks:
            return
        starts, keys, longest = entry
        lo = bisect.bisect_left(starts, blocks[0][0] - longest)
        hi = bisect.bisect_right(starts, blocks[-1][1] - 1)
        if lo >= hi:
            return
        multi = int(read.has_tag("NH") and read.get_tag("NH") > 1)
        for block_start, block_end in blocks:
            for k in range(lo, hi):
                key = keys[k]
                offset = key[1] - 1
                left = block_start if block_start > offset else offset
                right = block_end if block_end < key[2] else key[2]
                if left < right:
                    diff = self.blocks[key][multi]
                    diff[left - offset] += 1
                    diff[right - offset] -= 1

    def finalize(self):
        for key, (unique, multi) in self.blocks.items():
            self.results[key] = (summarize_depth(unique), summarize_depth(multi))
        self.blocks = {}

    def get(self, model):
        key = (model.chrom, tuple(model.cds))
        if key in self.cache:
            return self.cache[key]
        if not self.bams:
            return 0.0, 0.0, 0.0, 0.0
        covered = [0, 0]
        total_depth = [0, 0]
        for start, end in model.cds:
            result = self.results.get((model.chrom, start, end))
            if result is None:
                continue
            for klass in (0, 1):
                covered[klass] += result[klass][0]
                total_depth[klass] += result[klass][1]
        scale = max(model.cds_len, 1)
        result = (covered[0] / scale, total_depth[0] / scale, covered[1] / scale, total_depth[1] / scale)
        self.cache[key] = result
        return result


def summarize_depth(diff):
    current = 0
    covered = 0
    total = 0
    for value in diff[:-1]:
        current += value
        total += current
        covered += current > 0
    return covered, total
