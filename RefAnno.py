#!/usr/bin/env python3
"""RefAnno: multi-evidence gene structure refinement and representative annotation."""

import argparse
import bisect
import csv
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from urllib.parse import quote


TX_TYPES = {"mRNA", "transcript", "lncRNA", "ncRNA", "rRNA", "tRNA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def introns(blocks):
    blocks = sorted(blocks)
    return tuple((blocks[i - 1][1] + 1, blocks[i][0] - 1)
                 for i in range(1, len(blocks)) if blocks[i][0] > blocks[i - 1][1] + 1)


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def blocks_coverage(blocks, targets):
    return sum(overlap(c, e) for c in blocks for e in targets)


def attrs(text):
    text = text.strip()
    if '"' in text and "=" not in text.split(";")[0]:
        out = {}
        for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
            out[key] = value
        return out
    out = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def attr_text(items):
    return ";".join(f"{key}={quote(str(value), safe=':|+_,.-')}"
                    for key, value in items if value is not None and value != "")


def clean_id(value):
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)


@dataclass
class Model:
    set_name: str
    gene_id: str
    tx_id: str
    chrom: str
    strand: str
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)
    biotype: str = "protein_coding"
    source_weight: float = 1.0
    structure_coverage: float = 0.0
    structure_status: str = "untested"
    junctions_supported: int = 0
    canonical_junctions: int = 0
    complete_cds: bool = False
    rna_supported: bool = False
    protein_supported: bool = False
    depth: float = 0.0
    coverage: float = 0.0
    multi_depth: float = 0.0
    multi_coverage: float = 0.0
    family_support: str = ""
    iso_support: int = 0
    protein_coverage: float = 0.0
    protein_identity: float = 0.0
    protein_query_coverage: float = 0.0
    protein_target: str = ""
    output_id: str = ""
    utr_note: str = ""
    review: list = field(default_factory=list)

    @property
    def key(self):
        return f"{self.set_name}|{self.gene_id}|{self.tx_id}"

    @property
    def source_model(self):
        return f"{self.set_name}|{self.tx_id}"

    @property
    def cds_introns(self):
        return introns(self.cds)

    @property
    def cds_len(self):
        return sum(e - s + 1 for s, e in self.cds)

    @property
    def junctions_total(self):
        return len(self.cds_introns)

    def finish(self):
        self.exons = sorted(set(self.exons or self.cds))
        self.cds = sorted(set(self.cds))


@dataclass
class IsoChain:
    read_id: str
    chrom: str
    strand: str
    exons: list
    mapq: int = 60

    @property
    def introns(self):
        return introns(self.exons)

    @property
    def start(self):
        return self.exons[0][0]

    @property
    def end(self):
        return self.exons[-1][1]


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


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


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


def load_bam_junctions(paths, genome, sj, min_mapq=10, min_anchor=8):
    if not paths:
        return
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("读取 BAM 需要 pysam：conda install -c bioconda pysam") from exc
    for path in paths:
        with pysam.AlignmentFile(path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                    continue
                try:
                    xs = read.get_tag("XS")
                except KeyError:
                    xs = "."
                multi = read.get_tag("NH") > 1 if read.has_tag("NH") else False
                ref = read.reference_start + 1
                anchor = 0
                pending = None
                for op, length in read.cigartuples or []:
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


@dataclass
class ProteinHit:
    hit_id: str
    chrom: str
    strand: str
    start: int
    end: int
    identity: float
    target: str
    query_start: int
    query_end: int
    cds: list = field(default_factory=list)


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
    if os.path.exists(output) and os.path.getsize(output) > 0:
        print(f"miniprot: reuse {output}", flush=True)
        return output
    command = [executable, "-t", str(threads), "--gff", "--outn=1", genome, proteins]
    print("miniprot: " + " ".join(command), flush=True)
    try:
        with open(output, "w") as handle:
            subprocess.run(command, stdout=handle, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError("找不到 miniprot；请安装后用 --miniprot 指定可执行文件路径") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"miniprot 运行失败，退出码 {exc.returncode}") from exc
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
    def __init__(self, paths):
        self.cache = {}
        if not paths:
            self.bams = []
            return
        try:
            import pysam
        except ImportError as exc:
            raise RuntimeError("--bam 需要 pysam；也可改用 --sj 与 BED12，纯标准库运行") from exc
        self.bams = [pysam.AlignmentFile(path, "rb") for path in paths]

    def close(self):
        for bam in self.bams:
            bam.close()

    def get(self, model):
        key = (model.chrom, tuple(model.cds))
        if key in self.cache:
            return self.cache[key]
        if not self.bams:
            return 0.0, 0.0, 0.0, 0.0
        covered = [0, 0]
        total_depth = [0, 0]
        for start, end in model.cds:
            size = end - start + 1
            diffs = [[0] * (size + 1), [0] * (size + 1)]
            for bam in self.bams:
                try:
                    reads = bam.fetch(model.chrom, start - 1, end)
                except ValueError as exc:
                    raise RuntimeError(f"BAM 无法按坐标查询，请先 samtools index：{exc}") from exc
                for read in reads:
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    multi = read.has_tag("NH") and read.get_tag("NH") > 1
                    diff = diffs[int(multi)]
                    for block_start, block_end in read.get_blocks():
                        left = max(start - 1, block_start)
                        right = min(end, block_end)
                        if left < right:
                            diff[left - start + 1] += 1
                            diff[right - start + 1] -= 1
            for index, diff in enumerate(diffs):
                current = 0
                for value in diff[:size]:
                    current += value
                    total_depth[index] += current
                    covered[index] += current > 0
        result = (covered[0] / max(model.cds_len, 1), total_depth[0] / max(model.cds_len, 1),
                  covered[1] / max(model.cds_len, 1), total_depth[1] / max(model.cds_len, 1))
        self.cache[key] = result
        return result


def reciprocal_cds_overlap(a, b):
    value = blocks_coverage(a.cds, b.cds)
    return value / max(a.cds_len, 1), value / max(b.cds_len, 1)


def same_chain(a, b):
    if a.chrom != b.chrom or a.strand != b.strand or len(a.cds) != len(b.cds):
        return False
    if len(a.cds) == 1:
        return overlap(a.cds[0], b.cds[0]) / max(1, min(a.cds_len, b.cds_len)) >= 0.9
    return a.cds_introns == b.cds_introns


def cluster_loci(models):
    parent = {model.key: model.key for model in models}
    by_key = {model.key: model for model in models}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_chrom = defaultdict(list)
    for model in models:
        for block in model.cds:
            by_chrom[model.chrom].append((block[0], block[1], model.key))
    for intervals in by_chrom.values():
        active = []
        for start, end, key in sorted(intervals):
            active = [item for item in active if item[0] >= start]
            for _, other in active:
                union(key, other)
            active.append((end, key))
    groups = defaultdict(list)
    for key in parent:
        groups[find(key)].append(by_key[key])
    return list(groups.values())


def evaluate_models(models, sj, iso, proteins, depth, genome, weights, min_sj, min_iso, min_coverage, min_protein_identity):
    protein_cache = {}
    has_rna = bool(depth.bams or len(iso))
    for model in models:
        model.source_weight = weights.get(model.set_name.lower(), 1.0)
        model.coverage, model.depth, model.multi_coverage, model.multi_depth = depth.get(model)
        model.iso_support = iso.exact_support(model, min_coverage)
        for intron in model.cds_introns:
            unique, _ = sj.get(model, intron)
            iso_n = iso.junction_support(model, intron)
            if unique >= min_sj or iso_n >= min_iso:
                model.junctions_supported += 1
            if genome.splice_strand(model.chrom, [intron]) == model.strand:
                model.canonical_junctions += 1
        if model.cds_introns:
            model.structure_coverage = model.junctions_supported / model.junctions_total
        else:
            model.structure_coverage = 1.0 if model.iso_support >= min_iso else model.coverage
        model.rna_supported = bool(
            model.iso_support >= min_iso
            or (model.cds_introns and model.junctions_supported)
            or (not model.cds_introns and model.coverage >= 0.5 and model.depth >= 10)
        )

        protein_key = (model.chrom, model.strand, tuple(model.cds))
        if protein_key not in protein_cache:
            seq = genome.cds_sequence(model)
            protein = translate(seq)
            complete = len(seq) % 3 == 0 and seq.startswith("ATG") and protein.endswith("*") and "*" not in protein[:-1]
            protein_cache[protein_key] = (complete, "*" in protein[:-1])
        model.complete_cds, internal_stop = protein_cache[protein_key]
        if internal_stop:
            model.review.append("internal_stop")
        protein_hit = proteins.best(model, min_protein_identity)
        if protein_hit:
            model.protein_coverage, model.protein_query_coverage, model.protein_identity, hit = protein_hit
            model.protein_target = hit.target

        if model.iso_support >= min_iso:
            model.structure_status = "exact_isoseq"
        elif not has_rna:
            model.structure_status = "untested"
        elif model.cds_introns and model.junctions_supported == model.junctions_total:
            model.structure_status = "supported"
        elif not model.cds_introns and model.coverage >= 0.5 and model.depth >= 10:
            model.structure_status = "supported"
        elif model.junctions_supported or model.iso_support or model.depth > 0:
            model.structure_status = "partial"
        else:
            model.structure_status = "unsupported"


def model_has_expression(model, sj, iso, args):
    if model.iso_support >= args.min_isoform_molecules:
        return True
    if not model.cds_introns:
        return model.coverage >= 0.5 and model.depth >= 10
    for intron in model.cds_introns:
        unique, _ = sj.get(model, intron)
        if unique >= args.min_sj or iso.junction_support(model, intron) >= args.min_isoform_molecules:
            return True
    return False


def model_has_protein(model, proteins, args):
    if "internal_stop" in model.review:
        return False
    if model.protein_identity < args.min_protein_identity or model.protein_coverage < args.min_protein_coverage:
        return False
    return (not proteins.protein_lengths or
            model.protein_query_coverage >= args.min_protein_query_coverage)


def chain_groups(locus):
    groups = []
    for model in sorted(locus, key=lambda m: m.key):
        for group in groups:
            if same_chain(model, group[0]):
                group.append(model)
                break
        else:
            groups.append([model])
    return groups


def descending_rank(metrics, tie):
    return tuple(-float(value) for value in metrics) + (tie,)


def model_rank(model):
    structure_rank = {"exact_isoseq": 3, "supported": 2, "partial": 1}.get(model.structure_status, 0)
    metrics = (structure_rank, model.complete_cds, model.structure_coverage, model.iso_support,
               model.depth, model.coverage, model.protein_coverage, model.protein_query_coverage,
               model.protein_identity, model.source_weight)
    return descending_rank(metrics, model.key)


def best_in_group(group):
    return sorted(group, key=model_rank)[0]


def representative_basis(item, args):
    model, _, reason, _ = item
    if model.iso_support >= args.min_isoform_molecules:
        return 5, "isoseq_chain_abundance"
    if model.rna_supported:
        return 4, "rna_structure_depth"
    if model.protein_supported:
        return 3, "protein_coverage"
    if reason == "tandem_family_support":
        return 2, "family_support"
    if reason.endswith("set_chain_consensus"):
        return 1, "input_consensus"
    return 0, "trusted_source"


def representative_rank(item, args):
    model, _, _, sources = item
    evidence_rank, _ = representative_basis(item, args)
    protein_metrics = (model.protein_coverage, model.protein_query_coverage, model.protein_identity) if evidence_rank == 3 else (0, 0, 0)
    metrics = (evidence_rank, model.complete_cds, model.structure_coverage, model.iso_support,
               model.depth, model.coverage) + protein_metrics + (len(sources), model.source_weight)
    return descending_rank(metrics, model.key)


def choose_representatives(loci, sj, iso, proteins, args, trusted, tandem_info):
    winners = []
    dropped = {}
    for locus in loci:
        kept = []
        for group in chain_groups(locus):
            best = best_in_group(group)
            best.protein_supported = model_has_protein(best, proteins, args)
            family = tandem_info.get(best.key, {})
            sources = {m.set_name for m in group}
            if len(sources) >= 2:
                reason = f"{len(sources)}set_chain_consensus"
            elif best.set_name.lower() in trusted:
                reason = "trusted_unique"
            elif model_has_expression(best, sj, iso, args):
                reason = "unique_with_evidence"
            elif model_has_protein(best, proteins, args):
                reason = "unique_with_protein"
            elif family.get("family_support") and "internal_stop" not in best.review:
                reason = "tandem_family_support"
                best.family_support = family["family_support"]
                if "tandem_family_support" not in best.review:
                    best.review.append("tandem_family_support")
            else:
                for model in group:
                    dropped[model.key] = "unique_no_evidence"
                continue
            kept.append((best, group, reason, sources))
        kept.sort(key=lambda item: representative_rank(item, args))
        surviving = []
        for item in kept:
            best = item[0]
            if any(best.strand == other[0].strand and max(reciprocal_cds_overlap(best, other[0])) > args.conflict_overlap
                   for other in surviving):
                for model in item[1]:
                    dropped[model.key] = "locus_conflict_lower_rank"
            else:
                surviving.append(item)
        for best, group, reason, sources in surviving:
            for model in group:
                if model is not best:
                    dropped[model.key] = "same_chain_lower_rank"
            _, selection = representative_basis((best, group, reason, sources), args)
            winners.append({"model": best, "group": group, "reason": reason, "selection": selection,
                            "sources": sources, "locus": locus, "alternatives": []})
    return winners, dropped


def detect_duplicate_clusters(models, args):
    by_location = defaultdict(list)
    for model in models:
        by_location[(model.chrom, model.strand)].append(model)
    candidates = []
    for models in by_location.values():
        models.sort(key=lambda m: (m.cds[0][0], m.cds[-1][1], m.key))
        current = []
        for model in models:
            if current:
                previous = current[-1]
                gap = max(0, model.cds[0][0] - previous.cds[-1][1] - 1)
                length_ratio = max(model.cds_len, previous.cds_len) / max(1, min(model.cds_len, previous.cds_len))
                similar = (gap <= args.tandem_max_gap and length_ratio <= args.tandem_cds_ratio and
                           abs(len(model.cds) - len(previous.cds)) <= 1 and
                           max(reciprocal_cds_overlap(model, previous)) <= args.conflict_overlap)
                if similar:
                    current.append(model)
                else:
                    if len(current) >= 2:
                        candidates.append(current)
                    current = [model]
            else:
                current = [model]
        if len(current) >= 2:
            candidates.append(current)
    candidates.sort(key=lambda group: (natural_key(group[0].chrom), group[0].cds[0][0],
                                       group[0].strand, group[0].key))
    annotations = {}
    counters = Counter()
    for group in candidates:
        kind = "tandem" if len(group) >= args.tandem_min_members else "duplicate"
        counters[kind] += 1
        cluster_id = f"{kind.upper()}{counters[kind]:06d}"
        for model in group:
            annotations[model.key] = {"type": kind, "cluster": cluster_id, "size": len(group)}
    return annotations


def prepare_tandem_evidence(loci, sj, proteins, args):
    grouped_loci = [chain_groups(locus) for locus in loci]
    representatives = [best_in_group(group) for groups in grouped_loci for group in groups]
    annotations = detect_duplicate_clusters(representatives, args)
    cluster_members = defaultdict(list)
    for model in representatives:
        if model.key in annotations:
            cluster_members[annotations[model.key]["cluster"]].append(model)
    cluster_support = {}
    for cluster, members in cluster_members.items():
        support = []
        if any(model_has_protein(model, proteins, args) for model in members):
            support.append("protein")
        if any(model.multi_coverage >= 0.5 and model.multi_depth >= args.tandem_min_multi_depth
               for model in members):
            support.append("multi_depth")
        multi_junction = 0
        for model in members:
            for intron in model.cds_introns:
                _, multi = sj.get(model, intron)
                multi_junction = max(multi_junction, multi)
        if multi_junction >= 3 * args.min_sj:
            support.append("multi_junction")
        cluster_support[cluster] = "+".join(support)
    tandem_info = {}
    for groups in grouped_loci:
        for group in groups:
            annotation = annotations.get(best_in_group(group).key)
            if not annotation:
                continue
            annotation = dict(annotation)
            annotation["family_support"] = cluster_support[annotation["cluster"]]
            for model in group:
                tandem_info[model.key] = annotation
    return tandem_info


def consensus_exons(reads, node_window=3, end_window=50):
    clusters = []
    for read in sorted(reads, key=lambda r: (len(r.exons), r.exons, r.read_id)):
        for cluster in clusters:
            reference = cluster[0]
            if len(reference.exons) != len(read.exons):
                continue
            compatible = True
            for i, (left, right) in enumerate(zip(reference.exons, read.exons)):
                start_window = end_window if i == 0 else node_window
                stop_window = end_window if i == len(read.exons) - 1 else node_window
                if abs(left[0] - right[0]) > start_window or abs(left[1] - right[1]) > stop_window:
                    compatible = False
                    break
            if compatible:
                cluster.append(read)
                break
        else:
            clusters.append([read])
    if not clusters:
        return None, 0, 0.0
    best = max(clusters, key=lambda c: (len(c), -c[0].start, c[0].read_id))
    n = len(best[0].exons)
    medians = []
    for i in range(n):
        medians.append((sorted(r.exons[i][0] for r in best)[len(best) // 2],
                        sorted(r.exons[i][1] for r in best)[len(best) // 2]))
    medoid = min(best, key=lambda r: (sum(abs(r.exons[i][0] - medians[i][0]) +
                                             abs(r.exons[i][1] - medians[i][1]) for i in range(n)), r.read_id))
    return medoid, len(best), len(best) / len(reads)


def apply_utr_consensus(winners, iso, args):
    for winner in winners:
        model = winner["model"]
        reads = iso.utr_compatible(model, args.min_cds_coverage)
        if len(reads) < args.utr_min_molecules:
            continue
        best, support, dominance = consensus_exons(reads)
        if support < args.utr_min_molecules or dominance < args.utr_dominance:
            model.review.append("utr_consensus_ambiguous")
            continue
        if blocks_coverage(best.exons, model.cds) != model.cds_len:
            model.review.append("utr_consensus_incomplete")
            continue
        old = list(model.exons)
        model.exons = sorted(best.exons)
        model.utr_note = f"{best.read_id}:n={support}:frac={dominance:.2f}"
        winner["utr_changed"] = old != model.exons


def valid_alternative(model, sj, iso, genome, args):
    if not model.cds or not model.exons or model.cds_len % 3:
        return False
    if blocks_coverage(model.exons, model.cds) != model.cds_len:
        return False
    protein = translate(genome.cds_sequence(model))
    if "*" in protein[:-1]:
        return False
    for intron in model.cds_introns:
        if not 5 <= intron[1] - intron[0] + 1 <= 200000:
            return False
        motif = genome.splice_strand(model.chrom, [intron])
        unique, _ = sj.get(model, intron)
        if motif != model.strand and unique < args.min_sj and iso.junction_support(model, intron) < args.min_isoform_molecules:
            return False
    return True


def attach_isoforms(winners, sj, iso, genome, args):
    primary_keys = {w["model"].key for w in winners}
    for winner in winners:
        primary = winner["model"]
        candidates = []
        for model in winner["locus"]:
            if model.key in primary_keys or model is primary or model.strand != primary.strand:
                continue
            if model.cds_introns == primary.cds_introns:
                continue
            support = iso.exact_support(model, args.min_cds_coverage)
            if support < args.min_isoform_molecules:
                continue
            if not valid_alternative(model, sj, iso, genome, args):
                continue
            candidates.append((support, model))
        candidates.sort(key=lambda item: (-item[0],) + model_rank(item[1]))
        seen_chains = {primary.cds_introns}
        for support, model in candidates:
            if len(winner["alternatives"]) >= args.max_isoforms - 1:
                break
            if model.cds_introns in seen_chains:
                continue
            seen_chains.add(model.cds_introns)
            winner["alternatives"].append((model, support))


def utr_segments(exons, cds, strand):
    remaining = []
    for exon in sorted(exons):
        segments = [exon]
        for block in cds:
            next_segments = []
            for segment in segments:
                if overlap(segment, block) == 0:
                    next_segments.append(segment)
                    continue
                if segment[0] < block[0]:
                    next_segments.append((segment[0], block[0] - 1))
                if segment[1] > block[1]:
                    next_segments.append((block[1] + 1, segment[1]))
            segments = next_segments
        remaining.extend(segments)
    low, high = cds[0][0], cds[-1][1]
    five = [x for x in remaining if (x[1] < low if strand == "+" else x[0] > high)]
    three = [x for x in remaining if (x[0] > high if strand == "+" else x[1] < low)]
    return sorted(five), sorted(three)


def model_row(model, kind, start, end, attributes, phase="."):
    return "\t".join([model.chrom, "RefAnno", kind, str(start), str(end), ".",
                      model.strand, str(phase), attr_text(attributes)])


def write_transcript(handle, protein_handle, model, gene_id, rank, genome, role, iso_support=0):
    tx_id = f"{gene_id}-T{rank}"
    start, end = model.exons[0][0], model.exons[-1][1]
    items = [("ID", tx_id), ("Parent", gene_id), ("source_model", model.source_model),
             ("transcript_role", role)]
    if role == "alternative":
        items.extend([("isoform_molecules", iso_support), ("isoform_parent", f"{gene_id}-T1")])
    if model.utr_note:
        items.append(("utr_consensus", model.utr_note))
    handle.write(model_row(model, "mRNA", start, end, items) + "\n")
    for i, (start, end) in enumerate(model.exons, 1):
        handle.write(model_row(model, "exon", start, end,
                               [("ID", f"{tx_id}.exon.{i}"), ("Parent", tx_id)]) + "\n")
    ordered_cds = sorted(model.cds, reverse=model.strand == "-")
    elapsed = 0
    for i, (start, end) in enumerate(ordered_cds, 1):
        phase = (3 - elapsed % 3) % 3
        handle.write(model_row(model, "CDS", start, end,
                               [("ID", f"{tx_id}.CDS.{i}"), ("Parent", tx_id)], phase) + "\n")
        elapsed += end - start + 1
    five, three = utr_segments(model.exons, model.cds, model.strand)
    for kind, segments in (("five_prime_UTR", five), ("three_prime_UTR", three)):
        for i, (start, end) in enumerate(segments, 1):
            handle.write(model_row(model, kind, start, end,
                                   [("ID", f"{tx_id}.{kind}.{i}"), ("Parent", tx_id)]) + "\n")
    protein = translate(genome.cds_sequence(model))
    protein_handle.write(f">{tx_id} source_model={model.source_model} role={role}\n")
    for i in range(0, len(protein), 60):
        protein_handle.write(protein[i:i + 60] + "\n")


def natural_key(value):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", value)]


def write_outputs(winners, dropped, models, genome, args):
    winners.sort(key=lambda w: (natural_key(w["model"].chrom), w["model"].cds[0][0],
                                w["model"].cds[-1][1], w["model"].strand, w["model"].key))
    duplicate_info = detect_duplicate_clusters([w["model"] for w in winners], args)
    retained = set()
    transcript_count = 0
    iso_loci = 0
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "final.gff3"), "w") as gff, \
            open(os.path.join(args.outdir, "final.proteins.faa"), "w") as proteins, \
            open(os.path.join(args.outdir, "decision_report.tsv"), "w", newline="") as report:
        gff.write("##gff-version 3\n")
        fields = ["gene_id", "transcript_id", "role", "chrom", "start", "end", "strand", "source_model",
                  "structure_status", "structure_coverage", "junctions_supported", "junctions_total",
                  "canonical_junctions", "complete_cds", "cds_depth", "cds_coverage", "cds_multi_depth",
                  "cds_multi_coverage", "family_support", "isoform_molecules", "protein_target",
                  "protein_identity", "protein_coverage", "protein_query_coverage", "duplicate_type",
                  "duplicate_cluster", "duplicate_size", "support", "reason", "selection", "review"]
        writer = csv.DictWriter(report, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for index, winner in enumerate(winners, 1):
            model = winner["model"]
            gene_id = f"{args.id_prefix}{index:06d}"
            model.output_id = gene_id
            retained.add(model.key)
            alternatives = winner["alternatives"]
            if alternatives:
                iso_loci += 1
            gene_attrs = [("ID", gene_id), ("gene_biotype", "protein_coding"),
                          ("support", "+".join(sorted(winner["sources"]))),
                          ("source_model", model.source_model), ("reason", winner["reason"]),
                          ("selection", winner["selection"]), ("structure_status", model.structure_status),
                          ("structure_coverage", f"{model.structure_coverage:.3f}"),
                          ("junctions_supported", model.junctions_supported),
                          ("junctions_total", model.junctions_total),
                          ("canonical_junctions", model.canonical_junctions),
                          ("complete_cds", int(model.complete_cds)), ("cds_depth", f"{model.depth:.2f}"),
                          ("cds_coverage", f"{model.coverage:.3f}"),
                          ("cds_multi_depth", f"{model.multi_depth:.2f}"),
                          ("cds_multi_coverage", f"{model.multi_coverage:.3f}"),
                          ("family_support", model.family_support),
                          ("protein_coverage", f"{model.protein_coverage:.3f}"),
                          ("protein_identity", f"{model.protein_identity:.3f}"),
                          ("isoform_count", len(alternatives) + 1)]
            duplicate = duplicate_info.get(model.key)
            if duplicate:
                gene_attrs.extend([("duplicate_type", duplicate["type"]),
                                   ("duplicate_cluster", duplicate["cluster"]),
                                   ("duplicate_size", duplicate["size"])])
            if model.review:
                gene_attrs.append(("flag_review", ",".join(sorted(set(model.review)))))
            gff.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], gene_attrs) + "\n")
            write_transcript(gff, proteins, model, gene_id, 1, genome, "primary", model.iso_support)
            transcript_count += 1
            writer.writerow({"gene_id": gene_id, "transcript_id": f"{gene_id}-T1", "role": "primary",
                             "chrom": model.chrom, "start": model.exons[0][0], "end": model.exons[-1][1],
                             "strand": model.strand, "source_model": model.source_model,
                             "structure_status": model.structure_status,
                             "structure_coverage": f"{model.structure_coverage:.3f}",
                             "junctions_supported": model.junctions_supported,
                             "junctions_total": model.junctions_total,
                             "canonical_junctions": model.canonical_junctions,
                             "complete_cds": int(model.complete_cds),
                             "cds_depth": f"{model.depth:.3f}", "cds_coverage": f"{model.coverage:.3f}",
                             "cds_multi_depth": f"{model.multi_depth:.3f}",
                             "cds_multi_coverage": f"{model.multi_coverage:.3f}",
                             "family_support": model.family_support,
                             "isoform_molecules": model.iso_support, "protein_target": model.protein_target,
                             "protein_identity": f"{model.protein_identity:.3f}",
                             "protein_coverage": f"{model.protein_coverage:.3f}",
                             "protein_query_coverage": f"{model.protein_query_coverage:.3f}",
                             "duplicate_type": duplicate["type"] if duplicate else "",
                             "duplicate_cluster": duplicate["cluster"] if duplicate else "",
                             "duplicate_size": duplicate["size"] if duplicate else "",
                             "support": "+".join(sorted(winner["sources"])),
                             "reason": winner["reason"], "selection": winner["selection"],
                             "review": ",".join(sorted(set(model.review)))})
            for rank, (alternative, support) in enumerate(alternatives, 2):
                retained.add(alternative.key)
                alternative.output_id = gene_id
                write_transcript(gff, proteins, alternative, gene_id, rank, genome, "alternative", support)
                transcript_count += 1
                writer.writerow({"gene_id": gene_id, "transcript_id": f"{gene_id}-T{rank}", "role": "alternative",
                                 "chrom": alternative.chrom, "start": alternative.exons[0][0],
                                 "end": alternative.exons[-1][1], "strand": alternative.strand,
                                 "source_model": alternative.source_model,
                                 "structure_status": alternative.structure_status,
                                 "structure_coverage": f"{alternative.structure_coverage:.3f}",
                                 "junctions_supported": alternative.junctions_supported,
                                 "junctions_total": alternative.junctions_total,
                                 "canonical_junctions": alternative.canonical_junctions,
                                 "complete_cds": int(alternative.complete_cds),
                                 "cds_depth": f"{alternative.depth:.3f}", "cds_coverage": f"{alternative.coverage:.3f}",
                                 "cds_multi_depth": f"{alternative.multi_depth:.3f}",
                                 "cds_multi_coverage": f"{alternative.multi_coverage:.3f}",
                                 "family_support": alternative.family_support,
                                 "isoform_molecules": support, "protein_target": alternative.protein_target,
                                 "protein_identity": f"{alternative.protein_identity:.3f}",
                                 "protein_coverage": f"{alternative.protein_coverage:.3f}",
                                 "protein_query_coverage": f"{alternative.protein_query_coverage:.3f}",
                                 "duplicate_type": duplicate["type"] if duplicate else "",
                                 "duplicate_cluster": duplicate["cluster"] if duplicate else "",
                                 "duplicate_size": duplicate["size"] if duplicate else "",
                                 "support": alternative.set_name,
                                 "reason": "isoseq_exact_chain", "selection": "isoseq_chain_abundance",
                                 "review": ""})

    with open(os.path.join(args.outdir, "dropped.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        index = 0
        for model in sorted(models, key=lambda m: m.key):
            if model.key in retained:
                continue
            index += 1
            gid = f"D{index:06d}"
            reason = dropped.get(model.key, "not_selected")
            attrs0 = [("ID", gid), ("source_model", model.source_model), ("reason", reason)]
            out.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], attrs0) + "\n")
            txid = gid + "-T1"
            out.write(model_row(model, "mRNA", model.exons[0][0], model.exons[-1][1],
                                [("ID", txid), ("Parent", gid), ("source_model", model.source_model)]) + "\n")
            for i, (start, end) in enumerate(model.exons, 1):
                out.write(model_row(model, "exon", start, end, [("ID", f"{txid}.exon.{i}"), ("Parent", txid)]) + "\n")
            for i, (start, end) in enumerate(model.cds, 1):
                out.write(model_row(model, "CDS", start, end, [("ID", f"{txid}.CDS.{i}"), ("Parent", txid)]) + "\n")

    summary = {
        "input_transcripts": len(models), "genes": len(winners), "transcripts": transcript_count,
        "genes_with_alternatives": iso_loci, "alternative_transcripts": transcript_count - len(winners),
        "dropped_transcripts": len(models) - len(retained),
        "utr_consensus_applied": sum(bool(w["model"].utr_note) for w in winners),
        "utr_consensus_changed": sum(bool(w.get("utr_changed")) for w in winners),
        "protein_supported_transcripts": sum(m.protein_coverage > 0 for m in models),
        "protein_only_retained_genes": sum(w["reason"] == "unique_with_protein" for w in winners),
        "structure_supported_genes": sum(w["model"].structure_status in {"exact_isoseq", "supported"} for w in winners),
        "structure_partial_genes": sum(w["model"].structure_status == "partial" for w in winners),
        "structure_untested_genes": sum(w["model"].structure_status == "untested" for w in winners),
        "representative_basis": dict(Counter(w["selection"] for w in winners)),
        "tandem_clusters": len({d["cluster"] for d in duplicate_info.values() if d["type"] == "tandem"}),
        "tandem_genes": sum(d["type"] == "tandem" for d in duplicate_info.values()),
        "duplicate_clusters": len({d["cluster"] for d in duplicate_info.values() if d["type"] == "duplicate"}),
        "duplicate_genes": sum(d["type"] == "duplicate" for d in duplicate_info.values()),
        "tandem_family_retained_genes": sum(w["reason"] == "tandem_family_support" for w in winners),
        "max_isoforms": args.max_isoforms,
    }
    with open(os.path.join(args.outdir, "summary.json"), "w") as out:
        json.dump(summary, out, ensure_ascii=False, indent=2)
        out.write("\n")
    return summary


def named_paths(specs):
    out = []
    for spec in specs:
        if "=" not in spec:
            raise argparse.ArgumentTypeError(f"需要 NAME=PATH：{spec}")
        out.append(tuple(spec.split("=", 1)))
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="示例：\n  python3 RefAnno.py --set annevo=a.gff3 --set helixer=h.gff3 "
                                            "--genome genome.fa --proteins related.faa --bam RNA.bam "
                                            "--outdir RefAnno_out")
    parser.add_argument("--version", action="version", version="RefAnno 1.0.0")
    parser.add_argument("--set", action="append", required=True, metavar="NAME=GFF", help="输入注释，可多次")
    parser.add_argument("--genome", required=True, help="基因组 FASTA")
    parser.add_argument("--sj", action="append", default=[], help=argparse.SUPPRESS)
    parser.add_argument("--bam", action="append", default=[], help="二代坐标排序 BAM，自动提取 junction 并计算 CDS 深度，可多次")
    parser.add_argument("--proteins", help="近缘同源蛋白 FASTA；提供时自动运行 miniprot")
    parser.add_argument("--miniprot-gff", help=argparse.SUPPRESS)
    parser.add_argument("--miniprot", default="miniprot", help="miniprot 可执行文件")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-protein-identity", type=float, default=0.3)
    parser.add_argument("--min-protein-coverage", type=float, default=0.5, help="单条蛋白命中对模型 CDS 的最低覆盖")
    parser.add_argument("--min-protein-query-coverage", type=float, default=0.5, help="命中对查询蛋白的最低覆盖")
    parser.add_argument("--isoseq", action="append", default=[], help="可选：Iso-seq BAM 或 BED12 读链，可多次")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--id-prefix", default="REFA")
    parser.add_argument("--max-isoforms", type=int, default=3, help="每基因总转录本上限，含主代表")
    parser.add_argument("--min-isoform-molecules", type=int, default=2)
    parser.add_argument("--min-cds-coverage", type=float, default=0.9)
    parser.add_argument("--min-sj", type=int, default=3)
    parser.add_argument("--conflict-overlap", type=float, default=0.3)
    parser.add_argument("--tandem-max-gap", type=int, default=20000, help="串联重复相邻基因最大间距")
    parser.add_argument("--tandem-cds-ratio", type=float, default=1.8, help="串联重复相邻基因最大 CDS 长度比")
    parser.add_argument("--tandem-min-members", type=int, default=3, help="tandem 簇最少成员数；默认 2 成员标 duplicate")
    parser.add_argument("--tandem-min-multi-depth", type=float, default=5.0, help="串联家族多重比对 CDS 深度下限")
    parser.add_argument("--trusted-set", action="append", default=[], help="允许无 RNA 证据单独保留的集合名")
    parser.add_argument("--set-weight", action="append", default=[], metavar="NAME=FLOAT")
    parser.add_argument("--utr-consensus", action="store_true", help="用 Iso-seq 局部主导链替换 UTR 外显子，CDS 不变")
    parser.add_argument("--utr-min-molecules", type=int, default=2)
    parser.add_argument("--utr-dominance", type=float, default=0.6)
    args = parser.parse_args(argv)
    if args.max_isoforms < 1 or args.max_isoforms > 3:
        parser.error("--max-isoforms 必须为 1–3")
    thresholds = [args.min_cds_coverage, args.utr_dominance, args.conflict_overlap,
                  args.min_protein_identity, args.min_protein_coverage, args.min_protein_query_coverage]
    if any(not 0 <= value <= 1 for value in thresholds):
        parser.error("覆盖度、主导度、一致性和冲突重叠阈值必须在 0–1 之间")
    if args.min_isoform_molecules < 1 or args.utr_min_molecules < 1 or args.min_sj < 1 or args.threads < 1:
        parser.error("分子数、junction 支持阈值和线程数必须至少为 1")
    if args.tandem_max_gap < 0 or args.tandem_cds_ratio < 1 or args.tandem_min_members < 2 or args.tandem_min_multi_depth < 0:
        parser.error("串联间距和多重深度不能为负，CDS 长度比至少为 1，簇成员数至少为 2")
    return args


def main(argv=None):
    args = parse_args(argv)
    genome = Genome(args.genome)
    models = []
    for name, path in named_paths(args.set):
        loaded = parse_gff(path, name)
        print(f"[{name}] {len(loaded)} coding transcripts", flush=True)
        models.extend(loaded)
    if not models:
        raise RuntimeError("没有读到含 CDS 的转录本")
    sj = SJIndex()
    for path in args.sj:
        sj.load(path)
    load_bam_junctions(args.bam, genome, sj)
    chains = []
    for path in args.isoseq:
        if path.lower().endswith((".bam", ".cram")):
            chains.extend(load_iso_bam(path, genome))
        else:
            chains.extend(load_bed12(path))
    iso = IsoIndex(chains)
    protein_lengths = load_protein_lengths(args.proteins) if args.proteins else {}
    miniprot_gff = args.miniprot_gff
    if args.proteins and not miniprot_gff:
        miniprot_gff = run_miniprot(args.genome, args.proteins, args.miniprot, args.threads, args.outdir)
    protein_evidence = ProteinIndex(miniprot_gff or "", protein_lengths)
    print(f"evidence: {len(iso)} Iso-seq chains; {sum(len(v) for v in sj.data.values())} junctions; "
          f"{len(protein_evidence)} protein hits", flush=True)
    depth = DepthIndex(args.bam)
    try:
        weights = {"egapx": 3.0, "eviann": 3.0, "e": 3.0, "braker": 2.0, "helixer": 1.0, "annevo": 1.0}
        for name, value in named_paths(args.set_weight):
            weights[name.lower()] = float(value)
        evaluate_models(models, sj, iso, protein_evidence, depth, genome, weights, args.min_sj,
                        args.min_isoform_molecules, args.min_cds_coverage, args.min_protein_identity)
    finally:
        depth.close()
    loci = cluster_loci(models)
    trusted = {"egapx", "eviann", "e"} | {x.lower() for x in args.trusted_set}
    tandem_info = prepare_tandem_evidence(loci, sj, protein_evidence, args)
    winners, dropped = choose_representatives(loci, sj, iso, protein_evidence, args, trusted, tandem_info)
    print(f"selected {len(winners)} primary genes from {len(loci)} loci", flush=True)
    if args.utr_consensus:
        apply_utr_consensus(winners, iso, args)
    if args.max_isoforms > 1 and len(iso):
        attach_isoforms(winners, sj, iso, genome, args)
    summary = write_outputs(winners, dropped, models, genome, args)
    print(json.dumps(summary, ensure_ascii=False), flush=True)
    print(f"输出目录：{args.outdir}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
