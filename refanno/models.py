import re
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


def natural_key(value):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", value)]


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


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
    family_cluster: str = ""
    family_donors: list = field(default_factory=list)
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
