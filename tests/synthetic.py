import random

import pysam

CODONS = ["GCT", "GAA", "CCT", "GGT", "CTG", "TCT", "CAG", "GAT", "ACT", "GTT", "CAT", "CGT"]
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))


def write_fasta(path, genome):
    with path.open("w") as out:
        out.write(">chr1\n")
        for i in range(0, len(genome), 80):
            out.write("".join(genome[i:i + 80]) + "\n")


def write_bam(path, reads, length):
    reads.sort(key=lambda read: read.reference_start)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(str(path))


def make_read(name, start, cigar):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * sum(size for op, size in cigar if op == 0)
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start - 1
    read.mapping_quality = 60
    read.cigartuples = cigar
    read.set_tag("NH", 1)
    return read


def build_single_exon(root, related_duplicates=True):
    random.seed(17)
    length = 100000
    genome = list("A" * length)

    def set_orf(start, size, internal=None):
        internal = internal or [random.choice(CODONS) for _ in range(size // 3 - 2)]
        seq = "ATG" + "".join(internal) + "TAA"
        assert len(seq) == size
        genome[start - 1:start - 1 + size] = list(seq)
        return seq

    set_orf(1001, 300)
    set_orf(1201, 399)
    protein_seq = set_orf(30001, 303)
    dup_seq = set_orf(50001, 600)
    if not related_duplicates:
        set_orf(50001, 600)
    set_orf(51201, 600, [dup_seq[i:i + 3] for i in range(3, len(dup_seq) - 3, 3)])
    set_orf(70001, 300)
    write_fasta(root / "genome.fa", genome)

    models = [
        ("setA", "short", 1001, 1300), ("setB", "long", 1201, 1599),
        ("setA", "protein", 30001, 30303), ("setA", "dup1", 50001, 50600),
        ("setB", "dup2", 51201, 51800), ("setA", "isolated", 70001, 70300),
    ]
    for set_name in ("setA", "setB"):
        with (root / f"{set_name}.gff3").open("w") as out:
            out.write("##gff-version 3\n")
            for source, name, start, end in models:
                if source != set_name:
                    continue
                out.write(f"chr1\tsynthetic\tgene\t{start}\t{end}\t.\t+\t.\tID={name}\n")
                out.write(f"chr1\tsynthetic\tmRNA\t{start}\t{end}\t.\t+\t.\tID={name}.t1;Parent={name}\n")
                out.write(f"chr1\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\tParent={name}.t1\n")
                out.write(f"chr1\tsynthetic\tCDS\t{start}\t{end}\t.\t+\t0\tParent={name}.t1\n")

    reads = []
    counter = 0

    def add_coverage(start, end, depth):
        nonlocal counter
        for pos in range(start, end - 50 + 2, 50):
            for _ in range(depth):
                counter += 1
                read = make_read(f"read{counter}", pos, [(0, 50)])
                reads.append(read)

    add_coverage(1001, 1300, 20)
    add_coverage(1201, 1599, 5)
    write_bam(root / "rna.bam", reads, length)

    with (root / "proteins.faa").open("w") as out:
        out.write(f">protein\n{translate(protein_seq)[:-1]}\n>duplicate\n{translate(dup_seq)[:-1]}\n")
    with (root / "evidence.gff3").open("w") as out:
        out.write("##gff-version 3\n")
        out.write("chr1\tminiprot\tmRNA\t30001\t30303\t.\t+\t.\tID=proteinHit;Target=protein 1 100;Identity=1.0\n")
        out.write("chr1\tminiprot\tCDS\t30001\t30303\t.\t+\t0\tParent=proteinHit\n")
        out.write("chr1\tminiprot\tmRNA\t51201\t51800\t.\t+\t.\tID=dupHit;Target=duplicate 1 199;Identity=1.0\n")
        out.write("chr1\tminiprot\tCDS\t51201\t51800\t.\t+\t0\tParent=dupHit\n")
    with (root / "long_protein.gff3").open("w") as out:
        out.write("##gff-version 3\n")
        out.write("chr1\tminiprot\tmRNA\t1201\t1599\t.\t+\t.\tID=longHit;Target=longProtein 1 100;Identity=1.0\n")
        out.write("chr1\tminiprot\tCDS\t1201\t1599\t.\t+\t0\tParent=longHit\n")


def build_multiexon(root):
    random.seed(23)
    genome = list("A" * 6000)
    cds = "ATG" + "".join(random.choice(CODONS[:8]) for _ in range(98)) + "TAA"
    for start in (1001, 2001):
        genome[start - 1:start + 149] = list(cds[:150] if start == 1001 else cds[150:])
    genome[3000:3150] = list(cds[150:])
    genome[1150:1152] = list("GT")
    genome[1998:2000] = list("AG")
    genome[2998:3000] = list("AG")
    write_fasta(root / "genome.fa", genome)

    for set_name, second in (("setA", 2001), ("setB", 3001)):
        with (root / f"{set_name}.gff3").open("w") as out:
            out.write("##gff-version 3\n")
            out.write(f"chr1\tsynthetic\tgene\t1001\t{second + 149}\t.\t+\t.\tID=gene{set_name}\n")
            out.write(f"chr1\tsynthetic\tmRNA\t1001\t{second + 149}\t.\t+\t.\tID=tx{set_name};Parent=gene{set_name}\n")
            for start in (1001, second):
                out.write(f"chr1\tsynthetic\texon\t{start}\t{start + 149}\t.\t+\t.\tParent=tx{set_name}\n")
                out.write(f"chr1\tsynthetic\tCDS\t{start}\t{start + 149}\t.\t+\t0\tParent=tx{set_name}\n")

    reads = []
    counter = 0

    def add_read(start, cigar, name):
        nonlocal counter
        counter += 1
        read = make_read(f"{name}{counter}", start, cigar)
        read.set_tag("XS", "+")
        reads.append(read)

    for pos in list(range(1001, 1101, 50)) + list(range(2001, 2101, 50)):
        for _ in range(15):
            add_read(pos, [(0, 50)], "depth")
    for _ in range(6):
        add_read(1101, [(0, 50), (3, 850), (0, 50)], "splice")
    write_bam(root / "rna.bam", reads, len(genome))

    with (root / "iso.bed12").open("w") as out:
        for name, second in (("A1", 2001), ("A2", 2001), ("B1", 3001), ("B2", 3001)):
            offset = second - 901
            out.write(f"chr1\t900\t{second + 199}\t{name}\t60\t+\t900\t{second + 199}\t0\t2\t250,200,\t0,{offset},\n")
