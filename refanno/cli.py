import argparse
import json
import sys
from datetime import datetime, timezone

from . import DESCRIPTION, __version__
from .audit import write_run_manifest
from .adjudication import choose_representatives, cluster_loci, evaluate_models, prepare_tandem_evidence
from .consensus import apply_utr_consensus, attach_isoforms
from .evidence import (
    DepthIndex,
    Genome,
    IsoIndex,
    ProteinIndex,
    SJIndex,
    load_bam_junctions,
    load_bed12,
    load_iso_bam,
    load_protein_lengths,
    parse_gff,
    run_miniprot,
)
from .output import write_outputs
from .splice import SpliceIndex, pick_training_models, run_minisplice


def named_paths(specs):
    out = []
    for spec in specs:
        if "=" not in spec:
            raise argparse.ArgumentTypeError(f"需要 NAME=PATH：{spec}")
        out.append(tuple(spec.split("=", 1)))
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="示例：\n  refanno --set annevo=a.gff3 --set helixer=h.gff3 "
                                            "--genome genome.fa --proteins related.faa --bam RNA.bam "
                                            "--outdir RefAnno_out")
    parser.add_argument("--version", action="version", version=f"RefAnno {__version__}")
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
    parser.add_argument("--minisplice", help="minisplice 可执行文件路径（conda 无此包，需从源码编译）")
    parser.add_argument("--splice-model", help="预训练 minisplice 模型 .kan；缺省时用输入注释自动训练")
    parser.add_argument("--splice-scores", help="已有 minisplice 基因组分数 .tsv(.gz)，跳过训练与预测")
    parser.add_argument("--min-splice-score", type=float, default=2.0, help="剪接供体/受体位点通过校验的最低分数")
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
    parser.add_argument("--tandem-min-protein-similarity", type=float, default=0.7,
                        help="重复家族共享证据的最低蛋白序列相似度")
    parser.add_argument("--trusted-set", action="append", default=[], help="允许无 RNA 证据单独保留的集合名")
    parser.add_argument("--set-weight", action="append", default=[], metavar="NAME=FLOAT")
    parser.add_argument("--utr-consensus", action="store_true", help="用 Iso-seq 局部主导链替换 UTR 外显子，CDS 不变")
    parser.add_argument("--utr-min-molecules", type=int, default=2)
    parser.add_argument("--utr-dominance", type=float, default=0.6)
    args = parser.parse_args(argv)
    if args.max_isoforms < 1 or args.max_isoforms > 3:
        parser.error("--max-isoforms 必须为 1–3")
    thresholds = [args.min_cds_coverage, args.utr_dominance, args.conflict_overlap,
                  args.min_protein_identity, args.min_protein_coverage, args.min_protein_query_coverage,
                  args.tandem_min_protein_similarity]
    if any(not 0 <= value <= 1 for value in thresholds):
        parser.error("覆盖度、主导度、一致性和冲突重叠阈值必须在 0–1 之间")
    if args.min_isoform_molecules < 1 or args.utr_min_molecules < 1 or args.min_sj < 1 or args.threads < 1:
        parser.error("分子数、junction 支持阈值和线程数必须至少为 1")
    if args.tandem_max_gap < 0 or args.tandem_cds_ratio < 1 or args.tandem_min_members < 2 or args.tandem_min_multi_depth < 0:
        parser.error("串联间距和多重深度不能为负，CDS 长度比至少为 1，簇成员数至少为 2")
    if args.splice_scores and (args.minisplice or args.splice_model):
        parser.error("--splice-scores 不能与 --minisplice/--splice-model 同时使用")
    if args.splice_model and not args.minisplice:
        parser.error("--splice-model 需要同时提供 --minisplice 可执行文件路径")
    return args


def run(argv=None):
    started = datetime.now(timezone.utc)
    args = parse_args(argv)
    genome = Genome(args.genome)
    models = []
    for name, path in named_paths(args.set):
        loaded = parse_gff(path, name)
        print(f"[{name}] {len(loaded)} coding transcripts", flush=True)
        models.extend(loaded)
    if not models:
        raise RuntimeError("没有读到含 CDS 的转录本")
    trusted = {"egapx", "eviann", "e"} | {x.lower() for x in args.trusted_set}
    loci = cluster_loci(models)
    splice_scores = args.splice_scores
    if args.minisplice and not splice_scores:
        train_models, train_source = pick_training_models(models, loci, trusted)
        print(f"minisplice: training on {len(train_models)} transcripts ({train_source})", flush=True)
        splice_scores = run_minisplice(args.minisplice, args.genome, args.threads, args.outdir,
                                       train_models, models, model_path=args.splice_model)
    splice_index = None
    if splice_scores:
        splice_index = SpliceIndex()
        loaded = splice_index.load(splice_scores, models)
        print(f"splice sites: {loaded} model junction sites scored", flush=True)
    sj = SJIndex()
    for path in args.sj:
        sj.load(path)
    depth = DepthIndex(args.bam, models)
    bam_stats = load_bam_junctions(args.bam, genome, sj, depth)
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
        miniprot_gff = run_miniprot(args.genome, args.proteins, args.miniprot, args.threads, args.outdir,
                                    spsc=splice_scores)
    protein_evidence = ProteinIndex(miniprot_gff or "", protein_lengths)
    hits = [hit for values in protein_evidence.by_chrom.values() for hit in values]
    identical = sum(1 for hit in hits if hit.identity >= 0.999)
    if len(hits) >= 1000 and identical >= 0.3 * len(hits):
        print(f"警告: 蛋白比对中 {100 * identical / len(hits):.0f}% 的命中一致性为 1.0；--proteins 可能来自同一基因组自身的预测，"
              "蛋白支持会退化为自我确认，建议改用近缘物种或 Swiss-Prot 等独立来源的蛋白", file=sys.stderr, flush=True)
    print(f"evidence: {len(iso)} Iso-seq chains; {sum(len(v) for v in sj.data.values())} junctions; "
          f"{len(protein_evidence)} protein hits", flush=True)
    weights = {"egapx": 3.0, "eviann": 3.0, "e": 3.0, "braker": 2.0, "helixer": 1.0, "annevo": 1.0}
    for name, value in named_paths(args.set_weight):
        weights[name.lower()] = float(value)
    evaluate_models(models, sj, iso, protein_evidence, depth, genome, weights, args.min_sj,
                    args.min_isoform_molecules, args.min_cds_coverage, args.min_protein_identity,
                    splice=splice_index, min_splice_score=args.min_splice_score)
    tandem_info = prepare_tandem_evidence(loci, sj, protein_evidence, genome, args)
    winners, dropped, relations = choose_representatives(loci, sj, iso, protein_evidence, args, trusted, tandem_info)
    print(f"selected {len(winners)} primary genes from {len(loci)} loci", flush=True)
    if args.utr_consensus:
        apply_utr_consensus(winners, iso, args)
    if args.max_isoforms > 1 and len(iso):
        attach_isoforms(winners, sj, iso, genome, args)
    summary = write_outputs(winners, dropped, relations, models, sj, iso, genome, args,
                            bam_stats=bam_stats, protein_evidence=protein_evidence)
    write_run_manifest(args, argv, started, summary, bam_stats)
    print(json.dumps(summary, ensure_ascii=False), flush=True)
    print(f"输出目录：{args.outdir}", flush=True)
    return summary


def main(argv=None):
    try:
        run(argv)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0
