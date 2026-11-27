import csv
import html
import json
import os
from collections import Counter

from . import OUTPUT_SCHEMA_VERSION, __version__
from .adjudication import chain_groups, cluster_loci, detect_duplicate_clusters
from .models import attr_text, natural_key, overlap, translate


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
    if protein_handle is not None:
        protein = translate(genome.cds_sequence(model))
        protein_handle.write(f">{tx_id} source_model={model.source_model} role={role}\n")
        for i in range(0, len(protein), 60):
            protein_handle.write(protein[i:i + 60] + "\n")


def write_decision_evidence(outdir, models, records, relations, sj, iso, genome, args):
    event_fields = ["source_model", "gene_id", "transcript_id", "outcome", "evidence_class",
                    "reason", "selection", "rule", "related_source_model", "related_gene_id",
                    "cds_overlap_source", "cds_overlap_related"]
    junction_fields = ["gene_id", "transcript_id", "outcome", "source_model", "chrom",
                       "intron_start", "intron_end", "strand", "unique_reads", "multi_reads",
                       "isoseq_molecules", "motif", "supported", "source_bams"]
    model_by_key = {model.key: model for model in models}
    with open(os.path.join(outdir, "decision_events.tsv"), "w", newline="") as events, \
            open(os.path.join(outdir, "junction_evidence.tsv"), "w", newline="") as junctions, \
            open(os.path.join(outdir, "supporting_junctions.bed"), "w") as bed, \
            open(os.path.join(outdir, "family_rescues.tsv"), "w", newline="") as rescues:
        event_writer = csv.DictWriter(events, fieldnames=event_fields, delimiter="\t")
        junction_writer = csv.DictWriter(junctions, fieldnames=junction_fields, delimiter="\t")
        rescue_fields = ["rescued_gene_id", "rescued_source_model", "cluster_id", "donor_gene_id",
                         "donor_source_model", "evidence", "protein_target", "value"]
        rescue_writer = csv.DictWriter(rescues, fieldnames=rescue_fields, delimiter="\t")
        event_writer.writeheader()
        junction_writer.writeheader()
        rescue_writer.writeheader()
        for model in sorted(models, key=lambda item: item.key):
            record = records[model.key]
            relation = relations.get(model.key, {})
            related_key = record.get("related_model", relation.get("related_model", ""))
            related = model_by_key.get(related_key)
            related_record = records.get(related_key, {})
            overlap_values = relation.get("overlap", ())
            event_writer.writerow({
                "source_model": model.source_model, "gene_id": record["gene_id"],
                "transcript_id": record["transcript_id"], "outcome": record["outcome"],
                "evidence_class": record["evidence_class"], "reason": record["reason"],
                "selection": record["selection"],
                "rule": relation.get("rule", "isoseq_alternative" if related_key else "representative"),
                "related_source_model": related.source_model if related else "",
                "related_gene_id": related_record.get("gene_id", ""),
                "cds_overlap_source": f"{overlap_values[0]:.3f}" if overlap_values else "",
                "cds_overlap_related": f"{overlap_values[1]:.3f}" if overlap_values else "",
            })
            for donor in model.family_donors:
                donor_model = model_by_key[donor["donor"]]
                rescue_writer.writerow({
                    "rescued_gene_id": record["gene_id"], "rescued_source_model": model.source_model,
                    "cluster_id": model.family_cluster,
                    "donor_gene_id": records.get(donor["donor"], {}).get("gene_id", ""),
                    "donor_source_model": donor_model.source_model, "evidence": donor["evidence"],
                    "protein_target": donor["protein_target"], "value": donor["value"],
                })
            for intron in model.cds_introns:
                unique, multi = sj.get(model, intron)
                molecules = iso.junction_support(model, intron)
                supported = unique >= args.min_sj or molecules >= args.min_isoform_molecules
                motif = genome.splice_strand(model.chrom, [intron]) or "."
                junction_writer.writerow({
                    "gene_id": record["gene_id"], "transcript_id": record["transcript_id"],
                    "outcome": record["outcome"], "source_model": model.source_model,
                    "chrom": model.chrom, "intron_start": intron[0], "intron_end": intron[1],
                    "strand": model.strand, "unique_reads": unique, "multi_reads": multi,
                    "isoseq_molecules": molecules, "motif": motif, "supported": int(supported),
                    "source_bams": "+".join(os.path.basename(path) for path in args.bam),
                })
                if supported and record["outcome"] != "dropped":
                    bed.write("\t".join([
                        model.chrom, str(intron[0] - 1), str(intron[1]), record["transcript_id"],
                        str(min(1000, unique + multi + molecules)), model.strand,
                    ]) + "\n")


def format_blocks(blocks):
    return ",".join(f"{start}-{end}" for start, end in blocks)


def format_family_donors(donors):
    return ";".join(f"{donor['donor']}:{donor['evidence']}" for donor in donors)


def _h(value):
    return html.escape(str(value))


def _bar(value, total, color="#4c78a8"):
    pct = 100.0 * value / total if total else 0.0
    return (f'<td class="num">{value:,}</td>'
            f'<td class="barcell"><div class="track"><div class="fill" '
            f'style="width:{pct:.1f}%;background:{color}"></div></div><span>{pct:.1f}%</span></td>')


def _dist_table(counts, total=None, order=None, color="#4c78a8", label="Category"):
    total = sum(counts.values()) if total is None else total
    keys = [k for k in order if counts.get(k)] if order else sorted(counts, key=lambda k: (-counts[k], str(k)))
    keys += [k for k in sorted(counts) if k not in keys]
    rows = "".join(f"<tr><th>{_h(key)}</th>{_bar(counts[key], total, color)}</tr>" for key in keys if counts.get(key))
    return f'<table><tr><th>{_h(label)}</th><th>Count</th><th>Share</th></tr>{rows}</table>'


def _table(headers, rows, limit=None, note=""):
    shown = rows[:limit] if limit else rows
    body = "".join("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>" for row in shown)
    head = "".join(f"<th>{_h(h)}</th>" for h in headers)
    suffix = f'<p class="dim">仅列出前 {len(shown):,} 行（共 {len(rows):,} 行），全量见 {_h(note)}。</p>' if limit and len(rows) > limit else ""
    return f"<table><tr>{head}</tr>{body}</table>{suffix}"


EVIDENCE_ORDER = ["isoseq_exact", "short_read_complete", "short_read_partial",
                  "homology_only", "family_indirect", "sequence_only", "untested", "annotation_only"]
STATUS_ORDER = ["exact_isoseq", "supported", "sequence_supported", "partial", "unsupported", "untested"]


def _winner_junction_stats(winners, records, sj, iso, genome, args):
    totals = {"junctions": 0, "sr_supported": 0, "iso_supported": 0,
              "unsupported": 0, "noncanonical": 0}
    weak = []
    for winner in winners:
        model = winner["model"]
        for intron in model.cds_introns:
            unique, multi = sj.get(model, intron)
            molecules = iso.junction_support(model, intron)
            supported = unique >= args.min_sj or molecules >= args.min_isoform_molecules
            motif = genome.splice_strand(model.chrom, [intron])
            totals["junctions"] += 1
            totals["sr_supported"] += unique >= args.min_sj
            totals["iso_supported"] += molecules >= args.min_isoform_molecules
            totals["unsupported"] += not supported
            totals["noncanonical"] += motif not in {"+", "-"}
            if not supported and (unique or multi or molecules):
                weak.append((unique + multi, records[model.key]["gene_id"], model.source_model,
                             model.chrom, intron[0], intron[1], model.strand,
                             unique, multi, molecules, motif or "."))
    weak.sort(reverse=True)
    return totals, weak


def _bin(value, edges):
    for label, low, high in edges:
        if low <= value < high:
            return label
    return edges[-1][0]


def write_review_qc(outdir, winners, models, records, relations, genome, args,
                    sj=None, iso=None, bam_stats=None, protein_evidence=None, duplicate_info=None):
    review_winners = [winner for winner in winners
                      if winner["model"].review
                      or winner["model"].structure_status in {"partial", "unsupported", "untested"}
                      or winner["reason"] == "tandem_family_support"]
    changes = [winner for winner in winners if winner.get("change")]
    with open(os.path.join(outdir, "changes.tsv"), "w", newline="") as out:
        fields = ["gene_id", "transcript_id", "source_model", "change_type", "old_exons",
                  "new_exons", "cds_changed", "supporting_molecules", "molecule_count", "dominance"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for winner in changes:
            model = winner["model"]
            change = winner["change"]
            writer.writerow({
                "gene_id": records[model.key]["gene_id"], "transcript_id": records[model.key]["transcript_id"],
                "source_model": model.source_model, "change_type": change["type"],
                "old_exons": format_blocks(change["old_exons"]), "new_exons": format_blocks(change["new_exons"]),
                "cds_changed": 0, "supporting_molecules": "|".join(change["molecules"]),
                "molecule_count": change["support"], "dominance": f"{change['dominance']:.3f}",
            })

    with open(os.path.join(outdir, "review.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        for winner in review_winners:
            model = winner["model"]
            gene_id = records[model.key]["gene_id"]
            attributes = [("ID", gene_id), ("source_model", model.source_model),
                          ("reason", winner["reason"]), ("evidence_class", winner["evidence_class"]),
                          ("structure_status", model.structure_status)]
            if model.review:
                attributes.append(("flag_review", ",".join(sorted(set(model.review)))))
            out.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], attributes) + "\n")
            write_transcript(out, None, model, gene_id, 1, genome, "primary", model.iso_support)

    with open(os.path.join(outdir, "conflicting_models.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        for model in sorted(models, key=lambda item: item.key):
            record = records[model.key]
            if record["reason"] != "locus_conflict_lower_rank":
                continue
            related_key = relations.get(model.key, {}).get("related_model", "")
            related_gene = records.get(related_key, {}).get("gene_id", "")
            attributes = [("ID", record["gene_id"]), ("source_model", model.source_model),
                          ("reason", record["reason"]), ("related_gene_id", related_gene)]
            out.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], attributes) + "\n")
            write_transcript(out, None, model, record["gene_id"], 1, genome, "dropped")

    fates = Counter()
    for model in models:
        fates[(model.set_name, records[model.key]["outcome"])] += 1
    with open(os.path.join(outdir, "input_fates.tsv"), "w", newline="") as out:
        fields = ["set_name", "input_transcripts", "complete_cds", "structure_supported",
                  "structure_partial", "primary", "alternative", "dropped"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for set_name in sorted({model.set_name for model in models}):
            members = [model for model in models if model.set_name == set_name]
            writer.writerow({
                "set_name": set_name, "input_transcripts": len(members),
                "complete_cds": sum(model.complete_cds for model in members),
                "structure_supported": sum(model.structure_status in {"exact_isoseq", "supported", "sequence_supported"} for model in members),
                "structure_partial": sum(model.structure_status == "partial" for model in members),
                "primary": fates[(set_name, "primary")], "alternative": fates[(set_name, "alternative")],
                "dropped": fates[(set_name, "dropped")],
            })

    primary_models = [winner["model"] for winner in winners]
    total_genes = len(winners)
    evidence_counts = Counter(winner["evidence_class"] for winner in winners)
    status_counts = Counter(winner["model"].structure_status for winner in winners)
    reason_counts = Counter(winner["reason"] for winner in winners)
    selection_counts = Counter(winner["selection"] for winner in winners)
    tier_counts = Counter(output_tier(winner) for winner in winners)
    dropped_reasons = Counter(records[model.key]["reason"] for model in models
                              if records[model.key]["outcome"] == "dropped")

    sections = []

    set_counts = Counter(model.set_name for model in models)
    overview_rows = [
        ("RefAnno version", __version__),
        ("Output schema", OUTPUT_SCHEMA_VERSION),
        ("Input sets", "; ".join(f"{name} ({count:,} transcripts)" for name, count in sorted(set_counts.items()))),
        ("Input transcripts", len(models)),
        ("Primary genes", total_genes),
        ("Final (supported) genes", tier_counts.get("final", 0)),
        ("Weak (low-support) genes", tier_counts.get("weak", 0)),
        ("Transcripts", sum(1 + len(winner["alternatives"]) for winner in winners)),
        ("Structural changes (UTR consensus)", len(changes)),
        ("Genes requiring review", len(review_winners)),
        ("Key thresholds", f"min_sj={args.min_sj}; min_isoform_molecules={args.min_isoform_molecules}; "
                           f"max_isoforms={args.max_isoforms}; min_cds_coverage={args.min_cds_coverage}"),
    ]
    if bam_stats:
        overview_rows.append(("RNA-seq BAMs", "; ".join(os.path.basename(s["path"]) for s in bam_stats)))
    if iso is not None and len(iso):
        overview_rows.append(("Iso-seq chains", f"{len(iso):,}"))
    if protein_evidence is not None and len(protein_evidence):
        overview_rows.append(("Protein hits", f"{len(protein_evidence):,}"))
    if getattr(args, "splice_scores", None) or getattr(args, "minisplice", None):
        overview_rows.append(("Splice validation", "minisplice"))
    overview = "".join(
        f"<tr><th>{_h(k)}</th><td class='num'>{v:,}</td></tr>" if isinstance(v, int)
        else f"<tr><th>{_h(k)}</th><td>{_h(v)}</td></tr>" for k, v in overview_rows)
    sections.append(f"<h2>1. 运行概览</h2><table>{overview}</table>")

    cross = {}
    for winner in winners:
        model = winner["model"]
        bucket = cross.setdefault(winner["evidence_class"], [0, 0])
        bucket[0] += 1
        bucket[1] += len(model.cds) == 1
    cross_rows = []
    for key in [k for k in EVIDENCE_ORDER if k in cross] + sorted(set(cross) - set(EVIDENCE_ORDER)):
        total, single = cross[key]
        pct = f"{100.0 * single / total:.1f}%" if total else "0"
        cross_rows.append((f"<b>{_h(key)}</b>", f"{total:,}", f"{single:,}", pct))
    sections.append(
        "<h2>2. 证据链</h2>"
        "<h3>2.1 证据等级分布</h3>" + _dist_table(evidence_counts, order=EVIDENCE_ORDER, label="Evidence class") +
        "<h3>2.2 保留原因</h3>" + _dist_table(reason_counts, label="Reason") +
        "<h3>2.3 主代表选择依据</h3>" + _dist_table(selection_counts, label="Selection") +
        "<h3>2.4 结构状态</h3>" + _dist_table(status_counts, order=STATUS_ORDER, label="Structure status") +
        "<h3>2.5 证据等级 × 单外显子基因</h3>"
        '<p class="dim">单外显子比例高的证据等级（尤其 homology_only）通常意味着假阳性风险。</p>'
        + _table(["Evidence class", "Genes", "Single-exon", "Single-exon %"], cross_rows))

    tier_evidence = Counter((output_tier(winner), winner["evidence_class"]) for winner in winners)
    tier_rows = []
    for key in [k for k in EVIDENCE_ORDER if evidence_counts.get(k)] + sorted(set(evidence_counts) - set(EVIDENCE_ORDER)):
        final_n, weak_n = tier_evidence.get(("final", key), 0), tier_evidence.get(("weak", key), 0)
        tier_rows.append((f"<b>{_h(key)}</b>", f"{final_n:,}", f"{weak_n:,}"))
    sections.append(
        "<h3>2.6 final/weak 分层 × 证据等级</h3>"
        '<p class="dim">final 层：直接证据或多集合共识；weak 层：family_indirect、annotation_only、'
        'untested 且非共识，写入 weak.gff3。</p>'
        + _table(["Evidence class", "Final tier", "Weak tier"], tier_rows))

    splice_active = bool(getattr(args, "splice_scores", None) or getattr(args, "minisplice", None))
    if splice_active:
        multi = [m for m in primary_models if m.junctions_total]
        full = sum(m.splice_supported == m.junctions_total for m in multi)
        partial = sum(0 < m.splice_supported < m.junctions_total for m in multi)
        none = sum(m.splice_supported == 0 for m in multi)
        scored = [m.splice_min_score for m in multi if m.splice_min_score is not None]
        score_edges = [("<0", -1e9, 0), ("0–2", 0, 2), ("2–6", 2, 6), ("6–10", 6, 10), ("≥10", 10, 1e9)]
        score_bins = Counter(_bin(v, score_edges) for v in scored)
        splice_rows = [
            ("Multi-exon primary genes", len(multi)),
            (f"All junctions pass (donor & acceptor ≥ {args.min_splice_score:g})", full),
            ("Partially pass", partial),
            ("No junction passes", none),
            ("Flagged weak_splice_sites", sum("weak_splice_sites" in m.review for m in primary_models)),
        ]
        sections.append(
            "<h3>2.7 剪接位点序列校验（minisplice）</h3><table>"
            + "".join(f"<tr><th>{_h(k)}</th><td class='num'>{v:,}</td></tr>" for k, v in splice_rows)
            + "</table><h4>最弱位点分数分布</h4>"
            + _dist_table(score_bins, order=[e[0] for e in score_edges], label="Min site score")
            + '<p class="dim">序列校验是注释自洽性检查：ab initio 预测器本就会产生规范剪接位点，'
              '因此通过校验不代表独立证据，未通过则值得人工查看。</p>')

    fate_rows = []
    for set_name in sorted(set_counts):
        members = [model for model in models if model.set_name == set_name]
        kept = fates[(set_name, "primary")] + fates[(set_name, "alternative")]
        tier_final = sum(1 for winner in winners
                         if winner["model"].set_name == set_name and output_tier(winner) == "final")
        tier_weak = fates[(set_name, "primary")] - tier_final
        fate_rows.append(
            f"<tr><th>{_h(set_name)}</th><td class='num'>{len(members):,}</td>"
            f"<td class='num'>{sum(m.complete_cds for m in members):,}</td>"
            f"<td class='num'>{sum(m.structure_status in {'exact_isoseq', 'supported', 'sequence_supported'} for m in members):,}</td>"
            f"<td class='num'>{sum(m.structure_status == 'partial' for m in members):,}</td>"
            f"<td class='num'>{fates[(set_name, 'primary')]:,}</td>"
            f"<td class='num'>{tier_final:,}</td>"
            f"<td class='num'>{tier_weak:,}</td>"
            f"<td class='num'>{fates[(set_name, 'alternative')]:,}</td>"
            f"<td class='num'>{fates[(set_name, 'dropped')]:,}</td>"
            + _bar(kept, len(members), "#59a14f") + "</tr>")
    pair_counts = Counter()
    set_chains = Counter()
    for locus in cluster_loci(models):
        for group in chain_groups(locus):
            group_sets = sorted({model.set_name for model in group})
            for set_name in group_sets:
                set_chains[set_name] += 1
            for i, first in enumerate(group_sets):
                for second in group_sets[i + 1:]:
                    pair_counts[(first, second)] += 1
    pair_rows = []
    for (first, second), count in pair_counts.most_common(20):
        base = min(set_chains[first], set_chains[second]) or 1
        pair_rows.append((f"<b>{_h(first)}</b>", f"<b>{_h(second)}</b>", f"{count:,}",
                          f"{100.0 * count / base:.1f}%"))
    sections.append(
        "<h2>3. 输入集合命运</h2><table><tr><th>Set</th><th>Transcripts</th><th>Complete CDS</th>"
        "<th>Structure supported</th><th>Partial</th><th>Primary</th><th>Final tier</th><th>Weak tier</th>"
        "<th>Alternative</th><th>Dropped</th><th>Retention</th><th></th></tr>" + "".join(fate_rows) + "</table>"
        "<h3>3.1 淘汰原因分布</h3>" + _dist_table(dropped_reasons, label="Drop reason")
        + "<h3>3.2 集合间一致链</h3>"
        + '<p class="dim">同一座位内 CDS 剪接链完全一致（单外显子按 CDS 重叠 ≥90%）的集合对；'
          '百分比以两集合中较少的链数为基数。这些一致链是 2set/3set_chain_consensus 和 minisplice 自训练训练集的来源。</p>'
        + _table(["Set A", "Set B", "Shared chains", "Share"], pair_rows))

    exon_edges = [("1", 1, 2), ("2", 2, 3), ("3–5", 3, 6), ("6–10", 6, 11), ("11–20", 11, 21), (">20", 21, 1 << 62)]
    cds_edges = [("<300", 0, 300), ("300–600", 300, 600), ("600–1200", 600, 1200),
                 ("1200–2400", 1200, 2400), ("2400–6000", 2400, 6000), (">6000", 6000, 1 << 62)]
    exon_bins = Counter(_bin(len(m.exons), exon_edges) for m in primary_models)
    cds_bins = Counter(_bin(m.cds_len, cds_edges) for m in primary_models)
    isoform_bins = Counter(str(1 + len(winner["alternatives"])) for winner in winners)
    consensus3 = reason_counts.get("3set_chain_consensus", 0)
    consensus2 = reason_counts.get("2set_chain_consensus", 0)
    struct_rows = [
        ("Complete CDS", sum(m.complete_cds for m in primary_models)),
        ("Internal stop codon", sum("internal_stop" in m.review for m in primary_models)),
        ("Single-exon genes", sum(len(m.cds) == 1 for m in primary_models)),
        ("Mean exons per gene", f"{sum(len(m.exons) for m in primary_models) / total_genes:.2f}" if total_genes else 0),
        ("3-set chain consensus", consensus3),
        ("2-set chain consensus", consensus2),
        ("Multi-isoform loci", sum(bool(winner["alternatives"]) for winner in winners)),
    ]
    struct_overview = "".join(f"<tr><th>{_h(k)}</th><td class='num'>{_h(v)}</td></tr>" for k, v in struct_rows)
    sections.append(
        "<h2>4. 最终基因集结构</h2><table>" + struct_overview + "</table>"
        "<h3>4.1 外显子数分布</h3>" + _dist_table(exon_bins, order=[e[0] for e in exon_edges], label="Exons") +
        "<h3>4.2 CDS 长度分布 (bp)</h3>" + _dist_table(cds_bins, order=[e[0] for e in cds_edges], label="CDS length") +
        "<h3>4.3 每基因转录本数分布</h3>" + _dist_table(isoform_bins, order=["1", "2", "3"], label="Isoforms"))

    rna_sections = []
    if bam_stats:
        bam_rows = []
        for stat in bam_stats:
            pct = 100.0 * stat["spliced_reads"] / stat["mapq_pass"] if stat["mapq_pass"] else 0.0
            warn = " class='warn'" if stat["reads"] >= 10000 and pct < 1 else ""
            bam_rows.append((f"<b>{_h(os.path.basename(stat['path']))}</b>", f"{stat['reads']:,}",
                             f"{stat['mapq_pass']:,}", f"{stat['spliced_reads']:,}",
                             f"<span{warn}>{pct:.1f}%</span>"))
        rna_sections.append(
            "<h3>5.1 逐 BAM 比对统计</h3>"
            '<p class="dim">剪接读段占比 &lt;1% 的 BAM 会被标黄：它可能不是链特异性 RNA-seq（例如 Hi-C），对 junction 证据贡献有限。</p>'
            + _table(["BAM", "Reads", "MAPQ pass", "Spliced reads", "Spliced %"], bam_rows))
    if sj is not None and iso is not None:
        totals, weak = _winner_junction_stats(winners, records, sj, iso, genome, args)
        junction_rows = [
            ("CDS junctions (primary transcripts)", totals["junctions"]),
            (f"Supported by short reads (unique ≥ {args.min_sj})", totals["sr_supported"]),
            (f"Supported by Iso-seq (molecules ≥ {args.min_isoform_molecules})", totals["iso_supported"]),
            ("Unsupported", totals["unsupported"]),
            ("Non-canonical splice motif", totals["noncanonical"]),
        ]
        junction_overview = "".join(f"<tr><th>{_h(k)}</th><td class='num'>{v:,}</td></tr>" for k, v in junction_rows)
        weak_rows = [(f"<b>{_h(gid)}</b>", _h(src), _h(chrom), f"{start:,}–{end:,}", _h(strand),
                      f"{u:,}", f"{m:,}", f"{mol:,}", _h(motif))
                     for _, gid, src, chrom, start, end, strand, u, m, mol, motif in weak]
        rna_sections.append(
            "<h3>5.2 主代表 junction 支持</h3><table>" + junction_overview + "</table>"
            "<h3>5.3 弱支持 junction（Top 30）</h3>"
            + _table(["Gene", "Source model", "Chrom", "Intron", "Strand", "Unique", "Multi", "Iso-seq", "Motif"],
                     weak_rows, limit=30, note="junction_evidence.tsv"))
    if args.bam:
        depth_edges = [("0", -1e-9, 1e-9), ("(0, 1)", 1e-9, 1), ("1–10", 1, 10), ("≥10", 10, float("inf"))]
        cov_edges = [("≥0.9", 0.9, 2.0), ("0.5–0.9", 0.5, 0.9), ("<0.5", -1e-9, 0.5)]
        depth_bins = Counter(_bin(m.depth, depth_edges) for m in primary_models)
        cov_bins = Counter(_bin(m.coverage, cov_edges) for m in primary_models)
        dc = Counter((_bin(m.depth, depth_edges), _bin(m.coverage, cov_edges)) for m in primary_models)
        dc_rows = []
        for dlabel, _, _ in depth_edges:
            row = [f"<b>{_h(dlabel)}</b>"] + [f"{dc.get((dlabel, clabel), 0):,}" for clabel, _, _ in cov_edges]
            dc_rows.append(tuple(row))
        rna_sections.append(
            "<h3>5.4 主代表 CDS 深度与覆盖度</h3><div class='cols'><div><h4>Depth</h4>"
            + _dist_table(depth_bins, order=[e[0] for e in depth_edges], label="Depth")
            + "</div><div><h4>Coverage</h4>"
            + _dist_table(cov_bins, order=[e[0] for e in cov_edges], label="Coverage") + "</div></div>"
            + "<h3>5.5 深度 × 覆盖度</h3>"
            + '<p class="dim">高覆盖但零深度的基因通常是比对多重性或注释冗余造成的假象，值得留意。</p>'
            + _table(["Depth \\ Coverage"] + [e[0] for e in cov_edges], dc_rows))
    if rna_sections:
        sections.append("<h2>5. RNA 证据</h2>" + "".join(rna_sections))

    if iso is not None and len(iso):
        utr_applied = sum(bool(winner["model"].utr_note) for winner in winners)
        utr_changed = sum(bool(winner.get("utr_changed")) for winner in winners)
        utr_flags = Counter(flag for m in primary_models for flag in m.review if flag.startswith("utr_"))
        deltas = []
        for winner in winners:
            change = winner.get("change")
            if change and change.get("type") == "utr_consensus":
                old_len = sum(e - s + 1 for s, e in change["old_exons"])
                new_len = sum(e - s + 1 for s, e in change["new_exons"])
                deltas.append(new_len - old_len)
        delta_edges = [("≤-500", -1 << 60, -500), ("-500–-100", -500, -100), ("-100–100", -100, 100),
                       ("100–500", 100, 500), ("≥500", 500, 1 << 60)]
        delta_bins = Counter(_bin(d, delta_edges) for d in deltas)
        alt_total = sum(len(winner["alternatives"]) for winner in winners)
        iso_rows = [
            ("Iso-seq chains", len(iso)),
            ("Genes with exact Iso-seq support", evidence_counts.get("isoseq_exact", 0)),
            ("Alternative transcripts retained", alt_total),
            ("Multi-isoform loci", sum(bool(winner["alternatives"]) for winner in winners)),
            ("UTR consensus applied", utr_applied),
            ("UTR consensus changed exons", utr_changed),
            ("UTR consensus ambiguous", utr_flags.get("utr_consensus_ambiguous", 0)),
            ("UTR consensus incomplete", utr_flags.get("utr_consensus_incomplete", 0)),
        ]
        iso_html = "<table>" + "".join(f"<tr><th>{_h(k)}</th><td class='num'>{_h(v)}</td></tr>"
                                       for k, v in iso_rows) + "</table>"
        if deltas:
            mean_delta = sum(deltas) / len(deltas)
            iso_html += (f"<h3>6.1 UTR 修订长度变化（n={len(deltas):,}，均值 {mean_delta:+.0f} bp）</h3>"
                         + _dist_table(delta_bins, order=[e[0] for e in delta_edges], label="Δ exon length (bp)"))
        sections.append("<h2>6. Iso-seq 证据</h2>" + iso_html)

    if protein_evidence is not None and len(protein_evidence):
        hits = [hit for values in protein_evidence.by_chrom.values() for hit in values]
        id_edges = [("<0.5", -1e-9, 0.5), ("0.5–0.8", 0.5, 0.8), ("0.8–0.95", 0.8, 0.95),
                    ("0.95–<1.0", 0.95, 0.999), ("=1.0", 0.999, 2.0)]
        id_bins = Counter(_bin(hit.identity, id_edges) for hit in hits)
        identical = id_bins.get("=1.0", 0)
        warning = ""
        if len(hits) >= 1000 and identical >= 0.3 * len(hits):
            warning = (f"<p class='warnbox'>警告：{100.0 * identical / len(hits):.0f}% 的蛋白命中一致性为 1.0，"
                       "--proteins 很可能来自同一基因组自身的预测，蛋白支持退化为自我确认；"
                       "建议改用近缘物种或 Swiss-Prot 等独立来源。</p>")
        protein_only = [m for m in primary_models if records[m.key]["reason"] == "unique_with_protein"]
        risky = sum(len(m.cds) == 1 and m.depth == 0 for m in protein_only)
        protein_rows = [
            ("Protein hits", len(hits)),
            ("Genes with protein support", sum(m.protein_coverage > 0 for m in primary_models)),
            ("Protein-only retained genes", len(protein_only)),
            ("Protein-only + single-exon + zero depth (high risk)", risky),
        ]
        sections.append(
            "<h2>7. 蛋白同源证据</h2>" + warning
            + "<table>" + "".join(f"<tr><th>{_h(k)}</th><td class='num'>{v:,}</td></tr>"
                                  for k, v in protein_rows) + "</table>"
            + "<h3>7.1 命中一致性分布</h3>"
            + _dist_table(id_bins, order=[e[0] for e in id_edges], label="Identity"))

    if duplicate_info is None:
        duplicate_info = detect_duplicate_clusters(primary_models, args)
    tandem_clusters = len({d["cluster"] for d in duplicate_info.values() if d["type"] == "tandem"})
    tandem_genes = sum(d["type"] == "tandem" for d in duplicate_info.values())
    dup_clusters = len({d["cluster"] for d in duplicate_info.values() if d["type"] == "duplicate"})
    dup_genes = sum(d["type"] == "duplicate" for d in duplicate_info.values())
    rescues = []
    for winner in winners:
        model = winner["model"]
        for donor in model.family_donors:
            donor_record = records.get(donor["donor"], {})
            rescues.append((f"<b>{_h(records[model.key]['gene_id'])}</b>", _h(model.source_model),
                            _h(model.family_cluster), _h(donor_record.get("gene_id", "")),
                            _h(donor["evidence"]), _h(donor["protein_target"]),
                            f"{donor['value']:.3f}" if isinstance(donor["value"], float) else _h(donor["value"])))
    dup_rows = [("Tandem clusters", tandem_clusters), ("Tandem genes", tandem_genes),
                ("Duplicate clusters", dup_clusters), ("Duplicate genes", dup_genes),
                ("Family-rescued genes", sum(winner["reason"] == "tandem_family_support" for winner in winners))]
    dup_html = "<table>" + "".join(f"<tr><th>{_h(k)}</th><td class='num'>{v:,}</td></tr>" for k, v in dup_rows) + "</table>"
    if rescues:
        dup_html += ("<h3>8.1 家族救援明细</h3>"
                     + _table(["Rescued gene", "Source model", "Cluster", "Donor gene", "Evidence", "Protein", "Value"],
                              rescues, limit=100, note="family_rescues.tsv"))
    sections.append("<h2>8. 重复与基因家族</h2>" + dup_html)

    flag_counts = Counter(flag for m in primary_models for flag in set(m.review))
    flag_html = _dist_table(flag_counts, label="Review flag") if flag_counts else "<p class='dim'>无 review 标记。</p>"
    risk_rows = []
    for winner in winners:
        model = winner["model"]
        score = risk_score(model, records[model.key]["reason"], winner["evidence_class"])
        if score:
            risk_rows.append((score, model, winner))
    risk_rows.sort(key=lambda item: (-item[0], natural_key(item[1].chrom), item[1].cds[0][0], item[1].key))
    risk_table = []
    for score, model, winner in risk_rows[:50]:
        risk_table.append((
            str(score), f"<b>{_h(records[model.key]['gene_id'])}</b>", _h(model.source_model),
            _h(winner["evidence_class"]), _h(records[model.key]["reason"]),
            _h(",".join(sorted(set(model.review)))), str(len(model.exons)),
            f"{model.depth:.1f}", f"{model.coverage:.2f}",
            f"{model.junctions_supported}/{model.junctions_total}",
            f"{model.protein_identity:.2f}", f"{model.protein_coverage:.2f}"))
    sections.append(
        "<h2>9. 复查优先级</h2><h3>9.1 Review 标记计数</h3>" + flag_html
        + "<h3>9.2 高风险基因（Top 50）</h3>"
        + '<p class="dim">风险分：internal_stop +4；仅蛋白+单外显子+零深度 +3；family_indirect +2；'
          '存在不支持 junction +2；UTR 歧义/不完整 +1；weak_splice_sites +1；结构 partial/untested +1。全量见 decision_report.tsv 的 review 列。</p>'
        + _table(["Score", "Gene", "Source model", "Evidence class", "Reason", "Flags", "Exons",
                  "Depth", "Coverage", "Junctions", "Protein id", "Protein cov"], risk_table, limit=50))

    document = f"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><title>RefAnno QC report</title>
<style>
body{{font-family:-apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.45}}
h1{{border-bottom:3px solid #4c78a8;padding-bottom:.3rem}}
h2{{margin-top:2.2rem;border-bottom:1px solid #ddd;padding-bottom:.2rem;color:#1a3a5c}}
h3{{margin-top:1.4rem;color:#333}}h4{{margin:.5rem 0}}
table{{border-collapse:collapse;margin:.8rem 0;min-width:45%}}
th,td{{border:1px solid #ccc;padding:.4rem .6rem;text-align:left;vertical-align:middle}}
th{{background:#eef3f8}}td.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
.note{{padding:.8rem 1rem;background:#fff8dc;border-left:4px solid #e6c229}}
.warnbox{{padding:.8rem 1rem;background:#fdecea;border-left:4px solid #d64545}}
.warn{{background:#ffe08a;padding:.1rem .3rem}}
.dim{{color:#666;font-size:.92em}}
.track{{display:inline-block;width:220px;height:.8em;background:#eee;border-radius:3px;vertical-align:middle;margin-right:.5rem}}
.fill{{height:100%;border-radius:3px}}
.cols{{display:flex;gap:2.5rem;flex-wrap:wrap}}
</style></head><body>
<h1>RefAnno QC report</h1>
<p class="note">这些指标描述证据链与输入模型命运。基因数增加或 BUSCO 完整度本身并不能证明结构准确性；
请优先关注第 2、7、9 节中的风险信号。</p>
{''.join(sections)}
<p class="dim">详细证据：decision_report.tsv · decision_events.tsv · junction_evidence.tsv · changes.tsv ·
family_rescues.tsv · input_fates.tsv · run_manifest.json</p>
</body></html>
"""
    with open(os.path.join(outdir, "qc_report.html"), "w") as out:
        out.write(document)
    return {"review_genes": len(review_winners), "structural_changes": len(changes)}


CONFIDENCE_BASE = {
    "isoseq_exact": 900,
    "short_read_complete": 800,
    "short_read_partial": 550,
    "homology_only": 400,
    "family_indirect": 300,
    "sequence_only": 250,
    "untested": 250,
    "annotation_only": 200,
}
CONFIDENCE_TIERS = [
    (750, "26,152,80"),
    (500, "230,171,2"),
    (300, "244,109,67"),
    (0, "150,150,150"),
]


def risk_score(model, reason, evidence_class):
    score = 0
    if "internal_stop" in model.review:
        score += 4
    if reason == "unique_with_protein" and len(model.cds) == 1 and model.depth == 0:
        score += 3
    if evidence_class == "family_indirect":
        score += 2
    if model.junctions_total and model.junctions_supported < model.junctions_total:
        score += 2
    if any(flag.startswith("utr_") for flag in model.review):
        score += 1
    if "weak_splice_sites" in model.review:
        score += 1
    if model.structure_status in {"partial", "unsupported", "untested"}:
        score += 1
    return score


def confidence_score(model, evidence_class):
    if "internal_stop" in model.review:
        return 50
    score = CONFIDENCE_BASE.get(evidence_class, 100)
    score -= 30 * (model.junctions_total - model.junctions_supported)
    score -= 50 * len(set(model.review))
    return max(0, min(1000, score))


def confidence_rgb(model, score):
    if "internal_stop" in model.review:
        return "215,48,39"
    for threshold, rgb in CONFIDENCE_TIERS:
        if score >= threshold:
            return rgb
    return CONFIDENCE_TIERS[-1][1]


def _bed_row(model, name, score, rgb):
    return "\t".join([
        model.chrom, str(model.exons[0][0] - 1), str(model.exons[-1][1]), name, str(score),
        model.strand, str(model.cds[0][0] - 1), str(model.cds[-1][1]), rgb,
    ])


def write_curation_tracks(outdir, winners, records):
    with open(os.path.join(outdir, "gene_confidence.bed"), "w") as confidence, \
            open(os.path.join(outdir, "curation_sites.bed"), "w") as curation:
        confidence.write('track name="RefAnno gene confidence" itemRgb=On visibility=dense\n')
        curation.write('track name="RefAnno curation sites" itemRgb=On\n')
        sites = []
        for winner in winners:
            model = winner["model"]
            record = records[model.key]
            score = confidence_score(model, winner["evidence_class"])
            rgb = confidence_rgb(model, score)
            confidence.write(_bed_row(model, f"{record['gene_id']}|{winner['evidence_class']}", score, rgb) + "\n")
            needs_review = (model.review or winner["reason"] == "tandem_family_support"
                            or model.structure_status in {"partial", "unsupported", "untested"})
            if not needs_review:
                continue
            if "internal_stop" in model.review:
                issue = "internal_stop"
            elif model.review:
                issue = "+".join(sorted(set(model.review)))
            else:
                issue = winner["evidence_class"]
            risk = risk_score(model, record["reason"], winner["evidence_class"])
            sites.append((model, record, score, rgb, issue, risk))
        sites.sort(key=lambda item: (-item[5], natural_key(item[0].chrom), item[0].exons[0][0], item[0].key))
        for model, record, score, rgb, issue, _ in sites:
            curation.write(_bed_row(model, f"{record['gene_id']}|{issue}", score, rgb) + "\n")


def output_tier(winner):
    if winner["evidence_class"] in {"isoseq_exact", "short_read_complete", "short_read_partial",
                                    "homology_only", "sequence_only"}:
        return "final"
    if winner["reason"].endswith("set_chain_consensus"):
        return "final"
    return "weak"


def write_outputs(winners, dropped, relations, models, sj, iso, genome, args,
                  bam_stats=None, protein_evidence=None):
    winners.sort(key=lambda w: (natural_key(w["model"].chrom), w["model"].cds[0][0],
                                w["model"].cds[-1][1], w["model"].strand, w["model"].key))
    duplicate_info = detect_duplicate_clusters([w["model"] for w in winners], args)
    records = {}
    retained = set()
    transcript_count = 0
    iso_loci = 0
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "final.gff3"), "w") as gff, \
            open(os.path.join(args.outdir, "final.proteins.faa"), "w") as proteins, \
            open(os.path.join(args.outdir, "weak.gff3"), "w") as weak_gff, \
            open(os.path.join(args.outdir, "weak.proteins.faa"), "w") as weak_proteins, \
            open(os.path.join(args.outdir, "decision_report.tsv"), "w", newline="") as report:
        gff.write("##gff-version 3\n")
        weak_gff.write("##gff-version 3\n")
        fields = ["gene_id", "transcript_id", "role", "chrom", "start", "end", "strand", "source_model",
                  "structure_status", "evidence_class", "structure_coverage", "junctions_supported", "junctions_total",
                  "canonical_junctions", "splice_sites_supported", "splice_min_score", "complete_cds",
                  "cds_depth", "cds_coverage", "cds_multi_depth",
                  "cds_multi_coverage", "family_support", "family_donors", "isoform_molecules", "protein_target",
                  "protein_identity", "protein_coverage", "protein_query_coverage", "duplicate_type",
                  "duplicate_cluster", "duplicate_size", "support", "reason", "selection", "review", "output_tier"]
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
            tier = output_tier(winner)
            out_gff, out_proteins = (gff, proteins) if tier == "final" else (weak_gff, weak_proteins)
            gene_attrs = [("ID", gene_id), ("gene_biotype", "protein_coding"),
                          ("support", "+".join(sorted(winner["sources"]))),
                          ("source_model", model.source_model), ("reason", winner["reason"]),
                          ("selection", winner["selection"]), ("structure_status", model.structure_status),
                          ("evidence_class", winner["evidence_class"]),
                          ("structure_coverage", f"{model.structure_coverage:.3f}"),
                          ("junctions_supported", model.junctions_supported),
                          ("junctions_total", model.junctions_total),
                          ("canonical_junctions", model.canonical_junctions),
                          ("splice_sites_supported", model.splice_supported),
                          ("complete_cds", int(model.complete_cds)), ("cds_depth", f"{model.depth:.2f}"),
                          ("cds_coverage", f"{model.coverage:.3f}"),
                          ("cds_multi_depth", f"{model.multi_depth:.2f}"),
                          ("cds_multi_coverage", f"{model.multi_coverage:.3f}"),
                          ("family_support", model.family_support),
                          ("family_donors", format_family_donors(model.family_donors)),
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
            out_gff.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], gene_attrs) + "\n")
            write_transcript(out_gff, out_proteins, model, gene_id, 1, genome, "primary", model.iso_support)
            transcript_count += 1
            records[model.key] = {
                "gene_id": gene_id, "transcript_id": f"{gene_id}-T1", "outcome": "primary",
                "reason": winner["reason"], "selection": winner["selection"],
                "evidence_class": winner["evidence_class"],
            }
            writer.writerow({"gene_id": gene_id, "transcript_id": f"{gene_id}-T1", "role": "primary",
                             "chrom": model.chrom, "start": model.exons[0][0], "end": model.exons[-1][1],
                             "strand": model.strand, "source_model": model.source_model,
                             "structure_status": model.structure_status,
                             "evidence_class": winner["evidence_class"],
                             "structure_coverage": f"{model.structure_coverage:.3f}",
                             "junctions_supported": model.junctions_supported,
                             "junctions_total": model.junctions_total,
                             "canonical_junctions": model.canonical_junctions,
                             "splice_sites_supported": model.splice_supported,
                             "splice_min_score": "" if model.splice_min_score is None else f"{model.splice_min_score:.1f}",
                             "complete_cds": int(model.complete_cds),
                             "cds_depth": f"{model.depth:.3f}", "cds_coverage": f"{model.coverage:.3f}",
                             "cds_multi_depth": f"{model.multi_depth:.3f}",
                             "cds_multi_coverage": f"{model.multi_coverage:.3f}",
                             "family_support": model.family_support,
                             "family_donors": format_family_donors(model.family_donors),
                             "isoform_molecules": model.iso_support, "protein_target": model.protein_target,
                             "protein_identity": f"{model.protein_identity:.3f}",
                             "protein_coverage": f"{model.protein_coverage:.3f}",
                             "protein_query_coverage": f"{model.protein_query_coverage:.3f}",
                             "duplicate_type": duplicate["type"] if duplicate else "",
                             "duplicate_cluster": duplicate["cluster"] if duplicate else "",
                             "duplicate_size": duplicate["size"] if duplicate else "",
                             "support": "+".join(sorted(winner["sources"])),
                             "reason": winner["reason"], "selection": winner["selection"],
                             "review": ",".join(sorted(set(model.review))), "output_tier": tier})
            for rank, (alternative, support) in enumerate(alternatives, 2):
                retained.add(alternative.key)
                alternative.output_id = gene_id
                write_transcript(out_gff, out_proteins, alternative, gene_id, rank, genome, "alternative", support)
                transcript_count += 1
                records[alternative.key] = {
                    "gene_id": gene_id, "transcript_id": f"{gene_id}-T{rank}", "outcome": "alternative",
                    "reason": "isoseq_exact_chain", "selection": "isoseq_chain_abundance",
                    "evidence_class": "isoseq_exact", "related_model": model.key,
                }
                writer.writerow({"gene_id": gene_id, "transcript_id": f"{gene_id}-T{rank}", "role": "alternative",
                                 "chrom": alternative.chrom, "start": alternative.exons[0][0],
                                 "end": alternative.exons[-1][1], "strand": alternative.strand,
                                 "source_model": alternative.source_model,
                                 "structure_status": alternative.structure_status,
                                 "evidence_class": "isoseq_exact",
                                 "structure_coverage": f"{alternative.structure_coverage:.3f}",
                                 "junctions_supported": alternative.junctions_supported,
                                 "junctions_total": alternative.junctions_total,
                                 "canonical_junctions": alternative.canonical_junctions,
                                 "splice_sites_supported": alternative.splice_supported,
                                 "splice_min_score": "" if alternative.splice_min_score is None else f"{alternative.splice_min_score:.1f}",
                                 "complete_cds": int(alternative.complete_cds),
                                 "cds_depth": f"{alternative.depth:.3f}", "cds_coverage": f"{alternative.coverage:.3f}",
                                 "cds_multi_depth": f"{alternative.multi_depth:.3f}",
                                 "cds_multi_coverage": f"{alternative.multi_coverage:.3f}",
                                 "family_support": alternative.family_support,
                                 "family_donors": format_family_donors(alternative.family_donors),
                                 "isoform_molecules": support, "protein_target": alternative.protein_target,
                                 "protein_identity": f"{alternative.protein_identity:.3f}",
                                 "protein_coverage": f"{alternative.protein_coverage:.3f}",
                                 "protein_query_coverage": f"{alternative.protein_query_coverage:.3f}",
                                 "duplicate_type": duplicate["type"] if duplicate else "",
                                 "duplicate_cluster": duplicate["cluster"] if duplicate else "",
                                 "duplicate_size": duplicate["size"] if duplicate else "",
                                 "support": alternative.set_name,
                                 "reason": "isoseq_exact_chain", "selection": "isoseq_chain_abundance",
                                 "review": "", "output_tier": tier})

    with open(os.path.join(args.outdir, "dropped.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        index = 0
        for model in sorted(models, key=lambda m: m.key):
            if model.key in retained:
                continue
            index += 1
            gid = f"D{index:06d}"
            reason = dropped.get(model.key, "not_selected")
            records[model.key] = {
                "gene_id": gid, "transcript_id": f"{gid}-T1", "outcome": "dropped",
                "reason": reason, "selection": "not_selected", "evidence_class": "not_retained",
            }
            attrs0 = [("ID", gid), ("source_model", model.source_model), ("reason", reason)]
            out.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], attrs0) + "\n")
            txid = gid + "-T1"
            out.write(model_row(model, "mRNA", model.exons[0][0], model.exons[-1][1],
                                [("ID", txid), ("Parent", gid), ("source_model", model.source_model)]) + "\n")
            for i, (start, end) in enumerate(model.exons, 1):
                out.write(model_row(model, "exon", start, end, [("ID", f"{txid}.exon.{i}"), ("Parent", txid)]) + "\n")
            for i, (start, end) in enumerate(model.cds, 1):
                out.write(model_row(model, "CDS", start, end, [("ID", f"{txid}.CDS.{i}"), ("Parent", txid)]) + "\n")

    write_curation_tracks(args.outdir, winners, records)
    write_decision_evidence(args.outdir, models, records, relations, sj, iso, genome, args)
    qc = write_review_qc(args.outdir, winners, models, records, relations, genome, args,
                         sj=sj, iso=iso, bam_stats=bam_stats, protein_evidence=protein_evidence,
                         duplicate_info=duplicate_info)

    summary = {
        "refanno_version": __version__, "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "input_transcripts": len(models), "genes": len(winners), "transcripts": transcript_count,
        "final_genes": sum(output_tier(w) == "final" for w in winners),
        "weak_genes": sum(output_tier(w) == "weak" for w in winners),
        "genes_with_alternatives": iso_loci, "alternative_transcripts": transcript_count - len(winners),
        "dropped_transcripts": len(models) - len(retained),
        "utr_consensus_applied": sum(bool(w["model"].utr_note) for w in winners),
        "utr_consensus_changed": sum(bool(w.get("utr_changed")) for w in winners),
        "protein_supported_transcripts": sum(m.protein_coverage > 0 for m in models),
        "protein_only_retained_genes": sum(w["reason"] == "unique_with_protein" for w in winners),
        "structure_supported_genes": sum(w["model"].structure_status in {"exact_isoseq", "supported", "sequence_supported"} for w in winners),
        "structure_sequence_supported_genes": sum(w["model"].structure_status == "sequence_supported" for w in winners),
        "structure_partial_genes": sum(w["model"].structure_status == "partial" for w in winners),
        "structure_untested_genes": sum(w["model"].structure_status == "untested" for w in winners),
        "review_genes": qc["review_genes"], "structural_changes": qc["structural_changes"],
        "evidence_classes": dict(Counter(w["evidence_class"] for w in winners)),
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
