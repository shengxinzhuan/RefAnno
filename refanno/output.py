import csv
import html
import json
import os
from collections import Counter

from . import OUTPUT_SCHEMA_VERSION, __version__
from .adjudication import detect_duplicate_clusters
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


def write_review_qc(outdir, winners, models, records, relations, genome, args):
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
                "structure_supported": sum(model.structure_status in {"exact_isoseq", "supported"} for model in members),
                "structure_partial": sum(model.structure_status == "partial" for model in members),
                "primary": fates[(set_name, "primary")], "alternative": fates[(set_name, "alternative")],
                "dropped": fates[(set_name, "dropped")],
            })

    evidence_counts = Counter(winner["evidence_class"] for winner in winners)
    status_counts = Counter(winner["model"].structure_status for winner in winners)
    rows = [
        ("Primary genes", len(winners)),
        ("Transcripts", sum(1 + len(winner["alternatives"]) for winner in winners)),
        ("Structural changes", len(changes)),
        ("Genes requiring review", len(review_winners)),
        ("Evidence classes", ", ".join(f"{key}={value}" for key, value in sorted(evidence_counts.items()))),
        ("Structure status", ", ".join(f"{key}={value}" for key, value in sorted(status_counts.items()))),
    ]
    fate_rows = []
    for set_name in sorted({model.set_name for model in models}):
        fate_rows.append((set_name, fates[(set_name, "primary")], fates[(set_name, "alternative")],
                          fates[(set_name, "dropped")]))
    overview = "".join(f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>" for key, value in rows)
    fate_html = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
                        for row in fate_rows)
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>RefAnno QC report</title>
<style>body{{font-family:sans-serif;max-width:1100px;margin:2rem auto;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;min-width:45%}}th,td{{border:1px solid #ccc;padding:.45rem;text-align:left}}
th{{background:#eef3f8}}.note{{padding:1rem;background:#fff8dc}}</style></head><body>
<h1>RefAnno QC report</h1>
<p class="note">These metrics describe evidence and model fates. Increased gene counts or BUSCO completeness alone do not demonstrate structural accuracy.</p>
<h2>Run overview</h2><table>{overview}</table>
<h2>Input model fates</h2><table><tr><th>Set</th><th>Primary</th><th>Alternative</th><th>Dropped</th></tr>{fate_html}</table>
<p>Detailed evidence: decision_report.tsv, decision_events.tsv, junction_evidence.tsv, changes.tsv and run_manifest.json.</p>
</body></html>
"""
    with open(os.path.join(outdir, "qc_report.html"), "w") as out:
        out.write(document)
    return {"review_genes": len(review_winners), "structural_changes": len(changes)}


def write_outputs(winners, dropped, relations, models, sj, iso, genome, args):
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
            open(os.path.join(args.outdir, "decision_report.tsv"), "w", newline="") as report:
        gff.write("##gff-version 3\n")
        fields = ["gene_id", "transcript_id", "role", "chrom", "start", "end", "strand", "source_model",
                  "structure_status", "evidence_class", "structure_coverage", "junctions_supported", "junctions_total",
                  "canonical_junctions", "complete_cds", "cds_depth", "cds_coverage", "cds_multi_depth",
                  "cds_multi_coverage", "family_support", "family_donors", "isoform_molecules", "protein_target",
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
                          ("evidence_class", winner["evidence_class"]),
                          ("structure_coverage", f"{model.structure_coverage:.3f}"),
                          ("junctions_supported", model.junctions_supported),
                          ("junctions_total", model.junctions_total),
                          ("canonical_junctions", model.canonical_junctions),
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
            gff.write(model_row(model, "gene", model.exons[0][0], model.exons[-1][1], gene_attrs) + "\n")
            write_transcript(gff, proteins, model, gene_id, 1, genome, "primary", model.iso_support)
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
                             "review": ",".join(sorted(set(model.review)))})
            for rank, (alternative, support) in enumerate(alternatives, 2):
                retained.add(alternative.key)
                alternative.output_id = gene_id
                write_transcript(gff, proteins, alternative, gene_id, rank, genome, "alternative", support)
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

    write_decision_evidence(args.outdir, models, records, relations, sj, iso, genome, args)
    qc = write_review_qc(args.outdir, winners, models, records, relations, genome, args)

    summary = {
        "refanno_version": __version__, "output_schema_version": OUTPUT_SCHEMA_VERSION,
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
