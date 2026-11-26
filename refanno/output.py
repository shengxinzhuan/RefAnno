import csv
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
    protein = translate(genome.cds_sequence(model))
    protein_handle.write(f">{tx_id} source_model={model.source_model} role={role}\n")
    for i in range(0, len(protein), 60):
        protein_handle.write(protein[i:i + 60] + "\n")


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
