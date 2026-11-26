from .adjudication import model_rank
from .models import blocks_coverage, translate


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
        return None, 0, 0.0, []
    best = max(clusters, key=lambda c: (len(c), -c[0].start, c[0].read_id))
    n = len(best[0].exons)
    medians = []
    for i in range(n):
        medians.append((sorted(r.exons[i][0] for r in best)[len(best) // 2],
                        sorted(r.exons[i][1] for r in best)[len(best) // 2]))
    medoid = min(best, key=lambda r: (sum(abs(r.exons[i][0] - medians[i][0]) +
                                             abs(r.exons[i][1] - medians[i][1]) for i in range(n)), r.read_id))
    return medoid, len(best), len(best) / len(reads), best


def apply_utr_consensus(winners, iso, args):
    for winner in winners:
        model = winner["model"]
        reads = iso.utr_compatible(model, args.min_cds_coverage)
        if len(reads) < args.utr_min_molecules:
            continue
        best, support, dominance, members = consensus_exons(reads)
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
        if winner["utr_changed"]:
            winner["change"] = {
                "type": "UTR_consensus", "old_exons": old, "new_exons": list(model.exons),
                "molecules": [read.read_id for read in members], "support": support,
                "dominance": dominance,
            }


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
