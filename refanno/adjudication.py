from collections import Counter, defaultdict

from .models import blocks_coverage, natural_key, overlap, translate


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


def evidence_class(model, reason, args):
    if model.iso_support >= args.min_isoform_molecules:
        return "isoseq_exact"
    if model.rna_supported:
        return "short_read_complete" if model.structure_status == "supported" else "short_read_partial"
    if model.protein_supported:
        return "homology_only"
    if reason == "tandem_family_support":
        return "family_indirect"
    if model.structure_status == "untested":
        return "untested"
    return "annotation_only"


def choose_representatives(loci, sj, iso, proteins, args, trusted, tandem_info):
    winners = []
    dropped = {}
    relations = {}
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
                best.family_cluster = family["cluster"]
                best.family_donors = family["family_donors"]
                if "tandem_family_support" not in best.review:
                    best.review.append("tandem_family_support")
            else:
                for model in group:
                    dropped[model.key] = "unique_no_evidence"
                    relations[model.key] = {"rule": "retention"}
                continue
            kept.append((best, group, reason, sources))
        kept.sort(key=lambda item: representative_rank(item, args))
        surviving = []
        for item in kept:
            best = item[0]
            conflict = next((other for other in surviving
                             if best.strand == other[0].strand
                             and max(reciprocal_cds_overlap(best, other[0])) > args.conflict_overlap), None)
            if conflict is not None:
                for model in item[1]:
                    dropped[model.key] = "locus_conflict_lower_rank"
                    relations[model.key] = {
                        "rule": "locus_conflict", "related_model": conflict[0].key,
                        "overlap": reciprocal_cds_overlap(model, conflict[0]),
                    }
            else:
                surviving.append(item)
        for best, group, reason, sources in surviving:
            for model in group:
                if model is not best:
                    dropped[model.key] = "same_chain_lower_rank"
                    relations[model.key] = {"rule": "same_chain", "related_model": best.key}
            _, selection = representative_basis((best, group, reason, sources), args)
            winners.append({"model": best, "group": group, "reason": reason, "selection": selection,
                            "evidence_class": evidence_class(best, reason, args),
                            "sources": sources, "locus": locus, "alternatives": []})
    return winners, dropped, relations


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


def protein_similarity(a, b):
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, left in enumerate(a, 1):
        current = [i]
        for j, right in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (left != right)))
        previous = current
    return 1.0 - previous[-1] / max(len(a), len(b))


def related_family_members(model, members, genome, args, sequence_cache, similarity_cache):
    related = []
    for other in members:
        if other.key == model.key:
            continue
        pair = tuple(sorted((model.key, other.key)))
        if pair not in similarity_cache:
            same_target = bool(model.protein_target and model.protein_target == other.protein_target)
            if same_target:
                similarity_cache[pair] = True
            else:
                for item in (model, other):
                    if item.key not in sequence_cache:
                        sequence_cache[item.key] = translate(genome.cds_sequence(item)).rstrip("*")
                first, second = sequence_cache[model.key], sequence_cache[other.key]
                longer, shorter = max(len(first), len(second)), min(len(first), len(second))
                # 编辑距离至少是长度差，长度差已超出阈值允许的范围时不必算编辑距离
                if longer - shorter > (1.0 - args.tandem_min_protein_similarity) * longer:
                    similarity_cache[pair] = False
                else:
                    similarity = protein_similarity(first, second)
                    similarity_cache[pair] = similarity >= args.tandem_min_protein_similarity
        if similarity_cache[pair]:
            related.append(other)
    return related


def prepare_tandem_evidence(loci, sj, proteins, genome, args):
    grouped_loci = [chain_groups(locus) for locus in loci]
    representatives = [best_in_group(group) for groups in grouped_loci for group in groups]
    annotations = detect_duplicate_clusters(representatives, args)
    cluster_members = defaultdict(list)
    for model in representatives:
        if model.key in annotations:
            cluster_members[annotations[model.key]["cluster"]].append(model)
    support_by_model = {}
    donors_by_model = {}
    sequence_cache = {}
    similarity_cache = {}
    for cluster, members in cluster_members.items():
        for model in members:
            related = related_family_members(model, members, genome, args, sequence_cache, similarity_cache)
            if not related:
                continue
            support = []
            donors = []
            for donor in [model] + related:
                donor_support = []
                if donor.key != model.key and model_has_protein(donor, proteins, args):
                    donor_support.append(("protein", donor.protein_target, donor.protein_coverage))
                if donor.multi_coverage >= 0.5 and donor.multi_depth >= args.tandem_min_multi_depth:
                    donor_support.append(("multi_depth", "", donor.multi_depth))
                multi_junction = max((sj.get(donor, intron)[1] for intron in donor.cds_introns), default=0)
                if multi_junction >= 3 * args.min_sj:
                    donor_support.append(("multi_junction", "", multi_junction))
                for kind, target, value in donor_support:
                    if kind not in support:
                        support.append(kind)
                    donors.append({"donor": donor.key, "evidence": kind,
                                   "protein_target": target, "value": value})
            support_by_model[model.key] = "+".join(support)
            donors_by_model[model.key] = donors
    tandem_info = {}
    for groups in grouped_loci:
        for group in groups:
            representative = best_in_group(group)
            annotation = annotations.get(representative.key)
            if not annotation:
                continue
            annotation = dict(annotation)
            annotation["family_support"] = support_by_model.get(representative.key, "")
            annotation["family_donors"] = donors_by_model.get(representative.key, [])
            for model in group:
                tandem_info[model.key] = annotation
    return tandem_info
