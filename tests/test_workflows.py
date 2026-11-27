import contextlib
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.synthetic import build_multiexon, build_single_exon

REPORT_FIELDS = [
    "gene_id", "transcript_id", "role", "chrom", "start", "end", "strand", "source_model",
    "structure_status", "evidence_class", "structure_coverage", "junctions_supported", "junctions_total",
    "canonical_junctions", "splice_sites_supported", "splice_min_score", "complete_cds",
    "cds_depth", "cds_coverage", "cds_multi_depth",
    "cds_multi_coverage", "family_support", "family_donors", "isoform_molecules", "protein_target",
    "protein_identity", "protein_coverage", "protein_query_coverage", "duplicate_type",
    "duplicate_cluster", "duplicate_size", "support", "reason", "selection", "review",
]


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="refanno-test-", dir=Path(__file__).parent)
        self.root = Path(self.temp.name)
        self.repo = Path(__file__).resolve().parents[1]
        self.script = self.repo / "RefAnno.py"

    def tearDown(self):
        self.temp.cleanup()

    def command(self, root, outdir):
        return [
            sys.executable, str(self.script), "--set", f"setA={root / 'setA.gff3'}",
            "--set", f"setB={root / 'setB.gff3'}", "--genome", str(root / "genome.fa"),
            "--bam", str(root / "rna.bam"), "--threads", "2", "--outdir", str(outdir),
        ]

    def run_command(self, command):
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def load_result(self, outdir):
        with (outdir / "decision_report.tsv").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, REPORT_FIELDS)
            rows = list(reader)
        with (outdir / "summary.json").open() as handle:
            summary = json.load(handle)
        self.assertEqual(summary["refanno_version"], "1.2.0")
        self.assertEqual(summary["output_schema_version"], "1.2")
        return rows, summary

    def read_tsv(self, path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def assert_manifest(self, outdir):
        manifest = json.loads((outdir / "run_manifest.json").read_text())
        self.assertEqual(manifest["status"], "completed")
        self.assertEqual(manifest["refanno_version"], "1.2.0")
        self.assertEqual(manifest["output_schema_version"], "1.2")
        self.assertTrue(manifest["inputs"])
        self.assertTrue(all(len(record["sha256"]) == 64 for record in manifest["inputs"]))
        for record in manifest["inputs"]:
            if record["role"] == "short_read_bam":
                self.assertIn("reads", record)
                self.assertIn("mapq_pass_reads", record)
                self.assertIn("spliced_reads", record)
                self.assertGreaterEqual(record["mapq_pass_reads"], record["spliced_reads"])
        output_names = {record["path"] for record in manifest["outputs"]}
        self.assertIn("final.gff3", output_names)
        self.assertIn("decision_events.tsv", output_names)
        self.assertIn("junction_evidence.tsv", output_names)
        self.assertIn("qc_report.html", output_names)
        return manifest

    def assert_representative_and_family(self, rows, summary):
        by_source = {row["source_model"]: row for row in rows}
        self.assertEqual(summary["genes"], 4)
        self.assertIn("setA|short.t1", by_source)
        self.assertNotIn("setB|long.t1", by_source)
        short = by_source["setA|short.t1"]
        self.assertEqual(short["selection"], "rna_structure_depth")
        self.assertEqual(short["structure_status"], "supported")
        self.assertEqual(short["evidence_class"], "short_read_complete")
        protein = by_source["setA|protein.t1"]
        self.assertEqual(protein["reason"], "unique_with_protein")
        self.assertEqual(protein["selection"], "protein_coverage")
        self.assertEqual(protein["structure_status"], "unsupported")
        self.assertEqual(protein["evidence_class"], "homology_only")
        dup_reasons = {
            name: by_source[f"{source}|{name}.t1"]["reason"]
            for source, name in (("setA", "dup1"), ("setB", "dup2"))
        }
        self.assertEqual(sorted(dup_reasons.values()), ["tandem_family_support", "unique_with_protein"])
        self.assertEqual(summary["tandem_family_retained_genes"], 1)
        self.assertEqual(summary["duplicate_genes"], 2)
        self.assertNotIn("setA|isolated.t1", by_source)
        self.assertNotIn("score", rows[0])
        self.assertEqual(summary["representative_basis"]["rna_structure_depth"], 1)

    def test_representative_protein_retention_and_duplicate_rescue(self):
        build_single_exon(self.root)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + ["--miniprot-gff", str(self.root / "evidence.gff3")])
        rows, summary = self.load_result(outdir)
        self.assert_representative_and_family(rows, summary)
        gff = (outdir / "final.gff3").read_text()
        self.assertIn("\tRefAnno\tgene\t", gff)
        self.assertIn("ID=REFA000001", gff)
        self.assertIn("selection=rna_structure_depth", gff)
        self.assertIn("structure_status=", gff)
        self.assert_manifest(outdir)
        rescues = self.read_tsv(outdir / "family_rescues.tsv")
        self.assertEqual(len(rescues), 1)
        self.assertEqual(rescues[0]["rescued_source_model"], "setA|dup1.t1")
        self.assertEqual(rescues[0]["donor_source_model"], "setB|dup2.t1")
        self.assertEqual(rescues[0]["cluster_id"], "DUPLICATE000001")
        self.assertEqual(rescues[0]["evidence"], "protein")
        events = self.read_tsv(outdir / "decision_events.tsv")
        long_event = next(row for row in events if row["source_model"] == "setB|long.t1")
        self.assertEqual(long_event["outcome"], "dropped")
        self.assertEqual(long_event["reason"], "unique_no_evidence")
        review = (outdir / "review.gff3").read_text()
        self.assertIn("source_model=setA|dup1.t1", review)
        self.assertIn("source_model=setA|protein.t1", review)
        self.assertNotIn("source_model=setA|short.t1", review)
        self.assertIn("RefAnno QC report", (outdir / "qc_report.html").read_text())
        confidence = [line.split("\t") for line in (outdir / "gene_confidence.bed").read_text().splitlines()
                      if not line.startswith("track")]
        self.assertEqual(len(confidence), 4)
        for fields in confidence:
            self.assertEqual(len(fields), 9)
            self.assertTrue(0 <= int(fields[4]) <= 1000)
            self.assertRegex(fields[8], r"^\d+,\d+,\d+$")
            self.assertEqual(int(fields[1]), int(fields[6]))
        sites = [line.split("\t") for line in (outdir / "curation_sites.bed").read_text().splitlines()
                 if not line.startswith("track")]
        site_names = {fields[3].split("|")[1] for fields in sites}
        self.assertEqual(site_names, {"tandem_family_support", "homology_only"})
        self.assertNotIn("short.t1", {fields[3] for fields in sites})
        fates = self.read_tsv(outdir / "input_fates.tsv")
        self.assertEqual({row["set_name"]: row["primary"] for row in fates}, {"setA": "3", "setB": "1"})

    def test_unrelated_adjacent_models_do_not_share_family_evidence(self):
        build_single_exon(self.root, related_duplicates=False)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + ["--miniprot-gff", str(self.root / "evidence.gff3")])
        rows, summary = self.load_result(outdir)
        self.assertEqual(summary["genes"], 3)
        self.assertEqual(summary["tandem_family_retained_genes"], 0)
        self.assertNotIn("setA|dup1.t1", {row["source_model"] for row in rows})
        self.assertEqual(self.read_tsv(outdir / "family_rescues.tsv"), [])
        events = self.read_tsv(outdir / "decision_events.tsv")
        dup1 = next(row for row in events if row["source_model"] == "setA|dup1.t1")
        self.assertEqual(dup1["reason"], "unique_no_evidence")

    def test_rna_support_outranks_protein_for_representative(self):
        build_single_exon(self.root)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + ["--miniprot-gff", str(self.root / "long_protein.gff3")])
        rows, summary = self.load_result(outdir)
        self.assertEqual(summary["genes"], 1)
        self.assertEqual(rows[0]["source_model"], "setA|short.t1")
        self.assertEqual(rows[0]["selection"], "rna_structure_depth")
        events = self.read_tsv(outdir / "decision_events.tsv")
        long_event = next(row for row in events if row["source_model"] == "setB|long.t1")
        self.assertEqual(long_event["reason"], "locus_conflict_lower_rank")
        self.assertEqual(long_event["related_source_model"], "setA|short.t1")
        self.assertEqual(long_event["related_gene_id"], "REFA000001")
        self.assertTrue(long_event["cds_overlap_source"])

    def test_miniprot_fasta_workflow(self):
        miniprot = os.environ.get("REFANNO_MINIPROT") or shutil.which("miniprot")
        if not miniprot:
            self.skipTest("miniprot is not available")
        build_single_exon(self.root)
        outdir = self.root / "out"
        command = self.command(self.root, outdir) + [
            "--proteins", str(self.root / "proteins.faa"), "--miniprot", miniprot,
        ]
        self.run_command(command)
        rows, summary = self.load_result(outdir)
        self.assert_representative_and_family(rows, summary)
        cache_before = json.loads((outdir / "miniprot.cache.json").read_text())
        reused = self.run_command(command)
        self.assertIn("miniprot: reuse", reused.stdout)
        with (self.root / "proteins.faa").open("a") as out:
            out.write(">changed_input\nMKTAA\n")
        rerun = self.run_command(command)
        self.assertNotIn("miniprot: reuse", rerun.stdout)
        cache_after = json.loads((outdir / "miniprot.cache.json").read_text())
        self.assertNotEqual(cache_before["proteins_sha256"], cache_after["proteins_sha256"])

    def test_multiexon_junction_isoseq_and_utr_consensus(self):
        build_multiexon(self.root)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + [
            "--isoseq", str(self.root / "iso.bed12"), "--utr-consensus",
        ])
        rows, summary = self.load_result(outdir)
        self.assertEqual(summary["genes"], 1)
        self.assertEqual(summary["transcripts"], 2)
        self.assertEqual(summary["genes_with_alternatives"], 1)
        self.assertEqual(summary["utr_consensus_applied"], 1)
        self.assertEqual(summary["utr_consensus_changed"], 1)
        self.assertEqual(summary["structural_changes"], 1)
        primary = next(row for row in rows if row["role"] == "primary")
        alternative = next(row for row in rows if row["role"] == "alternative")
        self.assertEqual(primary["source_model"], "setA|txsetA")
        self.assertEqual(primary["structure_status"], "exact_isoseq")
        self.assertEqual(primary["junctions_supported"], "1")
        self.assertEqual(primary["junctions_total"], "1")
        self.assertEqual(primary["canonical_junctions"], "1")
        self.assertEqual(primary["selection"], "isoseq_chain_abundance")
        self.assertGreater(float(primary["cds_depth"]), float(alternative["cds_depth"]))
        self.assertEqual(alternative["source_model"], "setB|txsetB")
        self.assertEqual(alternative["isoform_molecules"], "2")
        changes = self.read_tsv(outdir / "changes.tsv")
        self.assertEqual(len(changes), 1)
        self.assertEqual(changes[0]["change_type"], "UTR_consensus")
        self.assertEqual(changes[0]["cds_changed"], "0")
        self.assertEqual(changes[0]["supporting_molecules"], "A1|A2")
        self.assertEqual(changes[0]["molecule_count"], "2")
        self.assertEqual(changes[0]["old_exons"], "1001-1150,2001-2150")
        self.assertEqual(changes[0]["new_exons"], "901-1150,2001-2200")
        junctions = self.read_tsv(outdir / "junction_evidence.tsv")
        self.assertEqual(len(junctions), 2)
        primary_junction = next(row for row in junctions if row["source_model"] == "setA|txsetA")
        self.assertEqual(primary_junction["unique_reads"], "6")
        self.assertEqual(primary_junction["isoseq_molecules"], "2")
        self.assertEqual(primary_junction["motif"], "+")
        self.assertEqual(primary_junction["supported"], "1")
        bed_lines = (outdir / "supporting_junctions.bed").read_text().strip().splitlines()
        self.assertEqual(len(bed_lines), 2)
        self.assert_manifest(outdir)

    def test_depth_index_accumulator(self):
        from refanno.evidence import DepthIndex, Genome, SJIndex, load_bam_junctions
        from refanno.models import Model
        from tests.synthetic import make_read, write_bam, write_fasta

        write_fasta(self.root / "genome.fa", "A" * 2000)
        multi = make_read("r2", 151, [(0, 200)])
        multi.set_tag("NH", 3)
        secondary = make_read("r4", 91, [(0, 160)])
        secondary.is_secondary = True
        supplementary = make_read("r5", 91, [(0, 160)])
        supplementary.is_supplementary = True
        reads = [make_read("r1", 91, [(0, 160)]), multi, make_read("r3", 291, [(0, 60), (3, 100), (0, 60)]),
                 secondary, supplementary]
        write_bam(self.root / "rna.bam", reads, 2000)

        def model(name, blocks):
            item = Model("setA", name, name, "chr1", "+")
            item.cds = list(blocks)
            item.finish()
            return item

        both = model("both", [(101, 200), (301, 400)])
        shared = model("shared", [(101, 200)])
        tail = model("tail", [(501, 600)])
        bam = str(self.root / "rna.bam")
        depth = DepthIndex([bam], [both, shared, tail])
        sj = SJIndex()
        stats = load_bam_junctions([bam], Genome(str(self.root / "genome.fa")), sj, depth)
        self.assertEqual(stats, [{"path": bam, "reads": 3, "mapq_pass": 3, "spliced_reads": 1}])
        self.assertEqual(depth.get(both), (0.75, 0.75, 0.5, 0.5))
        self.assertEqual(depth.get(shared), (1.0, 1.0, 0.5, 0.5))
        self.assertEqual(depth.get(tail), (0.1, 0.1, 0.0, 0.0))
        self.assertEqual(dict(sj.data["chr1"]), {(351, 450, "."): [1, 0]})

        flat = str(self.root / "hic.bam")
        write_bam(flat, [make_read(f"f{i}", 1 + i % 1800, [(0, 100)]) for i in range(10000)], 2000)
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            flat_stats = load_bam_junctions([flat], Genome(str(self.root / "genome.fa")), SJIndex())
        self.assertEqual(flat_stats[0]["reads"], 10000)
        self.assertEqual(flat_stats[0]["spliced_reads"], 0)
        self.assertIn("可能不是 RNA-seq", stderr.getvalue())

    def test_splice_scores_validation(self):
        build_multiexon(self.root)
        outdir = self.root / "out"
        with (self.root / "scores.tsv").open("w") as out:
            out.write("chr1\t1150\t+\tD\t9\nchr1\t3000\t+\tA\t11\nchr1\t2000\t+\tA\t0\n")
        self.run_command([
            sys.executable, str(self.script),
            "--set", f"setA={self.root / 'setA.gff3'}", "--set", f"setB={self.root / 'setB.gff3'}",
            "--genome", str(self.root / "genome.fa"),
            "--trusted-set", "setA", "--trusted-set", "setB",
            "--splice-scores", str(self.root / "scores.tsv"),
            "--threads", "2", "--outdir", str(outdir),
        ])
        rows, summary = self.load_result(outdir)
        self.assertEqual(summary["genes"], 1)
        self.assertEqual(summary["structure_sequence_supported_genes"], 1)
        winner = rows[0]
        self.assertEqual(winner["source_model"], "setB|txsetB")
        self.assertEqual(winner["structure_status"], "sequence_supported")
        self.assertEqual(winner["evidence_class"], "sequence_only")
        self.assertEqual(winner["splice_sites_supported"], "1")
        self.assertEqual(winner["splice_min_score"], "9.0")
        self.assertEqual(winner["reason"], "trusted_unique")
        events = self.read_tsv(outdir / "decision_events.tsv")
        weak = next(row for row in events if row["source_model"] == "setA|txsetA")
        self.assertEqual(weak["outcome"], "dropped")
        self.assertIn("2.6 剪接位点序列校验", (outdir / "qc_report.html").read_text())
        confidence = [line.split("\t") for line in (outdir / "gene_confidence.bed").read_text().splitlines()
                      if not line.startswith("track")]
        self.assertEqual(confidence[0][3], f"{winner['gene_id']}|sequence_only")

    def test_version(self):
        result = self.run_command([sys.executable, str(self.script), "--version"])
        self.assertEqual(result.stdout.strip(), "RefAnno 1.2.0")


if __name__ == "__main__":
    unittest.main()
