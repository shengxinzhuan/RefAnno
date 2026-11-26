import csv
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
    "canonical_junctions", "complete_cds", "cds_depth", "cds_coverage", "cds_multi_depth",
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
        self.assertEqual(summary["refanno_version"], "1.1.0")
        self.assertEqual(summary["output_schema_version"], "1.1")
        return rows, summary

    def read_tsv(self, path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def assert_manifest(self, outdir):
        manifest = json.loads((outdir / "run_manifest.json").read_text())
        self.assertEqual(manifest["status"], "completed")
        self.assertEqual(manifest["refanno_version"], "1.1.0")
        self.assertEqual(manifest["output_schema_version"], "1.1")
        self.assertTrue(manifest["inputs"])
        self.assertTrue(all(len(record["sha256"]) == 64 for record in manifest["inputs"]))
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

    def test_version(self):
        result = self.run_command([sys.executable, str(self.script), "--version"])
        self.assertEqual(result.stdout.strip(), "RefAnno 1.1.0")


if __name__ == "__main__":
    unittest.main()
