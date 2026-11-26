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
    "structure_status", "structure_coverage", "junctions_supported", "junctions_total",
    "canonical_junctions", "complete_cds", "cds_depth", "cds_coverage", "cds_multi_depth",
    "cds_multi_coverage", "family_support", "isoform_molecules", "protein_target",
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
        self.assertEqual(summary["refanno_version"], "1.0.0")
        self.assertEqual(summary["output_schema_version"], "1.0")
        return rows, summary

    def assert_representative_and_family(self, rows, summary):
        by_source = {row["source_model"]: row for row in rows}
        self.assertEqual(summary["genes"], 4)
        self.assertIn("setA|short.t1", by_source)
        self.assertNotIn("setB|long.t1", by_source)
        short = by_source["setA|short.t1"]
        self.assertEqual(short["selection"], "rna_structure_depth")
        self.assertEqual(short["structure_status"], "supported")
        protein = by_source["setA|protein.t1"]
        self.assertEqual(protein["reason"], "unique_with_protein")
        self.assertEqual(protein["selection"], "protein_coverage")
        self.assertEqual(protein["structure_status"], "unsupported")
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

    def test_rna_support_outranks_protein_for_representative(self):
        build_single_exon(self.root)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + ["--miniprot-gff", str(self.root / "long_protein.gff3")])
        rows, summary = self.load_result(outdir)
        self.assertEqual(summary["genes"], 1)
        self.assertEqual(rows[0]["source_model"], "setA|short.t1")
        self.assertEqual(rows[0]["selection"], "rna_structure_depth")

    def test_miniprot_fasta_workflow(self):
        miniprot = os.environ.get("REFANNO_MINIPROT") or shutil.which("miniprot")
        if not miniprot:
            self.skipTest("miniprot is not available")
        build_single_exon(self.root)
        outdir = self.root / "out"
        self.run_command(self.command(self.root, outdir) + [
            "--proteins", str(self.root / "proteins.faa"), "--miniprot", miniprot,
        ])
        rows, summary = self.load_result(outdir)
        self.assert_representative_and_family(rows, summary)

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

    def test_version(self):
        result = self.run_command([sys.executable, str(self.script), "--version"])
        self.assertEqual(result.stdout.strip(), "RefAnno 1.0.0")


if __name__ == "__main__":
    unittest.main()
