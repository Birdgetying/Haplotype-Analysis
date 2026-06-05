#!/usr/bin/env python3
"""Tests for lightweight star-gene data source metadata."""

from pathlib import Path
from contextlib import redirect_stdout
from io import StringIO
import tempfile
import unittest


class StarGeneDataTests(unittest.TestCase):
    def test_rice_figshare_is_large_and_skipped_by_default(self):
        from star_gene_data import build_download_commands, iter_data_files

        rice_files = list(iter_data_files(paper="rice2024"))
        rice_figshare = [f for f in rice_files if f.key == "rice2024_figshare_genotype_matrix"][0]

        self.assertTrue(rice_figshare.is_large)
        self.assertEqual(rice_figshare.default_action, "manual_large_download")
        self.assertIn("10.6084/m9.figshare.19166475", rice_figshare.source)
        self.assertEqual(build_download_commands(paper="rice2024", include_large=False), [])

    def test_maize_small_files_have_direct_download_commands(self):
        from star_gene_data import build_download_commands, iter_data_files

        maize_files = list(iter_data_files(paper="maize2019"))
        keys = {f.key for f in maize_files}

        self.assertIn("maize2019_sv_382254", keys)
        self.assertIn("maize2019_agronomic_blup", keys)

        commands = build_download_commands(paper="maize2019")
        command_text = "\n".join(commands)
        self.assertIn("SV.382254.zip", command_text)
        self.assertIn("blup_traits_final.csv", command_text)
        self.assertIn("external_data/maize_natgenet_2019/maizego", command_text)

    def test_wheat_sources_are_instruction_only_until_object_names_are_known(self):
        from star_gene_data import build_download_commands, iter_data_files

        wheat_files = list(iter_data_files(paper="wheat2024"))
        self.assertTrue(wheat_files)
        portal_files = [f for f in wheat_files if f.key in {"wheat2024_wwwg2b_portal", "wheat2024_earlham_watseq"}]
        self.assertEqual(len(portal_files), 2)
        self.assertTrue(all(f.default_action == "instruction_only" for f in portal_files))
        self.assertEqual(build_download_commands(paper="wheat2024"), [])
        self.assertTrue(any("wwwg2b.com" in f.source for f in wheat_files))

    def test_summarize_downloads_contains_status_and_target_paths(self):
        from star_gene_data import summarize_downloads

        summary = summarize_downloads(paper="maize2019")
        expected_path = str(Path("external_data/maize_natgenet_2019/maizego/SV.382254.zip")).replace("\\", "/")
        self.assertIn("maize2019_sv_382254", summary)
        self.assertIn("small_direct_download", summary)
        self.assertIn(expected_path, summary)

    def test_unknown_paper_filter_raises_clear_error(self):
        from star_gene_data import iter_data_files

        with self.assertRaisesRegex(ValueError, "No data files matched"):
            list(iter_data_files(paper="not_a_paper"))

    def test_cli_print_downloads_exits_before_validation(self):
        from star_gene_validation import main

        stdout = StringIO()
        with redirect_stdout(stdout):
            rc = main(["--print-downloads", "--paper", "maize2019"])
        self.assertEqual(rc, 0)
        self.assertIn("maize2019_sv_382254", stdout.getvalue())

    def test_build_database_from_wide_marker_matrix(self):
        from star_gene_data import build_database_from_marker_matrix

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            marker_path = tmp_path / "markers.tsv"
            phenotype_path = tmp_path / "phenotypes.tsv"
            out_root = tmp_path / "db"

            marker_path.write_text(
                "SampleID\tqHKW1_100\tqHKW1_200\n"
                "S1\tA\tT\n"
                "S2\tA\tT\n"
                "S3\tG\tT\n"
                "S4\tG\tC\n",
                encoding="utf-8",
            )
            phenotype_path.write_text(
                "SampleID\tHKW\n"
                "S1\t25.1\n"
                "S2\t25.3\n"
                "S3\t30.4\n"
                "S4\t31.0\n",
                encoding="utf-8",
            )

            db_dir = build_database_from_marker_matrix(
                marker_matrix=marker_path,
                phenotype_table=phenotype_path,
                output_root=out_root,
                target_id="qHKW1",
                chrom="chr1",
                start=100,
                end=200,
                phenotype_columns=["HKW"],
                marker_columns=["qHKW1_100", "qHKW1_200"],
                marker_positions={"qHKW1_100": 100, "qHKW1_200": 200},
                expected_direction="increases_trait",
            )

            self.assertTrue((db_dir / "gene_info.json").exists())
            self.assertTrue((db_dir / "haplotype_data.csv").exists())
            self.assertTrue((db_dir / "haplotype_samples.csv").exists())
            self.assertTrue((db_dir / "phenotype_data.csv").exists())
            self.assertTrue((db_dir / "variant_info.csv").exists())

            hap_samples = (db_dir / "haplotype_samples.csv").read_text(encoding="utf-8")
            self.assertIn("S1,A|T,Hap1", hap_samples)
            self.assertIn("S3,G|T,Hap2", hap_samples)
            phenotype_data = (db_dir / "phenotype_data.csv").read_text(encoding="utf-8")
            self.assertIn("HKW", phenotype_data)

    def test_marker_database_cli(self):
        from build_star_gene_database import main as build_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            marker_path = tmp_path / "markers.tsv"
            phenotype_path = tmp_path / "phenotypes.tsv"
            out_root = tmp_path / "db"
            marker_path.write_text(
                "SampleID\tm1\nS1\tA\nS2\tG\n",
                encoding="utf-8",
            )
            phenotype_path.write_text(
                "SampleID\tHKW\nS1\t25\nS2\t30\n",
                encoding="utf-8",
            )

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = build_main([
                    "--marker-matrix", str(marker_path),
                    "--phenotype-table", str(phenotype_path),
                    "--output-root", str(out_root),
                    "--target-id", "qHKW1",
                    "--chrom", "chr1",
                    "--start", "100",
                    "--end", "100",
                    "--phenotype-column", "HKW",
                    "--marker-column", "m1",
                    "--marker-position", "m1=100",
                ])

            self.assertEqual(rc, 0)
            self.assertTrue((out_root / "qHKW1" / "gene_info.json").exists())
            self.assertIn("Built star-gene database", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
