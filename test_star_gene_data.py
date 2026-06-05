#!/usr/bin/env python3
"""Tests for lightweight star-gene data source metadata."""

from pathlib import Path
from contextlib import redirect_stdout
from io import StringIO
import json
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

    def test_build_database_filters_configured_phenotype_missing_values(self):
        from star_gene_data import build_database_from_marker_matrix

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            marker_path = tmp_path / "markers.tsv"
            phenotype_path = tmp_path / "phenotypes.tsv"
            out_root = tmp_path / "db"

            marker_path.write_text(
                "SampleID\tm1\n"
                "S1\tNN\n"
                "S2\tAA\n"
                "S3\tAA\n",
                encoding="utf-8",
            )
            phenotype_path.write_text(
                "SampleID\tHKW\n"
                "S1\t-999.0\n"
                "S2\t25.1\n"
                "S3\t25.3\n",
                encoding="utf-8",
            )

            db_dir = build_database_from_marker_matrix(
                marker_matrix=marker_path,
                phenotype_table=phenotype_path,
                output_root=out_root,
                target_id="qHKW1",
                chrom="chr1",
                start=100,
                end=100,
                phenotype_columns=["HKW"],
                marker_columns=["m1"],
                marker_positions={"m1": 100},
                phenotype_missing_values=["-999"],
            )

            phenotype_data = (db_dir / "phenotype_data.csv").read_text(encoding="utf-8")
            hap_samples = (db_dir / "haplotype_samples.csv").read_text(encoding="utf-8")

            self.assertNotIn("S1", phenotype_data)
            self.assertNotIn("-999.0", phenotype_data)
            self.assertIn("S2", phenotype_data)
            self.assertNotIn("S1", hap_samples)

    def test_build_database_infers_sv_length_from_marker_id(self):
        from star_gene_data import build_database_from_marker_matrix

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            marker_path = tmp_path / "markers.tsv"
            phenotype_path = tmp_path / "phenotypes.tsv"
            out_root = tmp_path / "db"
            marker_name = "chr1_27976568_27994848_deletion_18280"

            marker_path.write_text(
                f"SampleID\t{marker_name}\n"
                "S1\tNN\n"
                "S2\tAA\n",
                encoding="utf-8",
            )
            phenotype_path.write_text(
                "SampleID\tHKW\n"
                "S1\t22\n"
                "S2\t25\n",
                encoding="utf-8",
            )

            db_dir = build_database_from_marker_matrix(
                marker_matrix=marker_path,
                phenotype_table=phenotype_path,
                output_root=out_root,
                target_id="qHKW1",
                chrom="1",
                start=27976568,
                end=27994848,
                phenotype_columns=["HKW"],
                marker_columns=[marker_name],
                marker_positions={marker_name: 28892142},
            )

            variant_info = (db_dir / "variant_info.csv").read_text(encoding="utf-8")

            self.assertIn("-18280", variant_info)
            self.assertIn("True", variant_info)
            self.assertIn("SV", variant_info)

    def test_scan_maizego_sv_candidates_filters_by_window_and_length(self):
        from star_gene_data import scan_maizego_sv_candidates

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "SV.386014.hmp.txt"
            matrix_path.write_text(
                "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tS1\tS2\tS3\tS4\n"
                "chr1_30450000_30458900_insertion_8900\tA/T\t1\t30454500\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tAA\tNN\tTT\n"
                "chr1_30452000_30452100_insertion_100\tA/T\t1\t30452000\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tAA\tAA\tAA\n"
                "chr1_40000000_40008900_insertion_8900\tA/T\t1\t40004500\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tNN\tNN\tTT\n",
                encoding="utf-8",
            )

            candidates = scan_maizego_sv_candidates(
                matrix_paths=[matrix_path],
                chrom="1",
                window_start=30440000,
                window_end=30540000,
                length_min=8500,
                length_max=9500,
            )

            self.assertEqual(len(candidates), 1)
            candidate = candidates[0]
            self.assertEqual(candidate["marker"], "chr1_30450000_30458900_insertion_8900")
            self.assertEqual(candidate["variant_type"], "insertion")
            self.assertEqual(candidate["sv_length"], 8900)
            self.assertEqual(candidate["counts"], "AA:2;NN:1;TT:1")
            self.assertEqual(candidate["source_path"], str(matrix_path))

    def test_extract_maizego_marker_matrix_transposes_selected_record(self):
        from star_gene_data import extract_maizego_marker_matrix

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "SV.386014.hmp.txt"
            output_path = tmp_path / "qHKW1_marker.tsv"
            marker = "chr1_30450000_30458900_insertion_8900"
            matrix_path.write_text(
                "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tS1\tS2\tS3\n"
                f"{marker}\tA/T\t1\t30454500\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tNN\tTT\n",
                encoding="utf-8",
            )

            written = extract_maizego_marker_matrix(
                matrix_path=matrix_path,
                marker_id=marker,
                output_path=output_path,
            )

            self.assertEqual(written, output_path)
            matrix_text = output_path.read_text(encoding="utf-8")
            self.assertEqual(
                matrix_text,
                f"SampleID\t{marker}\nS1\tAA\nS2\tNN\nS3\tTT\n",
            )

    def test_complete_database_can_supply_coordinates_and_phenotypes(self):
        from star_gene_validation import StarGeneValidator

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_root = tmp_path / "star_gene_database"
            results_root = tmp_path / "star_gene_results"
            db_dir = db_root / "paper1" / "qHKW1"
            db_dir.mkdir(parents=True)

            (db_dir / "gene_info.json").write_text(
                '{"gene_id":"qHKW1","chrom":"1","start":100,"end":200}',
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\nAA,2,Hap1,AA\nNN,2,Hap2,NN\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\nS1,AA,Hap1\nS2,AA,Hap1\nS3,NN,Hap2\nS4,NN,Hap2\n",
                encoding="utf-8",
            )
            (db_dir / "phenotype_data.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name,HKW\nS1,AA,Hap1,25\nS2,AA,Hap1,26\nS3,NN,Hap2,22\nS4,NN,Hap2,23\n",
                encoding="utf-8",
            )

            paper = {
                "paper_id": "paper1",
                "short_name": "p1",
                "local_expected_paths": [
                    {
                        "key": "target_variant_matrix",
                        "path": "missing/markers.tsv",
                        "format": "unknown",
                        "required": True,
                    },
                    {
                        "key": "phenotype_table",
                        "path": "missing/phenotype.tsv",
                        "format": "phenotype_table",
                        "required": True,
                    },
                ],
            }
            target = {
                "target_id": "qHKW1",
                "gene_or_locus": "qHKW1",
                "requires_coordinate_resolution": True,
                "coordinates": None,
                "traits": [{"name": "hundred-kernel weight", "local_columns": ["HKW"]}],
            }

            validator = StarGeneValidator(
                manifest={"papers": [paper]},
                database_root=db_root,
                results_root=results_root,
                check_only=False,
            )
            check = validator.check_target(paper, target)
            coords = validator.resolve_target_coordinates(target, check["database_path"])

            self.assertEqual(check["status"], "ready_for_analysis")
            self.assertIn("database_complete", check["notes"])
            self.assertIn("database_coordinates:qHKW1=1:100-200", check["notes"])
            self.assertIn("HKW", check["phenotype_columns"])
            self.assertEqual(coords, ("1", 100, 200))

    def test_required_database_source_blocks_substitute_marker_database(self):
        from star_gene_validation import StarGeneValidator

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_root = tmp_path / "star_gene_database"
            results_root = tmp_path / "star_gene_results"
            db_dir = db_root / "paper1" / "qHKW1"
            db_dir.mkdir(parents=True)

            (db_dir / "gene_info.json").write_text(
                '{"gene_id":"qHKW1","chrom":"1","start":100,"end":200,"source":"marker_matrix"}',
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\nAA,2,Hap1,AA\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\nS1,AA,Hap1\n",
                encoding="utf-8",
            )
            (db_dir / "phenotype_data.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name,HKW\nS1,AA,Hap1,25\n",
                encoding="utf-8",
            )

            paper = {"paper_id": "paper1", "short_name": "p1", "local_expected_paths": []}
            target = {
                "target_id": "qHKW1",
                "gene_or_locus": "qHKW1",
                "requires_coordinate_resolution": True,
                "coordinates": None,
                "required_database_source": "maizego_paper_marker",
                "traits": [{"name": "hundred-kernel weight", "local_columns": ["HKW"]}],
            }

            validator = StarGeneValidator(
                manifest={"papers": [paper]},
                database_root=db_root,
                results_root=results_root,
                check_only=True,
            )
            check = validator.check_target(paper, target)

            self.assertEqual(check["status"], "unsupported_input_format_for_analysis")
            self.assertTrue(any(note.startswith("database_source_mismatch:") for note in check["notes"]))

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

    def test_maize_qhkw1_prepare_cli_reports_missing_full_sv_package(self):
        from prepare_maize2019_qhkw1_paper_genotype import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main(["--sv-root", str(tmp_path)])

            self.assertEqual(rc, 2)
            output = stdout.getvalue()
            self.assertIn("SV.386014.zip", output)
            self.assertIn("paper genotype package is missing", output)

    def test_maize_qhkw1_prepare_ignores_non_matrix_text_files(self):
        from prepare_maize2019_qhkw1_paper_genotype import find_candidate_matrices

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            (tmp_path / "paper.txt").write_text("This is article text, not genotype data.\n", encoding="utf-8")
            matrix_path = tmp_path / "MS.step1.org.txt"
            matrix_path.write_text(
                "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tS1\n"
                "chr1_30450000_30458900_insertion_8900\tA/T\t1\t30454500\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\n",
                encoding="utf-8",
            )

            self.assertEqual(find_candidate_matrices(tmp_path), [matrix_path])

    def test_maize_qhkw1_prepare_builds_paper_marker_database(self):
        from prepare_maize2019_qhkw1_paper_genotype import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            sv_root = tmp_path / "maizego"
            sv_root.mkdir()
            marker = "chr1_30450000_30458900_insertion_8900"
            matrix_path = sv_root / "SV.386014.hmp.txt"
            matrix_path.write_text(
                "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tS1\tS2\tS3\n"
                f"{marker}\tA/T\t1\t30454500\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tNN\tTT\n",
                encoding="utf-8",
            )
            phenotype_path = tmp_path / "phenotypes.csv"
            phenotype_path.write_text(
                "<Trait>,100grainweight\n"
                "S1,25.1\n"
                "S2,29.3\n"
                "S3,-999\n",
                encoding="utf-8",
            )
            output_root = tmp_path / "db"
            candidate_output = tmp_path / "candidates.tsv"
            marker_output = tmp_path / "marker.tsv"

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main([
                    "--sv-root", str(sv_root),
                    "--phenotype-table", str(phenotype_path),
                    "--output-root", str(output_root),
                    "--candidate-output", str(candidate_output),
                    "--marker-output", str(marker_output),
                ])

            self.assertEqual(rc, 0)
            gene_info = json.loads((output_root / "qHKW1" / "gene_info.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_info["source"], "maizego_paper_marker")
            self.assertEqual(gene_info["source_marker"], marker)
            self.assertTrue(marker_output.exists())
            self.assertTrue(candidate_output.exists())


if __name__ == "__main__":
    unittest.main()
