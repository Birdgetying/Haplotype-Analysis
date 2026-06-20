#!/usr/bin/env python3
"""Tests for lightweight star-gene data source metadata."""

from pathlib import Path
from contextlib import redirect_stdout
from io import StringIO
import csv
import gzip
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
        from star_gene_data import iter_data_files

        wheat_files = list(iter_data_files(paper="wheat2024"))
        self.assertTrue(wheat_files)
        portal_files = [f for f in wheat_files if f.key in {"wheat2024_wwwg2b_portal", "wheat2024_earlham_watseq"}]
        self.assertEqual(len(portal_files), 2)
        self.assertTrue(all(f.default_action == "instruction_only" for f in portal_files))
        self.assertTrue(any("wwwg2b.com" in f.source for f in wheat_files))

    def test_wheat_q7b_small_files_have_wwwg2b_download_commands(self):
        from star_gene_data import build_download_commands, iter_data_files

        wheat_files = list(iter_data_files(paper="wheat2024"))
        keys = {f.key for f in wheat_files}

        self.assertIn("wheat2024_q7b_ph_figure3g", keys)
        self.assertIn("wheat2024_watkins_jic_phenotypes", keys)

        commands = build_download_commands(paper="wheat2024")
        command_text = "\n".join(commands)
        self.assertIn("get_download_url_form_onedrive", command_text)
        self.assertIn("Watseq_Figure_3g_NIL_Q7B-PH_field_data.xlsx", command_text)
        self.assertIn("external_data/wheat_nature_2024/wwwg2b/q7b_ph", command_text)
        self.assertIn("11032_2014_34_MOESM1_ESM.doc", command_text)
        self.assertIn("static-content.springer.com", command_text)
        self.assertIn("curl.exe -L -C -", command_text)

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

    def test_cli_accepts_robust_discovery_score_mode(self):
        from star_gene_validation import build_arg_parser

        args = build_arg_parser().parse_args([
            "--check-only",
            "--paper",
            "wheat2024",
            "--score-mode",
            "robust_discovery",
        ])

        self.assertEqual(args.score_mode, "robust_discovery")

    def test_direction_aware_top_haplotype_uses_expected_decrease(self):
        from star_gene_validation import StarGeneValidator

        runner = StarGeneValidator(
            manifest={"papers": []},
            database_root=Path("star_gene_database"),
            results_root=Path("star_gene_results"),
            check_only=True,
        )
        score_results = {
            "per_haplotype": {
                "TallHap": {"total": 2.0, "mean_phenotype": 100.0},
                "ShortHap": {"total": 1.0, "mean_phenotype": 80.0},
                "TinyShortHap": {"total": 3.0, "mean_phenotype": 70.0},
            },
            "per_sample": [
                {"sample_id": "S1", "haplotype": "TallHap", "score": 2.0, "phenotype": 100.0},
                {"sample_id": "S2", "haplotype": "TallHap", "score": 2.0, "phenotype": 102.0},
                {"sample_id": "S3", "haplotype": "ShortHap", "score": 1.0, "phenotype": 80.0},
                {"sample_id": "S4", "haplotype": "ShortHap", "score": 1.0, "phenotype": 78.0},
                {"sample_id": "S5", "haplotype": "TinyShortHap", "score": 3.0, "phenotype": 70.0},
            ],
        }

        self.assertEqual(runner.pick_directional_top_haplotype(score_results, "increases_trait")[0], "TallHap")
        self.assertEqual(runner.pick_directional_top_haplotype(score_results, "decreases_trait")[0], "ShortHap")

    def test_score_tooltips_use_viewport_coordinates(self):
        source = Path("haplotype_phenotype_analysis.py").read_text(encoding="utf-8")

        integrated_start = source.index("function drawHaplotypeScorePlot(scoreData)")
        integrated_end = source.index("function updateScoreLegend", integrated_start)
        integrated_block = source[integrated_start:integrated_end]

        standalone_start = source.index("def generate_haplotype_score_html")
        standalone_end = source.index("def generate_effect_boxplot_html", standalone_start)
        standalone_block = source[standalone_start:standalone_end]

        for block in (integrated_block, standalone_block):
            self.assertIn("clientX", block)
            self.assertIn("clientY", block)
            self.assertNotIn("pageX", block)
            self.assertNotIn("pageY", block)

    def test_haplotype_score_html_embeds_default_and_robust_modes(self):
        from haplotype_phenotype_analysis import ReportGenerator

        def score_bundle(mode, score):
            return {
                "score_mode": mode,
                "per_sample": [
                    {"sample_id": "S1", "haplotype": "Hap1", "score": score, "phenotype": 10.0},
                    {"sample_id": "S2", "haplotype": "Hap2", "score": score + 1.0, "phenotype": 12.0},
                ],
                "per_haplotype": {
                    "Hap1": {"total": score, "sample_reliability": 0.9},
                    "Hap2": {"total": score + 1.0, "sample_reliability": 0.8},
                },
                "component_weights": {"variant_effect": 1.0},
                "r_squared": 0.5,
                "regression_pvalue": 0.01,
                "slope": 1.0,
                "intercept": 9.0,
            }

        with tempfile.TemporaryDirectory() as tmp:
            robust_dir = Path(tmp) / "GeneA__robust_discovery"
            default_dir = Path(tmp) / "GeneA"
            robust_dir.mkdir()
            default_dir.mkdir()
            (default_dir / "haplotype_scores.json").write_text(
                json.dumps({"Trait": score_bundle("default", 1.0)}),
                encoding="utf-8",
            )

            generator = ReportGenerator(output_dir=str(robust_dir))
            generator._cached_all_haplotype_scores = {"Trait": score_bundle("robust_discovery", 2.0)}
            generator._cached_haplotype_score = generator._cached_all_haplotype_scores["Trait"]

            html_path = generator.generate_haplotype_score_html(gene_id="GeneA", phenotype_col="Trait")
            html = Path(html_path).read_text(encoding="utf-8")

            self.assertIn("var allScoreModeData", html)
            self.assertIn('"default"', html)
            self.assertIn('"robust_discovery"', html)
            self.assertIn('"current_mode": "robust_discovery"', html)
            self.assertIn('"score": 2.0', html)
            self.assertIn('"score": 1.0', html)
            self.assertIn("switchScoreMode('default')", html)
            self.assertIn("switchScoreMode('robust_discovery')", html)
            self.assertIn("mode-toggle-btn", html)

    def test_integrated_html_has_score_mode_toggle_wiring(self):
        source = Path("haplotype_phenotype_analysis.py").read_text(encoding="utf-8")

        integrated_start = source.index("def generate_integrated_html")
        integrated_end = source.index("def generate_haplotype_network_html", integrated_start)
        integrated_block = source[integrated_start:integrated_end]

        self.assertIn("var allScoreModeData", integrated_block)
        self.assertIn("function switchScoreMode", integrated_block)
        self.assertIn("switchScoreMode('default')", integrated_block)
        self.assertIn("switchScoreMode('robust_discovery')", integrated_block)
        self.assertIn("mode-toggle-btn", integrated_block)

    def test_integrated_html_has_local_post_gwas_evidence_panel(self):
        source = Path("haplotype_phenotype_analysis.py").read_text(encoding="utf-8")

        integrated_start = source.index("def generate_integrated_html")
        integrated_end = source.index("def generate_haplotype_network_html", integrated_start)
        integrated_block = source[integrated_start:integrated_end]

        self.assertIn("def _summarize_post_gwas_evidence", source)
        self.assertIn("post-gwas-evidence", integrated_block)
        self.assertIn("Post-GWAS Evidence", integrated_block)
        self.assertIn("Candidate rule", integrated_block)
        self.assertIn("Top local marker", integrated_block)
        self.assertIn("postGwasEvidence", integrated_block)

    def test_summarize_post_gwas_evidence_uses_local_variant_context(self):
        from haplotype_phenotype_analysis import _summarize_post_gwas_evidence

        evidence = _summarize_post_gwas_evidence(
            variant_positions=[90, 110, 206],
            variant_pvalues={90: 1e-3, 110: 1e-6, 206: 0.2},
            variant_info={
                90: {"ref": "A", "alt": "G", "maf": 0.20, "missing_rate": 0.01, "annotation": "promoter"},
                110: {"ref": "C", "alt": "T", "maf": 0.35, "missing_rate": 0.02, "annotation": "missense"},
                206: {"ref": "A", "alt": "<DEL>", "is_sv": True, "maf": 0.08, "missing_rate": 0.05},
            },
            snp_effects={90: "promoter", 110: "missense", 206: "SV"},
            gene_start=100,
            gene_end=200,
            promoter_start=50,
            promoter_end=99,
            flank_bp=5,
            gwas_threshold=1e-5,
        )

        self.assertEqual(evidence["top_marker"]["pos"], 110)
        self.assertEqual(evidence["top_marker"]["candidate_status"], "strict_local_gwas_window")
        self.assertEqual(evidence["candidate_rule"]["status"], "strict_local_gwas_window")
        self.assertEqual(evidence["candidate_rule"]["significant_in_window"], 1)
        self.assertEqual(evidence["variant_type_counts"]["SNP"], 2)
        self.assertEqual(evidence["variant_type_counts"]["SV"], 1)
        self.assertEqual(evidence["position_counts"]["promoter"], 1)
        self.assertEqual(evidence["position_counts"]["body"], 1)
        self.assertEqual(evidence["position_counts"]["outside"], 1)

    def test_score_mode_collection_prefers_score_json_over_stale_report_mode(self):
        from haplotype_phenotype_analysis import ReportGenerator

        with tempfile.TemporaryDirectory() as tmp:
            generator = ReportGenerator(output_dir=str(Path(tmp) / "GeneA__robust_discovery"))
            generator.score_mode = "default"

            collected = generator._collect_score_mode_data({
                "Trait": {
                    "score_mode": "robust_discovery",
                    "per_sample": [],
                    "per_haplotype": {},
                }
            })

            self.assertEqual(collected["current_mode"], "robust_discovery")
            self.assertIn("robust_discovery", collected["modes"])
            self.assertNotIn("default", collected["modes"])

    def test_robust_discovery_penalizes_tiny_ambiguous_haplotypes(self):
        import pandas as pd
        from haplotype_phenotype_analysis import HaplotypeScorer

        rows = []
        rows.extend(
            {
                "SampleID": f"H1_{i}",
                "Hap_Name": "Hap1",
                "Haplotype_Seq": "A|G|C",
                "TGW_mean": 41.0 + (0.1 if i % 2 else -0.1),
            }
            for i in range(40)
        )
        rows.extend(
            {
                "SampleID": f"H2_{i}",
                "Hap_Name": "Hap2",
                "Haplotype_Seq": "C|G|C",
                "TGW_mean": 39.0 + (0.1 if i % 2 else -0.1),
            }
            for i in range(35)
        )
        rows.extend(
            {
                "SampleID": f"H5_{i}",
                "Hap_Name": "Hap5",
                "Haplotype_Seq": "C/A|G|C",
                "TGW_mean": 40.6 + (0.1 if i % 2 else -0.1),
            }
            for i in range(3)
        )
        hap_sample_df = pd.DataFrame(rows)
        variant_positions = [291759689, 291760677, 291761315]
        variant_info = {
            291759689: {"ref": "C", "alt": "A", "maf": 0.48, "annotation": "promoter"},
            291760677: {"ref": "G", "alt": "A", "maf": 0.01, "annotation": "promoter"},
            291761315: {"ref": "T", "alt": "C", "maf": 0.47, "annotation": "promoter"},
        }
        gwas_data = [
            {"pos": 291759689, "pvalue": 1e-8},
            {"pos": 291760677, "pvalue": 1e-4},
            {"pos": 291761315, "pvalue": 1e-6},
        ]
        effect_results = {
            "haplotype_effects": [
                {"haplotype": "Hap1", "cohens_d": 0.9, "n_samples": 40},
                {"haplotype": "Hap2", "cohens_d": 0.4, "n_samples": 35},
                {"haplotype": "Hap5", "cohens_d": 2.5, "n_samples": 3},
            ]
        }

        robust = HaplotypeScorer(
            hap_sample_df,
            variant_positions,
            variant_info=variant_info,
            snp_effects={pos: "promoter" for pos in variant_positions},
            gwas_data=gwas_data,
            effect_results=effect_results,
            phenotype_col="TGW_mean",
            score_mode="robust_discovery",
        ).score_all()
        default = HaplotypeScorer(
            hap_sample_df,
            variant_positions,
            variant_info=variant_info,
            snp_effects={pos: "promoter" for pos in variant_positions},
            gwas_data=gwas_data,
            effect_results=effect_results,
            phenotype_col="TGW_mean",
        ).score_all()

        self.assertEqual(robust["score_mode"], "robust_discovery")
        self.assertNotIn("sample_reliability", default["per_haplotype"]["Hap5"])
        self.assertLess(robust["per_haplotype"]["Hap5"]["sample_reliability"], 0.2)
        self.assertLess(robust["per_haplotype"]["Hap5"]["ambiguity_factor"], 1.0)
        self.assertGreater(
            robust["per_haplotype"]["Hap1"]["total"],
            robust["per_haplotype"]["Hap5"]["total"],
        )

    def test_pca_plot_handles_two_variant_haplotypes(self):
        import pandas as pd
        from haplotype_phenotype_analysis import ReportGenerator

        with tempfile.TemporaryDirectory() as tmp:
            hap_sample_df = pd.DataFrame({
                "SampleID": ["S1", "S2", "S3", "S4"],
                "Hap_Name": ["Hap1", "Hap1", "Hap2", "Hap3"],
                "Haplotype_Seq": ["T|G", "T|G", "T|A", "C|G"],
                "TGW_mean": [39.7, 40.1, 31.4, 39.6],
            })
            generator = ReportGenerator(output_dir=tmp)
            stdout = StringIO()

            with redirect_stdout(stdout):
                html_path = generator.generate_pca_plot_html(
                    hap_sample_df,
                    phenotype_col="TGW_mean",
                )

            self.assertTrue(Path(html_path).exists())
            html = Path(html_path).read_text(encoding="utf-8")
            self.assertIn("Principal Component Analysis", html)
            self.assertIn("PC2", html)
            self.assertIn("PCA分析图已保存", stdout.getvalue())

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

    def test_build_database_from_marker_matrix_reads_long_tsv_with_truncated_sniffer_sample(self):
        from star_gene_data import build_database_from_marker_matrix

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            marker_path = tmp_path / "markers.tsv"
            phenotype_path = tmp_path / "phenotypes.tsv"
            out_root = tmp_path / "db"
            marker_cols = [f"chr6A_{237732717 + i}" for i in range(190)]

            with marker_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                writer.writerow(["SampleID", *marker_cols])
                writer.writerow(["S1", *(["A"] * len(marker_cols))])
                writer.writerow(["S2", *(["T"] * len(marker_cols))])
                writer.writerow(["S3", *(["A"] * len(marker_cols))])
                writer.writerow(["S4", *(["T"] * len(marker_cols))])
            with phenotype_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                writer.writerow(["SampleID", "Trait"])
                writer.writerow(["S1", 1.0])
                writer.writerow(["S2", 2.0])
                writer.writerow(["S3", 1.5])
                writer.writerow(["S4", 2.5])

            db_dir = build_database_from_marker_matrix(
                marker_matrix=marker_path,
                phenotype_table=phenotype_path,
                output_root=out_root,
                target_id="LongHeaderGene",
                chrom="1",
                start=1000,
                end=2000,
                phenotype_columns=["Trait"],
                min_haplotype_count=1,
            )

            self.assertTrue((db_dir / "haplotype_samples.csv").exists())
            with (db_dir / "variant_info.csv").open() as f:
                variant_rows = list(csv.DictReader(f))
            self.assertEqual(len(variant_rows), len(marker_cols))

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

    def test_scan_maizego_sv_catalogue_candidates_reports_records_without_samples(self):
        from star_gene_data import scan_maizego_sv_catalogue_candidates

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            catalogue_path = tmp_path / "svs.final.ms.txt"
            catalogue_path.write_text(
                "chr1\t30450000\t30458900\tinsertion\t8900\t+\tMo17.chr1.64\t29961002\t29969902\t295670\t181512\t324959\t0.558569\tACGT\n"
                "chr1\t30452000\t30452100\tinsertion\t100\t+\tMo17.chr1.64\t29961002\t29961102\t295670\t181512\t324959\t0.558569\tAC\n"
                "chr1\t40000000\t40008900\tdeletion\t8900\t+\tMo17.chr1.64\t39961002\t39961002\t295670\t181512\t324959\t0.558569\tACGT\n",
                encoding="utf-8",
            )

            candidates = scan_maizego_sv_catalogue_candidates(
                catalogue_paths=[catalogue_path],
                chrom="1",
                window_start=30440000,
                window_end=30540000,
                length_min=8500,
                length_max=9500,
            )

            self.assertEqual(len(candidates), 1)
            candidate = candidates[0]
            self.assertEqual(candidate["chrom"], "1")
            self.assertEqual(candidate["start"], 30450000)
            self.assertEqual(candidate["end"], 30458900)
            self.assertEqual(candidate["variant_type"], "insertion")
            self.assertEqual(candidate["sv_length"], 8900)
            self.assertEqual(candidate["sample_genotype_columns"], 0)
            self.assertEqual(candidate["has_sample_genotypes"], "False")
            self.assertEqual(candidate["source_path"], str(catalogue_path))

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

    def test_literature_audit_matches_top_haplotype_for_single_variant(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,A,G,0,False,0.5,0.0,promoter,m1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "G,2,Hap1,G\n"
                "A,1,Hap2,A\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,G,Hap1\n"
                "S2,G,Hap1\n"
                "S3,A,Hap2\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({
                    "Trait": {
                        "per_haplotype": {
                            "Hap1": {"total": 4.2},
                            "Hap2": {"total": 1.1},
                        },
                        "per_sample": [
                            {"sample_id": "S1", "haplotype": "Hap1"},
                            {"sample_id": "S2", "haplotype": "Hap1"},
                            {"sample_id": "S3", "haplotype": "Hap2"},
                        ],
                    }
                }),
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_variants": [
                        {
                            "variant_name": "KnownSNP",
                            "chrom": "chr1",
                            "position": 100,
                            "marker_id": "m1",
                            "expected_allele": "G",
                            "variant_class": "promoter_snp",
                            "source_note": "unit test",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            self.assertEqual(len(records), 1)
            row = records[0]
            self.assertEqual(row["record_type"], "variant")
            self.assertEqual(row["variant_name"], "KnownSNP")
            self.assertTrue(row["covered_in_database"])
            self.assertTrue(row["segregating_in_current_samples"])
            self.assertEqual(row["carrier_count"], 2)
            self.assertEqual(row["top_scored_haplotype"], "Hap1")
            self.assertEqual(row["top_haplotype_sample_count"], 2)
            self.assertEqual(row["top_haplotype_carrier_count"], 2)
            self.assertTrue(row["top_haplotype_contains_expected"])
            self.assertTrue(row["top_haplotype_exact_expected"])
            self.assertEqual(row["validation_status"], "matched_top_haplotype")
            self.assertTrue((result_dir / "literature_variant_audit.csv").exists())
            self.assertTrue((result_dir / "literature_variant_audit.json").exists())

    def test_literature_audit_reports_directional_top_match(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "1,D1a,D1b,0,False,0.5,0.0,functional_marker,Rht-D1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "D1a,20,HapTall,D1a\n"
                "D1b,30,HapShort,D1b\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,D1a,HapTall\n"
                "S2,D1b,HapShort\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({"Trait": {"per_haplotype": {"HapTall": {"total": 5.0}, "HapShort": {"total": 1.0}}}}),
                encoding="utf-8",
            )
            (result_dir.parent.parent / "validation_summary.csv").parent.mkdir(parents=True, exist_ok=True)
            (result_dir.parent.parent / "validation_summary.csv").write_text(
                "target_id,result_path,directional_top_haplotype,directional_top_haplotype_score,directional_top_haplotype_sample_count\n"
                f"RhtTest,{result_dir},HapShort,1.0,30\n",
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "RhtTest",
                    "literature_variants": [
                        {
                            "variant_name": "Rht-D1b",
                            "marker_id": "Rht-D1",
                            "expected_allele": "D1b",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            row = records[0]
            self.assertEqual(row["validation_status"], "present_but_not_top")
            self.assertEqual(row["directional_top_haplotype"], "HapShort")
            self.assertTrue(row["directional_top_haplotype_exact_expected"])
            self.assertEqual(row["directional_validation_status"], "matched_directional_top_haplotype")

    def test_literature_audit_distinguishes_heterozygous_top_variant_match(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,C,A,0,False,0.5,0.0,promoter,m1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "A,3,Hap1,A\n"
                "C/A,2,Hap2,C/A\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,A,Hap1\n"
                "S2,A,Hap1\n"
                "S3,A,Hap1\n"
                "S4,C/A,Hap2\n"
                "S5,C/A,Hap2\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({
                    "Trait": {
                        "per_haplotype": {
                            "Hap1": {"total": 2.0},
                            "Hap2": {"total": 5.0},
                        }
                    }
                }),
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_variants": [
                        {
                            "variant_name": "KnownSNP",
                            "marker_id": "m1",
                            "expected_allele": "A",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            row = records[0]
            self.assertTrue(row["top_haplotype_contains_expected"])
            self.assertFalse(row["top_haplotype_exact_expected"])
            self.assertEqual(row["exact_matching_haplotypes"], "Hap1:3")
            self.assertEqual(row["validation_status"], "contained_in_top_haplotype_not_exact")

    def test_literature_audit_reports_missing_and_nonsegregating_variants(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,A,G,0,False,0.0,0.0,promoter,m1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\nA,3,Hap1,A\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\nS1,A,Hap1\nS2,A,Hap1\nS3,A,Hap1\n",
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_variants": [
                        {
                            "variant_name": "ExpectedMissingAllele",
                            "chrom": "chr1",
                            "position": 100,
                            "marker_id": "m1",
                            "expected_allele": "G",
                        },
                        {
                            "variant_name": "AbsentMarker",
                            "chrom": "chr1",
                            "position": 200,
                            "marker_id": "m2",
                            "expected_allele": "T",
                        },
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            by_name = {row["variant_name"]: row for row in records}
            self.assertEqual(by_name["ExpectedMissingAllele"]["validation_status"], "present_not_segregating")
            self.assertTrue(by_name["ExpectedMissingAllele"]["covered_in_database"])
            self.assertFalse(by_name["ExpectedMissingAllele"]["segregating_in_current_samples"])
            self.assertEqual(by_name["ExpectedMissingAllele"]["carrier_count"], 0)
            self.assertEqual(by_name["ExpectedMissingAllele"]["allele_counts"], "A:3")
            self.assertEqual(by_name["AbsentMarker"]["validation_status"], "missing_from_database")
            self.assertFalse(by_name["AbsentMarker"]["covered_in_database"])

    def test_literature_audit_does_not_validate_monomorphic_expected_allele(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,G,A,0,False,0.0,0.0,promoter,m1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\nG,3,Hap1,G\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\nS1,G,Hap1\nS2,G,Hap1\nS3,G,Hap1\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({"Trait": {"per_haplotype": {"Hap1": {"total": 5.0}}}}),
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_variants": [
                        {
                            "variant_name": "MonomorphicExpected",
                            "marker_id": "m1",
                            "expected_allele": "G",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            row = records[0]
            self.assertEqual(row["carrier_count"], 3)
            self.assertFalse(row["segregating_in_current_samples"])
            self.assertTrue(row["top_haplotype_contains_expected"])
            self.assertEqual(row["validation_status"], "present_not_segregating")

    def test_literature_audit_reports_group_haplotype_status(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,A,G,0,False,0.5,0.0,promoter,m1\n"
                "200,C,T,0,False,0.5,0.0,promoter,m2\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "G|T,2,Hap1,G|T\n"
                "G|C,1,Hap2,G|C\n"
                "A|T,1,Hap3,A|T\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,G|T,Hap1\n"
                "S2,G|T,Hap1\n"
                "S3,G|C,Hap2\n"
                "S4,A|T,Hap3\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({
                    "Trait": {
                        "per_haplotype": {
                            "Hap1": {"total": 2.0},
                            "Hap2": {"total": 5.0},
                            "Hap3": {"total": 1.0},
                        }
                    }
                }),
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_haplotypes": [
                        {
                            "haplotype_name": "PublishedPair",
                            "expected_markers": {"m1": "G", "m2": "T"},
                            "source_note": "unit test pair",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            self.assertEqual(len(records), 1)
            row = records[0]
            self.assertEqual(row["record_type"], "haplotype")
            self.assertEqual(row["variant_name"], "PublishedPair")
            self.assertTrue(row["covered_in_database"])
            self.assertTrue(row["segregating_in_current_samples"])
            self.assertEqual(row["carrier_count"], 2)
            self.assertEqual(row["top_scored_haplotype"], "Hap2")
            self.assertFalse(row["top_haplotype_contains_expected"])
            self.assertEqual(row["validation_status"], "present_but_not_top")

    def test_literature_audit_requires_exact_top_haplotype_for_group_match(self):
        from star_gene_literature_audit import run_literature_audit

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db"
            result_dir = tmp_path / "result"
            db_dir.mkdir()
            result_dir.mkdir()
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,len_diff,is_sv,maf,missing_rate,annotation,marker_id\n"
                "100,C,A,0,False,0.5,0.0,promoter,m1\n"
                "200,G,A,0,False,0.5,0.0,promoter,m2\n"
                "300,T,C,0,False,0.5,0.0,promoter,m3\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "A|G|C,3,Hap1,A|G|C\n"
                "C/A|G|C,2,Hap2,C/A|G|C\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,A|G|C,Hap1\n"
                "S2,A|G|C,Hap1\n"
                "S3,A|G|C,Hap1\n"
                "S4,C/A|G|C,Hap2\n"
                "S5,C/A|G|C,Hap2\n",
                encoding="utf-8",
            )
            (result_dir / "haplotype_scores.json").write_text(
                json.dumps({
                    "Trait": {
                        "per_haplotype": {
                            "Hap1": {"total": 2.0},
                            "Hap2": {"total": 5.0},
                        }
                    }
                }),
                encoding="utf-8",
            )

            records = run_literature_audit(
                paper={"paper_id": "paper1"},
                target={
                    "target_id": "Gene1",
                    "literature_haplotypes": [
                        {
                            "haplotype_name": "PublishedAGC",
                            "expected_markers": {"m1": "A", "m2": "G", "m3": "C"},
                            "source_note": "unit test strict haplotype",
                        }
                    ],
                },
                database_path=db_dir,
                result_path=result_dir,
            )

            row = records[0]
            self.assertTrue(row["top_haplotype_contains_expected"])
            self.assertFalse(row["top_haplotype_exact_expected"])
            self.assertEqual(row["exact_matching_haplotypes"], "Hap1:3")
            self.assertEqual(row["validation_status"], "contained_in_top_haplotype_not_exact")

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

    def test_maize_qhkw1_prepare_reports_catalogue_candidates_but_keeps_blocked(self):
        from prepare_maize2019_qhkw1_paper_genotype import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            sv_root = tmp_path / "maizego"
            extracted_dir = sv_root / "SV.386014" / "SV.386014"
            extracted_dir.mkdir(parents=True)
            (extracted_dir / "svs.final.ms.txt").write_text(
                "chr1\t30450000\t30458900\tinsertion\t8900\t+\tMo17.chr1.64\t29961002\t29969902\t295670\t181512\t324959\t0.558569\tACGT\n",
                encoding="utf-8",
            )
            matrix_candidate_output = tmp_path / "matrix_candidates.tsv"
            catalogue_output = tmp_path / "catalogue_candidates.tsv"

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main([
                    "--sv-root", str(sv_root),
                    "--candidate-output", str(matrix_candidate_output),
                    "--catalogue-output", str(catalogue_output),
                ])

            self.assertEqual(rc, 3)
            output = stdout.getvalue()
            self.assertIn("Candidate matrix files: 0", output)
            self.assertIn("SV catalogue candidate records: 1", output)
            self.assertIn("not sample-level genotype matrices", output)
            self.assertTrue(matrix_candidate_output.exists())
            catalogue_text = catalogue_output.read_text(encoding="utf-8")
            self.assertIn("sample_genotype_columns", catalogue_text)
            self.assertIn("chr1", catalogue_text)

    def test_maize_qhkw1_prepare_does_not_reextract_existing_catalogues(self):
        from prepare_maize2019_qhkw1_paper_genotype import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            sv_root = tmp_path / "maizego"
            extracted_dir = sv_root / "SV.386014" / "SV.386014"
            extracted_dir.mkdir(parents=True)
            (extracted_dir / "svs.final.ms.txt").write_text(
                "chr1\t30450000\t30458900\tinsertion\t8900\t+\tMo17.chr1.64\t29961002\t29969902\t295670\t181512\t324959\t0.558569\tACGT\n",
                encoding="utf-8",
            )
            invalid_zip = sv_root / "SV.386014.zip"
            invalid_zip.write_text("not a zip file\n", encoding="utf-8")

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main([
                    "--sv-root", str(sv_root),
                    "--sv-zip", str(invalid_zip),
                    "--candidate-output", str(tmp_path / "matrix_candidates.tsv"),
                    "--catalogue-output", str(tmp_path / "catalogue_candidates.tsv"),
                ])

            self.assertEqual(rc, 3)
            self.assertIn("SV catalogue candidate records: 1", stdout.getvalue())

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

    def test_wheat_q7b_ph_prepare_builds_figure3g_database(self):
        from openpyxl import Workbook
        from prepare_wheat2024_q7b_ph_figure3g import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_path = tmp_path / "figure3g.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "Figure_3_NIL_Q7B-PH_W141_field_"
            ws.append([
                "site", "year", "siteyear", "Accession", "allele",
                "Replicate", "PH_M_cm", "GY_Calc_tha", "GW_M_g1000grn",
            ])
            ws.append(["JIC", "H15", "JIC_H15", "WL0019", "W", 1, 103, 8.34, 48.55])
            ws.append(["JIC", "H15", "JIC_H15", "WL0019", "W", 2, 101, 8.10, 49.00])
            ws.append(["JIC", "H15", "JIC_H15", "WL0025", "P", 1, 97, 8.52, 50.19])
            ws.append(["JIC", "H15", "JIC_H15", "Paragon", "P", 1, 74, 7.90, 46.00])
            ws.append(["JIC", "H15", "JIC_H15", "W10074", "Par", 1, 98, 8.10, None])
            wb.save(source_path)

            output_root = tmp_path / "db"
            marker_output = tmp_path / "q7b_marker.tsv"
            phenotype_output = tmp_path / "q7b_pheno.tsv"

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main([
                    "--source-xlsx", str(source_path),
                    "--output-root", str(output_root),
                    "--marker-output", str(marker_output),
                    "--phenotype-output", str(phenotype_output),
                ])

            self.assertEqual(rc, 0)
            db_dir = output_root / "Q7B-HT"
            gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_info["source"], "wwwg2b_q7b_ph_figure3g")
            self.assertEqual(gene_info["source_file"], str(source_path))
            self.assertEqual(gene_info["source_marker"], "Q7B-PH_allele")
            self.assertLess(gene_info["start"], gene_info["end"])
            self.assertTrue(marker_output.exists())
            self.assertTrue(phenotype_output.exists())
            hap_samples = (db_dir / "haplotype_samples.csv").read_text(encoding="utf-8")
            phenotype_data = (db_dir / "phenotype_data.csv").read_text(encoding="utf-8")
            self.assertIn("Q7B-PH_allele", (db_dir / "variant_info.csv").read_text(encoding="utf-8"))
            self.assertIn("WL0019__plot0001", hap_samples)
            self.assertIn("WL0019__plot0002", hap_samples)
            self.assertIn("WL0019__plot0001", phenotype_data)

    def test_wheat_tagw2_ad_prepare_builds_vcf_region_databases(self):
        from openpyxl import Workbook
        from prepare_wheat2024_tagw2_ad import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            vcf_a = tmp_path / "chr6A.vcf.gz"
            vcf_d = tmp_path / "chr6D.vcf.gz"

            vcf_header = (
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                "WATDE0001\tWATDE0002\tWATDE0003\tModernA\n"
            )
            with gzip.open(vcf_a, "wt", encoding="utf-8", newline="") as f:
                f.write(vcf_header)
                f.write("6A\t237732100\toutside\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n")
                f.write("6A\t237734700\tvarA1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1\t0/0\n")
                f.write("6A\t237735000\tvarA2\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t0/0\t0/0\n")
                f.write("6A\t237761000\tafter\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n")
            with gzip.open(vcf_d, "wt", encoding="utf-8", newline="") as f:
                f.write(vcf_header.replace("6A", "6D"))
                f.write("6D\t175712300\tvarD1\tG\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1\t0/0\n")
                f.write("6D\t175721000\tvarD2\tT\tC\t.\tPASS\t.\tGT\t1/1\t1/1\t0/0\t0/0\n")

            phenotype_xlsx = tmp_path / "watkins.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "WISP_Watkins_JIC_CFLN10"
            ws.append(["StoreCode", "GW_M_g1000grn-CFLN10"])
            ws.append(["WATDE0001", 40.0])
            ws.append(["WATDE0002", 44.0])
            ws.append(["WATDE0003", 42.0])
            ws2 = wb.create_sheet("DFW_Watkins_JI_CFLN14")
            ws2.append(["StoreCode", "GW_M_g1000grn-CFLN14"])
            ws2.append(["WATDE0001", 41.0])
            ws2.append(["WATDE0002", 45.0])
            ws2.append(["WATDE0003", None])
            wb.save(phenotype_xlsx)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "external"

            stdout = StringIO()
            with redirect_stdout(stdout):
                rc = prepare_main([
                    "--vcf-a", str(vcf_a),
                    "--vcf-d", str(vcf_d),
                    "--phenotype-xlsx", str(phenotype_xlsx),
                    "--output-root", str(output_root),
                    "--intermediate-root", str(intermediate_root),
                    "--min-haplotype-count", "1",
                ])

            self.assertEqual(rc, 0)
            for target_id in ["TaGW2-A1", "TaGW2-D1"]:
                db_dir = output_root / target_id
                self.assertTrue((db_dir / "gene_info.json").exists())
                gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
                self.assertEqual(gene_info["source"], "watseq_tagw2_ad_vcf")
                self.assertIn("TGW_mean", (db_dir / "phenotype_data.csv").read_text(encoding="utf-8"))
                self.assertIn("marker_id", (db_dir / "variant_info.csv").read_text(encoding="utf-8"))
                self.assertIn("WATDE0001", (db_dir / "haplotype_samples.csv").read_text(encoding="utf-8"))

    def test_wheat_rht1_prepare_builds_functional_snp_databases(self):
        import pandas as pd
        from openpyxl import Workbook
        from prepare_wheat2024_rht1_functional_snps import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            vcf_path = tmp_path / "rht1.vcf"
            vcf_path.write_text(
                "\n".join([
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWATDE0001\tWATDE0002\tWATDE0003\tWATDE0004",
                    "chr4B\t30861571\tchr4B_30861571\tC\tT\t.\t.\tANN=T|stop_gained|HIGH|TraesCS4B01G043100|TraesCS4B01G043100|transcript|TraesCS4B01G043100.1|protein_coding|1/1|c.190C>T|p.Gln64*|190/1866|190/1866|64/621||\tGT\t0/0\t1/1\t1/1\t0/0",
                    "chr4D\t18781242\tchr4D_18781242\tG\tT\t.\t.\tANN=T|stop_gained|HIGH|TraesCS4D01G040400|TraesCS4D01G040400|transcript|TraesCS4D01G040400.1|protein_coding|1/1|c.181G>T|p.Glu61*|181/1872|181/1872|61/623||\tGT\t0/0\t1/1\t0/0\t1/1",
                    "",
                ]),
                encoding="utf-8",
            )

            phenotype_xlsx = tmp_path / "phenotypes.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "WGIN_Watkins_JIC_CFLN06"
            ws.append(["StoreCode", "PH_M_cm-CFLN06", "GrwHabit_E_sw-CFLN06"])
            ws.append(["WATDE0001", 120, "w"])
            ws.append(["WATDE0002", 80, "s"])
            ws.append(["WATDE0003", 85, "s"])
            ws.append(["WATDE0004", 78, "w"])
            wb.save(phenotype_xlsx)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "intermediate"
            rc = prepare_main([
                "--vcf", str(vcf_path),
                "--phenotype-xlsx", str(phenotype_xlsx),
                "--output-root", str(output_root),
                "--intermediate-root", str(intermediate_root),
                "--target", "Rht-B1b",
                "--target", "Rht-D1b",
                "--min-haplotype-count", "1",
            ])

            self.assertEqual(rc, 0)
            for target_id, marker_id in [
                ("Rht-B1b", "chr4B_30861571"),
                ("Rht-D1b", "chr4D_18781242"),
            ]:
                db_dir = output_root / target_id
                self.assertTrue((db_dir / "gene_info.json").exists())
                gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
                self.assertEqual(gene_info["source"], "wheatomics_remote_rht1_functional_snp_vcf")
                self.assertIn("PlantHeight_CFLN06", (db_dir / "phenotype_data.csv").read_text(encoding="utf-8"))
                variant_info = pd.read_csv(db_dir / "variant_info.csv")
                self.assertEqual(variant_info.loc[0, "marker_id"], marker_id)
                self.assertEqual(variant_info.loc[0, "annotation"], "stop_gained")
                self.assertIn("T", (db_dir / "haplotype_data.csv").read_text(encoding="utf-8"))

    def test_wheat_rht1_prepare_continues_when_one_functional_snp_is_not_segregating(self):
        from openpyxl import Workbook
        from prepare_wheat2024_rht1_functional_snps import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            vcf_path = tmp_path / "rht1.vcf"
            vcf_path.write_text(
                "\n".join([
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWATDE0001\tWATDE0002\tWATDE0003\tWATDE0004",
                    "chr4B\t30861571\tchr4B_30861571\tC\tT\t.\t.\tANN=T|stop_gained|HIGH|TraesCS4B01G043100|TraesCS4B01G043100|transcript|TraesCS4B01G043100.1|protein_coding|1/1|c.190C>T|p.Gln64*|190/1866|190/1866|64/621||\tGT\t0/0\t0/0\t0/0\t0/0",
                    "chr4D\t18781242\tchr4D_18781242\tG\tT\t.\t.\tANN=T|stop_gained|HIGH|TraesCS4D01G040400|TraesCS4D01G040400|transcript|TraesCS4D01G040400.1|protein_coding|1/1|c.181G>T|p.Glu61*|181/1872|181/1872|61/623||\tGT\t0/0\t1/1\t0/0\t1/1",
                    "",
                ]),
                encoding="utf-8",
            )

            phenotype_xlsx = tmp_path / "phenotypes.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "WGIN_Watkins_JIC_CFLN06"
            ws.append(["StoreCode", "PH_M_cm-CFLN06", "GrwHabit_E_sw-CFLN06"])
            for sample, height in [
                ("WATDE0001", 120),
                ("WATDE0002", 80),
                ("WATDE0003", 118),
                ("WATDE0004", 78),
            ]:
                ws.append([sample, height, "s"])
            wb.save(phenotype_xlsx)

            output_root = tmp_path / "db"
            rc = prepare_main([
                "--vcf", str(vcf_path),
                "--phenotype-xlsx", str(phenotype_xlsx),
                "--output-root", str(output_root),
                "--intermediate-root", str(tmp_path / "intermediate"),
                "--target", "Rht-B1b",
                "--target", "Rht-D1b",
                "--min-haplotype-count", "1",
            ])

            self.assertEqual(rc, 0)
            self.assertFalse((output_root / "Rht-B1b" / "gene_info.json").exists())
            self.assertTrue((output_root / "Rht-D1b" / "gene_info.json").exists())

    def test_wheat_rht_zanke2014_prepare_builds_external_marker_database(self):
        import pandas as pd
        from openpyxl import Workbook
        from prepare_wheat2024_rht_zanke2014 import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_path = tmp_path / "zanke_table_s2.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "Table_S2"
            ws.append(["Table S2: List of varieties, genotyping data of candidate genes and phenotypic data."])
            ws.append([None, None, None, "Ellis et al., 2002, TAG", None, "Ellis et al., 2002, TAG"])
            ws.append([
                "Variety name", "GW no.", "habit", "Rht-B1", None, "Rht-D1", None,
                "Rht8", None, None, None, None, None, None, "Ppd-D1 wildtype",
                "Ppd-D1a mutant", "09.AND.PH", "09.SEL.PH", "09.WOH.PH",
                "10.AND.PH", "10.JAN.PH", "10.SAU.PH", "10.SEL.PH",
                "10.WOH.PH", "BLUES", "GAT_2012",
            ])
            ws.append([
                None, None, None, "mutant (dwarfed)", "wild type (tall)",
                "4D-mutant (DF-MR2)", "4D-wild type (DF2-WR2)",
            ])
            ws.append(["Tall1", "GW0001", "winter", None, 1, None, 1, None, None, None, None, None, None, None, 1, None, 110, 112, 108, 111, 109, 113, 107, 110, 110, 111])
            ws.append(["B1short", "GW0002", "winter", 1, None, None, 1, None, None, None, None, None, None, None, 1, None, 88, 90, 86, 89, 87, 91, 85, 88, 88, 89])
            ws.append(["D1short", "GW0003", "spring", None, 1, 1, None, None, None, None, None, None, None, None, 1, None, 82, 84, 80, 83, 81, 85, 79, 82, 82, 83])
            ws.append(["DoubleShort", "GW0004", "spring", 1, None, 1, None, None, None, None, None, None, None, None, 1, None, 70, 72, 68, 71, 69, 73, 67, 70, 70, 71])
            ws.append(["Sum", None, None, 2, 2, 2, 2])
            wb.save(source_path)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "intermediate"
            rc = prepare_main([
                "--source-workbook", str(source_path),
                "--output-root", str(output_root),
                "--intermediate-root", str(intermediate_root),
                "--min-haplotype-count", "1",
            ])

            self.assertEqual(rc, 0)
            db_dir = output_root / "Rht-Zanke2014"
            gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_info["source"], "zanke2014_rht_candidate_gene_table_s2")
            self.assertEqual(gene_info["marker_panel"], ["Rht-B1", "Rht-D1"])
            variant_info = pd.read_csv(db_dir / "variant_info.csv")
            self.assertEqual(variant_info["marker_id"].tolist(), ["Rht-B1", "Rht-D1"])
            self.assertTrue((variant_info["annotation"] == "functional_marker").all())
            self.assertIn("PlantHeight_BLUE", (db_dir / "phenotype_data.csv").read_text(encoding="utf-8"))
            self.assertIn("B1b|Rht-D1a", (db_dir / "haplotype_data.csv").read_text(encoding="utf-8"))

    def test_wheat_vrn_kiss2014_prepare_builds_diagnostic_marker_database(self):
        from openpyxl import Workbook
        from prepare_wheat2024_vrn_kiss2014 import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_path = tmp_path / "kiss2014.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "stepwise_reg2011_2012a"
            ws.append([
                "Pedigre", "VRN-A1", "VRN-B1", "VRN-D1", "PPD-D1", "PPDB1",
                "DEV49_2011", "DEV49_2012",
            ])
            ws.append(["WINTER-1", 0, 0, 0, 0, 0, 214, 216])
            ws.append(["VRND-SPRING", 0, 0, 1, 0, 0, 205, 207])
            ws.append(["VRNA-SPRING", 1, 0, 0, 0, 0, 212, 211])
            ws.append(["VRNB-SPRING", 0, 1, 0, 0, 0, 213, 212])
            ws2 = wb.create_sheet("stepwise_reg2011_2012b")
            ws2.append([
                "Pedigre", "VRN-A1", "VRN-B1", "VRN-D1", "PPD-D1", "PPDB1",
                "DEV59_2011", "DEV59_2012",
            ])
            ws2.append(["WINTER-1", 0, 0, 0, 0, 0, 222, 224])
            ws2.append(["VRND-SPRING", 0, 0, 1, 0, 0, 214, 216])
            ws2.append(["VRNA-SPRING", 1, 0, 0, 0, 0, 220, 221])
            ws2.append(["VRNB-SPRING", 0, 1, 0, 0, 0, 221, 220])
            wb.save(source_path)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "intermediate"
            rc = prepare_main([
                "--source-workbook", str(source_path),
                "--output-root", str(output_root),
                "--intermediate-root", str(intermediate_root),
                "--min-haplotype-count", "1",
            ])

            self.assertEqual(rc, 0)
            db_dir = output_root / "VRN-Kiss2014"
            gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_info["source"], "kiss2014_vrn_diagnostic_markers")
            self.assertEqual(gene_info["source_workbook"], str(source_path))
            self.assertEqual(gene_info["marker_panel"], ["VRN-A1", "VRN-B1", "VRN-D1"])
            self.assertTrue((intermediate_root / "VRN-Kiss2014" / "marker_matrix.tsv").exists())
            self.assertTrue((intermediate_root / "VRN-Kiss2014" / "phenotype.tsv").exists())

            variant_info = (db_dir / "variant_info.csv").read_text(encoding="utf-8")
            self.assertIn("VRN-D1", variant_info)
            self.assertIn("diagnostic_marker", variant_info)
            phenotype_data = (db_dir / "phenotype_data.csv").read_text(encoding="utf-8")
            self.assertIn("DEV49_mean", phenotype_data)
            self.assertIn("DEV59_mean", phenotype_data)
            hap_text = (db_dir / "haplotype_data.csv").read_text(encoding="utf-8")
            self.assertIn("0|0|1", hap_text)

    def test_wheat_vrn_kiss2014_prepare_builds_single_marker_targets(self):
        from openpyxl import Workbook
        from prepare_wheat2024_vrn_kiss2014 import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_path = tmp_path / "kiss2014.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "stepwise_reg2011_2012a"
            ws.append([
                "Pedigre", "VRN-A1", "VRN-B1", "VRN-D1", "PPD-D1", "PPDB1",
                "DEV49_2011", "DEV49_2012",
            ])
            ws.append(["WINTER-1", 0, 0, 0, 0, 0, 214, 216])
            ws.append(["VRND-SPRING", 0, 0, 1, 0, 0, 205, 207])
            ws.append(["VRNA-SPRING", 1, 0, 0, 0, 0, 212, 211])
            ws.append(["VRNB-SPRING", 0, 1, 0, 0, 0, 213, 212])
            ws2 = wb.create_sheet("stepwise_reg2011_2012b")
            ws2.append([
                "Pedigre", "VRN-A1", "VRN-B1", "VRN-D1", "PPD-D1", "PPDB1",
                "DEV59_2011", "DEV59_2012",
            ])
            ws2.append(["WINTER-1", 0, 0, 0, 0, 0, 222, 224])
            ws2.append(["VRND-SPRING", 0, 0, 1, 0, 0, 214, 216])
            ws2.append(["VRNA-SPRING", 1, 0, 0, 0, 0, 220, 221])
            ws2.append(["VRNB-SPRING", 0, 1, 0, 0, 0, 221, 220])
            wb.save(source_path)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "intermediate"
            rc = prepare_main([
                "--source-workbook", str(source_path),
                "--output-root", str(output_root),
                "--intermediate-root", str(intermediate_root),
                "--min-haplotype-count", "1",
                "--single-marker-targets",
            ])

            self.assertEqual(rc, 0)
            for marker in ["VRN-A1", "VRN-B1", "VRN-D1"]:
                db_dir = output_root / f"{marker}-Kiss2014"
                self.assertTrue((db_dir / "gene_info.json").exists())
                gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
                self.assertEqual(gene_info["gene_symbol"], marker)
                self.assertEqual(gene_info["marker_panel"], [marker])
                variant_info = (db_dir / "variant_info.csv").read_text(encoding="utf-8")
                self.assertIn(marker, variant_info)
                self.assertNotIn("|", (db_dir / "haplotype_data.csv").read_text(encoding="utf-8"))

    def test_precomputed_single_marker_haplotypes_load_as_strings(self):
        import pandas as pd
        from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "db" / "SingleMarker"
            db_dir.mkdir(parents=True)
            (db_dir / "variant_info.csv").write_text(
                "position,ref,alt,maf,missing_rate,len_diff,is_sv,annotation,marker_id\n"
                "1,0,1,0.25,0,0,False,diagnostic_marker,Marker1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_data.csv").write_text(
                "Haplotype_Seq,Count,Hap_Name,Alleles\n"
                "0,2,Hap1,0\n"
                "1,2,Hap2,1\n",
                encoding="utf-8",
            )
            (db_dir / "haplotype_samples.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name\n"
                "S1,0,Hap1\n"
                "S2,1,Hap2\n",
                encoding="utf-8",
            )
            (db_dir / "phenotype_data.csv").write_text(
                "SampleID,Haplotype_Seq,Hap_Name,Trait\n"
                "S1,0,Hap1,10\n"
                "S2,1,Hap2,20\n",
                encoding="utf-8",
            )
            (db_dir / "gene_info.json").write_text(
                json.dumps({"gene_id": "SingleMarker", "chrom": "diagnostic_marker_panel", "start": 1, "end": 1}),
                encoding="utf-8",
            )

            analyzer = HaplotypePhenotypeAnalyzer(
                vcf_file=str(tmp_path / "missing.vcf"),
                phenotype_file=str(db_dir / "phenotype_data.csv"),
                output_dir=str(tmp_path / "out"),
            )
            analyzer.analyze_gene(
                chrom="diagnostic_marker_panel",
                start=1,
                end=1,
                gene_id="SingleMarker",
                phenotype_cols=["Trait"],
                database_dir=str(tmp_path / "db"),
            )

            loaded = pd.read_csv(db_dir / "haplotype_samples.csv")
            self.assertIn(loaded["Haplotype_Seq"].dtype.kind, "iu")
            self.assertTrue((tmp_path / "out" / "SingleMarker.html").exists())
            scores = json.loads((tmp_path / "out" / "haplotype_scores.json").read_text(encoding="utf-8"))
            self.assertIn("Trait", scores)

    def test_wheat_tagw2_b1_prepare_builds_remote_snp_database(self):
        import pandas as pd
        from openpyxl import Workbook
        from prepare_wheat2024_tagw2_b1_remote_snp import main as prepare_main

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            vcf_path = tmp_path / "tagw2_b1.vcf"
            vcf_path.write_text(
                "\n".join([
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWATDE0001\tWATDE0002\tWATDE0003\tWATDE0004",
                    "chr6B\t291759689\tchr6B_291759689\tC\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t1/1\t0/0",
                    "chr6B\t291760677\tchr6B_291760677\tG\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1",
                    "chr6B\t291761315\tchr6B_291761315\tT\tC\t.\tPASS\t.\tGT\t1/1\t0/0\t1/1\t0/0",
                    "",
                ]),
                encoding="utf-8",
            )

            phenotype_xlsx = tmp_path / "watkins.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "WISP_Watkins_JIC_CFLN10"
            ws.append(["StoreCode", "GW_M_g1000grn-CFLN10"])
            ws.append(["WATDE0001", 48.0])
            ws.append(["WATDE0002", 40.0])
            ws.append(["WATDE0003", 47.0])
            ws.append(["WATDE0004", 41.0])
            ws2 = wb.create_sheet("DFW_Watkins_JI_CFLN14")
            ws2.append(["StoreCode", "GW_M_g1000grn-CFLN14"])
            ws2.append(["WATDE0001", 49.0])
            ws2.append(["WATDE0002", 39.0])
            ws2.append(["WATDE0003", 46.0])
            ws2.append(["WATDE0004", 42.0])
            wb.save(phenotype_xlsx)

            output_root = tmp_path / "db"
            intermediate_root = tmp_path / "intermediate"
            rc = prepare_main([
                "--vcf", str(vcf_path),
                "--phenotype-xlsx", str(phenotype_xlsx),
                "--output-root", str(output_root),
                "--intermediate-root", str(intermediate_root),
                "--min-haplotype-count", "1",
            ])

            self.assertEqual(rc, 0)
            db_dir = output_root / "TaGW2-B1-remoteSNP"
            gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_info["source"], "wheatomics_remote_tagw2_b1_promoter_snp_vcf")
            self.assertEqual(gene_info["atg_position"], 291761397)
            variant_info = pd.read_csv(db_dir / "variant_info.csv")
            self.assertEqual(
                variant_info["marker_id"].tolist(),
                ["chr6B_291759689", "chr6B_291760677", "chr6B_291761315"],
            )
            self.assertTrue((variant_info["annotation"] == "promoter_snp").all())
            hap_text = (db_dir / "haplotype_data.csv").read_text(encoding="utf-8")
            self.assertIn("A|G|C", hap_text)

    def test_tagw2_hap8_maps_published_promoter_offsets_to_refseq_positions(self):
        from analyze_tagw2_hap8_functional_groups import map_promoter_offset_to_position

        gene_start = 237_734_651
        cds_start_offset_1based = 185

        self.assertEqual(
            map_promoter_offset_to_position(gene_start, cds_start_offset_1based, -739),
            237_734_096,
        )
        self.assertEqual(
            map_promoter_offset_to_position(gene_start, cds_start_offset_1based, -593),
            237_734_242,
        )

    def test_tagw2_hap8_classifies_literature_and_background_patterns(self):
        import pandas as pd
        from analyze_tagw2_hap8_functional_groups import classify_literature_background_groups

        df = pd.DataFrame({
            "SampleID": ["S1", "S2", "S3", "S4"],
            "Hap_Name": ["Hap8", "Hap7", "Hap1", "Hap15"],
            "lit1": ["G", "G", "A", "G"],
            "lit2": ["A", "A", "G", "A"],
            "bg1": ["C", "T", "T", "C"],
            "bg2": ["T", "T", "T", "G"],
        })

        classified = classify_literature_background_groups(
            df,
            literature_alleles={"lit1": "G", "lit2": "A"},
            background_alleles={"bg1": "C", "bg2": "T"},
            full_haplotype_name="Hap8",
        )

        self.assertEqual(classified["PublishedPair"].tolist(), [True, True, False, True])
        self.assertEqual(classified["Hap8Background"].tolist(), [True, False, False, False])
        self.assertEqual(classified["FullHap8"].tolist(), [True, False, False, False])
        self.assertEqual(
            classified["Hap8SplitGroup"].tolist(),
            ["Published+Background", "PublishedOnly", "Neither", "PublishedOnly"],
        )

    def test_tagw2_hap8_binary_trait_summary_reports_effect_and_pvalue(self):
        import pandas as pd
        from analyze_tagw2_hap8_functional_groups import summarize_binary_trait

        df = pd.DataFrame({
            "SampleID": ["S1", "S2", "S3", "S4", "S5"],
            "group": [True, True, False, False, False],
            "TGW_mean": [30.0, 32.0, 40.0, 42.0, 44.0],
        })

        summary = summarize_binary_trait(df, "group", "TGW_mean", "case", "control")

        self.assertEqual(summary["case_label"], "case")
        self.assertEqual(summary["control_label"], "control")
        self.assertEqual(summary["case_n"], 2)
        self.assertEqual(summary["control_n"], 3)
        self.assertLess(summary["effect_case_minus_control"], 0)
        self.assertTrue(0 <= summary["welch_pvalue"] <= 1)


if __name__ == "__main__":
    unittest.main()
