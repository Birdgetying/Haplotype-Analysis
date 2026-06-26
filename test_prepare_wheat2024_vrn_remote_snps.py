#!/usr/bin/env python3
"""Tests for VRN WheatOmics SNP VCF database preparation."""

import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd


class WheatVrnRemoteSnpPrepareTests(unittest.TestCase):
    def test_build_target_database_from_single_gene_snp_vcf(self):
        from prepare_wheat2024_vrn_remote_snps import TARGETS, build_target_database

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            vcf_path = tmp_path / "vrn-a1.test.vcf.gz"
            with gzip.open(vcf_path, "wt", encoding="utf-8", newline="") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWAT001\tWAT002\tWAT003\n")
                f.write("chr5A\t587409581\tchr5A_587409581\tC\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/0\n")
                f.write("chr5A\t587410000\tchr5A_587410000\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t1/1\n")
                f.write("chr5A\t587500000\toutside\tG\tA\t.\t.\t.\tGT\t0/0\t1/1\t0/0\n")

            phenotype_df = pd.DataFrame({
                "SampleID": ["WAT001", "WAT002", "WAT003"],
                "GrowthHabitSpringScore_CFLN06": [1.0, 0.0, 0.75],
                "PlantHeight_CFLN06": [91.0, 85.0, 88.0],
            })

            db_dir = build_target_database(
                target_id="VRN-A1-remoteSNP",
                target=TARGETS["VRN-A1-remoteSNP"],
                vcf_path=vcf_path,
                phenotype_df=phenotype_df,
                phenotype_xlsx=tmp_path / "phenotype.xlsx",
                output_root=tmp_path / "db",
                intermediate_root=tmp_path / "intermediate",
                min_haplotype_count=1,
                max_missing_rate=0.2,
            )

            variant_info = pd.read_csv(db_dir / "variant_info.csv")
            phenotype_data = pd.read_csv(db_dir / "phenotype_data.csv")
            gene_info = pd.read_json(db_dir / "gene_info.json", typ="series").to_dict()

            self.assertEqual(len(variant_info), 2)
            self.assertEqual(len(phenotype_data), 3)
            self.assertEqual(gene_info["gene_id"], "VRN-A1-remoteSNP")
            self.assertEqual(gene_info["traes_id"], "TraesCS5A01G391700")
            self.assertEqual(gene_info["source"], "wheatomics_remote_vrn_snp_vcf")
            self.assertIn("GrowthHabitSpringScore_CFLN06", phenotype_data.columns)


if __name__ == "__main__":
    unittest.main()
