#!/usr/bin/env python3
"""
水稻单倍型分析测试脚本
使用 D:\Desktop\data\水稻 中的数据，分析3个测试基因
"""

import sys
import os

# Fix Windows encoding issue
if sys.platform == 'win32':
    import io
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper):
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper):
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

sys.path.insert(0, r'd:\Desktop\project1')
import pandas as pd
from genome_wide_haplotype_scan import (
    ScanConfig, parse_gff3_genes, process_single_gene,
    HAPLOTYPE_MODULE_AVAILABLE
)

# ========== 覆盖 ScanConfig ==========
ScanConfig.VCF_FILE = r'D:\Desktop\data\水稻\VCF format\rice4k_geno_add_del.vcf.gz'
ScanConfig.GFF_FILE = r'D:\Desktop\data\水稻\rice_test_genes.gff3'
ScanConfig.PHENO_FILE = r'D:\Desktop\data\水稻\Phenotypes\phenos_modified.tsv'
ScanConfig.FASTA_FILE = ''  # 无FASTA，跳过missense注释

# 测试模式
ScanConfig.TARGET_GENES = [
    'LOC_Os01g02660',
    'LOC_Os01g25140',
    'LOC_Os01g57340',
]

# 输出目录
ScanConfig.DATABASE_DIR = r'D:\Desktop\project1\rice_database'
ScanConfig.RESULTS_DIR = r'D:\Desktop\project1\rice_results'

# 分析参数
ScanConfig.PROMOTER_LENGTH = 2000
ScanConfig.PROMOTER_EXTENDED_LENGTH = 5000
ScanConfig.SNP_ONLY = False
ScanConfig.MIN_SAMPLES = 1
ScanConfig.PVALUE_THRESHOLD = 0.05
ScanConfig.TEST_REGION_LENGTH = 0

from genome_wide_haplotype_scan import run_genome_scan

if __name__ == '__main__':
    print("=" * 60)
    print("水稻单倍型分析 - 测试运行")
    print("=" * 60)
    print(f"VCF: {ScanConfig.VCF_FILE}")
    print(f"GFF: {ScanConfig.GFF_FILE}")
    print(f"Phenotype: {ScanConfig.PHENO_FILE}")
    print(f"Target genes: {ScanConfig.TARGET_GENES}")
    print(f"Database dir: {ScanConfig.DATABASE_DIR}")
    print(f"Results dir: {ScanConfig.RESULTS_DIR}")
    print(f"HAPLOTYPE_MODULE_AVAILABLE: {HAPLOTYPE_MODULE_AVAILABLE}")
    print()

    # 先测试GFF解析
    genes_df = parse_gff3_genes(ScanConfig.GFF_FILE)
    print(f"\nParsed {len(genes_df)} genes from GFF:")
    print(genes_df.to_string())

    # 运行
    results = run_genome_scan(
        vcf_file=ScanConfig.VCF_FILE,
        gff_file=ScanConfig.GFF_FILE,
        pheno_file=ScanConfig.PHENO_FILE,
        database_dir=ScanConfig.DATABASE_DIR,
        results_dir=ScanConfig.RESULTS_DIR,
        gene_filter=ScanConfig.TARGET_GENES,
        min_samples=ScanConfig.MIN_SAMPLES,
        run_analysis=True,
        pvalue_threshold=ScanConfig.PVALUE_THRESHOLD,
        test_region_length=ScanConfig.TEST_REGION_LENGTH,
        cluster_haplotypes=False,
        fasta_file=ScanConfig.FASTA_FILE if ScanConfig.FASTA_FILE else None,
    )

    print("\nDone!")
    print(results.to_string())
