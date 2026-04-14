#!/usr/bin/env python3
"""快速测试单个基因"""
import os
import sys
import json
import pandas as pd

sys.path.insert(0, 'd:/Desktop/project1')
from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

GENE_ID = 'HORVU.MOREX.r3.4HG0387570'
DATABASE_DIR = "d:/Desktop/project1/database"
OUTPUT_DIR = f"d:/Desktop/project1/results2/{GENE_ID}"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 加载基因信息
gene_info_path = os.path.join(DATABASE_DIR, GENE_ID, 'gene_info.json')
with open(gene_info_path, 'r') as f:
    gene_info = json.load(f)

# 创建分析器
analyzer = HaplotypePhenotypeAnalyzer(
    vcf_file=None,
    phenotype_file=None,
    output_dir=OUTPUT_DIR,
    gtf_file=None
)

analyzer.max_missing_rate = 0.1
analyzer.min_maf = 0.1

# 加载表型数据
phenotype_data_path = os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')
analyzer.phenotype_df = pd.read_csv(phenotype_data_path)

# 删除 VCF 标记文件
vcf_mark_path = os.path.join(DATABASE_DIR, GENE_ID, 'variants.vcf.gz')
if os.path.exists(vcf_mark_path) and os.path.getsize(vcf_mark_path) < 100:
    os.remove(vcf_mark_path)

# 获取表型列名
temp_pheno_df = pd.read_csv(phenotype_data_path)
pheno_cols = [c for c in temp_pheno_df.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
actual_pheno_col = pheno_cols[0] if pheno_cols else 'Phenotype_1'

print(f"分析基因: {GENE_ID}")
print(f"表型列: {actual_pheno_col}")

result = analyzer.analyze_gene(
    chrom=gene_info['chrom'],
    start=gene_info['start'],
    end=gene_info['end'],
    gene_id=GENE_ID,
    phenotype_cols=[actual_pheno_col],
    database_dir=DATABASE_DIR
)

print(f"\n✓ 分析完成!")
print(f"  输出: {OUTPUT_DIR}/integrated_analysis.html")

# 验证
html_file = f"{OUTPUT_DIR}/integrated_analysis.html"
with open(html_file, 'r', encoding='utf-8') as f:
    content = f.read()

print(f"\n验证结果:")
print(f"  Copy All 按钮: {'✓' if 'Copy All' in content else '✗'}")
print(f"  Copy 按钮 70px: {'✓' if 'width:70px' in content else '✗'}")
print(f"  font-weight:600: {'✓' if 'font-weight:600' in content else '✗'}")
print(f"  font-weight:500: {'✓' if 'font-weight:500' in content else '✗'}")
print(f"  无 Copied!: {'✓' if 'Copied!' not in content else '✗'}")
