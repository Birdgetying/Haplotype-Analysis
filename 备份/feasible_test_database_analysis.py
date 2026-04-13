#!/usr/bin/env python3
"""
测试从数据库加载数据运行分析
"""

import os
import sys
import json
import pandas as pd
import numpy as np

sys.path.insert(0, 'd:/Desktop/project1')

from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

# 配置
GENE_ID = "CSIAAS1BG1157200HC"
DATABASE_DIR = "d:/Desktop/project1/database"
OUTPUT_DIR = "d:/Desktop/project1/results2"

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 加载基因信息
gene_info_path = os.path.join(DATABASE_DIR, GENE_ID, 'gene_info.json')
with open(gene_info_path, 'r') as f:
    gene_info = json.load(f)

print(f"[INFO] 基因信息:")
print(f"  - ID: {gene_info['gene_id']}")
print(f"  - 染色体: {gene_info['chrom']}")
print(f"  - 区域: {gene_info['start']:,}-{gene_info['end']:,}")
print(f"  - 链方向: {gene_info['strand']}")
print(f"  - 外显子数: {len(gene_info.get('exons', []))}")
print(f"  - CDS数: {len(gene_info.get('cds', []))}")

# 创建分析器（使用数据库目录）
# 注意：VCF和GTF文件路径在本地不可用，但数据库已有预计算数据
analyzer = HaplotypePhenotypeAnalyzer(
    vcf_file=gene_info['vcf_file'],  # 数据库中已保存数据，不需要实际访问
    phenotype_file=None,  # 使用数据库中的phenotype_data.csv
    output_dir=OUTPUT_DIR,
    gtf_file=None  # 使用数据库中的gene_info.json
)

# **设置过滤参数**（可以在这里调整）
analyzer.max_missing_rate = 0.1  # 缺失率阈值（> 此值的变异被过滤）
analyzer.min_maf = 0.1  # MAF阈值（< 此值的变异被过滤）

print(f"[INFO] 过滤参数:")
print(f"  - max_missing_rate: {analyzer.max_missing_rate}")
print(f"  - min_maf: {analyzer.min_maf}")

# 从数据库加载表型数据
phenotype_data_path = os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')
if os.path.exists(phenotype_data_path):
    analyzer.phenotype_df = pd.read_csv(phenotype_data_path)
    print(f"[INFO] 从数据库加载表型数据: {len(analyzer.phenotype_df)} 个样本")
    print(f"[DEBUG] 表型数据列: {list(analyzer.phenotype_df.columns)}")
else:
    print(f"[WARNING] 数据库中缺少表型数据: {phenotype_data_path}")

# 检查 hap_sample_df 的列名
print(f"[DEBUG] hap_sample_df 列: {list(analyzer.hap_sample_df.columns) if analyzer.hap_sample_df is not None else 'None'}")

# 运行分析，传入database_dir以使用预计算数据
# 注意：使用数据库中实际存在的表型列名
result = analyzer.analyze_gene(
    chrom=gene_info['chrom'],
    start=gene_info['start'],
    end=gene_info['end'],
    gene_id=GENE_ID,
    phenotype_cols=['TFW_DSI'],  # 数据库中的实际列名
    database_dir=DATABASE_DIR
)

print(f"\n[INFO] 分析完成!")
print(f"  - 输出目录: {OUTPUT_DIR}")
print(f"  - 结果: {result}")
