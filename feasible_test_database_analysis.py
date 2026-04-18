#!/usr/bin/env python3
"""
测试从数据库加载数据运行分析（本地版）
"""

import os
import sys
import json
import pandas as pd
import numpy as np

# 设置UTF-8输出编码
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr.encoding != 'utf-8':
    sys.stderr.reconfigure(encoding='utf-8')

sys.path.insert(0, 'd:/Desktop/project1')

from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

# 配置
DATABASE_DIR = "d:/Desktop/project1/database"
BASE_OUTPUT_DIR = "d:/Desktop/project1/results3"
VCF_FILE = "d:/Desktop/project1/chrALL.impute.vcf.gz"  # 需要VCF文件用于HTML生成
GFF_FILE = "d:/Desktop/project1/barley_morex_v3.chr.gff3"  # 需要GFF文件用于HTML生成

# 协变量文件配置（如果需要）
COPHE_FILES = [
    ("d:/Desktop/project1/WT_TPM.txt", "WT_TPM"),
    ("d:/Desktop/project1/SA_TPM.txt", "SA_TPM"),
]

# 老师指定的30个基因（从 Excel 提取）
TEACHER_GENES = [
    # K-SNP_Indel_CandidateGene (16个)
    "HORVU.MOREX.r3.1HG0081480", "HORVU.MOREX.r3.2HG0213430", "HORVU.MOREX.r3.2HG0214000",
    "HORVU.MOREX.r3.3HG0289370", "HORVU.MOREX.r3.4HG0341890", "HORVU.MOREX.r3.4HG0356420",
    "HORVU.MOREX.r3.4HG0356370", "HORVU.MOREX.r3.4HG0398010", "HORVU.MOREX.r3.4HG0415770",
    "HORVU.MOREX.r3.6HG0625930", "HORVU.MOREX.r3.6HG0626440", "HORVU.MOREX.r3.7HG0705960",
    "HORVU.MOREX.r3.7HG0733200", "HORVU.MOREX.r3.4HG0415690", "HORVU.MOREX.r3.4HG0415780",
    "HORVU.MOREX.r3.5HG0492640",
    # K-SV（自己结果）(15个，其中 HORVU.MOREX.r3.4HG0356420 重复)
    "HORVU.MOREX.r3.6HG0542030", "HORVU.MOREX.r3.6HG0542250", "HORVU.MOREX.r3.4HG0356430",
    "HORVU.MOREX.r3.3HG0299830", "HORVU.MOREX.r3.3HG0299890", "HORVU.MOREX.r3.1HG0080620",
    "HORVU.MOREX.r3.1HG0080220", "HORVU.MOREX.r3.2HG0096990", "HORVU.MOREX.r3.7HG0639230",
    "HORVU.MOREX.r3.5HG0432290", "HORVU.MOREX.r3.5HG0490420", "HORVU.MOREX.r3.5HG0490030",
    "HORVU.MOREX.r3.5HG0490190", "HORVU.MOREX.r3.6HG0626350"
]

# 优先使用 database 目录，检查老师指定的基因是否存在
if os.path.exists(DATABASE_DIR):
    gene_folders = [d for d in os.listdir(DATABASE_DIR) 
                    if os.path.isdir(os.path.join(DATABASE_DIR, d)) and 
                    os.path.exists(os.path.join(DATABASE_DIR, d, 'gene_info.json'))]
else:
    gene_folders = []

# 如果 database 目录为空，使用 results2 中的现有基因
if not gene_folders:
    print("[WARNING] database 目录为空，使用 results2 中的现有基因")
    if os.path.exists(BASE_OUTPUT_DIR):
        gene_folders = [d for d in os.listdir(BASE_OUTPUT_DIR) 
                        if os.path.isdir(os.path.join(BASE_OUTPUT_DIR, d))]
    else:
        print(f"[ERROR] results2 目录也不存在: {BASE_OUTPUT_DIR}")
        exit(1)

print(f"[INFO] 数据库中找到 {len(gene_folders)} 个基因")

# 检查老师指定的30个基因是否在数据库中
teacher_genes_found = [g for g in TEACHER_GENES if g in gene_folders]
teacher_genes_missing = [g for g in TEACHER_GENES if g not in gene_folders]

print(f"\n{'='*70}")
print(f"老师指定的30个基因检查")
print(f"{'='*70}")
print(f"✅ 找到: {len(teacher_genes_found)}/{len(TEACHER_GENES)} 个")
if teacher_genes_missing:
    print(f"❌ 缺失: {len(teacher_genes_missing)} 个")
    for g in teacher_genes_missing:
        print(f"  - {g}")
else:
    print(f"🎉 所有30个基因都在数据库中！")

print()

# 使用老师指定的30个基因进行测试
ANALYZE_ALL = False  # 只分析老师指定的基因
# TARGET_GENES = teacher_genes_found  # 使用找到的老师基因
TARGET_GENES = ["HORVU.MOREX.r3.1HG0080220"]  # 只测试这一个基因

if ANALYZE_ALL:
    genes_to_analyze = gene_folders
else:
    genes_to_analyze = [g for g in TARGET_GENES if g in gene_folders]

print(f"[INFO] 将分析 {len(genes_to_analyze)} 个基因（老师指定的30个基因）")
print(f"[INFO] 结果文件夹将添加[Teacher]前缀便于查找\n")

# 批量分析所有基因
for GENE_ID in genes_to_analyze:
    print(f"\n{'='*70}")
    print(f"[INFO] 开始分析: {GENE_ID}")
    print(f"{'='*70}\n")
    
    # 为每个基因创建独立输出目录（添加[Teacher]前缀）
    OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, f"[Teacher]_{GENE_ID}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"[INFO] 输出目录: {OUTPUT_DIR}")
    
    # 加载基因信息
    gene_info_path = os.path.join(DATABASE_DIR, GENE_ID, 'gene_info.json')
    if os.path.exists(gene_info_path):
        with open(gene_info_path, 'r') as f:
            gene_info = json.load(f)
        print(f"[INFO] 从 database 加载基因信息")
    else:
        # 从 results2 中使用，尝试从文件名提取信息
        print(f"[INFO] 使用 results2 中的现有数据")
        gene_info = {
            'gene_id': GENE_ID,
            'chrom': GENE_ID.split('.')[3] if len(GENE_ID.split('.')) > 3 else 'Unknown',
            'start': 0,
            'end': 1000000,
            'strand': '+'
        }

    print(f"[INFO] 基因信息:")
    print(f"  - ID: {gene_info['gene_id']}")
    print(f"  - 染色体: {gene_info['chrom']}")
    print(f"  - 区域: {gene_info['start']:,}-{gene_info['end']:,}")
    print(f"  - 链方向: {gene_info['strand']}")
    print(f"  - 外显子数: {len(gene_info.get('exons', []))}")
    print(f"  - CDS数: {len(gene_info.get('cds', []))}")

    # 创建分析器（使用数据库目录）
    # 注意：VCF文件是标记文件，需要从CSV加载数据
    analyzer = HaplotypePhenotypeAnalyzer(
        vcf_file=VCF_FILE,  # 使用真实VCF文件用于HTML生成
        phenotype_file=None,  # 使用数据库中的phenotype_data.csv
        output_dir=OUTPUT_DIR,
        gtf_file=GFF_FILE  # 使用GFF文件用于HTML生成
    )

    # **设置过滤参数**（可以在这里调整）
    analyzer.max_missing_rate = 0.1  # 缺失率阈值（> 此值的变异被过滤）
    analyzer.min_maf = 0.1  # MAF阈值（< 此值的变异被过滤）

    print(f"[INFO] 过滤参数:")
    print(f"  - max_missing_rate: {analyzer.max_missing_rate}")
    print(f"  - min_maf: {analyzer.min_maf}")

    # 从数据库加载表型数据
    phenotype_data_path = os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')
    if not os.path.exists(phenotype_data_path):
        # 尝试从 results2 加载
        phenotype_data_path = os.path.join(BASE_OUTPUT_DIR, GENE_ID, 'sample_haplotypes.csv')
    
    if os.path.exists(phenotype_data_path):
        analyzer.phenotype_df = pd.read_csv(phenotype_data_path)
        print(f"[INFO] 加载表型数据: {phenotype_data_path}")
        print(f"[INFO] {len(analyzer.phenotype_df)} 个样本")
        print(f"[DATA FORMAT] 表型数据格式")
        print(f"  列名: {list(analyzer.phenotype_df.columns)}")
        print(f"  数据类型: {dict(analyzer.phenotype_df.dtypes)}")
        print(f"  样本ID示例 (前5个): {list(analyzer.phenotype_df['SampleID'].head())}")
        print(f"  样本ID数据类型: {analyzer.phenotype_df['SampleID'].dtype}")
        # 查找表型列（排除 SampleID, Hap_Name, Haplotype_Seq）
        pheno_cols = [c for c in analyzer.phenotype_df.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
        print(f"  表型列: {pheno_cols}")
        if pheno_cols:
            print(f"  表型列 '{pheno_cols[0]}' 详情:")
            print(f"    数据类型: {analyzer.phenotype_df[pheno_cols[0]].dtype}")
            print(f"    非空值数: {analyzer.phenotype_df[pheno_cols[0]].notna().sum()}")
            print(f"    统计: mean={analyzer.phenotype_df[pheno_cols[0]].mean():.4f}, std={analyzer.phenotype_df[pheno_cols[0]].std():.4f}")
    else:
        print(f"[WARNING] 缺少表型数据: {phenotype_data_path}")
        continue  # 跳过这个基因，继续下一个

    # 检查并删除标记文件（避免分析器误读）
    vcf_mark_path = os.path.join(DATABASE_DIR, GENE_ID, 'variants.vcf.gz')
    if os.path.exists(vcf_mark_path):
        vcf_size = os.path.getsize(vcf_mark_path)
        if vcf_size < 100:  # 标记文件通常很小
            print(f"[INFO] 检测到 VCF 标记文件 ({vcf_size} bytes)，删除后使用 CSV 数据")
            os.remove(vcf_mark_path)
            print(f"[INFO] 已删除: {vcf_mark_path}")

    # 运行分析，传入database_dir以使用预计算数据
    # 注意：使用数据库中实际的表型列名
    if os.path.exists(os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')):
        phenotype_data_path = os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')
    else:
        phenotype_data_path = os.path.join(BASE_OUTPUT_DIR, GENE_ID, 'sample_haplotypes.csv')
    
    temp_pheno_df = pd.read_csv(phenotype_data_path)
    pheno_cols = [c for c in temp_pheno_df.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
    actual_pheno_col = pheno_cols[0] if pheno_cols else 'Phenotype_1'
    
    # 过滤掉协变量列（只保留第一个表型列）
    covariate_cols = [c for c in pheno_cols if c != actual_pheno_col and ('WT' in c or 'SA' in c or 'TPM' in c or 'expr' in c.lower())]
    if covariate_cols:
        print(f"[INFO] 检测到协变量列: {covariate_cols}")
        print(f"[INFO] 数据库已包含协变量数据，无需额外加载")

    print(f"\n[INFO] 使用表型列: {actual_pheno_col}")

    # 注意：数据库中的 phenotype_data.csv 已包含协变量列，analyze_gene 会自动读取
    result = analyzer.analyze_gene(
        chrom=gene_info['chrom'],
        start=gene_info['start'],
        end=gene_info['end'],
        gene_id=GENE_ID,
        phenotype_cols=[actual_pheno_col],  # 使用实际的列名
        database_dir=DATABASE_DIR
    )

    print(f"\n[INFO] {GENE_ID} 分析完成!")
    print(f"  - 输出目录: {OUTPUT_DIR}")
    print(f"  - 结果: {result}\n")

print(f"\n{'='*70}")
print(f"[INFO] 所有基因分析完成！")
print(f"{'='*70}")
