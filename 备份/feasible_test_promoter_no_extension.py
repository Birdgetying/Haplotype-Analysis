#!/usr/bin/env python3
"""
测试无扩展启动子（2000bp内就有变异）的分析逻辑

测试目标：
1. 启动子区域2000bp内有变异 → promoter_expansion_status = 'none'
2. 画图时使用2000bp启动子长度
3. 验证启动子区域正确绘制（橙色虚线框）
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
DATABASE_DIR = "d:/Desktop/project1/database_text"
BASE_OUTPUT_DIR = "d:/Desktop/project1/result_test"
VCF_FILE = "d:/Desktop/project1/chrALL.impute.vcf.gz"
GFF_FILE = "d:/Desktop/project1/barley_morex_v3.chr.gff3"

# 从 database_text 读取全部30个基因
gene_dirs = [d for d in os.listdir(DATABASE_DIR) if os.path.isdir(os.path.join(DATABASE_DIR, d))]
TEST_GENES = sorted([d for d in gene_dirs if os.path.exists(os.path.join(DATABASE_DIR, d, 'gene_info.json'))])

print(f"\n{'='*70}")
print(f"启动子无扩展测试（2000bp内有变异）")
print(f"{'='*70}")
print(f"[INFO] 测试基因数: {len(TEST_GENES)}")
print(f"[INFO] 预期行为: 启动子区域使用2000bp长度，不扩展\n")

# 批量分析所有测试基因
for GENE_ID in TEST_GENES:
    print(f"\n{'='*70}")
    print(f"[TEST] 测试基因: {GENE_ID}")
    print(f"{'='*70}\n")
    
    # 为每个基因创建独立输出目录（直接用基因ID命名）
    OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, GENE_ID)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"[INFO] 输出目录: {OUTPUT_DIR}")
    
    # 加载基因信息
    gene_info_path = os.path.join(DATABASE_DIR, GENE_ID, 'gene_info.json')
    if os.path.exists(gene_info_path):
        with open(gene_info_path, 'r') as f:
            gene_info = json.load(f)
        print(f"[INFO] 从 database 加载基因信息")
    else:
        print(f"[ERROR] 基因信息文件不存在: {gene_info_path}")
        continue

    print(f"[INFO] 基因信息:")
    print(f"  - ID: {gene_info['gene_id']}")
    print(f"  - 染色体: {gene_info['chrom']}")
    print(f"  - 区域: {gene_info['start']:,}-{gene_info['end']:,}")
    print(f"  - 基因体: {gene_info.get('gene_start', 'N/A'):,}-{gene_info.get('gene_end', 'N/A'):,}")
    print(f"  - 链方向: {gene_info['strand']}")
    print(f"  - 启动子长度: {gene_info.get('promoter_length', 'N/A')}")
    print(f"  - 扩展状态: {gene_info.get('promoter_expansion_status', 'N/A')}")
    print(f"  - 实际长度: {gene_info.get('promoter_actual_length', 'N/A')}")

    # 创建分析器
    analyzer = HaplotypePhenotypeAnalyzer(
        vcf_file=VCF_FILE,
        phenotype_file=None,
        output_dir=OUTPUT_DIR,
        gtf_file=GFF_FILE
    )

    # 设置过滤参数
    analyzer.max_missing_rate = 0.1
    analyzer.min_maf = 0.1

    # 从数据库加载表型数据
    phenotype_data_path = os.path.join(DATABASE_DIR, GENE_ID, 'phenotype_data.csv')
    if not os.path.exists(phenotype_data_path):
        print(f"[WARNING] 缺少表型数据: {phenotype_data_path}")
        continue
    
    analyzer.phenotype_df = pd.read_csv(phenotype_data_path)
    print(f"[INFO] 加载表型数据: {len(analyzer.phenotype_df)} 个样本")
    
    # 查找表型列
    pheno_cols = [c for c in analyzer.phenotype_df.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
    actual_pheno_col = pheno_cols[0] if pheno_cols else 'Phenotype_1'
    print(f"[INFO] 使用表型列: {actual_pheno_col}")

    # 运行分析
    print(f"\n[TEST] 开始分析...")
    result = analyzer.analyze_gene(
        chrom=gene_info['chrom'],
        start=gene_info['start'],
        end=gene_info['end'],
        gene_id=GENE_ID,
        phenotype_cols=[actual_pheno_col],
        database_dir=DATABASE_DIR
    )

    print(f"\n[TEST] {GENE_ID} 分析完成!")
    print(f"  - 输出目录: {OUTPUT_DIR}")
    print(f"  - 请检查 integrated_analysis.html 中的启动子区域")
    print(f"  - 预期: 橙色虚线框，长度2000bp\n")

print(f"\n{'='*70}")
print(f"[INFO] 所有测试完成！")
print(f"{'='*70}")
print(f"\n[验证清单]")
print(f"  1. 打开 result_test/{{基因ID}}/{{基因ID}}.html")
print(f"  2. 检查基因结构图上的启动子区域（橙色虚线框）")
print(f"  3. 验证启动子长度 = 2000bp（不是5000或其他值）")
print(f"  4. 检查日志中的 'promoter_actual_length=2000'")
print(f"\n")
