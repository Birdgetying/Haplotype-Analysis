#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检查 database 目录是否包含老师所需的30个基因
"""

import pandas as pd
import os

# Excel 文件路径
excel_file = r'd:\Desktop\project1\database\大麦备份\大麦GWAS及SV筛选基因-2026.03.22(1).xlsx'

# 提取老师所需的30个基因
df1 = pd.read_excel(excel_file, sheet_name='K-SNP_Indel_CandidateGene')
df2 = pd.read_excel(excel_file, sheet_name='K-SV（自己结果）')

genes1 = df1.iloc[:, 0].dropna().tolist()
genes2 = df2.iloc[:, 0].dropna().tolist()
all_genes = set(genes1 + genes2)

print("=" * 80)
print("老师所需的30个基因")
print("=" * 80)
for i, g in enumerate(sorted(all_genes), 1):
    print(f"  {i:2d}. {g}")

# 检查 database 目录
db_dir = r'd:\Desktop\project1\database'
print(f"\n{'=' * 80}")
print(f"检查数据库目录: {db_dir}")
print(f"{'=' * 80}")

# 获取所有基因文件夹
gene_dirs = [d for d in os.listdir(db_dir) 
             if os.path.isdir(os.path.join(db_dir, d)) and d.startswith('HORVU')]

print(f"\n数据库中的基因文件夹: {len(gene_dirs)} 个")

# 检查匹配
found_genes = []
missing_genes = []
extra_genes = []

for gene in all_genes:
    if gene in gene_dirs:
        found_genes.append(gene)
    else:
        missing_genes.append(gene)

for gene in gene_dirs:
    if gene not in all_genes:
        extra_genes.append(gene)

print(f"\n{'=' * 80}")
print(f"检查结果")
print(f"{'=' * 80}")

print(f"\n✅ 找到的基因: {len(found_genes)}/{len(all_genes)}")
for g in sorted(found_genes):
    print(f"  ✓ {g}")

if missing_genes:
    print(f"\n❌ 缺失的基因: {len(missing_genes)}")
    for g in sorted(missing_genes):
        print(f"  ✗ {g}")
else:
    print(f"\n🎉 所有 {len(all_genes)} 个基因都在数据库中！")

if extra_genes:
    print(f"\n📁 额外的基因文件夹 (不在30个基因列表中): {len(extra_genes)}")
    for g in sorted(extra_genes):
        print(f"  + {g}")

# 统计
print(f"\n{'=' * 80}")
print(f"统计")
print(f"{'=' * 80}")
print(f"老师所需基因数: {len(all_genes)}")
print(f"数据库中已有: {len(found_genes)}")
print(f"缺失: {len(missing_genes)}")
print(f"数据库中的总基因数: {len(gene_dirs)}")
print(f"额外基因: {len(extra_genes)}")

if len(found_genes) == len(all_genes):
    print(f"\n✅✅✅ 完美！新数据库包含老师所需的全部 {len(all_genes)} 个基因！")
else:
    print(f"\n⚠️  还需要添加 {len(missing_genes)} 个基因")
