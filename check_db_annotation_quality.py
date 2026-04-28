#!/usr/bin/env python
"""测试数据库中的注释质量"""
import os
import json
import pandas as pd
from pathlib import Path

db_dir = Path("database_text")
summary = pd.read_csv(db_dir / "summary.csv")

print("=" * 80)
print("数据库注释质量检查")
print("=" * 80)

# 检查每个基因的annotation分布
for _, row in summary.iterrows():
    gene_id = row['gene_id']
    var_file = db_dir / gene_id / "variant_info.csv"
    gene_info_file = db_dir / gene_id / "gene_info.json"
    
    if not var_file.exists():
        continue
    
    df = pd.read_csv(var_file)
    
    # 统计annotation分布
    ann_counts = df['annotation'].value_counts()
    
    # 检查是否有CDS
    with open(gene_info_file) as f:
        gene_info = json.load(f)
    
    has_cds = bool(gene_info.get('cds', []))
    n_cds_intervals = len(gene_info.get('cds', []))
    
    # 计算CDS范围
    if has_cds:
        cds_all = [pos for interval in gene_info['cds'] for pos in range(interval[0], interval[1]+1)]
        cds_variants = df[df['position'].isin(cds_all)]
        cds_ann = cds_variants['annotation'].value_counts() if len(cds_variants) > 0 else pd.Series()
    else:
        cds_variants = pd.DataFrame()
        cds_ann = pd.Series()
    
    # 计算启动子范围
    promoter_start = gene_info['gene_start'] - gene_info.get('promoter_length', 2000)
    promoter_end = gene_info['gene_start'] - 1
    promoter_variants = df[(df['position'] >= promoter_start) & (df['position'] <= promoter_end)]
    promoter_ann = promoter_variants['annotation'].value_counts() if len(promoter_variants) > 0 else pd.Series()
    
    print(f"\n{gene_id} (chr: {row['chrom']}, strand: {row['strand']})")
    print(f"  总变异数: {len(df)}, CDS区间数: {n_cds_intervals}")
    print(f"  全部annotation分布:")
    for ann, cnt in ann_counts.items():
        print(f"    {ann}: {cnt}")
    
    if len(cds_variants) > 0:
        print(f"  CDS内变异 ({len(cds_variants)}个):")
        for ann, cnt in cds_ann.items():
            print(f"    {ann}: {cnt}")
        
        # 检查问题
        if 'missense' not in cds_ann and 'synonymous' not in cds_ann:
            print(f"    ⚠️ 警告: CDS内变异没有missense/synonymous，全是: {list(cds_ann.index)}")
    
    if len(promoter_variants) > 0:
        print(f"  启动子区变异 ({len(promoter_variants)}个):")
        for ann, cnt in promoter_ann.items():
            print(f"    {ann}: {cnt}")

print("\n" + "=" * 80)
print("检查完成")
print("=" * 80)
