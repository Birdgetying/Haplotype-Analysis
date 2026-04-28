#!/usr/bin/env python
"""验证单个基因的CDS和exon坐标"""
import json
import pandas as pd

gene_id = "HORVU.MOREX.r3.1HG0080220"

with open(f"database_text/{gene_id}/gene_info.json") as f:
    gene_info = json.load(f)

df = pd.read_csv(f"database_text/{gene_id}/variant_info.csv")

print("=" * 80)
print(f"基因: {gene_id}")
print(f"基因体: {gene_info['gene_start']}-{gene_info['gene_end']} ({gene_info['strand']})")
print(f"Exon区间 ({len(gene_info['exons'])}个):")
for i, (s, e) in enumerate(gene_info['exons']):
    print(f"  Exon {i+1}: {s}-{e}")

print(f"\nCDS区间 ({len(gene_info['cds'])}个):")
for i, (s, e) in enumerate(gene_info['cds']):
    print(f"  CDS {i+1}: {s}-{e}")

# 合并所有CDS区域
cds_all = []
for s, e in gene_info['cds']:
    cds_all.extend(range(s, e+1))
cds_set = set(cds_all)

print(f"\nCDS总碱基数: {len(cds_set)}")
print(f"CDS范围: {min(cds_all)}-{max(cds_all)}")

# 检查每个变异是否在CDS内
print("\n变异位置检查:")
for _, row in df.iterrows():
    pos = row['position']
    in_cds = pos in cds_set
    in_gene = gene_info['gene_start'] <= pos <= gene_info['gene_end']
    
    # 检查在哪个exon/cds
    in_which_exon = []
    in_which_cds = []
    for i, (s, e) in enumerate(gene_info['exons']):
        if s <= pos <= e:
            in_which_exon.append(i+1)
    for i, (s, e) in enumerate(gene_info['cds']):
        if s <= pos <= e:
            in_which_cds.append(i+1)
    
    print(f"  Pos {pos}: 在基因内={in_gene}, 在CDS内={in_cds}, Exon={in_which_exon}, CDS={in_which_cds}, 注释={row['annotation']}")
