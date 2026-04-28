#!/usr/bin/env python
"""重新注释database_text中所有基因的variant_info.csv"""
import os
import json
import pandas as pd
from pathlib import Path

db_dir = Path("database_text")
summary = pd.read_csv(db_dir / "summary.csv")

print("=" * 80)
print("重新注释数据库中的变异")
print("=" * 80)

# 手动定义判断函数
def _pos_in_any_interval(pos, intervals):
    for s, e in intervals:
        if s <= pos <= e:
            return True
    return False

def classify_variant(ref, alt):
    if isinstance(alt, str) and alt.startswith('<') and alt.endswith('>'):
        return 'SV'
    len_diff = abs(len(ref) - len(alt))
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len_diff >= 50:
        return 'SV'
    elif len(alt) > len(ref):
        return 'INS'
    elif len(ref) > len(alt):
        return 'DEL'
    else:
        return 'SNP'

# 对每个基因重新注释
for _, row in summary.iterrows():
    gene_id = row['gene_id']
    gene_folder = db_dir / gene_id
    var_file = gene_folder / "variant_info.csv"
    gene_info_file = gene_folder / "gene_info.json"
    
    if not var_file.exists() or not gene_info_file.exists():
        print(f"\n[SKIP] {gene_id}: 文件不存在")
        continue
    
    # 加载数据
    var_df = pd.read_csv(var_file)
    with open(gene_info_file) as f:
        gene_info = json.load(f)
    
    cds_intervals = gene_info.get('cds', [])
    exon_intervals = gene_info.get('exons', [])
    strand = gene_info['strand']
    gene_start = gene_info['gene_start']
    gene_end = gene_info['gene_end']
    
    # 计算promoter
    if strand == '+':
        promoter_start = max(1, gene_start - gene_info.get('promoter_actual_length', 2000))
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + gene_info.get('promoter_actual_length', 2000)
    
    # 重新注释每个变异
    new_annotations = []
    for _, var in var_df.iterrows():
        pos = var['position']
        ref = var['ref']
        alt = var['alt']
        
        var_type = classify_variant(str(ref), str(alt))
        in_cds = _pos_in_any_interval(pos, cds_intervals)
        in_exon = _pos_in_any_interval(pos, exon_intervals)
        
        # 判断逻辑（无FASTA，无法判断missense/synonymous）
        if var_type == 'SV':
            ann = 'SV'
        elif var_type in ('INS', 'DEL'):
            ann = var_type
        elif not in_exon:
            if promoter_start <= pos <= promoter_end:
                ann = 'promoter'
            elif gene_start <= pos <= gene_end:
                ann = 'intron'
            else:
                ann = 'other'
        elif not in_cds:
            ann = 'UTR'
        else:
            # CDS内的SNP，无FASTA时标注为other
            ann = 'other'
        
        new_annotations.append(ann)
    
    # 更新DataFrame
    var_df['annotation'] = new_annotations
    
    # 保存
    var_df.to_csv(var_file, index=False)
    
    # 统计
    ann_counts = var_df['annotation'].value_counts()
    print(f"\n{gene_id}:")
    for ann, cnt in ann_counts.items():
        print(f"  {ann}: {cnt}")

print("\n" + "=" * 80)
print("重新注释完成")
print("=" * 80)
