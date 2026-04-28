#!/usr/bin/env python
"""深度检查数据库注释质量和逻辑一致性"""
import os
import json
import pandas as pd
from pathlib import Path

db_dir = Path("database_text")
summary = pd.read_csv(db_dir / "summary.csv")

print("=" * 80)
print("数据库注释深度检查")
print("=" * 80)

errors = []

for _, row in summary.iterrows():
    gene_id = row['gene_id']
    gene_folder = db_dir / gene_id
    var_file = gene_folder / "variant_info.csv"
    gene_info_file = gene_folder / "gene_info.json"
    
    if not var_file.exists() or not gene_info_file.exists():
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
    promoter_actual_length = gene_info.get('promoter_actual_length', 2000)
    
    # 计算promoter
    if strand == '+':
        promoter_start = max(1, gene_start - promoter_actual_length)
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + promoter_actual_length
    
    print(f"\n{'='*60}")
    print(f"基因: {gene_id} (strand={strand})")
    print(f"基因体: {gene_start}-{gene_end}")
    print(f"Promoter: {promoter_start}-{promoter_end}")
    print(f"CDS区间数: {len(cds_intervals)}, Exon区间数: {len(exon_intervals)}")
    
    # 合并CDS和Exon区域
    cds_all = set()
    for s, e in cds_intervals:
        cds_all.update(range(s, e+1))
    
    exon_all = set()
    for s, e in exon_intervals:
        exon_all.update(range(s, e+1))
    
    # 逐行检查每个变异
    for idx, var in var_df.iterrows():
        pos = var['position']
        ref = str(var['ref'])
        alt = str(var['alt'])
        ann = var['annotation']
        
        # 计算实际位置
        in_cds = pos in cds_all
        in_exon = pos in exon_all
        in_promoter = promoter_start <= pos <= promoter_end
        in_gene = gene_start <= pos <= gene_end
        
        # 计算变异类型
        len_diff = abs(len(ref) - len(alt))
        is_symbolic = ref.startswith('<') and ref.endswith('>') or alt.startswith('<') and alt.endswith('>')
        
        if is_symbolic:
            expected_type = 'SV'
        elif len_diff >= 50:
            expected_type = 'SV'
        elif len(alt) > len(ref) and len_diff > 0:
            expected_type = 'INS'
        elif len(ref) > len(alt) and len_diff > 0:
            expected_type = 'DEL'
        else:
            expected_type = 'SNP'
        
        # 根据逻辑判断期望的注释
        if expected_type == 'SV':
            expected_ann = 'SV'
        elif expected_type in ('INS', 'DEL'):
            # INS/DEL应该保留类型，不区分位置
            expected_ann = expected_type
        elif not in_exon:
            if in_promoter:
                expected_ann = 'promoter'
            elif in_gene:
                expected_ann = 'intron'
            else:
                expected_ann = 'other'
        elif not in_cds:
            expected_ann = 'UTR'
        else:
            # CDS内的SNP，无FASTA时应该是other
            expected_ann = 'other'
        
        # 检查是否匹配
        if ann != expected_ann:
            error_msg = f"  ❌ {gene_id} Pos {pos}: 期望={expected_ann}, 实际={ann} (in_cds={in_cds}, in_exon={in_exon}, in_promoter={in_promoter}, var_type={expected_type})"
            print(error_msg)
            errors.append(error_msg)
    
    # 统计
    ann_counts = var_df['annotation'].value_counts()
    print(f"  注释分布: {dict(ann_counts)}")

print(f"\n{'='*80}")
print(f"检查完成")
print(f"发现 {len(errors)} 个错误")
if errors:
    print("\n错误详情:")
    for err in errors[:20]:  # 只显示前20个
        print(err)
    if len(errors) > 20:
        print(f"  ... 还有 {len(errors)-20} 个错误")
print("=" * 80)
