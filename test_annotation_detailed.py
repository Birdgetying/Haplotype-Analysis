# -*- coding: utf-8 -*-
"""详细测试同义/错义突变注释功能"""
import os
import sys

sys.path.insert(0, 'd:/Desktop/project1')

def test_database_variant_info():
    """测试数据库中的variant_info"""
    
    test_gene = 'HORVU.MOREX.r3.2HG0213430'
    test_db_dir = 'd:/Desktop/project1/database_text'
    
    import pandas as pd
    import json
    
    gene_db_dir = os.path.join(test_db_dir, test_gene)
    variant_info_path = os.path.join(gene_db_dir, 'variant_info.csv')
    
    print("=" * 60)
    print("测试数据库 variant_info")
    print("=" * 60)
    
    # 读取variant_info
    df = pd.read_csv(variant_info_path)
    print(f"\n[INFO] variant_info.csv 统计:")
    print(f"  - 总变异数: {len(df)}")
    print(f"  - annotation分布:")
    for ann, count in df['annotation'].value_counts().items():
        print(f"    {ann}: {count}")
    
    # 读取gene_info
    gene_info_path = os.path.join(gene_db_dir, 'gene_info.json')
    with open(gene_info_path, 'r') as f:
        gene_info = json.load(f)
    
    print(f"\n[INFO] gene_info.json:")
    print(f"  - 染色体: {gene_info.get('chrom')}")
    print(f"  - 基因位置: {gene_info.get('gene_start')}-{gene_info.get('gene_end')}")
    print(f"  - CDS数量: {len(gene_info.get('cds', []))}")
    print(f"  - 外显子数量: {len(gene_info.get('exons', []))}")
    print(f"  - FASTA: {gene_info.get('fasta_file', 'N/A')}")
    
    # 分析CDS区域内有len_diff=0的变异（潜在SNP）
    cds_intervals = gene_info.get('cds', [])
    print(f"\n[INFO] CDS区间:")
    for cds in cds_intervals:
        print(f"  - {cds}")
    
    # 检查有多少变异位于CDS区域内
    cds_vars = []
    for _, row in df.iterrows():
        pos = row['position']
        for cds_start, cds_end in cds_intervals:
            if cds_start <= pos <= cds_end:
                cds_vars.append(row)
                break
    
    print(f"\n[INFO] 位于CDS区域的变异: {len(cds_vars)}")
    if cds_vars:
        cds_df = pd.DataFrame(cds_vars)
        print(f"  - annotation分布:")
        for ann, count in cds_df['annotation'].value_counts().items():
            print(f"    {ann}: {count}")
    
    return True


def test_annotate_function_directly():
    """直接测试annotate_snp_effects_for_region函数"""
    
    from haplotype_phenotype_analysis import annotate_snp_effects_for_region, _build_coding_context
    
    test_gene = 'HORVU.MOREX.r3.2HG0213430'
    test_db_dir = 'd:/Desktop/project1/database_text'
    
    import json
    import pandas as pd
    
    print("\n" + "=" * 60)
    print("直接测试 annotate_snp_effects_for_region 函数")
    print("=" * 60)
    
    gene_db_dir = os.path.join(test_db_dir, test_gene)
    
    # 读取gene_info
    with open(os.path.join(gene_db_dir, 'gene_info.json'), 'r') as f:
        gene_info = json.load(f)
    
    # 读取variant_info
    df = pd.read_csv(os.path.join(gene_db_dir, 'variant_info.csv'))
    
    chrom = gene_info['chrom']
    gene_start = gene_info['gene_start']
    gene_end = gene_info['gene_end']
    strand = gene_info['strand']
    cds_intervals = gene_info['cds']
    exon_intervals = gene_info['exons']
    positions = df['position'].tolist()
    
    print(f"\n[INFO] 基因信息:")
    print(f"  - 染色体: {chrom}")
    print(f"  - 基因位置: {gene_start}-{gene_end}")
    print(f"  - CDS区间: {cds_intervals}")
    print(f"  - 待注释位点数: {len(positions)}")
    
    # 计算promoter区间
    promoter_length = gene_info.get('promoter_length', 2000)
    if strand == '+':
        promoter_start = max(1, gene_start - promoter_length)
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + promoter_length
    
    print(f"  - promoter区间: {promoter_start}-{promoter_end}")
    
    # 直接调用annotate_snp_effects_for_region（不使用FASTA）
    print(f"\n[TEST] 调用 annotate_snp_effects_for_region (无FASTA)...")
    effects = annotate_snp_effects_for_region(
        vcf_file=None,
        fasta_path=None,
        gene_chrom=chrom,
        cds_intervals=cds_intervals,
        exon_intervals=exon_intervals,
        gene_strand=strand,
        positions=positions,
        gene_start=gene_start,
        gene_end=gene_end,
        promoter_start=promoter_start,
        promoter_end=promoter_end
    )
    
    print(f"\n[INFO] 注释结果分布:")
    effect_counts = {}
    for eff in effects.values():
        effect_counts[eff] = effect_counts.get(eff, 0) + 1
    for eff, count in sorted(effect_counts.items(), key=lambda x: -x[1]):
        print(f"  - {eff}: {count}")
    
    return True


if __name__ == '__main__':
    test_database_variant_info()
    test_annotate_function_directly()
    print("\n" + "=" * 60)
    print("测试完成")
    print("=" * 60)