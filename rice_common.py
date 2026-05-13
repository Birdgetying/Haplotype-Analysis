#!/usr/bin/env python3
"""水稻分析共享工具：SnpEff注释映射、H5注释加载、变异分类"""

import os
import pandas as pd

# SnpEff注释 → 内部分类映射
SNPEFF_TO_CATEGORY = {
    'missense_variant': 'missense',
    'stop_gained': 'missense',
    'stop_lost': 'missense',
    'start_lost': 'missense',
    'synonymous_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',
    'intron_variant': 'intron',
    '3_prime_UTR_variant': 'UTR',
    '5_prime_UTR_variant': 'UTR',
    'frameshift_variant': 'frameshift',
    'inframe_deletion': 'inframe_del',
    'inframe_insertion': 'inframe_ins',
    'disruptive_inframe_deletion': 'inframe_del',
    'disruptive_inframe_insertion': 'inframe_ins',
    'splice_donor_variant': 'splice_site',
    'splice_acceptor_variant': 'splice_site',
    'splice_region_variant': 'splice_region',
}

IMPACT_PRIORITY = {'HIGH': 5, 'MODERATE': 4, 'LOW': 3, 'MODIFIER': 2}


def load_h5_annotations(h5_dir, chrom, target_positions):
    """从H5文件加载指定位置的SnpEff注释

    Returns:
        dict: {pos: snpeff_anno} — 每个位置取最高impact的注释
    """
    chr_num = chrom.replace('chr', '')
    h5_file = os.path.join(h5_dir,
                           f'chr{chr_num}_snpeff_coovar_polyphen_sift_merge_anno.h5')
    if not os.path.exists(h5_file):
        print(f"[INFO] H5注释文件不存在: {h5_file}，使用位置回退")
        return {}

    try:
        with pd.HDFStore(h5_file, 'r') as store:
            key = f'/chr{chr_num}'
            if key not in store.keys():
                print(f"[WARN] H5中找不到key: {key}")
                return {}

            target_set = set(target_positions)
            # 只读3个需要的列（var + snpeff_anno + snpeff_impact），避免加载20列×6.5M行
            df = store.select(key, columns=['var', 'snpeff_anno', 'snpeff_impact'])
            positions_series = df['var'].str[4:].astype(int)
            mask = positions_series.isin(target_set)

            if not mask.any():
                return {}

            df_f = df.loc[mask].copy()
            df_f['pos'] = positions_series[mask].values
            df_f['impact_score'] = (
                df_f['snpeff_impact'].map(IMPACT_PRIORITY).fillna(0))

            idx = df_f.groupby('pos')['impact_score'].idxmax()
            best = df_f.loc[idx]

            result = dict(zip(best['pos'], best['snpeff_anno']))
            print(f"[INFO] 从H5加载了 {len(result)}/{len(target_positions)} 个位点的注释")
            return result
    except Exception as e:
        print(f"[WARN] 加载H5注释失败: {e}")
        return {}


def classify_variant(pos, ref, alt, gene_start, gene_end,
                     promoter_start, promoter_end, h5_anno=None, len_diff=None):
    """综合SnpEff注释、len_diff和基因组位置确定变异分类"""
    if len_diff is None:
        len_diff = len(alt) - len(ref) if alt else 0

    if h5_anno and h5_anno != 'N':
        for snpeff_key, category in SNPEFF_TO_CATEGORY.items():
            if snpeff_key in h5_anno:
                return category

        if 'upstream_gene_variant' in h5_anno or 'downstream_gene_variant' in h5_anno:
            if promoter_start <= pos <= promoter_end:
                return 'promoter'
            if gene_start <= pos <= gene_end:
                return 'intron'
            return 'other'

        if 'intergenic_region' in h5_anno:
            return 'other'

    if abs(len_diff) >= 50:
        return 'SV'
    if len_diff > 0:
        return 'INS'
    if len_diff < 0:
        return 'DEL'

    if promoter_start <= pos <= promoter_end:
        return 'promoter'
    if gene_start <= pos <= gene_end:
        return 'intron'
    return 'other'
