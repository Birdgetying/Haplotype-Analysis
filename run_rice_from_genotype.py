#!/usr/bin/env python3
"""
水稻单倍型分析 - 从tab-separated基因型文件直接构建
基因型文件包含所有4726样本的填充基因型，与表型样本(C001-C202, W001-W327)匹配
"""

import sys, os, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, r'd:\Desktop\project1')

import pandas as pd
import numpy as np
import gzip, csv, json, re, time
from datetime import datetime
from collections import Counter

# ========== 配置 ==========
GENO_DIR = r'D:\Desktop\data\水稻\Tab-separated genotype files'
PHENO_PATH = r'D:\Desktop\data\水稻\Phenotypes\phenos.csv'
H5_DIR = r'D:\Desktop\data\水稻\Results of variant effect annotation'
OUT_DB = r'D:\Desktop\project1\rice_database2'

# 3个测试基因 (chr01)
GENES = [
    {'gene_id': 'LOC_Os01g02660', 'chrom': 'chr01', 'start': 900450, 'end': 906961, 'strand': '+'},
    {'gene_id': 'LOC_Os01g25140', 'chrom': 'chr01', 'start': 14188041, 'end': 14193462, 'strand': '+'},
    {'gene_id': 'LOC_Os01g57340', 'chrom': 'chr01', 'start': 33143983, 'end': 33150573, 'strand': '-'},
]

PROMOTER_LENGTH = 2000
MIN_HAP_COUNT = 2  # 单倍型最少样本数，低于此数量的不参与打分

from rice_common import SNPEFF_TO_CATEGORY, load_h5_annotations, classify_variant


def load_phenotype():
    """加载并标准化表型数据"""
    df = pd.read_csv(PHENO_PATH)
    # 使用 id_name 作为 SampleID
    df = df.drop(columns=['id']).rename(columns={'id_name': 'SampleID'})
    return df


def build_haplotype_for_gene(gene_info, pheno_df):
    """直接从tab-separated基因型文件构建单倍型数据"""
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    strand = gene_info['strand']

    # 计算启动子扩展区域
    gene_start = gene_info['start']
    gene_end = gene_info['end']
    if strand == '+':
        region_start = max(1, gene_start - PROMOTER_LENGTH)
        region_end = gene_end
    else:
        region_start = gene_start
        region_end = gene_end + PROMOTER_LENGTH

    # 确定染色体编号
    chr_num = chrom.replace('chr', '')
    chr_idx = int(chr_num)

    geno_file = os.path.join(GENO_DIR, f'rice4k_chr{chr_num}_geno.txt.gz')
    if not os.path.exists(geno_file):
        print(f"[ERROR] Genotype file not found: {geno_file}")
        return None

    print(f"[INFO] {gene_id}: 区域 {chrom}:{region_start}-{region_end}")

    # 读取表型样本
    pheno_samples = set(pheno_df['SampleID'].tolist())

    # 扫描基因型文件，提取目标区域的变异
    print(f"[INFO] 扫描基因型文件 {os.path.basename(geno_file)}...")
    t0 = time.time()

    variant_data = {}  # {var_id: {sample: allele}}
    header_samples = []

    with gzip.open(geno_file, 'rt', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        header_samples = header[1:]  # Skip first column

        for row in reader:
            var_id = row[0]
            # Parse position from var ID: e.g., "0100001151" -> pos=1151
            # Format: chr(2 digits) + position(8 digits)
            try:
                v_chr = int(var_id[:2])
                v_pos = int(var_id[2:])
            except ValueError:
                continue

            if v_chr != chr_idx:
                if v_chr > chr_idx:
                    break  # 已经超出目标染色体
                continue

            if v_pos < region_start:
                continue
            if v_pos > region_end:
                break  # 已超出目标区域

            # 收集等位基因
            genotypes = row[1:]
            allele_dict = {}
            for sample, allele in zip(header_samples, genotypes):
                if sample in pheno_samples:
                    allele = allele.strip('"').strip()
                    if allele and allele != 'N' and allele != '':
                        allele_dict[sample] = allele.upper()

            if allele_dict:
                variant_data[var_id] = allele_dict

    print(f"[INFO] 找到 {len(variant_data)} 个变异的基因型数据 (耗时 {time.time()-t0:.1f}s)")

    if not variant_data:
        print(f"[WARNING] 目标区域无变异数据")
        return None

    # 过滤：只保留有 >1 种等位基因的位点
    filtered_variants = {}
    for var_id, alleles in variant_data.items():
        unique_alleles = set(alleles.values())
        if len(unique_alleles) > 1:
            filtered_variants[var_id] = alleles

    print(f"[INFO] 过滤后有 {len(filtered_variants)} 个多态性位点")

    if not filtered_variants:
        return None

    # 确定样本集（所有有基因型数据的表型样本）
    all_typed_samples = set()
    for alleles in filtered_variants.values():
        all_typed_samples.update(alleles.keys())
    all_typed_samples = sorted(all_typed_samples)
    print(f"[INFO] 有基因型数据的表型样本: {len(all_typed_samples)}")

    # 构建位置列表（按位置排序）
    var_ids = sorted(filtered_variants.keys(), key=lambda x: int(x[2:]))
    positions = [int(v[2:]) for v in var_ids]

    # 构建单倍型序列
    sample_alleles = {s: [] for s in all_typed_samples}
    for var_id in var_ids:
        alleles = filtered_variants[var_id]
        for sample in all_typed_samples:
            sample_alleles[sample].append(alleles.get(sample, 'N'))

    # 过滤掉有N的样本
    valid_samples = []
    for sample, alleles in sample_alleles.items():
        if 'N' not in alleles:
            valid_samples.append(sample)

    print(f"[INFO] 完整基因型样本: {len(valid_samples)}")

    if len(valid_samples) < 2:
        print(f"[WARNING] 有效样本过少")
        return None

    # 构建单倍型
    hap_seqs = {}
    for sample in valid_samples:
        seq = '|'.join(sample_alleles[sample])
        hap_seqs[sample] = seq

    # 命名单倍型（过滤数量极少的）
    seq_counts = Counter(hap_seqs.values())
    valid_seqs = {seq for seq, count in seq_counts.items() if count >= MIN_HAP_COUNT}
    n_filtered = len(seq_counts) - len(valid_seqs)
    if n_filtered > 0:
        print(f"[INFO] {gene_id}: 过滤 {n_filtered} 个低样本数单倍型 (min_count={MIN_HAP_COUNT})")

    sorted_seqs = seq_counts.most_common()
    seq_to_name = {}
    for i, (seq, count) in enumerate(sorted_seqs):
        if seq not in valid_seqs:
            continue
        seq_to_name[seq] = f'Hap{i+1}'

    # 构建 hap_sample_df
    sample_haps = []
    filtered_samples = []
    for sample in valid_samples:
        hap_seq = hap_seqs[sample]
        if hap_seq not in valid_seqs:
            filtered_samples.append(sample)
            continue
        sample_haps.append({
            'SampleID': sample,
            'Haplotype_Seq': hap_seq,
            'Hap_Name': seq_to_name.get(hap_seq, 'Other')
        })
    hap_sample_df = pd.DataFrame(sample_haps)
    if filtered_samples:
        print(f"[INFO] {gene_id}: 因单倍型过滤跳过的样本: {len(filtered_samples)}")

    # 构建 hap_df
    hap_list = []
    for seq, name in seq_to_name.items():
        alleles = seq.split('|')
        hap_list.append({
            'Haplotype_Seq': seq,
            'Count': seq_counts[seq],
            'Hap_Name': name,
            'Alleles': alleles
        })
    hap_df = pd.DataFrame(hap_list)

    # 计算启动子区域（用于注释分类）
    if strand == '+':
        promoter_start = max(1, gene_start - PROMOTER_LENGTH)
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + PROMOTER_LENGTH

    # 加载H5注释
    h5_annotations = load_h5_annotations(H5_DIR, chrom, positions)

    # 构建 variant_info
    variant_info = {}
    for i, var_id in enumerate(var_ids):
        pos = positions[i]
        alleles_at_pos = [sample_alleles[s][i] for s in valid_samples]
        allele_counts = Counter(alleles_at_pos)
        # 按频率降序：最常见等位基因视为ref，第二常见视为alt
        sorted_by_freq = sorted(allele_counts.keys(), key=lambda a: allele_counts[a], reverse=True)
        sorted_counts = sorted(allele_counts.values(), reverse=True)
        maf = sorted_counts[1] / len(alleles_at_pos) if len(sorted_counts) >= 2 else 0.5

        ref = sorted_by_freq[0] if sorted_by_freq else ''
        alt = sorted_by_freq[1] if len(sorted_by_freq) > 1 else ''
        len_diff = len(alt) - len(ref) if alt else 0
        is_sv = abs(len_diff) >= 50

        h5_anno = h5_annotations.get(pos)
        annotation = classify_variant(
            pos, ref, alt, gene_start, gene_end,
            promoter_start, promoter_end, h5_anno=h5_anno, len_diff=len_diff
        )

        variant_info[pos] = {
            'ref': ref,
            'alt': alt,
            'len_diff': len_diff,
            'is_sv': is_sv,
            'maf': maf,
            'missing_rate': 0.0,
            'annotation': annotation
        }

    return {
        'gene_id': gene_id,
        'chrom': chrom,
        'start': region_start,
        'end': region_end,
        'gene_start': gene_start,
        'gene_end': gene_end,
        'strand': strand,
        'positions': positions,
        'hap_df': hap_df,
        'hap_sample_df': hap_sample_df,
        'variant_info': variant_info,
        'n_variants': len(positions),
        'n_haplotypes': len(hap_df),
        'n_samples': len(valid_samples),
    }


def save_database(data, out_dir):
    """保存数据库文件"""
    gene_id = data['gene_id']
    gene_dir = os.path.join(out_dir, gene_id)
    os.makedirs(gene_dir, exist_ok=True)

    # gene_info.json
    gene_info = {
        'gene_id': gene_id,
        'chrom': data['chrom'],
        'start': data['start'],
        'end': data['end'],
        'gene_start': data['gene_start'],
        'gene_end': data['gene_end'],
        'strand': data['strand'],
        'length': data['end'] - data['start'] + 1,
        'promoter_length': PROMOTER_LENGTH,
        'exons': [],
        'cds': []
    }
    with open(os.path.join(gene_dir, 'gene_info.json'), 'w') as f:
        json.dump(gene_info, f, indent=2)

    # haplotype_data.csv
    data['hap_df'].to_csv(os.path.join(gene_dir, 'haplotype_data.csv'), index=False)

    # haplotype_samples.csv
    data['hap_sample_df'].to_csv(os.path.join(gene_dir, 'haplotype_samples.csv'), index=False)

    # variant_info.csv
    vi = data['variant_info']
    variant_df = pd.DataFrame([
        {
            'position': pos,
            'ref': info.get('ref', ''),
            'alt': info.get('alt', ''),
            'len_diff': info.get('len_diff', 0),
            'is_sv': info.get('is_sv', False),
            'maf': info.get('maf', 0.5),
            'missing_rate': info.get('missing_rate', 0.0),
            'annotation': info.get('annotation', 'other')
        }
        for pos, info in vi.items()
    ])
    variant_df.to_csv(os.path.join(gene_dir, 'variant_info.csv'), index=False)

    print(f"[INFO] 数据库已保存: {gene_dir}")


def merge_phenotype(data, pheno_df):
    """合并表型数据"""
    merged = pd.merge(data['hap_sample_df'], pheno_df, on='SampleID', how='inner')
    merged.to_csv(os.path.join(OUT_DB, data['gene_id'], 'phenotype_data.csv'), index=False)

    # 统计单倍型
    hap_col = 'Hap_Name'
    hap_stats_list = []
    for hap_name in merged[hap_col].unique():
        hap_data = merged[merged[hap_col] == hap_name]
        pheno_col = 'Heading_date'  # 第一个表型
        pheno_vals = hap_data[pheno_col].dropna()
        hap_stats_list.append({
            'haplotype_name': hap_name,
            'haplotype_count': len(hap_data),
            'haplotype_freq': round(len(hap_data) / len(merged), 4),
            'phenotype_mean': round(pheno_vals.mean(), 4) if len(pheno_vals) > 0 else None,
            'phenotype_sd': round(pheno_vals.std(), 4) if len(pheno_vals) > 1 else None,
        })

    hap_stats_df = pd.DataFrame(hap_stats_list)
    hap_stats_df.to_csv(os.path.join(OUT_DB, data['gene_id'], 'haplotype_stats.csv'), index=False)

    return merged


def run_anova_analysis(merged, gene_id):
    """ANOVA关联分析"""
    from scipy import stats

    hap_col = 'Hap_Name'
    pheno_cols = [c for c in merged.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]

    results = {}
    for pheno_col in pheno_cols:
        analysis_df = merged[merged[hap_col] != 'Other'].copy()
        if len(analysis_df) < 3:
            continue

        groups = [analysis_df[analysis_df[hap_col] == h][pheno_col].dropna().values
                  for h in analysis_df[hap_col].unique()
                  if len(analysis_df[analysis_df[hap_col] == h]) >= 2]

        if len(groups) < 2:
            continue

        try:
            f_stat, p_value = stats.f_oneway(*groups)
            grand_mean = analysis_df[pheno_col].mean()
            ss_between = sum(len(analysis_df[analysis_df[hap_col] == h]) *
                            (analysis_df[analysis_df[hap_col] == h][pheno_col].mean() - grand_mean)**2
                            for h in analysis_df[hap_col].unique())
            ss_total = ((analysis_df[pheno_col] - grand_mean)**2).sum()
            pve = (ss_between / ss_total) * 100 if ss_total > 0 else 0

            results[pheno_col] = {
                'f_statistic': f_stat, 'p_value': p_value,
                'pve_percent': pve,
                'significant': p_value < 0.05
            }
        except Exception as e:
            print(f"  [WARNING] ANOVA failed for {pheno_col}: {e}")

    return results


# ========== 主流程 ==========
if __name__ == '__main__':
    print("=" * 60)
    print("水稻单倍型分析 - 从基因型文件构建")
    print("=" * 60)

    os.makedirs(OUT_DB, exist_ok=True)

    # 1. 加载表型
    print("\n[1] 加载表型数据...")
    pheno_df = load_phenotype()
    print(f"  表型样本: {len(pheno_df)}")
    print(f"  表型列: {[c for c in pheno_df.columns if c != 'SampleID']}")

    # 2. 处理每个基因
    all_results = []
    for gene_info in GENES:
        gene_id = gene_info['gene_id']
        print(f"\n{'='*60}")
        print(f"[2] 处理基因: {gene_id}")
        print(f"{'='*60}")

        # 构建单倍型数据
        data = build_haplotype_for_gene(gene_info, pheno_df)
        if data is None:
            print(f"[SKIP] {gene_id}: 无有效数据")
            continue

        # 保存数据库
        save_database(data, OUT_DB)

        # 合并表型
        merged = merge_phenotype(data, pheno_df)
        print(f"[INFO] 表型合并: {len(merged)} 个样本匹配")
        print(f"  单倍型分布: {dict(merged['Hap_Name'].value_counts())}")

        # ANOVA分析
        print(f"\n[3] 关联分析...")
        anova_results = run_anova_analysis(merged, gene_id)

        for pheno_col, res in anova_results.items():
            sig_mark = '***' if res['p_value'] < 0.001 else ('**' if res['p_value'] < 0.01 else ('*' if res['p_value'] < 0.05 else ''))
            print(f"  {pheno_col}: P={res['p_value']:.4f} {sig_mark}  PVE={res['pve_percent']:.2f}%")

        # 保存关联结果
        if anova_results:
            best_pheno = min(anova_results.items(), key=lambda x: x[1]['p_value'])
            result_summary = {
                'gene_id': gene_id,
                'chrom': data['chrom'],
                'start': data['start'],
                'end': data['end'],
                'strand': data['strand'],
                'n_variants': data['n_variants'],
                'n_haplotypes': data['n_haplotypes'],
                'n_samples': data['n_samples'],
                'best_phenotype': best_pheno[0],
                'best_pvalue': best_pheno[1]['p_value'],
                'best_pve': best_pheno[1]['pve_percent'],
                'significant': best_pheno[1]['significant'],
            }
            all_results.append(result_summary)

            # 保存
            assoc_df = pd.DataFrame([{
                'gene_id': gene_id,
                'phenotype': pheno_col,
                'f_statistic': res['f_statistic'],
                'p_value': res['p_value'],
                'pve_percent': res['pve_percent'],
                'significant': res['significant'],
            } for pheno_col, res in anova_results.items()])
            assoc_df.to_csv(os.path.join(OUT_DB, gene_id, 'association_result.csv'), index=False)

    # 4. 汇总
    print(f"\n{'='*60}")
    print(f"[4] 分析汇总")
    print(f"{'='*60}")
    if all_results:
        summary_df = pd.DataFrame(all_results)
        summary_df.to_csv(os.path.join(OUT_DB, 'summary.csv'), index=False)
        print(summary_df.to_string())

        # 找出显著基因
        sig = summary_df[summary_df['significant']]
        print(f"\n显著基因 ({len(sig)}/{len(summary_df)}):")
        for _, r in sig.iterrows():
            print(f"  {r['gene_id']}: P={r['best_pvalue']:.6f}, PVE={r['best_pve']:.2f}%, phenotype={r['best_phenotype']}")
    else:
        print("无有效分析结果")

    print(f"\n数据库保存于: {OUT_DB}")
    print("Done!")
