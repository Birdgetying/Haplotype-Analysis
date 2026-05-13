#!/usr/bin/env python3
"""
水稻单倍型分析 - 谢为博2021论文关键基因
从tab-separated基因型文件构建数据库，生成HTML报告
"""
import sys, os, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, r'd:\Desktop\project1')

import pandas as pd
import numpy as np
import gzip, csv, json, time
from collections import Counter

from rice_common import SNPEFF_TO_CATEGORY, classify_variant

# ========== 配置 ==========
GENO_DIR = r'D:\Desktop\data\水稻\Tab-separated genotype files'
PHENO_PATH = r'D:\Desktop\data\水稻\Phenotypes\phenos.csv'
H5_DIR = r'D:\Desktop\data\水稻\Results of variant effect annotation'
OUT_DB = r'D:\Desktop\project1\rice_database'
OUT_RESULTS = r'D:\Desktop\project1\rice_results'

PROMOTER_LENGTH = 2000

GENES = [
    {'gene_id': 'DEP1', 'msu_id': 'LOC_Os09g26999', 'chrom': 'chr09',
     'start': 16406348, 'end': 16420832, 'strand': '+'},
    {'gene_id': 'OsMADS25', 'msu_id': 'LOC_Os04g23910', 'chrom': 'chr04',
     'start': 13667728, 'end': 13680873, 'strand': '-'},
    {'gene_id': 'GW7', 'msu_id': 'LOC_Os07g41200', 'chrom': 'chr07',
     'start': 24659182, 'end': 24674279, 'strand': '+'},
]


def load_phenotype():
    df = pd.read_csv(PHENO_PATH)
    df = df.drop(columns=['id']).rename(columns={'id_name': 'SampleID'})
    return df


def get_gene_info_from_h5(gene_info):
    """从H5文件中提取基因的CDS/外显子信息"""
    msu_id = gene_info['msu_id']
    chr_id = gene_info['chrom'].replace('chr', '')

    h5_file = os.path.join(H5_DIR, f'chr{chr_id}_snpeff_coovar_polyphen_sift_merge_anno.h5')
    if not os.path.exists(h5_file):
        return None

    df = pd.read_hdf(h5_file, key=f'/chr{chr_id}')
    gene_df = df[df['gene'] == msu_id]
    del df

    if len(gene_df) == 0:
        return None

    # Parse positions and annotations
    positions = []
    annotations = {}
    for _, row in gene_df.iterrows():
        try:
            pos = int(str(row['var'])[4:])  # vgXXPPPPPPPP -> position
            positions.append(pos)
            impact = row.get('snpeff_impact', 'MODIFIER')
            anno = row.get('snpeff_anno', '')
            annotations[pos] = {'impact': impact, 'anno': anno}
        except:
            pass

    if not positions:
        return None

    positions = sorted(set(positions))
    gene_start = gene_info['start']
    gene_end = gene_info['end']

    # Classify positions into CDS/exons based on impact
    cds_positions = [p for p in positions
                     if annotations.get(p, {}).get('impact') in ('MODERATE', 'HIGH')]
    utr_positions = [p for p in positions
                     if annotations.get(p, {}).get('impact') == 'LOW']

    # Build CDS intervals
    cds = []
    if cds_positions:
        cds_positions = sorted(set(cds_positions))
        cds_start = cds_positions[0]
        cds_end = cds_positions[-1]
        cds.append([cds_start, cds_end])

    # Build exon intervals (CDS + UTR)
    exons = []
    all_func_positions = sorted(set(cds_positions + utr_positions))
    if all_func_positions:
        exons.append([all_func_positions[0], all_func_positions[-1]])
    elif cds:
        exons = cds[:]

    return {
        'positions': positions,
        'exons': exons,
        'cds': cds,
        'annotations': annotations,
    }


def build_haplotype_for_gene(gene_info, pheno_df):
    """从基因型文件构建单倍型数据"""
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    strand = gene_info['strand']
    gene_start = gene_info['start']
    gene_end = gene_info['end']

    if strand == '+':
        region_start = max(1, gene_start - PROMOTER_LENGTH)
        region_end = gene_end
    else:
        region_start = gene_start
        region_end = gene_end + PROMOTER_LENGTH

    chr_str = chrom.replace('chr', '')  # keep leading zero: '09', '04', '07'
    geno_file = os.path.join(GENO_DIR, f'rice4k_chr{chr_str}_geno.txt.gz')

    print(f"[INFO] {gene_id}: region {chrom}:{region_start}-{region_end}")
    print(f"[INFO] Scanning {os.path.basename(geno_file)}...")
    t0 = time.time()

    pheno_samples = set(pheno_df['SampleID'].tolist())
    variant_data = {}
    header_samples = []

    with gzip.open(geno_file, 'rt', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        header_samples = header[1:]

        for row in reader:
            var_id = row[0]
            try:
                v_chr = int(var_id[:2])
                v_pos = int(var_id[2:])
            except ValueError:
                continue

            if v_chr > int(chr_str):
                break
            if v_chr < int(chr_str):
                continue
            if v_pos < region_start:
                continue
            if v_pos > region_end:
                break

            genotypes = row[1:]
            allele_dict = {}
            for sample, allele in zip(header_samples, genotypes):
                if sample in pheno_samples:
                    allele = allele.strip('"').strip()
                    if allele and allele != 'N' and allele != '':
                        allele_dict[sample] = allele.upper()

            if allele_dict:
                variant_data[var_id] = allele_dict

    print(f"[INFO] {len(variant_data)} variants found ({time.time()-t0:.1f}s)")

    if not variant_data:
        return None

    # Filter to polymorphic positions
    filtered_variants = {}
    for var_id, alleles in variant_data.items():
        if len(set(alleles.values())) > 1:
            filtered_variants[var_id] = alleles

    print(f"[INFO] {len(filtered_variants)} polymorphic positions")

    if not filtered_variants:
        return None

    # Get all typed samples
    all_typed = set()
    for alleles in filtered_variants.values():
        all_typed.update(alleles.keys())
    all_typed = sorted(all_typed)
    print(f"[INFO] {len(all_typed)} samples with genotype data")

    # Build position-sorted sequences
    var_ids = sorted(filtered_variants.keys(), key=lambda x: int(x[2:]))
    positions = [int(v[2:]) for v in var_ids]

    sample_alleles = {s: [] for s in all_typed}
    for var_id in var_ids:
        alleles = filtered_variants[var_id]
        for sample in all_typed:
            sample_alleles[sample].append(alleles.get(sample, 'N'))

    # Filter samples with complete data
    valid_samples = [s for s in all_typed if 'N' not in sample_alleles[s]]
    print(f"[INFO] {len(valid_samples)} complete samples")

    if len(valid_samples) < 2:
        return None

    # Build haplotypes
    hap_seqs = {}
    for sample in valid_samples:
        hap_seqs[sample] = '|'.join(sample_alleles[sample])

    seq_counts = Counter(hap_seqs.values())
    seq_to_name = {}
    for i, (seq, count) in enumerate(seq_counts.most_common()):
        seq_to_name[seq] = f'Hap{i+1}'

    sample_haps = []
    for sample in valid_samples:
        seq = hap_seqs[sample]
        sample_haps.append({
            'SampleID': sample,
            'Haplotype_Seq': seq,
            'Hap_Name': seq_to_name.get(seq, 'Other')
        })
    hap_sample_df = pd.DataFrame(sample_haps)

    hap_list = []
    for seq, name in seq_to_name.items():
        hap_list.append({
            'Haplotype_Seq': seq,
            'Count': seq_counts[seq],
            'Hap_Name': name,
            'Alleles': seq.split('|')
        })
    hap_df = pd.DataFrame(hap_list)

    # Build variant_info
    variant_info = {}
    # 计算启动子区域
    if strand == '+':
        promoter_start = max(1, gene_start - PROMOTER_LENGTH)
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + PROMOTER_LENGTH

    for i, var_id in enumerate(var_ids):
        pos = int(var_id[2:])
        alleles_at_pos = [sample_alleles[s][i] for s in valid_samples]
        allele_counts = Counter(alleles_at_pos)
        # 按频率降序：最常见=ref，第二常见=alt
        sorted_by_freq = sorted(allele_counts.keys(), key=lambda a: allele_counts[a], reverse=True)
        sorted_counts = sorted(allele_counts.values(), reverse=True)
        maf = sorted_counts[1] / len(alleles_at_pos) if len(sorted_counts) >= 2 else 0.5

        ref = sorted_by_freq[0] if sorted_by_freq else ''
        alt = sorted_by_freq[1] if len(sorted_by_freq) > 1 else ''
        len_diff = len(alt) - len(ref) if alt else 0
        is_sv = abs(len_diff) >= 50

        variant_info[pos] = {
            'ref': ref,
            'alt': alt,
            'len_diff': len_diff,
            'is_sv': is_sv,
            'maf': maf,
            'missing_rate': 0.0,
            'promoter_start': promoter_start,
            'promoter_end': promoter_end,
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


def save_database(data, out_dir, h5_info=None):
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
        'exons': h5_info.get('exons', []) if h5_info else [],
        'cds': h5_info.get('cds', []) if h5_info else [],
    }
    with open(os.path.join(gene_dir, 'gene_info.json'), 'w', encoding='utf-8') as f:
        json.dump(gene_info, f, indent=2, ensure_ascii=False)

    data['hap_df'].to_csv(os.path.join(gene_dir, 'haplotype_data.csv'), index=False)
    data['hap_sample_df'].to_csv(os.path.join(gene_dir, 'haplotype_samples.csv'), index=False)

    vi = data['variant_info']
    # 从h5_info获取SnpEff注释
    h5_annotations = h5_info.get('annotations', {}) if h5_info else {}
    gene_start = data.get('gene_start', 0)
    gene_end = data.get('gene_end', 0)
    strand = data.get('strand', '+')

    variant_rows = []
    for pos, info in vi.items():
        ref = info.get('ref', '')
        alt = info.get('alt', '')
        len_diff = info.get('len_diff', 0)
        is_sv = info.get('is_sv', False)
        promoter_start = info.get('promoter_start', 0)
        promoter_end = info.get('promoter_end', 0)

        # 从H5获取SnpEff注释
        h5_entry = h5_annotations.get(pos, {})
        snpeff_anno = h5_entry.get('anno', '') if isinstance(h5_entry, dict) else ''

        annotation = classify_variant(
            pos, ref, alt, gene_start, gene_end,
            promoter_start, promoter_end, h5_anno=snpeff_anno, len_diff=len_diff
        )

        variant_rows.append({
            'position': pos,
            'ref': ref,
            'alt': alt,
            'len_diff': len_diff,
            'is_sv': is_sv,
            'maf': info.get('maf', 0.5),
            'missing_rate': info.get('missing_rate', 0.0),
            'annotation': annotation
        })
    variant_df = pd.DataFrame(variant_rows)
    variant_df.to_csv(os.path.join(gene_dir, 'variant_info.csv'), index=False)

    # 识别启动子区域变异并保存（供后续分析使用）
    strand = data['strand']
    gene_start = data['gene_start']
    gene_end = data['gene_end']
    if strand == '+':
        promoter_start = max(1, gene_start - PROMOTER_LENGTH)
        promoter_end = gene_start - 1
    else:
        promoter_start = gene_end + 1
        promoter_end = gene_end + PROMOTER_LENGTH

    promoter_variants = []
    for pos, info in vi.items():
        if promoter_start <= pos <= promoter_end:
            distance_to_tss = abs(pos - gene_start) if strand == '+' else abs(pos - gene_end)
            promoter_variants.append({
                'position': pos,
                'ref': info.get('ref', ''),
                'alt': info.get('alt', ''),
                'distance_to_tss': distance_to_tss,
                'is_sv': info.get('is_sv', False),
                'maf': info.get('maf', 0.5),
                'overlaps_cds': False,
                'overlapping_genes': ''
            })

    if promoter_variants:
        promoter_df = pd.DataFrame(promoter_variants)
        promoter_df.to_csv(os.path.join(gene_dir, 'promoter_variants.csv'), index=False)
        # 详细文本报告
        with open(os.path.join(gene_dir, 'promoter_variants_detail.txt'), 'w', encoding='utf-8') as f:
            f.write(f"启动子区域变异分析报告\n")
            f.write(f"{'=' * 60}\n")
            f.write(f"基因: {gene_id}\n")
            f.write(f"染色体: {data['chrom']}\n")
            f.write(f"启动子区域: {promoter_start:,}-{promoter_end:,}\n")
            f.write(f"变异总数: {len(promoter_variants)}\n\n")
            for v in promoter_variants:
                sv_marker = " [SV]" if v['is_sv'] else ""
                f.write(f"  位置: {v['position']:,} | 距离TSS: {v['distance_to_tss']:,}bp | "
                       f"REF: {v['ref']} | ALT: {v['alt']} | MAF: {v['maf']:.3f}{sv_marker}\n")
        print(f"[INFO] 启动子变异: {len(promoter_variants)} 个 -> promoter_variants.csv")
    else:
        print(f"[INFO] 启动子区域未发现变异")

    print(f"[INFO] Database saved: {gene_dir}")


def merge_phenotype(data, pheno_df, out_dir):
    gene_id = data['gene_id']
    merged = pd.merge(data['hap_sample_df'], pheno_df, on='SampleID', how='inner')
    merged.to_csv(os.path.join(out_dir, gene_id, 'phenotype_data.csv'), index=False)

    hap_col = 'Hap_Name'
    pheno_cols = [c for c in merged.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
    hap_stats = []
    for hap_name in merged[hap_col].unique():
        hap_data = merged[merged[hap_col] == hap_name]
        pheno = pheno_cols[0]
        vals = hap_data[pheno].dropna()
        hap_stats.append({
            'haplotype_name': hap_name,
            'haplotype_count': len(hap_data),
            'haplotype_freq': round(len(hap_data) / len(merged), 4),
            'phenotype_mean': round(vals.mean(), 4) if len(vals) > 0 else None,
            'phenotype_sd': round(vals.std(), 4) if len(vals) > 1 else None,
        })
    pd.DataFrame(hap_stats).to_csv(os.path.join(out_dir, gene_id, 'haplotype_stats.csv'), index=False)

    return merged, pheno_cols


def generate_html(gene_id, gene_info_dict):
    """调用 HaplotypePhenotypeAnalyzer 生成 HTML"""
    from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

    db_dir = OUT_DB
    results_dir = os.path.join(OUT_RESULTS, gene_id)
    os.makedirs(results_dir, exist_ok=True)

    gene_dir = os.path.join(db_dir, gene_id)
    pheno_data = pd.read_csv(os.path.join(gene_dir, 'phenotype_data.csv'))
    pheno_cols = [c for c in pheno_data.columns
                  if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]

    vcf_file = r'D:\Desktop\data\水稻\VCF format\rice4k_geno_add_del.vcf.gz'
    gff_file = r'D:\Desktop\data\水稻\rice_test_genes_paper.gff3'

    analyzer = HaplotypePhenotypeAnalyzer(
        vcf_file=vcf_file if os.path.exists(vcf_file) else None,
        phenotype_file=None,
        output_dir=results_dir,
        gtf_file=gff_file if os.path.exists(gff_file) else None
    )
    analyzer.phenotype_df = pheno_data

    try:
        result = analyzer.analyze_gene(
            chrom=gene_info_dict['chrom'],
            start=gene_info_dict['start'],
            end=gene_info_dict['end'],
            gene_id=gene_id,
            phenotype_cols=pheno_cols,
            cluster_haplotypes=False,
            database_dir=db_dir
        )

        for check_name in [f'{gene_id}.html', 'integrated_analysis.html']:
            html_path = os.path.join(results_dir, check_name)
            if os.path.exists(html_path):
                size_kb = os.path.getsize(html_path) / 1024
                print(f"[OK] {gene_id}: HTML generated ({size_kb:.0f} KB) -> {check_name}")
                return True

        print(f"[WARN] {gene_id}: HTML not found")
        return False
    except Exception as e:
        import traceback
        print(f"[ERROR] {gene_id}: {e}")
        traceback.print_exc()
        return False


# ========== 主流程 ==========
if __name__ == '__main__':
    print("=" * 60)
    print("水稻单倍型分析 - 谢为博2021论文关键基因")
    print("=" * 60)

    os.makedirs(OUT_DB, exist_ok=True)
    os.makedirs(OUT_RESULTS, exist_ok=True)

    # Load phenotype
    print("\n[1] Loading phenotype data...")
    pheno_df = load_phenotype()
    print(f"  {len(pheno_df)} samples, cols: {[c for c in pheno_df.columns if c != 'SampleID']}")

    # Process each gene
    for gene_info in GENES:
        gene_id = gene_info['gene_id']
        print(f"\n{'='*60}")
        print(f"[{gene_id}] {gene_info['msu_id']}")
        print(f"{'='*60}")

        # Get H5 annotation info
        h5_info = get_gene_info_from_h5(gene_info)
        if h5_info:
            print(f"[INFO] H5: {len(h5_info['positions'])} annotated positions, "
                  f"CDS={h5_info['cds']}, exons={h5_info['exons']}")
        else:
            print(f"[WARN] No H5 annotation found")

        # Build haplotypes
        data = build_haplotype_for_gene(gene_info, pheno_df)
        if data is None:
            print(f"[SKIP] {gene_id}: no valid data")
            continue

        # Save database
        save_database(data, OUT_DB, h5_info)

        # Merge phenotype
        merged, pheno_cols = merge_phenotype(data, pheno_df, OUT_DB)
        print(f"[INFO] Merged: {len(merged)} samples, haplotypes: {dict(merged['Hap_Name'].value_counts())}")

        # Update gene_info.json with H5 annotation data
        gene_dir = os.path.join(OUT_DB, gene_id)
        with open(os.path.join(gene_dir, 'gene_info.json')) as f:
            gi = json.load(f)
        if h5_info:
            gi['exons'] = h5_info.get('exons', [])
            gi['cds'] = h5_info.get('cds', [])

        # Update gene_info for HTML generation
        gi.setdefault('promoter_length', PROMOTER_LENGTH)
        gi.setdefault('promoter_actual_length', PROMOTER_LENGTH)
        gi.setdefault('promoter_expansion_status', 'none')
        if 'promoter_start' not in gi:
            if gi['strand'] == '+':
                gi['promoter_start'] = max(1, gi['gene_start'] - PROMOTER_LENGTH)
                gi['promoter_end'] = gi['gene_start'] - 1
            else:
                gi['promoter_start'] = gi['gene_end'] + 1
                gi['promoter_end'] = gi['gene_end'] + PROMOTER_LENGTH

        with open(os.path.join(gene_dir, 'gene_info.json'), 'w', encoding='utf-8') as f:
            json.dump(gi, f, indent=2, ensure_ascii=False)

    # Build GFF3 for paper genes
    print("\n[2] Building GFF3...")
    gff_lines = ['##gff-version 3']
    for gene_info in GENES:
        gene_id = gene_info['gene_id']
        msu = gene_info['msu_id']
        chrom = gene_info['chrom']
        start = gene_info['start']
        end = gene_info['end']
        strand = gene_info['strand']
        gene_dir = os.path.join(OUT_DB, gene_id)

        with open(os.path.join(gene_dir, 'gene_info.json')) as f:
            gi = json.load(f)

        gff_lines.append(f'{chrom}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};Name={msu}')

        for i, (es, ee) in enumerate(gi.get('exons', [])):
            gff_lines.append(f'{chrom}\t.\tmRNA\t{es}\t{ee}\t.\t{strand}\t.\tID={gene_id}.mRNA;Parent={gene_id}')
            gff_lines.append(f'{chrom}\t.\texon\t{es}\t{ee}\t.\t{strand}\t.\tID={gene_id}.exon{i+1};Parent={gene_id}.mRNA')
            break  # One mRNA/exon

        for i, (cs, ce) in enumerate(gi.get('cds', [])):
            gff_lines.append(f'{chrom}\t.\tCDS\t{cs}\t{ce}\t.\t{strand}\t.\tID={gene_id}.cds{i+1};Parent={gene_id}.mRNA')

    gff_path = r'D:\Desktop\data\水稻\rice_test_genes_paper.gff3'
    with open(gff_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(gff_lines) + '\n')
    print(f"GFF3 saved: {gff_path}")

    # Generate HTML for each gene
    print("\n[3] Generating HTML reports...")
    success = 0
    for gene_info in GENES:
        gene_id = gene_info['gene_id']
        gene_dir = os.path.join(OUT_DB, gene_id)

        with open(os.path.join(gene_dir, 'gene_info.json')) as f:
            gi = json.load(f)

        print(f"\n--- {gene_id} ---")
        if generate_html(gene_id, gi):
            success += 1

    print(f"\n{'='*60}")
    print(f"Done: {success}/{len(GENES)} HTML reports generated")
    print(f"DB: {OUT_DB}")
    print(f"Results: {OUT_RESULTS}")
