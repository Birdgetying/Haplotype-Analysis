#!/usr/bin/env python3
"""
小麦 Rht 基因 - 单倍型-株高关联分析 (最终版)
==============================================
数据: Core819Samples SNP+INDEL VCF, Phe.raw2 (PH=株高), CS-IAAS GFF3
"""
import sys, os, io, gzip, json, time, re
from collections import Counter
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, 'd:/Desktop/project1')
from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer, DataConfig

# Override HPC paths with local data
DataConfig.GENOME_BASE = 'd:/Desktop/wheat_data/'
DataConfig.FASTA_PATH = DataConfig.GENOME_BASE + 'CS-IAAS_v1.1.fasta'

DESKTOP = 'd:/Desktop'
SNP_VCF = f'{DESKTOP}/wheat_data/Core819Samples_snp.filter.final.id_gt.813m.vcf.gz'
INDEL_VCF = f'{DESKTOP}/wheat_data/Core819Samples_indel.filter.final.id_gt.813m.vcf.gz'
SV_VCF = f'{DESKTOP}/wheat_data/SV.new.vcf.gz'
GFF_FILE = f'{DESKTOP}/wheat_data/CS-IAAS_v1.1_HC.gene.gff3'
FASTA_FILE = f'{DESKTOP}/wheat_data/CS-IAAS_v1.1.fasta'
PHENO_FILE = f'{DESKTOP}/wheat_data/Phe.raw2'
OUTPUT_DIR = f'{DESKTOP}/project1/wheat_rht_results'
DB_DIR = f'{DESKTOP}/project1/wheat_rht_database'
VCF_OUT = f'{DESKTOP}/wheat_rht_vcf'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(DB_DIR, exist_ok=True)
os.makedirs(VCF_OUT, exist_ok=True)

PROMOTER = 2000
FLANK = 5000

# Rht基因 - 使用GFF3里最接近的基因坐标
# Rht-A1: GCA=4A:492,561,193-492,563,420, GFF3最近基因=CSIAAS4AG0494000HC(492,432,167-492,436,930)和CSIAAS4AG0494200HC(492,629,767-492,632,880)
#   直接用GCA坐标+扩展
# Rht-B1: GCA=4B:23,156,849-23,159,304, 但Chr4B在该区域有356Kb空白, 最近标记在23,237,232(80Kb外)
#   扩展搜索范围到±100Kb
# Rht-D1: GFF3候选 CSIAAS4DG0091400HC(18,775,201-18,781,609) + CSIAAS4DG0091600HC(18,782,933-18,787,353)
#   用这个区间

RHT = {
    'TraesCS4A02G271000': {
        'name': 'Rht-A1', 'chrom': 'Chr4A',
        'start': 492559193, 'end': 492565420,
        'strand': '-', 'note': 'GCA坐标±2Kb promoter',
    },
    'TraesCS4B02G043100': {
        'name': 'Rht-B1', 'chrom': 'Chr4B',
        'start': 22880000, 'end': 23240000,
        'strand': '+', 'note': '标记空白区,扩展±180Kb',
    },
    'TraesCS4D02G040400': {
        'name': 'Rht-D1', 'chrom': 'Chr4D',
        'start': 18775201, 'end': 18787353,
        'strand': '+', 'note': 'GFF3相邻基因区间',
    },
}


def load_phenotype():
    print("=" * 60)
    print("Load phenotype (Phe.raw2)")
    print("=" * 60)
    df = pd.read_csv(PHENO_FILE, sep='\t', encoding='utf-8')
    pheno = df[['VCFID', 'PH']].copy()
    pheno.columns = ['SampleID', 'PlantHeight']
    pheno['SampleID'] = pheno['SampleID'].astype(str).str.strip()
    pheno['PlantHeight'] = pd.to_numeric(pheno['PlantHeight'], errors='coerce')
    pheno = pheno.dropna(subset=['PlantHeight'])
    print(f"Records: {len(pheno)}, PH range: {pheno['PlantHeight'].min():.1f}-{pheno['PlantHeight'].max():.1f}")
    return pheno


def extract_regions():
    """扫描SNP+INDEL VCF,提取Rht基因区域"""
    print("\n" + "=" * 60)
    print("Extract Rht gene regions from VCFs")
    print("=" * 60)

    gene_vcfs = {}
    for gene_id, info in RHT.items():
        chrom = info['chrom']
        r_start = info['start'] - PROMOTER - FLANK
        r_end = info['end'] + FLANK

        out_path = f'{VCF_OUT}/{gene_id}.vcf'
        f_out = open(out_path, 'w', encoding='utf-8')

        gene_vcfs[gene_id] = {'f': f_out, 'chrom': chrom,
                               'start': r_start, 'end': r_end, 'count': 0}

    for label, vcf_path in [('SNP', SNP_VCF), ('INDEL', INDEL_VCF), ('SV', SV_VCF)]:
        print(f"  Scanning {label} VCF...")
        t0 = time.time()
        with gzip.open(vcf_path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('##'):
                    for gw in gene_vcfs.values():
                        gw['f'].write(line)
                    continue
                if line.startswith('#CHROM'):
                    for gw in gene_vcfs.values():
                        gw['f'].write(line)
                    continue
                parts = line.split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                for gene_id, gw in gene_vcfs.items():
                    if chrom == gw['chrom'] and gw['start'] <= pos <= gw['end']:
                        gw['f'].write(line)
                        gw['count'] += 1
        print(f"    {time.time()-t0:.0f}s")

    for gene_id, gw in gene_vcfs.items():
        gw['f'].close()
        out_path = f'{VCF_OUT}/{gene_id}.vcf'
        size_kb = os.path.getsize(out_path) / 1024
        print(f"  {RHT[gene_id]['name']}: {gw['count']} variants, {size_kb:.1f} KB")
        gene_vcfs[gene_id] = out_path

    return gene_vcfs


def run_analysis(gene_vcfs, pheno_df):
    print("\n" + "=" * 60)
    print("Haplotype Analysis + Scoring")
    print("=" * 60)

    for gene_id in RHT:
        vcf_path = gene_vcfs[gene_id]
        info = RHT[gene_id]
        name, chrom, start, end = info['name'], info['chrom'], info['start'], info['end']

        with open(vcf_path, 'r') as f:
            n_var = sum(1 for l in f if not l.startswith('#'))

        if n_var == 0:
            print(f"\n  {name}: [SKIP] No variants ({info['note']})")
            continue

        print(f"\n{'='*50}")
        print(f"  {name} ({gene_id}) — {info['note']}")
        print(f"  {chrom}:{start:,}-{end:,} | {n_var} variants")
        print(f"{'='*50}")

        results_dir = f'{OUTPUT_DIR}/{gene_id}'
        os.makedirs(results_dir, exist_ok=True)

        try:
            analyzer = HaplotypePhenotypeAnalyzer(
                vcf_file=vcf_path, phenotype_file=None,
                output_dir=results_dir, gtf_file=GFF_FILE
            )
            analyzer.phenotype_df = pheno_df

            result = analyzer.analyze_gene(
                chrom=chrom, start=start, end=end,
                gene_id=gene_id, phenotype_cols=['PlantHeight'],
                cluster_haplotypes=False, database_dir=DB_DIR,
            )

            htmls = sorted([f for f in os.listdir(results_dir) if f.endswith('.html')])
            if htmls:
                print(f"  Generated: {len(htmls)} HTML files")
                for h in htmls:
                    sz = os.path.getsize(os.path.join(results_dir, h)) / 1024
                    print(f"    {h} ({sz:.1f} KB)")
            else:
                print(f"  No HTML generated")
                for f in sorted(os.listdir(results_dir))[:10]:
                    print(f"    {f}")
        except Exception as e:
            import traceback
            print(f"  [ERROR] {e}")
            traceback.print_exc()


if __name__ == '__main__':
    t0 = time.time()
    print("=" * 60)
    print("  Wheat Rht Haplotype-PlantHeight Analysis")
    print("=" * 60)

    pheno_df = load_phenotype()
    # 检查是否已有提取好的VCF
    gene_vcfs = {}
    all_exist = True
    for gene_id in RHT:
        vcf_path = f'{VCF_OUT}/{gene_id}.vcf'
        if os.path.exists(vcf_path):
            gene_vcfs[gene_id] = vcf_path
        else:
            all_exist = False

    if not all_exist:
        gene_vcfs = extract_regions()
    else:
        print(f'\nUsing existing VCF extracts (skip re-scan)')

    run_analysis(gene_vcfs, pheno_df)

    print(f"\nTotal: {int((time.time()-t0)/60)} min")
    print(f"Results: {OUTPUT_DIR}")
