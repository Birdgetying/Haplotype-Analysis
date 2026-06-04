#!/usr/bin/env python3
"""
生成水稻单倍型分析HTML报告
从已有的数据库(CSV)生成，无需VCF

Usage: 先运行 run_rice_from_genotype.py 建立数据库，再运行本脚本生成HTML
"""

import sys, os, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, r'd:\Desktop\project1')

import pandas as pd
import json

from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

# ========== 配置 ==========
DATABASE_DIR = r'D:\Desktop\project1\rice_database'
RESULTS_DIR = r'D:\Desktop\project1\rice_results'

GENES = ['LOC_Os01g02660', 'LOC_Os01g25140', 'LOC_Os01g57340']

VCF_FILE = r'D:\Desktop\data\水稻\VCF format\rice4k_geno_add_del.vcf.gz'
GFF_FILE = r'D:\Desktop\data\水稻\rice_test_genes.gff3'


def update_gene_info(gene_dir):
    """更新gene_info.json，补充HTML生成需要的字段"""
    info_path = os.path.join(gene_dir, 'gene_info.json')
    with open(info_path, 'r', encoding='utf-8') as f:
        info = json.load(f)

    # 确保有必要的字段
    info.setdefault('promoter_length', 2000)
    info.setdefault('promoter_actual_length', 2000)
    info.setdefault('promoter_expansion_status', 'none')
    info.setdefault('vcf_mtime', 0)
    info.setdefault('vcf_file', VCF_FILE)
    info.setdefault('extraction_time', '2026-05-09 00:00:00')
    if 'promoter_start' not in info:
        if info['strand'] == '+':
            info['promoter_start'] = max(1, info['gene_start'] - 2000)
            info['promoter_end'] = info['gene_start'] - 1
        else:
            info['promoter_start'] = info['gene_end'] + 1
            info['promoter_end'] = info['gene_end'] + 2000

    with open(info_path, 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2)
    return info


def generate_html(gene_id):
    """为单个基因生成HTML报告"""
    gene_dir = os.path.join(DATABASE_DIR, gene_id)
    results_dir = os.path.join(RESULTS_DIR, gene_id)

    required_files = ['gene_info.json', 'phenotype_data.csv', 'haplotype_data.csv', 'haplotype_samples.csv']
    missing = [f for f in required_files if not os.path.exists(os.path.join(gene_dir, f))]
    if missing:
        print(f"\n[WARN] {gene_id}: 数据库不完整，跳过HTML生成，缺少 {missing}")
        return False

    os.makedirs(results_dir, exist_ok=True)

    # 更新 gene_info
    gene_info = update_gene_info(gene_dir)
    print(f"\n[INFO] {gene_id}: gene_info updated")

    # 读取表型数据
    pheno_path = os.path.join(gene_dir, 'phenotype_data.csv')
    pheno_data = pd.read_csv(pheno_path, encoding='utf-8')
    print(f"[INFO] {gene_id}: {len(pheno_data)} 个样本有表型数据")

    # 获取表型列
    pheno_cols = [c for c in pheno_data.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]

    # 创建分析器
    analyzer = HaplotypePhenotypeAnalyzer(
        vcf_file=VCF_FILE if os.path.exists(VCF_FILE) else None,
        phenotype_file=None,
        output_dir=results_dir,
        gtf_file=GFF_FILE
    )

    # 设置表型数据
    analyzer.phenotype_df = pheno_data

    # 运行分析 (从数据库加载)
    try:
        result = analyzer.analyze_gene(
            chrom=gene_info['chrom'],
            start=gene_info['start'],
            end=gene_info['end'],
            gene_id=gene_id,
            phenotype_cols=pheno_cols,
            cluster_haplotypes=False,
            database_dir=DATABASE_DIR
        )

        # 检查HTML是否生成
        html_path1 = os.path.join(results_dir, f'{gene_id}.html')
        html_path2 = os.path.join(results_dir, 'integrated_analysis.html')
        if os.path.exists(html_path1):
            size_kb = os.path.getsize(html_path1) / 1024
            print(f"[OK] {gene_id}: HTML generated ({size_kb:.1f} KB)")
            return True
        elif os.path.exists(html_path2):
            size_kb = os.path.getsize(html_path2) / 1024
            print(f"[OK] {gene_id}: HTML generated ({size_kb:.1f} KB)")
            return True
        else:
            print(f"[WARN] {gene_id}: HTML not generated")
            return False

    except Exception as e:
        import traceback
        print(f"[ERROR] {gene_id}: {e}")
        traceback.print_exc()
        return False


if __name__ == '__main__':
    print("=" * 60)
    print("水稻单倍型分析 - HTML报告生成")
    print("=" * 60)

    success = 0
    for gene_id in GENES:
        print(f"\n{'='*60}")
        print(f"[{gene_id}]")
        print(f"{'='*60}")
        if generate_html(gene_id):
            success += 1

    print(f"\n{'='*60}")
    print(f"完成: {success}/{len(GENES)} 个基因生成HTML成功")
    print(f"HTML输出目录: {RESULTS_DIR}")
    print("Done!")
