#!/usr/bin/env python3
"""
全基因组单倍型扫描脚本
========================
扫描所有基因，提取单倍型统计信息，生成数据库总表

输出表格包含：
- gene_id: 基因ID
- chrom: 染色体
- start, end: 基因区间
- strand: 链方向
- n_variants: 变异位点数
- n_haplotypes: 单倍型数量
- haplotype_name: 单倍型名称
- haplotype_freq: 单倍型频率
- haplotype_count: 单倍型样本数
- phenotype_mean: 表型均值
- phenotype_sd: 表型标准差
"""

import os
import sys
import argparse
import gzip
import time
import json
import logging
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

# 尝试导入 psutil 用于资源监控
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# 导入主分析模块
try:
    from haplotype_phenotype_analysis import (
        HaplotypeExtractor, DataConfig, setup_logging, PerformanceMonitor
    )
except ImportError:
    print("[ERROR] 无法导入 haplotype_phenotype_analysis 模块")
    print("请确保 haplotype_phenotype_analysis.py 在当前目录")
    sys.exit(1)

# 尝试导入 pysam
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    print("[WARNING] pysam 不可用，单倍型提取功能受限")


# ============================================================================
# 配置
# ============================================================================

class ScanConfig:
    """扫描配置"""
    # 数据路径 (CS-IAAS T2T v1.1)
    VCF_FILE = "/storage/public/home/2024110093/data/Variation/CSIAAS/Core819Samples_ALL.vcf.gz"
    GFF_FILE = "/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1_HC.gff3"
    PHENO_FILE = "/storage/public/home/2024110093/data/Variation/CSIAAS/Phe.txt"
    FASTA_FILE = "/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1.fasta"
    
    # 输出路径
    OUTPUT_DIR = "./results/genome_scan"
    
    # 并行参数
    N_WORKERS = 8  # 并行进程数
    BATCH_SIZE = 100  # 每批处理的基因数
    
    # 过滤参数
    MIN_VARIANTS = 1  # 最小变异数
    MIN_SAMPLES = 5   # 最小样本数


# ============================================================================
# GFF3 解析
# ============================================================================

def parse_gff3_genes(gff_file: str) -> pd.DataFrame:
    """
    从 GFF3 文件解析所有基因信息
    
    Returns:
        DataFrame: gene_id, chrom, start, end, strand
    """
    genes = []
    
    if not os.path.exists(gff_file):
        print(f"[ERROR] GFF3 文件不存在: {gff_file}")
        return pd.DataFrame()
    
    open_func = gzip.open if gff_file.endswith('.gz') else open
    
    print(f"[INFO] 解析 GFF3 文件: {gff_file}")
    start_time = time.time()
    
    with open_func(gff_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            # 只取 gene 行
            if parts[2] != 'gene':
                continue
            
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            # 解析属性
            attr = parts[8]
            gene_id = None
            
            # GFF3 格式: ID=GENEID;...
            for item in attr.split(';'):
                if item.startswith('ID='):
                    gene_id = item[3:]
                    break
            
            if gene_id:
                genes.append({
                    'gene_id': gene_id,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                })
    
    elapsed = time.time() - start_time
    df = pd.DataFrame(genes)
    print(f"[INFO] 解析完成: {len(df)} 个基因 ({elapsed:.1f}s)")
    
    return df


# ============================================================================
# 单基因处理函数（用于并行）
# ============================================================================

def process_single_gene(gene_info: dict, vcf_file: str, pheno_df: pd.DataFrame, 
                        min_samples: int = 5) -> list:
    """
    处理单个基因，返回单倍型统计信息
    
    Returns:
        list of dict: 每个单倍型一行
    """
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    start = gene_info['start']
    end = gene_info['end']
    strand = gene_info['strand']
    
    results = []
    
    try:
        # 创建单倍型提取器
        extractor = HaplotypeExtractor(vcf_file)
        
        # 提取单倍型
        positions, hap_df, hap_sample_df = extractor.extract_region(
            chrom, start, end, min_samples=min_samples
        )
        
        if hap_df is None or len(hap_df) == 0:
            # 没有变异或单倍型
            return [{
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'n_variants': 0,
                'n_haplotypes': 0,
                'haplotype_name': None,
                'haplotype_freq': None,
                'haplotype_count': None,
                'phenotype_mean': None,
                'phenotype_sd': None,
                'status': 'no_variants'
            }]
        
        n_variants = len(positions) if positions else 0
        
        # 合并表型数据
        if 'SampleID' not in hap_sample_df.columns:
            return [{
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'n_variants': n_variants,
                'n_haplotypes': 0,
                'haplotype_name': None,
                'haplotype_freq': None,
                'haplotype_count': None,
                'phenotype_mean': None,
                'phenotype_sd': None,
                'status': 'no_sample_id'
            }]
        
        merged = pd.merge(hap_sample_df, pheno_df, on='SampleID', how='inner')
        
        if len(merged) == 0:
            return [{
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'n_variants': n_variants,
                'n_haplotypes': len(hap_df),
                'haplotype_name': None,
                'haplotype_freq': None,
                'haplotype_count': None,
                'phenotype_mean': None,
                'phenotype_sd': None,
                'status': 'no_phenotype_match'
            }]
        
        # 获取表型列（假设第一个数值列是表型）
        pheno_cols = [c for c in pheno_df.columns if c != 'SampleID' and pheno_df[c].dtype in ['float64', 'int64']]
        if not pheno_cols:
            pheno_cols = [c for c in pheno_df.columns if c != 'SampleID']
        
        pheno_col = pheno_cols[0] if pheno_cols else None
        
        # 统计每个单倍型
        hap_col = 'Hap_Name' if 'Hap_Name' in merged.columns else 'Haplotype'
        hap_stats = merged.groupby(hap_col).agg({
            'SampleID': 'count'
        }).reset_index()
        hap_stats.columns = [hap_col, 'count']
        
        total_samples = len(merged)
        n_haplotypes = len(hap_stats)
        
        for _, row in hap_stats.iterrows():
            hap_name = row[hap_col]
            hap_count = row['count']
            hap_freq = hap_count / total_samples if total_samples > 0 else 0
            
            # 计算表型统计
            if pheno_col:
                hap_data = merged[merged[hap_col] == hap_name][pheno_col].dropna()
                pheno_mean = hap_data.mean() if len(hap_data) > 0 else None
                pheno_sd = hap_data.std() if len(hap_data) > 1 else None
            else:
                pheno_mean = None
                pheno_sd = None
            
            results.append({
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'n_variants': n_variants,
                'n_haplotypes': n_haplotypes,
                'haplotype_name': hap_name,
                'haplotype_freq': round(hap_freq, 4) if hap_freq else None,
                'haplotype_count': hap_count,
                'phenotype_mean': round(pheno_mean, 4) if pheno_mean is not None else None,
                'phenotype_sd': round(pheno_sd, 4) if pheno_sd is not None else None,
                'status': 'success'
            })
        
        return results
        
    except Exception as e:
        import traceback
        error_msg = str(e) + " | " + traceback.format_exc().replace('\n', ' ')
        return [{
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'n_variants': None,
            'n_haplotypes': None,
            'haplotype_name': None,
            'haplotype_freq': None,
            'haplotype_count': None,
            'phenotype_mean': None,
            'phenotype_sd': None,
            'status': f'error: {error_msg[:200]}'
        }]


# ============================================================================
# 主扫描函数
# ============================================================================

def run_genome_scan(vcf_file: str, gff_file: str, pheno_file: str,
                    output_dir: str, n_workers: int = 4,
                    chrom_filter: str = None, gene_filter: list = None,
                    batch_size: int = 100, min_samples: int = 5) -> pd.DataFrame:
    """
    运行全基因组扫描
    
    Args:
        vcf_file: VCF 文件路径
        gff_file: GFF3 文件路径
        pheno_file: 表型文件路径
        output_dir: 输出目录
        n_workers: 并行进程数
        chrom_filter: 染色体过滤（如 'chr1A'）
        gene_filter: 基因ID过滤列表
        batch_size: 批处理大小
        min_samples: 最小样本数
    
    Returns:
        DataFrame: 扫描结果总表
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 性能监控
    perf = PerformanceMonitor()
    perf.start()
    
    print("=" * 60)
    print("全基因组单倍型扫描")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"VCF: {vcf_file}")
    print(f"GFF: {gff_file}")
    print(f"表型: {pheno_file}")
    print(f"输出: {output_dir}")
    print(f"并行进程: {n_workers}")
    print("=" * 60)
    
    # 1. 解析基因列表
    perf.step_start("parse_gff")
    genes_df = parse_gff3_genes(gff_file)
    perf.step_end("parse_gff")
    
    if len(genes_df) == 0:
        print("[ERROR] 未找到任何基因")
        return pd.DataFrame()
    
    # 应用过滤
    if chrom_filter:
        genes_df = genes_df[genes_df['chrom'] == chrom_filter]
        print(f"[INFO] 过滤后基因数 (chrom={chrom_filter}): {len(genes_df)}")
    
    if gene_filter:
        genes_df = genes_df[genes_df['gene_id'].isin(gene_filter)]
        print(f"[INFO] 过滤后基因数 (gene_filter): {len(genes_df)}")
    
    total_genes = len(genes_df)
    print(f"\n[INFO] 待扫描基因数: {total_genes}")
    
    # 2. 加载表型数据
    perf.step_start("load_phenotype")
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    # 标准化列名
    if 'SampleID' not in pheno_df.columns:
        if pheno_df.columns[0] not in ['SampleID', 'sample_id', 'ID']:
            pheno_df.columns = ['SampleID'] + list(pheno_df.columns[1:])
        else:
            pheno_df = pheno_df.rename(columns={pheno_df.columns[0]: 'SampleID'})
    print(f"[INFO] 表型样本数: {len(pheno_df)}, 列: {list(pheno_df.columns)}")
    perf.step_end("load_phenotype")
    
    # 3. 批量处理
    perf.step_start("scan_genes")
    all_results = []
    processed = 0
    errors = 0
    
    gene_list = genes_df.to_dict('records')
    
    # 使用进度条
    print(f"\n[INFO] 开始扫描...")
    
    # 串行处理（更稳定）或并行处理
    if n_workers == 1:
        # 串行处理
        for i, gene_info in enumerate(gene_list):
            results = process_single_gene(gene_info, vcf_file, pheno_df, min_samples)
            all_results.extend(results)
            processed += 1
            
            if processed % 100 == 0:
                print(f"  进度: {processed}/{total_genes} ({processed*100/total_genes:.1f}%)")
            
            if results[0].get('status', '').startswith('error'):
                errors += 1
    else:
        # 并行处理（批量）
        for batch_start in range(0, total_genes, batch_size):
            batch_end = min(batch_start + batch_size, total_genes)
            batch = gene_list[batch_start:batch_end]
            
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(process_single_gene, gene, vcf_file, pheno_df, min_samples): gene
                    for gene in batch
                }
                
                for future in as_completed(futures):
                    try:
                        results = future.result()
                        all_results.extend(results)
                        processed += 1
                        
                        if results[0].get('status', '').startswith('error'):
                            errors += 1
                    except Exception as e:
                        errors += 1
                        print(f"  [ERROR] 处理失败: {e}")
            
            print(f"  进度: {processed}/{total_genes} ({processed*100/total_genes:.1f}%)")
    
    perf.step_end("scan_genes")
    
    # 4. 生成结果表
    perf.step_start("save_results")
    results_df = pd.DataFrame(all_results)
    
    # 保存完整结果
    output_file = os.path.join(output_dir, "haplotype_database.csv")
    results_df.to_csv(output_file, index=False)
    print(f"\n[INFO] 完整结果已保存: {output_file}")
    
    # 生成基因级汇总
    gene_summary = results_df.groupby('gene_id').agg({
        'chrom': 'first',
        'start': 'first',
        'end': 'first',
        'strand': 'first',
        'n_variants': 'first',
        'n_haplotypes': 'first',
        'status': 'first'
    }).reset_index()
    
    summary_file = os.path.join(output_dir, "gene_summary.csv")
    gene_summary.to_csv(summary_file, index=False)
    print(f"[INFO] 基因汇总已保存: {summary_file}")
    
    perf.step_end("save_results")
    
    # 5. 打印统计
    print("\n" + "=" * 60)
    print("扫描完成统计")
    print("=" * 60)
    print(f"总基因数: {total_genes}")
    print(f"成功处理: {processed - errors}")
    print(f"处理失败: {errors}")
    print(f"结果行数: {len(results_df)}")
    
    # 状态统计
    if 'status' in results_df.columns:
        status_counts = results_df.groupby('gene_id')['status'].first().value_counts()
        print(f"\n状态分布:")
        for status, count in status_counts.items():
            print(f"  {status}: {count}")
    
    # 性能报告
    perf.report_performance()
    
    return results_df


# ============================================================================
# 命令行入口
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="全基因组单倍型扫描 - 生成单倍型数据库",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--vcf", default=ScanConfig.VCF_FILE,
                        help="VCF 文件路径")
    parser.add_argument("--gff", default=ScanConfig.GFF_FILE,
                        help="GFF3 注释文件路径")
    parser.add_argument("--phenotype", default=ScanConfig.PHENO_FILE,
                        help="表型文件路径")
    parser.add_argument("--output-dir", default=ScanConfig.OUTPUT_DIR,
                        help="输出目录")
    parser.add_argument("--workers", type=int, default=1,
                        help="并行进程数 (1=串行)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="批处理大小")
    parser.add_argument("--chrom", default=None,
                        help="只扫描指定染色体")
    parser.add_argument("--genes", nargs="+", default=None,
                        help="只扫描指定基因列表")
    parser.add_argument("--min-samples", type=int, default=5,
                        help="最小样本数阈值")
    
    args = parser.parse_args()
    
    # 运行扫描
    run_genome_scan(
        vcf_file=args.vcf,
        gff_file=args.gff,
        pheno_file=args.phenotype,
        output_dir=args.output_dir,
        n_workers=args.workers,
        chrom_filter=args.chrom,
        gene_filter=args.genes,
        batch_size=args.batch_size,
        min_samples=args.min_samples
    )


if __name__ == "__main__":
    main()
