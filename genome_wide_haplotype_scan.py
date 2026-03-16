#!/usr/bin/env python3
"""
全基因组单倍型扫描脚本
========================
扫描所有基因，提取单倍型统计信息，生成数据库总表，
并进行单倍型-表型关联分析，生成HTML可视化报告

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
- p_value: 关联分析P值（可选）
- significance: 是否显著（可选）
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
        HaplotypeExtractor, DataConfig, setup_logging, PerformanceMonitor,
        HaplotypePhenotypeAnalyzer, HaplotypeReporter
    )
    HAPLOTYPE_MODULE_AVAILABLE = True
except ImportError as e:
    print(f"[WARNING] 无法完全导入 haplotype_phenotype_analysis 模块: {e}")
    print("[WARNING] 将使用基础功能运行")
    HAPLOTYPE_MODULE_AVAILABLE = False
    # 定义占位符类
    class HaplotypeExtractor:
        def __init__(self, vcf_file):
            raise NotImplementedError("haplotype_phenotype_analysis 模块不可用")
    class HaplotypePhenotypeAnalyzer:
        pass
    class HaplotypeReporter:
        pass
    class DataConfig:
        pass
    def setup_logging(*args, **kwargs):
        pass
    class PerformanceMonitor:
        def start(self): pass
        def step_start(self, name): pass
        def step_end(self, name): return 0
        def report_performance(self): pass

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
    
    # 分析参数
    PVALUE_THRESHOLD = 0.05  # P值显著性阈值
    CORRECTION_METHOD = 'bonferroni'  # 多重检验校正方法


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


def analyze_gene_association(gene_info: dict, vcf_file: str, pheno_df: pd.DataFrame,
                            min_samples: int = 5, pvalue_threshold: float = 0.05,
                            gff_file: str = None) -> dict:
    """
    对单个基因进行单倍型-表型关联分析，并生成综合HTML图
    
    Returns:
        dict: 包含关联分析结果和HTML路径的字典
    """
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    start = gene_info['start']
    end = gene_info['end']
    strand = gene_info['strand']
    
    result = {
        'gene_id': gene_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'n_variants': 0,
        'n_haplotypes': 0,
        'p_value': None,
        'pve_percent': None,
        'significant': False,
        'status': 'pending',
        'html_file': None
    }
    
    try:
        # 创建单倍型提取器
        extractor = HaplotypeExtractor(vcf_file)
        
        # 提取单倍型
        positions, hap_df, hap_sample_df = extractor.extract_region(
            chrom, start, end, min_samples=min_samples
        )
        
        if hap_df is None or len(hap_df) == 0 or len(hap_sample_df) < min_samples:
            result['status'] = 'insufficient_data'
            return result
        
        result['n_variants'] = len(positions) if positions else 0
        result['n_haplotypes'] = len(hap_df)
        
        # 合并表型数据
        if 'SampleID' not in hap_sample_df.columns:
            result['status'] = 'no_sample_id'
            return result
        
        merged = pd.merge(hap_sample_df, pheno_df, on='SampleID', how='inner')
        
        if len(merged) < min_samples:
            result['status'] = 'insufficient_phenotype'
            return result
        
        # 获取表型列
        pheno_cols = [c for c in pheno_df.columns if c != 'SampleID' and 
                      pheno_df[c].dtype in ['float64', 'int64']]
        if not pheno_cols:
            pheno_cols = [c for c in pheno_df.columns if c != 'SampleID']
        
        if not pheno_cols:
            result['status'] = 'no_phenotype'
            return result
        
        pheno_col = pheno_cols[0]
        
        # 进行关联分析
        hap_col = 'Hap_Name' if 'Hap_Name' in merged.columns else 'Haplotype'
        
        # 过滤掉Other单倍型
        analysis_df = merged[merged[hap_col] != 'Other'].copy()
        
        if len(analysis_df) < min_samples:
            result['status'] = 'insufficient_haplotypes'
            return result
        
        # ANOVA检验
        from scipy import stats
        
        hap_means = analysis_df.groupby(hap_col)[pheno_col].mean()
        groups = [analysis_df[analysis_df[hap_col] == h][pheno_col].dropna().values 
                  for h in hap_means.index if len(analysis_df[analysis_df[hap_col] == h]) >= 2]
        
        if len(groups) < 2:
            result['status'] = 'single_haplotype'
            return result
        
        f_stat, p_value = stats.f_oneway(*groups)
        
        result['p_value'] = p_value
        result['significant'] = p_value < pvalue_threshold
        result['status'] = 'success'
        
        # 计算PVE
        grand_mean = analysis_df[pheno_col].mean()
        ss_between = sum(len(analysis_df[analysis_df[hap_col] == h]) * 
                        (analysis_df[analysis_df[hap_col] == h][pheno_col].mean() - grand_mean)**2 
                        for h in hap_means.index)
        ss_total = ((analysis_df[pheno_col] - grand_mean)**2).sum()
        
        if ss_total > 0:
            result['pve_percent'] = (ss_between / ss_total) * 100
        
        # 如果显著，生成综合HTML图
        if result['significant'] and gff_file:
            try:
                # 创建输出目录
                gene_output_dir = os.path.join("./results/genome_scan/significant_genes", gene_id)
                os.makedirs(gene_output_dir, exist_ok=True)
                
                # 使用 HaplotypePhenotypeAnalyzer 生成完整分析
                analyzer = HaplotypePhenotypeAnalyzer(
                    vcf_file=vcf_file,
                    phenotype_file=None,
                    output_dir=gene_output_dir,
                    gtf_file=gff_file
                )
                
                # 设置已加载的数据
                analyzer.positions = positions
                analyzer.hap_df = hap_df
                analyzer.hap_sample_df = hap_sample_df
                analyzer.pheno_df = pheno_df
                
                # 运行完整分析流程（会生成综合HTML）
                analysis_result = analyzer.analyze_gene(
                    gene_id=gene_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    phenotypes=[pheno_col]
                )
                
                # 查找生成的HTML文件
                integrated_html = os.path.join(gene_output_dir, "integrated_analysis.html")
                if os.path.exists(integrated_html):
                    result['html_file'] = integrated_html
                    print(f"[INFO] 生成综合HTML图: {integrated_html}")
                
            except Exception as e:
                print(f"[WARNING] 生成综合HTML图失败 ({gene_id}): {e}")
        
        return result
        
    except Exception as e:
        import traceback
        error_msg = str(e) + " | " + traceback.format_exc().replace('\n', ' ')
        result['status'] = f'error: {error_msg[:200]}'
        return result


# ============================================================================
# 主扫描函数
# ============================================================================

def run_genome_scan(vcf_file: str, gff_file: str, pheno_file: str,
                    output_dir: str, n_workers: int = 4,
                    chrom_filter: str = None, gene_filter: list = None,
                    batch_size: int = 100, min_samples: int = 5,
                    run_analysis: bool = True, pvalue_threshold: float = 0.05,
                    generate_html: bool = True) -> pd.DataFrame:
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
        run_analysis: 是否进行关联分析
        pvalue_threshold: P值显著性阈值
        generate_html: 是否生成HTML报告
    
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
    
    # 4. 生成结果表（在关联分析前创建）
    results_df = pd.DataFrame(all_results)
    
    # 5. 关联分析（可选）
    analysis_results = []
    if run_analysis:
        print("\n[INFO] 开始关联分析...")
        perf.step_start("association_analysis")
        
        # 对每个基因进行关联分析
        analysis_genes = results_df[results_df['status'] == 'success']['gene_id'].unique()
        
        if n_workers == 1:
            # 串行分析
            for i, gene_id in enumerate(analysis_genes):
                gene_info = genes_df[genes_df['gene_id'] == gene_id].iloc[0].to_dict()
                assoc_result = analyze_gene_association(
                    gene_info, vcf_file, pheno_df, min_samples, pvalue_threshold, gff_file
                )
                analysis_results.append(assoc_result)
                
                if (i + 1) % 100 == 0:
                    print(f"  分析进度: {i+1}/{len(analysis_genes)}")
        else:
            # 并行分析
            analysis_genes_list = genes_df[genes_df['gene_id'].isin(analysis_genes)].to_dict('records')
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(analyze_gene_association, g, vcf_file, pheno_df, min_samples, pvalue_threshold, gff_file): g
                    for g in analysis_genes_list
                }
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        analysis_results.append(result)
                    except Exception as e:
                        print(f"  [ERROR] 分析失败: {e}")
        
        # 保存关联分析结果
        assoc_df = pd.DataFrame(analysis_results)
        assoc_file = os.path.join(output_dir, "association_analysis.csv")
        assoc_df.to_csv(assoc_file, index=False)
        print(f"[INFO] 关联分析结果已保存: {assoc_file}")
        
        # 显著基因
        sig_genes = assoc_df[assoc_df['significant'] == True]
        print(f"[INFO] 显著关联基因数 (P < {pvalue_threshold}): {len(sig_genes)}")
        
        # 保存显著基因列表
        if len(sig_genes) > 0:
            sig_file = os.path.join(output_dir, "significant_genes.csv")
            sig_genes.to_csv(sig_file, index=False)
            print(f"[INFO] 显著基因列表已保存: {sig_file}")
        
        perf.step_end("association_analysis")
    
    # 6. 生成HTML报告
    if generate_html and run_analysis and len(analysis_results) > 0:
        print("\n[INFO] 生成HTML报告...")
        perf.step_start("generate_html")
        generate_genome_scan_html_report(results_df, assoc_df, output_dir, pheno_file)
        perf.step_end("generate_html")
    
    # 7. 保存结果表
    perf.step_start("save_results")
    
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
    parser.add_argument("--no-analysis", action="store_true",
                        help="跳过关联分析（只生成单倍型数据库）")
    parser.add_argument("--pvalue-threshold", type=float, default=0.05,
                        help="P值显著性阈值")
    parser.add_argument("--no-html", action="store_true",
                        help="不生成HTML报告")
    
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
        min_samples=args.min_samples,
        run_analysis=not args.no_analysis,
        pvalue_threshold=args.pvalue_threshold,
        generate_html=not args.no_html
    )


# ============================================================================
# HTML报告生成
# ============================================================================

def generate_genome_scan_html_report(haplo_db_df: pd.DataFrame, assoc_df: pd.DataFrame,
                                     output_dir: str, pheno_file: str) -> str:
    """
    生成全基因组扫描HTML报告
    
    Returns:
        str: 生成的HTML文件路径
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import io
    import base64
    
    # 安全格式化函数
    def fmt(val, default='NA'):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return default
        if isinstance(val, float):
            return f"{val:.4f}"
        return str(val)
    
    # 读取表型信息
    pheno_name = "Phenotype"
    try:
        pheno_df = pd.read_csv(pheno_file, sep='\t')
        pheno_cols = [c for c in pheno_df.columns if c != 'SampleID']
        if pheno_cols:
            pheno_name = pheno_cols[0]
    except:
        pass
    
    # 统计信息
    total_genes = len(haplo_db_df['gene_id'].unique())
    success_genes = len(haplo_db_df[haplo_db_df['status'] == 'success']['gene_id'].unique())
    
    # 关联分析统计
    if assoc_df is not None and len(assoc_df) > 0:
        analyzed = len(assoc_df[assoc_df['status'] == 'success'])
        sig_count = len(assoc_df[assoc_df['significant'] == True])
        total_tests = len(assoc_df[assoc_df['p_value'].notna()])
    else:
        analyzed = 0
        sig_count = 0
        total_tests = 0
    
    # 生成P值分布图
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # P值直方图
    if assoc_df is not None and 'p_value' in assoc_df.columns:
        pvals = assoc_df['p_value'].dropna()
        if len(pvals) > 0:
            axes[0].hist(pvals, bins=50, color='#3498db', edgecolor='white', alpha=0.7)
            axes[0].axvline(x=0.05, color='red', linestyle='--', label='P=0.05')
            axes[0].set_xlabel('P-value')
            axes[0].set_ylabel('Count')
            axes[0].set_title('P-value Distribution')
            axes[0].legend()
    
    # 曼哈顿图风格的P值图
    if assoc_df is not None and len(assoc_df) > 0:
        plot_df = assoc_df[assoc_df['p_value'].notna()].copy()
        if len(plot_df) > 0:
            # 按染色体排序
            chrom_order = sorted(plot_df['chrom'].unique(), 
                                key=lambda x: (len(x), x))
            plot_df['chrom_order'] = plot_df['chrom'].map({c: i for i, c in enumerate(chrom_order)})
            plot_df = plot_df.sort_values(['chrom_order', 'start'])
            
            # 转换P值为 -log10
            plot_df['-log10p'] = -np.log10(plot_df['p_value'])
            
            # 散点图
            colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']
            for i, chrom in enumerate(plot_df['chrom'].unique()):
                chrom_data = plot_df[plot_df['chrom'] == chrom]
                axes[1].scatter(chrom_data['start'], chrom_data['-log10p'], 
                               c=colors[i % len(colors)], alpha=0.6, s=10, label=chrom)
            
            axes[1].axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
            axes[1].set_xlabel('Genomic Position')
            axes[1].set_ylabel('-log10(P-value)')
            axes[1].set_title('Genome-wide Association Results')
            axes[1].legend(fontsize=8, loc='upper right')
    
    plt.tight_layout()
    
    # 保存图表为base64
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', dpi=150, bbox_inches='tight')
    img_buffer.seek(0)
    img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
    plt.close()
    
    # 生成显著基因表格（包含综合HTML图链接）
    sig_table_html = ""
    sig_genes_with_html = []
    if assoc_df is not None and len(assoc_df[assoc_df['significant'] == True]) > 0:
        sig_genes = assoc_df[assoc_df['significant'] == True].sort_values('p_value').head(50)
        for _, row in sig_genes.iterrows():
            # 检查是否有综合HTML图
            html_link = ""
            if 'html_file' in row and row['html_file'] and os.path.exists(row['html_file']):
                rel_path = os.path.relpath(row['html_file'], output_dir)
                html_link = f'<a href="{rel_path}" target="_blank" style="color:#3498db;text-decoration:none;">📊 查看</a>'
                sig_genes_with_html.append({
                    'gene_id': row['gene_id'],
                    'p_value': row['p_value'],
                    'html_file': row['html_file']
                })
            
            sig_table_html += f"""
            <tr>
                <td>{row['gene_id']}</td>
                <td>{row['chrom']}</td>
                <td>{row['start']:,}</td>
                <td>{row['end']:,}</td>
                <td>{fmt(row['n_variants'])}</td>
                <td>{fmt(row['n_haplotypes'])}</td>
                <td style="color:red;font-weight:bold;">{fmt(row['p_value'], 'NA')}</td>
                <td>{fmt(row['pve_percent'], 'NA')}%</td>
                <td>{html_link}</td>
            </tr>
            """
    
    # HTML内容
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome-Wide Haplotype Scan Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f7fa;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .summary-box {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-card.green {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }}
        .stat-card.orange {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .stat-card.yellow {{
            background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        }}
        .stat-value {{
            font-size: 32px;
            font-weight: bold;
        }}
        .stat-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: 600;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .chart-container {{
            margin: 20px 0;
            text-align: center;
        }}
        .chart-container img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            color: #7f8c8d;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>全基因组单倍型扫描分析报告</h1>
        <p><strong>表型:</strong> {pheno_name}</p>
        <p><strong>生成时间:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>统计摘要</h2>
        <div class="summary-box">
            <div class="stat-card">
                <div class="stat-value">{total_genes:,}</div>
                <div class="stat-label">总基因数</div>
            </div>
            <div class="stat-card green">
                <div class="stat-value">{success_genes:,}</div>
                <div class="stat-label">成功处理</div>
            </div>
            <div class="stat-card orange">
                <div class="stat-value">{analyzed:,}</div>
                <div class="stat-label">关联分析</div>
            </div>
            <div class="stat-card yellow">
                <div class="stat-value">{sig_count:,}</div>
                <div class="stat-label">显著基因 (P&lt;0.05)</div>
            </div>
        </div>
        
        <h2>关联分析结果</h2>
        <div class="chart-container">
            <img src="data:image/png;base64,{img_base64}" alt="Association Analysis Plot">
        </div>
        
        <h2>显著关联基因列表</h2>
        {"<table><thead><tr><th>Gene ID</th><th>Chrom</th><th>Start</th><th>End</th><th>Variants</th><th>Haplotypes</th><th>P-value</th><th>PVE %</th><th>综合图</th></tr></thead><tbody>" + sig_table_html + "</tbody></table>" if sig_table_html else "<p>未发现显著关联基因</p>"}
        
        {f'<h2>综合可视化图（基因结构+单倍型序列+效应图+箱线图）</h2><div class="summary-box">' + ''.join([f'<div style="margin:10px 0;padding:15px;background:#f8f9fa;border-radius:8px;"><h3>{g["gene_id"]}</h3><p>P-value: {g["p_value"]:.2e}</p><iframe src="{os.path.relpath(g["html_file"], output_dir)}" width="100%" height="600px" style="border:1px solid #ddd;border-radius:4px;"></iframe></div>' for g in sig_genes_with_html[:5]]) + '</div>' if sig_genes_with_html else ''}
        
        <div class="footer">
            <p>全基因组单倍型扫描分析系统 | Genome-Wide Haplotype Scan Analysis</p>
        </div>
    </div>
</body>
</html>
"""
    
    # 保存HTML文件
    html_file = os.path.join(output_dir, "genome_scan_report.html")
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"[INFO] HTML报告已生成: {html_file}")
    return html_file


if __name__ == "__main__":
    main()
