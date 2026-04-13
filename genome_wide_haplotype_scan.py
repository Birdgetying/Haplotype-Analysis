#!/usr/bin/env python3
"""
单倍型数据集构建脚本
====================
针对指定基因构建单倍型数据集

目录结构:
    database/                           # 数据库文件夹
    ├── {gene_id}/
    │   ├── gene_info.json              # 基因基本信息
    │   ├── haplotype_data.csv          # 单倍型数据
    │   ├── haplotype_samples.csv       # 样本-单倍型对应
    │   ├── haplotype_stats.csv         # 单倍型统计
    │   ├── phenotype_data.csv          # 表型数据
    │   └── association_result.csv      # 关联分析结果
    └── summary.csv                      # 所有基因汇总

    results/                            # 结果文件夹
    └── {gene_id}/
        └── integrated_analysis.html    # 综合分析HTML图

默认处理基因:
- CSIAAS1BG1157200HC (di19)
- CSIAAS4BG0701800HC (hox)
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
        HaplotypePhenotypeAnalyzer
    )
    HAPLOTYPE_MODULE_AVAILABLE = True
    print("[INFO] haplotype_phenotype_analysis 模块导入成功")
except ImportError as e:
    import traceback
    print(f"[WARNING] 无法导入 haplotype_phenotype_analysis 模块")
    print(f"[WARNING] 错误详情: {e}")
    traceback.print_exc()
    print("[WARNING] 将使用内置基础功能运行")
    HAPLOTYPE_MODULE_AVAILABLE = False
    
    # 定义基础的 DataConfig
    class DataConfig:
        pass
    
    def setup_logging(*args, **kwargs):
        pass
    
    class PerformanceMonitor:
        def start(self): pass
        def step_start(self, name): pass
        def step_end(self, name): return 0
        def report_performance(self): pass
    
    class HaplotypePhenotypeAnalyzer:
        pass

# 尝试导入 pysam
try:
    import pysam
    PYSAM_AVAILABLE = True
    print("[INFO] pysam 导入成功")
except ImportError:
    PYSAM_AVAILABLE = False
    print("[WARNING] pysam 不可用，单倍型提取功能受限")


# ============================================================================
# VCF 工具函数
# ============================================================================

def create_subset_vcf(input_vcf: str, chrom: str, start: int, end: int, 
                      output_vcf: str, sample_ids: list = None):
    """
    从原始VCF中提取指定区域和样本的子集VCF（纯Python实现，不依赖pysam）
    
    Args:
        input_vcf: 输入VCF文件路径
        chrom: 染色体
        start: 起始位置
        end: 结束位置
        output_vcf: 输出VCF文件路径
        sample_ids: 要保留的样本ID列表（None表示保留所有样本）
    """
    import gzip
    
    print(f"[DEBUG] 开始提取VCF子集: {chrom}:{start}-{end}")
    
    # 打开输入VCF
    if input_vcf.endswith('.gz'):
        f_in = gzip.open(input_vcf, 'rt', encoding='utf-8', errors='ignore')
    else:
        f_in = open(input_vcf, 'r', encoding='utf-8', errors='ignore')
    
    # 读取头部
    header_lines = []
    sample_line = None
    all_samples = []
    
    for line in f_in:
        if line.startswith('##'):
            header_lines.append(line.rstrip('\n'))
        elif line.startswith('#CHROM'):
            sample_line = line.rstrip('\n')
            parts = line.rstrip('\n').split('\t')
            all_samples = parts[9:] if len(parts) > 9 else []
            break
    
    # 确定要保留的样本
    if sample_ids:
        available_samples = set(all_samples)
        filtered_samples = [s for s in sample_ids if s in available_samples]
        if len(filtered_samples) != len(sample_ids):
            missing = set(sample_ids) - available_samples
            print(f"[WARNING] 以下样本在VCF中不存在: {missing}")
    else:
        filtered_samples = all_samples
    
    # 如果需要过滤样本，重新构建样本行
    if sample_ids and filtered_samples != all_samples:
        # 找到样本在原始样本列表中的索引
        sample_indices = [all_samples.index(s) for s in filtered_samples if s in all_samples]
        # 重新构建样本行
        parts = sample_line.split('\t')[:9]  # 前9个字段
        sample_fields = filtered_samples
        new_sample_line = '\t'.join(parts + sample_fields)
    else:
        new_sample_line = sample_line
        sample_indices = list(range(len(all_samples)))
    
    # 打开输出VCF
    if output_vcf.endswith('.gz'):
        f_out = gzip.open(output_vcf, 'wt', encoding='utf-8', compresslevel=6)
    else:
        f_out = open(output_vcf, 'w', encoding='utf-8')
    
    # 写入头部
    for line in header_lines:
        f_out.write(line + '\n')
    f_out.write(new_sample_line + '\n')
    
    # 提取指定区域的记录
    record_count = 0
    for line in f_in:
        if line.startswith('#'):
            continue
        
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 8:
            continue
        
        rec_chrom = parts[0]
        rec_pos = int(parts[1])
        
        # 检查染色体
        if rec_chrom != chrom and rec_chrom != chrom.replace('chr', ''):
            continue
        
        # 检查位置
        if rec_pos < start:
            continue
        if rec_pos > end:
            break  # 超出区域，停止（假设文件按位置排序）
        
        # 如果需要过滤样本
        if sample_ids and sample_indices != list(range(len(all_samples))):
            # 保留前9个字段 + 指定样本的数据
            selected_parts = parts[:9] + [parts[9 + i] for i in sample_indices if 9 + i < len(parts)]
            f_out.write('\t'.join(selected_parts) + '\n')
        else:
            # 保留所有样本
            f_out.write(line)
        
        record_count += 1
    
    f_in.close()
    f_out.close()
    
    print(f"[INFO] 子集VCF创建完成: {record_count} 个变异, {len(filtered_samples)} 个样本")


# ============================================================================
# 内置单倍型提取器（当主模块不可用时使用）
# ============================================================================

class BuiltinHaplotypeExtractor:
    """内置的单倍型提取器，用于当 haplotype_phenotype_analysis 模块不可用时"""
    
    def __init__(self, vcf_file: str):
        if not PYSAM_AVAILABLE:
            raise ImportError("pysam 不可用，无法提取单倍型")
        self.vcf_file = vcf_file
        self.vcf = pysam.VariantFile(vcf_file)
        self.samples = list(self.vcf.header.samples)
        print(f"[INFO] 内置提取器初始化成功，样本数: {len(self.samples)}")
    
    def extract_region(self, chrom: str, start: int, end: int, 
                       min_samples: int = 5, snp_only: bool = True) -> tuple:
        """
        提取指定区域的单倍型
        
        Args:
            snp_only: 是否只保留SNP（默认True，过滤indel和SV）
        
        Returns:
            tuple: (positions, hap_df, hap_sample_df)
        """
        # SNP判断函数
        def is_snp(ref_allele, alt_allele):
            valid_bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
            return (len(ref_allele) == 1 and len(alt_allele) == 1 and
                    ref_allele.upper() in valid_bases and alt_allele.upper() in valid_bases)
        
        positions = []
        genotypes_by_sample = {s: [] for s in self.samples}
        filtered_count = 0
        
        try:
            for record in self.vcf.fetch(chrom, start, end):
                ref = record.ref
                alts = record.alts if record.alts else []
                alt0 = alts[0] if alts else ""
                
                # SNP过滤
                if snp_only and not is_snp(ref, alt0):
                    filtered_count += 1
                    continue
                
                pos = record.pos
                positions.append(pos)
                
                for sample in self.samples:
                    gt = record.samples[sample]['GT']
                    if gt is None or None in gt:
                        allele = 'N'
                    else:
                        # 简化：取第一个等位基因
                        allele_idx = gt[0]
                        if allele_idx == 0:
                            allele = ref
                        elif allele_idx <= len(alts):
                            allele = alts[allele_idx - 1]
                        else:
                            allele = 'N'
                    genotypes_by_sample[sample].append(allele)
        except Exception as e:
            print(f"[WARNING] 提取区域 {chrom}:{start}-{end} 失败: {e}")
            return None, None, None
        
        if filtered_count > 0:
            print(f"[INFO] 过滤掉的非SNP变异数: {filtered_count}")
        
        if len(positions) == 0:
            print(f"[INFO] 区域 {chrom}:{start}-{end} 无变异")
            return positions, None, None
        
        # 过滤无变异位点（所有样本在该位点的等位基因相同）
        invariant_positions = []
        for i in range(len(positions)):
            alleles_at_pos = set()
            for sample in self.samples:
                if genotypes_by_sample[sample] and i < len(genotypes_by_sample[sample]):
                    allele = genotypes_by_sample[sample][i]
                    if allele != 'N':
                        alleles_at_pos.add(allele)
            if len(alleles_at_pos) <= 1:
                invariant_positions.append(i)
        
        if invariant_positions:
            print(f"[INFO] 过滤掉的无变异位点数: {len(invariant_positions)}")
            keep_indices = [i for i in range(len(positions)) if i not in invariant_positions]
            positions = [positions[i] for i in keep_indices]
            for sample in self.samples:
                if genotypes_by_sample[sample]:
                    genotypes_by_sample[sample] = [genotypes_by_sample[sample][i] for i in keep_indices]
            print(f"[INFO] 过滤后有效变异位点数: {len(positions)}")
        
        if len(positions) == 0:
            print(f"[INFO] 区域 {chrom}:{start}-{end} 过滤后无有效变异")
            return positions, None, None
        
        # 构建单倍型序列
        hap_seqs = {}
        sample_haps = []
        
        for sample in self.samples:
            alleles = genotypes_by_sample[sample]
            hap_seq = '|'.join(alleles)
            hap_seqs[sample] = hap_seq
            sample_haps.append({'SampleID': sample, 'Haplotype_Seq': hap_seq})
        
        # 统计单倍型频率并命名
        from collections import Counter
        seq_counts = Counter(hap_seqs.values())
        
        # 给频率排名前 N 的单倍型命名
        sorted_seqs = seq_counts.most_common()
        seq_to_name = {}
        for i, (seq, count) in enumerate(sorted_seqs):
            if count >= min_samples:
                seq_to_name[seq] = f"Hap{i+1}"
            else:
                seq_to_name[seq] = "Other"
        
        # 添加单倍型名称
        for item in sample_haps:
            item['Hap_Name'] = seq_to_name.get(item['Haplotype_Seq'], 'Other')
        
        hap_sample_df = pd.DataFrame(sample_haps)
        
        # 构建单倍型汇总表
        hap_list = []
        for seq, name in seq_to_name.items():
            if name != 'Other':
                alleles = seq.split('|')
                hap_list.append({
                    'Hap_Name': name,
                    'Haplotype_Seq': seq,
                    'Count': seq_counts[seq],
                    'Alleles': alleles
                })
        
        hap_df = pd.DataFrame(hap_list) if hap_list else None
        
        # 注意：不再基于主要单倍型进行第二次过滤
        # 因为次要单倍型中的变异可能对表型有显著影响
        
        print(f"[INFO] 提取完成: {len(positions)} 个变异, {len(hap_list)} 个单倍型")
        return positions, hap_df, hap_sample_df


# 根据模块可用性选择提取器
if HAPLOTYPE_MODULE_AVAILABLE:
    # 使用完整模块
    pass
else:
    # 使用内置提取器
    HaplotypeExtractor = BuiltinHaplotypeExtractor
    print("[INFO] 使用内置 HaplotypeExtractor")


# ============================================================================
# 配置
# ============================================================================

class ScanConfig:
    """扫描配置"""
    # ========== 原有数据路径（CS-IAAS T2T v1.1）备份 ==========
    # VCF_FILE = "/storage/public/home/2024110093/data/Variation/CSIAAS/Core819Samples_ALL.vcf.gz"
    # GFF_FILE = "/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1_HC.gff3"
    # PHENO_FILE = "/storage/public/home/2024110093/data/Variation/CSIAAS/Phe.txt"
    # FASTA_FILE = "/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1.fasta"
    
    # ========== 新数据路径（大麦 Barley Morex_v3）==========
    VCF_FILE = "/home/qinz/project/tmp_Proj/07.GWAS/20260103_Barley_salt/00.Rawdata/chrALL.impute.vcf.gz"
    GFF_FILE = "/home/qinz/data/genomes/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.gff3"
    PHENO_FILE = "/home/qinz/project/tmp_Proj/07.GWAS/20260103_Barley_salt/04.Gemma/00.pheno/K.rep1.SA.gemma"
    FASTA_FILE = "/home/qinz/data/genomes/Morex_v3/MorexV3_MtPt.fasta"
    GWAS_DIR = "/home/qinz/project/tmp_Proj/07.GWAS/20260103_Barley_salt"
    GWAS_RESULT_FILE = "/home/qinz/project/tmp_Proj/07.GWAS/20260103_Barley_salt/04.Gemma/04.Gemma/output/K.rep1.SA.assoc.txt"
    
    # 输出路径
    DATABASE_DIR = "./database"           # 数据库文件夹（存放提取的数据）
    RESULTS_DIR = "./results"             # 结果文件夹（存放HTML图）
    
    # 指定要分析的基因列表（大麦 GWAS 及 SV 筛选基因）
    TARGET_GENES = [
        "HORVU.MOREX.r3.4HG0405190",  # myosin-like protein XIF
        "HORVU.MOREX.r3.4HG0405240",  # Jasmonate ZIM domain protein
        "HORVU.MOREX.r3.4HG0405310",  # Protein ROOT PRIMORDIUM DEFECTIVE 1
        "HORVU.MOREX.r3.4HG0405230",  # Jasmonate ZIM domain protein
        "HORVU.MOREX.r3.5HG0426750",  # Cation/H(+) antiporter
        "HORVU.MOREX.r3.3HG0281790",  # Laccase
        "HORVU.MOREX.r3.3HG0282160",  # Serine/threonine-protein kinase
        "HORVU.MOREX.r3.3HG0282170",  # Methyltransferase
        "HORVU.MOREX.r3.4HG0387570",  # 70 kDa heat shock protein
        "HORVU.MOREX.r3.7HG0634850"   # Acyl-coenzyme A oxidase
    ]
    
    # 启动子区域参数
    PROMOTER_LENGTH = 2000  # 启动子长度（bp），在TSS上游扩展
    
    # 变异类型过滤
    SNP_ONLY = False  # False=包含所有变异类型(SNP/indel/SV)
    
    # 过滤参数
    MIN_VARIANTS = 1  # 最小变异数
    MIN_SAMPLES = 1   # 最小样本数（保留所有单倍型，不过滤）
    
    # 分析参数
    PVALUE_THRESHOLD = 0.05  # P值显著性阈值
    
    # 测试模式：限制区间长度（设为0表示不限制，使用完整基因区间）
    TEST_REGION_LENGTH = 0  # 0=使用完整基因区间


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
    
    # 记录数据格式信息
    if len(df) > 0:
        print("\n" + "=" * 60)
        print("[DATA FORMAT] GFF3 基因格式")
        print("=" * 60)
        print(f"  基因ID示例 (前10个):")
        for i, gene_id in enumerate(df['gene_id'].head(10)):
            print(f"    [{i+1}] {gene_id}")
        print(f"\n  染色体列表 (前10个): {list(df['chrom'].unique()[:10])}")
        print(f"  总染色体数: {df['chrom'].nunique()}")
        print(f"  基因数量: {len(df)}")
        print(f"  ID格式特征: 长度范围={df['gene_id'].str.len().min()}-{df['gene_id'].str.len().max()}, 前缀示例={df['gene_id'].iloc[0][:5] if len(df) > 0 else 'N/A'}")
        print("=" * 60 + "\n")
    
    return df


# ============================================================================
# 单基因处理函数（用于并行）
# ============================================================================

def parse_gtf_for_gene(gtf_file: str, target_gene_id: str) -> dict:
    """
    从 GTF/GFF3 文件解析目标基因的外显子和 CDS 位置
    简化版，用于数据库建立时保存基因结构信息
    """
    if not gtf_file or not os.path.exists(gtf_file):
        return {'exons': [], 'cds': [], 'strand': '+', 'chrom': None, 'gene_start': None, 'gene_end': None}
    
    exons = []
    cds = []
    strand = '+'
    chrom = None
    gene_start = None
    gene_end = None
    
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    is_gff3 = gtf_file.endswith('.gff3') or gtf_file.endswith('.gff')
    
    try:
        with open_func(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                chrom = parts[0]
                attr = parts[8]
                
                # 检查是否匹配目标基因
                gene_match = False
                if is_gff3:
                    # GFF3: ID=GENEID 或 Parent=GENEID
                    if f'ID={target_gene_id}' in attr or f'Parent={target_gene_id}' in attr:
                        gene_match = True
                else:
                    # GTF: gene_id "GENEID"
                    if f'gene_id "{target_gene_id}"' in attr:
                        gene_match = True
                
                if not gene_match:
                    continue
                
                # 记录基因坐标
                if feature_type == 'gene':
                    gene_start = start
                    gene_end = end
                elif feature_type == 'exon':
                    exons.append((start, end))
                elif feature_type == 'CDS':
                    cds.append((start, end))
        
        return {
            'exons': sorted(exons),
            'cds': sorted(cds),
            'strand': strand,
            'chrom': chrom,
            'gene_start': gene_start,
            'gene_end': gene_end
        }
    except Exception as e:
        print(f"[WARNING] 解析GTF失败: {e}")
        return {'exons': [], 'cds': [], 'strand': '+', 'chrom': None, 'gene_start': None, 'gene_end': None}


def process_single_gene(gene_info: dict, vcf_file: str, pheno_df: pd.DataFrame, 
                        database_dir: str, results_dir: str = None,
                        min_samples: int = 5, test_region_length: int = 0,
                        gff_file: str = None) -> dict:
    """
    处理单个基因，提取单倍型并保存到专属文件夹
    
    文件夹结构: 
        database/{gene_id}/     # 数据文件
        results/{gene_id}/      # HTML图表
    
    Args:
        test_region_length: 测试模式下的区间长度，0表示使用完整区间
    
    Returns:
        dict: 基因处理结果摘要
    """
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    gene_start = gene_info['start']  # GFF中的基因体起始
    gene_end = gene_info['end']      # GFF中的基因体终止
    strand = gene_info['strand']
    
    # 扩展启动子区域
    promoter_length = ScanConfig.PROMOTER_LENGTH
    if strand == '+':  # +链：TSS在gene_start，启动子在上游(更小坐标)
        start = max(1, gene_start - promoter_length)  # 向上游扩展
        end = gene_end
    else:  # -链：TSS在gene_end，启动子在下游(更大坐标)
        start = gene_start
        end = gene_end + promoter_length  # 向下游扩展
    
    print(f"[INFO] {gene_id}: 基因体={gene_start}-{gene_end}({strand}), 启动子扩展后={start}-{end}")
    
    # 测试模式：限制区间长度
    original_end = end
    if test_region_length > 0:
        end = min(start + test_region_length, original_end)
        print(f"[TEST] 测试模式: 区间限制为 {start}-{end} ({test_region_length}bp)")
    
    # 创建基因专属文件夹（数据库目录）
    gene_folder_name = gene_id  # 简化文件夹名，只用gene_id
    gene_data_dir = os.path.join(database_dir, gene_folder_name)
    os.makedirs(gene_data_dir, exist_ok=True)
    
    # 如果指定了results_dir，也创建对应目录
    gene_results_dir = None
    if results_dir:
        gene_results_dir = os.path.join(results_dir, gene_folder_name)
        os.makedirs(gene_results_dir, exist_ok=True)
    
    result_summary = {
        'gene_id': gene_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'n_variants': 0,
        'n_haplotypes': 0,
        'n_samples': 0,
        'data_folder': gene_folder_name,
        'results_folder': gene_folder_name if results_dir else None,
        'status': 'pending'
    }
    
    try:
        # 创建单倍型提取器
        extractor = HaplotypeExtractor(vcf_file)
        
        # 提取单倍型
        positions, hap_df, hap_sample_df = extractor.extract_region(
            chrom, start, end, min_samples=min_samples, snp_only=ScanConfig.SNP_ONLY
        )
        
        # 解析GTF获取基因结构信息
        gtf_data = parse_gtf_for_gene(gff_file, gene_id) if gff_file else {'exons': [], 'cds': [], 'strand': strand}
        
        # 保存基因基本信息
        gene_info_dict = {
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,  # 启动子扩展后的起始
            'end': end,      # 启动子扩展后的终止
            'gene_start': gene_start,  # 原始基因体起始
            'gene_end': gene_end,      # 原始基因体终止
            'strand': strand,
            'length': end - start + 1,
            'promoter_length': promoter_length,
            'vcf_file': vcf_file,
            'vcf_mtime': os.path.getmtime(vcf_file),  # 关键：保存VCF修改时间用于缓存判断
            'extraction_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'exons': gtf_data.get('exons', []),
            'cds': gtf_data.get('cds', [])
        }
        with open(os.path.join(gene_data_dir, 'gene_info.json'), 'w') as f:
            json.dump(gene_info_dict, f, indent=2)
        
        if hap_df is None or len(hap_df) == 0:
            print(f"[WARNING] hap_df 为空，但 hap_sample_df 有 {len(hap_sample_df) if hap_sample_df is not None else 0} 行")
            result_summary['status'] = 'no_variants'
            return result_summary
        
        # 关键检查：如果 hap_sample_df 为空，直接返回
        if hap_sample_df is None or len(hap_sample_df) == 0:
            print(f"[WARNING] hap_sample_df 为空，无法继续")
            result_summary['status'] = 'no_valid_haplotypes'
            return result_summary
        
        n_variants = len(positions) if positions else 0
        result_summary['n_variants'] = n_variants
        
        # 保存单倍型数据
        hap_df.to_csv(os.path.join(gene_data_dir, 'haplotype_data.csv'), index=False)
        
        # [DATA FORMAT LOG] 单倍型数据格式
        print(f"\n[DATA FORMAT] haplotype_data.csv")
        print(f"  行数: {len(hap_df)}, 列数: {len(hap_df.columns)}")
        print(f"  列名: {list(hap_df.columns)}")
        print(f"  列数据类型: {dict(hap_df.dtypes)}")
        print(f"  示例数据 (前3行):")
        print(hap_df.head(3).to_string(index=False))
        print()
        
        # 保存变异注释信息（从extractor获取）
        if hasattr(extractor, 'variant_info') and extractor.variant_info:
            variant_info_df = pd.DataFrame([
                {
                    'position': pos,
                    'ref': info.get('ref', ''),
                    'alt': info.get('alt', ''),
                    'len_diff': info.get('len_diff', 0),
                    'is_sv': info.get('is_sv', False),
                    'maf': info.get('maf', 0.5),
                    'missing_rate': info.get('missing_rate', 0.0)
                }
                for pos, info in extractor.variant_info.items()
            ])
            variant_info_df.to_csv(os.path.join(gene_data_dir, 'variant_info.csv'), index=False)
            print(f"[INFO] 变异注释信息已保存: {len(variant_info_df)} 个位点")
        
        # **新增**: 分析并保存启动子区域变异信息
        # 计算启动子区域（基因体外的上游区域）
        if strand == '+':
            promoter_start = max(1, gene_start - promoter_length)
            promoter_end = gene_start - 1
        else:
            promoter_start = gene_end + 1
            promoter_end = gene_end + promoter_length
        
        # **关键修复**：加载所有基因的CDS信息，用于检查启动子变异是否与其他基因重叠
        all_genes_cds = []  # [(chrom, cds_start, cds_end, gene_id), ...]
        if gff_file and os.path.exists(gff_file):
            open_func_gff = gzip.open if gff_file.endswith('.gz') else open
            is_gff3_format = gff_file.endswith('.gff3') or gff_file.endswith('.gff')
            current_gene_id = None
            
            try:
                with open_func_gff(gff_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) < 9:
                            continue
                        
                        feature_type = parts[2]
                        gff_chrom = parts[0]
                        gff_start = int(parts[3])
                        gff_end = int(parts[4])
                        attr = parts[8]
                        
                        # 解析基因ID
                        if is_gff3_format:
                            if feature_type == 'gene':
                                for item in attr.split(';'):
                                    if item.startswith('ID='):
                                        current_gene_id = item[3:]
                                        break
                            elif feature_type == 'CDS' and current_gene_id:
                                for item in attr.split(';'):
                                    if item.startswith('Parent='):
                                        parent_id = item[7:]
                                        # 如果Parent是当前关注的基因，跳过（不检查自己）
                                        if parent_id == gene_id:
                                            break
                                        all_genes_cds.append((gff_chrom, gff_start, gff_end, parent_id))
                                        break
                        else:
                            # GTF格式
                            if feature_type == 'gene':
                                for item in attr.split(';'):
                                    item = item.strip()
                                    if item.startswith('gene_id "'):
                                        current_gene_id = item[9:-1]
                                        break
                            elif feature_type == 'CDS':
                                for item in attr.split(';'):
                                    item = item.strip()
                                    if item.startswith('gene_id "'):
                                        cds_gene_id = item[9:-1]
                                        if cds_gene_id != gene_id:  # 不检查自己
                                            all_genes_cds.append((gff_chrom, gff_start, gff_end, cds_gene_id))
                                        break
            except Exception as e:
                print(f"[WARNING] 加载其他基因CDS信息失败: {e}")
            
            print(f"[INFO] 加载了 {len(all_genes_cds)} 个其他基因的CDS区域（用于重叠检查）")
        
        # 提取启动子区域的变异
        promoter_variants = []
        if positions:
            for pos in positions:
                if promoter_start <= pos <= promoter_end:
                    info = extractor.variant_info.get(pos, {}) if hasattr(extractor, 'variant_info') else {}
                    
                    # **关键检查**：该位置是否在其他基因的CDS区域内
                    overlaps_cds = False
                    overlapping_genes = []
                    
                    for cds_chrom, cds_start, cds_end, cds_gene_id in all_genes_cds:
                        if cds_chrom == chrom and cds_start <= pos <= cds_end:
                            overlaps_cds = True
                            overlapping_genes.append(cds_gene_id)
                    
                    # 如果与CDS重叠，标记但保留（让用户知道）
                    if overlaps_cds:
                        promoter_variants.append({
                            'position': pos,
                            'ref': info.get('ref', ''),
                            'alt': info.get('alt', ''),
                            'distance_to_tss': abs(pos - gene_start) if strand == '+' else abs(pos - gene_end),
                            'is_sv': info.get('is_sv', False),
                            'maf': info.get('maf', 0.5),
                            'overlaps_cds': True,
                            'overlapping_genes': ';'.join(overlapping_genes)
                        })
                    else:
                        promoter_variants.append({
                            'position': pos,
                            'ref': info.get('ref', ''),
                            'alt': info.get('alt', ''),
                            'distance_to_tss': abs(pos - gene_start) if strand == '+' else abs(pos - gene_end),
                            'is_sv': info.get('is_sv', False),
                            'maf': info.get('maf', 0.5),
                            'overlaps_cds': False,
                            'overlapping_genes': ''
                        })
        
        # 保存启动子变异信息
        if promoter_variants:
            promoter_df = pd.DataFrame(promoter_variants)
            promoter_df.to_csv(os.path.join(gene_data_dir, 'promoter_variants.csv'), index=False)
            print(f"[INFO] 启动子变异信息已保存: {len(promoter_df)} 个变异")
            
            # 统计CDS重叠的变异数量
            n_overlaps_cds = sum(1 for v in promoter_variants if v.get('overlaps_cds', False))
            n_pure_promoter = len(promoter_variants) - n_overlaps_cds
            
            # 同时保存详细文本报告
            with open(os.path.join(gene_data_dir, 'promoter_variants_detail.txt'), 'w') as f:
                f.write(f"启动子区域变异分析报告\n")
                f.write(f"=" * 60 + "\n")
                f.write(f"基因: {gene_id}\n")
                f.write(f"染色体: {chrom}\n")
                f.write(f"链方向: {strand}\n")
                f.write(f"启动子区域: {promoter_start:,}-{promoter_end:,}\n")
                f.write(f"TSS位置: {gene_start:,} (正链) 或 {gene_end:,} (负链)\n")
                f.write(f"变异总数: {len(promoter_variants)}\n")
                f.write(f"  - 纯启动子变异: {n_pure_promoter}\n")
                f.write(f"  - 与其他基因CDS重叠: {n_overlaps_cds}\n")
                f.write(f"\n变异详情:\n")
                f.write(f"-" * 60 + "\n")
                for v in promoter_variants:
                    sv_marker = " [SV]" if v['is_sv'] else ""
                    cds_marker = ""
                    if v.get('overlaps_cds', False):
                        cds_marker = f" [CDS重叠: {v.get('overlapping_genes', 'N/A')}]"
                    f.write(f"  位置: {v['position']:,} | 距离TSS: {v['distance_to_tss']:,}bp | "
                           f"REF: {v['ref']} | ALT: {v['alt']} | MAF: {v['maf']:.3f}{sv_marker}{cds_marker}\n")
        else:
            print(f"[INFO] 启动子区域未发现变异")
        
        # 合并表型数据
        print(f"\n[DATA FORMAT] 样本-单倍型映射表 (hap_sample_df)")
        print(f"  行数: {len(hap_sample_df)}, 列数: {len(hap_sample_df.columns)}")
        print(f"  列名: {list(hap_sample_df.columns)}")
        if len(hap_sample_df) > 0 and 'SampleID' in hap_sample_df.columns:
            print(f"  SampleID列数据类型: {hap_sample_df['SampleID'].dtype}")
            print(f"  SampleID示例 (前5个): {hap_sample_df['SampleID'].head(5).tolist()}")
        print()
        
        print(f"[DATA FORMAT] 表型数据 (pheno_df)")
        print(f"  行数: {len(pheno_df)}, 列数: {len(pheno_df.columns)}")
        print(f"  列名: {list(pheno_df.columns)}")
        if len(pheno_df) > 0:
            print(f"  SampleID列数据类型: {pheno_df['SampleID'].dtype if 'SampleID' in pheno_df.columns else 'N/A'}")
            print(f"  SampleID示例 (前5个): {pheno_df['SampleID'].head(5).tolist() if 'SampleID' in pheno_df.columns else 'N/A'}")
        print()
        
        if 'SampleID' not in hap_sample_df.columns:
            print(f"[ERROR] hap_sample_df 缺少 'SampleID' 列！")
            print(f"[ERROR] 实际列名: {list(hap_sample_df.columns)}")
            result_summary['status'] = 'no_sample_id'
            return result_summary
        
        # 保存样本-单倍型对应表
        hap_sample_df.to_csv(os.path.join(gene_data_dir, 'haplotype_samples.csv'), index=False)
        
        # **保存完整VCF子集**（与原始数据完全一致）
        # 注意：无索引时会全文件扫描，耗时较长但数据完整
        subset_vcf_path = os.path.join(gene_data_dir, 'variants.vcf.gz')
        if not os.path.exists(subset_vcf_path) or os.path.getsize(subset_vcf_path) < 1024:
            try:
                # 获取所有在hap_sample_df中的样本（这些样本在VCF中有基因型数据）
                vcf_sample_ids = hap_sample_df['SampleID'].tolist()
                print(f"[INFO] 开始保存VCF子集: {subset_vcf_path}")
                print(f"[INFO] VCF区域: {chrom}:{start}-{end}, 样本数: {len(vcf_sample_ids)}")
                
                # 保存VCF子集（包含所有样本的基因型数据）
                if not vcf_sample_ids:
                    print(f"[WARNING] 样本ID列表为空，保存所有样本")
                    create_subset_vcf(vcf_file, chrom, start, end, subset_vcf_path, sample_ids=None)
                else:
                    create_subset_vcf(vcf_file, chrom, start, end, subset_vcf_path, 
                                    sample_ids=vcf_sample_ids)
                
                # 验证 VCF 文件是否成功生成
                vcf_size = os.path.getsize(subset_vcf_path)
                if vcf_size > 1024:
                    print(f"[INFO] VCF子集已保存: {subset_vcf_path} ({vcf_size/1024:.1f} KB)")
                else:
                    print(f"[WARNING] VCF子集文件过小 ({vcf_size} bytes)，可能提取失败")
            except Exception as e:
                import traceback
                print(f"[WARNING] 保存VCF子集失败: {e}")
                print(f"[WARNING] 错误详情: {traceback.format_exc()}")
        else:
            print(f"[INFO] VCF 子集已存在: {subset_vcf_path} ({os.path.getsize(subset_vcf_path)/1024:.1f} KB)")
        
        merged = pd.merge(hap_sample_df, pheno_df, on='SampleID', how='inner')
        print(f"[INFO] 表型合并完成: merged 行数={len(merged)}")
        
        print(f"[DATA FORMAT] 合并后的表型数据 (merged)")
        print(f"  行数: {len(merged)} (hap_sample_df: {len(hap_sample_df)}, pheno_df: {len(pheno_df)})")
        print(f"  列数: {len(merged.columns)}")
        print(f"  列名: {list(merged.columns)}")
        if len(merged) > 0:
            print(f"  合并成功！样本匹配率: {len(merged)/len(hap_sample_df)*100:.1f}%")
            print(f"  示例数据 (前3行):")
            print(merged.head(3).to_string(index=False))
        else:
            print(f"[WARNING] 合并后为空！检查 SampleID 匹配情况：")
            hap_samples = set(hap_sample_df['SampleID'].astype(str))
            pheno_samples = set(pheno_df['SampleID'].astype(str))
            print(f"  hap_sample_df 样本数: {len(hap_samples)}")
            print(f"  pheno_df 样本数: {len(pheno_samples)}")
            print(f"  交集样本数: {len(hap_samples & pheno_samples)}")
            if len(hap_samples & pheno_samples) == 0:
                print(f"  [ERROR] 完全没有匹配的样本！")
                print(f"  hap_sample_df SampleID 示例: {list(hap_samples)[:5]}")
                print(f"  pheno_df SampleID 示例: {list(pheno_samples)[:5]}")
        print()
        
        if len(merged) == 0:
            print(f"[WARNING] 表型数据匹配失败，保存已有数据后退出")
            # 即使表型不匹配，也保存单倍型统计（不含表型）
            print(f"[INFO] 开始保存 haplotype_stats.csv（无表型）...")
            hap_stats_list = []
            hap_sample_counts = hap_sample_df['Hap_Name'].value_counts()
            for hap_name, count in hap_sample_counts.items():
                hap_stats_list.append({
                    'haplotype_name': hap_name,
                    'haplotype_count': count,
                    'haplotype_freq': round(count / len(hap_sample_df), 4),
                    'phenotype_mean': None,
                    'phenotype_sd': None
                })
            hap_stats_df = pd.DataFrame(hap_stats_list)
            hap_stats_df.to_csv(os.path.join(gene_data_dir, 'haplotype_stats.csv'), index=False)
            print(f"[INFO] 已保存 haplotype_stats.csv（无表型数据）")
            print(f"[INFO] 返回结果: no_phenotype_match")
            
            result_summary['status'] = 'no_phenotype_match'
            return result_summary
        
        # 保存合并后的表型数据
        print(f"[INFO] 保存 phenotype_data.csv...")
        merged.to_csv(os.path.join(gene_data_dir, 'phenotype_data.csv'), index=False)
        print(f"[INFO] phenotype_data.csv 已保存")
        
        # 获取表型列
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
        
        result_summary['n_haplotypes'] = n_haplotypes
        result_summary['n_samples'] = total_samples
        
        # 保存单倍型统计
        hap_stats_list = []
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
            
            hap_stats_list.append({
                'haplotype_name': hap_name,
                'haplotype_count': hap_count,
                'haplotype_freq': round(hap_freq, 4),
                'phenotype_mean': round(pheno_mean, 4) if pheno_mean is not None else None,
                'phenotype_sd': round(pheno_sd, 4) if pheno_sd is not None else None
            })
        
        hap_stats_df = pd.DataFrame(hap_stats_list)
        hap_stats_df.to_csv(os.path.join(gene_data_dir, 'haplotype_stats.csv'), index=False)
        print(f"[INFO] haplotype_stats.csv 已保存 ({n_haplotypes} 个单倍型)")
        
        result_summary['status'] = 'success'
        print(f"[INFO] 基因处理完成: {gene_id}, status=success")
        return result_summary
        
    except Exception as e:
        import traceback
        error_msg = str(e) + " | " + traceback.format_exc().replace('\n', ' ')
        result_summary['status'] = f'error: {error_msg[:200]}'
        return result_summary


def analyze_gene_association(gene_info: dict, vcf_file: str, pheno_df: pd.DataFrame,
                            min_samples: int = 5, pvalue_threshold: float = 0.05,
                            gff_file: str = None, results_dir: str = None,
                            cluster_haplotypes: bool = False) -> dict:
    """
    对单个基因进行单倍型-表型关联分析，并生成综HTML图
    
    Args:
        results_dir: HTML图表输出目录
        cluster_haplotypes: 是否按序列相似度聚类排序单倍型
    
    Returns:
        dict: 包含关联分析结果和HTML路径的字典
    """
    gene_id = gene_info['gene_id']
    chrom = gene_info['chrom']
    gene_start = gene_info['start']  # GFF中的基因体起始
    gene_end = gene_info['end']      # GFF中的基因体终止
    strand = gene_info['strand']
    
    # 扩展启动子区域
    promoter_length = ScanConfig.PROMOTER_LENGTH
    if strand == '+':  # +链：TSS在gene_start，启动子在上游
        start = max(1, gene_start - promoter_length)
        end = gene_end
    else:  # -链：TSS在gene_end，启动子在下游
        start = gene_start
        end = gene_end + promoter_length
    
    result = {
        'gene_id': gene_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'gene_start': gene_start,  # 保留原始基因体坐标
        'gene_end': gene_end,
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
            chrom, start, end, min_samples=min_samples, snp_only=ScanConfig.SNP_ONLY
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
        
        # 如果显著，生成综HTML图
        if result['significant'] and gff_file and results_dir:
            try:
                # 创建输出目录（results文件夹下）
                gene_output_dir = os.path.join(results_dir, gene_id)
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
                    chrom=chrom,
                    start=start,
                    end=end,
                    gene_id=gene_id,
                    phenotype_cols=[pheno_col],
                    cluster_haplotypes=cluster_haplotypes
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
                    database_dir: str, results_dir: str = None,
                    n_workers: int = 1, chrom_filter: str = None, 
                    gene_filter: list = None, batch_size: int = 100, 
                    min_samples: int = 5, run_analysis: bool = True, 
                    pvalue_threshold: float = 0.05,
                    generate_html: bool = True,
                    test_region_length: int = 0,
                    cluster_haplotypes: bool = False) -> pd.DataFrame:
    """
    运行指定基因的单倍型数据集构建
    
    目录结构:
        database/                           # 数据库文件夹
        ├── {gene_id}/
        │   ├── gene_info.json              # 基因基本信息
        │   ├── haplotype_data.csv          # 单倍型数据
        │   ├── haplotype_samples.csv       # 样本-单倍型对应
        │   ├── haplotype_stats.csv         # 单倍型统计
        │   ├── phenotype_data.csv          # 表型数据
        │   └── association_result.csv      # 关联分析结果
        └── summary.csv                      # 所有基因汇总
        
        results/                            # 结果文件夹
        └── {gene_id}/
            └── integrated_analysis.html    # 综合分析HTML图
    
    Args:
        vcf_file: VCF 文件路径
        gff_file: GFF3 文件路径
        pheno_file: 表型文件路径
        database_dir: 数据库输出目录
        results_dir: 结果(图表)输出目录
        n_workers: 并行进程数
        chrom_filter: 染色体过滤（如 'chr1A'）
        gene_filter: 基因ID过滤列表（默认使用TARGET_GENES）
        batch_size: 批处理大小
        min_samples: 最小样本数
        run_analysis: 是否进行关联分析
        pvalue_threshold: P值显著性阈值
        generate_html: 是否生成HTML报告
        cluster_haplotypes: 是否按序列相似度聚类排序单倍型
    
    Returns:
        DataFrame: 扫描结果总表
    """
    # 如果未指定基因过滤，使用默认的目标基因
    if gene_filter is None:
        gene_filter = ScanConfig.TARGET_GENES
    
    os.makedirs(database_dir, exist_ok=True)
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
    
    # 性能监控
    perf = PerformanceMonitor()
    perf.start()
    
    print("=" * 60)
    print("单倍型数据集构建")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"VCF: {vcf_file}")
    print(f"GFF: {gff_file}")
    print(f"表型: {pheno_file}")
    print(f"数据库目录: {database_dir}")
    print(f"结果目录: {results_dir}")
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
    
    # 自动检测表型文件格式
    print("\n[INFO] 检测表型文件格式...")
    with open(pheno_file, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
    
    # 检测分隔符
    if '\t' in first_line:
        separator = '\t'
        sep_name = '制表符(\\t)'
    elif ',' in first_line:
        separator = ','
        sep_name = '逗号(,)'
    elif ' ' in first_line:
        separator = r'\s+'
        sep_name = '空格'
    else:
        separator = '\t'  # 默认
        sep_name = '制表符(\\t)（默认）'
    
    print(f"  检测到分隔符: {sep_name}")
    
    # 尝试读取表型文件
    try:
        pheno_df = pd.read_csv(pheno_file, sep=separator)
        print(f"  列数: {len(pheno_df.columns)}")
        print(f"  行数: {len(pheno_df)}")
        print(f"  列名: {list(pheno_df.columns)}")
        
        # 检测第一列是否为样本 ID
        first_col = pheno_df.columns[0]
        first_col_sample = pheno_df[first_col].iloc[0]
        
        # 判断第一列是字符串（样本ID）还是数值（表型值）
        try:
            float(first_col_sample)
            is_numeric = True
        except (ValueError, TypeError):
            is_numeric = False
        
        if is_numeric:
            # 第一列是数值，说明没有样本 ID 列，需要从 VCF 读取真实样本ID
            print(f"\n  ⚠️  检测到表型文件无样本ID列（第一列为数值）")
            print(f"  第一列示例值: {first_col_sample}")
            
            # 关键修复：从 VCF 读取真实样本ID
            vcf_samples = []
            if vcf_file and os.path.exists(vcf_file):
                try:
                    import pysam
                    vcf = pysam.VariantFile(vcf_file)
                    vcf_samples = list(vcf.header.samples)
                    vcf.close()
                    print(f"  从 VCF 读取到 {len(vcf_samples)} 个真实样本ID")
                    print(f"  VCF 样本ID示例: {vcf_samples[:5]}")
                except Exception as e:
                    print(f"  [WARNING] 读取 VCF 样本失败: {e}")
            
            # 重新读取，所有列都作为表型数据
            pheno_df = pd.read_csv(pheno_file, sep=separator, header=None)
            
            # 如果第一行是列名（尝试转换为浮点数失败）
            try:
                float(pheno_df.iloc[0, 0])
                # 第一行就是数据，没有表头
                pheno_df.columns = [f'Phenotype_{i+1}' for i in range(len(pheno_df.columns))]
                print(f"  检测到：无表头，所有列为表型数据")
            except (ValueError, TypeError):
                # 第一行是表头
                pheno_df.columns = pheno_df.iloc[0]
                pheno_df = pheno_df.iloc[1:].reset_index(drop=True)
                print(f"  检测到：有表头，列名: {list(pheno_df.columns)}")
            
            # 关键修复：使用 VCF 的真实样本ID
            if len(vcf_samples) >= len(pheno_df):
                pheno_df.insert(0, 'SampleID', vcf_samples[:len(pheno_df)])
                print(f"  已添加 SampleID 列（从 VCF 读取，示例: {vcf_samples[:3]}）")
            elif len(vcf_samples) > 0:
                print(f"  [WARNING] VCF 样本数({len(vcf_samples)}) < 表型行数({len(pheno_df)})")
                pheno_df.insert(0, 'SampleID', vcf_samples + [f'Unknown_{i}' for i in range(len(pheno_df) - len(vcf_samples))])
            else:
                print(f"  [WARNING] 无法读取 VCF 样本，使用自动生成ID")
                pheno_df.insert(0, 'SampleID', [f'Sample_{i+1:03d}' for i in range(len(pheno_df))])
            
            
        else:
            # 第一列是字符串，可能是样本 ID
            print(f"\n  ✓ 检测到样本ID列: '{first_col}'")
            print(f"  样本ID示例: {first_col_sample}")
            
            # 标准化列名
            if 'SampleID' not in pheno_df.columns:
                if pheno_df.columns[0] not in ['SampleID', 'sample_id', 'ID']:
                    pheno_df.columns = ['SampleID'] + list(pheno_df.columns[1:])
                    print(f"  已将第一列重命名为 'SampleID'")
                else:
                    pheno_df = pheno_df.rename(columns={pheno_df.columns[0]: 'SampleID'})
    
    except Exception as e:
        print(f"\n  ✗ 表型文件读取失败: {e}")
        print(f"  尝试使用空白分隔符读取...")
        pheno_df = pd.read_csv(pheno_file, sep=r'\s+', engine='python')
        if 'SampleID' not in pheno_df.columns:
            if pheno_df.columns[0] not in ['SampleID', 'sample_id', 'ID']:
                pheno_df.columns = ['SampleID'] + list(pheno_df.columns[1:])
            else:
                pheno_df = pheno_df.rename(columns={pheno_df.columns[0]: 'SampleID'})
    
    print(f"\n[INFO] 表型样本数: {len(pheno_df)}, 列: {list(pheno_df.columns)}")
    
    # 记录表型数据格式
    print("\n" + "=" * 60)
    print("[DATA FORMAT] 表型数据格式")
    print("=" * 60)
    print(f"  样本ID列名: 'SampleID'")
    print(f"  样本ID示例 (前5个): {list(pheno_df['SampleID'].head(5))}")
    print(f"  样本ID数据类型: {pheno_df['SampleID'].dtype}")
    print(f"  表型列: {[c for c in pheno_df.columns if c != 'SampleID']}")
    print(f"  总样本数: {len(pheno_df)}")
    
    # 显示表型列的统计信息
    for col in pheno_df.columns:
        if col != 'SampleID':
            print(f"\n  表型列 '{col}' 详情:")
            print(f"    数据类型: {pheno_df[col].dtype}")
            print(f"    缺失值: {pheno_df[col].isna().sum()}")
            if pheno_df[col].dtype in ['float64', 'int64']:
                print(f"    范围: {pheno_df[col].min():.4f} - {pheno_df[col].max():.4f}")
                print(f"    均值: {pheno_df[col].mean():.4f}")
                print(f"    标准差: {pheno_df[col].std():.4f}")
            # 显示前5个值作为示例
            print(f"    示例值 (前5个): {list(pheno_df[col].head(5))}")
    print("=" * 60 + "\n")
    perf.step_end("load_phenotype")
    
    # 3. 批量处理
    perf.step_start("scan_genes")
    all_results = []
    processed = 0
    errors = 0
    
    gene_list = genes_df.to_dict('records')
    
    # 使用进度条
    print(f"\n[INFO] 开始处理 {total_genes} 个基因...")
    print(f"[INFO] 数据库输出目录: {database_dir}")
    if results_dir:
        print(f"[INFO] 结果输出目录: {results_dir}")
    if test_region_length > 0:
        print(f"[TEST] 测试模式已启用，区间限制为 {test_region_length}bp")
    
    # 串行处理（更稳定）
    for i, gene_info in enumerate(gene_list):
        gene_id = gene_info['gene_id']
        gene_data_dir = os.path.join(database_dir, gene_id)
        
        # 检查数据库缓存是否存在
        if os.path.exists(gene_data_dir) and os.path.exists(os.path.join(gene_data_dir, 'gene_info.json')):
            # 检查 VCF 是否被修改
            vcf_mtime = os.path.getmtime(vcf_file)
            cached_gene_info_file = os.path.join(gene_data_dir, 'gene_info.json')
            try:
                with open(cached_gene_info_file, 'r') as f:
                    cached_info = json.load(f)
                cached_vcf_mtime = cached_info.get('vcf_mtime', 0)
                
                # 关键修复：如果旧缓存没有 vcf_mtime，自动补充
                if cached_vcf_mtime == 0:
                    print(f"  [{processed+1}/{total_genes}] {gene_id}: [CACHE FIX] 旧缓存缺少 vcf_mtime，自动补充")
                    cached_info['vcf_mtime'] = vcf_mtime
                    with open(cached_gene_info_file, 'w') as f:
                        json.dump(cached_info, f, indent=2)
                    cached_vcf_mtime = vcf_mtime  # 设置为当前时间，下次就能命中缓存
                
                if vcf_mtime <= cached_vcf_mtime:
                    # 关键：验证缓存文件完整性
                    required_files = ['gene_info.json', 'haplotype_data.csv', 'haplotype_samples.csv', 'variant_info.csv']
                    missing_files = [f for f in required_files if not os.path.exists(os.path.join(gene_data_dir, f))]
                    
                    if missing_files:
                        print(f"  [{processed+1}/{total_genes}] {gene_id}: [CACHE INCOMPLETE] 缺少文件 {missing_files}，重新处理")
                    else:
                        print(f"  [{processed+1}/{total_genes}] {gene_id}: [CACHED] 从数据库加载")
                        # 从缓存加载统计信息
                        hap_stats_file = os.path.join(gene_data_dir, 'haplotype_stats.csv')
                        if os.path.exists(hap_stats_file):
                            hap_stats = pd.read_csv(hap_stats_file)
                            # 安全检查：确保列存在
                            n_samples = int(hap_stats['haplotype_count'].sum()) if 'haplotype_count' in hap_stats.columns else 0
                            result_summary = {
                                'gene_id': gene_id,
                                'chrom': cached_info.get('chrom'),
                                'start': cached_info.get('start'),
                                'end': cached_info.get('end'),
                                'strand': cached_info.get('strand'),
                                'n_variants': len(pd.read_csv(os.path.join(gene_data_dir, 'variant_info.csv'))) if os.path.exists(os.path.join(gene_data_dir, 'variant_info.csv')) else 0,
                                'n_haplotypes': len(hap_stats),
                                'n_samples': n_samples,
                                'data_folder': gene_id,
                                'results_folder': gene_id if results_dir else None,
                                'status': 'success'
                            }
                        else:
                            result_summary = {
                                'gene_id': gene_id,
                                'chrom': cached_info.get('chrom'),
                                'start': cached_info.get('start'),
                                'end': cached_info.get('end'),
                                'strand': cached_info.get('strand'),
                                'n_variants': 0,
                                'n_haplotypes': 0,
                                'n_samples': 0,
                                'data_folder': gene_id,
                                'results_folder': gene_id if results_dir else None,
                                'status': 'cached_no_stats'
                            }
                        all_results.append(result_summary)
                        processed += 1
                        continue
                else:
                    print(f"  [{processed+1}/{total_genes}] {gene_id}: VCF已更新，重新扫描")
            except Exception as e:
                print(f"  [{processed+1}/{total_genes}] {gene_id}: 缓存检测失败 ({e})，重新扫描")
        
        result = process_single_gene(gene_info, vcf_file, pheno_df, 
                                     database_dir, results_dir, min_samples,
                                     test_region_length=test_region_length,
                                     gff_file=gff_file)
        all_results.append(result)
        processed += 1
        
        # 关键修复：只要处理过（无论成功与否），都保存 VCF 修改时间，避免重复扫描
        if result.get('data_folder'):
            gene_info_file = os.path.join(database_dir, result['data_folder'], 'gene_info.json')
            if os.path.exists(gene_info_file):
                try:
                    with open(gene_info_file, 'r') as f:
                        info = json.load(f)
                    info['vcf_mtime'] = os.path.getmtime(vcf_file)
                    with open(gene_info_file, 'w') as f:
                        json.dump(info, f, indent=2)
                except Exception as e:
                    print(f"  [WARNING] 保存 VCF 修改时间失败: {e}")
        
        print(f"  [{processed}/{total_genes}] {result['gene_id']}: {result['status']}")
        
        # 记录第一个成功基因的数据格式
        if result['status'] == 'success' and processed == 1:
            gene_data_dir = os.path.join(database_dir, result['data_folder'])
            print("\n" + "=" * 60)
            print(f"[DATA FORMAT] 单倍型数据库格式 (基因: {result['gene_id']})")
            print("=" * 60)
            
            # 读取并显示 gene_info.json
            gene_info_file = os.path.join(gene_data_dir, 'gene_info.json')
            if os.path.exists(gene_info_file):
                with open(gene_info_file, 'r') as f:
                    gene_info = json.load(f)
                print(f"\n  gene_info.json:")
                print(f"    gene_id: {gene_info.get('gene_id')}")
                print(f"    chrom: {gene_info.get('chrom')}")
                print(f"    start: {gene_info.get('start'):,}")
                print(f"    end: {gene_info.get('end'):,}")
                print(f"    strand: {gene_info.get('strand')}")
                print(f"    length: {gene_info.get('length'):,}")
                print(f"    promoter_length: {gene_info.get('promoter_length')}")
                print(f"    exons: {len(gene_info.get('exons', []))} 个")
                print(f"    cds: {len(gene_info.get('cds', []))} 个")
            
            # 读取并显示 haplotype_data.csv
            hap_data_file = os.path.join(gene_data_dir, 'haplotype_data.csv')
            if os.path.exists(hap_data_file):
                hap_data = pd.read_csv(hap_data_file)
                print(f"\n  haplotype_data.csv:")
                print(f"    列: {list(hap_data.columns)}")
                print(f"    单倍型数: {len(hap_data)}")
                if len(hap_data) > 0:
                    print(f"    示例 (前3行):")
                    for i, row in hap_data.head(3).iterrows():
                        hap_name = row.get('Hap_Name', 'N/A')
                        count = row.get('Count', 'N/A')
                        print(f"      [{i+1}] {hap_name}: count={count}")
            
            # 读取并显示 haplotype_samples.csv
            hap_samples_file = os.path.join(gene_data_dir, 'haplotype_samples.csv')
            if os.path.exists(hap_samples_file):
                hap_samples = pd.read_csv(hap_samples_file)
                print(f"\n  haplotype_samples.csv:")
                print(f"    列: {list(hap_samples.columns)}")
                print(f"    样本数: {len(hap_samples)}")
                if len(hap_samples) > 0:
                    print(f"    示例 (前5行):")
                    for i, row in hap_samples.head(5).iterrows():
                        sample_id = row.get('SampleID', 'N/A')
                        hap_name = row.get('Hap_Name', 'N/A')
                        print(f"      [{i+1}] SampleID={sample_id}, Hap_Name={hap_name}")
            
            # 读取并显示 variant_info.csv
            variant_info_file = os.path.join(gene_data_dir, 'variant_info.csv')
            if os.path.exists(variant_info_file):
                variant_info = pd.read_csv(variant_info_file)
                print(f"\n  variant_info.csv:")
                print(f"    列: {list(variant_info.columns)}")
                print(f"    变异数: {len(variant_info)}")
                if len(variant_info) > 0:
                    print(f"    示例 (前5行):")
                    for i, row in variant_info.head(5).iterrows():
                        pos = row.get('position', 'N/A')
                        ref = row.get('ref', 'N/A')
                        alt = row.get('alt', 'N/A')
                        maf = row.get('maf', 'N/A')
                        is_sv = row.get('is_sv', 'N/A')
                        print(f"      [{i+1}] pos={pos}, ref={ref}, alt={alt}, maf={maf}, is_sv={is_sv}")
            
            # 读取并显示 phenotype_data.csv
            pheno_data_file = os.path.join(gene_data_dir, 'phenotype_data.csv')
            if os.path.exists(pheno_data_file):
                pheno_data = pd.read_csv(pheno_data_file)
                print(f"\n  phenotype_data.csv:")
                print(f"    列: {list(pheno_data.columns)}")
                print(f"    样本数: {len(pheno_data)}")
                if len(pheno_data) > 0:
                    print(f"    示例 (前5行):")
                    for i, row in pheno_data.head(5).iterrows():
                        sample_id = row.get('SampleID', 'N/A')
                        hap_name = row.get('Hap_Name', 'N/A')
                        pheno_val = row.iloc[-1] if len(row) > 2 else 'N/A'
                        print(f"      [{i+1}] SampleID={sample_id}, Hap_Name={hap_name}, 表型值={pheno_val}")
            
            # 读取并显示 haplotype_stats.csv
            hap_stats_file = os.path.join(gene_data_dir, 'haplotype_stats.csv')
            if os.path.exists(hap_stats_file):
                hap_stats = pd.read_csv(hap_stats_file)
                print(f"\n  haplotype_stats.csv:")
                print(f"    列: {list(hap_stats.columns)}")
                print(f"    单倍型数: {len(hap_stats)}")
                if len(hap_stats) > 0:
                    print(f"    完整内容:")
                    for i, row in hap_stats.iterrows():
                        hap_name = row.get('haplotype_name', 'N/A')
                        count = row.get('haplotype_count', 'N/A')
                        freq = row.get('haplotype_freq', 'N/A')
                        mean = row.get('phenotype_mean', 'N/A')
                        sd = row.get('phenotype_sd', 'N/A')
                        print(f"      [{i+1}] {hap_name}: count={count}, freq={freq}, mean={mean}, sd={sd}")
            
            print("=" * 60 + "\n")
        
        if result['status'].startswith('error'):
            errors += 1
    
    perf.step_end("scan_genes")
    
    # 4. 生成结果表
    results_df = pd.DataFrame(all_results)
    
    # 5. 为每个基因进行关联分析并保存到对应文件夹
    if run_analysis:
        print("\n[INFO] 开始关联分析...")
        perf.step_start("association_analysis")
        
        # 对每个成功处理的基因进行关联分析
        for result in all_results:
            if result['status'] != 'success':
                continue
            
            gene_id = result['gene_id']
            gene_data_dir = os.path.join(database_dir, result['data_folder'])
            
            try:
                # 读取该基因的数据
                hap_sample_df = pd.read_csv(os.path.join(gene_data_dir, 'haplotype_samples.csv'))
                pheno_data = pd.read_csv(os.path.join(gene_data_dir, 'phenotype_data.csv'))
                
                # 获取表型列
                pheno_cols = [c for c in pheno_df.columns if c != 'SampleID' and pheno_df[c].dtype in ['float64', 'int64']]
                if not pheno_cols:
                    pheno_cols = [c for c in pheno_df.columns if c != 'SampleID']
                pheno_col = pheno_cols[0] if pheno_cols else None
                
                if pheno_col is None:
                    continue
                
                # 进行ANOVA分析
                from scipy import stats
                hap_col = 'Hap_Name' if 'Hap_Name' in pheno_data.columns else 'Haplotype'
                
                # 过滤掉Other单倍型
                analysis_df = pheno_data[pheno_data[hap_col] != 'Other'].copy()
                
                if len(analysis_df) < min_samples:
                    continue
                
                hap_means = analysis_df.groupby(hap_col)[pheno_col].mean()
                groups = [analysis_df[analysis_df[hap_col] == h][pheno_col].dropna().values 
                          for h in hap_means.index if len(analysis_df[analysis_df[hap_col] == h]) >= 2]
                
                if len(groups) < 2:
                    continue
                
                f_stat, p_value = stats.f_oneway(*groups)
                
                # 计算PVE
                grand_mean = analysis_df[pheno_col].mean()
                ss_between = sum(len(analysis_df[analysis_df[hap_col] == h]) * 
                                (analysis_df[analysis_df[hap_col] == h][pheno_col].mean() - grand_mean)**2 
                                for h in hap_means.index)
                ss_total = ((analysis_df[pheno_col] - grand_mean)**2).sum()
                pve = (ss_between / ss_total) * 100 if ss_total > 0 else 0
                
                # 保存关联分析结果
                assoc_result = {
                    'gene_id': gene_id,
                    'chrom': result['chrom'],
                    'start': result['start'],
                    'end': result['end'],
                    'f_statistic': f_stat,
                    'p_value': p_value,
                    'significant': p_value < pvalue_threshold,
                    'pve_percent': pve,
                    'n_haplotypes': result['n_haplotypes'],
                    'n_samples': result['n_samples']
                }
                
                assoc_df = pd.DataFrame([assoc_result])
                assoc_df.to_csv(os.path.join(gene_data_dir, 'association_result.csv'), index=False)
                print(f"  [关联分析] {gene_id}: P={p_value:.4f}, PVE={pve:.2f}%")
                
                # 如果显著，生成HTML图
                if p_value < pvalue_threshold and results_dir and HAPLOTYPE_MODULE_AVAILABLE:
                    try:
                        gene_results_dir = os.path.join(results_dir, gene_id)
                        os.makedirs(gene_results_dir, exist_ok=True)
                        
                        # 读取单倍型数据
                        hap_df = pd.read_csv(os.path.join(gene_data_dir, 'haplotype_data.csv'))
                        
                        # 创建分析器并生成HTML
                        analyzer = HaplotypePhenotypeAnalyzer(
                            vcf_file=vcf_file,
                            phenotype_file=pheno_file,
                            output_dir=gene_results_dir,
                            gtf_file=gff_file
                        )
                        
                        # 运行分析生成HTML（传入database_dir以使用预计算数据）
                        analyzer.analyze_gene(
                            chrom=result['chrom'],
                            start=result['start'],
                            end=result['end'],
                            gene_id=gene_id,
                            phenotype_cols=[pheno_col],
                            cluster_haplotypes=cluster_haplotypes,
                            database_dir=database_dir
                        )
                        print(f"  [HTML] {gene_id}: 综合分析图已生成")
                        
                        # 保存分析中间数据到数据库，供后续直接使用
                        # 1. 保存 variant_info
                        if hasattr(analyzer, 'extractor') and hasattr(analyzer.extractor, 'variant_info'):
                            variant_info = analyzer.extractor.variant_info
                            if variant_info:
                                variant_info_df = pd.DataFrame([
                                    {
                                        'position': pos,
                                        'ref': info.get('ref', ''),
                                        'alt': info.get('alt', ''),
                                        'len_diff': info.get('len_diff', 0),
                                        'is_sv': info.get('is_sv', False),
                                        'maf': info.get('maf', 0.5),
                                        'missing_rate': info.get('missing_rate', 0.0)
                                    }
                                    for pos, info in variant_info.items()
                                ])
                                variant_info_df.to_csv(os.path.join(gene_data_dir, 'variant_info.csv'), index=False)
                                print(f"  [DB] {gene_id}: variant_info 已保存")
                        
                        # 2. 保存 positions
                        if hasattr(analyzer, 'positions') and analyzer.positions:
                            positions_df = pd.DataFrame({'position': analyzer.positions})
                            positions_df.to_csv(os.path.join(gene_data_dir, 'variant_positions_detailed.csv'), index=False)
                    except Exception as html_e:
                        print(f"  [WARNING] {gene_id} HTML生成失败: {html_e}")
                
            except Exception as e:
                print(f"  [WARNING] {gene_id} 关联分析失败: {e}")
        
        perf.step_end("association_analysis")
    
    # 6. 生成数据集汇总
    perf.step_start("save_results")
    
    # 保存汇总表
    summary_file = os.path.join(database_dir, "summary.csv")
    results_df.to_csv(summary_file, index=False)
    print(f"\n[INFO] 数据集汇总已保存: {summary_file}")
    
    # 打印数据集结构
    print("\n" + "=" * 60)
    print("单倍型数据集结构")
    print("=" * 60)
    print(f"\n数据库目录: {database_dir}/")
    for result in all_results:
        gene_data_dir = os.path.join(database_dir, result['data_folder'])
        print(f"  {result['data_folder']}/")
        if os.path.exists(gene_data_dir):
            for f in sorted(os.listdir(gene_data_dir)):
                print(f"    ├── {f}")
    
    if results_dir:
        print(f"\n结果目录: {results_dir}/")
        for result in all_results:
            if result['results_folder']:
                gene_results_dir = os.path.join(results_dir, result['results_folder'])
                if os.path.exists(gene_results_dir):
                    print(f"  {result['results_folder']}/")
                    for f in sorted(os.listdir(gene_results_dir)):
                        print(f"    ├── {f}")
    
    perf.step_end("save_results")
    
    # 7. 打印统计
    print("\n" + "=" * 60)
    print("数据集构建完成统计")
    print("=" * 60)
    print(f"总基因数: {total_genes}")
    print(f"成功处理: {processed - errors}")
    print(f"处理失败: {errors}")
    print(f"数据库目录: {database_dir}")
    if results_dir:
        print(f"结果目录: {results_dir}")
    
    # 状态统计
    if 'status' in results_df.columns:
        print(f"\n状态分布:")
        for status, count in results_df['status'].value_counts().items():
            print(f"  {status}: {count}")
    
    # 性能报告
    perf.report_performance()
    
    return results_df


# ============================================================================
# 命令行入口
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="构建单倍型数据集 - 针对指定基因",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--vcf", default=ScanConfig.VCF_FILE,
                        help="VCF 文件路径")
    parser.add_argument("--gff", default=ScanConfig.GFF_FILE,
                        help="GFF3 注释文件路径")
    parser.add_argument("--phenotype", default=ScanConfig.PHENO_FILE,
                        help="表型文件路径")
    parser.add_argument("--database-dir", default=ScanConfig.DATABASE_DIR,
                        help="数据库输出目录")
    parser.add_argument("--results-dir", default=ScanConfig.RESULTS_DIR,
                        help="结果(图表)输出目录")
    parser.add_argument("--genes", nargs="+", default=None,
                        help=f"指定基因列表（默认: {ScanConfig.TARGET_GENES})")
    parser.add_argument("--min-samples", type=int, default=5,
                        help="最小样本数阈值")
    parser.add_argument("--no-analysis", action="store_true",
                        help="跳过关联分析")
    parser.add_argument("--pvalue-threshold", type=float, default=0.05,
                        help="P值显著性阈值")
    parser.add_argument("--test-region", type=int, default=ScanConfig.TEST_REGION_LENGTH,
                        help="测试模式：限制区间长度(bp)，0表示使用完整基因区间")
    parser.add_argument("--cluster", action="store_true",
                        help="按序列相似度聚类排序单倍型（默认按数量排序）")
    
    args = parser.parse_args()
    
    # 运行扫描
    run_genome_scan(
        vcf_file=args.vcf,
        gff_file=args.gff,
        pheno_file=args.phenotype,
        database_dir=args.database_dir,
        results_dir=args.results_dir,
        gene_filter=args.genes,
        min_samples=args.min_samples,
        run_analysis=not args.no_analysis,
        pvalue_threshold=args.pvalue_threshold,
        test_region_length=args.test_region,
        cluster_haplotypes=args.cluster
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
        analyzed = len(assoc_df[assoc_df['status'] == 'success']) if 'status' in assoc_df.columns else 0
        sig_count = len(assoc_df[assoc_df['significant'] == True]) if 'significant' in assoc_df.columns else 0
        total_tests = len(assoc_df[assoc_df['p_value'].notna()]) if 'p_value' in assoc_df.columns else 0
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
    if assoc_df is not None and 'significant' in assoc_df.columns and len(assoc_df[assoc_df['significant'] == True]) > 0:
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
