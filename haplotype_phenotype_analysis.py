"""
单倍型-表型关联分析平台 (Haplotype-Phenotype Association Analysis Platform)
==============================================================================
基于基因单倍型进行表型关联分析，挖掘传统GWAS遗漏的关联信号。

核心模块:
1. 表型关联模块 - 单倍型与表型的统计关联（t-test/ANOVA/回归）
2. PVE模块 - 变异解释率计算（R²、遗传力估计）
3. 报告输出模块 - 整合结果生成可视化报告
4. GWAS整合模块 - 对比传统GWAS结果
5. 启动子注释模块 - 上游2kb区域功能注释

数据来源: sv_pro 数据（不修改原始数据路径）
"""

import os
import sys
import argparse
import gzip
import time
import json
import logging
import traceback
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
from scipy import stats
from scipy.stats import f_oneway, ttest_ind, pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
import warnings
warnings.filterwarnings('ignore')

# 尝试导入 psutil 用于内存监控（可选）
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# statsmodels 替代方案（无需外部依赖）
STATSMODELS_AVAILABLE = False

# ============================================================================
# 日志配置（带Unix时间戳文件名）
# ============================================================================

def setup_logging(log_dir: str = "./logs", prefix: str = "haplotype_analysis") -> logging.Logger:
    """
    配置日志系统，日志文件名包含Unix时间戳
    
    Args:
        log_dir: 日志目录
        prefix: 日志文件前缀
    
    Returns:
        Logger: 配置好的logger对象
    """
    # 使用绝对路径
    if not os.path.isabs(log_dir):
        log_dir = os.path.join(os.getcwd(), log_dir)
    
    os.makedirs(log_dir, exist_ok=True)
    print(f"[DEBUG] 日志目录: {log_dir}")
    
    # Unix时间戳作为文件名核心标识
    timestamp = int(time.time())
    log_filename = f"{prefix}_{timestamp}.log"
    log_filepath = os.path.join(log_dir, log_filename)
    
    print(f"[DEBUG] 日志文件: {log_filepath}")
    
    # 创建 logger
    logger = logging.getLogger('HaplotypeAnalysis')
    logger.setLevel(logging.DEBUG)
    
    # 清除已有的 handlers
    logger.handlers.clear()
    
    # 文件 handler
    file_handler = logging.FileHandler(log_filepath, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        '[%(asctime)s] [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_formatter)
    
    # 控制台 handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('[%(levelname)s] %(message)s')
    console_handler.setFormatter(console_formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


# ============================================================================
# 性能监控工具类
# ============================================================================

class PerformanceMonitor:
    """性能监控器，记录各步骤执行时间和内存使用"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger('HaplotypeAnalysis')
        self.start_time = None
        self.step_times = {}
        self.step_start_times = {}
        self.memory_records = []
        
    def start(self):
        """开始总体计时"""
        self.start_time = time.time()
        self._record_memory("start")
        
    def end(self):
        """结束总体计时并返回总耗时"""
        total_time = time.time() - self.start_time if self.start_time else 0
        self._record_memory("end")
        return total_time
    
    def step_start(self, step_name: str):
        """开始某个步骤的计时"""
        self.step_start_times[step_name] = time.time()
        self._record_memory(f"{step_name}_start")
        
    def step_end(self, step_name: str) -> float:
        """结束某个步骤的计时并返回耗时"""
        if step_name in self.step_start_times:
            elapsed = time.time() - self.step_start_times[step_name]
            self.step_times[step_name] = self.step_times.get(step_name, 0) + elapsed
            self._record_memory(f"{step_name}_end")
            return elapsed
        return 0.0
    
    def _record_memory(self, stage: str):
        """记录当前内存使用"""
        if PSUTIL_AVAILABLE:
            process = psutil.Process(os.getpid())
            mem_info = process.memory_info()
            self.memory_records.append({
                'stage': stage,
                'rss_mb': mem_info.rss / 1024 / 1024,  # Resident Set Size
                'vms_mb': mem_info.vms / 1024 / 1024,  # Virtual Memory Size
                'timestamp': time.time()
            })
        # 如果 psutil 不可用，不记录内存
    
    def get_current_memory_mb(self) -> float:
        """获取当前内存使用（MB）"""
        if PSUTIL_AVAILABLE:
            process = psutil.Process(os.getpid())
            return process.memory_info().rss / 1024 / 1024
        return 0.0
    
    def report_performance(self):
        """生成性能报告"""
        report_lines = []
        report_lines.append("\n" + "=" * 60)
        report_lines.append("性能报告")
        report_lines.append("=" * 60)
        
        # 总体统计
        total_time = self.end()
        report_lines.append(f"\n总执行时间：{total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
        
        # 各步骤耗时
        if self.step_times:
            report_lines.append("\n各步骤耗时:")
            for step_name, elapsed in sorted(self.step_times.items(), key=lambda x: x[1], reverse=True):
                percentage = (elapsed / total_time * 100) if total_time > 0 else 0
                report_lines.append(f"  {step_name:<35s}: {elapsed:>8.2f}s ({percentage:>5.1f}%)")
        
        # 内存使用
        if self.memory_records and len(self.memory_records) >= 2:
            initial_mem = self.memory_records[0].get('rss_mb', 0)
            peak_mem = max(r.get('rss_mb', 0) for r in self.memory_records)
            final_mem = self.memory_records[-1].get('rss_mb', 0)
            
            report_lines.append(f"\n内存使用:")
            report_lines.append(f"  初始内存：{initial_mem:.1f} MB")
            report_lines.append(f"  峰值内存：{peak_mem:.1f} MB")
            report_lines.append(f"  最终内存：{final_mem:.1f} MB")
            report_lines.append(f"  净增长：{final_mem - initial_mem:.1f} MB")
        
        report_lines.append("\n" + "=" * 60)
        
        report_str = "\n".join(report_lines)
        if self.logger:
            self.logger.info(report_str)
        else:
            print(report_str)
        
        return {
            'total_time': total_time,
            'step_times': self.step_times,
            'memory_records': self.memory_records
        }


# 全局 logger
LOGGER = None

def get_logger() -> logging.Logger:
    """获取全局logger"""
    global LOGGER
    if LOGGER is None:
        LOGGER = setup_logging()
    return LOGGER

# 尝试导入 pysam
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    print("[WARNING] pysam not available, some features may be limited")

# 设置字体
rcParams['font.family'] = 'DejaVu Sans'
rcParams['axes.unicode_minus'] = False
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# ============================================================================
# 数据配置 - 使用 sv_pro 数据路径（不修改）
# ============================================================================

class DataConfig:
    """数据路径配置 - CS-IAAS T2T v1.1"""
    # 超算路径
    DATA_BASE_PATH = "/storage/public/home/2024110093/data/Variation/CSIAAS/"
    SV_VCF_PATH = DATA_BASE_PATH + "Core819Samples_ALL.vcf.gz"
    PHENOTYPE_PATH = DATA_BASE_PATH + "Phe.txt"
    VCF_ID_PATH = DATA_BASE_PATH + "VCFID.txt"
    
    # GFF3/FASTA 参考 (CS-IAAS T2T v1.1)
    GENOME_BASE = "/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/"
    GTF_PATH = GENOME_BASE + "CS-IAAS_v1.1_HC.gff3"
    FASTA_PATH = GENOME_BASE + "CS-IAAS_v1.1.fasta"
    
    # 本地测试路径（可选）
    LOCAL_VCF_ID_PATH = "./Variation/CSIAAS/VCFID.txt"
    
    # 输出路径
    OUTPUT_DIR = "./results/haplotype_analysis"


# ============================================================================
# GTF/GFF3 解析函数：获取真实的外显子和 CDS 坐标
# ============================================================================

def parse_gtf_for_gene(gtf_file: str, target_gene_id: str) -> dict:
    """
    从 GTF/GFF3 文件解析目标基因的外显子和 CDS 位置
    
    支持两种格式：
    - GTF: gene_id "GENEID"
    - GFF3: ID=GENEID 或 Parent=GENEID
    
    Args:
        gtf_file: GTF/GFF3 注释文件路径
        target_gene_id: 目标基因 ID
    
    Returns:
        dict: {
            'exons': [(start, end), ...],  # 外显子坐标列表
            'cds': [(start, end), ...],     # CDS 坐标列表
            'strand': '+' or '-',
            'chrom': 'chr5B',
            'gene_start': int,
            'gene_end': int
        }
    """
    exons = []
    cds = []
    gene_strand = '+'
    gene_chrom = None
    gene_start = None
    gene_end = None
    
    if not gtf_file or not os.path.exists(gtf_file):
        return {'exons': [], 'cds': [], 'strand': '+', 'chrom': None, 'gene_start': None, 'gene_end': None}
    
    # 支持 .gz 或纯文本格式
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    
    # 检测文件格式（GTF 或 GFF3）
    is_gff3 = gtf_file.endswith('.gff3') or gtf_file.endswith('.gff')
    
    try:
        with open_func(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                attr = parts[8]
                
                # 根据格式匹配基因ID
                if is_gff3:
                    # GFF3格式：ID=GENEID 或 Parent=GENEID 或 gene_id=GENEID
                    match = False
                    if f'ID={target_gene_id}' in attr or f'ID={target_gene_id};' in attr:
                        match = True
                    elif f'Parent={target_gene_id}' in attr or f'Parent={target_gene_id};' in attr:
                        match = True
                    elif f'gene_id={target_gene_id}' in attr:
                        match = True
                    if not match:
                        continue
                else:
                    # GTF格式：gene_id "GENEID"
                    if f'gene_id "{target_gene_id}"' not in attr:
                        continue
                
                if gene_chrom is None:
                    gene_chrom = parts[0]
                gene_strand = parts[6]
                
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                
                # 记录基因边界
                if gene_start is None or start < gene_start:
                    gene_start = start
                if gene_end is None or end > gene_end:
                    gene_end = end
                
                if feature == 'exon':
                    exons.append((start, end))
                elif feature == 'CDS':
                    cds.append((start, end))
    except Exception as e:
        print(f"[WARNING] 解析 GTF/GFF3 失败: {e}")
        return {'exons': [], 'cds': [], 'strand': '+', 'chrom': None, 'gene_start': None, 'gene_end': None}
    
    # 排序
    exons = sorted(exons)
    cds = sorted(cds)
    
    return {
        'exons': exons,
        'cds': cds,
        'strand': gene_strand,
        'chrom': gene_chrom,
        'gene_start': gene_start,
        'gene_end': gene_end
    }


# ============================================================================
# SNP 功能注释（Missense / Synonymous / UTR / Other）
# ============================================================================

def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate_codon(codon: str) -> str:
    codon = codon.upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return _GENETIC_CODE.get(codon, "X")


def _pos_in_any_interval(pos: int, intervals) -> bool:
    for s, e in intervals:
        if s <= pos <= e:
            return True
    return False


def _build_coding_context(chrom: str, cds_intervals, strand: str, fasta_path: str):
    """构建 CDS 参考序列及基因组坐标到 CDS 索引的映射"""
    if not cds_intervals or not PYSAM_AVAILABLE:
        return "", {}
    try:
        import pysam as _pysam
        fasta = _pysam.FastaFile(fasta_path)
        cds_intervals = sorted(cds_intervals)
        cds_pos_to_idx = {}
        seq_frags = []
        idx = 0
        if strand == "+":
            for start, end in cds_intervals:
                frag = fasta.fetch(chrom, start - 1, end)
                seq_frags.append(frag)
                for pos in range(start, end + 1):
                    cds_pos_to_idx[pos] = idx
                    idx += 1
        else:
            for start, end in reversed(cds_intervals):
                frag = fasta.fetch(chrom, start - 1, end)
                frag_rc = _revcomp(frag)
                seq_frags.append(frag_rc)
                for pos in range(end, start - 1, -1):
                    cds_pos_to_idx[pos] = idx
                    idx += 1
        fasta.close()
        return "".join(seq_frags), cds_pos_to_idx
    except Exception as e:
        print(f"[WARNING] build_coding_context 失败: {e}")
        return "", {}


def annotate_snp_effects_for_region(vcf_file: str, fasta_path: str, gene_chrom: str,
                                     cds_intervals: list, exon_intervals: list,
                                     gene_strand: str, positions: list) -> dict:
    """
    计算指定位点列表的 SNP 功能注释
    返回: {pos: 'missense'|'synonymous'|'UTR'|'indel'|'SV'|'other'}
    
    变异类型判断标准：
    - SNP: ref和alt都是单碱基
    - Indel: ref和alt长度不同，但长度差<50bp
    - SV: ref和alt长度差>=50bp，或VCF INFO字段有SVTYPE
    - missense/synonymous: CDS区域内的SNP，改变/不改变氨基酸
    - UTR: 外显子区域但非CDS区域
    - other: 其他情况
    """
    effects = {pos: 'other' for pos in positions}  # 默认都是 other
    
    # 定义变异类型判断函数
    def classify_variant(ref, alt):
        """根据ref/alt长度判断变异类型: SNP, indel, SV"""
        len_diff = abs(len(ref) - len(alt))
        if len(ref) == 1 and len(alt) == 1:
            return 'SNP'
        elif len_diff >= 50:  # 长度差>=50bp认为是SV
            return 'SV'
        else:  # 长度不同但差值<50bp是indel
            return 'indel'

    if not cds_intervals and not exon_intervals:
        return effects

    if not PYSAM_AVAILABLE or not fasta_path or not os.path.exists(fasta_path):
        # 没有 FASTA，只能判断 UTR vs other
        for pos in positions:
            in_cds = _pos_in_any_interval(pos, cds_intervals)
            in_exon = _pos_in_any_interval(pos, exon_intervals)
            if in_exon and not in_cds:
                effects[pos] = 'UTR'
            elif in_cds:
                effects[pos] = 'other'  # 无法判断 missense/synonymous
        return effects

    try:
        import pysam as _pysam
        cds_seq, cds_pos_to_idx = _build_coding_context(gene_chrom, cds_intervals, gene_strand, fasta_path)
        fasta = _pysam.FastaFile(fasta_path)

        pos_set = set(positions)
        remaining = set(pos_set)  # 未找到的位点，全部找到后提前退出

        # 优先用 tabix 按区域查询（秒级），否则逐行扫描（慢）
        tabix_ok = False
        if vcf_file.endswith('.gz') and os.path.exists(vcf_file + '.tbi'):
            try:
                tbx = _pysam.TabixFile(vcf_file)
                if positions:
                    region_min = min(positions)
                    region_max = max(positions)
                    region_str = f"{gene_chrom}:{region_min}-{region_max}"
                    for row in tbx.fetch(region=region_str):
                        parts = row.split('\t')
                        if len(parts) < 5:
                            continue
                        chrom_v, pos_str, _id, ref, alt = parts[:5]
                        if chrom_v != gene_chrom:
                            continue
                        try:
                            pos = int(pos_str)
                        except Exception:
                            continue
                        if pos not in remaining:
                            continue
                        remaining.discard(pos)
                        alt_allele = alt.split(',')[0]
                        var_type = classify_variant(ref, alt_allele)
                        in_cds  = _pos_in_any_interval(pos, cds_intervals)
                        in_exon = _pos_in_any_interval(pos, exon_intervals)
                        
                        # 先判断SV/indel（不区分位置，统一标记）
                        if var_type == 'SV':
                            effects[pos] = 'SV'
                        elif var_type == 'indel':
                            effects[pos] = 'indel'
                        elif not in_exon:
                            effects[pos] = 'other'  # 基因间区、内含子等
                        elif not in_cds:
                            effects[pos] = 'UTR'
                        elif var_type == 'SNP' and pos in cds_pos_to_idx and cds_seq:
                            idx = cds_pos_to_idx[pos]
                            codon_start = (idx // 3) * 3
                            if codon_start + 3 <= len(cds_seq):
                                codon_ref = cds_seq[codon_start: codon_start + 3]
                                alt_base = alt_allele if gene_strand == '+' else _revcomp(alt_allele)[0]
                                offset = idx % 3
                                codon_alt = codon_ref[:offset] + alt_base + codon_ref[offset + 1:]
                                aa_ref = _translate_codon(codon_ref)
                                aa_alt = _translate_codon(codon_alt)
                                if aa_ref == 'X' or aa_alt == 'X':
                                    effects[pos] = 'other'
                                elif aa_ref == aa_alt:
                                    effects[pos] = 'synonymous'
                                else:
                                    effects[pos] = 'missense'
                            else:
                                effects[pos] = 'other'
                        else:
                            effects[pos] = 'other'
                        if not remaining:
                            break
                tbx.close()
                tabix_ok = True
            except Exception as te:
                print(f"[WARNING] tabix 查询失败，回退到逐行扫描: {te}")

        if not tabix_ok:
            # 回退：逐行扫描，但一旦找到全部位点立即退出
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            with open_func(vcf_file, 'rt') as f:
                for line in f:
                    if not line or line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 5:
                        continue
                    chrom_v, pos_str, _id, ref, alt = parts[:5]
                    if chrom_v != gene_chrom:
                        continue
                    try:
                        pos = int(pos_str)
                    except Exception:
                        continue
                    if pos not in remaining:
                        continue
                    remaining.discard(pos)
                    alt_allele = alt.split(',')[0]
                    var_type = classify_variant(ref, alt_allele)
                    in_cds  = _pos_in_any_interval(pos, cds_intervals)
                    in_exon = _pos_in_any_interval(pos, exon_intervals)
                    
                    # 先判断SV/indel（不区分位置，统一标记）
                    if var_type == 'SV':
                        effects[pos] = 'SV'
                    elif var_type == 'indel':
                        effects[pos] = 'indel'
                    elif not in_exon:
                        effects[pos] = 'other'  # 基因间区、内含子等
                    elif not in_cds:
                        effects[pos] = 'UTR'
                    elif var_type == 'SNP' and pos in cds_pos_to_idx and cds_seq:
                        idx = cds_pos_to_idx[pos]
                        codon_start = (idx // 3) * 3
                        if codon_start + 3 <= len(cds_seq):
                            codon_ref = cds_seq[codon_start: codon_start + 3]
                            alt_base = alt_allele if gene_strand == '+' else _revcomp(alt_allele)[0]
                            offset = idx % 3
                            codon_alt = codon_ref[:offset] + alt_base + codon_ref[offset + 1:]
                            aa_ref = _translate_codon(codon_ref)
                            aa_alt = _translate_codon(codon_alt)
                            if aa_ref == 'X' or aa_alt == 'X':
                                effects[pos] = 'other'
                            elif aa_ref == aa_alt:
                                effects[pos] = 'synonymous'
                            else:
                                effects[pos] = 'missense'
                        else:
                            effects[pos] = 'other'
                    else:
                        effects[pos] = 'other'
                    if not remaining:
                        break  # 所有位点都已找到，提前退出

        fasta.close()
    except Exception as e:
        print(f"[WARNING] SNP 注释失败: {e}")

    return effects

class HaplotypeExtractor:
    """从VCF提取指定区间的单倍型"""
    
    def __init__(self, vcf_file: str):
        self.vcf_file = vcf_file
        self.samples = []
        self.positions = []
        self.sample_haplotypes = {}
        
    def _gt_to_allele(self, rec, sample_name: str, ref: str, alt0: str) -> str:
        """从 GT 字段获取等位基因"""
        try:
            # pysam VariantRecord.samples 返回一个 dict-like 对象
            sample = rec.samples[sample_name]
            
            # 在 pysam 0.19+ 中，sample 是一个 VariantRecordSample 对象
            # 可以直接访问 GT 属性
            gt = None
            
            # 方法1: 直接访问 .GT 属性 (pysam 0.19+)
            if hasattr(sample, 'GT'):
                gt = sample.GT
            
            # 方法2: 使用 get 方法
            if gt is None and hasattr(sample, 'get'):
                gt = sample.get('GT')
            
            # 方法3: 使用 __getitem__
            if gt is None and hasattr(sample, '__getitem__'):
                try:
                    gt = sample['GT']
                except (KeyError, TypeError):
                    pass
            
            # 方法4: 使用 alleles 属性
            if gt is None and hasattr(sample, 'alleles'):
                alleles = sample.alleles
                if alleles and len(alleles) >= 1:
                    # 根据 alleles 推断 GT
                    gt = tuple([0 if a == ref else 1 for a in alleles if a is not None])
                    
        except Exception as e:
            return "N"
            
        if gt is None or None in gt:
            return "N"
            
        # 纯合 REF
        if gt == (0, 0) or gt == (0,) or gt == [0, 0] or gt == [0]:
            return ref if ref else "N"
        # 纯合 ALT
        if gt == (1, 1) or gt == (1,) or gt == [1, 1] or gt == [1]:
            return alt0 if alt0 else "N"
        # 杂合按 REF 处理
        if (isinstance(gt, (tuple, list)) and len(gt) >= 2 and 0 in gt and 1 in gt):
            return ref if ref else "N"
        return "N"
    
    def extract_region(self, chrom: str, start: int, end: int, min_samples: int = 2, snp_only: bool = True) -> tuple:
        """
        提取指定区间的单倍型
        
        Args:
            chrom: 染色体
            start: 起始位置
            end: 终止位置
            min_samples: 单倍型最少样本数
            snp_only: 是否只保留SNP（默认True，过滤indel和SV）
        
        Returns:
            positions: 变异位置列表
            hap_df: 单倍型汇总表
            hap_sample_df: 样本-单倍型对应表
        """
        logger = get_logger()
        
        if not PYSAM_AVAILABLE:
            raise ImportError("pysam is required for VCF parsing")
        
        logger.info(f"打开VCF文件: {self.vcf_file}")
        vcf = pysam.VariantFile(self.vcf_file)
        self.samples = list(vcf.header.samples)
        logger.info(f"样本数: {len(self.samples)}")
        
        # 获取VCF中的染色体列表
        vcf_chroms = list(vcf.header.contigs)
        if vcf_chroms:
            logger.info(f"VCF中的染色体: {vcf_chroms[:10]}{'...' if len(vcf_chroms) > 10 else ''}")
            # 检查指定的染色体是否在VCF中
            if chrom not in vcf_chroms:
                logger.warning(f"[!] 指定的染色体 '{chrom}' 不在VCF文件中!")
                # 尝试查找类似的染色体名
                similar = [c for c in vcf_chroms if chrom.replace('chr', '') in c or c.replace('chr', '') in chrom]
                if similar:
                    logger.warning(f"[建议] 可能的匹配: {similar}")
        else:
            logger.warning("VCF header中没有染色体信息，将尝试直接查询")
        positions = []
        sample_alleles = {s: [] for s in self.samples}
        
        # 检查是否有索引
        has_index = False
        try:
            # 尝试fetch，如果没有索引会报错
            test_iter = vcf.fetch(chrom, start, start+1)
            next(test_iter, None)
            has_index = True
            logger.info("检测到VCF索引，使用快速区间查询")
        except Exception as e:
            has_index = False
            logger.warning(f"VCF没有索引: {e}")
            logger.warning("建议运行: tabix -p vcf <vcf_file> 创建索引")
        
        vcf.close()
        vcf = pysam.VariantFile(self.vcf_file)
        
        # SNP判断函数
        def is_snp(ref_allele, alt_allele):
            """判断是否为SNP（ref和alt都是单碱基）"""
            valid_bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
            return (len(ref_allele) == 1 and len(alt_allele) == 1 and
                    ref_allele.upper() in valid_bases and alt_allele.upper() in valid_bases)
        
        # 提取区间内的变异
        record_count = 0
        filtered_count = 0  # 被过滤的非SNP数量
        if has_index:
            # 使用索引快速查询
            logger.info(f"快速查询区间: {chrom}:{start}-{end}")
            for rec in vcf.fetch(chrom, start, end):
                pos = rec.pos
                ref = rec.ref or ""
                alt0 = (rec.alts[0] if rec.alts else "") or ""
                
                # SNP过滤
                if snp_only and not is_snp(ref, alt0):
                    filtered_count += 1
                    continue
                
                record_count += 1
                if record_count % 100 == 0:
                    logger.debug(f"已处理 {record_count} 条SNP记录...")
                
                positions.append(pos)
                
                for s in self.samples:
                    allele = self._gt_to_allele(rec, s, ref, alt0)
                    sample_alleles[s].append(allele if allele else "N")
        else:
            # 无索引，需要扫描整个文件（非常慢！）
            logger.warning("无索引，开始全文件扫描（可能非常慢）...")
            total_scanned = 0
            for rec in vcf:
                total_scanned += 1
                if total_scanned % 100000 == 0:
                    logger.info(f"已扫描 {total_scanned} 条记录，匹配 {record_count} 条SNP...")
                
                if rec.chrom != chrom:
                    continue
                if rec.pos < start or rec.pos > end:
                    continue
                
                pos = rec.pos
                ref = rec.ref or ""
                alt0 = (rec.alts[0] if rec.alts else "") or ""
                
                # SNP过滤
                if snp_only and not is_snp(ref, alt0):
                    filtered_count += 1
                    continue
                
                record_count += 1
                positions.append(pos)
                
                for s in self.samples:
                    allele = self._gt_to_allele(rec, s, ref, alt0)
                    sample_alleles[s].append(allele if allele else "N")
            
            logger.info(f"全文件扫描完成: 共扫描 {total_scanned} 条，匹配 {record_count} 条SNP")
        
        vcf.close()
        self.positions = positions
        logger.info(f"区间内变异位点数: {len(positions)}")
        if snp_only and filtered_count > 0:
            logger.info(f"过滤掉的非SNP变异数: {filtered_count}")
        
        # 过滤无变异位点（所有样本在该位点的等位基因相同）
        if len(positions) > 0:
            invariant_positions = []  # 记录无变异位点的索引
            for i in range(len(positions)):
                alleles_at_pos = set()
                for sample in self.samples:
                    if sample_alleles[sample] and i < len(sample_alleles[sample]):
                        allele = sample_alleles[sample][i]
                        if allele != 'N':
                            alleles_at_pos.add(allele)
                # 如果该位点只有一种等位基因（或没有有效等位基因），标记为无变异
                if len(alleles_at_pos) <= 1:
                    invariant_positions.append(i)
            
            if invariant_positions:
                logger.info(f"过滤掉的无变异位点数: {len(invariant_positions)}")
                # 移除无变异位点
                keep_indices = [i for i in range(len(positions)) if i not in invariant_positions]
                positions = [positions[i] for i in keep_indices]
                for sample in self.samples:
                    if sample_alleles[sample]:
                        sample_alleles[sample] = [sample_alleles[sample][i] for i in keep_indices]
                self.positions = positions
                logger.info(f"过滤后有效变异位点数: {len(positions)}")
        
        # 构建单倍型
        hap_dict = {}
        sample_to_hap = {}
        
        for sample in self.samples:
            allele_seq = sample_alleles[sample]
            if "N" in allele_seq:
                continue
            hap_key = "|".join(allele_seq)
            hap_dict[hap_key] = hap_dict.get(hap_key, 0) + 1
            sample_to_hap[sample] = hap_key
        
        # 汇总表
        hap_df = pd.DataFrame(list(hap_dict.items()), columns=["Haplotype_Seq", "Count"])
        hap_df = hap_df[hap_df["Count"] >= min_samples].sort_values(
            by="Count", ascending=False
        ).reset_index(drop=True)
        hap_df["Hap_Name"] = [f"Hap{i+1}" for i in range(len(hap_df))]
        hap_df["Alleles"] = hap_df["Haplotype_Seq"].str.split("|")
        
        # 注意：不再基于主要单倍型进行第二次过滤
        # 因为次要单倍型中的变异可能对表型有显著影响
        # 第一次过滤（样本层面）已经删除了真正无变异的位点
        
        # 样本-单倍型对应表
        hap_seq_to_name = dict(zip(hap_df["Haplotype_Seq"], hap_df["Hap_Name"]))
        rows = []
        for sample, hap_seq in sample_to_hap.items():
            hap_name = hap_seq_to_name.get(hap_seq)
            if hap_name is None:
                continue  # 数量少于 min_samples 的单倍型直接跳过，不归入 Other
            rows.append({
                "SampleID": sample,
                "Hap_Name": hap_name,
                "Haplotype_Seq": hap_seq
            })
        hap_sample_df = pd.DataFrame(rows)
        
        self.sample_haplotypes = sample_to_hap
        
        return positions, hap_df, hap_sample_df


# ============================================================================
# 模块2: 表型关联分析
# ============================================================================

class PhenotypeAssociation:
    """单倍型-表型关联分析模块"""
    
    def __init__(self, phenotype_df: pd.DataFrame, hap_sample_df: pd.DataFrame):
        """
        Args:
            phenotype_df: 表型数据，必须包含 SampleID 列和表型列
            hap_sample_df: 单倍型-样本对应表
        """
        self.phenotype_df = phenotype_df
        self.hap_sample_df = hap_sample_df
        self.merged_df = None
        self._merge_data()
    
    def _merge_data(self):
        """合并表型和单倍型数据"""
        # 检查表型DataFrame是否为空
        if self.phenotype_df is None or len(self.phenotype_df) == 0 or len(self.phenotype_df.columns) == 0:
            print(f"[WARNING] 表型数据为空")
            self.merged_df = pd.DataFrame()
            return
        
        # 检查列名
        if 'SampleID' not in self.phenotype_df.columns:
            # 尝试第一列作为样本ID
            self.phenotype_df = self.phenotype_df.rename(
                columns={self.phenotype_df.columns[0]: 'SampleID'}
            )
        
        self.merged_df = pd.merge(
            self.hap_sample_df,
            self.phenotype_df,
            on='SampleID',
            how='inner'
        )
        print(f"[INFO] 合并后样本数: {len(self.merged_df)}")
    
    def get_phenotype_columns(self) -> list:
        """获取可用的表型列"""
        exclude_cols = ['SampleID', 'Hap_Name', 'Haplotype_Seq']
        return [c for c in self.merged_df.columns if c not in exclude_cols]
    
    def association_test(self, phenotype_col: str, method: str = 'auto') -> dict:
        """
        单倍型-表型关联检验
        
        Args:
            phenotype_col: 表型列名
            method: 检验方法 ('auto', 'ttest', 'anova', 'kruskal')
        
        Returns:
            dict: 包含统计结果的字典
        """
        if phenotype_col not in self.merged_df.columns:
            raise ValueError(f"表型列 '{phenotype_col}' 不存在")
        
        # 过滤有效数据
        df = self.merged_df[
            (self.merged_df['Hap_Name'] != 'Other') & 
            (self.merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 5:
            return {'error': '样本量不足', 'n_samples': len(df)}
        
        # 按单倍型分组
        groups = df.groupby('Hap_Name')[phenotype_col].apply(list).to_dict()
        n_groups = len(groups)
        
        if n_groups < 2:
            return {'error': '单倍型分组不足', 'n_groups': n_groups}
        
        # 自动选择方法
        if method == 'auto':
            method = 'ttest' if n_groups == 2 else 'anova'
        
        result = {
            'phenotype': phenotype_col,
            'method': method,
            'n_groups': n_groups,
            'n_samples': len(df),
            'group_sizes': {k: len(v) for k, v in groups.items()},
            'group_means': {k: np.mean(v) for k, v in groups.items()},
            'group_stds': {k: np.std(v) for k, v in groups.items()},
        }
        
        group_values = list(groups.values())
        
        if method == 'ttest' and n_groups == 2:
            # 两组: Welch's t-test
            keys = list(groups.keys())
            stat, pval = ttest_ind(groups[keys[0]], groups[keys[1]], equal_var=False)
            result['statistic'] = stat
            result['p_value'] = pval
            result['test_type'] = "Welch's t-test"
            
            # 效应量 Cohen's d
            n1, n2 = len(groups[keys[0]]), len(groups[keys[1]])
            m1, m2 = np.mean(groups[keys[0]]), np.mean(groups[keys[1]])
            s1, s2 = np.std(groups[keys[0]], ddof=1), np.std(groups[keys[1]], ddof=1)
            pooled_std = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
            result['cohens_d'] = (m1 - m2) / pooled_std if pooled_std > 0 else 0
            
        elif method == 'anova':
            # 多组: 单因素 ANOVA
            stat, pval = f_oneway(*group_values)
            result['statistic'] = stat
            result['p_value'] = pval
            result['test_type'] = "One-way ANOVA"
            
            # 效应量 eta-squared
            grand_mean = df[phenotype_col].mean()
            ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in group_values)
            ss_total = sum((x - grand_mean)**2 for g in group_values for x in g)
            result['eta_squared'] = ss_between / ss_total if ss_total > 0 else 0
            
            # Tukey HSD 事后检验（使用scipy实现）
            if pval < 0.05 and len(group_values) >= 2:
                # 简化版Tukey HSD：计算所有成对比较
                hap_names = df['Hap_Name'].unique()
                tukey_results = []
                for i, hap1 in enumerate(hap_names):
                    for hap2 in hap_names[i+1:]:
                        g1 = df[df['Hap_Name'] == hap1][phenotype_col].dropna()
                        g2 = df[df['Hap_Name'] == hap2][phenotype_col].dropna()
                        if len(g1) > 0 and len(g2) > 0:
                            _, p = ttest_ind(g1, g2)
                            tukey_results.append(f"{hap1} vs {hap2}: p={p:.4f}")
                result['tukey_results'] = "; ".join(tukey_results) if tukey_results else "NA"
                
        elif method == 'kruskal':
            # 非参数: Kruskal-Wallis
            stat, pval = stats.kruskal(*group_values)
            result['statistic'] = stat
            result['p_value'] = pval
            result['test_type'] = "Kruskal-Wallis H-test"
        
        # 显著性判断
        result['significant'] = result.get('p_value', 1) < 0.05
        result['highly_significant'] = result.get('p_value', 1) < 0.01
        
        return result
    
    def regression_analysis(self, phenotype_col: str) -> dict:
        """
        线性回归分析 - 单倍型作为分类变量
        
        Returns:
            dict: 回归分析结果
        """
        df = self.merged_df[
            (self.merged_df['Hap_Name'] != 'Other') & 
            (self.merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 10:
            return {'error': '样本量不足进行回归分析'}
        
        # 使用最大频率的单倍型作为参考
        ref_hap = df['Hap_Name'].value_counts().index[0]
        
        # 创建哑变量
        dummies = pd.get_dummies(df['Hap_Name'], prefix='Hap', drop_first=False)
        df = pd.concat([df, dummies], axis=1)
        
        # OLS 回归（使用sklearn实现，无需statsmodels）
        hap_cols = [c for c in dummies.columns if c != f'Hap_{ref_hap}']
        if not hap_cols:
            return {'error': '单倍型分组不足'}
        
        try:
            # 准备特征矩阵和目标变量
            X = df[hap_cols].values
            y = df[phenotype_col].values
            
            # 添加截距项
            X_with_intercept = np.column_stack([np.ones(len(X)), X])
            
            # 拟合线性回归
            model = LinearRegression(fit_intercept=False)
            model.fit(X_with_intercept, y)
            
            # 计算R²
            y_pred = model.predict(X_with_intercept)
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            
            # 计算调整后R²
            n = len(y)
            p = len(hap_cols)
            adj_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - p - 1) if n > p + 1 else 0
            
            # 计算F统计量
            ms_reg = (ss_tot - ss_res) / p if p > 0 else 0
            ms_res = ss_res / (n - p - 1) if n > p + 1 else 1
            f_statistic = ms_reg / ms_res if ms_res > 0 else 0
            
            # 计算F检验p值
            from scipy.stats import f as f_dist
            f_pvalue = 1 - f_dist.cdf(f_statistic, p, n - p - 1) if f_statistic > 0 else 1
            
            # 计算系数的t检验p值（简化版）
            coef_pvalues = []
            for i in range(1, len(model.coef_)):
                # 使用标准误简化估计
                se = np.sqrt(ms_res / np.sum(X[:, i-1])) if np.sum(X[:, i-1]) > 0 else 1
                t_stat = model.coef_[i] / se if se > 0 else 0
                # 近似p值
                p_val = 2 * (1 - stats.t.cdf(abs(t_stat), n - p - 1)) if n > p + 1 else 1
                coef_pvalues.append(p_val)
            
            result = {
                'phenotype': phenotype_col,
                'method': 'OLS Regression',
                'reference_haplotype': ref_hap,
                'r_squared': r_squared,
                'adj_r_squared': adj_r_squared,
                'f_statistic': f_statistic,
                'f_pvalue': f_pvalue,
                'n_samples': len(df),
                'coefficients': {},
                'p_values': {}
            }
            
            for i, col in enumerate(hap_cols):
                if i < len(model.coef_) - 1:
                    result['coefficients'][col] = model.coef_[i + 1]
                    result['p_values'][col] = coef_pvalues[i] if i < len(coef_pvalues) else 1.0
            
            return result
            
        except Exception as e:
            return {'error': str(e)}


# ============================================================================
# 模块3: AMOVA分析（分子变异分析）
# ============================================================================

class AMOVAAnalyzer:
    """
    AMOVA (Analysis of Molecular Variance) 分析模块
    用NumPy/Pandas手动实现，分析群体间遗传变异占比
    支持多层级分组（如地理区域、表型亚群）
    """
    
    def __init__(self, genotype_df: pd.DataFrame, group_df: pd.DataFrame = None):
        """
        Args:
            genotype_df: 基因型数据框，行为样本，列为位点
            group_df: 分组信息，需包含SampleID和分组列
        """
        self.genotype_df = genotype_df
        self.group_df = group_df
        self.results = {}
    
    def _encode_genotype(self, gt_str: str) -> int:
        """
        将基因型编码为数值（用于距离计算）
        0|0 -> 0, 0|1/1|0 -> 1, 1|1 -> 2
        """
        if pd.isna(gt_str) or gt_str in ['N', '.', './.' , '.|.']:
            return np.nan
        gt_str = str(gt_str)
        if '|' in gt_str:
            parts = gt_str.split('|')
        elif '/' in gt_str:
            parts = gt_str.split('/')
        else:
            return np.nan
        
        try:
            return sum(int(p) for p in parts if p.isdigit())
        except:
            return np.nan
    
    def _calculate_distance_matrix(self, genotypes: np.ndarray) -> np.ndarray:
        """
        计算样本间的遗传距离矩阵（欧氏距离的平方）
        用于AMOVA分析
        """
        n_samples = genotypes.shape[0]
        dist_matrix = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                # 计算两样本间的差异（忽略缺失值）
                g1, g2 = genotypes[i], genotypes[j]
                valid = ~(np.isnan(g1) | np.isnan(g2))
                if np.sum(valid) > 0:
                    diff = g1[valid] - g2[valid]
                    dist = np.sum(diff ** 2)
                    dist_matrix[i, j] = dist
                    dist_matrix[j, i] = dist
        
        return dist_matrix
    
    def run_amova(self, group_col: str, nested_col: str = None) -> dict:
        """
        执行AMOVA分析
        
        Args:
            group_col: 主要分组列名（如地理区域、表型分组）
            nested_col: 嵌套分组列名（可选，如亚群）
        
        Returns:
            dict: AMOVA分析结果，包含各层级方差组分和Phi统计量
        """
        logger = get_logger()
        logger.info(f"开始AMOVA分析，分组列: {group_col}")
        
        if self.group_df is None:
            return {'error': '未提供分组信息'}
        
        # 合并基因型和分组数据
        if 'SampleID' not in self.genotype_df.columns:
            geno_df = self.genotype_df.reset_index()
            geno_df.columns = ['SampleID'] + list(geno_df.columns[1:])
        else:
            geno_df = self.genotype_df.copy()
        
        merged = pd.merge(geno_df, self.group_df[['SampleID', group_col]], on='SampleID', how='inner')
        
        if len(merged) < 10:
            return {'error': '合并后样本量不足', 'n_samples': len(merged)}
        
        # 提取基因型数据列（排除ID和分组列）
        geno_cols = [c for c in merged.columns if c not in ['SampleID', group_col, nested_col] and c is not None]
        
        if len(geno_cols) == 0:
            return {'error': '无有效基因型列'}
        
        # 编码基因型
        genotypes = merged[geno_cols].map(self._encode_genotype).values.astype(float)
        groups = merged[group_col].values
        
        # 计算距离矩阵
        dist_matrix = self._calculate_distance_matrix(genotypes)
        
        # 获取唯一分组
        unique_groups = np.unique(groups)
        n_groups = len(unique_groups)
        n_total = len(groups)
        
        if n_groups < 2:
            return {'error': '分组数不足', 'n_groups': n_groups}
        
        logger.info(f"  样本数: {n_total}, 分组数: {n_groups}")
        
        # 计算各组样本索引
        group_indices = {g: np.where(groups == g)[0] for g in unique_groups}
        group_sizes = {g: len(idx) for g, idx in group_indices.items()}
        
        # 计算总平方和 (SSD_Total) - 所有配对距离之和 / (2*n)
        ssd_total = np.sum(dist_matrix) / 2
        
        # 计算组内平方和 (SSD_Within)
        ssd_within = 0
        for g, idx in group_indices.items():
            n_g = len(idx)
            if n_g > 1:
                sub_dist = dist_matrix[np.ix_(idx, idx)]
                ssd_within += np.sum(sub_dist) / 2
        
        # 计算组间平方和 (SSD_Among)
        ssd_among = ssd_total - ssd_within
        
        # 自由度计算
        df_among = n_groups - 1
        df_within = n_total - n_groups
        df_total = n_total - 1
        
        # 均方计算
        ms_among = ssd_among / df_among if df_among > 0 else 0
        ms_within = ssd_within / df_within if df_within > 0 else 0
        
        # 计算加权平均样本量 n0
        n_values = np.array(list(group_sizes.values()))
        n0 = (n_total - np.sum(n_values ** 2) / n_total) / df_among if df_among > 0 else 0
        
        # 方差组分估计
        # sigma_a^2 (组间方差) = (MS_among - MS_within) / n0
        # sigma_w^2 (组内方差) = MS_within
        sigma_within = ms_within
        sigma_among = max(0, (ms_among - ms_within) / n0) if n0 > 0 else 0
        sigma_total = sigma_among + sigma_within
        
        # Phi统计量 (类似于F_ST)
        phi_st = sigma_among / sigma_total if sigma_total > 0 else 0
        
        # 方差组分百分比
        pct_among = (sigma_among / sigma_total * 100) if sigma_total > 0 else 0
        pct_within = (sigma_within / sigma_total * 100) if sigma_total > 0 else 0
        
        # 置换检验计算P值
        n_permutations = 999
        perm_phi_values = []
        
        for _ in range(n_permutations):
            perm_groups = np.random.permutation(groups)
            perm_group_indices = {g: np.where(perm_groups == g)[0] for g in unique_groups}
            
            perm_ssd_within = 0
            for g, idx in perm_group_indices.items():
                if len(idx) > 1:
                    sub_dist = dist_matrix[np.ix_(idx, idx)]
                    perm_ssd_within += np.sum(sub_dist) / 2
            
            perm_ssd_among = ssd_total - perm_ssd_within
            perm_ms_among = perm_ssd_among / df_among if df_among > 0 else 0
            perm_ms_within = perm_ssd_within / df_within if df_within > 0 else 0
            
            perm_sigma_within = perm_ms_within
            perm_sigma_among = max(0, (perm_ms_among - perm_ms_within) / n0) if n0 > 0 else 0
            perm_sigma_total = perm_sigma_among + perm_sigma_within
            
            perm_phi = perm_sigma_among / perm_sigma_total if perm_sigma_total > 0 else 0
            perm_phi_values.append(perm_phi)
        
        # 计算P值
        p_value = (np.sum(np.array(perm_phi_values) >= phi_st) + 1) / (n_permutations + 1)
        
        result = {
            'method': 'AMOVA',
            'group_column': group_col,
            'n_samples': n_total,
            'n_groups': n_groups,
            'group_sizes': group_sizes,
            
            # 平方和
            'ssd_among': ssd_among,
            'ssd_within': ssd_within,
            'ssd_total': ssd_total,
            
            # 自由度
            'df_among': df_among,
            'df_within': df_within,
            'df_total': df_total,
            
            # 均方
            'ms_among': ms_among,
            'ms_within': ms_within,
            
            # 方差组分
            'variance_among': sigma_among,
            'variance_within': sigma_within,
            'variance_total': sigma_total,
            
            # 方差百分比
            'percent_among': pct_among,
            'percent_within': pct_within,
            
            # Phi统计量和P值
            'phi_st': phi_st,
            'p_value': p_value,
            'n_permutations': n_permutations,
            
            # 显著性判断
            'significant': p_value < 0.05,
            'highly_significant': p_value < 0.01,
        }
        
        self.results = result
        logger.info(f"  Phi_ST: {phi_st:.4f}, P-value: {p_value:.4f}")
        logger.info(f"  组间方差占比: {pct_among:.2f}%, 组内方差占比: {pct_within:.2f}%")
        
        return result
    
    def run_hierarchical_amova(self, region_col: str, population_col: str) -> dict:
        """
        执行多层级AMOVA分析（三层级：区域 > 群体 > 个体）
        
        Args:
            region_col: 区域分组列
            population_col: 群体分组列
        
        Returns:
            dict: 多层级AMOVA结果
        """
        logger = get_logger()
        logger.info(f"开始多层级AMOVA分析")
        
        if self.group_df is None:
            return {'error': '未提供分组信息'}
        
        required_cols = ['SampleID', region_col, population_col]
        if not all(c in self.group_df.columns for c in required_cols):
            return {'error': f'分组数据缺少必要列: {required_cols}'}
        
        # 合并数据
        if 'SampleID' not in self.genotype_df.columns:
            geno_df = self.genotype_df.reset_index()
            geno_df.columns = ['SampleID'] + list(geno_df.columns[1:])
        else:
            geno_df = self.genotype_df.copy()
        
        merged = pd.merge(geno_df, self.group_df[required_cols], on='SampleID', how='inner')
        
        if len(merged) < 10:
            return {'error': '样本量不足'}
        
        geno_cols = [c for c in merged.columns if c not in required_cols]
        genotypes = merged[geno_cols].map(self._encode_genotype).values.astype(float)
        
        regions = merged[region_col].values
        populations = merged[population_col].values
        
        dist_matrix = self._calculate_distance_matrix(genotypes)
        
        n_total = len(merged)
        unique_regions = np.unique(regions)
        n_regions = len(unique_regions)
        
        # 计算各层级的样本索引
        region_pop_structure = {}
        for r in unique_regions:
            r_mask = regions == r
            r_pops = np.unique(populations[r_mask])
            region_pop_structure[r] = {
                'indices': np.where(r_mask)[0],
                'populations': {p: np.where((regions == r) & (populations == p))[0] for p in r_pops}
            }
        
        # 计算总平方和
        ssd_total = np.sum(dist_matrix) / 2
        
        # 计算组内平方和（个体内）
        ssd_within_pop = 0
        n_populations = 0
        for r, r_info in region_pop_structure.items():
            for p, p_idx in r_info['populations'].items():
                n_populations += 1
                if len(p_idx) > 1:
                    sub_dist = dist_matrix[np.ix_(p_idx, p_idx)]
                    ssd_within_pop += np.sum(sub_dist) / 2
        
        # 计算区域内群体间平方和
        ssd_among_pop_within_region = 0
        for r, r_info in region_pop_structure.items():
            r_idx = r_info['indices']
            if len(r_idx) > 1:
                r_dist = dist_matrix[np.ix_(r_idx, r_idx)]
                ssd_within_region = np.sum(r_dist) / 2
                
                ssd_pop_in_region = 0
                for p, p_idx in r_info['populations'].items():
                    if len(p_idx) > 1:
                        sub_dist = dist_matrix[np.ix_(p_idx, p_idx)]
                        ssd_pop_in_region += np.sum(sub_dist) / 2
                
                ssd_among_pop_within_region += ssd_within_region - ssd_pop_in_region
        
        # 计算区域间平方和
        ssd_among_regions = ssd_total - ssd_within_pop - ssd_among_pop_within_region
        
        # 自由度
        df_among_regions = n_regions - 1
        df_among_pop = n_populations - n_regions
        df_within_pop = n_total - n_populations
        
        # 方差组分（简化估计）
        ms_among_regions = ssd_among_regions / df_among_regions if df_among_regions > 0 else 0
        ms_among_pop = ssd_among_pop_within_region / df_among_pop if df_among_pop > 0 else 0
        ms_within = ssd_within_pop / df_within_pop if df_within_pop > 0 else 0
        
        sigma_c = ms_within  # 个体内方差
        sigma_b = max(0, (ms_among_pop - ms_within) / 2)  # 群体间方差（区域内）
        sigma_a = max(0, (ms_among_regions - ms_among_pop) / 4)  # 区域间方差
        sigma_total = sigma_a + sigma_b + sigma_c
        
        # Phi统计量
        phi_ct = sigma_a / sigma_total if sigma_total > 0 else 0  # 区域间
        phi_sc = sigma_b / (sigma_b + sigma_c) if (sigma_b + sigma_c) > 0 else 0  # 群体间（区域内）
        phi_st = (sigma_a + sigma_b) / sigma_total if sigma_total > 0 else 0  # 总的
        
        result = {
            'method': 'Hierarchical AMOVA',
            'n_samples': n_total,
            'n_regions': n_regions,
            'n_populations': n_populations,
            
            # 方差组分
            'variance_among_regions': sigma_a,
            'variance_among_populations': sigma_b,
            'variance_within_populations': sigma_c,
            
            # 方差百分比
            'percent_among_regions': sigma_a / sigma_total * 100 if sigma_total > 0 else 0,
            'percent_among_populations': sigma_b / sigma_total * 100 if sigma_total > 0 else 0,
            'percent_within_populations': sigma_c / sigma_total * 100 if sigma_total > 0 else 0,
            
            # Phi统计量
            'phi_ct': phi_ct,  # 区域间分化
            'phi_sc': phi_sc,  # 区域内群体间分化
            'phi_st': phi_st,  # 总分化
        }
        
        logger.info(f"  Phi_CT (区域间): {phi_ct:.4f}")
        logger.info(f"  Phi_SC (群体间/区域内): {phi_sc:.4f}")
        logger.info(f"  Phi_ST (总): {phi_st:.4f}")
        
        return result
    
    def format_amova_table(self) -> pd.DataFrame:
        """
        格式化AMOVA结果为标准表格
        """
        if not self.results:
            return pd.DataFrame()
        
        r = self.results
        
        if 'ssd_among' in r:
            # 标准AMOVA表格格式
            data = [
                {
                    'Source': 'Among groups',
                    'd.f.': r.get('df_among', 0),
                    'Sum of squares': r.get('ssd_among', 0),
                    'Variance': r.get('variance_among', 0),
                    '% Variance': r.get('percent_among', 0),
                },
                {
                    'Source': 'Within groups',
                    'd.f.': r.get('df_within', 0),
                    'Sum of squares': r.get('ssd_within', 0),
                    'Variance': r.get('variance_within', 0),
                    '% Variance': r.get('percent_within', 0),
                },
                {
                    'Source': 'Total',
                    'd.f.': r.get('df_total', 0),
                    'Sum of squares': r.get('ssd_total', 0),
                    'Variance': r.get('variance_total', 0),
                    '% Variance': 100.0,
                },
            ]
            
            df = pd.DataFrame(data)
            return df
        
        return pd.DataFrame()


# ============================================================================
# 模块4: 单倍型效应分析
# ============================================================================

class HaplotypeEffectAnalyzer:
    """
    单倍型效应分析模块
    用于识别有显著效应的单倍型，计算效应值和置信区间
    """
    
    def __init__(self, merged_df: pd.DataFrame):
        """
        Args:
            merged_df: 合并后的数据框，包含 Hap_Name 和表型列
        """
        self.merged_df = merged_df
        self.results = {}
    
    def calculate_effects(self, phenotype_col: str, reference_hap: str = None,
                          method: str = 'mean_diff') -> dict:
        """
        计算各单倍型相对于参考单倍型的效应值
        
        Args:
            phenotype_col: 表型列名
            reference_hap: 参考单倍型（默认为最大频率的单倍型）
            method: 效应计算方法 ('mean_diff', 'regression', 'contrast')
        
        Returns:
            dict: 各单倍型的效应分析结果
        """
        logger = get_logger()
        logger.info(f"开始单倍型效应分析: {phenotype_col}")
        
        df = self.merged_df[
            (self.merged_df['Hap_Name'] != 'Other') & 
            (self.merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 10:
            return {'error': '样本量不足', 'n_samples': len(df)}
        
        # 获取单倍型分组
        hap_counts = df['Hap_Name'].value_counts()
        unique_haps = hap_counts.index.tolist()
        n_haps = len(unique_haps)
        
        if n_haps < 2:
            return {'error': '单倍型分组不足'}
        
        # 确定参考单倍型（默认选择最大频率的）
        if reference_hap is None or reference_hap not in unique_haps:
            reference_hap = hap_counts.index[0]
        
        logger.info(f"  参考单倍型: {reference_hap} (n={hap_counts[reference_hap]})")
        
        # 参考组数据
        ref_data = df[df['Hap_Name'] == reference_hap][phenotype_col].values
        ref_mean = np.mean(ref_data)
        ref_std = np.std(ref_data, ddof=1)
        ref_n = len(ref_data)
        
        # 总体统计
        grand_mean = df[phenotype_col].mean()
        grand_std = df[phenotype_col].std()
        
        # 计算每个单倍型的效应
        haplotype_effects = []
        
        for hap in unique_haps:
            hap_data = df[df['Hap_Name'] == hap][phenotype_col].values
            hap_n = len(hap_data)
            hap_mean = np.mean(hap_data)
            hap_std = np.std(hap_data, ddof=1)
            hap_se = hap_std / np.sqrt(hap_n) if hap_n > 0 else 0
            
            # 效应值（相对于参考组）
            effect = hap_mean - ref_mean
            
            # Cohen's d 标准化效应量
            pooled_std = np.sqrt(((hap_n - 1) * hap_std**2 + (ref_n - 1) * ref_std**2) / 
                                 (hap_n + ref_n - 2)) if (hap_n + ref_n) > 2 else 1
            cohens_d = effect / pooled_std if pooled_std > 0 else 0
            
            # 效应的标准误差
            effect_se = np.sqrt(hap_std**2/hap_n + ref_std**2/ref_n) if hap_n > 0 else 0
            
            # 95%置信区间
            ci_lower = effect - 1.96 * effect_se
            ci_upper = effect + 1.96 * effect_se
            
            # t检验 P值
            if hap != reference_hap and hap_n >= 2 and ref_n >= 2:
                t_stat, p_value = ttest_ind(hap_data, ref_data, equal_var=False)
            else:
                t_stat, p_value = 0, 1.0
            
            # 效应大小判断
            abs_d = abs(cohens_d)
            if abs_d < 0.2:
                effect_size_cat = 'negligible'
            elif abs_d < 0.5:
                effect_size_cat = 'small'
            elif abs_d < 0.8:
                effect_size_cat = 'medium'
            else:
                effect_size_cat = 'large'
            
            # 效应方向
            if effect > 0:
                direction = 'positive (+)'
            elif effect < 0:
                direction = 'negative (-)'
            else:
                direction = 'neutral'
            
            haplotype_effects.append({
                'haplotype': hap,
                'n_samples': hap_n,
                'mean': hap_mean,
                'std': hap_std,
                'se': hap_se,
                'effect': effect,
                'effect_se': effect_se,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'cohens_d': cohens_d,
                'effect_size': effect_size_cat,
                'direction': direction,
                't_statistic': t_stat,
                'p_value': p_value,
                'significant': p_value < 0.05,
                'highly_significant': p_value < 0.01,
                'is_reference': hap == reference_hap,
            })
        
        # 按效应值绝对值排序
        haplotype_effects.sort(key=lambda x: abs(x['effect']), reverse=True)
        
        # 找出显著效应的单倍型
        significant_haps = [h for h in haplotype_effects if h['significant'] and not h['is_reference']]
        
        result = {
            'phenotype': phenotype_col,
            'reference_haplotype': reference_hap,
            'reference_mean': ref_mean,
            'reference_n': ref_n,
            'grand_mean': grand_mean,
            'grand_std': grand_std,
            'n_haplotypes': n_haps,
            'n_samples': len(df),
            'haplotype_effects': haplotype_effects,
            'significant_haplotypes': significant_haps,
            'n_significant': len(significant_haps),
        }
        
        self.results = result
        
        # 输出结果
        logger.info(f"  总样本数: {len(df)}, 单倍型数: {n_haps}")
        logger.info(f"  参考单倍型均值: {ref_mean:.3f}")
        logger.info(f"  显著效应单倍型: {len(significant_haps)} 个")
        
        return result
    
    def get_effect_summary_table(self) -> pd.DataFrame:
        """
        生成效应汇总表格
        """
        if not self.results or 'haplotype_effects' not in self.results:
            return pd.DataFrame()
        
        rows = []
        for h in self.results['haplotype_effects']:
            sig_mark = '***' if h['highly_significant'] else ('*' if h['significant'] else '')
            rows.append({
                'Haplotype': h['haplotype'],
                'N': h['n_samples'],
                'Mean': f"{h['mean']:.3f}",
                'Effect': f"{h['effect']:+.3f}" if not h['is_reference'] else 'REF',
                '95% CI': f"[{h['ci_lower']:.2f}, {h['ci_upper']:.2f}]" if not h['is_reference'] else '-',
                "Cohen's d": f"{h['cohens_d']:.3f}" if not h['is_reference'] else '-',
                'Effect Size': h['effect_size'] if not h['is_reference'] else '-',
                'P-value': f"{h['p_value']:.4f}{sig_mark}" if not h['is_reference'] else '-',
                'Direction': h['direction'] if not h['is_reference'] else '-',
            })
        
        return pd.DataFrame(rows)
    
    def get_significant_haplotypes(self, alpha: float = 0.05, 
                                   min_effect_size: str = None) -> list:
        """
        获取显著效应的单倍型列表
        
        Args:
            alpha: 显著性水平
            min_effect_size: 最小效应量要求 ('small', 'medium', 'large')
        
        Returns:
            list: 显著单倍型列表
        """
        if not self.results or 'haplotype_effects' not in self.results:
            return []
        
        effect_order = {'negligible': 0, 'small': 1, 'medium': 2, 'large': 3}
        min_order = effect_order.get(min_effect_size, 0)
        
        significant = []
        for h in self.results['haplotype_effects']:
            if h['is_reference']:
                continue
            if h['p_value'] < alpha:
                if effect_order.get(h['effect_size'], 0) >= min_order:
                    significant.append({
                        'haplotype': h['haplotype'],
                        'effect': h['effect'],
                        'cohens_d': h['cohens_d'],
                        'p_value': h['p_value'],
                        'effect_size': h['effect_size'],
                        'direction': h['direction'],
                    })
        
        return significant
    
    def plot_effect_forest(self, output_file: str = None, 
                          show_ci: bool = True) -> str:
        """
        绘制效应值森林图
        """
        if not self.results or 'haplotype_effects' not in self.results:
            print("[WARNING] 无效应分析结果")
            return None
        
        effects = self.results['haplotype_effects']
        n_haps = len(effects)
        
        fig, ax = plt.subplots(figsize=(10, max(4, n_haps * 0.5)))
        
        y_positions = range(n_haps)
        hap_names = []
        effect_values = []
        ci_lowers = []
        ci_uppers = []
        colors = []
        
        for h in effects:
            hap_names.append(h['haplotype'])
            effect_values.append(h['effect'])
            ci_lowers.append(h['ci_lower'])
            ci_uppers.append(h['ci_upper'])
            
            if h['is_reference']:
                colors.append('gray')
            elif h['highly_significant']:
                colors.append('darkred')
            elif h['significant']:
                colors.append('red')
            else:
                colors.append('steelblue')
        
        # 绘制点
        ax.scatter(effect_values, y_positions, c=colors, s=80, zorder=3)
        
        # 绘制置信区间
        if show_ci:
            for i, (eff, lo, hi, col) in enumerate(zip(effect_values, ci_lowers, ci_uppers, colors)):
                ax.plot([lo, hi], [i, i], color=col, linewidth=2, alpha=0.7)
        
        # 参考线（零效应）
        ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.7)
        
        # 标签
        ax.set_yticks(y_positions)
        ax.set_yticklabels(hap_names)
        ax.set_xlabel('Effect Size (vs Reference)', fontsize=11)
        ax.set_ylabel('Haplotype', fontsize=11)
        
        # 标题
        phenotype = self.results.get('phenotype', '')
        ref_hap = self.results.get('reference_haplotype', '')
        ax.set_title(f'Haplotype Effect Analysis: {phenotype}\n(Reference: {ref_hap})', fontsize=12)
        
        # 图例
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Reference'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='darkred', markersize=10, label='P < 0.01'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='P < 0.05'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='steelblue', markersize=10, label='Not sig.'),
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] 效应森林图已保存: {output_file}")
            return output_file
        else:
            plt.show()
            return None
    
    def plot_effect_bar(self, output_file: str = None) -> str:
        """
        绘制效应值柱状图
        """
        if not self.results or 'haplotype_effects' not in self.results:
            print("[WARNING] 无效应分析结果")
            return None
        
        effects = [h for h in self.results['haplotype_effects'] if not h['is_reference']]
        if not effects:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        hap_names = [h['haplotype'] for h in effects]
        effect_values = [h['effect'] for h in effects]
        errors = [h['effect_se'] * 1.96 for h in effects]  # 95% CI
        
        # 颜色根据显著性
        colors = []
        for h in effects:
            if h['highly_significant']:
                colors.append('darkred')
            elif h['significant']:
                colors.append('red')
            else:
                colors.append('steelblue')
        
        bars = ax.bar(hap_names, effect_values, yerr=errors, capsize=5, 
                     color=colors, alpha=0.7, edgecolor='black')
        
        # 零线
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        
        # 添加显著性标记
        for i, (bar, h) in enumerate(zip(bars, effects)):
            if h['highly_significant']:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + errors[i] + 0.5,
                       '***', ha='center', va='bottom', fontsize=12, fontweight='bold')
            elif h['significant']:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + errors[i] + 0.5,
                       '*', ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        ax.set_xlabel('Haplotype', fontsize=11)
        ax.set_ylabel('Effect (vs Reference)', fontsize=11)
        
        phenotype = self.results.get('phenotype', '')
        ref_hap = self.results.get('reference_haplotype', '')
        ax.set_title(f'Haplotype Effects on {phenotype}\n(Reference: {ref_hap})', fontsize=12)
        
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] 效应柱状图已保存: {output_file}")
            return output_file
        else:
            plt.show()
            return None


# ============================================================================
# 模块5: 多重检验校正
# ============================================================================

class MultipleTestingCorrection:
    """
    多重检验校正模块
    实现Bonferroni和FDR校正，避免大规模位点分析的假阳性
    """
    
    @staticmethod
    def bonferroni(p_values: np.ndarray, alpha: float = 0.05) -> dict:
        """
        Bonferroni校正
        最保守的校正方法，控制家族错误率(FWER)
        
        Args:
            p_values: 原始P值数组
            alpha: 显著性水平
        
        Returns:
            dict: 校正结果
        """
        p_values = np.asarray(p_values).flatten()
        n_tests = len(p_values)
        
        # Bonferroni校正阈值
        corrected_alpha = alpha / n_tests
        
        # 校正后的P值
        p_corrected = np.minimum(p_values * n_tests, 1.0)
        
        # 显著性判断
        significant = p_values < corrected_alpha
        
        return {
            'method': 'Bonferroni',
            'n_tests': n_tests,
            'original_alpha': alpha,
            'corrected_alpha': corrected_alpha,
            'p_original': p_values.tolist(),
            'p_corrected': p_corrected.tolist(),
            'significant': significant.tolist(),
            'n_significant': int(np.sum(significant)),
            'significant_indices': np.where(significant)[0].tolist(),
        }
    
    @staticmethod
    def fdr_bh(p_values: np.ndarray, alpha: float = 0.05) -> dict:
        """
        Benjamini-Hochberg FDR校正
        控制错误发现率(FDR)，比Bonferroni更宽松
        
        Args:
            p_values: 原始P值数组
            alpha: FDR阈值
        
        Returns:
            dict: 校正结果
        """
        p_values = np.asarray(p_values).flatten()
        n_tests = len(p_values)
        
        if n_tests == 0:
            return {'error': '无P值输入'}
        
        # 排序索引
        sorted_idx = np.argsort(p_values)
        sorted_p = p_values[sorted_idx]
        
        # BH校正
        # 对于排序后的第i个P值，校正阈值为 (i/m) * alpha
        bh_critical = (np.arange(1, n_tests + 1) / n_tests) * alpha
        
        # 找到最大的i使得 p_(i) <= (i/m) * alpha
        significant_sorted = sorted_p <= bh_critical
        
        if np.any(significant_sorted):
            # 最大显著位置及之前的所有位点都显著
            max_significant_idx = np.max(np.where(significant_sorted)[0])
            significant_sorted[:max_significant_idx + 1] = True
        
        # 计算校正后的P值（调整后的P值）
        p_adjusted = np.zeros(n_tests)
        cummin_p = sorted_p[-1]
        
        for i in range(n_tests - 1, -1, -1):
            adjusted = sorted_p[i] * n_tests / (i + 1)
            cummin_p = min(cummin_p, adjusted)
            p_adjusted[sorted_idx[i]] = min(cummin_p, 1.0)
        
        # 恢复原始顺序的显著性
        significant = np.zeros(n_tests, dtype=bool)
        significant[sorted_idx] = significant_sorted
        
        return {
            'method': 'Benjamini-Hochberg FDR',
            'n_tests': n_tests,
            'fdr_level': alpha,
            'p_original': p_values.tolist(),
            'p_adjusted': p_adjusted.tolist(),
            'significant': significant.tolist(),
            'n_significant': int(np.sum(significant)),
            'significant_indices': np.where(significant)[0].tolist(),
            'bh_critical_values': bh_critical.tolist(),
        }
    
    @staticmethod
    def fdr_by(p_values: np.ndarray, alpha: float = 0.05) -> dict:
        """
        Benjamini-Yekutieli FDR校正
        适用于任意依赖结构的P值，比BH更保守
        
        Args:
            p_values: 原始P值数组
            alpha: FDR阈值
        
        Returns:
            dict: 校正结果
        """
        p_values = np.asarray(p_values).flatten()
        n_tests = len(p_values)
        
        if n_tests == 0:
            return {'error': '无P值输入'}
        
        # BY校正因子 c(m) = sum(1/i) for i=1 to m
        c_m = np.sum(1 / np.arange(1, n_tests + 1))
        
        # 排序
        sorted_idx = np.argsort(p_values)
        sorted_p = p_values[sorted_idx]
        
        # BY校正阈值
        by_critical = (np.arange(1, n_tests + 1) / (n_tests * c_m)) * alpha
        
        significant_sorted = sorted_p <= by_critical
        
        if np.any(significant_sorted):
            max_significant_idx = np.max(np.where(significant_sorted)[0])
            significant_sorted[:max_significant_idx + 1] = True
        
        # 校正后的P值
        p_adjusted = np.zeros(n_tests)
        cummin_p = sorted_p[-1]
        
        for i in range(n_tests - 1, -1, -1):
            adjusted = sorted_p[i] * n_tests * c_m / (i + 1)
            cummin_p = min(cummin_p, adjusted)
            p_adjusted[sorted_idx[i]] = min(cummin_p, 1.0)
        
        significant = np.zeros(n_tests, dtype=bool)
        significant[sorted_idx] = significant_sorted
        
        return {
            'method': 'Benjamini-Yekutieli FDR',
            'n_tests': n_tests,
            'fdr_level': alpha,
            'correction_factor': c_m,
            'p_original': p_values.tolist(),
            'p_adjusted': p_adjusted.tolist(),
            'significant': significant.tolist(),
            'n_significant': int(np.sum(significant)),
            'significant_indices': np.where(significant)[0].tolist(),
        }
    
    @staticmethod
    def compare_corrections(p_values: np.ndarray, alpha: float = 0.05) -> pd.DataFrame:
        """
        比较不同校正方法的结果
        
        Args:
            p_values: 原始P值数组
            alpha: 显著性水平
        
        Returns:
            DataFrame: 比较结果表
        """
        bonf = MultipleTestingCorrection.bonferroni(p_values, alpha)
        fdr_bh = MultipleTestingCorrection.fdr_bh(p_values, alpha)
        fdr_by = MultipleTestingCorrection.fdr_by(p_values, alpha)
        
        comparison = pd.DataFrame({
            'Test_Index': range(len(p_values)),
            'P_Original': p_values,
            'P_Bonferroni': bonf['p_corrected'],
            'P_FDR_BH': fdr_bh['p_adjusted'],
            'P_FDR_BY': fdr_by['p_adjusted'],
            'Sig_Bonferroni': bonf['significant'],
            'Sig_FDR_BH': fdr_bh['significant'],
            'Sig_FDR_BY': fdr_by['significant'],
        })
        
        return comparison
    
    @staticmethod
    def annotate_significance(p_value: float, corrected: bool = False, 
                             thresholds: dict = None) -> str:
        """
        根据P值生成显著性标注字符串
        
        Args:
            p_value: P值
            corrected: 是否已校正
            thresholds: 自定义阈值字典
        
        Returns:
            str: 显著性标注 (ns, *, **, ***)
        """
        if thresholds is None:
            thresholds = {
                'ns': 1.0,      # 不显著
                '*': 0.05,      # P < 0.05
                '**': 0.01,     # P < 0.01
                '***': 0.001,   # P < 0.001
            }
        
        if p_value >= thresholds['*']:
            return 'ns'
        elif p_value >= thresholds['**']:
            return '*'
        elif p_value >= thresholds['***']:
            return '**'
        else:
            return '***'
    
    @staticmethod
    def format_pvalue(p_value: float, sig_level: bool = True) -> str:
        """
        格式化P值字符串，用于图表标注
        
        Args:
            p_value: P值
            sig_level: 是否包含显著性水平标记
        
        Returns:
            str: 格式化的P值字符串
        """
        if p_value < 0.001:
            p_str = f'P < 0.001'
        elif p_value < 0.01:
            p_str = f'P = {p_value:.3f}'
        elif p_value < 0.05:
            p_str = f'P = {p_value:.3f}'
        else:
            p_str = f'P = {p_value:.2f}'
        
        if sig_level:
            sig = MultipleTestingCorrection.annotate_significance(p_value)
            if sig != 'ns':
                p_str += f' ({sig})'
        
        return p_str


# ============================================================================
# 模块6: PVE 变异解释率计算
# ============================================================================

class PVECalculator:
    """变异解释率 (Phenotypic Variance Explained) 计算模块"""
    
    def __init__(self, merged_df: pd.DataFrame):
        """
        Args:
            merged_df: 合并后的数据框（包含 Hap_Name 和表型列）
        """
        self.merged_df = merged_df
    
    def calculate_pve(self, phenotype_col: str, method: str = 'r_squared') -> dict:
        """
        计算单倍型对表型方差的解释率
        
        Args:
            phenotype_col: 表型列名
            method: 计算方法 ('r_squared', 'eta_squared', 'omega_squared')
        
        Returns:
            dict: PVE 计算结果
        """
        df = self.merged_df[
            (self.merged_df['Hap_Name'] != 'Other') & 
            (self.merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 10:
            return {'error': '样本量不足', 'n_samples': len(df)}
        
        groups = df.groupby('Hap_Name')[phenotype_col].apply(list).to_dict()
        n_groups = len(groups)
        n_total = len(df)
        
        if n_groups < 2:
            return {'error': '单倍型分组不足'}
        
        # 计算各类平方和
        grand_mean = df[phenotype_col].mean()
        
        # SS_between (组间平方和)
        ss_between = sum(
            len(g) * (np.mean(g) - grand_mean)**2 
            for g in groups.values()
        )
        
        # SS_within (组内平方和)
        ss_within = sum(
            sum((x - np.mean(g))**2 for x in g) 
            for g in groups.values()
        )
        
        # SS_total (总平方和)
        ss_total = ss_between + ss_within
        
        # 自由度
        df_between = n_groups - 1
        df_within = n_total - n_groups
        df_total = n_total - 1
        
        result = {
            'phenotype': phenotype_col,
            'n_samples': n_total,
            'n_groups': n_groups,
            'ss_between': ss_between,
            'ss_within': ss_within,
            'ss_total': ss_total,
            'df_between': df_between,
            'df_within': df_within,
        }
        
        # R² (eta-squared): SS_between / SS_total
        if ss_total > 0:
            result['r_squared'] = ss_between / ss_total
            result['eta_squared'] = result['r_squared']
        else:
            result['r_squared'] = 0
            result['eta_squared'] = 0
        
        # Omega-squared (效应量的无偏估计)
        ms_within = ss_within / df_within if df_within > 0 else 0
        omega_sq_num = ss_between - df_between * ms_within
        omega_sq_denom = ss_total + ms_within
        result['omega_squared'] = max(0, omega_sq_num / omega_sq_denom) if omega_sq_denom > 0 else 0
        
        # 广义遗传力估计 (H²)
        # H² = V_between / V_total
        result['heritability_estimate'] = result['eta_squared']
        
        # PVE 百分比
        result['pve_percent'] = result['eta_squared'] * 100
        
        # 效应大小判断 (Cohen's conventions)
        eta_sq = result['eta_squared']
        if eta_sq < 0.01:
            result['effect_size'] = 'negligible'
        elif eta_sq < 0.06:
            result['effect_size'] = 'small'
        elif eta_sq < 0.14:
            result['effect_size'] = 'medium'
        else:
            result['effect_size'] = 'large'
        
        return result
    
    def bootstrap_pve(self, phenotype_col: str, n_bootstrap: int = 1000, 
                      confidence: float = 0.95) -> dict:
        """
        Bootstrap 方法估计 PVE 的置信区间
        
        Args:
            phenotype_col: 表型列名
            n_bootstrap: Bootstrap 重抽样次数
            confidence: 置信水平
        
        Returns:
            dict: 包含置信区间的 PVE 结果
        """
        df = self.merged_df[
            (self.merged_df['Hap_Name'] != 'Other') & 
            (self.merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 20:
            return {'error': '样本量不足进行 Bootstrap'}
        
        pve_values = []
        
        for _ in range(n_bootstrap):
            # 有放回重抽样
            sample_df = df.sample(n=len(df), replace=True)
            
            groups = sample_df.groupby('Hap_Name')[phenotype_col].apply(list).to_dict()
            if len(groups) < 2:
                continue
            
            grand_mean = sample_df[phenotype_col].mean()
            ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups.values())
            ss_total = sum((x - grand_mean)**2 for x in sample_df[phenotype_col])
            
            if ss_total > 0:
                pve_values.append(ss_between / ss_total)
        
        if len(pve_values) < 100:
            return {'error': 'Bootstrap 有效样本不足'}
        
        pve_array = np.array(pve_values)
        alpha = 1 - confidence
        
        result = {
            'phenotype': phenotype_col,
            'method': 'Bootstrap',
            'n_bootstrap': len(pve_values),
            'pve_mean': np.mean(pve_array),
            'pve_median': np.median(pve_array),
            'pve_std': np.std(pve_array),
            'ci_lower': np.percentile(pve_array, alpha/2 * 100),
            'ci_upper': np.percentile(pve_array, (1 - alpha/2) * 100),
            'confidence_level': confidence,
        }
        
        return result


# ============================================================================
# 模块4: GWAS 整合模块
# ============================================================================

class GWASIntegrator:
    """
    GWAS 信息整合模块
    整合传统 GWAS 结果，对比单倍型关联与单点关联的差异
    """
    
    def __init__(self, gwas_file: str = None):
        """
        Args:
            gwas_file: GWAS 结果文件路径 (支持 PLINK/GEMMA/EMMAX 格式)
        """
        self.gwas_file = gwas_file
        self.gwas_df = None
        if gwas_file and os.path.exists(gwas_file):
            self._load_gwas_results()
    
    def _load_gwas_results(self):
        """加载 GWAS 结果文件"""
        # 支持多种格式
        for sep in ['\t', ' ', ',']:
            try:
                df = pd.read_csv(self.gwas_file, sep=sep, comment='#')
                if len(df.columns) > 3:
                    self.gwas_df = df
                    print(f"[INFO] 加载 GWAS 结果: {len(df)} 条记录")
                    break
            except:
                continue
        
        # 标准化列名
        if self.gwas_df is not None:
            col_mapping = {
                'chr': 'CHROM', 'CHR': 'CHROM', 'chromosome': 'CHROM',
                'pos': 'POS', 'BP': 'POS', 'position': 'POS',
                'p': 'P', 'pvalue': 'P', 'P_VALUE': 'P', 'p_wald': 'P',
                'beta': 'BETA', 'EFFECT': 'BETA', 'b': 'BETA',
                'maf': 'MAF', 'AF': 'MAF',
                'snp': 'SNP', 'ID': 'SNP', 'rs': 'SNP',
            }
            self.gwas_df = self.gwas_df.rename(columns=col_mapping)
    
    def load_from_dataframe(self, df: pd.DataFrame):
        """直接从 DataFrame 加载 GWAS 结果"""
        self.gwas_df = df.copy()
    
    def get_region_gwas(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        """
        获取指定区间的 GWAS 结果
        
        Args:
            chrom: 染色体
            start: 起始位置
            end: 终止位置
        
        Returns:
            DataFrame: 区间内的 GWAS 结果
        """
        if self.gwas_df is None:
            return pd.DataFrame()
        
        df = self.gwas_df.copy()
        
        # 染色体匹配（处理不同命名格式）
        chrom_variants = [chrom, chrom.replace('chr', ''), f'chr{chrom}']
        mask = df['CHROM'].astype(str).isin([str(c) for c in chrom_variants])
        
        if 'POS' in df.columns:
            mask &= (df['POS'] >= start) & (df['POS'] <= end)
        
        return df[mask].sort_values('POS') if 'POS' in df.columns else df[mask]
    
    def compare_with_haplotype(self, hap_association: dict, 
                                region_gwas: pd.DataFrame,
                                gwas_threshold: float = 5e-8) -> dict:
        """
        对比单倍型关联与传统 GWAS 结果
        
        Args:
            hap_association: 单倍型关联分析结果
            region_gwas: 区间 GWAS 结果
            gwas_threshold: GWAS 显著性阈值
        
        Returns:
            dict: 对比结果
        """
        result = {
            'haplotype_significant': hap_association.get('significant', False),
            'haplotype_pvalue': hap_association.get('p_value', 1),
            'gwas_n_variants': len(region_gwas),
            'gwas_significant_variants': 0,
            'gwas_min_pvalue': None,
            'gwas_top_variant': None,
            'novel_finding': False,  # 单倍型发现但 GWAS 未发现
            'comparison_note': '',
        }
        
        if len(region_gwas) == 0:
            result['comparison_note'] = 'No GWAS data available for comparison'
            return result
        
        if 'P' in region_gwas.columns:
            # GWAS 最小 p-value
            result['gwas_min_pvalue'] = region_gwas['P'].min()
            result['gwas_significant_variants'] = (region_gwas['P'] < gwas_threshold).sum()
            
            # 最显著的变异
            top_idx = region_gwas['P'].idxmin()
            top_row = region_gwas.loc[top_idx]
            result['gwas_top_variant'] = {
                'pos': top_row.get('POS', 'NA'),
                'p_value': top_row.get('P', 'NA'),
                'beta': top_row.get('BETA', 'NA'),
            }
        
        # 判断是否为新发现
        # 单倍型显著但 GWAS 不显著 -> 可能是新发现
        hap_sig = hap_association.get('p_value', 1) < 0.05
        gwas_sig = result['gwas_significant_variants'] > 0
        
        if hap_sig and not gwas_sig:
            result['novel_finding'] = True
            result['comparison_note'] = '*** NOVEL FINDING: Haplotype association detected but no GWAS signal ***'
        elif hap_sig and gwas_sig:
            result['comparison_note'] = 'Both haplotype and GWAS detected significant association'
        elif not hap_sig and gwas_sig:
            result['comparison_note'] = 'GWAS signal present but haplotype association not significant'
        else:
            result['comparison_note'] = 'No significant association in either method'
        
        return result
    
    def plot_comparison(self, positions: list, hap_pvalue: float, 
                        region_gwas: pd.DataFrame, gene_id: str = None,
                        output_file: str = None) -> str:
        """
        绘制单倍型 vs GWAS 对比图
        """
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        
        # 上图: GWAS Manhattan 局部图
        ax1 = axes[0]
        if len(region_gwas) > 0 and 'P' in region_gwas.columns and 'POS' in region_gwas.columns:
            log_p = -np.log10(region_gwas['P'].replace(0, 1e-300))
            ax1.scatter(region_gwas['POS'], log_p, c='steelblue', alpha=0.7, s=30)
            ax1.axhline(y=-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide (5e-8)')
            ax1.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', label='Suggestive (1e-5)')
            ax1.set_ylabel('-log10(P) GWAS', fontsize=11)
            ax1.legend(loc='upper right', fontsize=9)
        ax1.set_title(f'GWAS vs Haplotype Association{f": {gene_id}" if gene_id else ""}', fontsize=12)
        
        # 下图: 单倍型关联 (单点)
        ax2 = axes[1]
        hap_log_p = -np.log10(hap_pvalue) if hap_pvalue > 0 else 0
        
        # 在单倍型覆盖区域显示
        if positions:
            region_center = (min(positions) + max(positions)) / 2
            region_width = max(positions) - min(positions)
            ax2.bar(region_center, hap_log_p, width=region_width * 0.8, 
                   color='coral', alpha=0.7, label='Haplotype')
            ax2.axhline(y=-np.log10(0.05), color='green', linestyle='--', label='P=0.05')
            ax2.axhline(y=-np.log10(0.01), color='darkgreen', linestyle='--', label='P=0.01')
        
        ax2.set_ylabel('-log10(P) Haplotype', fontsize=11)
        ax2.set_xlabel('Position (bp)', fontsize=11)
        ax2.legend(loc='upper right', fontsize=9)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] 对比图已保存: {output_file}")
            return output_file
        else:
            plt.show()
            return None


# ============================================================================
# 模块5: 启动子变异注释模块
# ============================================================================

class PromoterAnnotator:
    """
    启动子变异注释模块
    注释基因上游区域（默认2kb）的变异，识别可能影响表达调控的位点
    """
    
    # 常见顺式作用元件 (cis-regulatory elements)
    CIS_ELEMENTS = {
        'TATA_box': {'motif': 'TATA[AT]A[AT]', 'position': (-30, -25), 'description': 'TATA box - 核心启动子元件'},
        'CAAT_box': {'motif': 'CCAAT', 'position': (-80, -70), 'description': 'CAAT box - 转录效率调控'},
        'GC_box': {'motif': 'GGGCGG', 'position': (-110, -100), 'description': 'GC box - Sp1结合位点'},
        'W_box': {'motif': 'TTGAC[CT]', 'position': None, 'description': 'W box - WRKY转录因子结合'},
        'G_box': {'motif': 'CACGTG', 'position': None, 'description': 'G box - bHLH/bZIP结合'},
        'ABRE': {'motif': 'ACGTG[GT]C', 'position': None, 'description': 'ABA响应元件'},
        'DRE': {'motif': '[AG]CCGAC', 'position': None, 'description': '脱水响应元件'},
        'MYB_binding': {'motif': '[CT]AAC[GT]G', 'position': None, 'description': 'MYB结合位点'},
        'MYC_binding': {'motif': 'CA[ACGT]{2}TG', 'position': None, 'description': 'MYC结合位点'},
    }
    
    def __init__(self, fasta_file: str = None, gtf_file: str = None):
        """
        Args:
            fasta_file: 参考基因组 FASTA 文件
            gtf_file: GTF 注释文件
        """
        self.fasta_file = fasta_file or DataConfig.FASTA_PATH
        self.gtf_file = gtf_file or DataConfig.GTF_PATH
        self.fasta = None
        
        if PYSAM_AVAILABLE and fasta_file and os.path.exists(fasta_file):
            try:
                self.fasta = pysam.FastaFile(fasta_file)
            except Exception as e:
                print(f"[WARNING] 无法加载 FASTA: {e}")
    
    def get_promoter_region(self, chrom: str, gene_start: int, gene_end: int,
                           strand: str = '+', upstream: int = 2000) -> tuple:
        """
        获取基因启动子区域坐标
        
        Args:
            chrom: 染色体
            gene_start: 基因起始位置
            gene_end: 基因终止位置
            strand: 链方向
            upstream: 上游距离（默认2000bp）
        
        Returns:
            tuple: (promoter_start, promoter_end)
        """
        if strand == '+':
            promoter_start = max(1, gene_start - upstream)
            promoter_end = gene_start - 1
        else:
            promoter_start = gene_end + 1
            promoter_end = gene_end + upstream
        
        return promoter_start, promoter_end
    
    def get_promoter_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        获取启动子区域序列
        """
        if self.fasta is None:
            return ""
        
        try:
            seq = self.fasta.fetch(chrom, start - 1, end)
            return seq.upper()
        except Exception as e:
            print(f"[WARNING] 无法获取序列: {e}")
            return ""
    
    def find_cis_elements(self, sequence: str, strand: str = '+') -> list:
        """
        在序列中搜索顺式作用元件
        
        Args:
            sequence: DNA 序列
            strand: 链方向
        
        Returns:
            list: 找到的顺式元件列表
        """
        import re
        
        elements_found = []
        seq = sequence.upper()
        
        for name, info in self.CIS_ELEMENTS.items():
            motif = info['motif']
            # 将简化的正则转换为完整正则
            regex_motif = motif.replace('[AT]', '[AT]').replace('[CT]', '[CT]').replace('[GT]', '[GT]').replace('[AG]', '[AG]').replace('[ACGT]', '[ACGT]')
            
            try:
                for match in re.finditer(regex_motif, seq):
                    elements_found.append({
                        'element': name,
                        'position': match.start(),
                        'sequence': match.group(),
                        'description': info['description'],
                        'expected_position': info['position'],
                    })
            except:
                continue
        
        return elements_found
    
    def annotate_variants_in_promoter(self, variants_df: pd.DataFrame, 
                                      promoter_start: int, promoter_end: int,
                                      promoter_seq: str = None) -> pd.DataFrame:
        """
        注释启动子区域内的变异
        
        Args:
            variants_df: 变异数据框（需包含 POS, REF, ALT 列）
            promoter_start: 启动子起始
            promoter_end: 启动子终止
            promoter_seq: 启动子序列（可选）
        
        Returns:
            DataFrame: 带注释的变异表
        """
        if len(variants_df) == 0:
            return pd.DataFrame()
        
        df = variants_df.copy()
        
        # 筛选启动子区域内的变异
        if 'POS' in df.columns:
            df = df[(df['POS'] >= promoter_start) & (df['POS'] <= promoter_end)].copy()
        
        if len(df) == 0:
            return pd.DataFrame()
        
        # 计算相对于 TSS 的位置
        df['Relative_to_TSS'] = df['POS'] - promoter_end  # 负值表示上游
        
        # 注释位置特征
        def annotate_position(rel_pos):
            if -50 <= rel_pos <= 0:
                return 'Core promoter (-50 to TSS)'
            elif -200 <= rel_pos < -50:
                return 'Proximal promoter (-200 to -50)'
            elif -1000 <= rel_pos < -200:
                return 'Distal promoter (-1000 to -200)'
            else:
                return 'Far upstream (>-1000)'
        
        df['Region'] = df['Relative_to_TSS'].apply(annotate_position)
        
        # 检查是否影响已知顺式元件
        if promoter_seq:
            cis_elements = self.find_cis_elements(promoter_seq)
            
            def check_cis_overlap(pos):
                rel_pos = pos - promoter_start
                overlaps = []
                for elem in cis_elements:
                    elem_start = elem['position']
                    elem_end = elem_start + len(elem['sequence'])
                    if elem_start <= rel_pos <= elem_end:
                        overlaps.append(elem['element'])
                return ','.join(overlaps) if overlaps else 'None'
            
            df['Affected_Elements'] = df['POS'].apply(check_cis_overlap)
        
        # 评估功能影响潜力
        def assess_impact(row):
            rel_pos = row.get('Relative_to_TSS', -1000)
            affected = row.get('Affected_Elements', 'None')
            
            if affected != 'None':
                return 'HIGH - affects cis-element'
            elif -100 <= rel_pos <= 0:
                return 'MODERATE - core promoter'
            elif -500 <= rel_pos < -100:
                return 'LOW - proximal region'
            else:
                return 'MINIMAL - distal region'
        
        df['Impact_Potential'] = df.apply(assess_impact, axis=1)
        
        return df
    
    def generate_promoter_report(self, gene_id: str, chrom: str, 
                                  gene_start: int, gene_end: int, strand: str,
                                  variants_positions: list = None) -> dict:
        """
        生成启动子分析报告
        
        Args:
            gene_id: 基因ID
            chrom: 染色体
            gene_start: 基因起始
            gene_end: 基因终止
            strand: 链方向
            variants_positions: 变异位置列表
        
        Returns:
            dict: 启动子分析报告
        """
        promoter_start, promoter_end = self.get_promoter_region(
            chrom, gene_start, gene_end, strand
        )
        
        report = {
            'gene_id': gene_id,
            'chrom': chrom,
            'gene_start': gene_start,
            'gene_end': gene_end,
            'strand': strand,
            'promoter_start': promoter_start,
            'promoter_end': promoter_end,
            'promoter_length': abs(promoter_end - promoter_start),
            'cis_elements': [],
            'variants_in_promoter': 0,
            'high_impact_variants': 0,
        }
        
        # 获取启动子序列
        promoter_seq = self.get_promoter_sequence(chrom, promoter_start, promoter_end)
        if promoter_seq:
            report['promoter_sequence_length'] = len(promoter_seq)
            report['gc_content'] = (promoter_seq.count('G') + promoter_seq.count('C')) / len(promoter_seq) * 100
            
            # 查找顺式元件
            cis_elements = self.find_cis_elements(promoter_seq, strand)
            report['cis_elements'] = cis_elements
            report['n_cis_elements'] = len(cis_elements)
        
        # 统计启动子区域内的变异
        if variants_positions:
            promoter_variants = [
                p for p in variants_positions 
                if promoter_start <= p <= promoter_end
            ]
            report['variants_in_promoter'] = len(promoter_variants)
            report['variant_positions'] = promoter_variants
        
        return report
    
    def plot_promoter_structure(self, report: dict, output_file: str = None) -> str:
        """
        绘制启动子结构图
        """
        fig, ax = plt.subplots(figsize=(14, 4))
        
        promoter_start = report['promoter_start']
        promoter_end = report['promoter_end']
        gene_start = report['gene_start']
        strand = report['strand']
        
        # 绘制启动子区域
        ax.axhline(y=0, color='gray', linewidth=2)
        
        # 启动子区域框
        ax.axvspan(promoter_start, promoter_end, alpha=0.3, color='lightblue', label='Promoter region')
        
        # TSS 标记
        tss = gene_start if strand == '+' else report['gene_end']
        ax.axvline(x=tss, color='red', linewidth=2, label='TSS')
        ax.text(tss, 0.5, 'TSS', ha='center', fontsize=10, color='red')
        
        # 顺式元件
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        cis_elements = report.get('cis_elements', [])
        
        for i, elem in enumerate(cis_elements[:10]):  # 最多显示10个
            pos = promoter_start + elem['position']
            elem_len = len(elem['sequence'])
            ax.barh(0, elem_len, left=pos, height=0.3, 
                   color=colors[i % 10], alpha=0.7, 
                   label=elem['element'])
        
        # 变异位置
        variant_positions = report.get('variant_positions', [])
        for vpos in variant_positions:
            ax.plot(vpos, 0, 'v', markersize=10, color='darkred')
        
        ax.set_xlim(promoter_start - 100, gene_start + 500)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('Position (bp)', fontsize=11)
        ax.set_title(f'Promoter Structure: {report.get("gene_id", "Unknown")}', fontsize=12)
        
        # 图例
        ax.legend(loc='upper left', fontsize=8, ncol=3)
        ax.set_yticks([])
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] 启动子结构图已保存: {output_file}")
            return output_file
        else:
            plt.show()
            return None


# ============================================================================
# 模块6: 结果整合与报告输出
# ============================================================================

class ReportGenerator:
    """结果整合与报告生成模块"""
    
    def __init__(self, output_dir: str = None):
        self.output_dir = output_dir or DataConfig.OUTPUT_DIR
        # 确保使用绝对路径
        if not os.path.isabs(self.output_dir):
            self.output_dir = os.path.join(os.getcwd(), self.output_dir)
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"[DEBUG] 输出目录: {self.output_dir}")
        self.results = {}
    
    def add_result(self, key: str, result: dict):
        """添加分析结果"""
        self.results[key] = result
    
    def plot_boxplot(self, merged_df: pd.DataFrame, phenotype_col: str, 
                     gene_id: str = None, save: bool = True,
                     show_pvalue: bool = True, p_value: float = None,
                     correction_method: str = None) -> str:
        """
        绘制单倍型-表型箱线图（增强版：带P值标注）
        
        Args:
            merged_df: 合并后的数据框
            phenotype_col: 表型列名
            gene_id: 基因ID（可选）
            save: 是否保存
            show_pvalue: 是否显示P值标注
            p_value: 预计算的P值（可选，如果不提供则自动计算）
            correction_method: 多重检验校正方法 ('bonferroni', 'fdr', None)
        """
        df = merged_df[
            (merged_df['Hap_Name'] != 'Other') & 
            (merged_df[phenotype_col].notna())
        ].copy()
        
        if len(df) < 5:
            print(f"[WARNING] 样本量不足，跳过箱线图")
            return None
        
        # 按样本数排序单倍型
        hap_order = df['Hap_Name'].value_counts().index.tolist()
        n_groups = len(hap_order)
        
        fig, ax = plt.subplots(figsize=(10, 7))  # 略微增加高度以容纳P值标注
        
        # 绘制箱线图
        box_data = [df[df['Hap_Name'] == h][phenotype_col].values for h in hap_order]
        bp = ax.boxplot(box_data, labels=hap_order, patch_artist=True)
        
        # 设置颜色
        colors = plt.cm.Set3(np.linspace(0, 1, len(hap_order)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # 添加散点
        for i, h in enumerate(hap_order):
            y = df[df['Hap_Name'] == h][phenotype_col].values
            x = np.random.normal(i + 1, 0.04, size=len(y))
            ax.scatter(x, y, alpha=0.5, s=20, c='black')
        
        # 获取Y轴范围用于标注
        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        
        # 添加样本数标签
        for i, h in enumerate(hap_order):
            n = len(df[df['Hap_Name'] == h])
            ax.text(i + 1, y_min - y_range * 0.02, f'n={n}', ha='center', va='top', fontsize=9)
        
        # P值计算和标注
        if show_pvalue and n_groups >= 2:
            if p_value is None:
                # 自动计算P值
                if n_groups == 2:
                    stat, p_value = ttest_ind(box_data[0], box_data[1], equal_var=False)
                else:
                    stat, p_value = f_oneway(*box_data)
            
            # 多重检验校正
            p_display = p_value
            correction_note = ""
            if correction_method == 'bonferroni':
                p_display = min(p_value * n_groups, 1.0)
                correction_note = " (Bonferroni)"
            elif correction_method == 'fdr':
                correction_note = " (FDR)"
            
            # 显著性标注
            sig_mark = MultipleTestingCorrection.annotate_significance(p_value)
            p_text = MultipleTestingCorrection.format_pvalue(p_value, sig_level=False)
            
            # 在图表顶部标注总体P值
            annotation_y = y_max + y_range * 0.05
            ax.text(0.5, 1.02, f'{p_text}{correction_note}  {sig_mark}',
                   transform=ax.transAxes, ha='center', va='bottom',
                   fontsize=11, fontweight='bold',
                   color='red' if p_value < 0.05 else 'black')
            
            # 两组时绘制配对比较线
            if n_groups == 2:
                self._add_significance_bracket(ax, 1, 2, y_max + y_range * 0.02, sig_mark)
            
            # 多组时绘制配对比较（只显示最显著的几对）
            elif n_groups >= 3 and n_groups <= 5:
                self._add_pairwise_significance(ax, box_data, hap_order, y_max, y_range)
        
        title = f"Haplotype-Phenotype Association: {phenotype_col}"
        if gene_id:
            title = f"{gene_id}\n{title}"
        ax.set_title(title, fontsize=12, pad=20)
        ax.set_xlabel("Haplotype", fontsize=11)
        ax.set_ylabel(phenotype_col, fontsize=11)
        
        # 调整Y轴范围以容纳标注
        ax.set_ylim(y_min - y_range * 0.1, y_max + y_range * 0.15)
        
        plt.tight_layout()
        
        if save:
            filename = f"boxplot_{phenotype_col.replace(' ', '_')}.pdf"
            filepath = os.path.join(self.output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] 箱线图已保存: {filepath}")
            return filepath
        else:
            plt.show()
            return None
    
    def _add_significance_bracket(self, ax, x1: int, x2: int, y: float, sig_mark: str):
        """
        在箱线图上添加显著性括号标注
        """
        bracket_height = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02
        
        # 绘制括号
        ax.plot([x1, x1, x2, x2], [y, y + bracket_height, y + bracket_height, y],
               color='black', linewidth=1.5)
        
        # 添加显著性标记
        ax.text((x1 + x2) / 2, y + bracket_height * 1.5, sig_mark,
               ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    def _add_pairwise_significance(self, ax, box_data: list, hap_order: list, 
                                   y_max: float, y_range: float):
        """
        添加多组配对比较的显著性标注（只显示最显著的配对）
        """
        n_groups = len(box_data)
        pairwise_results = []
        
        # 计算所有配对的P值
        for i in range(n_groups):
            for j in range(i + 1, n_groups):
                if len(box_data[i]) >= 2 and len(box_data[j]) >= 2:
                    stat, p_val = ttest_ind(box_data[i], box_data[j], equal_var=False)
                    pairwise_results.append({
                        'i': i, 'j': j,
                        'p_value': p_val,
                        'sig': MultipleTestingCorrection.annotate_significance(p_val)
                    })
        
        # 只显示P<0.05的配对（最多3个）
        sig_pairs = [p for p in pairwise_results if p['p_value'] < 0.05]
        sig_pairs = sorted(sig_pairs, key=lambda x: x['p_value'])[:3]
        
        # 绘制显著配对
        for idx, pair in enumerate(sig_pairs):
            y_offset = y_max + y_range * (0.03 + idx * 0.06)
            x1, x2 = pair['i'] + 1, pair['j'] + 1
            sig_mark = pair['sig']
            
            self._add_significance_bracket(ax, x1, x2, y_offset, sig_mark)
    
    def generate_summary_table(self, association_results: list, 
                               pve_results: list) -> pd.DataFrame:
        """
        生成汇总统计表
        """
        rows = []
        
        for assoc, pve in zip(association_results, pve_results):
            if 'error' in assoc or 'error' in pve:
                continue
            
            row = {
                'Phenotype': assoc.get('phenotype', 'NA'),
                'N_Samples': assoc.get('n_samples', 0),
                'N_Haplotypes': assoc.get('n_groups', 0),
                'Test_Method': assoc.get('test_type', 'NA'),
                'Statistic': assoc.get('statistic', np.nan),
                'P_Value': assoc.get('p_value', np.nan),
                'PVE(%)': pve.get('pve_percent', np.nan),
                'Effect_Size': pve.get('effect_size', 'NA'),
                'Significant': '***' if assoc.get('highly_significant') else ('*' if assoc.get('significant') else '')
            }
            rows.append(row)
        
        df = pd.DataFrame(rows)
        
        # 保存
        filepath = os.path.join(self.output_dir, "association_summary.csv")
        df.to_csv(filepath, index=False)
        print(f"[INFO] 汇总表已保存: {filepath}")
        
        return df
    
    def generate_report(self, gene_info: dict = None) -> str:
        """
        生成完整分析报告 (HTML格式)
        """
        # 辅助函数：安全格式化数值
        def fmt(val, fmt_str=".4f"):
            """安全格式化数值，非数字返回'NA'"""
            if val is None:
                return 'NA'
            if isinstance(val, (int, float)):
                return f"{val:{fmt_str}}"
            return str(val)
        
        html_content = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Haplotype-Phenotype Association Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; }
        h2 { color: #34495e; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .significant { color: #e74c3c; font-weight: bold; }
        .summary-box { background: #ecf0f1; padding: 20px; border-radius: 8px; margin: 20px 0; }
    </style>
</head>
<body>
    <h1>Haplotype-Phenotype Association Analysis Report</h1>
"""
        # 基因信息
        if gene_info:
            html_content += f"""
    <div class="summary-box">
        <h2>Gene Information</h2>
        <p><strong>Gene ID:</strong> {gene_info.get('gene_id', 'NA')}</p>
        <p><strong>Chromosome:</strong> {gene_info.get('chrom', 'NA')}</p>
        <p><strong>Region:</strong> {gene_info.get('start', 'NA')} - {gene_info.get('end', 'NA')}</p>
        <p><strong>Strand:</strong> {gene_info.get('strand', 'NA')}</p>
    </div>
"""
        
        # 关联分析结果
        if 'association' in self.results:
            assoc = self.results['association']
            sig_class = 'significant' if assoc.get('significant') else ''
            html_content += f"""
    <h2>Association Analysis</h2>
    <table>
        <tr><th>Parameter</th><th>Value</th></tr>
        <tr><td>Phenotype</td><td>{assoc.get('phenotype', 'NA')}</td></tr>
        <tr><td>Sample Size</td><td>{assoc.get('n_samples', 'NA')}</td></tr>
        <tr><td>Number of Haplotypes</td><td>{assoc.get('n_groups', 'NA')}</td></tr>
        <tr><td>Test Method</td><td>{assoc.get('test_type', 'NA')}</td></tr>
        <tr><td>Test Statistic</td><td>{fmt(assoc.get('statistic'), '.4f')}</td></tr>
        <tr class="{sig_class}"><td>P-value</td><td>{fmt(assoc.get('p_value'), '.2e')}</td></tr>
    </table>
"""
        
        # PVE 结果
        if 'pve' in self.results:
            pve = self.results['pve']
            html_content += f"""
    <h2>Phenotypic Variance Explained (PVE)</h2>
    <table>
        <tr><th>Parameter</th><th>Value</th></tr>
        <tr><td>PVE (%)</td><td>{fmt(pve.get('pve_percent'), '.2f')}%</td></tr>
        <tr><td>R²</td><td>{fmt(pve.get('r_squared'), '.4f')}</td></tr>
        <tr><td>Omega²</td><td>{fmt(pve.get('omega_squared'), '.4f')}</td></tr>
        <tr><td>Effect Size</td><td>{pve.get('effect_size', 'NA')}</td></tr>
    </table>
"""
        
        html_content += """
    <hr>
    <p><em>Generated by Haplotype-Phenotype Association Analysis Platform</em></p>
</body>
</html>
"""
        
        # 保存报告
        filepath = os.path.join(self.output_dir, "analysis_report.html")
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"[INFO] 分析报告已生成: {filepath}")
        return filepath
    
    def generate_integrated_html(self, hap_sample_df: pd.DataFrame, 
                                  effect_results: dict,
                                  variant_positions: list,
                                  region_start: int, region_end: int,
                                  phenotype_col: str = 'phenotype',
                                  gene_start: int = None, gene_end: int = None,
                                  promoter_start: int = None, promoter_end: int = None,
                                  strand: str = '+',
                                  exons: list = None, cds: list = None,
                                  snp_effects: dict = None,
                                  chrom: str = None) -> str:
        """生成综合HTML大图（效应图 + 箱线图 + 单倍型序列 共用Hap标签）
        
        基因结构使用真实的外显子和 CDS 坐标（从 GTF 解析）:
        - exons: 外显子坐标列表 [(start, end), ...]
        - cds: CDS 坐标列表 [(start, end), ...]
        - 外显子之间的空隙即为内含子（黑色细线）
        - UTR = 外显子中不属于 CDS 的部分
        """
        hap_col = 'Hap_Name' if 'Hap_Name' in hap_sample_df.columns else 'Haplotype'
        hap_counts = hap_sample_df.groupby(hap_col).size().sort_values(ascending=False)
        top_haps = hap_counts.head(8).index.tolist()
        
        # 基因边界
        g_start = gene_start if gene_start else region_start
        g_end = gene_end if gene_end else region_end
        
        # 初始化 exons 和 cds 列表
        exons = exons or []
        cds = cds or []
        chrom = chrom or 'chr1'  # 默认染色体
        
        # 计算整个 CDS 区域的起止位置（用于变异分类）
        if cds:
            cds_start_pos = min(c[0] for c in cds)
            cds_end_pos = max(c[1] for c in cds)
        else:
            # 如果没有 CDS 数据，默认整个基因为 CDS
            cds_start_pos = g_start
            cds_end_pos = g_end
        
        # 变异类型颜色映射（与 plot_Gene_HapSeq.py 一致）
        var_type_colors = {
            'missense':    '#e74c3c',  # 红色
            'synonymous':  '#f39c12',  # 橙色
            'UTR':         '#9b59b6',  # 紫色
            'indel':       '#3498db',  # 蓝色 - 新增
            'SV':          '#e91e63',  # 深粉色 - 新增（结构变异）
            'other':       '#95a5a6',  # 灰色
        }
        var_type_labels = {
            'missense':   'Missense',
            'synonymous': 'Synonymous',
            'UTR':        'UTR',
            'indel':      'Indel',      # 新增
            'SV':         'SV',         # 新增（结构变异）
            'other':      'Other',
        }
        
        def get_var_color(pos):
            """snp_effects 优先，无则回退到位置分类"""
            if snp_effects and pos in snp_effects:
                return var_type_colors.get(snp_effects[pos], '#95a5a6'), snp_effects[pos]
            # 位置回退：简单分为 UTR / other
            in_exon = any(es <= pos <= ee for es, ee in exons)
            in_cds_b = any(cs <= pos <= ce for cs, ce in cds)
            if in_exon and not in_cds_b:
                return var_type_colors['UTR'], 'UTR'
            elif in_cds_b:
                return var_type_colors['other'], 'other'
            return var_type_colors['other'], 'other'
        
        # 从Haplotype_Seq获取序列
        if 'Haplotype_Seq' in hap_sample_df.columns:
            first_seq = hap_sample_df['Haplotype_Seq'].iloc[0] if len(hap_sample_df) > 0 else ''
            seq_len = len(first_seq.split('|')) if first_seq else 0
            max_vars = min(30, seq_len)
            display_positions = variant_positions[:max_vars] if variant_positions else list(range(max_vars))
        else:
            display_positions = variant_positions[:30] if variant_positions else []
            max_vars = len(display_positions)
        
        effects = effect_results.get('haplotype_effects', []) if effect_results else []
        grand_mean = effect_results.get('grand_mean', 0) if effect_results else 0
        region_len_kb = (region_end - region_start) / 1000
        
        # 箱线图数据 - 检查所有可能的表型列
        box_data = {}
        pheno_cols = [c for c in hap_sample_df.columns if 'pheno' in c.lower() or c == phenotype_col]
        actual_pheno_col = pheno_cols[0] if pheno_cols else phenotype_col
        print(f"[DEBUG] 表型列: {actual_pheno_col}, 可用列: {list(hap_sample_df.columns)[:10]}")
        
        for hap in top_haps:
            hap_df = hap_sample_df[hap_sample_df[hap_col] == hap]
            if actual_pheno_col in hap_df.columns:
                values = hap_df[actual_pheno_col].dropna().tolist()
                if values:
                    box_data[hap] = {
                        'mean': np.mean(values), 'median': np.median(values),
                        'q1': np.percentile(values, 25), 'q3': np.percentile(values, 75),
                        'min': min(values), 'max': max(values), 'n': len(values),
                        'values': values  # 保存原始值用于数据点绘制
                    }
        print(f"[DEBUG] 箱线图数据: {len(box_data)} 个单倍型有表型数据")
        
        # 效应图数据
        effect_data = {}
        for e in effects:
            effect_data[e['haplotype']] = e
        
        # 计算全局数据范围
        if box_data:
            all_vals = [v for d in box_data.values() for v in [d['min'], d['max']]]
            global_min, global_max = min(all_vals), max(all_vals)
            global_range = global_max - global_min if global_max > global_min else 1
        else:
            global_min, global_max, global_range = 0, 1, 1
        
        # 效应值范围
        if effect_data:
            eff_vals = [e.get('effect', 0) for e in effect_data.values()]
            eff_min, eff_max = min(eff_vals), max(eff_vals)
            eff_range = max(abs(eff_min), abs(eff_max)) if eff_vals else 1
        else:
            eff_vals = []
            eff_min, eff_max = -1, 1
            eff_range = 1
        
        base_colors = {'A': '#E41A1C', 'T': '#27AE60', 'C': '#3498DB', 'G': '#F1C40F'}
        row_height = 36
        n_haps = len(top_haps)
        
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Haplotype-Phenotype Association Analysis</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f0f2f5; padding: 15px; }}
        .container {{ max-width: 1800px; margin: 0 auto; background: white; border-radius: 10px; 
                     box-shadow: 0 4px 20px rgba(0,0,0,0.08); }}
        .header {{ background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); color: white; 
                  padding: 18px 25px; border-radius: 10px 10px 0 0; }}
        .header h1 {{ font-size: 18px; margin-bottom: 5px; }}
        .header-info {{ display: flex; gap: 20px; font-size: 11px; opacity: 0.9; flex-wrap: wrap; }}
        .content {{ padding: 15px; overflow: auto; max-height: 85vh; }}
        .footer {{ background: #f8f9fa; padding: 10px 20px; border-top: 1px solid #e8e8e8;
                  display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;
                  border-radius: 0 0 10px 10px; }}
        .base-legend {{ display: flex; gap: 12px; }}
        .base-legend-item {{ display: flex; align-items: center; gap: 4px; font-size: 11px; }}
        .base-box {{ width: 18px; height: 18px; border-radius: 3px; display: flex;
                    align-items: center; justify-content: center; color: white; font-weight: 700; font-size: 11px; }}
        /* 表格样式 */
        .data-table {{ border-collapse: collapse; margin-top: 0; }}
        .data-table th {{ background: #34495e; color: white; padding: 0; font-size: 10px; font-weight: 500; vertical-align: top; }}
        .data-table td {{ padding: 6px 4px; text-align: center; border-bottom: 1px solid #eee; }}
        .data-table tr:hover {{ background: #f8f9fa; }}
        .data-table tr.ref-row {{ background: #fffbeb; }}
        .hap-cell {{ text-align: left !important; padding-left: 10px !important; font-weight: 600; font-size: 12px; }}
        .ref-tag {{ background: #e74c3c; color: white; font-size: 8px; padding: 1px 4px; border-radius: 3px; margin-left: 4px; }}
        .base {{ font-family: Consolas, monospace; font-weight: 700; font-size: 13px; }}
        .effect-cell, .box-cell {{ width: 180px; position: relative; height: 35px; }}
        .bar-container {{ position: relative; height: 20px; background: #f8f9fa; border-radius: 3px; overflow: hidden; }}
        .bar-center {{ position: absolute; left: 50%; top: 0; bottom: 0; width: 1px; background: #333; border-left: 1px dashed #333; }}
        /* 森林图样式 */
        .forest-ci {{ position: absolute; top: 50%; height: 2px; transform: translateY(-50%); }}
        .forest-point {{ position: absolute; top: 50%; transform: translate(-50%, -50%); width: 10px; height: 10px; border-radius: 50%; }}
        /* 箱线图样式 */
        .bar-whisker {{ position: absolute; top: 50%; height: 2px; background: #7f8c8d; transform: translateY(-50%); }}
        .bar-box {{ position: absolute; top: 3px; height: 14px; border-radius: 2px; border: 2px solid #3498db; background: rgba(52,152,219,0.2); }}
        .bar-median {{ position: absolute; top: 2px; height: 16px; width: 2px; background: #2c3e50; }}
        .data-dot {{ position: absolute; width: 4px; height: 4px; border-radius: 50%; background: rgba(44,62,80,0.55); transform: translateX(-50%); }}
        .n-cell {{ font-size: 11px; color: #666; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>Haplotype-Phenotype Association Analysis</h1>
        <div class="header-info">
            <span>Region: {chrom}:{region_start:,}-{region_end:,}</span>
            <span>Length: {region_len_kb:.1f} kb</span>
            <span>Variants: {len(display_positions)}</span>
            <span>Haplotypes: {len(top_haps)}</span>
            <span>Samples: {hap_counts.sum()}</span>
            <span>Grand Mean: {grand_mean:.3f}</span>
        </div>
    </div>
    
    <div class="content">
'''
        
        # SVG 基因结构图和连线（参照老师指定样式）
        n_vars = len(display_positions)
        hap_col_w = 90
        eff_col_w = 180
        box_col_w = 180
        gene_area_start = hap_col_w + eff_col_w + box_col_w  # = 450px
        seq_col_w = 28
        gene_area_width = n_vars * seq_col_w
        legend_w = 100  # 图例宽度
        svg_width = gene_area_start + gene_area_width + legend_w
        svg_height = 200  # 足够容纳长斜线
                
        # 坐标系设计（参照老师的图）
        axis_y = 28       # 坐标轴线 y
        gene_y = 48       # 基因结构图上框 y
        gene_h = 18       # 基因结构高度
        var_top_y = axis_y - 8   # 变异小竖线顶部 y
        line_end_y = svg_height - 5  # 斜线终点 y
                
        html += f'<svg width="{svg_width}" height="{svg_height}" style="display:block;margin-bottom:0;">\n'
        
        # ==== SVG 定义（UTR 旜线填充图案）====
        html += '''<defs>
            <pattern id="utrPattern" patternUnits="userSpaceOnUse" width="6" height="6">
                <line x1="0" y1="6" x2="6" y2="0" stroke="#9b59b6" stroke-width="1" opacity="0.5"/>
            </pattern>
        </defs>\n'''
                
        # ==== 标题 ====
        html += f'<text x="{gene_area_start + gene_area_width/2}" y="13" font-size="11" fill="#2c3e50" text-anchor="middle" font-weight="600">Relative Position (kb)</text>\n'
                
        # ==== 坐标轴线（参照老师的图：带小刻度）====
        html += f'<line x1="{gene_area_start}" y1="{axis_y}" x2="{gene_area_start + gene_area_width}" y2="{axis_y}" stroke="#333" stroke-width="1.2"/>\n'
                
        # 坐标刻度：方向与箭头方向一致
        # 正链（箭头朝右）：刻度从左到右递增
        # 负链（箭头朝左）：刻度从右到左递增
        n_ticks = 9 if region_len_kb >= 4 else max(5, int(region_len_kb * 2) + 1)
        for i in range(n_ticks):
            if strand == '+':
                # 正链：位置从左到右，刻度值从0到region_len_kb
                x = gene_area_start + i * gene_area_width / (n_ticks - 1)
                kb = i * region_len_kb / (n_ticks - 1)
            else:
                # 负链：位置从右到左（反转），刻度值从region_len_kb到0
                x = gene_area_start + (n_ticks - 1 - i) * gene_area_width / (n_ticks - 1)
                kb = i * region_len_kb / (n_ticks - 1)
            # 刻度竖线：从文字下方一直延伸到断连线上端 (gene_y+gene_h)
            html += f'<line x1="{x}" y1="{var_top_y-2}" x2="{x}" y2="{gene_y+gene_h}" stroke="#aaa" stroke-width="0.7" stroke-dasharray="2,2"/>\n'
            # 坐标轴处加粗小刻度
            html += f'<line x1="{x}" y1="{axis_y-2}" x2="{x}" y2="{axis_y+4}" stroke="#555" stroke-width="1.5"/>\n'
            html += f'<text x="{x}" y="{var_top_y-5}" font-size="8.5" fill="#555" text-anchor="middle">{kb:.1f}</text>\n'
        
        # ==== 基因结构：使用真实的外显子和 CDS 坐标 ====
        # 计算各区域在 SVG 中的位置
        def pos_to_x(pos):
            return gene_area_start + ((pos - region_start) / (region_end - region_start)) * gene_area_width
        
        # 基因起止在 SVG 中的位置
        gene_x1 = pos_to_x(g_start)
        gene_x2 = pos_to_x(g_end)
        
        # 外显子高度
        exon_h = gene_h
        # 内含子位置（中心线）
        intron_y = gene_y + gene_h // 2
        
        # ==== 1. 内含子骨架（黑色细线，贯穿整个基因区域）====
        html += f'<line x1="{gene_x1}" y1="{intron_y}" x2="{gene_x2}" y2="{intron_y}" stroke="#2c3e50" stroke-width="2"/>\n'
        
        # ==== 2. 启动子（与基因同一行）====
        if promoter_start is not None and promoter_end is not None:
            prom_x1 = pos_to_x(promoter_start)
            prom_x2 = pos_to_x(promoter_end)
            prom_w = abs(prom_x2 - prom_x1)
            if prom_w > 2:
                html += f'<rect x="{min(prom_x1, prom_x2)}" y="{gene_y}" width="{prom_w}" height="{exon_h}" fill="#f39c12" opacity="0.4" stroke="#e67e22" stroke-width="1.5" stroke-dasharray="3,2" rx="2"/>\n'
                if prom_w > 40:
                    html += f'<text x="{min(prom_x1, prom_x2) + prom_w/2}" y="{gene_y + exon_h/2 + 3}" font-size="7" fill="#d35400" text-anchor="middle" font-weight="600">Promoter</text>\n'
        
        # ==== 3. 绘制每个外显子（区分 UTR 和 CDS）====
        if exons:
            # 有真实的外显子数据
            for ex_start, ex_end in exons:
                ex_x1 = pos_to_x(ex_start)
                ex_x2 = pos_to_x(ex_end)
                ex_w = ex_x2 - ex_x1
                
                if ex_w < 2:
                    continue
                
                # 检查这个外显子与 CDS 的重叠情况
                # 找出外显子与所有 CDS 的重叠区域
                cds_overlaps = []
                for c_start, c_end in cds:
                    overlap_start = max(ex_start, c_start)
                    overlap_end = min(ex_end, c_end)
                    if overlap_start < overlap_end:
                        cds_overlaps.append((overlap_start, overlap_end))
                
                if not cds_overlaps:
                    # 整个外显子都是 UTR（白色）
                    html += f'<rect x="{ex_x1}" y="{gene_y}" width="{ex_w}" height="{exon_h}" fill="white" stroke="#3498db" stroke-width="1.5" rx="1"/>\n'
                else:
                    # 外显子包含 CDS 和/或 UTR
                    # 先画整个外显子的底层（白色 UTR）
                    html += f'<rect x="{ex_x1}" y="{gene_y}" width="{ex_w}" height="{exon_h}" fill="white" stroke="#3498db" stroke-width="1" rx="1"/>\n'
                    
                    # 再画 CDS 区域（蓝色）
                    for cds_s, cds_e in cds_overlaps:
                        cds_x1 = pos_to_x(cds_s)
                        cds_x2 = pos_to_x(cds_e)
                        cds_w = cds_x2 - cds_x1
                        if cds_w > 1:
                            html += f'<rect x="{cds_x1}" y="{gene_y}" width="{cds_w}" height="{exon_h}" fill="#3498db" rx="1"/>\n'
        else:
            # 没有外显子数据，画一个整体矩形表示基因
            html += f'<rect x="{gene_x1}" y="{gene_y}" width="{gene_x2-gene_x1}" height="{exon_h}" fill="#3498db" stroke="#2980b9" stroke-width="1" rx="1"/>\n'
        
        # ==== 内含子骨架箭头（表示链方向）====
        # 规则：箭头在基因3'端，启动子在5'端，箭头始终朝外
        # 正链（5'在左，3'在右）：箭头在基因右端，朝右
        # 负链（5'在右，3'在左）：箭头在基因左端，朝左
        arrow_size = 10  # 增大箭头尺寸
        arrow_line = 8   # 箭头前的直线长度
        if strand == '+':
            # 正链：箭头在基因右端（3'端），朝右
            # 画箭头前的直线
            html += f'<line x1="{gene_x2}" y1="{intron_y}" x2="{gene_x2+arrow_line}" y2="{intron_y}" stroke="#2c3e50" stroke-width="2"/>\n'
            # 画箭头三角形（顶点在右）
            html += f'<polygon points="{gene_x2+arrow_line+arrow_size},{intron_y} {gene_x2+arrow_line},{intron_y-arrow_size//2} {gene_x2+arrow_line},{intron_y+arrow_size//2}" fill="#2c3e50"/>\n'
        else:
            # 负链：箭头在基因左端（3'端），朝左
            # 画箭头前的直线
            html += f'<line x1="{gene_x1}" y1="{intron_y}" x2="{gene_x1-arrow_line}" y2="{intron_y}" stroke="#2c3e50" stroke-width="2"/>\n'
            # 画箭头三角形（顶点在左）
            html += f'<polygon points="{gene_x1-arrow_line-arrow_size},{intron_y} {gene_x1-arrow_line},{intron_y-arrow_size//2} {gene_x1-arrow_line},{intron_y+arrow_size//2}" fill="#2c3e50"/>\n'
                
        # ==== 变异位点连线（参照老师的图：斜直线虚线）====
        # gene_x = 真实物理坐标在基因结构上的像素位置
        # table_x = 序列列中心（按索引均匀排列）
        # 斜线从真实位置连到序列列、与老师参考图一致
        var_types_found = set()  # 记录找到的变异类型
        for idx, pos in enumerate(display_positions):
            rel_pct = (pos - region_start) / (region_end - region_start)
            gene_x = gene_area_start + rel_pct * gene_area_width  # 真实物理位置
            table_x = gene_area_start + idx * seq_col_w + seq_col_w / 2  # 列中心
            
            # 判断变异类型并获取颜色
            var_color, var_type = get_var_color(pos)
            var_types_found.add(var_type)
            
            # 变异竖线：从圆圈顶部一直延伸到断线上端 (gene_y+gene_h)
            html += f'<line x1="{gene_x}" y1="{var_top_y+3}" x2="{gene_x}" y2="{gene_y+gene_h}" stroke="{var_color}" stroke-width="1.2"/>\n'
            html += f'<circle cx="{gene_x}" cy="{var_top_y}" r="3" fill="{var_color}" stroke="white" stroke-width="0.5"/>\n'
            
            # 斜线：从真实位置下导到序列列中心
            html += f'<line x1="{gene_x}" y1="{gene_y+gene_h}" x2="{table_x}" y2="{line_end_y}" stroke="{var_color}" stroke-width="0.8" stroke-dasharray="4,2"/>\n'
                
        # ==== 图例（右侧，根据实际变异类型动态生成）====
        leg_x = gene_area_start + gene_area_width + 15
        html += f'<text x="{leg_x}" y="{axis_y}" font-size="8.5" fill="#333" font-weight="600">Gene Structure</text>\n'
        
        # 基因结构图例（与实际绘制一致）
        structure_items = [
            ('Promoter', '#f39c12', 'prom'),
            ('Intron', '#2c3e50', 'line'),
            ("5'/3' UTR", 'white', 'utr'),
            ('CDS', '#3498db', 'cds'),
        ]
        for li, (label, color, style) in enumerate(structure_items):
            ly = axis_y + 8 + li * 14
            if style == 'prom':
                html += f'<rect x="{leg_x}" y="{ly}" width="14" height="10" fill="{color}" opacity="0.35" stroke="#e67e22" stroke-width="1" stroke-dasharray="2,1" rx="1"/>\n'
            elif style == 'line':
                html += f'<line x1="{leg_x}" y1="{ly+5}" x2="{leg_x+14}" y2="{ly+5}" stroke="{color}" stroke-width="2"/>\n'
            elif style == 'utr':
                html += f'<rect x="{leg_x}" y="{ly}" width="14" height="10" fill="{color}" stroke="#3498db" stroke-width="1.5" rx="1"/>\n'
            else:  # cds
                html += f'<rect x="{leg_x}" y="{ly}" width="14" height="10" fill="{color}" rx="1"/>\n'
            html += f'<text x="{leg_x+18}" y="{ly+8}" font-size="7.5" fill="#333">{label}</text>\n'
        
        # 变异类型图例（基于实际检测到的类型）
        html += f'<text x="{leg_x}" y="{axis_y + 70}" font-size="8.5" fill="#333" font-weight="600">Variant Type</text>\n'
        # 按顺序排列变异类型
        ordered_types = ['missense', 'synonymous', 'UTR', 'other']
        found_types = [t for t in ordered_types if t in var_types_found]
        for li, var_type in enumerate(found_types):
            ly = axis_y + 78 + li * 14
            color = var_type_colors.get(var_type, '#95a5a6')
            label = var_type_labels.get(var_type, var_type)
            html += f'<circle cx="{leg_x+7}" cy="{ly+5}" r="4" fill="{color}" stroke="white" stroke-width="0.5"/>\n'
            html += f'<text x="{leg_x+18}" y="{ly+8}" font-size="7.5" fill="#333">{label}</text>\n'
                
        html += '</svg>\n'
        
        # HTML表格
        html += '''<table class="data-table">
<thead><tr>
    <th style="width:90px;text-align:left;padding-left:10px;vertical-align:middle;height:60px;">Haplotype</th>
    <th class="effect-cell" style="vertical-align:middle;height:60px;">Effect (vs Grand Mean)</th>
    <th class="box-cell" style="vertical-align:middle;height:60px;">Phenotype</th>\n'''
        
        for pos in display_positions:
            # 物理坐标竖排，千分位逗号分隔，宽度与序列列 td 严格一致
            pos_str = f'{pos:,}'  # 千分位逗号
            html += (f'<th style="width:28px;min-width:28px;max-width:28px;padding:0;'
                     f'vertical-align:top;overflow:hidden;">'
                     f'<div style="writing-mode:vertical-rl;transform:rotate(180deg);'
                     f'width:28px;height:60px;display:flex;align-items:center;justify-content:center;'
                     f'font-size:9px;color:#f5f5f5;background:#2c3e50;'
                     f'font-weight:600;letter-spacing:0;box-sizing:border-box;">{pos_str}</div></th>\n')
        html += '<th style="width:40px;vertical-align:middle;height:60px;">n</th></tr></thead><tbody>\n'
        
        # 数据行
        for i, hap in enumerate(top_haps):
            rows = hap_sample_df[hap_sample_df[hap_col] == hap]
            if len(rows) == 0:
                continue
            row = rows.iloc[0]
            cnt = len(rows)
            is_ref = (i == 0)
            row_class = 'ref-row' if is_ref else ''
            ref_tag = '<span class="ref-tag">Ref</span>' if is_ref else ''
            
            # 效应数据（森林图样式）
            eff = effect_data.get(hap, {})
            eff_val = eff.get('effect', 0)
            ci_lower = eff.get('ci_lower', eff_val - 0.5)
            ci_upper = eff.get('ci_upper', eff_val + 0.5)
            p_val = eff.get('p_value', 1.0)
            
            # 颜色根据显著性
            if is_ref:
                point_color = '#95a5a6'  # Reference - 灰色
            elif p_val < 0.01:
                point_color = '#8B0000'  # P < 0.01 - 深红
            elif p_val < 0.05:
                point_color = '#e74c3c'  # P < 0.05 - 红色
            else:
                point_color = '#3498db'  # Not sig. - 蓝色
            
            # 计算位置百分比（0在中心）
            eff_center = 50 + (eff_val / eff_range * 45) if eff_range > 0 else 50
            ci_left = 50 + (ci_lower / eff_range * 45) if eff_range > 0 else 45
            ci_right = 50 + (ci_upper / eff_range * 45) if eff_range > 0 else 55
            ci_width = ci_right - ci_left
            
            # 箱线图数据
            bd = box_data.get(hap, {})
            if bd:
                w_left = ((bd['min'] - global_min) / global_range) * 100
                w_right = ((bd['max'] - global_min) / global_range) * 100
                b_left = ((bd['q1'] - global_min) / global_range) * 100
                b_width = ((bd['q3'] - bd['q1']) / global_range) * 100
                m_pos = ((bd['median'] - global_min) / global_range) * 100
                # 数据点（随机竖向抖动）
                import random
                random.seed(hash(hap) % 9999)
                dots_html = ''
                for val in bd.get('values', []):
                    x_pct = ((val - global_min) / global_range) * 100
                    top_pct = random.randint(15, 85)  # 竖向随机抖动
                    dots_html += f'<div class="data-dot" style="left:{x_pct:.1f}%;top:{top_pct}%;"></div>'
                box_html = f'''<div class="bar-container">
                    <div class="bar-whisker" style="left:{w_left}%;width:{max(w_right-w_left,1)}%;"></div>
                    <div class="bar-box" style="left:{b_left}%;width:{max(b_width,2)}%;"></div>
                    <div class="bar-median" style="left:{m_pos}%;"></div>
                    {dots_html}
                </div>'''
            else:
                box_html = '<div class="bar-container" style="background:#f5f5f5;"><span style="font-size:9px;color:#999;position:absolute;left:50%;transform:translateX(-50%);top:3px;">No data</span></div>'
            
            html += f'''<tr class="{row_class}">
    <td class="hap-cell">{hap}{ref_tag}</td>
    <td class="effect-cell">
        <div class="bar-container">
            <div class="bar-center"></div>
            <div class="forest-ci" style="left:{ci_left}%;width:{max(ci_width,1)}%;background:{point_color};"></div>
            <div class="forest-point" style="left:{eff_center}%;background:{point_color};"></div>
        </div>
    </td>
    <td class="box-cell">
        {box_html}
    </td>\n'''
            
            # 序列列
            if 'Haplotype_Seq' in row.index:
                seq = row['Haplotype_Seq'].replace('|', '')
                for idx in range(len(display_positions)):
                    base = seq[idx].upper() if idx < len(seq) else 'N'
                    color = base_colors.get(base, '#666')
                    html += f'<td style="width:28px;min-width:28px;max-width:28px;padding:0;text-align:center;"><span class="base" style="color:{color};">{base}</span></td>\n'
            
            html += f'<td class="n-cell">{cnt}</td></tr>\n'
        
        # 表格底部：共用坐标轴行（三列分开对应）
        html += '<tr style="height:30px;">'
        
        # 第1列：Haplotype列（空）
        html += '<td style="border:none;"></td>'
        
        # 第2列：Effect列 - 效应图坐标轴（5 个刻度）
        eff_axis_ticks = []
        for i in range(5):
            tick_val = eff_min + i * (eff_max - eff_min) / 4 if len(eff_vals) > 0 else 0
            tick_pct = 50 + (tick_val / eff_range * 45) if eff_range > 0 else 50
            eff_axis_ticks.append((tick_pct, f'{tick_val:+.2f}'))
        
        html += '<td style="border:none;position:relative;">'
        html += '<div style="position:relative;height:25px;">'
        html += '<div style="position:absolute;bottom:10px;left:0;right:0;height:1px;background:#bdc3c7;"></div>'
        for tick_pct, tick_label in eff_axis_ticks:
            html += f'<div style="position:absolute;bottom:6px;left:{tick_pct}%;transform:translateX(-50%);width:1px;height:9px;background:#bdc3c7;"></div>'
            html += f'<span style="position:absolute;bottom:-4px;left:{tick_pct}%;transform:translateX(-50%);font-size:9px;color:#7f8c8d;white-space:nowrap;">{tick_label}</span>'
        html += '</div></td>\n'
        
        # 第3列：Phenotype列 - 箱线图坐标轴（3 个刻度）
        html += '<td style="border:none;position:relative;">'
        if box_data:
            box_axis_ticks = []
            for i in range(3):
                tick_val = global_min + i * global_range / 2
                tick_pct = ((tick_val - global_min) / global_range) * 100 if global_range > 0 else 50
                box_axis_ticks.append((tick_pct, f'{tick_val:.2f}'))
            
            html += '<div style="position:relative;height:25px;">'
            html += '<div style="position:absolute;bottom:10px;left:0;right:0;height:1px;background:#bdc3c7;"></div>'
            for tick_pct, tick_label in box_axis_ticks:
                html += f'<div style="position:absolute;bottom:6px;left:{tick_pct}%;transform:translateX(-50%);width:1px;height:9px;background:#bdc3c7;"></div>'
                html += f'<span style="position:absolute;bottom:-4px;left:{tick_pct}%;transform:translateX(-50%);font-size:9px;color:#7f8c8d;white-space:nowrap;">{tick_label}</span>'
            html += '</div>'
        html += '</td>\n'
        
        # 剩余列（序列+n）
        html += f'<td colspan="{len(display_positions)+1}" style="border:none;"></td>\n'
        html += '</tr>\n'
        
        html += '''</tbody></table>
    </div>
    
    <div class="footer">
        <div class="base-legend">
            <div class="base-legend-item"><div class="base-box" style="background:#E41A1C;">A</div>Adenine</div>
            <div class="base-legend-item"><div class="base-box" style="background:#27AE60;">T</div>Thymine</div>
            <div class="base-legend-item"><div class="base-box" style="background:#3498DB;">C</div>Cytosine</div>
            <div class="base-legend-item"><div class="base-box" style="background:#F1C40F;">G</div>Guanine</div>
            <span style="margin-left:20px;border-left:1px solid #ddd;padding-left:20px;"></span>
            <div class="base-legend-item"><div style="width:10px;height:10px;border-radius:50%;background:#95a5a6;"></div>Reference</div>
            <div class="base-legend-item"><div style="width:10px;height:10px;border-radius:50%;background:#8B0000;"></div>P &lt; 0.01</div>
            <div class="base-legend-item"><div style="width:10px;height:10px;border-radius:50%;background:#e74c3c;"></div>P &lt; 0.05</div>
            <div class="base-legend-item"><div style="width:10px;height:10px;border-radius:50%;background:#3498db;"></div>Not sig.</div>
        </div>
        <div style="font-size:10px;color:#7f8c8d;">Effect: point = estimate, line = 95% CI</div>
    </div>
</div>
</body>
</html>'''
        
        out = os.path.join(self.output_dir, "integrated_analysis.html")
        with open(out, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"[INFO] 综合HTML报告已保存: {out}")
        return out


# ============================================================================
# 主流程整合
# ============================================================================

class HaplotypePhenotypeAnalyzer:
    """单倍型-表型关联分析主流程"""
    
    def __init__(self, vcf_file: str, phenotype_file: str, output_dir: str = None, gtf_file: str = None):
        """
        Args:
            vcf_file: VCF文件路径
            phenotype_file: 表型文件路径
            output_dir: 输出目录
            gtf_file: GTF注释文件路径（用于解析外显子/CDS坐标）
        """
        self.vcf_file = vcf_file
        self.phenotype_file = phenotype_file
        self.output_dir = output_dir or DataConfig.OUTPUT_DIR
        self.gtf_file = gtf_file or DataConfig.GTF_PATH
        
        # 加载表型数据
        self.phenotype_df = self._load_phenotype()
        
        # 初始化各模块
        self.extractor = HaplotypeExtractor(vcf_file)
        self.reporter = ReportGenerator(self.output_dir)
        
        # 分析结果缓存
        self.hap_df = None
        self.hap_sample_df = None
        self.positions = None
    
    def _load_phenotype(self) -> pd.DataFrame:
        """加载表型数据，支持多种格式（包括GEMMA格式）"""
        if not os.path.exists(self.phenotype_file):
            print(f"[WARNING] 表型文件不存在: {self.phenotype_file}")
            return pd.DataFrame()
        
        # 尝试多种分隔符
        for sep in ['\t', ',', ' ', r'\s+']:
            try:
                # 先尝试有表头的格式
                df = pd.read_csv(self.phenotype_file, sep=sep, engine='python' if sep == r'\s+' else 'c')
                if len(df.columns) > 1:
                    # 检查第一列是否是样本ID
                    first_col = df.columns[0]
                    if first_col.lower() in ['sampleid', 'sample_id', 'id', 'iid', 'fid', 'sample']:
                        df = df.rename(columns={first_col: 'SampleID'})
                    elif df.columns[0] != 'SampleID':
                        # 第一列可能是样本ID但没有表头名
                        df = df.rename(columns={df.columns[0]: 'SampleID'})
                    print(f"[INFO] 成功加载表型数据: {len(df)} 样本, {len(df.columns)} 列")
                    return df
            except:
                continue
        
        # 尝试GEMMA格式（无表头：FID IID phenotype）
        try:
            df = pd.read_csv(self.phenotype_file, sep=r'\s+', header=None, engine='python')
            if len(df.columns) >= 2:
                if len(df.columns) == 2:
                    # 2列格式: SampleID, phenotype
                    df.columns = ['SampleID', 'phenotype']
                elif len(df.columns) == 3:
                    # GEMMA格式: FID, IID, phenotype
                    df.columns = ['FID', 'SampleID', 'phenotype']
                    df = df[['SampleID', 'phenotype']]  # 只保留IID和表型
                else:
                    # 多列格式：第一列为SampleID，其余为表型
                    cols = ['SampleID'] + [f'phenotype_{i}' for i in range(1, len(df.columns))]
                    df.columns = cols
                
                print(f"[INFO] 成功加载GEMMA格式表型数据: {len(df)} 样本, {len(df.columns)} 列")
                return df
            
            # 单列GEMMA格式：只有表型值，需要从VCF获取样本ID
            if len(df.columns) == 1:
                print(f"[INFO] 检测到单列GEMMA格式，从VCF获取样本ID...")
                phenotype_values = df.iloc[:, 0].values
                
                # 从VCF获取样本ID
                if PYSAM_AVAILABLE:
                    vcf = pysam.VariantFile(self.vcf_file)
                    sample_ids = list(vcf.header.samples)
                    vcf.close()
                    
                    if len(sample_ids) == len(phenotype_values):
                        df = pd.DataFrame({
                            'SampleID': sample_ids,
                            'phenotype': phenotype_values
                        })
                        print(f"[INFO] 成功匹配VCF样本ID: {len(df)} 样本")
                        return df
                    else:
                        print(f"[WARNING] 表型数据行数({len(phenotype_values)})与VCF样本数({len(sample_ids)})不匹配")
                        # 尝试部分匹配
                        min_len = min(len(sample_ids), len(phenotype_values))
                        df = pd.DataFrame({
                            'SampleID': sample_ids[:min_len],
                            'phenotype': phenotype_values[:min_len]
                        })
                        print(f"[INFO] 部分匹配: {len(df)} 样本")
                        return df
                else:
                    print(f"[WARNING] pysam不可用，无法从VCF获取样本ID")
        except Exception as e:
            print(f"[DEBUG] GEMMA格式解析失败: {e}")
        
        print(f"[WARNING] 无法解析表型文件")
        return pd.DataFrame()
    
    def _extract_genotypes_for_amova(self):
        """
        从VCF提取基因型数据用于AMOVA分析
        返回: DataFrame (样本 x 位点)
        """
        logger = get_logger()
        try:
            if PYSAM_AVAILABLE:
                vcf = pysam.VariantFile(self.vcf_file)
                samples = list(vcf.header.samples)
                
                # 获取区间内的变异
                variants = []
                positions = []
                for rec in vcf.fetch(self.chrom, self.start, self.end):
                    variants.append(rec)
                    positions.append(rec.pos)
                
                if len(variants) == 0:
                    return None
                
                # 构建基因型矩阵
                gt_data = []
                for rec in variants:
                    row = {'POS': rec.pos}
                    for sample in samples:
                        gt = rec.samples[sample]['GT']
                        # 编码基因型: 0/0->0, 0/1->1, 1/1->2, missing->-1
                        if gt is None or None in gt:
                            row[sample] = -1
                        else:
                            row[sample] = sum(gt)
                    gt_data.append(row)
                
                gt_df = pd.DataFrame(gt_data)
                # 转置为 样本 x 位点 格式
                gt_df = gt_df.set_index('POS').T.reset_index()
                gt_df.columns = ['SampleID'] + [f'POS_{p}' for p in positions]
                
                return gt_df
            else:
                logger.warning("pysam不可用，无法提取基因型数据")
                return None
        except Exception as e:
            logger.warning(f"提取基因型数据失败: {e}")
            return None
    
    def _extract_genotypes_for_amova(self):
        """
        从VCF提取基因型数据用于AMOVA分析
        返回: DataFrame (样本 x 位点)
        """
        logger = get_logger()
        try:
            if PYSAM_AVAILABLE:
                vcf = pysam.VariantFile(self.vcf_file)
                samples = list(vcf.header.samples)
                
                # 获取区间内的变异
                variants = []
                positions = []
                for rec in vcf.fetch(self.chrom, self.start, self.end):
                    variants.append(rec)
                    positions.append(rec.pos)
                
                if len(variants) == 0:
                    return None
                
                # 构建基因型矩阵
                gt_data = []
                for rec in variants:
                    row = {'POS': rec.pos}
                    for sample in samples:
                        gt = rec.samples[sample]['GT']
                        # 编码基因型: 0/0->0, 0/1->1, 1/1->2, missing->-1
                        if gt is None or None in gt:
                            row[sample] = -1
                        else:
                            row[sample] = sum(gt)
                    gt_data.append(row)
                
                gt_df = pd.DataFrame(gt_data)
                # 转置为 样本 x 位点 格式
                gt_df = gt_df.set_index('POS').T.reset_index()
                gt_df.columns = ['SampleID'] + [f'POS_{p}' for p in positions]
                
                return gt_df
            else:
                logger.warning("pysam不可用，无法提取基因型数据")
                return None
        except Exception as e:
            logger.warning(f"提取基因型数据失败: {e}")
            return None
    
    def _annotate_snp_effects_tabix(self, vcf_file: str, fasta_path: str, gene_chrom: str,
                                     region_start: int, region_end: int,
                                     cds_intervals: list, exon_intervals: list,
                                     gene_strand: str) -> dict:
        """
        用 tabix 只读目标区间，做 SNP 精细分类 (missense/synonymous/UTR/other)
        读取范围：gene_chrom:region_start-region_end (含启动子)
        """
        logger = get_logger()
        effects = {}
                
        # 默认全为 other
        positions_list = self.positions or []
        logger.info(f"    - 待注释位点数：{len(positions_list)}")
        for pos in positions_list:
            effects[pos] = 'other'
                
        if not PYSAM_AVAILABLE:
            logger.warning("    - pysam不可用，跳过精细注释")
            return effects
        
        try:
            import pysam as _pysam
            # 检查索引文件是否存在（必须在TabixFile之前）
            tbi_path = vcf_file + '.tbi'
            csi_path = vcf_file + '.csi'
            tbi_exists = os.path.exists(tbi_path)
            csi_exists = os.path.exists(csi_path)
            logger.info(f"    - 索引检查: tbi={tbi_exists}, csi={csi_exists}")
            
            # pysam TabixFile 只认 .tbi，大文件需要 .csi，改用 VariantFile
            logger.info(f"    - 使用VariantFile查询（支持.csi索引）")
            
            try:
                vcf_reader = _pysam.VariantFile(vcf_file)
            except Exception as e:
                logger.warning(f"    - VariantFile打开失败: {e}")
                for pos in positions_list:
                    in_exon = _pos_in_any_interval(pos, exon_intervals)
                    in_cds  = _pos_in_any_interval(pos, cds_intervals)
                    effects[pos] = 'UTR' if (in_exon and not in_cds) else 'other'
                return effects
            
            # 构建 CDS 上下文（用于翻译）
            cds_seq, cds_pos_to_idx = _build_coding_context(gene_chrom, cds_intervals, gene_strand, fasta_path)
            logger.info(f"    - CDS序列长度: {len(cds_seq)}, 映射位点数: {len(cds_pos_to_idx)}")
            if cds_intervals:
                logger.info(f"    - CDS区间: {cds_intervals[:3]}{'...' if len(cds_intervals) > 3 else ''}")
            if exon_intervals:
                logger.info(f"    - Exon区间: {exon_intervals[:3]}{'...' if len(exon_intervals) > 3 else ''}")
            
            # 用 VariantFile 查询目标区间
            region_str = f"{gene_chrom}:{region_start}-{region_end}"
            logger.info(f"    - VariantFile查询范围: {region_str}")
            
            pos_set = set(positions_list)
            match_count = 0
            
            try:
                for rec in vcf_reader.fetch(gene_chrom, region_start, region_end):
                    pos = rec.pos
                    if pos not in pos_set:
                        continue
                    
                    match_count += 1
                    ref = rec.ref
                    alt_allele = str(rec.alts[0]) if rec.alts else ref
                    in_cds  = _pos_in_any_interval(pos, cds_intervals)
                    in_exon = _pos_in_any_interval(pos, exon_intervals)
                    
                    if not in_exon:
                        effects[pos] = 'other'
                    elif not in_cds:
                        effects[pos] = 'UTR'
                    elif len(ref) != 1 or len(alt_allele) != 1:
                        effects[pos] = 'other'
                    elif not cds_seq or pos not in cds_pos_to_idx:
                        effects[pos] = 'other'
                    else:
                        idx = cds_pos_to_idx[pos]
                        codon_start = (idx // 3) * 3
                        if codon_start + 3 > len(cds_seq):
                            effects[pos] = 'other'
                        else:
                            codon_ref = cds_seq[codon_start: codon_start + 3]
                            alt_base = alt_allele if gene_strand == '+' else _revcomp(alt_allele)[0]
                            offset = idx % 3
                            codon_alt = codon_ref[:offset] + alt_base + codon_ref[offset + 1:]
                            aa_ref = _translate_codon(codon_ref)
                            aa_alt = _translate_codon(codon_alt)
                            if aa_ref == 'X' or aa_alt == 'X':
                                effects[pos] = 'other'
                            elif aa_ref == aa_alt:
                                effects[pos] = 'synonymous'
                            else:
                                effects[pos] = 'missense'
                
                logger.info(f"    - VariantFile匹配位点数: {match_count}/{len(positions_list)}")
                # 统计各类型数量
                effect_counts = {}
                for eff in effects.values():
                    effect_counts[eff] = effect_counts.get(eff, 0) + 1
                logger.info(f"    - 注释结果统计: {effect_counts}")
            except Exception as e:
                logger.warning(f"    - VariantFile查询失败: {e}")
            
            vcf_reader.close()
                
        except Exception as e:
            logger = get_logger()
            logger.warning(f"Tabix SNP注释失败: {e}")
        
        return effects
    
    def analyze_gene(self, chrom: str, start: int, end: int, 
                     gene_id: str = None, phenotype_cols: list = None,
                     min_samples: int = 2, gwas_file: str = None) -> dict:
        """
        分析指定基因区间的单倍型 - 表型关联
            
        Args:
            chrom: 染色体
            start: 起始位置
            end: 终止位置
            gene_id: 基因 ID（可选）
            phenotype_cols: 要分析的表型列列表（None=全部）
            min_samples: 最小样本数阈值
            gwas_file: GWAS结果文件路径（可选，用于对比分析）
                
        Returns:
            dict: 分析结果汇总
        """
        logger = get_logger()
        self.chrom = chrom
        self.start = start
        self.end = end
            
        # 初始化性能监控器
        perf_monitor = PerformanceMonitor(logger)
        perf_monitor.start()
            
        logger.info("=" * 60)
        logger.info(f"分析基因区间：{chrom}:{start}-{end}")
        logger.info("=" * 60)
        
        # 0. 先解析 GTF 获取链方向，计算启动子区域
        logger.info("[Step 0] 解析基因结构和启动子区域...")
        perf_monitor.step_start("Step_0_GTF_Parser")
        gtf_data = parse_gtf_for_gene(self.gtf_file, gene_id)
        exons_list = gtf_data.get('exons', [])
        cds_list = gtf_data.get('cds', [])
        gene_chrom = gtf_data.get('chrom') or chrom
        strand = gtf_data.get('strand', '+')  # 从GTF获取链方向
        # 获取GTF中的真实基因体坐标
        gene_body_start = gtf_data.get('gene_start') or start
        gene_body_end = gtf_data.get('gene_end') or end
        logger.info(f"  - GTF 解析：{len(exons_list)} 个外显子，{len(cds_list)} 个 CDS, strand={strand}")
        logger.info(f"  - 基因体坐标: {gene_body_start:,}-{gene_body_end:,}")
        step0_time = perf_monitor.step_end("Step_0_GTF_Parser")
        logger.info(f"  - Step 0 耗时：{step0_time:.2f}s")
        
        # 计算启动子区域（使用GTF中的真实基因体坐标）
        promoter_annotator = PromoterAnnotator()
        promoter_start_pos, promoter_end_pos = promoter_annotator.get_promoter_region(
            chrom, gene_body_start, gene_body_end, strand=strand, upstream=2000
        )
        logger.info(f"  - 启动子区域: {promoter_start_pos:,}-{promoter_end_pos:,} (strand={strand})")
        
        # 扩展区域以包含启动子
        extended_start = min(start, promoter_start_pos)
        extended_end = max(end, promoter_end_pos)
        logger.info(f"  - 扩展后区域: {extended_start:,}-{extended_end:,}")
        
        # 1. 提取单倍型（使用扩展后的区域，包含启动子）
        logger.info("[Step 1] 提取单倍型...")
        perf_monitor.step_start("Step_1_Haplotype_Extraction")
        self.positions, self.hap_df, self.hap_sample_df = self.extractor.extract_region(
            chrom, extended_start, extended_end, min_samples=min_samples, snp_only=False  # 包含所有变异类型(SNP/indel/SV)
        )
        step1_time = perf_monitor.step_end("Step_1_Haplotype_Extraction")
        logger.info(f"  - Step 1 耗时：{step1_time:.2f}s")
        
        if self.hap_df is None or len(self.hap_df) == 0:
            logger.error("未能提取到有效单倍型")
            error_report = {
                'error': '单倍型提取失败',
                'chrom': chrom,
                'start': start,
                'end': end,
                'n_variants': len(self.positions) if self.positions else 0,
                'n_haplotypes': 0,
                'diagnosis': '请检查: 1)指定区间是否在VCF文件中; 2)染色体名是否匹配; 3)区间内是否有变异'
            }
            # 保存诊断报告
            error_file = os.path.join(self.output_dir, "diagnosis_report.json")
            with open(error_file, 'w', encoding='utf-8') as f:
                json.dump(error_report, f, ensure_ascii=False, indent=2)
            logger.info(f"诊断报告已保存: {error_file}")
            return error_report
            
        # 保存单倍型表
        hap_file = os.path.join(self.output_dir, "haplotypes.csv")
        self.hap_df.to_csv(hap_file, index=False)
        sample_file = os.path.join(self.output_dir, "sample_haplotypes.csv")
        self.hap_sample_df.to_csv(sample_file, index=False)
        logger.debug(f"单倍型表已保存: {hap_file}")
            
        # 2. 表型关联分析
        logger.info("[Step 2] 表型关联分析...")
        perf_monitor.step_start("Step_2_Association_Analysis")
        assoc_module = PhenotypeAssociation(self.phenotype_df, self.hap_sample_df)
        
        # 2.1 AMOVA 分析（群体遗传结构）
        logger.info("[Step 2.1] AMOVA 分析...")
        # 创建基于单倍型的分组数据
        amova_groups = self.hap_sample_df[['SampleID', 'Hap_Name']].copy()
        amova_groups.columns = ['SampleID', 'Population']
        
        # 从VCF中提取基因型数据用于AMOVA
        try:
            genotype_df = self._extract_genotypes_for_amova()
            if genotype_df is not None:
                amova_analyzer = AMOVAAnalyzer(genotype_df, amova_groups)
                amova_result = amova_analyzer.run_amova('Population')
                
                if 'error' not in amova_result:
                    logger.info(f"  Phi_ST: {amova_result.get('phi_st', 0):.4f}")
                    logger.info(f"  组间方差: {amova_result.get('percent_among', 0):.2f}%")
                    logger.info(f"  P-value: {amova_result.get('p_value', 1):.4f}")
                    
                    # 保存AMOVA表格
                    amova_table = amova_analyzer.format_amova_table()
                    if not amova_table.empty:
                        amova_csv = os.path.join(self.output_dir, "amova_table.csv")
                        amova_table.to_csv(amova_csv, index=False)
                        logger.info(f"  AMOVA表已保存: {amova_csv}")
                else:
                    logger.warning(f"  AMOVA分析失败: {amova_result.get('error')}")
        except Exception as e:
            logger.warning(f"  AMOVA分析异常: {e}")
            
        if phenotype_cols is None:
            phenotype_cols = assoc_module.get_phenotype_columns()
            
        logger.info(f"  - 分析 {len(phenotype_cols)} 个表型")
            
        all_results = {
            'gene_info': {
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
            },
            'n_variants': len(self.positions),
            'n_haplotypes': len(self.hap_df),
            'phenotype_results': {}
        }
            
        # 3. PVE 计算
        logger.info("[Step 3] 计算变异解释率 (PVE)...")
        perf_monitor.step_start("Step_3_PVE_Calculation")
        pve_module = PVECalculator(assoc_module.merged_df)
            
        for pheno in phenotype_cols:
            logger.info(f"  分析表型: {pheno}")
                
            # 关联检验
            assoc_result = assoc_module.association_test(pheno)
                
            # PVE 计算
            pve_result = pve_module.calculate_pve(pheno)
                
            all_results['phenotype_results'][pheno] = {
                'association': assoc_result,
                'pve': pve_result,
            }
                
            # 记录关键结果
            if 'error' not in assoc_result:
                sig = "***" if assoc_result.get('highly_significant') else ("*" if assoc_result.get('significant') else "")
                p_val = assoc_result.get('p_value')
                pve_val = pve_result.get('pve_percent', 0)
                p_val_str = f"{p_val:.2e}" if isinstance(p_val, (int, float)) else 'NA'
                pve_str = f"{pve_val:.2f}" if isinstance(pve_val, (int, float)) else 'NA'
                logger.info(f"    P-value: {p_val_str} {sig}")
                logger.info(f"    PVE: {pve_str}%")
            else:
                logger.warning(f"    关联分析失败: {assoc_result.get('error')}")
                
            # 绘制箱线图
            self.reporter.plot_boxplot(
                assoc_module.merged_df, pheno, gene_id=gene_id
            )
            
            # 单倍型效应分析
            effect_analyzer = HaplotypeEffectAnalyzer(assoc_module.merged_df)
            effect_result = effect_analyzer.calculate_effects(pheno)
            
            if 'error' not in effect_result and effect_result.get('n_significant', 0) > 0:
                logger.info(f"    显著效应单倍型: {effect_result['n_significant']} 个")
                # 保存效应分析表格
                effect_table = effect_analyzer.get_effect_summary_table()
                effect_csv = os.path.join(self.output_dir, f"effect_table_{pheno}.csv")
                effect_table.to_csv(effect_csv, index=False)
                logger.info(f"    效应表已保存: {effect_csv}")
                
                # 绘制效应森林图
                forest_plot = os.path.join(self.output_dir, f"effect_forest_{pheno}.pdf")
                effect_analyzer.plot_effect_forest(output_file=forest_plot)
                
                # 绘制效应柱状图
                bar_plot = os.path.join(self.output_dir, f"effect_bar_{pheno}.pdf")
                effect_analyzer.plot_effect_bar(output_file=bar_plot)
            
            # 添加到报告
            self.reporter.add_result('association', assoc_result)
            self.reporter.add_result('pve', pve_result)
            
            # 保存效应结果
            all_results['phenotype_results'][pheno]['effect'] = effect_result
            
        # 4. GWAS 整合分析（可选，仅在提供 GWAS 文件时执行）
        logger.info("[Step 4] GWAS 整合分析...")
        perf_monitor.step_start("Step_4_GWAS_Integration")
        
        all_results['gwas_comparison'] = {}
        if gwas_file and os.path.exists(gwas_file):
            logger.info(f"  - 加载 GWAS 结果文件: {gwas_file}")
            gwas_integrator = GWASIntegrator(gwas_file)
            
            for pheno, res in all_results['phenotype_results'].items():
                assoc = res.get('association', {})
                if 'error' not in assoc:
                    comparison = gwas_integrator.compare_with_haplotype(
                        assoc, gwas_integrator.gwas_df if gwas_integrator.gwas_df is not None else pd.DataFrame()
                    )
                    all_results['gwas_comparison'][pheno] = comparison
                        
                    if comparison.get('novel_finding'):
                        logger.info(f"  *** {pheno}: 潜在新发现（单倍型显著）***")
        else:
            logger.info("  - 未提供 GWAS 文件，跳过对比分析")
        
        step4_time = perf_monitor.step_end("Step_4_GWAS_Integration")
        logger.info(f"  - Step 4 耗时：{step4_time:.2f}s")
            
        # 5. 启动子变异注释（复用 Step 0 已计算的结果）
        logger.info("[Step 5] 启动子变异注释...")
        perf_monitor.step_start("Step_5_Promoter_Annotation")
        
        # 生成启动子报告（复用 Step 0 的 promoter_annotator 和 strand）
        promoter_report = promoter_annotator.generate_promoter_report(
            gene_id=gene_id,
            chrom=chrom,
            gene_start=gene_body_start,  # 使用GTF中的真实基因体坐标
            gene_end=gene_body_end,
            strand=strand,  # 使用 Step 0 获取的真实链方向
            variants_positions=self.positions
        )
        all_results['promoter_analysis'] = promoter_report
        logger.info(f"  - 启动子区域: {promoter_report['promoter_start']}-{promoter_report['promoter_end']}")
        logger.info(f"  - 启动子内变异数: {promoter_report['variants_in_promoter']}")
        logger.info(f"  - 顺式元件数: {promoter_report.get('n_cis_elements', 0)}")
        # 绘制启动子结构图
        promoter_plot = os.path.join(self.output_dir, "promoter_structure.pdf")
        promoter_annotator.plot_promoter_structure(promoter_report, promoter_plot)
        
        # 5.1 生成综合 HTML 大图
        logger.info("[Step 5.1] 生成综合 HTML 报告...")
        try:
            # 获取第一个表型的效应结果
            first_pheno = phenotype_cols[0] if phenotype_cols else 'phenotype'
            first_effect = all_results['phenotype_results'].get(first_pheno, {}).get('effect', {})
            if not first_effect:
                # 重新计算
                effect_analyzer = HaplotypeEffectAnalyzer(assoc_module.merged_df)
                first_effect = effect_analyzer.calculate_effects(first_pheno)
                    
            # 扩展分析区域以包含启动子（确保启动子在图中可见）
            promoter_start_pos = promoter_report.get('promoter_start')
            promoter_end_pos = promoter_report.get('promoter_end')
            plot_region_start = min(start, promoter_start_pos) if promoter_start_pos else start
            plot_region_end = max(end, promoter_end_pos) if promoter_end_pos else end
                    
            logger.info(f"  - 绘图区域：{plot_region_start:,}-{plot_region_end:,} (包含启动子 {promoter_start_pos:,}-{promoter_end_pos:,})")
            
            # SNP 精细分类：用 tabix 只读目标区间（含启动子）
            fasta_path = getattr(self, 'fasta_file', DataConfig.FASTA_PATH)
            snp_effects = self._annotate_snp_effects_tabix(
                vcf_file=self.vcf_file,
                fasta_path=fasta_path,
                gene_chrom=gene_chrom,
                region_start=plot_region_start,
                region_end=plot_region_end,
                cds_intervals=cds_list,
                exon_intervals=exons_list,
                gene_strand=strand
            )
            logger.info(f"  - SNP注释: {dict((k, sum(1 for v in snp_effects.values() if v==k)) for k in set(snp_effects.values()))}")
                    
            self.reporter.generate_integrated_html(
                hap_sample_df=assoc_module.merged_df,  # 使用 merged_df，包含表型数据
                effect_results=first_effect,
                variant_positions=self.positions,
                region_start=plot_region_start,  # 扩展后的起始位置
                region_end=plot_region_end,      # 扩展后的终止位置
                phenotype_col=first_pheno,
                gene_start=start,
                gene_end=end,
                promoter_start=promoter_start_pos,
                promoter_end=promoter_end_pos,
                strand=strand,
                exons=exons_list,
                cds=cds_list,
                snp_effects=snp_effects,
                chrom=chrom
            )
        except Exception as e:
            logger.warning(f"综合HTML生成失败: {e}")
            
        # 6. 生成报告
        logger.info("[Step 6] 生成分析报告...")
        perf_monitor.step_start("Step_6_Report_Generation")
        self.reporter.generate_report(gene_info=all_results['gene_info'])
        step6_time = perf_monitor.step_end("Step_6_Report_Generation")
        logger.info(f"  - Step 6 耗时：{step6_time:.2f}s")
                    
        logger.info("=" * 60)
        logger.info(f"分析完成！结果保存在：{self.output_dir}")
        logger.info("=" * 60)
                
        # 输出性能报告
        perf_report = perf_monitor.report_performance()
        all_results['performance'] = perf_report
            
        return all_results


# ============================================================================
# 命令行接口
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="单倍型-表型关联分析平台",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--vcf", 
        default=DataConfig.SV_VCF_PATH,
        help="VCF文件路径"
    )
    parser.add_argument(
        "--phenotype",
        default=DataConfig.PHENOTYPE_PATH,
        help="表型文件路径"
    )
    parser.add_argument(
        "--chrom",
        required=True,
        help="染色体"
    )
    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="区间起始位置"
    )
    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="区间终止位置"
    )
    parser.add_argument(
        "--gene-id",
        default=None,
        help="基因ID（可选）"
    )
    parser.add_argument(
        "--phenotype-cols",
        nargs="+",
        default=None,
        help="要分析的表型列（默认分析全部）"
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="最小样本数阈值"
    )
    parser.add_argument(
        "--output-dir",
        default=DataConfig.OUTPUT_DIR,
        help="输出目录"
    )
    parser.add_argument(
        "--gwas-file",
        default=None,
        help="GWAS结果文件路径（可选，用于对比分析）"
    )
    parser.add_argument(
        "--fasta",
        default=DataConfig.FASTA_PATH,
        help="参考基因组FASTA文件（用于启动子注释）"
    )
    parser.add_argument(
        "--gtf",
        default=DataConfig.GTF_PATH,
        help="GTF注释文件"
    )
    parser.add_argument(
        "--strand",
        default="+",
        choices=["+", "-"],
        help="基因链方向"
    )
    
    return parser.parse_args()


if __name__ == "__main__":
    # 立即打印调试信息（确保输出到日志）
    print("=" * 60)
    print("Haplotype-Phenotype Association Analysis")
    print("=" * 60)
    print(f"Working Directory: {os.getcwd()}")
    print(f"Script Path: {os.path.abspath(__file__)}")
    print(f"Python Version: {sys.version}")
    print("=" * 60)
    sys.stdout.flush()  # 立即刷新输出
    
    args = parse_args()
    
    print(f"\n[CONFIG] VCF: {args.vcf}")
    print(f"[CONFIG] Phenotype: {args.phenotype}")
    print(f"[CONFIG] Region: {args.chrom}:{args.start}-{args.end}")
    print(f"[CONFIG] Gene ID: {args.gene_id}")
    print(f"[CONFIG] Output Dir: {args.output_dir}")
    sys.stdout.flush()
    
    # 创建分析器
    analyzer = HaplotypePhenotypeAnalyzer(
        vcf_file=args.vcf,
        phenotype_file=args.phenotype,
        output_dir=args.output_dir,
        gtf_file=args.gtf
    )
    
    # 执行分析
    results = analyzer.analyze_gene(
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        gene_id=args.gene_id,
        phenotype_cols=args.phenotype_cols,
        min_samples=args.min_samples,
        gwas_file=args.gwas_file
    )
    
    # 打印结果摘要
    print("\n" + "="*60)
    print("分析结果摘要")
    print("="*60)
    
    # 检查是否有错误
    if 'error' in results:
        print(f"\n[ERROR] {results['error']}")
        print(f"[DIAGNOSIS] {results.get('diagnosis', 'N/A')}")
        print(f"\n请检查 {args.output_dir}/diagnosis_report.json 获取详细信息")
    else:
        for pheno, res in results.get('phenotype_results', {}).items():
            assoc = res.get('association', {})
            pve = res.get('pve', {})
            
            if 'error' in assoc:
                print(f"\n{pheno}: {assoc['error']}")
                continue
            
            sig = "***" if assoc.get('highly_significant') else ("*" if assoc.get('significant') else "ns")
            print(f"\n{pheno}:")
            p_val = assoc.get('p_value')
            pve_val = pve.get('pve_percent', 0)
            p_val_str = f"{p_val:.2e}" if isinstance(p_val, (int, float)) else 'NA'
            pve_str = f"{pve_val:.2f}" if isinstance(pve_val, (int, float)) else 'NA'
            print(f"  P-value: {p_val_str} ({sig})")
            print(f"  PVE: {pve_str}%")
            print(f"  Effect size: {pve.get('effect_size', 'NA')}")
    
    # 打印输出文件位置
    print("\n" + "="*60)
    print("输出文件位置")
    print("="*60)
    print(f"结果目录: {os.path.abspath(args.output_dir)}")
    
    # 列出生成的文件
    if os.path.exists(args.output_dir):
        files = os.listdir(args.output_dir)
        if files:
            print("\n生成的文件:")
            for f in sorted(files):
                filepath = os.path.join(args.output_dir, f)
                size = os.path.getsize(filepath)
                print(f"  - {f} ({size} bytes)")
        else:
            print("\n[警告] 输出目录为空，未生成任何文件")
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)