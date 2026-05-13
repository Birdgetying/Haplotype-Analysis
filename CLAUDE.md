# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

程序修改结束后提交并推送至github仓库，若网络连接有问题则不必一直尝试。
总结用中文。
修改完之后需要运行代码测试。
不要臆测，不要隐藏困惑，把权衡讲出来。
可以提出你认为合理的但与要求不同的建议让我决定。

## 代码架构

项目是一个**单倍型-表型关联分析平台**，用于植物基因组学（小麦/大麦/水稻）。

### 核心模块（两个大脚本，约15k行）

**`haplotype_phenotype_analysis.py`** (12k行) — 核心分析引擎，被其他脚本 import 使用：
- `HaplotypeExtractor` — 从VCF提取单倍型
- `PhenotypeAssociation` — 表型关联（t-test/ANOVA/回归）
- `HaplotypeScorer` — 单倍型打分
- `PopulationStructureCorrector` — PCA+ANCOVA群体结构矫正
- `GWASIntegrator` — GWAS结果对比整合
- `PromoterAnnotator` — 启动子区域功能注释
- `ReportGenerator` — 生成 integrated_analysis.html 综合报告
- `HaplotypePhenotypeAnalyzer` — 主流程编排器，`analyze_gene()` 是核心入口
- 大量变异注释函数：`annotate_snp_effects_for_region()` 等

**`genome_wide_haplotype_scan.py`** (3k行) — 全基因组批量扫描：
- `ScanConfig` — 批量扫描配置类（数据路径、阈值等，通过修改类属性覆盖）
- 从 `haplotype_phenotype_analysis` import 核心类
- `run_genome_scan()` — 从GFF3解析所有基因，对每个基因调用 `process_single_gene()`，支持多进程并行
- `generate_genome_scan_html_report()` — 全基因组汇总HTML

### 数据流

```
输入: VCF + GFF3 + 表型文件 + FASTA
  ↓ 提取基因区间变异
  ↓ 构建单倍型
  ↓ 关联分析 + 群体矫正 + GWAS对比
  ↓ 输出数据库 + HTML报告 + PDF图表
```

### 输出目录结构

```
rice_database/{gene_id}/   — 水稻分析数据库 (gene_info.json, haplotype_data.csv, phenotype_data.csv, variant_info.csv)
rice_results/{gene_id}/     — 水稻分析结果 (GW7.html + boxplot/effect_forest/pca_plot 等 PDF)
results/{gene_id}/          — 大麦/小麦分析结果
```

水稻数据源: `D:\Desktop\data\水稻\`

### 水稻脚本

| 脚本 | 用途 | 数据源 |
|------|------|--------|
| `run_rice_paper_genes.py` | **主要脚本** — 论文基因(GW7/DEP1/OsMADS25)从基因型构建数据库+生成HTML | tab-separated genotype files |
| `run_rice_html.py` | 从已有数据库重新生成HTML报告 | rice_database/ |
| `run_rice_test.py` | 快速测试 (3个chr01基因, VCF提取) | rice4k_geno_add_del.vcf.gz |
| `run_rice_from_genotype.py` | 从tab-separated基因型构建数据库(通用) | tab-separated genotype files |

### 小麦/大麦脚本

| 脚本 | 用途 |
|------|------|
| `run_genome_scan.sh` | 大麦 Morex_v3 全基因组批量扫描 (集群) |
| `run_haplotype_analysis.sh` | 小麦 CS-IAAS 单基因分析 (集群) |
| `plot_Gene_HapSeq.py` | 基因单倍型序列可视化 |

### 依赖

Python 包: `numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-learn`, `pysam` (可选，用于tabix快速VCF查询)

Windows 编码修复: 两个主脚本均在开头有 `sys.stdout/stderr` UTF-8 wrapper。

### 运行单基因分析

```bash
python haplotype_phenotype_analysis.py \
    --vcf <vcf_file> --phenotype <pheno_file> \
    --gtf <gff3_file> --fasta <fasta_file> \
    --chrom <chr> --start <pos> --end <pos> \
    --gene-id <id> --strand <+/-> \
    --output-dir <dir>
```

### 运行全基因组扫描

```bash
python genome_wide_haplotype_scan.py \
    --vcf <vcf> --sv-vcf <sv_vcf> --gff <gff3> \
    --phenotype <pheno> --fasta <fasta> \
    --database-dir <dir> --results-dir <dir>
```

### 代码测试方法

修改核心脚本后，用 `run_rice_test.py` 做快速验证（3个小基因，运行快）：

```bash
cd d:/Desktop/project1 && python run_rice_test.py
```
