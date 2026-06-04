## 代码架构

项目是一个**单倍型-表型关联分析平台**，用于植物基因组学（小麦/大麦/水稻/玉米）。

### 核心模块（两个大脚本，约15k行）

**`haplotype_phenotype_analysis.py`** (12k行) — 核心分析引擎，被其他脚本 import 使用：
- `HaplotypeExtractor` — 从VCF提取单倍型
- `PhenotypeAssociation` — 表型关联（t-test/ANOVA/回归）
- `HaplotypeScorer` — 单倍型打分（含EB收缩 + LD软加权）
- `PopulationStructureCorrector` — PCA+ANCOVA群体结构矫正
- `GWASIntegrator` — GWAS结果对比整合
- `PromoterAnnotator` — 启动子区域功能注释
- `ReportGenerator` — 生成 integrated_analysis.html 综合报告
- `HaplotypePhenotypeAnalyzer` — 主流程编排器，`analyze_gene()` 是核心入口
- 大量变异注释函数：`annotate_snp_effects_for_region()` 等
- `DataConfig` — 全局数据路径配置（FASTA/GFF3/VCF等）

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
rice_database/{gene_id}/   — 水稻分析数据库
rice_results/{gene_id}/    — 水稻分析结果
maize_database/{gene_id}/  — 玉米分析数据库
maize_results/{gene_id}/   — 玉米分析结果
wheat_rht_database/{gene_id}/ — 小麦Rht基因数据库
wheat_rht_results/{gene_id}/  — 小麦Rht基因分析结果
database/{gene_id}/        — 大麦/小麦(CSIAAS)分析数据库
results/{gene_id}/         — 大麦/小麦(CSIAAS)分析结果
```

### 小麦数据路径 (本地, d:/Desktop/wheat_data/)

| 文件 | 大小 | 用途 |
|------|------|------|
| `Core819Samples_snp.filter.final.id_gt.813m.vcf.gz` | 10.7 GB | 小麦SNP (813样本, CS-IAAS T2T v1.1) |
| `Core819Samples_indel.filter.final.id_gt.813m.vcf.gz` | 721 MB | 小麦INDEL |
| `Core819Samples_ALL.vcf.gz` | 12.1 GB | 完整VCF (SNP+INDEL+SV) |
| `SV.new.vcf.gz` | 1.2 GB | 结构变异VCF |
| `CS-IAAS_v1.1_HC.gene.gff3` | 14.8 MB | 基因注释 |
| `CS-IAAS_v1.1.fasta` | 14 GB | 参考基因组 |
| `Phe.raw2` | 48 KB | 株高(PH)表型 (835样本, VCFID=WMGxxxA) |
| `Phe.txt` | 13 KB | TFW_DSI表型 |
| `812VCF` / `vcfID` / `VCFID.txt` | — | 样本列表/ID映射 |

小麦VCF染色体命名: Chr1A-Chr7D (CS-IAAS T2T v1.1), 样本命名: WMGxxxA/WATDExxxx

### 水稻数据路径

`D:\Desktop\data\水稻\`
- VCF: `VCF format/rice4k_geno_add_del.vcf.gz` (14 GB, 4726样本, IRGSP 1.0)
- 表型: `Phenotypes/phenos_modified.tsv` (529样本, C*/W*命名, 10个性状含Plant_height)
- GFF3: `rice_test_genes_paper.gff3` (仅GW7/DEP1/OsMADS25)
- 注意: 水稻VCF样本ID与表型ID部分匹配(C*/W*共529个交集)，IRIS_313/B*/CX*样本无表型

### 玉米数据路径 (CSV转VCF)

`d:\Desktop\论文复现\data2\`
- `Iranian_Samples.csv` → 转换得 `iranian_maize.vcf.gz` (47K标记, 2442样本, chr1)
- `phenotype_iranian.csv` → `maize_phenotype.tsv` (Heat_dtm/dth, Drought_dtm/dth)
- 玉米GFF3: `D:\Desktop\data\玉米\Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3`
- 注意: 玉米是SNP芯片数据(CloneID/CallRate/PIC)，非重测序；仅chr1上40/5892个基因有标记覆盖

### 分析脚本

| 脚本 | 用途 |
|------|------|
| `run_rice_test.py` | 水稻快速测试 (3个chr01基因) |
| `run_rice_html.py` | 水稻从数据库重新生成HTML |
| `run_rice_paper_genes.py` | 水稻论文基因(GW7/DEP1/OsMADS25) |
| `run_wheat_rht_final.py` | **小麦Rht基因×株高** 单倍型分析 (SNP+INDEL+SV合并) |
| `setup_wheat_data.py` | 小麦VCF染色体子集提取 |
| `run_genome_scan.sh` | 大麦 Morex_v3 全基因组批量扫描 (集群) |
| `run_haplotype_analysis.sh` | 小麦 CS-IAAS 单基因分析 (集群) |
| `plot_Gene_HapSeq.py` | 基因单倍型序列可视化 |