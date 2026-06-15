# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

程序修改结束后提交并推送至github仓库，若网络连接有问题则不必一直尝试。
总结用中文。
修改完之后需要运行代码测试。
不要臆测，不要隐藏困惑，把权衡讲出来。
可以提出你认为合理的但与要求不同的建议让我决定。

### 明星基因验证记录

进行明星基因/文献功能变异/单倍型打分验证时，必须实时更新 `star_gene_validation_record.md`。
每次新增或重跑目标都要记录：目标基因、文献功能变异或单倍型、数据来源、score mode、运行命令、输出目录、top-scored haplotype、是否匹配文献功能单倍型、样本数/可靠性限制，以及 blocked 原因。
文献功能变异只用于事后验证，不得作为 discovery scoring 的输入特征。

### 项目结构见DOCS.md

### 依赖

Python 包: `numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-learn`, `pysam` (可选，Windows编译困难)

Windows 编码修复: 所有脚本均在开头有 `sys.stdout/stderr` UTF-8 wrapper。

### 已知Bug

- `SimpleVCFParser.fetch()` 返回 dict 但部分变异注释代码期望 pysam 对象(`rec.pos`), 已修复 `annotate_snp_effects_for_region()` 中一处，其他地方可能仍有
- 无pysam时大VCF(>10GB)线性扫描极慢(~65min/次)，建议创建基因区域微型VCF子集
- `DataConfig.FASTA_PATH` 默认指向HPC路径，本地运行时需在导入后覆盖:
  ```python
  from haplotype_phenotype_analysis import DataConfig
  DataConfig.GENOME_BASE = 'd:/Desktop/wheat_data/'
  ```

### 代码测试方法

修改核心脚本后，用 `run_rice_test.py` 做快速验证（3个小基因，运行快）：

```bash
cd d:/Desktop/project1 && python run_rice_test.py
```

小麦Rht基因分析：
```bash
cd d:/Desktop/project1 && python run_wheat_rht_final.py
```
