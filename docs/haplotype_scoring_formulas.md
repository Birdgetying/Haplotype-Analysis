# 单倍型打分模型公式详解 (v3)

## 总体框架

对每个单倍型 \(h\)，综合得分由 **7 个加权组件**求和：

\[
S(h) = \sum_{c \in \mathcal{C}} w_c \cdot \text{Norm}\big(S_c(h)\big)
\]

组件集合及默认权重与 `HaplotypeScorer.DEFAULT_COMPONENT_WEIGHTS` 保持一致：

| 组件 \(c\) | 默认权重 \(w_c\) | 含义 |
|------------|-----------------|------|
| `variant_effect` | 1.0 | LD硬剪枝后的变异功能严重度 |
| `burden` | 1.0 | LD硬剪枝后的稀有变异负荷 |
| `eb_effect` | 0.9 | GeneBayes 启发的 EB 收缩 + LD 软加权效应 |
| `multi_omics` | 0.8 | 多组学/协变量标准化偏差 |
| `fine_mapping` | 0.8 | 近似贝叶斯因子 + LD剪枝，整合 GWAS 信号 |
| `effect_size` | 0.9 | 贝叶斯收缩后的单倍型效应量 |
| `genetic_distinct` | 0.5 | 功能加权的遗传分化距离 |

权重通过 `component_weights` 参数可配置。`gwas` 不再是独立评分组件；GWAS 信号已合并到 `fine_mapping`。

---

## 前置过滤：稀有单倍型

单倍型打分前，`HaplotypeExtractor.extract_region()` 会按 `min_samples` 过滤样本数不足的单倍型：

\[
\text{valid}(h) = \mathbb{1}[n_h \geq \text{min\_samples}]
\]

`Count < min_samples` 的单倍型及其样本不进入后续 `hap_sample_df`，因此不会参与关联分析、效应量估计和打分。默认快速测试/批量扫描中 `min_samples=2`，避免单样本单倍型造成不稳定高分。

---

## 归一化：自适应 Winsorize + MinMax

归一化函数 \(\text{Norm}(\cdot)\) 将各组件原始得分映射到 \([0, 1]\)：

```text
若 n_haplotypes < 10:  不剪裁，直接 MinMax
若 n_haplotypes < 30:  5%/95% 分位数剪裁
若 n_haplotypes ≥ 30:  2.5%/97.5% 标准剪裁
```

\[
x' = \text{clip}(x,\ p_{\text{lower}},\ p_{\text{upper}}), \qquad
\text{Norm}(x) = \frac{x' - p_{\text{lower}}}{p_{\text{upper}} - p_{\text{lower}}}
\]

当 \(p_{\text{upper}} = p_{\text{lower}}\)（所有得分相等），归一化为全零。最终总分不再二次归一化，保留组件权重累加后的绝对尺度。

---

## LD 处理策略

### LD block 构建

位点 \(i\) 和 \(j\) 若满足：

\[
r^2(i,j) \geq 0.5
\]

则归入同一 LD block（Union-Find 算法）。

### 两类剪枝/降权

| 策略 | 方法 | 使用组件 |
|------|------|----------|
| **硬剪枝 (BlockMax)** | block 内只保留贡献最大的位点，其余置零 | `variant_effect`, `burden` |
| **PIP软剪枝** | 连锁位点的 PIP 乘以 \((1-r^2)\) | `fine_mapping` |
| **EB软加权** | 若与更高贡献位点 \(r^2>0.6\)，贡献乘以 \(\max(1-\sqrt{r^2}, 0.01)\) | `eb_effect` |

---

## 组件 1：VariantEffectScore

### 公式

\[
S_{\text{ve}}(h) = \sum_{i \in \text{ALT}(h)} \text{BlockMax}\big(i,\ f_{\text{func}}(pos_i)\big)
\]

### 位点功能权重 \(f_{\text{func}}\)

功能分类按以下优先级确定：

1. `snp_effects` 精细注释（如 missense、synonymous、UTR、promoter）。
2. `variant_info` 中的大型 SV / 符号等位基因 / 长度差 ≥50bp。
3. 位置推断：CDS > UTR > intron > promoter > intergenic。

小 indel 不作为单独 raw variant type 评分项强行覆盖位置注释；它可在后验富集中解释。启动子按距 TSS 距离分级：core(≤50bp) > proximal(≤200bp) > distal(≤1000bp) > promoter(>1000bp)。

| 功能类别 | 权重 |
|----------|------|
| missense_non_conservative | 5.0 |
| frameshift | 5.0 |
| missense_semi_conservative | 4.0 |
| INS | 3.5 |
| DEL | 3.5 |
| missense_conservative | 3.0 |
| stop_gain | 3.0 |
| SV | 3.0 |
| promoter_core | 2.5 |
| UTR | 2.0 |
| splice_region | 2.0 |
| promoter_proximal | 2.0 |
| promoter_distal | 1.5 |
| promoter | 1.0 |
| intron | 0.5 |
| synonymous | 0.5 |
| intergenic | 0.1 |
| other | 0.1 |

---

## 组件 2：BurdenScore

### 公式

\[
S_{\text{burden}}(h) = \sum_{i \in \text{ALT}(h)} \text{BlockMax}\Big(i,\ f_{\text{func}}(pos_i) \times \mathbb{1}[\text{MAF}_i < 0.05] \times (-\log_{10} \text{MAF}_i)\Big)
\]

### 说明

- 只有 MAF < 5% 的稀有变异参与计算。
- 稀有度权重 \(-\log_{10}(\text{MAF})\)：MAF=0.01 → 2.0；MAF=0.001 → 3.0。
- 常见变异（MAF ≥ 5%）贡献为零。
- LD block 内只保留贡献最大的位点，避免连锁位点重复累计。

---

## 组件 3：EBEffectScore

### 原理

`eb_effect` 是 v3 新增组件，用于识别“携带稀有且表型效应较强变异”的功能单倍型，同时通过经验贝叶斯收缩抑制低 MAF 位点的噪声。

### Step 1 — 从单倍型重建 ALT 矩阵

从 `hap_sample_df` 的单倍型序列和 `variant_info` 位点列表重建：

\[
X_{s,i} = \mathbb{1}[\text{sample }s\text{ 在位点 }i\text{ 携带 ALT}]
\]

其中 `+`/`-` 视为 indel ALT，普通碱基与 `variant_info[pos]['ref']` 比较。

### Step 2 — 单变量效应回收

对每个位点计算单变量线性回归斜率的绝对值：

\[
\hat{\beta}_i = \left|\frac{\text{Cov}(X_i, y)}{\text{Var}(X_i)}\right|
\]

方差下限截断为 \(10^{-8}\)。

### Step 3 — 按 MAF 分箱的 EB 收缩

将 MAF 按 log-scale 分成最多 20 个箱。对每个箱 \(b\)：

\[
\lambda_i = \frac{\text{noise}_i}{\text{noise}_i + \text{signal}_b + 10^{-12}}
\]

\[
\hat{\beta}^{EB}_i = (1-\lambda_i)\hat{\beta}_i + \lambda_i\bar{\beta}_b
\]

其中：

\[
\text{noise}_i = \frac{1}{\text{MAF}_i(1-\text{MAF}_i) + 10^{-8}}
\]

\[
\text{signal}_b = \max\big(\text{Var}(\hat{\beta}_{i \in b}) - \overline{\text{noise}}_b,\ 0\big)
\]

当某个 MAF 箱内位点数 <5 时，代码不做收缩，保留原始效应。

### Step 4 — 稀有度 × 效应贡献

每个位点的原始贡献为：

\[
C_i = \big(-\log_{10}\max(\text{MAF}_i, 10^{-4})\big) \times \left(1 + \hat{\beta}^{EB}_i\right)
\]

### Step 5 — LD 软加权

若位点 \(i\) 与另一个更高贡献位点 \(j\) 满足 \(r^2(i,j) > 0.6\)，则：

\[
C^{LD}_i = C_i \times \max(1-\sqrt{r^2(i,j)},\ 0.01)
\]

若不存在更高贡献的强 LD 位点，则 \(C^{LD}_i = C_i\)。

### Step 6 — 按单倍型汇总

\[
S_{\text{eb}}(h) = \sum_{i \in \text{ALT}(h)} C^{LD}_i
\]

---

## 组件 4：MultiOmicsScore

### 公式

\[
S_{\text{mo}}(h) = \frac{\sum_k \nu_k \cdot D_k(h)}{\sum_k \nu_k}
\]

\(D_k(h)\) 为第 \(k\) 个协变量组件的偏差得分，\(\nu_k\) 为其信息量权重。

### 4a. 数值型协变量 — PCA 偏差

所有数值型协变量经标准化后 PCA，保留解释 ≥95% 方差的主成分。当前代码对每条 PC 计算标准化偏差后使用：

\[
D_{\text{PCA}}(h) = \frac{\sqrt{\sum_{m=1}^{M}\left(\frac{\bar{x}_{h,m} - \bar{x}_{\cdot,m}}{\sigma_{\cdot,m}}\right)^2}}{M}
\]

\[
\nu_{\text{PCA}} = n_{\text{features}} \times r^2_{\text{retained}}
\]

其中 \(n_{\text{features}}\) 为原始数值列数，\(r^2_{\text{retained}}\) 为保留主成分的累积方差比。若 sklearn 不可用，回退为逐列标准化偏差的均值。

### 4b. 二分类协变量 — 比例标准化偏差

\[
D_{\text{bin}}(h) = \frac{|p_h - p_{\text{all}}|}{\sqrt{p_{\text{all}}(1 - p_{\text{all}})}}
\]

\[
\nu_{\text{bin}} = H_{\text{norm}} = \frac{-p\log_2 p - (1-p)\log_2(1-p)}{\log_2 2}
\]

### 4c. 多分类协变量 — one-hot 哑变量偏差均值

\[
D_{\text{cat}}(h) = \frac{1}{K}\sum_{k=1}^{K}\frac{|p_{h,k} - p_{\text{all},k}|}{\sqrt{p_{\text{all},k}(1 - p_{\text{all},k})}}
\]

\[
\nu_{\text{cat}} = \frac{-\sum_{i=1}^{K} p_i \log_2 p_i}{\log_2 K}
\]

样本数 <2 的类别哑变量跳过。

---

## 组件 5：FineMappingScore

### 原理

参考 Wakefield (2009) ABF 近似，将 GWAS 的 \(p\)-value 转换为后验包含概率 (PIP)，再按单倍型携带的 ALT 等位基因汇总。

### Step 1 — p → Z² → logBF

\[
p \rightarrow Z^2 = \chi^2_{1,\text{isf}(p)}
\]

\[
\log\text{BF}(pos) = \frac{Z^2}{2}
\]

无 GWAS 信号（\(-\log_{10}p = 0\)）的位点记为无效信号。

### Step 2 — Log-sum-exp → PIP

\[
\text{PIP}(pos) = \frac{\exp\big(\log\text{BF}(pos) - \max\text{BF}\big)}{\sum_j \exp\big(\log\text{BF}(pos_j) - \max\text{BF}\big)}
\]

### Step 3 — LD 软剪枝

按 PIP 降序遍历所有位点。若与已选中位点的最大 \(r^2 \geq 0.5\)：

\[
\text{PIP}_{\text{pruned}}(pos) = \text{PIP}(pos) \times (1-r^2_{\max})
\]

否则保留完整 PIP。PIP < 10⁻⁶ 的位点先过滤，剪枝后 ≤10⁻⁴ 的位点丢弃。

### Step 4 — 按单倍型汇总

\[
S_{\text{fm}}(h) = \sum_{i \in \text{ALT}(h)} \text{PIP}_{\text{pruned}}(pos_i)
\]

若所有位点均无有效 GWAS 信号，直接返回全零。

---

## 组件 6：EffectSizeScore

### 原理

使用 James-Stein 经验贝叶斯收缩，向精度加权均值收缩小样本组的 Cohen's d，防止不稳定估计。

### Step 1 — 精确方差

\[
\text{var}(d_i) = \frac{n_1 + n_2}{n_1 \cdot n_2} + \frac{d_i^2}{2(n_1 + n_2)}
\]

其中 \(n_1\) = 该单倍型组样本数，\(n_2 = N_{\text{total}} - n_1\) = 其余样本数。方差下限截断为 10⁻⁸。

### Step 2 — 精度加权均值

\[
\text{precision}_i = \frac{1}{\text{var}(d_i)}, \qquad
\bar{d}_{pw} = \frac{\sum_i \text{precision}_i \cdot d_i}{\sum_i \text{precision}_i}
\]

### Step 3 — James-Stein 收缩因子

\[
z^2_{\text{sum}} = \sum_i \text{precision}_i \cdot (d_i - \bar{d}_{pw})^2
\]

\[
\lambda = \max\left(0,\ 1 - \frac{k-2}{z^2_{\text{sum}}}\right)
\]

其中 \(k\) 为有效单倍型数。\(k \leq 2\) 或 \(z^2_{\text{sum}} = 0\) 时 \(\lambda = 1\)（不收缩）。

### Step 4 — 收缩估计

\[
d_i^{\text{shrunk}} = \bar{d}_{pw} + \lambda \cdot (d_i - \bar{d}_{pw})
\]

\[
S_{\text{es}}(h) = |d_i^{\text{shrunk}}|
\]

---

## 组件 7：GeneticDistinctivenessScore

### 原理

功能加权的 Hamming 距离，反映每个单倍型在遗传上与其他单倍型的差异程度。LD block 内位点先压缩权重，避免连锁位点重复计数。

### Step 1 — 位点权重（LD 压缩）

\[
w(pos) = \frac{\max_{j \in \text{block}(pos)} f_{\text{func}}(pos_j)}{|\text{block}(pos)|}
\]

然后全局归一化：

\[
\tilde{w}(pos) = \frac{w(pos)}{\sum_{p} w(p)}
\]

### Step 2 — 加权 Hamming 距离

\[
D(h_i, h_j) = \sum_{k:\ \text{seq}_i[k] \neq \text{seq}_j[k]} \tilde{w}(pos_k) + \frac{|\text{len}_i - \text{len}_j|}{n_{\text{positions}}}
\]

长度差异项仅在序列长度不同时生效。

### Step 3 — 样本量加权平均

\[
S_{\text{gd}}(h_i) = \frac{\sum_{j \neq i} n_j \cdot D(h_i, h_j)}{\sum_{j \neq i} n_j}
\]

单倍型数 <2 时返回全零。

---

## 变异类型后验富集

变异类型（SNP / INDEL / SV）不作为独立 raw scoring component 加入总分。代码中 `variant_type_enrichment()` 用于后验解释：比较被选中/高贡献位点与全区域背景的变异类型分布，帮助解释评分结果，而不是直接决定分值。

---

## 置信度校准

| PVE 范围 | confidence_level | low_confidence |
|----------|-----------------|----------------|
| ≥ 20% | high | False |
| 5% – 20% | moderate | False |
| < 5% | low | True |
| 未知 (无 PVE) | unknown | False |

---

## 循环论证检测

FineMapping 组件使用 GWAS 的 \(-\log_{10}p\) 作为输入，而最终 score-phenotype 回归使用相同表型作为因变量。当以下条件同时满足时触发警告：

- PVE < 10%
- score-phenotype R² > 0.5

输出中 `circularity_warning = True` 提示用户注意此项风险。

---

## 最终输出结构

```python
{
    'per_haplotype': {
        'Hap1': {
            'variant_effect': 0.85,
            'burden': 0.32,
            'eb_effect': 0.41,
            'multi_omics': 0.20,
            'fine_mapping': 0.05,
            'effect_size': 0.66,
            'genetic_distinct': 0.18,
            'total': 2.14,
        },
        ...
    },
    'per_sample': [
        {'sample_id': 'S1', 'haplotype': 'Hap1', 'score': 2.14, 'phenotype': 10.5},
        ...
    ],
    'r_squared': 0.42,
    'regression_pvalue': 0.003,
    'slope': 1.25,
    'intercept': 0.05,
    'component_weights': {
        'variant_effect': 1.0,
        'burden': 1.0,
        'eb_effect': 0.9,
        'multi_omics': 0.8,
        'fine_mapping': 0.8,
        'effect_size': 0.9,
        'genetic_distinct': 0.5,
    },
    'confidence_level': 'moderate',
    'low_confidence': False,
    'pve': 0.12,
    'circularity_warning': False,
}
```
