# 单倍型打分模型公式详解 (v2)

## 总体框架

对每个单倍型 \(h\)，综合得分由 **6 个加权组件**求和：

\[
S(h) = \sum_{c \in \mathcal{C}} w_c \cdot \text{Norm}\big(S_c(h)\big)
\]

组件集合及默认权重：

| 组件 \(c\) | 默认权重 \(w_c\) | 含义 |
|------------|-----------------|------|
| `variant_effect` | 1.0 | LD剪枝后的变异功能严重度 |
| `burden` | 1.0 | LD剪枝后的稀有变异负荷 |
| `multi_omics` | 0.8 | 多组学协变量标准化偏差 |
| `fine_mapping` | 0.8 | 近似贝叶斯因子 + LD剪枝 |
| `effect_size` | 0.9 | 贝叶斯收缩后的效应量 |
| `genetic_distinct` | 0.5 | 功能加权的遗传分化距离 |

权重通过 `component_weights` 参数可配置。

---

## 归一化：自适应 Winsorize + MinMax

归一化函数 \(\text{Norm}(\cdot)\) 将各组件的原始得分映射到 \([0, 1]\)：

```
若 n_haplotypes < 10:  不剪裁，直接 MinMax
若 n_haplotypes < 30:  5%/95% 分位数剪裁
若 n_haplotypes ≥ 30:  2.5%/97.5% 标准剪裁
```

\[
x' = \text{clip}(x,\ p_{\text{lower}},\ p_{\text{upper}}), \qquad
\text{Norm}(x) = \frac{x' - p_{\text{lower}}}{p_{\text{upper}} - p_{\text{lower}}}
\]

当 \(p_{\text{upper}} = p_{\text{lower}}\)（所有得分相等），归一化为全零。

---

## LD Block 构建

位点 \(i\) 和 \(j\) 若满足 \(r^2(i,j) \geq 0.5\) 则归入同一 LD block（Union-Find 算法）。

两种剪枝策略：

| 策略 | 方法 | 使用组件 |
|------|------|----------|
| **硬剪枝 (BlockMax)** | block 内只保留贡献最大的位点，其余置零 | variant_effect, burden |
| **软剪枝 (SoftPrune)** | LD 冗余位点的 PIP 乘以 \((1-r^2)\) 降权 | fine_mapping |

---

## 组件 1：VariantEffectScore

### 公式

\[
S_{\text{ve}}(h) = \sum_{i \in \text{ALT}(h)} \text{BlockMax}\big(i,\ f_{\text{func}}(pos_i)\big)
\]

### 位点功能权重 \(f_{\text{func}}\)

功能分类按以下优先级确定：`snp_effects` 注释 → `variant_info` (SV/indel检测) → 位置推断 (CDS > UTR > intron > promoter > intergenic)。

启动子按距 TSS 距离分级：core(≤50bp) > proximal(≤200bp) > distal(≤1000bp) > promoter(>1000bp)。

| 功能类别 | 权重 |
|----------|------|
| missense_non_conservative | 5.0 |
| frameshift | 5.0 |
| missense_semi_conservative | 4.0 |
| INS (小indel) | 3.5 |
| DEL (小indel) | 3.5 |
| missense_conservative | 3.0 |
| stop_gain | 3.0 |
| SV (结构变异) | 3.0 |
| promoter_core (≤50bp) | 2.5 |
| UTR | 2.0 |
| splice_region | 2.0 |
| promoter_proximal (≤200bp) | 2.0 |
| promoter_distal (≤1000bp) | 1.5 |
| promoter (>1000bp) | 1.0 |
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

- 只有 MAF < 5% 的稀有变异参与计算
- 稀有度权重 \(-\log_{10}(\text{MAF})\)：MAF=0.01 → 2.0; MAF=0.001 → 3.0
- 常见变异（MAF ≥ 5%）贡献为零

---

## 组件 3：MultiOmicsScore

### 公式

\[
S_{\text{mo}}(h) = \frac{\sum_k \nu_k \cdot D_k(h)}{\sum_k \nu_k}
\]

\(D_k(h)\) 为第 \(k\) 个协变量组件的偏差得分，\(\nu_k\) 为其信息量权重。

### 3a. 数值型协变量 — PCA + Mahalanobis 距离

所有数值型协变量经标准化后 PCA，保留解释 ≥95% 方差的主成分。对每条 PC 计算标准化偏差：

\[
D_{\text{PCA}}(h) = \sqrt{\frac{1}{M}\sum_{m=1}^{M}\left(\frac{\bar{x}_{h,m} - \bar{x}_{\cdot,m}}{\sigma_{\cdot,m}}\right)^2}
\]

\[
\nu_{\text{PCA}} = n_{\text{features}} \times r^2_{\text{retained}}
\]

其中 \(n_{\text{features}}\) 为原始数值列数，\(r^2_{\text{retained}}\) 为保留主成分的累积方差比。

若 sklearn 不可用，回退为逐列标准化偏差的均值（不处理协变量间相关性）。

### 3b. 二分类协变量 — 比例标准化偏差

\[
D_{\text{bin}}(h) = \frac{|p_h - p_{\text{all}}|}{\sqrt{p_{\text{all}}(1 - p_{\text{all}})}}
\]

\[
\nu_{\text{bin}} = H_{\text{norm}} = \frac{-p\log_2 p - (1-p)\log_2(1-p)}{\log_2 2}
\]

信息量权重 = 归一化香农熵。平衡类别（50%/50%）→ 1.0；极端分布（99%/1%）→ ≈0。

### 3c. 多分类协变量 — one-hot 哑变量偏差均值

\[
D_{\text{cat}}(h) = \frac{1}{K}\sum_{k=1}^{K}\frac{|p_{h,k} - p_{\text{all},k}|}{\sqrt{p_{\text{all},k}(1 - p_{\text{all},k})}}
\]

\[
\nu_{\text{cat}} = \frac{-\sum_{i=1}^{K} p_i \log_2 p_i}{\log_2 K}
\]

跳过样本数 <2 的类别哑变量。

---

## 组件 4：FineMappingScore

### 原理

参考 Wakefield (2009) ABF 近似，将 GWAS 的 \(p\)-value 转换为后验包含概率 (PIP)，再按单倍型携带的 ALT 等位基因汇总。

### Step 1 — p → Z → logBF

\[
p \rightarrow Z^2 = \chi^2_{1,1-p} \quad\text{(使用 chi2.isf 避免 } 1-p \text{ 精度丢失)}
\]

\[
\log\text{BF}(pos) = \frac{Z^2}{2}
\]

无 GWAS 信号 (\(\log_{10}p = 0\)) 的位点：\(\log\text{BF} = -\infty\)。

### Step 2 — Log-sum-exp → PIP

\[
\text{PIP}(pos) = \frac{\exp\big(\log\text{BF}(pos) - \max\text{BF}\big)}{\sum_j \exp\big(\log\text{BF}(pos_j) - \max\text{BF}\big)}
\]

### Step 3 — LD 软剪枝

按 PIP 降序遍历所有位点：若与已选中位点的最大 \(r^2 \geq 0.5\)：

\[
\text{PIP}_{\text{pruned}}(pos) = \text{PIP}(pos) \times (1 - r^2_{\max})
\]

否则保留完整 PIP。PIP < 10⁻⁶ 的位点被过滤，剪枝后 < 10⁻⁴ 的位点被丢弃。

### Step 4 — 按单倍型汇总

\[
S_{\text{fm}}(h) = \sum_{i \in \text{ALT}(h)} \text{PIP}_{\text{pruned}}(pos_i)
\]

### 提前退出

若所有位点的 GWAS \(-\log_{10}p = 0\)（无有效信号），直接返回全零，跳过所有计算。

---

## 组件 5：EffectSizeScore

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

当 \(z^2_{\text{sum}} < k-2\) 时 \(\lambda = 0\)，所有效应量完全收缩至精度加权均值。

---

## 组件 6：GeneticDistinctivenessScore

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

单倍型数 < 2 时返回全零。

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
        'Hap1': {'variant_effect': 0.85, 'burden': 0.32, ... 'total': 2.14},
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
    'component_weights': {...},
    'confidence_level': 'moderate',
    'low_confidence': False,
    'pve': 0.12,
    'circularity_warning': False,
}
```
