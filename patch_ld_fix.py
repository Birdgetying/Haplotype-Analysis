"""
直接修复haplotype_phenotype_analysis.py中的LD r²矩阵计算逻辑
将错误的geno_vec编码改为正确的all_bases两步编码
"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 旧代码（错误的编码方式）
old_ld_code = '''            try:
                # 构建每个样本在每个显示位置的基因型向量（0=参考，1=变异）
                # 参考碱基：取最多的等位基因作为参考
                # 使用display_orig_indices从序列中提取正确的索引
                geno_matrix = []  # shape: (n_samples, n_positions)
                for _, sample_row in hap_sample_df.iterrows():
                    seq = str(sample_row.get('Haplotype_Seq', '')).replace('|', '')
                    geno_vec = []
                    for pos_i, orig_idx in enumerate(display_orig_indices):
                        if orig_idx < len(seq):
                            base = seq[orig_idx].upper()
                            # 非参考碱基 = 变异（1），参考碱基 = 0
                            # 用'N'和'R'（参考）作为参考状态
                            geno_vec.append(0 if base in ('N',) else 1)
                        else:
                            geno_vec.append(0)
                    geno_matrix.append(geno_vec)
                
                # 修正：确定每列的参考等位基因（出现最多的碱基 = 参考）
                # 参考碱基频率最高的作为0
                geno_array = np.array(geno_matrix, dtype=float)  # (n_samples, n_positions)
                
                # 重新确定每列的0/1编码：让列均值<=0.5（即使最多的碱基为参考）
                for col_i in range(geno_array.shape[1]):
                    col_mean = np.nanmean(geno_array[:, col_i])
                    if col_mean > 0.5:
                        geno_array[:, col_i] = 1 - geno_array[:, col_i]'''

# 新代码（正确的两步编码）
new_ld_code = '''            try:
                # ---- 第一步：收集每个样本在每个位置的实际碱基 ----
                # hap_sample_df 每行是一个样本（不是单倍型类型）
                all_bases = []  # shape: (n_samples, n_positions)
                for _, sample_row in hap_sample_df.iterrows():
                    seq = str(sample_row.get('Haplotype_Seq', '')).replace('|', '')
                    row_bases = []
                    for orig_idx in display_orig_indices:
                        if orig_idx < len(seq):
                            base = seq[orig_idx].upper()
                            # N/-/? 视为缺失
                            row_bases.append(base if base not in ('N', '-', '?') else None)
                        else:
                            row_bases.append(None)
                    all_bases.append(row_bases)

                # ---- 第二步：对每列确定多数碱基=参考(0)，其余=1，缺失=NaN ----
                n_samples_ld = len(all_bases)
                n_cols_ld = len(display_orig_indices)
                geno_array = np.full((n_samples_ld, n_cols_ld), np.nan)
                for col_i in range(n_cols_ld):
                    col_bases = [all_bases[si][col_i] for si in range(n_samples_ld)]
                    counts = {}
                    for b in col_bases:
                        if b is not None:
                            counts[b] = counts.get(b, 0) + 1
                    if not counts:
                        continue
                    ref_base = max(counts, key=counts.get)  # 最多碱基 = 参考
                    for si in range(n_samples_ld):
                        b = all_bases[si][col_i]
                        if b is None:
                            geno_array[si, col_i] = np.nan
                        elif b == ref_base:
                            geno_array[si, col_i] = 0.0
                        else:
                            geno_array[si, col_i] = 1.0'''

if old_ld_code in content:
    content = content.replace(old_ld_code, new_ld_code, 1)
    print("成功替换LD r²编码逻辑")
else:
    print("未找到旧代码！检查内容...")
    idx = content.find('geno_matrix = []  # shape: (n_samples, n_positions)')
    print(f"geno_matrix=[] 在字符 {idx}")
    if idx >= 0:
        print(repr(content[idx-200:idx+400]))

with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("文件保存完成")

# 验证
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    verify = f.read()
print("验证 all_bases:", 'all_bases = []' in verify)
print("验证 row_bases:", 'row_bases = []' in verify)
print("验证 geno_matrix已删除:", 'geno_matrix = []  # shape: (n_samples, n_positions)' not in verify)

