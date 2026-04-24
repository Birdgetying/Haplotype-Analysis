import pandas as pd, numpy as np

df = pd.read_csv(r'results3/[NoExtend]_HORVU.MOREX.r3.1HG0081480/sample_haplotypes.csv')
seqs_clean = [s.replace('|','') for s in df['Haplotype_Seq']]
print(f'样本数: {len(seqs_clean)}, 序列长: {len(seqs_clean[0])}')

# 只取前31个位点进行LD计算（模拟代码行为）
display_orig_indices = list(range(31))

all_bases = []
for seq in seqs_clean:
    row_bases = []
    for orig_idx in display_orig_indices:
        if orig_idx < len(seq):
            base = seq[orig_idx].upper()
            row_bases.append(base if base not in ('N', '-', '?') else None)
        else:
            row_bases.append(None)
    all_bases.append(row_bases)

# 编码
n_samples = len(all_bases)
n_cols = len(display_orig_indices)
geno_array = np.full((n_samples, n_cols), np.nan)
for col_i in range(n_cols):
    col_bases = [all_bases[si][col_i] for si in range(n_samples)]
    counts = {}
    for b in col_bases:
        if b is not None:
            counts[b] = counts.get(b, 0) + 1
    if not counts:
        continue
    ref_base = max(counts, key=counts.get)
    for si in range(n_samples):
        b = all_bases[si][col_i]
        if b is None: geno_array[si,col_i] = np.nan
        elif b == ref_base: geno_array[si,col_i] = 0.0
        else: geno_array[si,col_i] = 1.0

print(f'geno_array shape: {geno_array.shape}')
for ci in range(min(5, n_cols)):
    u = np.unique(geno_array[:,ci][~np.isnan(geno_array[:,ci])])
    print(f'  列{ci} unique: {u}')

# 计算r2矩阵前5个
nonzero = 0
for i in range(min(5, n_cols)):
    for j in range(i+1, min(5, n_cols)):
        xi = geno_array[:,i]; xj = geno_array[:,j]
        mask = (~np.isnan(xi)) & (~np.isnan(xj))
        xi_m, xj_m = xi[mask], xj[mask]
        pi, pj = np.mean(xi_m), np.mean(xj_m)
        if pi > 0 and pi < 1 and pj > 0 and pj < 1:
            cov = np.mean(xi_m * xj_m) - pi * pj
            denom = (pi * (1 - pi) * pj * (1 - pj)) ** 0.5
            r2 = (cov / denom) ** 2 if denom > 0 else 0.0
            r2 = min(1.0, max(0.0, r2))
        else:
            r2 = 0.0
        if r2 > 0:
            nonzero += 1
        print(f'  r2({i},{j}): pi={pi:.3f}, pj={pj:.3f}, r2={r2:.4f}')

print(f'非零r2数: {nonzero}')
