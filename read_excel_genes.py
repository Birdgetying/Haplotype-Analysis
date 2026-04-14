import pandas as pd

xl = pd.ExcelFile('database/大麦备份/大麦GWAS及SV筛选基因-2026.03.22(1).xlsx')

df1 = pd.read_excel(xl, sheet_name='K-SNP_Indel_CandidateGene')
df2 = pd.read_excel(xl, sheet_name='K-SV（自己结果）')

print('=' * 80)
print('K-SNP_Indel_CandidateGene 基因列表 (16个):')
print('=' * 80)
gene_list1 = df1.iloc[:, 0].dropna().tolist()
for i, gene in enumerate(gene_list1, 1):
    print(f"{i:2d}. {gene}")

print('\n' + '=' * 80)
print('K-SV（自己结果）基因列表 (15个):')
print('=' * 80)
gene_list2 = df2.iloc[:, 0].dropna().tolist()
for i, gene in enumerate(gene_list2, 1):
    print(f"{i:2d}. {gene}")

print('\n' + '=' * 80)
print('合并后的基因列表 (去重):')
print('=' * 80)
all_genes = list(set(gene_list1 + gene_list2))
all_genes.sort()
for i, gene in enumerate(all_genes, 1):
    print(f"{i:2d}. {gene}")

print(f'\n总计: {len(all_genes)} 个唯一基因')
