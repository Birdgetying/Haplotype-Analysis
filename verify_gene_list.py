# 验证基因列表
TARGET_GENES = [
    # K-SNP_Indel_CandidateGene
    "HORVU.MOREX.r3.1HG0081480",
    "HORVU.MOREX.r3.2HG0213430",
    "HORVU.MOREX.r3.2HG0214000",
    "HORVU.MOREX.r3.3HG0289370",
    "HORVU.MOREX.r3.4HG0341890",
    "HORVU.MOREX.r3.4HG0356420",
    "HORVU.MOREX.r3.4HG0356370",
    "HORVU.MOREX.r3.4HG0398010",
    "HORVU.MOREX.r3.4HG0415770",
    "HORVU.MOREX.r3.6HG0625930",
    "HORVU.MOREX.r3.6HG0626440",
    "HORVU.MOREX.r3.7HG0705960",
    "HORVU.MOREX.r3.7HG0733200",
    "HORVU.MOREX.r3.4HG0415690",
    "HORVU.MOREX.r3.4HG0415780",
    "HORVU.MOREX.r3.5HG0492640",
    # K-SV（自己结果）
    "HORVU.MOREX.r3.6HG0542030",
    "HORVU.MOREX.r3.6HG0542250",
    # "HORVU.MOREX.r3.4HG0356420",  # 重复，已在上
    "HORVU.MOREX.r3.4HG0356430",
    "HORVU.MOREX.r3.3HG0299830",
    "HORVU.MOREX.r3.3HG0299890",
    "HORVU.MOREX.r3.1HG0080620",
    "HORVU.MOREX.r3.1HG0080220",
    "HORVU.MOREX.r3.2HG0096990",
    "HORVU.MOREX.r3.7HG0639230",
    "HORVU.MOREX.r3.5HG0432290",
    "HORVU.MOREX.r3.5HG0490420",
    "HORVU.MOREX.r3.5HG0490030",
    "HORVU.MOREX.r3.5HG0490190",
    "HORVU.MOREX.r3.6HG0626350",
]

print('=' * 80)
print('genome_wide_haplotype_scan.py 中的基因列表验证')
print('=' * 80)
print(f'\n基因总数: {len(TARGET_GENES)} 个')

print('\n' + '=' * 80)
print('基因列表详情:')
print('=' * 80)
for i, gene in enumerate(TARGET_GENES, 1):
    print(f'{i:2d}. {gene}')

# 检查是否有重复
unique_genes = set(TARGET_GENES)
if len(unique_genes) == len(TARGET_GENES):
    print(f'\n✓ 无重复基因')
else:
    print(f'\n⚠ 发现重复基因！')
    from collections import Counter
    counts = Counter(TARGET_GENES)
    for gene, count in counts.items():
        if count > 1:
            print(f'  - {gene}: 出现 {count} 次')

# 按染色体排序显示
print('\n' + '=' * 80)
print('按染色体排序:')
print('=' * 80)
sorted_genes = sorted(TARGET_GENES)
for i, gene in enumerate(sorted_genes, 1):
    print(f'{i:2d}. {gene}')
