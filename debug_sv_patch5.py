with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 检查 sv_variants 周围的上下文
idx = c.find("sv_variants.vcf.gz")
print("sv_variants at:", idx)
print(repr(c[max(0,idx-200):idx+100]))

print()

# 精确查找 required_files 中的 variants.vcf.gz
idx2 = c.find("'variants.vcf.gz'")
if idx2 >= 0:
    print("Found 'variants.vcf.gz' at:", idx2)
    print("Context:", repr(c[idx2-30:idx2+100]))
