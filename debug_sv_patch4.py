with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 在 required_files 区域查找
idx = c.find("required_files = [")
print("required_files at:", idx)
if idx >= 0:
    snippet = c[idx:idx+600]
    print("CONTENT:", repr(snippet[:500]))

print()

# 检查是否已有 sv_variants
idx2 = c.find("sv_variants.vcf.gz")
print("sv_variants in content:", idx2)

# 检查 sv_vcf 相关关键字
for kw in ["sv_vcf_file", "sv_extractor", "sv_positions", "--sv-vcf"]:
    idx = c.find(kw)
    print(f"'{kw}': {'FOUND at ' + str(idx) if idx >= 0 else 'NOT FOUND'}")
