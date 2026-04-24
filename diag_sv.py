"""诊断：检查各模式在文件中的实际文本"""
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 1. cophe3-title
idx = content.find('--cophe3-title')
if idx >= 0:
    print("cophe3-title found at:", idx)
    print(repr(content[idx:idx+200]))
else:
    print("cophe3-title NOT FOUND")

print()

# 2. required_files
idx = content.find("'variants.vcf.gz'")
if idx >= 0:
    print("variants.vcf.gz found at:", idx)
    print(repr(content[idx-50:idx+200]))
else:
    print("variants.vcf.gz NOT FOUND")

print()

# 3. extractor pattern - find the right anchor
idx = content.find("positions, hap_df, hap_sample_df = extractor.extract_region(")
if idx >= 0:
    print("extract_region found at:", idx)
    print(repr(content[idx:idx+500]))
else:
    print("extract_region NOT FOUND")
