with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 找 vcf_mtime 周围内容
idx = c.find('vcf_mtime')
print('=== vcf_mtime context ===')
print('Position:', idx)
print(repr(c[max(0,idx-100):idx+100]))

print()

# 找 variants.vcf.gz 周围内容
idx2 = c.find('variants.vcf.gz')
print('=== variants.vcf.gz context ===')
print('Position:', idx2)
print(repr(c[max(0,idx2-50):idx2+150]))
