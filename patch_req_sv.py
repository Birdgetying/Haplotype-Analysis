"""
直接在 required_files 中添加 sv_variants.vcf.gz
"""
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 找到 required_files 区域（必须是数组内容中的）
idx = c.find("required_files = [")
print("required_files at:", idx)

# 在该位置之后的2000字节内查找 variants.vcf.gz
snippet = c[idx:idx+2000]
idx2 = snippet.find("'variants.vcf.gz'  # 关键")
print("In-snippet at:", idx2)

if idx2 >= 0:
    # 在该行后插入 sv_variants.vcf.gz
    target = "'variants.vcf.gz'  # 关键：VCF子集文件\n        ]"
    replacement = "'variants.vcf.gz',  # 关键：VCF子集文件\n            'sv_variants.vcf.gz',  # 结构变异VCF子集\n        ]"
    if target in c:
        c = c.replace(target, replacement, 1)
        with open('genome_wide_haplotype_scan.py', 'w', encoding='utf-8') as f:
            f.write(c)
        print("OK: sv_variants.vcf.gz added to required_files")
    else:
        print("FAIL: exact pattern not found")
        # 打印实际内容
        abs_idx2 = idx + idx2
        print("Actual at", abs_idx2, ":", repr(c[abs_idx2:abs_idx2+100]))
        # 尝试找到 ] 后的位置
        bracket_idx = c.find('\n        ]', abs_idx2)
        print("Next ] at:", bracket_idx)
        print("Between:", repr(c[abs_idx2:bracket_idx]))
else:
    print("FAIL: 'variants.vcf.gz  # 关键' not found in required_files snippet")
    print("Snippet:", repr(snippet[:300]))
