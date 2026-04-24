"""精确检查 required_files 的内容"""
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 找到 required_files 开始
idx = c.find("required_files = [")
end_idx = idx + c[idx:].find("        if os.path.exists(gene_data_dir):")

req_content = c[idx:end_idx]
print("Required_files content:")
print("=" * 60)
for i, line in enumerate(req_content.split('\n')):
    if line.strip():
        # 计算空格数
        spaces = len(line) - len(line.lstrip())
        print(f"Line {i}: spaces={spaces:2d} | {repr(line)}")
print("=" * 60)

# 现在找到正确的插入位置
# 在 variants.vcf.gz 行之后添加 sv_variants.vcf.gz
# variants.vcf.gz 行的格式: "            'variants.vcf.gz'  # 关键：VCF子集文件"
lines = req_content.split('\n')
vcf_gz_idx = -1
for i, line in enumerate(lines):
    if "'variants.vcf.gz'" in line and "# 关键" in line:
        vcf_gz_idx = i
        print(f"Found variants.vcf.gz at line {i}: {repr(line)}")
        print(f"Spaces before: {len(line) - len(line.lstrip())}")
        break

if vcf_gz_idx >= 0:
    # 下一行应该是 ] - 但需要确认没有 sv_variants.vcf.gz
    if vcf_gz_idx + 1 < len(lines):
        next_line = lines[vcf_gz_idx + 1]
        print(f"Next line: {repr(next_line)}")
        if "'sv_variants.vcf.gz'" in next_line:
            print("sv_variants.vcf.gz already present!")
        else:
            print("Need to add sv_variants.vcf.gz before the ]")
            # 计算需要多少空格（和variants.vcf.gz行一样）
            vcf_gz_spaces = len(lines[vcf_gz_idx]) - len(lines[vcf_gz_idx].lstrip())
            indent = ' ' * vcf_gz_spaces
            new_line = indent + "'sv_variants.vcf.gz',  # 结构变异VCF子集"
            # 替换
            old_text = lines[vcf_gz_idx] + '\n' + next_line
            new_text = lines[vcf_gz_idx] + ',\n' + new_line + '\n' + next_line
            full_old = req_content.replace(old_text, new_text, 1)
            if old_text in req_content:
                c = c[:idx] + full_old + c[idx+len(req_content):]
                with open('genome_wide_haplotype_scan.py', 'w', encoding='utf-8') as f:
                    f.write(c)
                print("SUCCESS: sv_variants.vcf.gz added!")
            else:
                print("FAIL: couldn't find exact text to replace")
                print("Trying old_text:", repr(old_text[:100]))
    else:
        print("No next line")
