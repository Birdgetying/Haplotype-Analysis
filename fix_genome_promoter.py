"""修复 genome_wide_haplotype_scan.py: gene_info_dict 添加 promoter_start/end"""
import py_compile

with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "'promoter_actual_length'" in line and 'vcf_file' not in line:
        indent = '            '
        new_lines = [
            indent + "'promoter_start': (max(1, gene_start - promoter_actual_length) if strand == '+' else gene_end + 1),\n",
            indent + "'promoter_end': (gene_start - 1 if strand == '+' else gene_end + promoter_actual_length),\n",
        ]
        lines = lines[:i+1] + new_lines + lines[i+1:]
        print(f'在行{i+2}后插入promoter_start/end')
        break

with open('genome_wide_haplotype_scan.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)
print('写入完成')

# 验证
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    content = f.read()
if "'promoter_start'" in content and "gene_end + 1" in content:
    print('验证: promoter_start 已添加')
    py_compile.compile('genome_wide_haplotype_scan.py', doraise=True)
    print('语法检查通过')
else:
    print('验证失败')
