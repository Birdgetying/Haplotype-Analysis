# -*- coding: utf-8 -*-
"""修复 annotate_snp_effects_for_region 函数签名"""
import re

with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 旧的函数签名模式
old_pattern = r'(def annotate_snp_effects_for_region\(vcf_file: str, fasta_path: str, gene_chrom: str,\n\s+cds_intervals: list, exon_intervals: list,\n\s+gene_strand: str, positions: list\) -> dict:)'

# 新的函数签名
new_signature = '''def annotate_snp_effects_for_region(vcf_file: str, fasta_path: str, gene_chrom: str,
                                     cds_intervals: list, exon_intervals: list,
                                     gene_strand: str, positions: list,
                                     gene_start: int = None, gene_end: int = None,
                                     promoter_start: int = None, promoter_end: int = None) -> dict:'''

# 检查旧模式是否存在
if re.search(old_pattern, content):
    print(f'[INFO] 找到旧的函数签名，准备替换...')
    new_content = re.sub(old_pattern, new_signature, content)
    
    with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print('[OK] 函数签名已更新')
else:
    print('[WARNING] 未找到旧的函数签名，尝试其他模式...')
    
    # 检查当前函数签名
    match = re.search(r'def annotate_snp_effects_for_region\([^)]+\) -> dict:', content)
    if match:
        current_sig = match.group(0)
        print(f'当前签名: {current_sig[:100]}...')
        
        # 检查是否已包含新参数
        if 'gene_start' in current_sig:
            print('[OK] 函数签名已包含 gene_start 参数')
        else:
            print('[ERROR] 函数签名需要修复')
    else:
        print('[ERROR] 找不到函数定义')