# -*- coding: utf-8 -*-
"""修复 CDS 区域 indel 被错误归类为 missense 的问题"""
import re

target_file = 'd:/Desktop/project1/haplotype_phenotype_analysis.py'

with open(target_file, 'r', encoding='utf-8') as f:
    content = f.read()

# 找到问题代码块
old_block = """            # 构建功能注释信息（用于indel分类）
            functional_ann = ann
            # 如果是indel，尝试获取更详细的功能注释
            if ann in ('indel', 'INS', 'DEL'):
                # 检查位置是否在CDS区域
                in_cds_check = any(cs <= ip <= ce for cs, ce in cds) if cds else False
                in_exon_check = any(es <= ip <= ee for es, ee in exons) if exons else False
                if in_cds_check:
                    functional_ann = 'missense'  # 在CDS区域的indel归类为错义
                elif in_exon_check:
                    functional_ann = 'UTR'
                elif g_start and g_end and g_start <= ip <= g_end:
                    functional_ann = 'intron'
                else:
                    functional_ann = 'other'
            
            gwas_data.append({"""

# 替换为正确逻辑：indel 保持原始类型
new_block = """            # 构建功能注释信息
            # **修复**：indel 应该保持原始类型（INS/DEL），只有 SNP 才会被分为 missense/synonymous
            functional_ann = ann
            
            gwas_data.append({"""

if old_block in content:
    content = content.replace(old_block, new_block)
    with open(target_file, 'w', encoding='utf-8') as f:
        f.write(content)
    print("[OK] 已修复 CDS 区域 indel 归类问题")
else:
    print("[ERROR] 未找到目标代码块")
    # 尝试找到附近的代码
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if 'functional_ann = ' in line and 'missense' not in line:
            print(f"行{i+1}: {line}")