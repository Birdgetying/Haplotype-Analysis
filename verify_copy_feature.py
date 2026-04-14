#!/usr/bin/env python3
"""验证所有基因的HTML文件是否包含复制样本功能"""
import os

results_dir = 'results2'
genes = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]

print(f"总计 {len(genes)} 个基因\n")
print("="*70)

count_with_copy = 0
for i, g in enumerate(sorted(genes), 1):
    html = f'{results_dir}/{g}/integrated_analysis.html'
    if os.path.exists(html):
        content = open(html, 'r', encoding='utf-8').read()
        has_copy = 'function copySamples' in content
        if has_copy:
            count_with_copy += 1
            status = '✓ 有复制功能'
        else:
            status = '✗ 无复制功能'
    else:
        status = '✗ HTML不存在'
    
    print(f'{i:2d}. {g:40s} {status}')

print("="*70)
print(f"\n统计: {count_with_copy}/{len(genes)} 个基因已添加复制样本功能")
