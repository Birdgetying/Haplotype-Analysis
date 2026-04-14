#!/usr/bin/env python3
"""
测试复制样本功能
从现有 HTML 提取数据，生成新的 HTML 来验证功能
"""

import os
import json
import pandas as pd

# 配置
GENE_ID = 'HORVU.MOREX.r3.4HG0387570'
RESULTS_DIR = f'results2/{GENE_ID}'
HTML_FILE = f'{RESULTS_DIR}/integrated_analysis.html'

# 读取旧 HTML 提取关键信息
print(f"正在处理: {GENE_ID}")
print(f"HTML 文件: {HTML_FILE}")

# 检查文件是否存在
if not os.path.exists(HTML_FILE):
    print(f"错误: 文件不存在")
    exit(1)

# 读取 HTML 文件
with open(HTML_FILE, 'r', encoding='utf-8') as f:
    html_content = f.read()

print(f"✓ HTML 文件读取成功 ({len(html_content)} 字符)")

# 检查是否已经有复制功能
if 'copySamples' in html_content:
    print("✓ HTML 文件中已有复制功能")
else:
    print("✗ HTML 文件中没有复制功能")
    
# 提取关键数据用于验证
import re

# 提取基因信息
gene_match = re.search(r'var geneLabelText = "([^"]+)";', html_content)
if gene_match:
    print(f"\n基因ID: {gene_match.group(1)}")

# 提取单倍型信息
hap_match = re.findall(r'data-hap="([^"]+)"', html_content)
if hap_match:
    print(f"\n单倍型数量: {len(set(hap_match))}")
    print(f"单倍型列表: {sorted(set(hap_match))[:10]}...")

# 提取变异位置
pos_match = re.findall(r'<th[^>]*>[^<]*?<div[^>]*>([^<]+)</div>', html_content)
if pos_match:
    print(f"\n变异位点数量: {len(pos_match)}")
    print(f"前5个位置: {pos_match[:5]}")

# 统计表格行数
row_match = re.findall(r'<tr class="[^"]*" data-hap="[^"]+">', html_content)
if row_match:
    print(f"\n表格行数: {len(row_match)}")

print("\n" + "="*70)
print("测试总结:")
print("="*70)
print("由于 database 目录为空，无法重新生成 HTML")
print("新数据库建立完成后，请使用以下命令重新生成:")
print(f"  python feasible_test_database_analysis.py")
print("\n新生成的 HTML 将包含:")
print("  1. 复制样本按钮列")
print("  2. Copy All 按钮（表头）")
print("  3. 每行的 Copy 按钮")
print("  4. JavaScript 复制功能")
