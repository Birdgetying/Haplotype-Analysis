#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检查 Excel 文件的 sheet 名称和基因列表
"""

import pandas as pd

# Excel 文件路径
excel_file = r'd:\Desktop\project1\database\大麦备份\大麦GWAS及SV筛选基因-2026.03.22(1).xlsx'

# 查看所有 sheet 名称
xls = pd.ExcelFile(excel_file)
print("Excel 文件中的 sheet 名称:")
for name in xls.sheet_names:
    print(f"  - '{name}'")

# 读取所有 sheet
for sheet_name in xls.sheet_names:
    print(f"\n{'=' * 80}")
    print(f"Sheet: {sheet_name}")
    print(f"{'=' * 80}")
    
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    print(f"行数: {len(df)}")
    print(f"列名: {list(df.columns)}")
    print(f"\n前3行数据:")
    print(df.head(3))
