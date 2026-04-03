#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""验证exportSVG功能是否正确生成"""

import haplotype_phenotype_analysis as hpa
import pandas as pd

# 创建最小测试数据 - 使用正确的列名
hap_df = pd.DataFrame({
    'Hap_Name': ['H1', 'H2'],
    'Sample': ['S1', 'S2'],
    'phenotype': [1.0, 2.0]
})

# 创建ReportGenerator实例
generator = hpa.ReportGenerator()

# 生成HTML - 注意参数顺序
html = generator.generate_integrated_html(
    hap_sample_df=hap_df,
    effect_results={},
    variant_positions=[],
    region_start=1,
    region_end=2000,
    phenotype_col='TEST'
)

# 验证内容
print("=" * 50)
print("验证结果:")
print("=" * 50)
print(f"1. exportSVG函数: {'✓ 存在' if 'function exportSVG()' in html else '✗ 不存在'}")
print(f"2. SVG按钮: {'✓ 存在' if 'SVG</button>' in html else '✗ 不存在'}")
print(f"3. Print按钮: {'✓ 存在' if 'Print</button>' in html else '✗ 不存在'}")
print(f"4. Export标签: {'✓ 存在' if 'Export:' in html else '✗ 不存在'}")
print("=" * 50)

# 保存测试文件
with open('test_export_verify.html', 'w', encoding='utf-8') as f:
    f.write(html)
print("测试文件已保存: test_export_verify.html")
print("=" * 50)
