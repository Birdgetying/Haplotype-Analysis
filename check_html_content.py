import os

html_path = 'result_test/[NoExtend]_HORVU.MOREX.r3.1HG0080220/integrated_analysis.html'
with open(html_path, 'r', encoding='utf-8') as f:
    content = f.read()

print(f'文件大小: {len(content)} bytes ({len(content)//1024} KB)')
print(f'文件总行数: {content.count(chr(10))}')
print()

checks = [
    ('过滤面板CSS', 'filter-panel'),
    ('变异类型勾选框class', 'ann-cb'),
    ('INS勾选', 'value="INS"'),
    ('DEL勾选', 'value="DEL"'),
    ('SV勾选', 'value="SV"'),
    ('Promoter勾选', 'value="promoter"'),
    ('intron勾选', 'value="intron"'),
    ('UTR勾选', 'value="UTR"'),
    ('启动子文字标注', 'Promoter'),
    ('基因ID', '1HG0080220'),
    ('applyFilters函数', 'applyFilters'),
    ('exportSVG按钮', 'exportSVG'),
    ('MAF滑块', 'mafSlider'),
    ('缩放控件', 'zoom-controls'),
    ('LD矩阵', 'LD'),
    ('箱线图', 'boxplot'),
    ('单倍型表格', 'haplotype'),
]

for name, key in checks:
    status = 'OK' if key in content else 'MISS'
    print(f'  [{status}] {name}')

print()
# 检查标题里的基因ID格式
import re
title_match = re.search(r'<title>(.*?)</title>', content)
if title_match:
    print(f'页面标题: {title_match.group(1)}')

h1_match = re.search(r'<h1[^>]*>(.*?)</h1>', content, re.DOTALL)
if h1_match:
    print(f'H1标题: {h1_match.group(1)[:100]}')

# 找 gene_id 或 gene 名称显示
gene_id_matches = re.findall(r'HORVU[^\s<"\']{5,30}', content)
if gene_id_matches:
    print(f'基因ID出现: {list(set(gene_id_matches))[:5]}')
