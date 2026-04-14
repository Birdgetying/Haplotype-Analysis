#!/usr/bin/env python3
"""修改 haplotype_phenotype_analysis.py 的表格排版"""

file_path = 'haplotype_phenotype_analysis.py'

with open(file_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# 修改表头（删除 height:60px，调整 Copy All 按钮）
lines[5172] = lines[5172].replace('height:60px;', '')
lines[5173] = lines[5173].replace('height:60px;', '')
lines[5174] = lines[5174].replace('height:60px;', '')
lines[5175] = lines[5175].replace('width:80px', 'width:70px').replace('height:60px;', '')
lines[5176] = lines[5176].replace('padding:4px 8px', 'padding:5px 10px').replace('font-size:10px', 'font-size:11px;font-weight:600')

# 修改 n 列表头
lines[5188] = lines[5188].replace('height:60px;', '')

# 修改 Copy 按钮（在 n 列之前）
lines[5312] = lines[5312].replace('width:80px', 'width:70px').replace('padding:4px 8px', 'padding:4px 10px').replace('font-size:9px', 'font-size:10px;font-weight:500')

# 修改 JavaScript 函数，移除"已复制"提示
# copySamples 函数
for i in range(5430, 5470):
    if i < len(lines):
        if 'btn.innerText' in lines[i] or 'btn.style.background' in lines[i] or 'setTimeout' in lines[i]:
            lines[i] = lines[i].replace('btn.innerText', '// btn.innerText').replace('btn.style.background', '// btn.style.background').replace('setTimeout', '// setTimeout')

# copyAllSamples 函数
for i in range(5470, 5515):
    if i < len(lines):
        if 'btn.innerText' in lines[i] or 'btn.style.background' in lines[i] or 'setTimeout' in lines[i]:
            lines[i] = lines[i].replace('btn.innerText', '// btn.innerText').replace('btn.style.background', '// btn.style.background').replace('setTimeout', '// setTimeout')

with open(file_path, 'w', encoding='utf-8') as f:
    f.writelines(lines)

print('✓ 修改完成')
print('  - 删除表头 height:60px')
print('  - Copy All 按钮: 70px 宽，11px 字体，加粗')
print('  - Copy 按钮: 70px 宽，10px 字体，加粗')
print('  - JavaScript: 注释掉"已复制"提示')
