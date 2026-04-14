#!/usr/bin/env python3
"""修复：Copy 按钮移到序列右边，Copy All 在表头与 n 同行"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

file_path = 'haplotype_phenotype_analysis.py'

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# 1. 修改表头：删除第4列的 Copy All，移到最后（n 列前面）
old_thead = '''        html += f\'\'\'<table class="data-table" style="width:{table_width}px;">
<thead><tr style="height:45px;">
    <th style="width:90px;min-width:90px;max-width:90px;text-align:left;padding-left:10px;vertical-align:middle;">Haplotype</th>
    <th class="effect-cell" style="vertical-align:middle;">Effect (vs Grand Mean)</th>
    <th class="box-cell" style="vertical-align:middle;">Phenotype</th>
    <th style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;">
        <button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:12px;font-weight:700;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy All</button>
    </th>\\n\'\'\''''

new_thead = '''        html += f\'\'\'<table class="data-table" style="width:{table_width}px;">
<thead><tr style="height:45px;">
    <th style="width:90px;min-width:90px;max-width:90px;text-align:left;padding-left:10px;vertical-align:middle;">Haplotype</th>
    <th class="effect-cell" style="vertical-align:middle;">Effect (vs Grand Mean)</th>
    <th class="box-cell" style="vertical-align:middle;">Phenotype</th>\\n\'\'\''''

content = content.replace(old_thead, new_thead)

# 2. 在序列列之后、n 列之前添加 Copy All 表头
old_n_header = '''        html += '<th style="width:40px;vertical-align:middle;">n</th></tr></thead><tbody>\\n\''''

new_n_header = '''        html += f'<th style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;"><button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:12px;font-weight:700;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy All</button></th>\\n'
        html += '<th style="width:40px;vertical-align:middle;">n</th></tr></thead><tbody>\\n\''''

content = content.replace(old_n_header, new_n_header)

# 3. 删除原来的 Copy 按钮列（在数据行开头的位置）
old_copy_in_rows = '''            # 添加复制样本按钮列（在 n 列之前）
            samples_for_hap = hap_samples_map.get(hap, [])
            samples_str = ','.join(samples_for_hap) if samples_for_hap else ''
            html += f'<td style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;"><button class="copy-samples-btn" onclick="copySamples(\\'{hap}\\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:5px 12px;border-radius:4px;cursor:pointer;font-size:11px;font-weight:600;white-space:nowrap;transition:all 0.2s;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy</button></td>\\n'
            
            html += f'<td class="n-cell" style="padding-left:15px;">{cnt}</td></tr>\\n\''''

new_copy_in_rows = '''            # 添加复制样本按钮列（在序列列之后，n 列之前）
            samples_for_hap = hap_samples_map.get(hap, [])
            samples_str = ','.join(samples_for_hap) if samples_for_hap else ''
            html += f'<td style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;"><button class="copy-samples-btn" onclick="copySamples(\\'{hap}\\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:5px 12px;border-radius:4px;cursor:pointer;font-size:11px;font-weight:600;white-space:nowrap;transition:all 0.2s;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy</button></td>\\n'
            
            html += f'<td class="n-cell">{cnt}</td></tr>\\n\''''

content = content.replace(old_copy_in_rows, new_copy_in_rows)

# 4. 修改底部坐标轴行：Copy 列和 n 列分开
old_bottom = '''        # 第4列：复制样本列（空）
        html += '<td style="width:80px;min-width:80px;max-width:80px;border:none;"></td>\\n'
        
        # 第5列：n 列（空）
        html += '<td style="width:40px;min-width:40px;max-width:40px;border:none;"></td>\\n'
        
        # 剩余列（序列）
        html += f'<td colspan="{len(display_positions)}" style="border:none;"></td>\\n\''''

new_bottom = '''        # 序列列（空）
        html += f'<td colspan="{len(display_positions)}" style="border:none;"></td>\\n'
        
        # Copy 列（空）
        html += '<td style="width:80px;min-width:80px;max-width:80px;border:none;"></td>\\n'
        
        # n 列（空）
        html += '<td style="width:40px;min-width:40px;max-width:40px;border:none;"></td>\\n\''''

content = content.replace(old_bottom, new_bottom)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)

print('✓ 修复完成')
print('  列顺序: Haplotype | Effect | Phenotype | 序列列... | Copy All/Copy | n')
print('  Copy All 在表头，与 n 列相邻')
print('  Copy 按钮在每行的序列列右边、n 列左边')
