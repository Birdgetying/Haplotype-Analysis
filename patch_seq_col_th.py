"""修复序列列th添加seq-col-th class"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

old = '''        for pos in display_positions:
            # 物理坐标竖排，千分位逗号分隔，宽度与序列列 td 严格一致
            pos_str = f'{pos:,}'  # 千分位逗号
            html += (f'<th style="width:20px;min-width:20px;max-width:20px;padding:0;'
                     f'vertical-align:top;overflow:hidden;">'
                     f'<div style="writing-mode:vertical-rl;transform:rotate(180deg);'
                     f'width:20px;height:60px;display:flex;align-items:center;justify-content:center;'
                     f'font-size:9px;color:#f5f5f5;background:#2c3e50;'
                     f'font-weight:600;letter-spacing:0;box-sizing:border-box;">{pos_str}</div></th>\\n')'''

new = '''        for pos in display_positions:
            # 物理坐标竖排，千分位逗号分隔，宽度与序列列 td 严格一致
            pos_str = f'{pos:,}'  # 千分位逗号
            html += (f'<th class="seq-col-th" data-pos="{pos}" style="width:20px;min-width:20px;max-width:20px;padding:0;'
                     f'vertical-align:top;overflow:hidden;">'
                     f'<div style="writing-mode:vertical-rl;transform:rotate(180deg);'
                     f'width:20px;height:60px;display:flex;align-items:center;justify-content:center;'
                     f'font-size:9px;color:#f5f5f5;background:#2c3e50;'
                     f'font-weight:600;letter-spacing:0;box-sizing:border-box;">{pos_str}</div></th>\\n')'''

if old in content:
    content = content.replace(old, new, 1)
    with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
        f.write(content)
    print("SUCCESS: seq-col-th class added to sequence column th")
else:
    print("FAIL: old string not found")
    # 找一下实际文本
    idx = content.find('for pos in display_positions')
    if idx >= 0:
        print("Found 'for pos in display_positions' at:", idx)
        print(repr(content[idx:idx+400]))
