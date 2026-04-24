"""
修复JS端的序列列识别（用seq-col-th class代替idx>=3）和colorbar颜色
"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 1. 修复JS中的序列列识别 - 旧代码
old_js_filter = '''    var allThs = Array.from(table.querySelectorAll('thead th'));
    var visibleVarThs = allThs.filter(function(th, idx) {{
        return idx >= 3 && idx < allThs.length - 1 && th.style.display !== 'none';
    }});'''

new_js_filter = '''    // 用 seq-col-th class 识别序列列，而非硬编码 idx>=3
    // 这样即使前面有更多协变量列也能正确识别
    var allThs = Array.from(table.querySelectorAll('thead th'));
    var visibleVarThs = allThs.filter(function(th) {{
        return th.classList.contains('seq-col-th') && th.style.display !== 'none';
    }});'''

if old_js_filter in content:
    content = content.replace(old_js_filter, new_js_filter, 1)
    print("成功修复JS序列列识别")
else:
    print("未找到JS序列列识别旧代码")
    idx = content.find('idx >= 3 && idx < allThs.length - 1')
    print(f"idx >= 3 在字符 {idx}")

# 2. 修复colorbar颜色 - 旧代码（蓝->红）
old_colorbar = '''<div id="ld-colorbar" style="display:flex;align-items:center;margin-top:4px;padding-left:0px;">
        <span style="font-size:9px;color:#555;margin-right:4px;">r²:</span>
        <div style="background:linear-gradient(to right,#313695,#4575b4,#74add1,#abd9e9,#e0f3f8,#fee090,#fdae61,#f46d43,#d73027,#a50026);width:80px;height:10px;border-radius:2px;"></div>
        <span style="font-size:9px;color:#555;margin-left:4px;">0</span>
        <span style="font-size:9px;color:#555;margin-left:2px;">→</span>
        <span style="font-size:9px;color:#555;margin-left:2px;">1</span>'''

new_colorbar = '''<div id="ld-colorbar" style="display:flex;align-items:center;margin-top:4px;padding-left:0px;">
        <span style="font-size:9px;color:#555;margin-right:4px;">r²:</span>
        <div style="background:linear-gradient(to right,#ffffff,#ffffcc,#fdcc8a,#fc8d59,#e34a33,#b30000);width:80px;height:10px;border-radius:2px;border:1px solid #ccc;"></div>
        <span style="font-size:9px;color:#555;margin-left:4px;">0</span>
        <span style="font-size:9px;color:#555;margin-left:2px;">→</span>
        <span style="font-size:9px;color:#555;margin-left:2px;">1</span>'''

if old_colorbar in content:
    content = content.replace(old_colorbar, new_colorbar, 1)
    print("成功修复colorbar颜色")
else:
    print("未找到colorbar旧代码")
    idx2 = content.find('313695')
    print(f"313695 在字符 {idx2}")

with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("文件保存完成")

# 验证
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    verify = f.read()
print("验证 seq-col-th:", 'seq-col-th' in verify)
print("验证 idx>=3已删除:", 'idx >= 3 && idx < allThs.length - 1' not in verify)
print("验证新colorbar:", 'ffffff,#ffffcc,#fdcc8a' in verify)
print("验证旧colorbar已删除:", '313695' not in verify)

