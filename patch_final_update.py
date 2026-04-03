#!/usr/bin/env python3
"""补丁脚本：更新现有 integrated_analysis.html 文件
1. 删除 gene-gwas-panel CSS 和 HTML
2. 扩大 network-panel 到 100% 宽度
3. 删除 drawGWASPlot 函数和调用
4. 更新网络图使用容器实际尺寸
"""
import re
import os
import glob

def patch_html(html_path):
    print(f"\n{'='*60}")
    print(f"正在处理: {html_path}")
    
    with open(html_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_len = len(content)
    
    # ========== 1. CSS: 修改 network-panel 宽度 ==========
    content = re.sub(
        r'\.network-panel\s*\{\s*width:\s*350px;\s*min-width:\s*350px;\s*height:\s*280px;',
        '.network-panel { width: 100%; height: 350px;',
        content
    )
    print("  [1] network-panel CSS 已更新为 100% 宽度")
    
    # ========== 2. CSS: 删除 gene-gwas-panel 和 gene-gwas-title 样式 ==========
    content = re.sub(
        r'\.gene-gwas-panel\s*\{[^}]+\}\s*',
        '',
        content
    )
    content = re.sub(
        r'\.gene-gwas-title\s*\{[^}]+\}\s*',
        '',
        content
    )
    print("  [2] gene-gwas-panel/title CSS 已删除")
    
    # ========== 3. HTML: 删除 gene-gwas-panel div ==========
    # 匹配 <div class="gene-gwas-panel">...</div> 整个块
    content = re.sub(
        r'\s*<div class="gene-gwas-panel">.*?</div>\s*</div>',
        '',
        content,
        flags=re.DOTALL,
        count=1
    )
    print("  [3] gene-gwas-panel HTML 已删除")
    
    # ========== 4. JS: 网络图使用容器实际尺寸 ==========
    content = re.sub(
        r'var W = 350,\s*H = 280;',
        'var W = container.clientWidth || 800;\n    var H = container.clientHeight || 350;',
        content
    )
    print("  [4] 网络图尺寸已改为容器实际尺寸")
    
    # ========== 5. JS: 删除 drawGWASPlot(filtered) 调用 ==========
    content = content.replace('drawGWASPlot(filtered);\n', '')
    content = content.replace('drawGWASPlot(filtered);', '')
    print("  [5] drawGWASPlot(filtered) 调用已删除")
    
    # ========== 6. JS: 删除整个 drawGWASPlot 函数 ==========
    content = re.sub(
        r'// =+ 基因结构图[^=]*=+\s*\nfunction drawGWASPlot\(data\)\s*\{[\s\S]*?\n\}\s*\n(?=// =+ 页面初始化)',
        '',
        content
    )
    print("  [6] drawGWASPlot 函数已删除")
    
    # ========== 7. JS: 删除 drawGWASPlot(gwasData) 初始化调用 ==========
    content = content.replace('    drawGWASPlot(gwasData);\n', '')
    content = content.replace('drawGWASPlot(gwasData);\n', '')
    print("  [7] drawGWASPlot(gwasData) 初始化调用已删除")
    
    # ========== 写入文件 ==========
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    new_len = len(content)
    print(f"\n  文件大小: {original_len:,} → {new_len:,} ({new_len - original_len:+,} 字符)")
    
    # ========== 验证 ==========
    checks = {
        'drawNetworkPlot': 'drawNetworkPlot' in content,
        'drawGWASPlot 已删除': 'drawGWASPlot' not in content,
        'gene-gwas-panel 已删除': 'gene-gwas-panel' not in content,
        'network-panel width:100%': 'width: 100%' in content,
    }
    print("\n  验证结果:")
    all_ok = True
    for name, ok in checks.items():
        status = '✓' if ok else '✗'
        print(f"    {status} {name}")
        if not ok:
            all_ok = False
    
    if all_ok:
        print("\n  ✓ 所有检查通过!")
    else:
        print("\n  ⚠ 部分检查未通过，请手动确认")
    
    return all_ok


# 查找所有 integrated_analysis.html 文件
html_files = glob.glob(os.path.join('results', '**', 'integrated_analysis.html'), recursive=True)

# 过滤：只处理包含 drawGWASPlot 或 network-panel 的文件
target_files = []
for f in html_files:
    with open(f, 'r', encoding='utf-8') as fh:
        text = fh.read()
    if 'drawGWASPlot' in text or 'gene-gwas-panel' in text:
        target_files.append(f)

if not target_files:
    print("没有找到需要更新的 HTML 文件")
else:
    print(f"找到 {len(target_files)} 个需要更新的 HTML 文件:")
    for f in target_files:
        print(f"  - {f}")
    
    results = []
    for f in target_files:
        ok = patch_html(f)
        results.append((f, ok))
    
    print(f"\n{'='*60}")
    print("处理完成:")
    for f, ok in results:
        print(f"  {'✓' if ok else '✗'} {f}")
