#!/usr/bin/env python3
"""补丁脚本：将现有HTML的过滤功能从包含模式改为排除模式"""
import re
import glob
import os

def patch_filter(html_path):
    print(f"\n{'='*60}")
    print(f"处理: {html_path}")
    
    with open(html_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_len = len(content)
    modified = False
    
    # 1. 替换Annotation下拉框为多选框
    old_select = r'<select id="annotationFilter"[^>]*>[\s\S]*?</select>'
    new_checkboxes = '''<div class="annotation-checkboxes" id="annotationExclude">
                <label class="ann-checkbox"><input type="checkbox" value="missense" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#e74c3c"></span>Missense</label>
                <label class="ann-checkbox"><input type="checkbox" value="synonymous" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#f39c12"></span>Synonymous</label>
                <label class="ann-checkbox"><input type="checkbox" value="UTR" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#9b59b6"></span>UTR</label>
                <label class="ann-checkbox"><input type="checkbox" value="indel" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#3498db"></span>Indel</label>
                <label class="ann-checkbox"><input type="checkbox" value="SV" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#e67e22"></span>SV</label>
                <label class="ann-checkbox"><input type="checkbox" value="promoter" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#27ae60"></span>Promoter</label>
                <label class="ann-checkbox"><input type="checkbox" value="intron" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#95a5a6"></span>Intron</label>
                <label class="ann-checkbox"><input type="checkbox" value="other" onchange="updateAnnotationExclude()"> <span class="ann-color" style="background:#bdc3c7"></span>Other</label>
            </div>'''
    
    if re.search(old_select, content):
        content = re.sub(old_select, new_checkboxes, content)
        print("  [1] Annotation下拉框已替换为多选框")
        modified = True
    
    # 2. 更新filter-group标签
    content = re.sub(
        r'<div class="filter-group">\s*<label>Annotation:</label>',
        '<div class="filter-group annotation-exclude">\n            <label>Exclude Types:</label>',
        content
    )
    
    # 3. 添加CSS样式（如果不存在）
    if '.annotation-exclude' not in content:
        css_addition = '''        /* 排除模式多选框样式 */
        .annotation-exclude { align-items: flex-start; }
        .annotation-checkboxes { display: flex; flex-wrap: wrap; gap: 6px 12px; max-width: 400px; }
        .ann-checkbox { display: flex; align-items: center; gap: 4px; font-size: 11px; cursor: pointer;
                        padding: 2px 6px; background: #fff; border: 1px solid #e0e0e0; border-radius: 3px; }
        .ann-checkbox:hover { background: #f0f0f0; }
        .ann-checkbox input[type="checkbox"] { margin: 0; cursor: pointer; }
        .ann-color { display: inline-block; width: 10px; height: 10px; border-radius: 2px; }
'''
        # 插入到 </style> 之前
        content = re.sub(r'(</style>)', css_addition + r'\1', content, count=1)
        print("  [2] CSS样式已添加")
        modified = True
    
    # 4. 更新JS过滤逻辑
    # 替换 currentFilter 初始化
    content = re.sub(
        r"currentFilter\s*=\s*\{\s*maf:\s*0,\s*missingRate:\s*1\.0,\s*annotation:\s*'all'\s*\}",
        "currentFilter = { maf: 0, missingRate: 1.0, excludedAnnotations: [] }",
        content
    )
    
    # 删除 annotation 相关的 updateFilter 逻辑
    content = re.sub(
        r"\} else if \(type === 'annotation'\) \{\s*currentFilter\.annotation = value;\s*\}",
        "}",
        content
    )
    
    # 添加 updateAnnotationExclude 函数
    if 'function updateAnnotationExclude' not in content:
        new_func = '''function updateAnnotationExclude() {
    var checkboxes = document.querySelectorAll('#annotationExclude input[type="checkbox"]');
    currentFilter.excludedAnnotations = [];
    checkboxes.forEach(function(cb) {
        if (cb.checked) {
            currentFilter.excludedAnnotations.push(cb.value);
        }
    });
}
'''
        # 在 resetFilters 之前插入
        content = re.sub(
            r'(function resetFilters\(\))',
            new_func + r'\1',
            content
        )
        print("  [3] updateAnnotationExclude函数已添加")
        modified = True
    
    # 5. 更新 resetFilters 函数
    content = re.sub(
        r"document\.getElementById\('annotationFilter'\)\.value = 'all';",
        "// 清除所有排除勾选\n    document.querySelectorAll('#annotationExclude input[type=\"checkbox\"]').forEach(function(cb) {\n        cb.checked = false;\n    });",
        content
    )
    content = re.sub(
        r"currentFilter = \{ maf: 0, missingRate: 1\.0, annotation: 'all' \}",
        "currentFilter = { maf: 0, missingRate: 1.0, excludedAnnotations: [] }",
        content
    )
    
    # 6. 更新 applyFilters 函数
    old_filter = r"&& \(currentFilter\.annotation === 'all' \|\| d\.annotation === currentFilter\.annotation\)"
    new_filter = "&& currentFilter.excludedAnnotations.indexOf(d.annotation || 'other') === -1"
    content = re.sub(old_filter, new_filter, content)
    
    # 添加 updateAnnotationExclude() 调用到 applyFilters 开头
    if 'updateAnnotationExclude()' not in content:
        content = re.sub(
            r'(function applyFilters\(\) \{)',
            r'\1\n    updateAnnotationExclude();',
            content
        )
        print("  [4] applyFilters函数已更新为排除模式")
        modified = True
    
    # 写入文件
    if modified:
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(content)
        new_len = len(content)
        print(f"\n  文件大小: {original_len:,} → {new_len:,} ({new_len - original_len:+,} 字符)")
        print("  ✓ 文件已更新")
    else:
        print("  - 无需修改")
    
    return modified

# 查找所有 integrated_analysis.html
html_files = glob.glob(os.path.join('results', '**', 'integrated_analysis.html'), recursive=True)

if not html_files:
    print("未找到 HTML 文件")
else:
    print(f"找到 {len(html_files)} 个 HTML 文件")
    results = []
    for f in html_files:
        ok = patch_filter(f)
        results.append((f, ok))
    
    print(f"\n{'='*60}")
    print("处理完成:")
    for f, ok in results:
        print(f"  {'✓' if ok else '-'} {f}")
