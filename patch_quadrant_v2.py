#!/usr/bin/env python3
"""
补丁脚本 v2：精确更新HTML为四象限布局
"""
import os
import re
import glob
import json

def patch_html_v2(filepath):
    """精确更新HTML文件"""
    print(f"处理: {filepath}")
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    changes = []
    
    # 1. 在quadrant-top-right内，SVG之前添加P值散点图容器
    if 'id="pvalue-scatter-plot"' not in content:
        # 查找quadrant-top-right开始后的位置
        pattern = r'(<div class="quadrant-top-right">[\s\S]*?<div class="gene-structure-area">)\s*(<svg)'
        if re.search(pattern, content):
            replacement = r'''\1
                    <div class="pvalue-scatter" id="pvalue-scatter-plot"></div>
                \2'''
            content = re.sub(pattern, replacement, content)
            changes.append("[1] P值散点图容器已添加")
        else:
            # 尝试另一种模式
            pattern2 = r'(<div class="quadrant-top-right">\s*<div class="gene-structure-area">)\s*(<svg)'
            if re.search(pattern2, content):
                replacement2 = r'''<div class="quadrant-top-right">
                <div class="gene-pvalue-container">
                    <div class="pvalue-scatter" id="pvalue-scatter-plot"></div>
                    <div class="gene-structure-area">
                \2'''
                content = re.sub(pattern2, replacement2, content)
                changes.append("[1] P值散点图容器和包装层已添加")
    
    # 2. 检查并添加左下和右下象限
    if 'class="quadrant-bottom-left"' not in content:
        # 找到表格结束位置
        table_end_pattern = r'(</tbody></table>\s*)(</div><!-- gene-structure-area -->|</div><!-- main-data-section -->)'
        
        bottom_sections = '''
            </div><!-- gene-structure-area -->
                </div><!-- gene-pvalue-container -->
            </div><!-- quadrant-top-right -->
            
            <!-- 左下：效应图摘要 -->
            <div class="quadrant-bottom-left">
                <div class="effect-box-panel">
                    <div class="effect-box-panel-title">Effect Size Summary</div>
                    <div id="effect-summary-chart" style="width:100%;height:200px;"></div>
                </div>
            </div>
            
            <!-- 右下：图例和信息 -->
            <div class="quadrant-bottom-right">
                <div class="sequence-panel">
                    <div class="sequence-panel-title">Legend & Information</div>
                    <div style="font-size:11px;color:#666;padding:10px;">
                        <p style="margin-bottom:8px;"><strong>Variant Types:</strong></p>
                        <div style="display:flex;flex-wrap:wrap;gap:8px;margin-bottom:10px;">
                            <span><span style="display:inline-block;width:12px;height:12px;background:#e74c3c;border-radius:50%;"></span> Missense</span>
                            <span><span style="display:inline-block;width:12px;height:12px;background:#3498db;border-radius:50%;"></span> Synonymous</span>
                            <span><span style="display:inline-block;width:12px;height:12px;background:#2ecc71;border-radius:50%;"></span> UTR</span>
                            <span><span style="display:inline-block;width:12px;height:12px;background:#9b59b6;border-radius:50%;"></span> InDel</span>
                            <span><span style="display:inline-block;width:12px;height:12px;background:#f39c12;border-radius:50%;"></span> SV</span>
                            <span><span style="display:inline-block;width:12px;height:12px;background:#1abc9c;border-radius:50%;"></span> Promoter</span>
                        </div>
                        <p><em>Hover over scatter points for details. Red dashed line indicates p=0.05 significance threshold.</em></p>
                    </div>
                </div>
            </div>
        </div><!-- quadrant-layout -->
    </div>
    </div>'''
        
        if re.search(table_end_pattern, content):
            content = re.sub(table_end_pattern, r'\1' + bottom_sections, content)
            changes.append("[2] 左下和右下象限已添加")
        else:
            # 尝试找到integrated-view或quadrant-layout的结束
            alt_pattern = r'(</tbody></table>[\s\S]*?)(</div><!-- integrated-view -->|</div><!-- quadrant-layout -->)'
            if re.search(alt_pattern, content):
                content = re.sub(alt_pattern, 
                    '</tbody></table>\n                    </div><!-- gene-structure-area -->\n                </div><!-- gene-pvalue-container -->' + bottom_sections.replace('</div><!-- quadrant-layout -->\n    </div>\n    </div>', ''),
                    content, count=1)
                changes.append("[2] 左下和右下象限已添加(alt)")
    
    # 3. 确保CSS中有gene-pvalue-container样式
    if '.gene-pvalue-container' not in content:
        css_add = '''
        .gene-pvalue-container { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .pvalue-scatter { width: 100%; height: 150px; border-bottom: 1px solid #eee; margin-bottom: 5px; }'''
        content = re.sub(r'(\.quadrant-top-right\s*\{[^}]*\})', r'\1' + css_add, content)
        changes.append("[3] gene-pvalue-container CSS已添加")
    
    # 保存
    if content != original_content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        for c in changes:
            print(f"  {c}")
        return True
    else:
        print("  无需更新")
        return False

def main():
    results_dir = os.path.join(os.path.dirname(__file__), 'results')
    html_files = glob.glob(os.path.join(results_dir, '**/integrated_analysis.html'), recursive=True)
    
    print(f"找到 {len(html_files)} 个HTML文件\n")
    
    updated = 0
    for f in html_files:
        if patch_html_v2(f):
            updated += 1
        print()
    
    print(f"[完成] 更新 {updated}/{len(html_files)} 个文件")

if __name__ == '__main__':
    main()
