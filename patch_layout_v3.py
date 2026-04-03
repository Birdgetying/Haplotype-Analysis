#!/usr/bin/env python3
"""
补丁脚本v3：更新为正确的四象限布局
根据示意图：
- 左上: 单倍型网络 (180px宽)
- 右上: P值散点图 + 基因结构图 (跨2列)
- 左下: 效应森林图 (180px宽)
- 左中: 箱线图 (120px宽)
- 右下: 单倍型序列 (剩余宽度)
"""
import os
import re
import glob
import json

def patch_html_v3(filepath):
    """更新HTML文件为正确的四象限布局"""
    print(f"处理: {filepath}")
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    changes = []
    
    # 1. 更新CSS为3列布局
    old_grid = r'grid-template-columns:\s*380px\s+1fr'
    new_grid = 'grid-template-columns: 180px 120px 1fr'
    if re.search(old_grid, content):
        content = re.sub(old_grid, new_grid, content)
        changes.append("[1] Grid布局更新为3列 (180px, 120px, 1fr)")
    
    # 2. 更新quadrant-top-right跨2列
    content = re.sub(r'\.quadrant-top-right\s*\{\s*grid-column:\s*2;',
                     '.quadrant-top-right { grid-column: 2 / 4;', content)
    changes.append("[2] quadrant-top-right跨2列")
    
    # 3. 更新网络图高度
    content = re.sub(r'\.network-panel\s*\{[^}]*height:\s*300px',
                     '.network-panel { width: 100%; height: 250px', content)
    changes.append("[3] 网络图高度调整为250px")
    
    # 4. 更新P值散点图高度
    content = re.sub(r'\.pvalue-scatter\s*\{[^}]*height:\s*150px',
                     '.pvalue-scatter { width: 100%; height: 120px', content)
    changes.append("[4] P值散点图高度调整为120px")
    
    # 5. 添加或更新左中和右下象限CSS
    if '.quadrant-bottom-center' not in content:
        css_add = '''
        /* 左中：箱线图 */
        .quadrant-bottom-center { grid-column: 2; grid-row: 2; }
        .boxplot-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; height: 100%; }
        
        /* 右下：单倍型序列 */
        .quadrant-bottom-right { grid-column: 3; grid-row: 2; }
        .sequence-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; overflow-x: auto; height: 100%; }'''
        content = re.sub(r'(\.quadrant-bottom-left\s*\{[^}]*\})', r'\1' + css_add, content)
        changes.append("[5] 添加quadrant-bottom-center和更新quadrant-bottom-right CSS")
    
    # 6. 更新左下象限为效应森林图
    content = re.sub(r'class="effect-box-panel"', 'class="effect-forest-panel"', content)
    content = re.sub(r'id="effect-summary-chart"', 'id="effect-forest-chart"', content)
    changes.append("[6] 左下象限改为effect-forest-panel")
    
    # 7. 更新HTML结构 - 添加左中象限
    if 'quadrant-bottom-center' not in content:
        # 找到左下象限后的位置
        pattern = r'(<!-- 左下：效应[^>]*-->.*?<div id="effect-forest-chart"[^>]*></div>\s*</div>\s*</div>)'
        replacement = r'''\1
            
            <!-- 左中：箱线图 -->
            <div class="quadrant-bottom-center">
                <div class="boxplot-panel">
                    <div id="boxplot-chart" style="width:100%;height:100%;"></div>
                </div>
            </div>'''
        if re.search(pattern, content, re.DOTALL):
            content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            changes.append("[7] 添加左中象限HTML")
    
    # 8. 更新右下象限为单倍型序列
    if 'sequence-table-container' not in content:
        pattern = r'(<!-- 右下：[^>]*-->.*?<div class="sequence-panel"[^>]*>)'
        replacement = r'<div class="quadrant-bottom-right">\n                <div class="sequence-panel" id="sequence-table-container">'
        if re.search(pattern, content, re.DOTALL):
            content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            changes.append("[8] 更新右下象限HTML")
    
    # 9. 更新JavaScript函数
    # 更新drawEffectSummary为drawEffectForest
    content = re.sub(r'function drawEffectSummary\(\)', 'function drawEffectForest()', content)
    content = re.sub(r'getElementById\([\'"]effect-summary-chart[\'"]\)', 'getElementById("effect-forest-chart")', content)
    changes.append("[9] JS函数名更新")
    
    # 10. 添加drawBoxplot函数（如果不存在）
    if 'function drawBoxplot()' not in content:
        boxplot_func = '''
// ==================== 箱线图绘制 ====================
function drawBoxplot() {
    var container = document.getElementById('boxplot-chart');
    if (!container || typeof haplotypeData === 'undefined' || haplotypeData.length === 0) return;
    
    var W = container.clientWidth || 100;
    var H = container.clientHeight || 200;
    var margin = {top: 10, right: 10, bottom: 25, left: 30};
    var width = W - margin.left - margin.right;
    var height = H - margin.top - margin.bottom;
    
    container.innerHTML = '';
    
    var svg = d3.select(container).append('svg').attr('width', W).attr('height', H);
    var g = svg.append('g').attr('transform', 'translate('+margin.left+','+margin.top+')');
    
    var boxData = haplotypeData.map(function(h) {
        return {
            name: h.name,
            min: h.mean ? h.mean - 1 : 0,
            q1: h.mean ? h.mean - 0.5 : 0.25,
            median: h.median || h.mean || 0.5,
            q3: h.mean ? h.mean + 0.5 : 0.75,
            max: h.mean ? h.mean + 1 : 1,
            color: h.color
        };
    }).sort(function(a, b) { return a.median - b.median; });
    
    var yScale = d3.scaleBand().domain(boxData.map(function(d){return d.name;})).range([0, height]).padding(0.3);
    var xScale = d3.scaleLinear().domain([0, 1]).range([0, width]);
    
    g.append('g').attr('transform', 'translate(0,'+height+')').call(d3.axisBottom(xScale).ticks(3));
    g.append('g').call(d3.axisLeft(yScale).tickSize(0));
    
    boxData.forEach(function(d) {
        var y = yScale(d.name);
        var bw = yScale.bandwidth();
        
        g.append('rect').attr('x', xScale(d.q1)).attr('y', y + 2)
            .attr('width', xScale(d.q3) - xScale(d.q1)).attr('height', bw - 4)
            .attr('fill', d.color || '#3498db').attr('opacity', 0.3).attr('stroke', d.color || '#3498db');
        
        g.append('line').attr('x1', xScale(d.median)).attr('x2', xScale(d.median))
            .attr('y1', y).attr('y2', y + bw).attr('stroke', '#333').attr('stroke-width', 2);
    });
    
    g.append('text').attr('x', width/2).attr('y', height+20).attr('text-anchor', 'middle')
        .style('font-size', '9px').text('Value');
}
'''
        # 在drawEffectForest后添加
        content = re.sub(r'(function drawEffectForest\(\)[\s\S]*?\n\})', r'\1' + boxplot_func, content)
        changes.append("[10] 添加drawBoxplot函数")
    
    # 11. 更新DOMContentLoaded
    if 'drawBoxplot();' not in content:
        content = re.sub(r"(drawEffectForest\(\);)", r"\1\n    drawBoxplot();", content)
        changes.append("[11] DOMContentLoaded添加drawBoxplot调用")
    
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
        if patch_html_v3(f):
            updated += 1
        print()
    
    print(f"[完成] 更新 {updated}/{len(html_files)} 个文件")

if __name__ == '__main__':
    main()
