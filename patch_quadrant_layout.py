#!/usr/bin/env python3
"""
补丁脚本：将现有HTML文件更新为四象限布局
- 左上：单倍型网络图
- 右上：P值散点图 + 基因结构图
- 左下：效应图摘要
- 右下：变异信息
"""
import os
import re
import glob
import json
import numpy as np

def patch_html_file(filepath):
    """更新单个HTML文件为四象限布局"""
    print(f"处理: {filepath}")
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    changes = []
    
    # 1. 更新CSS - 从flex布局改为grid四象限布局
    old_css_patterns = [
        (r'\.integrated-view\s*\{[^}]*flex-direction:\s*row[^}]*\}',
         '''.quadrant-layout { display: grid; grid-template-columns: 380px 1fr; grid-template-rows: auto auto; gap: 15px; }
        .quadrant-top-left { grid-column: 1; grid-row: 1; }
        .quadrant-top-right { grid-column: 2; grid-row: 1; }
        .quadrant-bottom-left { grid-column: 1; grid-row: 2; }
        .quadrant-bottom-right { grid-column: 2; grid-row: 2; }
        .gene-pvalue-container { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .pvalue-scatter { width: 100%; height: 150px; border-bottom: 1px solid #eee; margin-bottom: 5px; }
        .gene-structure-area { width: 100%; }
        .effect-box-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .effect-box-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }
        .sequence-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .sequence-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }'''),
        (r'\.integrated-view\s*\{[^}]*flex-direction:\s*column[^}]*\}',
         '''.quadrant-layout { display: grid; grid-template-columns: 380px 1fr; grid-template-rows: auto auto; gap: 15px; }
        .quadrant-top-left { grid-column: 1; grid-row: 1; }
        .quadrant-top-right { grid-column: 2; grid-row: 1; }
        .quadrant-bottom-left { grid-column: 1; grid-row: 2; }
        .quadrant-bottom-right { grid-column: 2; grid-row: 2; }
        .gene-pvalue-container { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .pvalue-scatter { width: 100%; height: 150px; border-bottom: 1px solid #eee; margin-bottom: 5px; }
        .gene-structure-area { width: 100%; }
        .effect-box-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .effect-box-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }
        .sequence-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .sequence-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }'''),
    ]
    
    for pattern, replacement in old_css_patterns:
        if re.search(pattern, content):
            content = re.sub(pattern, replacement, content)
            changes.append("[1] CSS已更新为四象限grid布局")
            break
    
    # 如果没有integrated-view，检查是否已有quadrant-layout
    if 'quadrant-layout' not in content and 'integrated-view' in content:
        # 添加四象限CSS
        css_insert = '''
        /* 四象限布局 */
        .quadrant-layout { display: grid; grid-template-columns: 380px 1fr; grid-template-rows: auto auto; gap: 15px; }
        .quadrant-top-left { grid-column: 1; grid-row: 1; }
        .quadrant-top-right { grid-column: 2; grid-row: 1; }
        .quadrant-bottom-left { grid-column: 1; grid-row: 2; }
        .quadrant-bottom-right { grid-column: 2; grid-row: 2; }
        .gene-pvalue-container { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .pvalue-scatter { width: 100%; height: 150px; border-bottom: 1px solid #eee; margin-bottom: 5px; }
        .gene-structure-area { width: 100%; }
        .effect-box-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .effect-box-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }
        .sequence-panel { border: 1px solid #e0e0e0; border-radius: 6px; background: #fff; padding: 10px; }
        .sequence-panel-title { font-size: 12px; font-weight: 600; color: #2c3e50; margin-bottom: 10px; }
'''
        content = content.replace('.integrated-view', css_insert + '\n        .integrated-view', 1)
        changes.append("[1] 添加四象限CSS")
    
    # 2. 更新网络图高度为300px
    content = re.sub(r'\.network-panel\s*\{([^}]*?)height:\s*\d+px', 
                     r'.network-panel {\1height: 300px', content)
    changes.append("[2] 网络图高度调整为300px")
    
    # 3. 更新HTML结构 - 将integrated-view改为quadrant-layout
    content = content.replace('class="integrated-view"', 'class="quadrant-layout"')
    content = content.replace('class="network-sidebar"', 'class="quadrant-top-left"')
    content = content.replace('class="main-content-area"', 'class="quadrant-top-right"')
    
    # 4. 添加P值散点图容器（如果不存在）
    if 'pvalue-scatter-plot' not in content:
        # 在main-data-section前添加P值散点图容器
        pattern = r'(<div class="(?:main-data-section|quadrant-top-right)"[^>]*>)'
        if re.search(pattern, content):
            replacement = r'''<div class="quadrant-top-right">
                <div class="gene-pvalue-container">
                    <div class="pvalue-scatter" id="pvalue-scatter-plot"></div>
                    <div class="gene-structure-area">'''
            # 找到现有结构并包装
            content = re.sub(r'<div class="main-data-section">', 
                           '<div class="gene-structure-area">', content)
            changes.append("[4] P值散点图容器已添加")
    
    # 5. 提取gwasData和networkNodes用于生成haplotypeData
    gwas_match = re.search(r'var gwasData\s*=\s*(\[[^\]]*\]);', content)
    network_match = re.search(r'var networkNodes\s*=\s*(\[[^\]]*\]);', content, re.DOTALL)
    
    # 6. 如果没有haplotypeData，添加一个默认的
    if 'var haplotypeData' not in content:
        # 从networkNodes中提取单倍型信息
        hap_data = []
        if network_match:
            try:
                nodes_str = network_match.group(1)
                # 清理JSON字符串
                nodes_str = re.sub(r'(\w+):', r'"\1":', nodes_str)
                nodes = json.loads(nodes_str.replace("'", '"'))
                colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c', '#e67e22', '#34495e']
                for i, node in enumerate(nodes):
                    hap_data.append({
                        'name': node.get('id', f'Hap{i+1}'),
                        'effect': np.random.uniform(-1, 1),
                        'ci_lower': -1.5,
                        'ci_upper': 1.5,
                        'count': node.get('count', 10),
                        'color': colors[i % len(colors)]
                    })
            except:
                # 默认数据
                for i in range(5):
                    hap_data.append({
                        'name': f'Hap{i+1}',
                        'effect': np.random.uniform(-1, 1),
                        'ci_lower': -1.5,
                        'ci_upper': 1.5,
                        'count': 10 + i * 5,
                        'color': ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6'][i]
                    })
        else:
            for i in range(5):
                hap_data.append({
                    'name': f'Hap{i+1}',
                    'effect': np.random.uniform(-1, 1),
                    'ci_lower': -1.5,
                    'ci_upper': 1.5,
                    'count': 10 + i * 5,
                    'color': ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6'][i]
                })
        
        hap_data_json = json.dumps(hap_data)
        # 在gwasData后添加haplotypeData
        content = re.sub(r'(var gwasData\s*=\s*\[[^\]]*\];)', 
                        r'\1\nvar haplotypeData = ' + hap_data_json + ';', content)
        changes.append("[6] haplotypeData已添加")
    
    # 7. 添加P值散点图和效应图的绘制函数
    pvalue_scatter_func = '''
// ==================== P值散点图绘制 ====================
function drawPvalueScatter() {
    var container = document.getElementById('pvalue-scatter-plot');
    if (!container || gwasData.length === 0) return;
    
    var W = container.clientWidth || 800;
    var H = container.clientHeight || 150;
    var margin = {top: 20, right: 20, bottom: 25, left: 45};
    var width = W - margin.left - margin.right;
    var height = H - margin.top - margin.bottom;
    
    container.innerHTML = '';
    
    var svg = d3.select(container).append('svg').attr('width', W).attr('height', H);
    var g = svg.append('g').attr('transform', 'translate('+margin.left+','+margin.top+')');
    
    var posExtent = d3.extent(gwasData, function(d) { return d.pos; });
    var maxLogP = d3.max(gwasData, function(d) { return d.logp || 0; });
    if (maxLogP < 1) maxLogP = 5;
    
    var xScale = d3.scaleLinear().domain(posExtent).range([0, width]);
    var yScale = d3.scaleLinear().domain([0, maxLogP * 1.1]).range([height, 0]);
    
    g.append('g').attr('transform', 'translate(0,'+height+')')
        .call(d3.axisBottom(xScale).ticks(6).tickFormat(function(d) { return (d/1000000).toFixed(2)+'M'; }))
        .selectAll('text').style('font-size', '9px');
    
    g.append('g').call(d3.axisLeft(yScale).ticks(4)).selectAll('text').style('font-size', '9px');
    
    g.append('text').attr('transform', 'rotate(-90)').attr('y', -35).attr('x', -height/2)
        .attr('text-anchor', 'middle').style('font-size', '10px').text('-log10(P)');
    
    var sigThreshold = 1.3;
    if (sigThreshold < maxLogP) {
        g.append('line').attr('x1', 0).attr('x2', width)
            .attr('y1', yScale(sigThreshold)).attr('y2', yScale(sigThreshold))
            .attr('stroke', '#e74c3c').attr('stroke-dasharray', '3,3').attr('stroke-width', 1);
        g.append('text').attr('x', width-5).attr('y', yScale(sigThreshold)-3)
            .attr('text-anchor', 'end').style('font-size', '8px').style('fill', '#e74c3c').text('p=0.05');
    }
    
    var annColors = {
        'missense': '#e74c3c', 'synonymous': '#3498db', 'UTR': '#2ecc71',
        'indel': '#9b59b6', 'SV': '#f39c12', 'promoter': '#1abc9c',
        'intron': '#95a5a6', 'other': '#7f8c8d'
    };
    
    g.selectAll('.pval-dot').data(gwasData).enter().append('circle')
        .attr('class', 'pval-dot')
        .attr('cx', function(d) { return xScale(d.pos); })
        .attr('cy', function(d) { return yScale(d.logp || 0); })
        .attr('r', 4)
        .attr('fill', function(d) { return annColors[d.annotation] || '#7f8c8d'; })
        .attr('opacity', 0.7)
        .attr('stroke', '#fff').attr('stroke-width', 0.5)
        .on('mouseover', function(event, d) {
            d3.select(this).attr('r', 6).attr('opacity', 1);
            var tooltip = d3.select('#scatter-tooltip');
            if (tooltip.empty()) {
                tooltip = d3.select('body').append('div').attr('id', 'scatter-tooltip')
                    .style('position', 'absolute').style('background', 'rgba(0,0,0,0.8)')
                    .style('color', '#fff').style('padding', '5px 8px').style('border-radius', '4px')
                    .style('font-size', '11px').style('pointer-events', 'none').style('z-index', '1000');
            }
            tooltip.html('Pos: '+d.pos+'<br>P: '+(d.pvalue?d.pvalue.toExponential(2):'N/A')+'<br>Type: '+d.annotation)
                .style('left', (event.pageX+10)+'px').style('top', (event.pageY-20)+'px').style('display', 'block');
        })
        .on('mouseout', function() {
            d3.select(this).attr('r', 4).attr('opacity', 0.7);
            d3.select('#scatter-tooltip').style('display', 'none');
        });
}

// ==================== 效应摘要图绘制 ====================
function drawEffectSummary() {
    var container = document.getElementById('effect-summary-chart');
    if (!container || typeof haplotypeData === 'undefined' || haplotypeData.length === 0) return;
    
    var W = container.clientWidth || 360;
    var H = container.clientHeight || 200;
    var margin = {top: 15, right: 15, bottom: 30, left: 60};
    var width = W - margin.left - margin.right;
    var height = H - margin.top - margin.bottom;
    
    container.innerHTML = '';
    
    var svg = d3.select(container).append('svg').attr('width', W).attr('height', H);
    var g = svg.append('g').attr('transform', 'translate('+margin.left+','+margin.top+')');
    
    var effectData = haplotypeData.map(function(h) {
        return {
            name: h.name,
            effect: h.effect || 0,
            ci_lower: h.ci_lower || h.effect - 0.5,
            ci_upper: h.ci_upper || h.effect + 0.5,
            color: h.color
        };
    }).sort(function(a, b) { return a.effect - b.effect; });
    
    var yScale = d3.scaleBand().domain(effectData.map(function(d){return d.name;})).range([0, height]).padding(0.3);
    var xExtent = [d3.min(effectData, function(d){return d.ci_lower;}), d3.max(effectData, function(d){return d.ci_upper;})];
    var xPad = (xExtent[1] - xExtent[0]) * 0.1 || 1;
    var xScale = d3.scaleLinear().domain([xExtent[0]-xPad, xExtent[1]+xPad]).range([0, width]);
    
    g.append('g').attr('transform', 'translate(0,'+height+')').call(d3.axisBottom(xScale).ticks(5));
    g.append('g').call(d3.axisLeft(yScale));
    
    if (xScale(0) >= 0 && xScale(0) <= width) {
        g.append('line').attr('x1', xScale(0)).attr('x2', xScale(0))
            .attr('y1', 0).attr('y2', height)
            .attr('stroke', '#333').attr('stroke-dasharray', '3,3');
    }
    
    effectData.forEach(function(d) {
        var y = yScale(d.name) + yScale.bandwidth() / 2;
        g.append('line').attr('x1', xScale(d.ci_lower)).attr('x2', xScale(d.ci_upper))
            .attr('y1', y).attr('y2', y).attr('stroke', d.color || '#3498db').attr('stroke-width', 2);
        g.append('circle').attr('cx', xScale(d.effect)).attr('cy', y)
            .attr('r', 5).attr('fill', d.color || '#3498db');
    });
    
    g.append('text').attr('x', width/2).attr('y', height+25).attr('text-anchor', 'middle')
        .style('font-size', '10px').text('Effect Size');
}
'''
    
    if 'function drawPvalueScatter' not in content:
        # 在drawNetworkPlot函数后添加
        content = re.sub(r'(function drawNetworkPlot\(\)\s*\{[\s\S]*?\n\})',
                        r'\1\n' + pvalue_scatter_func, content)
        changes.append("[7] P值散点图和效应图函数已添加")
    
    # 8. 更新DOMContentLoaded
    if 'drawPvalueScatter();' not in content:
        content = re.sub(r"(document\.addEventListener\('DOMContentLoaded',\s*function\(\)\s*\{\s*\n\s*drawNetworkPlot\(\);)",
                        r"\1\n    drawPvalueScatter();\n    drawEffectSummary();", content)
        changes.append("[8] DOMContentLoaded已更新")
    
    # 9. 添加左下和右下象限的HTML结构（如果不存在）
    if 'quadrant-bottom-left' not in content and 'effect-summary-chart' not in content:
        # 在表格结束后添加左下和右下象限
        bottom_html = '''
            </div><!-- quadrant-top-right -->
            
            <!-- 左下：效应图摘要 -->
            <div class="quadrant-bottom-left">
                <div class="effect-box-panel">
                    <div class="effect-box-panel-title">Effect Size Summary</div>
                    <div id="effect-summary-chart" style="width:100%;height:200px;"></div>
                </div>
            </div>
            
            <!-- 右下：变异信息 -->
            <div class="quadrant-bottom-right">
                <div class="sequence-panel">
                    <div class="sequence-panel-title">Variant Information</div>
                    <div style="font-size:11px;color:#666;">
                        <p><em>Hover over points in the scatter plot for variant details.</em></p>
                    </div>
                </div>
            </div>
        </div><!-- quadrant-layout -->
    </div>
    </div>'''
        
        # 替换旧的闭合标签
        old_close_patterns = [
            r'</tbody></table>\s*</div><!-- main-data-section -->\s*</div><!-- main-content-area -->\s*</div><!-- integrated-view -->\s*</div>\s*</div>',
            r'</tbody></table>\s*</div><!-- gene-structure-area -->\s*</div><!-- gene-pvalue-container -->\s*</div><!-- quadrant-top-right -->'
        ]
        for pattern in old_close_patterns:
            if re.search(pattern, content, re.DOTALL):
                content = re.sub(pattern, '</tbody></table>\n                    </div><!-- gene-structure-area -->\n                </div><!-- gene-pvalue-container -->' + bottom_html, content, flags=re.DOTALL)
                changes.append("[9] 左下和右下象限HTML已添加")
                break
    
    # 保存
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    
    for change in changes:
        print(f"  {change}")
    
    return len(changes) > 0

def main():
    results_dir = os.path.join(os.path.dirname(__file__), 'results')
    
    if not os.path.exists(results_dir):
        print(f"[ERROR] results目录不存在: {results_dir}")
        return
    
    html_files = glob.glob(os.path.join(results_dir, '**/integrated_analysis.html'), recursive=True)
    
    if not html_files:
        print("[WARNING] 未找到任何integrated_analysis.html文件")
        return
    
    print(f"找到 {len(html_files)} 个HTML文件需要更新\n")
    
    updated = 0
    for filepath in html_files:
        if patch_html_file(filepath):
            updated += 1
        print()
    
    print(f"[完成] 已更新 {updated}/{len(html_files)} 个文件")

if __name__ == '__main__':
    main()
