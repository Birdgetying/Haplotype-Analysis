#!/usr/bin/env python3
"""
修复 integrated_analysis.html 中的问题：
1. 修复网络图clipPath
2. 修复GWAS图竖线位置
3. 修复颜色映射
"""
import re
import json

html_path = r'd:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html'

with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

# 1. 修复网络图clipPath
old_network = '''var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 添加clipPath限制显示范围
    svg.append('defs').append('clipPath')
        .attr('id', 'network-clip')
        .append('rect')
        .attr('x', margin.left).attr('y', margin.top)
        .attr('width', innerW).attr('height', innerH);

    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');
    svg.call(d3.zoom().scaleExtent([0.2, 6]).on('zoom', function(e) { g.attr('transform', e.transform); }));'''

new_network = '''var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 添加clipPath限制显示范围 - 使用SVG坐标系
    svg.append('defs').append('clipPath')
        .attr('id', 'network-clip')
        .append('rect')
        .attr('x', 0).attr('y', 0)
        .attr('width', W).attr('height', H);

    // 主绘图组 - 应用clipPath
    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');
    
    // 缩放行为 - 限制平移范围
    var zoom = d3.zoom()
        .scaleExtent([0.3, 5])
        .extent([[0, 0], [W, H]])
        .translateExtent([[-W, -H], [2*W, 2*H]])
        .on('zoom', function(e) { g.attr('transform', e.transform); });
    
    svg.call(zoom);'''

if old_network in html:
    html = html.replace(old_network, new_network)
    print("✓ 已修复网络图clipPath")
else:
    print("⚠ 网络图代码可能已更新或格式不同")

# 2. 修复GWAS图 - 添加X轴线并修复颜色映射
old_gwas_stem = '''// 棒棒糖竖线 - 按annotation着色
    var stems = g.append('g').attr('class','gwas-stems').selectAll('.stem').data(data).join('line')
        .attr('class','stem')
        .attr('x1',function(d){return xSc(d.pos);}).attr('x2',function(d){return xSc(d.pos);})
        .attr('y1', gwasH).attr('y2', function(d){return ySc(d.logp);})
        .attr('stroke', function(d){ return annColor[d.annotation] || annColor.other; })
        .attr('stroke-width', 1.5).attr('opacity', 0.8);

    // 棒棒糖圆点 - 按annotation着色（与图例一致）
    var dots = g.append('g').attr('class','gwas-dots').selectAll('.dot').data(data).join('circle')
        .attr('class','dot')
        .attr('cx',function(d){return xSc(d.pos);})
        .attr('cy',function(d){return ySc(d.logp);})
        .attr('r', 4.5)
        .attr('fill',function(d){ return annColor[d.annotation] || annColor.other; })'''

new_gwas_stem = '''// 绘制X轴线（作为棒棒糖的基础线）
    g.append('line')
        .attr('x1', 0).attr('x2', iW)
        .attr('y1', gwasH).attr('y2', gwasH)
        .attr('stroke', '#ccc').attr('stroke-width', 1);

    // 棒棒糖竖线 - 从X轴(y=gwasH)向上延伸到数据点
    var stems = g.append('g').attr('class','gwas-stems').selectAll('.stem').data(data).join('line')
        .attr('class','stem')
        .attr('x1', function(d){return xSc(d.pos);})
        .attr('x2', function(d){return xSc(d.pos);})
        .attr('y1', gwasH)
        .attr('y2', function(d){return ySc(d.logp);})
        .attr('stroke', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })
        .attr('stroke-width', 1.5).attr('opacity', 0.8);

    // 棒棒糖圆点 - 按annotation着色
    var dots = g.append('g').attr('class','gwas-dots').selectAll('.dot').data(data).join('circle')
        .attr('class','dot')
        .attr('cx', function(d){return xSc(d.pos);})
        .attr('cy', function(d){return ySc(d.logp);})
        .attr('r', 4.5)
        .attr('fill', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })'''

if old_gwas_stem in html:
    html = html.replace(old_gwas_stem, new_gwas_stem)
    print("✓ 已修复GWAS图竖线和颜色映射")
else:
    print("⚠ GWAS图代码可能已更新或格式不同")

# 3. 修复tooltip中的annotation显示
old_tooltip = "+'Ann: '+d.annotation;"
new_tooltip = "+'Ann: '+(d.annotation || 'other');"

if old_tooltip in html:
    html = html.replace(old_tooltip, new_tooltip)
    print("✓ 已修复tooltip显示")

# 4. 修复CSS边框
html = html.replace('border: 1px solid #e0e0e0;', 'border: none;')
print("✓ 已移除边框")

# 保存
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n✅ HTML文件已更新: {html_path}")
print(f"   文件大小: {len(html)} 字节")

# 验证关键函数存在
print("\n验证关键组件:")
print(f"  - drawNetworkPlot: {'✓' if 'drawNetworkPlot' in html else '✗'}")
print(f"  - drawGWASPlot: {'✓' if 'drawGWASPlot' in html else '✗'}")
print(f"  - clipPath: {'✓' if 'network-clip' in html else '✗'}")
print(f"  - translateExtent: {'✓' if 'translateExtent' in html else '✗'}")
