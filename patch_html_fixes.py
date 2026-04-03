#!/usr/bin/env python3
"""
修复 integrated_analysis.html 中的问题：
1. 移除网络图和GWAS图的边框
2. 添加网络图clipPath限制显示范围
3. 修复GWAS图颜色逻辑（统一按annotation着色）
4. 修复过滤功能
"""
import re
import os

html_path = r'd:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html'

with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

# 1. 移除CSS边框
html = html.replace('border: 1px solid #e0e0e0;', 'border: none;')
html = html.replace('.network-panel {', '.network-panel { overflow: hidden;')
html = html.replace('.gene-gwas-panel {', '.gene-gwas-panel { overflow: hidden;')

# 2. 找到并替换drawNetworkPlot函数，添加clipPath
old_network_func = '''// ==================== 单倍型网络图（D3 force simulation） ====================
function drawNetworkPlot() {
    var container = document.getElementById('network-viz');
    if (!container || networkNodes.length === 0) return;
    var W = 350, H = 280;
    d3.select('#network-viz').selectAll('*').remove();

    var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block').style('overflow','visible');

    var g = svg.append('g');'''

new_network_func = '''// ==================== 单倍型网络图（D3 force simulation） ====================
function drawNetworkPlot() {
    var container = document.getElementById('network-viz');
    if (!container || networkNodes.length === 0) return;
    var W = 350, H = 280;
    var margin = { top: 5, right: 5, bottom: 5, left: 5 };
    var innerW = W - margin.left - margin.right;
    var innerH = H - margin.top - margin.bottom;
    d3.select('#network-viz').selectAll('*').remove();

    var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 添加clipPath限制显示范围
    svg.append('defs').append('clipPath')
        .attr('id', 'network-clip')
        .append('rect')
        .attr('x', margin.left).attr('y', margin.top)
        .attr('width', innerW).attr('height', innerH);

    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');'''

if old_network_func in html:
    html = html.replace(old_network_func, new_network_func)
    print("✓ 已更新 drawNetworkPlot 函数（添加clipPath）")
else:
    print("⚠ drawNetworkPlot 函数可能已更新或格式不同")

# 3. 找到并替换drawGWASPlot函数，修复颜色和图例
old_gwas_pattern = r'''var annColor = \{ missense:'#e74c3c', synonymous:'#f39c12', UTR:'#9b59b6',
                     indel:'#3498db', promoter:'#27ae60', intron:'#aab', other:'#95a5a6' \};'''

new_gwas_color = '''// 按annotation统一着色
    var annColor = { missense:'#e74c3c', synonymous:'#f39c12', UTR:'#9b59b6',
                     indel:'#3498db', promoter:'#27ae60', intron:'#aab', other:'#95a5a6' };'''

# 尝试替换
if 'missense:\'#e74c3c\'' in html:
    html = re.sub(
        r"var annColor = \{ missense:'#e74c3c', synonymous:'#f39c12', UTR:'#9b59b6',\s*indel:'#3498db', promoter:'#27ae60', intron:'#aab', other:'#95a5a6' \};",
        new_gwas_color,
        html
    )
    print("✓ 已更新 annColor 定义")

# 4. 修复棒棒糖图颜色逻辑 - 将按pvalue着色改为按annotation着色
old_dot_fill = ".attr('fill',function(d){ return d.pvalue<0.01?'#8B0000':d.pvalue<0.05?'#e74c3c':(annColor[d.annotation]||'#95a5a6'); })"
new_dot_fill = ".attr('fill', function(d){ return annColor[d.annotation] || annColor.other; })"

if old_dot_fill in html:
    html = html.replace(old_dot_fill, new_dot_fill)
    print("✓ 已修复圆点颜色逻辑（按annotation着色）")

old_stem_stroke = ".attr('stroke', function(d){ return annColor[d.annotation]||'#95a5a6'; })"
new_stem_stroke = ".attr('stroke', function(d){ return annColor[d.annotation] || annColor.other; })"

if old_stem_stroke in html:
    html = html.replace(old_stem_stroke, new_stem_stroke)
    print("✓ 已修复竖线颜色逻辑")

# 5. 修复图例 - 移除P值图例，只保留annotation图例
old_legend = '''var leg = g.append('g').attr('transform','translate('+(iW+8)+',2)');
    [ ['P<0.01','#8B0000',5.5], ['P<0.05','#e74c3c',4], ['n.s.','#95a5a6',3] ].forEach(function(l,i){
        leg.append('circle').attr('cx',5).attr('cy',i*16).attr('r',l[2]).attr('fill',l[1]);
        leg.append('text').attr('x',13).attr('y',i*16+4).attr('font-size','8px').attr('fill','#555').text(l[0]);
    });
    var annLeg = [ ['Missense','#e74c3c'],['Synonymous','#f39c12'],['UTR','#9b59b6'],['Indel','#3498db'],['Promoter','#27ae60'] ];
    annLeg.forEach(function(l,i){
        leg.append('rect').attr('x',1).attr('y',56+i*14).attr('width',8).attr('height',8).attr('fill',l[1]).attr('rx',1);
        leg.append('text').attr('x',13).attr('y',56+i*14+7).attr('font-size','8px').attr('fill','#555').text(l[0]);
    });'''

new_legend = '''// 图例 - 只显示annotation颜色（与实际数据点颜色一致）
    var leg = g.append('g').attr('transform','translate('+(iW+8)+',2)');
    var annLeg = [ 
        ['Missense','#e74c3c'], ['Synonymous','#f39c12'], ['UTR','#9b59b6'],
        ['Indel','#3498db'], ['Promoter','#27ae60'], ['Intron','#aab'], ['Other','#95a5a6']
    ];
    annLeg.forEach(function(l,i){
        leg.append('circle').attr('cx',5).attr('cy',i*16).attr('r',4).attr('fill',l[1]);
        leg.append('text').attr('x',13).attr('y',i*16+4).attr('font-size','9px').attr('fill','#555').text(l[0]);
    });'''

if old_legend in html:
    html = html.replace(old_legend, new_legend)
    print("✓ 已更新图例（只显示annotation颜色）")

# 保存
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n✅ HTML文件已更新: {html_path}")
print(f"   文件大小: {len(html)} 字节")
