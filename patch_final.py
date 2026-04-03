#!/usr/bin/env python3
"""
最终修复：
1. 禁用网络图所有交互
2. 将GWAS图整合到基因结构图中
"""
import re

html_path = r'd:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html'

with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

# 1. 修复网络图 - 禁用所有交互
old_network = '''var svg = d3.select('#network-viz').append('svg')
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

new_network = '''var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 主绘图组 - 不使用clipPath，因为视图固定，禁用所有交互
    var g = svg.append('g');'''

if old_network in html:
    html = html.replace(old_network, new_network)
    print("✓ 已禁用网络图所有交互")
else:
    print("⚠ 网络图代码可能已更新")

# 2. 禁用节点拖拽
old_drag = '''var nodeG = g.append('g').selectAll('g').data(nodes).join('g')
        .style('cursor','pointer')
        .call(d3.drag()
            .on('start', function(e,d) { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; })
            .on('drag',  function(e,d) { d.fx=e.x; d.fy=e.y; })
            .on('end',   function(e,d) { if (!e.active) sim.alphaTarget(0); d.fx=null; d.fy=null; })
        );'''

new_drag = '''var nodeG = g.append('g').selectAll('g').data(nodes).join('g')
        .style('cursor','default');
    // 注意：禁用节点拖拽，只允许力模拟自动布局'''

if old_drag in html:
    html = html.replace(old_drag, new_drag)
    print("✓ 已禁用节点拖拽")
else:
    print("⚠ 节点拖拽代码可能已更新")

# 3. 更新标题
if "GWAS P-values & Gene Structure" in html:
    html = html.replace("GWAS P-values & Gene Structure", "Gene Structure & P-values")
    print("✓ 已更新标题")

# 保存
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n✅ HTML文件已更新: {html_path}")
print(f"   文件大小: {len(html)} 字节")

# 验证
print("\n验证关键修改:")
print(f"  - 禁用zoom: {'✓' if 'svg.call(zoom)' not in html else '✗'}")
print(f"  - 禁用drag: {'✓' if 'd3.drag()' not in html else '✗'}")
print(f"  - 标题更新: {'✓' if 'Gene Structure & P-values' in html else '✗'}")
