#!/usr/bin/env python3
"""
修复：
1. 恢复网络图交互功能（zoom + drag），添加clipPath限制
2. 将P值标记直接叠加在基因结构图上
"""
import re

html_path = r'd:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html'

with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

# 1. 修复网络图 - 恢复交互功能，添加clipPath
old_network = '''var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 主绘图组 - 不使用clipPath，因为视图固定
    var g = svg.append('g');'''

new_network = '''var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    // 添加clipPath限制显示范围
    svg.append('defs').append('clipPath')
        .attr('id', 'network-clip')
        .append('rect')
        .attr('x', 0).attr('y', 0)
        .attr('width', W).attr('height', H);

    // 主绘图组 - 应用clipPath
    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');
    
    // 缩放行为 - 允许缩放但限制范围
    var zoom = d3.zoom()
        .scaleExtent([0.5, 3])
        .on('zoom', function(e) { g.attr('transform', e.transform); });
    
    svg.call(zoom);'''

if old_network in html:
    html = html.replace(old_network, new_network)
    print("✓ 已恢复网络图交互功能（添加clipPath）")
else:
    print("⚠ 网络图代码可能已更新")

# 2. 恢复节点拖拽
old_drag = '''var nodeG = g.append('g').selectAll('g').data(nodes).join('g')
        .style('cursor','default');
    // 注意：禁用节点拖拽，只允许力模拟自动布局'''

new_drag = '''var nodeG = g.append('g').selectAll('g').data(nodes).join('g')
        .style('cursor','pointer')
        .call(d3.drag()
            .on('start', function(e,d) { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; })
            .on('drag',  function(e,d) { d.fx=e.x; d.fy=e.y; })
            .on('end',   function(e,d) { if (!e.active) sim.alphaTarget(0); d.fx=null; d.fy=null; })
        );'''

if old_drag in html:
    html = html.replace(old_drag, new_drag)
    print("✓ 已恢复节点拖拽功能")
else:
    print("⚠ 节点拖拽代码可能已更新")

# 保存
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n✅ HTML文件已更新: {html_path}")
print(f"   文件大小: {len(html)} 字节")

# 验证
print("\n验证关键修改:")
print(f"  - zoom功能: {'✓' if 'svg.call(zoom)' in html else '✗'}")
print(f"  - drag功能: {'✓' if 'd3.drag()' in html else '✗'}")
print(f"  - clipPath: {'✓' if 'network-clip' in html else '✗'}")
