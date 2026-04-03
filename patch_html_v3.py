#!/usr/bin/env python3
"""
直接修复 integrated_analysis.html 中的GWAS图代码
"""
import re

html_path = r'd:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html'

with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

# 1. 在棒棒糖竖线前添加X轴线
old_before_stem = '''var tip = document.getElementById('d3-tooltip');

    // 棒棒糖竖线 - 按annotation着色
    var stems = g.append('g').attr('class','gwas-stems').selectAll('.stem').data(data).join('line')
        .attr('class','stem')
        .attr('x1',function(d){return xSc(d.pos);}).attr('x2',function(d){return xSc(d.pos);})
        .attr('y1', gwasH).attr('y2', function(d){return ySc(d.logp);})
        .attr('stroke', function(d){ return annColor[d.annotation] || annColor.other; })
        .attr('stroke-width', 1.5).attr('opacity', 0.8);'''

new_before_stem = '''var tip = document.getElementById('d3-tooltip');

    // 绘制X轴线（作为棒棒糖的基础线）
    g.append('line')
        .attr('x1', 0).attr('x2', iW)
        .attr('y1', gwasH).attr('y2', gwasH)
        .attr('stroke', '#ccc').attr('stroke-width', 1);

    // 棒棒糖竖线 - 从X轴向上延伸到数据点
    var stems = g.append('g').attr('class','gwas-stems').selectAll('.stem').data(data).join('line')
        .attr('class','stem')
        .attr('x1',function(d){return xSc(d.pos);}).attr('x2',function(d){return xSc(d.pos);})
        .attr('y1', gwasH).attr('y2', function(d){return ySc(d.logp);})
        .attr('stroke', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })
        .attr('stroke-width', 1.5).attr('opacity', 0.8);'''

if old_before_stem in html:
    html = html.replace(old_before_stem, new_before_stem)
    print("✓ 已添加X轴线并修复竖线颜色")
else:
    print("⚠ 未找到匹配代码（棒棒糖竖线部分）")

# 2. 修复圆点颜色
old_dot = '''// 棒棒糖圆点 - 按annotation着色（与图例一致）
    var dots = g.append('g').attr('class','gwas-dots').selectAll('.dot').data(data).join('circle')
        .attr('class','dot')
        .attr('cx',function(d){return xSc(d.pos);})
        .attr('cy',function(d){return ySc(d.logp);})
        .attr('r', 4.5)
        .attr('fill', function(d){ return annColor[d.annotation] || annColor.other; })'''

new_dot = '''// 棒棒糖圆点 - 按annotation着色
    var dots = g.append('g').attr('class','gwas-dots').selectAll('.dot').data(data).join('circle')
        .attr('class','dot')
        .attr('cx',function(d){return xSc(d.pos);})
        .attr('cy',function(d){return ySc(d.logp);})
        .attr('r', 4.5)
        .attr('fill', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })'''

if old_dot in html:
    html = html.replace(old_dot, new_dot)
    print("✓ 已修复圆点颜色")
else:
    print("⚠ 未找到匹配代码（圆点部分）")

# 3. 修复tooltip
if "+'Ann: '+d.annotation;" in html:
    html = html.replace("+'Ann: '+d.annotation;", "+'Ann: '+(d.annotation || 'other');")
    print("✓ 已修复tooltip")

# 保存
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n✅ HTML文件已更新: {html_path}")
