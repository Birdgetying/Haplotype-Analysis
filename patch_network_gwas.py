#!/usr/bin/env python3
"""
补丁脚本：修复网络图zoom/pan禁用 + GWAS图直接在基因上显示P值
"""
import os
import re

def patch_html(html_path):
    if not os.path.exists(html_path):
        print(f"[ERROR] 文件不存在: {html_path}")
        return False
    
    with open(html_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # ========== 1. 修复网络图：移除zoom/pan功能 ==========
    # 查找并替换网络图中的zoom相关代码
    old_network_zoom = r'''// 主绘图组 - 应用clipPath\s*
\s*var g = svg\.append\('g'\)\.attr\('clip-path', 'url\(#network-clip\)'\);\s*
\s*
\s*// 缩放行为 - 允许缩放但限制范围\s*
\s*var zoom = d3\.zoom\(\)\s*
\s*\.scaleExtent\(\[0\.5, 3\]\)\s*
\s*\.on\('zoom', function\(e\) \{ g\.attr\('transform', e\.transform\); \}\);\s*
\s*
\s*svg\.call\(zoom\);'''
    
    new_network_zoom = '''// 主绘图组 - 应用clipPath（禁用zoom/pan，只使用clipPath限制）
    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');
    
    // 注意：完全禁用D3 zoom/pan功能，不添加任何zoom行为'''
    
    content = re.sub(old_network_zoom, new_network_zoom, content)
    
    # 简化替换：直接替换关键行
    if 'svg.call(zoom);' in content:
        # 移除zoom调用
        content = re.sub(
            r"var zoom = d3\.zoom\(\)[^;]+;\s*svg\.call\(zoom\);",
            "// zoom功能已禁用",
            content,
            flags=re.DOTALL
        )
    
    # ========== 2. 修复节点拖拽：添加边界限制 ==========
    old_drag = r"\.on\('drag',\s*function\(e,d\)\s*\{\s*d\.fx=e\.x;\s*d\.fy=e\.y;\s*\}\)"
    new_drag = """.on('drag',  function(e,d) {
                // 限制拖拽范围在容器内
                var r = d.size * 0.5 + 5;
                d.fx = Math.max(r, Math.min(W - r, e.x));
                d.fy = Math.max(r, Math.min(H - r, e.y));
            })"""
    
    content = re.sub(old_drag, new_drag, content)
    
    # ========== 3. 修复GWAS图：P值直接在基因上显示 ==========
    # 替换整个drawGWASPlot函数
    old_gwas_pattern = r"// =+ 基因结构图（P值叠加显示） =+\s*function drawGWASPlot\(data\) \{[\s\S]+?// =+ 页面初始化 =+"
    
    new_gwas_code = '''// ==================== 基因结构图（P值直接在基因上显示） ====================
function drawGWASPlot(data) {
    var container = document.getElementById('gwas-gene-viz');
    if (!container) return;
    var W = container.clientWidth || 680;
    var H = 220;  // 减小高度，因为P值直接在基因上显示
    var ml = 52, mr = 85, mt = 25, mb = 35;
    var iW = W - ml - mr, iH = H - mt - mb;

    d3.select('#gwas-gene-viz').selectAll('*').remove();

    var svg = d3.select('#gwas-gene-viz').append('svg')
        .attr('width', W).attr('height', H).style('display','block');
    var g = svg.append('g').attr('transform','translate('+ml+','+mt+')');

    // X轴比例尺 - 使用regionStart/regionEnd
    var xSc = d3.scaleLinear().domain([regionStart, regionEnd]).range([0, iW]);

    // 按annotation统一着色
    var annColor = { missense:'#e74c3c', synonymous:'#f39c12', UTR:'#9b59b6',
                     indel:'#3498db', promoter:'#27ae60', intron:'#aab', other:'#95a5a6' };

    var tip = document.getElementById('d3-tooltip');

    // ---------- 基因结构图（主体层） ----------
    var geneY = iH * 0.5;   // 基因结构图Y位置（居中）
    var geneH = iH * 0.35;  // 基因结构图高度
    
    var geneGroup = g.append('g').attr('class', 'gene-structure');
    var gX1 = xSc(geneStart), gX2 = xSc(geneEnd);
    var midY = geneY + geneH / 2;

    // 基因方向线
    geneGroup.append('line')
        .attr('x1', gX1).attr('x2', gX2)
        .attr('y1', midY).attr('y2', midY)
        .attr('stroke', '#2c3e50').attr('stroke-width', 2.5);

    // 基因主体（矩形）
    geneGroup.append('rect')
        .attr('x', gX1).attr('y', midY-12)
        .attr('width', Math.max(4, gX2-gX1)).attr('height', 24)
        .attr('fill', '#3498db').attr('rx', 4).attr('opacity', 0.85);

    // 方向箭头
    geneGroup.append('text')
        .attr('x', (gX1+gX2)/2).attr('y', midY+1)
        .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
        .attr('font-size', '14px').attr('fill', '#fff').attr('font-weight', 'bold')
        .text('→');

    // 基因起始/终止位置标签
    geneGroup.append('text')
        .attr('x', gX1).attr('y', midY+30)
        .attr('font-size', '9px').attr('fill', '#666').attr('text-anchor', 'middle')
        .text((geneStart/1e6).toFixed(3)+'M');
    geneGroup.append('text')
        .attr('x', gX2).attr('y', midY+30)
        .attr('font-size', '9px').attr('fill', '#666').attr('text-anchor', 'middle')
        .text((geneEnd/1e6).toFixed(3)+'M');

    // ---------- P值标记（直接在基因主体上显示） ----------
    // P值圆点直接放在基因主体的上边缘，大小和透明度反映显著性
    var pvalGroup = g.append('g').attr('class','pval-markers');
    
    // 计算最大logP用于调整圆点大小
    var maxLogP = data.length > 0 ? Math.max(d3.max(data, function(d){ return d.logp; }), 1) : 2;
    var sizeScale = d3.scaleLinear().domain([0, maxLogP]).range([3, 10]);  // 圆点大小范围

    // P值圆点（直接在基因主体上方边缘位置）
    pvalGroup.selectAll('.pval-dot').data(data).join('circle')
        .attr('class', 'pval-dot')
        .attr('cx', function(d){ return xSc(d.pos); })
        .attr('cy', midY - 14)  // 固定在基因主体上边缘
        .attr('r', function(d){ return sizeScale(d.logp); })  // 大小反映P值
        .attr('fill', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })
        .attr('stroke', '#fff').attr('stroke-width', 1.5)
        .attr('opacity', function(d){ return d.logp > 1.301 ? 0.95 : 0.7; })  // 显著的点更亮
        .style('cursor', 'pointer')
        .on('mouseover', function(e,d) {
            d3.select(this).attr('r', sizeScale(d.logp) + 3).attr('stroke-width', 2);
            var sig = d.logp > 1.301 ? ' **' : '';
            tip.innerHTML = '<b>Pos: '+d.pos.toLocaleString()+'</b>'+sig+'<br>'
                +'P = '+d.pvalue.toExponential(2)+'<br>'
                +'-log₁₀P = '+d.logp.toFixed(2)+'<br>'
                +'MAF = '+d.maf.toFixed(3)+'<br>'
                +'Ann: '+(d.annotation || 'other');
            tip.style.display = 'block';
        }).on('mousemove', function(e){
            tip.style.left = (e.clientX+14)+'px'; 
            tip.style.top = (e.clientY-10)+'px';
        }).on('mouseout', function(e,d){
            d3.select(this).attr('r', sizeScale(d.logp)).attr('stroke-width', 1.5);
            tip.style.display = 'none';
        });

    // X轴
    g.append('g').attr('transform','translate(0,'+iH+')')
        .call(d3.axisBottom(xSc).ticks(6).tickFormat(function(d){ return (d/1e6).toFixed(2)+'M'; }))
        .selectAll('text').attr('font-size', '9px').attr('transform', 'rotate(-20)').attr('text-anchor', 'end');

    // 图例
    var leg = g.append('g').attr('transform','translate('+(iW+8)+',5)');
    
    // Annotation颜色图例
    var annLeg = [ 
        ['Missense','#e74c3c'], ['Synonymous','#f39c12'], ['UTR','#9b59b6'],
        ['Indel','#3498db'], ['Promoter','#27ae60'], ['Intron','#aab'], ['Other','#95a5a6']
    ];
    annLeg.forEach(function(l,i){
        leg.append('circle').attr('cx', 4).attr('cy', i*15).attr('r', 4).attr('fill', l[1]);
        leg.append('text').attr('x', 12).attr('y', i*15+4).attr('font-size', '9px').attr('fill', '#555').text(l[0]);
    });
    
    // P值大小说明
    leg.append('text').attr('x', 0).attr('y', 115).attr('font-size', '9px').attr('fill', '#666').text('Size = -logP');
    leg.append('circle').attr('cx', 4).attr('cy', 128).attr('r', 6).attr('fill', '#999');
    leg.append('text').attr('x', 14).attr('y', 131).attr('font-size', '8px').attr('fill', '#888').text('sig.');
    leg.append('circle').attr('cx', 4).attr('cy', 143).attr('r', 3).attr('fill', '#ccc');
    leg.append('text').attr('x', 14).attr('y', 146).attr('font-size', '8px').attr('fill', '#888').text('n.s.');
}

// ==================== 页面初始化 =='''
    
    content = re.sub(old_gwas_pattern, new_gwas_code, content)
    
    # 备份并保存
    backup_path = html_path + '.bak3'
    with open(backup_path, 'w', encoding='utf-8') as f:
        with open(html_path, 'r', encoding='utf-8') as orig:
            f.write(orig.read())
    
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"[OK] 已更新: {html_path}")
    print(f"    备份: {backup_path}")
    return True

if __name__ == '__main__':
    # 更新现有的HTML文件
    html_files = [
        r"d:\Desktop\project1\results\CSIAAS1BG1157200HC\integrated_analysis.html",
    ]
    
    for html_file in html_files:
        if os.path.exists(html_file):
            patch_html(html_file)
        else:
            print(f"[SKIP] 文件不存在: {html_file}")
    
    print("\n[DONE] 补丁应用完成")
    print("修复内容:")
    print("  1. 网络图：完全禁用D3 zoom/pan功能")
    print("  2. 网络图：节点拖拽限制在容器内")
    print("  3. GWAS图：P值直接在基因结构图上显示（圆点在基因上边缘）")
