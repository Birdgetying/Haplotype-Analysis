#!/usr/bin/env python3
"""
补丁脚本：
1. GWAS图：用竖线长度表示P值（删除圆点）
2. 网络图：扩大显示范围，允许滚轮缩放内部节点
"""
import os
import re

def patch_html(html_path):
    if not os.path.exists(html_path):
        print(f"[ERROR] 文件不存在: {html_path}")
        return False
    
    with open(html_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # ========== 1. 替换整个drawGWASPlot函数 ==========
    gwas_pattern = r"// =+ 基因结构图[^=]+==+\s*function drawGWASPlot\(data\) \{[\s\S]+?(?=// =+ 页面初始化)"
    
    new_gwas_code = """// ==================== 基因结构图（P值用竖线长度表示） ====================
function drawGWASPlot(data) {
    var container = document.getElementById('gwas-gene-viz');
    if (!container) return;
    var W = container.clientWidth || 680;
    var H = 180;
    var ml = 52, mr = 85, mt = 20, mb = 35;
    var iW = W - ml - mr, iH = H - mt - mb;

    d3.select('#gwas-gene-viz').selectAll('*').remove();

    var svg = d3.select('#gwas-gene-viz').append('svg')
        .attr('width', W).attr('height', H).style('display','block');
    var g = svg.append('g').attr('transform','translate('+ml+','+mt+')');

    var xSc = d3.scaleLinear().domain([regionStart, regionEnd]).range([0, iW]);

    var annColor = { missense:'#e74c3c', synonymous:'#f39c12', UTR:'#9b59b6',
                     indel:'#3498db', promoter:'#27ae60', intron:'#aab', other:'#95a5a6' };

    var tip = document.getElementById('d3-tooltip');

    // 基因结构图
    var geneY = iH * 0.65;
    var geneH = iH * 0.25;
    
    var geneGroup = g.append('g').attr('class', 'gene-structure');
    var gX1 = xSc(geneStart), gX2 = xSc(geneEnd);
    var midY = geneY + geneH / 2;

    geneGroup.append('line')
        .attr('x1', gX1).attr('x2', gX2)
        .attr('y1', midY).attr('y2', midY)
        .attr('stroke', '#2c3e50').attr('stroke-width', 2.5);

    geneGroup.append('rect')
        .attr('x', gX1).attr('y', midY-10)
        .attr('width', Math.max(4, gX2-gX1)).attr('height', 20)
        .attr('fill', '#3498db').attr('rx', 4).attr('opacity', 0.85);

    geneGroup.append('text')
        .attr('x', (gX1+gX2)/2).attr('y', midY+1)
        .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
        .attr('font-size', '12px').attr('fill', '#fff').attr('font-weight', 'bold')
        .text('→');

    geneGroup.append('text')
        .attr('x', gX1).attr('y', midY+26)
        .attr('font-size', '9px').attr('fill', '#666').attr('text-anchor', 'middle')
        .text((geneStart/1e6).toFixed(3)+'M');
    geneGroup.append('text')
        .attr('x', gX2).attr('y', midY+26)
        .attr('font-size', '9px').attr('fill', '#666').attr('text-anchor', 'middle')
        .text((geneEnd/1e6).toFixed(3)+'M');

    // P值竖线（长度反映显著性）
    var pvalGroup = g.append('g').attr('class','pval-lines');
    
    var maxLogP = data.length > 0 ? Math.max(d3.max(data, function(d){ return d.logp; }), 1.5) : 2;
    var lineBaseY = midY - 12;
    var maxLineH = iH * 0.5;
    var lineScale = d3.scaleLinear().domain([0, maxLogP * 1.1]).range([3, maxLineH]);

    pvalGroup.selectAll('.pval-line').data(data).join('line')
        .attr('class', 'pval-line')
        .attr('x1', function(d){ return xSc(d.pos); })
        .attr('x2', function(d){ return xSc(d.pos); })
        .attr('y1', lineBaseY)
        .attr('y2', function(d){ return lineBaseY - lineScale(d.logp); })
        .attr('stroke', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })
        .attr('stroke-width', 2)
        .attr('opacity', function(d){ return d.logp > 1.301 ? 0.9 : 0.6; })
        .style('cursor', 'pointer')
        .on('mouseover', function(e,d) {
            d3.select(this).attr('stroke-width', 4).attr('opacity', 1);
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
            d3.select(this).attr('stroke-width', 2).attr('opacity', d.logp > 1.301 ? 0.9 : 0.6);
            tip.style.display = 'none';
        });

    // 竖线顶端小圆点
    pvalGroup.selectAll('.pval-cap').data(data).join('circle')
        .attr('class', 'pval-cap')
        .attr('cx', function(d){ return xSc(d.pos); })
        .attr('cy', function(d){ return lineBaseY - lineScale(d.logp); })
        .attr('r', 3)
        .attr('fill', function(d){ 
            var ann = d.annotation || 'other';
            return annColor[ann] || annColor.other; 
        })
        .attr('stroke', '#fff').attr('stroke-width', 1)
        .style('pointer-events', 'none');

    // X轴
    g.append('g').attr('transform','translate(0,'+iH+')')
        .call(d3.axisBottom(xSc).ticks(6).tickFormat(function(d){ return (d/1e6).toFixed(2)+'M'; }))
        .selectAll('text').attr('font-size', '9px').attr('transform', 'rotate(-20)').attr('text-anchor', 'end');

    // 图例
    var leg = g.append('g').attr('transform','translate('+(iW+8)+',0)');
    
    var annLeg = [ 
        ['Missense','#e74c3c'], ['Synonymous','#f39c12'], ['UTR','#9b59b6'],
        ['Indel','#3498db'], ['Promoter','#27ae60'], ['Intron','#aab'], ['Other','#95a5a6']
    ];
    annLeg.forEach(function(l,i){
        leg.append('line').attr('x1', 0).attr('x2', 12).attr('y1', i*14).attr('y2', i*14)
           .attr('stroke', l[1]).attr('stroke-width', 2);
        leg.append('text').attr('x', 16).attr('y', i*14+3).attr('font-size', '9px').attr('fill', '#555').text(l[0]);
    });
    
    leg.append('text').attr('x', 0).attr('y', 108).attr('font-size', '8px').attr('fill', '#666').text('Line = -logP');
}

"""
    
    content = re.sub(gwas_pattern, new_gwas_code, content)
    
    # ========== 2. 替换整个drawNetworkPlot函数 ==========
    network_pattern = r"// =+ 单倍型网络图[^=]+==+\s*function drawNetworkPlot\(\) \{[\s\S]+?(?=// =+ 基因结构图)"
    
    new_network_code = """// ==================== 单倍型网络图（D3 force simulation） ====================
function drawNetworkPlot() {
    var container = document.getElementById('network-viz');
    if (!container || networkNodes.length === 0) return;
    var W = 350, H = 280;
    d3.select('#network-viz').selectAll('*').remove();

    var svg = d3.select('#network-viz').append('svg')
        .attr('width', W).attr('height', H)
        .style('display','block');

    svg.append('defs').append('clipPath')
        .attr('id', 'network-clip')
        .append('rect')
        .attr('x', 0).attr('y', 0)
        .attr('width', W).attr('height', H);

    var g = svg.append('g').attr('clip-path', 'url(#network-clip)');
    
    var zoomG = g.append('g').attr('class', 'zoom-group');
    
    var currentScale = 1;
    var minScale = 0.5, maxScale = 2.5;
    
    // 使用原生事件确保滚轮缩放正常工作
    container.addEventListener('wheel', function(e) {
        e.preventDefault();
        e.stopPropagation();
        var delta = e.deltaY > 0 ? 0.9 : 1.1;
        var newScale = Math.max(minScale, Math.min(maxScale, currentScale * delta));
        if (newScale !== currentScale) {
            currentScale = newScale;
            zoomG.attr('transform', 'translate('+W/2+','+H/2+') scale('+currentScale+') translate('+-W/2+','+-H/2+')');
        }
    }, { passive: false });

    var nodes = networkNodes.map(function(d) { return Object.assign({}, d); });
    var nodeById = {};
    nodes.forEach(function(d) { nodeById[d.id] = d; });

    var links = networkEdges
        .filter(function(e) { return nodeById[e.source] && nodeById[e.target]; })
        .map(function(e) { return { source: nodeById[e.source], target: nodeById[e.target], distance: e.distance || 1 }; });

    if (links.length === 0 && nodes.length > 1) {
        for (var i = 0; i < nodes.length - 1; i++) {
            links.push({ source: nodes[i], target: nodes[i+1], distance: 1 });
        }
    }

    var sim = d3.forceSimulation(nodes)
        .force('link',      d3.forceLink(links).distance(function(d) { return Math.min(80, 35 + d.distance * 12); }).strength(0.8))
        .force('charge',    d3.forceManyBody().strength(-180))
        .force('center',    d3.forceCenter(W / 2, H / 2))
        .force('collision', d3.forceCollide().radius(function(d) { return d.size * 0.45 + 4; }))
        .force('x', d3.forceX(W / 2).strength(0.05))
        .force('y', d3.forceY(H / 2).strength(0.05));

    var linkSel = zoomG.append('g').selectAll('line').data(links).join('line')
        .attr('stroke', '#99aabb').attr('stroke-opacity', 0.75)
        .attr('stroke-width', function(d) { return Math.max(1, 3.5 - d.distance * 0.4); });

    var linkLabelSel = zoomG.append('g').selectAll('text').data(links).join('text')
        .attr('font-size','8px').attr('fill','#778899').attr('text-anchor','middle')
        .text(function(d) { return d.distance > 1 ? d.distance : ''; });

    var nodeG = zoomG.append('g').selectAll('g').data(nodes).join('g')
        .style('cursor','pointer')
        .call(d3.drag()
            .on('start', function(e,d) { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; })
            .on('drag',  function(e,d) {
                var r = d.size * 0.5 + 5;
                d.fx = Math.max(r, Math.min(W - r, e.x));
                d.fy = Math.max(r, Math.min(H - r, e.y));
            })
            .on('end',   function(e,d) { if (!e.active) sim.alphaTarget(0); d.fx=null; d.fy=null; })
        );

    nodeG.append('circle')
        .attr('r',      function(d) { return d.size * 0.5; })
        .attr('fill',   function(d) { return d.color; })
        .attr('stroke', '#fff').attr('stroke-width', 2.5);

    nodeG.append('text')
        .attr('text-anchor','middle').attr('dy','0.35em')
        .attr('font-size', function(d) { return Math.min(11, Math.max(7, d.size * 0.28)) + 'px'; })
        .attr('fill','#fff').attr('font-weight','bold').attr('pointer-events','none')
        .text(function(d) { return d.id.replace('Hap','H'); });

    nodeG.append('text')
        .attr('text-anchor','middle').attr('dy', function(d) { return d.size * 0.5 + 13; })
        .attr('font-size','9px').attr('fill','#445566').attr('pointer-events','none')
        .text(function(d) { return 'n='+d.count; });

    var tip = document.getElementById('d3-tooltip');
    nodeG.on('mouseover', function(e,d) {
        d3.select(this).select('circle').attr('stroke','#f39c12').attr('stroke-width',3.5);
        tip.innerHTML = '<b>'+d.id+'</b><br>Count: '+d.count+'<br>Mean: '+d.phenoMean;
        tip.style.display = 'block';
    }).on('mousemove', function(e) {
        tip.style.left = (e.clientX+14)+'px'; tip.style.top = (e.clientY-10)+'px';
    }).on('mouseout', function(e) {
        d3.select(this).select('circle').attr('stroke','#fff').attr('stroke-width',2.5);
        tip.style.display = 'none';
    });

    sim.on('tick', function() {
        linkSel
            .attr('x1', function(d) { return d.source.x; }).attr('y1', function(d) { return d.source.y; })
            .attr('x2', function(d) { return d.target.x; }).attr('y2', function(d) { return d.target.y; });
        linkLabelSel
            .attr('x', function(d) { return (d.source.x + d.target.x)/2; })
            .attr('y', function(d) { return (d.source.y + d.target.y)/2 - 3; });
        nodeG.attr('transform', function(d) { return 'translate('+d.x+','+d.y+')'; });
    });
}

"""
    
    content = re.sub(network_pattern, new_network_code, content)
    
    # 备份并保存
    backup_path = html_path + '.bak4'
    with open(backup_path, 'w', encoding='utf-8') as f:
        with open(html_path, 'r', encoding='utf-8') as orig:
            f.write(orig.read())
    
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"[OK] 已更新: {html_path}")
    print(f"    备份: {backup_path}")
    return True

if __name__ == '__main__':
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
    print("  1. GWAS图：用竖线长度表示P值（删除圆点）")
    print("  2. 网络图：扩大显示范围，允许滚轮缩放内部节点")
    print("  3. 容器边界固定，只能缩放内部内容")
