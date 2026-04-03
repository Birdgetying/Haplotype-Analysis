#!/usr/bin/env python3
"""检查单倍型网络边连接情况"""
import re
import json
import glob
import os

def analyze_network(html_path):
    print(f"\n{'='*60}")
    print(f"分析: {html_path}")
    
    with open(html_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 提取 networkEdges 数据
    edges_match = re.search(r'var networkEdges = (\[.*?\]);', content, re.DOTALL)
    nodes_match = re.search(r'var networkNodes = (\[.*?\]);', content, re.DOTALL)
    
    if not edges_match or not nodes_match:
        print("  ✗ 无法提取网络数据")
        return
    
    try:
        edges = json.loads(edges_match.group(1))
        nodes = json.loads(nodes_match.group(1))
    except json.JSONDecodeError as e:
        print(f"  ✗ JSON解析错误: {e}")
        return
    
    node_ids = {n['id'] for n in nodes}
    print(f"\n  节点数量: {len(nodes)}")
    print(f"  边数量: {len(edges)}")
    
    # 构建连接图
    connections = {nid: set() for nid in node_ids}
    for e in edges:
        s, t = e['source'], e['target']
        if s in connections and t in connections:
            connections[s].add(t)
            connections[t].add(s)
    
    # 查找连通分量
    visited = set()
    components = []
    
    def dfs(node, component):
        if node in visited:
            return
        visited.add(node)
        component.append(node)
        for neighbor in connections[node]:
            dfs(neighbor, component)
    
    for nid in node_ids:
        if nid not in visited:
            component = []
            dfs(nid, component)
            components.append(component)
    
    print(f"\n  连通分量: {len(components)} 个")
    for i, comp in enumerate(components, 1):
        print(f"    组{i}: {', '.join(sorted(comp))}")
    
    # 检查缺失的边
    if len(components) > 1:
        print(f"\n  ⚠ 警告: 网络被分割成 {len(components)} 个不连通的组!")
        print("  建议: 放宽边生成阈值或添加最小生成树连接")
    else:
        print(f"\n  ✓ 网络完全连通")
    
    # 显示边详情
    print(f"\n  边详情 (source → target, distance):")
    for e in sorted(edges, key=lambda x: (x['source'], x['target'])):
        print(f"    {e['source']} → {e['target']}, dist={e.get('distance', 'N/A')}")

# 查找所有 integrated_analysis.html
html_files = glob.glob(os.path.join('results', '**', 'integrated_analysis.html'), recursive=True)

if not html_files:
    print("未找到 HTML 文件")
else:
    for f in html_files:
        analyze_network(f)
