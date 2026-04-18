#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检查新数据库是否包含老师所需的30个基因
"""

import pandas as pd
import sqlite3
import os

# Excel 文件路径
excel_file = r'd:\Desktop\project1\database\大麦备份\大麦GWAS及SV筛选基因-2026.03.22(1).xlsx'

print("=" * 80)
print("提取老师所需的基因列表")
print("=" * 80)

# 读取两个 sheet
df1 = pd.read_excel(excel_file, sheet_name='K-SNP_Indel_CandidateGene')
df2 = pd.read_excel(excel_file, sheet_name='K-SV（自己结果）')

# 基因名在第一列（Unnamed: 0）
genes1 = df1.iloc[:, 0].dropna().tolist()
genes2 = df2.iloc[:, 0].dropna().tolist()

print(f"\nK-SNP_Indel_CandidateGene: {len(genes1)} 个基因")
for g in genes1:
    print(f"  - {g}")

print(f"\nK-SV（自己结果）: {len(genes2)} 个基因")
for g in genes2:
    print(f"  - {g}")

# 合并去重
all_genes = list(set(genes1 + genes2))
print(f"\n总计（去重后）: {len(all_genes)} 个基因")
print(f"\n完整列表:")
for i, g in enumerate(sorted(all_genes), 1):
    print(f"  {i}. {g}")

# 检查数据库
print("\n" + "=" * 80)
print("检查数据库")
print("=" * 80)

db_dir = r'd:\Desktop\project1\database'
db_files = [f for f in os.listdir(db_dir) if f.endswith(('.db', '.sqlite', '.sqlite3'))]

if not db_files:
    print("未找到数据库文件！")
    print("数据库目录中的文件:")
    for f in os.listdir(db_dir):
        print(f"  - {f}")
else:
    for db_file in db_files:
        db_path = os.path.join(db_dir, db_file)
        print(f"\n{'=' * 80}")
        print(f"数据库: {db_file}")
        print(f"{'=' * 80}")
        
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # 获取所有表
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [t[0] for t in cursor.fetchall()]
        print(f"\n表列表: {tables}")
        
        # 检查每个表
        for table_name in tables:
            print(f"\n--- 表: {table_name} ---")
            
            # 获取列名
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [col[1] for col in cursor.fetchall()]
            print(f"列: {columns}")
            
            # 查找基因列
            gene_col = None
            for col in columns:
                if 'gene' in col.lower() or 'Gene' in col or '基因' in col or 'gene_id' in col.lower():
                    gene_col = col
                    break
            
            if not gene_col and 'GeneID' in columns:
                gene_col = 'GeneID'
            
            if gene_col:
                # 获取所有基因
                cursor.execute(f"SELECT DISTINCT {gene_col} FROM {table_name}")
                db_genes = [row[0] for row in cursor.fetchall() if row[0]]
                print(f"\n数据库中的基因数: {len(db_genes)}")
                
                # 检查匹配
                found = [g for g in all_genes if g in db_genes]
                missing = [g for g in all_genes if g not in db_genes]
                
                print(f"\n✅ 找到的基因 ({len(found)}/{len(all_genes)}):")
                for g in found:
                    print(f"  ✓ {g}")
                
                if missing:
                    print(f"\n❌ 缺失的基因 ({len(missing)}):")
                    for g in missing:
                        print(f"  ✗ {g}")
                
                if len(found) == len(all_genes):
                    print(f"\n🎉 所有 {len(all_genes)} 个基因都在数据库中！")
                else:
                    print(f"\n⚠️  缺失 {len(missing)}/{len(all_genes)} 个基因")
            else:
                print("未找到基因列")
        
        conn.close()
