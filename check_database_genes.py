#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检查新数据库是否包含老师所需的30个基因
"""

import pandas as pd
import os

# Excel 文件路径
excel_file = r'd:\Desktop\project1\database\大麦备份\大麦GWAS及SV筛选基因-2026.03.22(1).xlsx'

# 读取两个 sheet
print("=" * 80)
print("读取 Excel 文件中的基因列表")
print("=" * 80)

df1 = pd.read_excel(excel_file, sheet_name='K-SNP_Indel_CandidateGene')
df2 = pd.read_excel(excel_file, sheet_name='K-SV(自己结果)')

print(f"\nSheet 'K-SNP_Indel_CandidateGene': {len(df1)} 行")
print(f"Sheet 'K-SV(自己结果)': {len(df2)} 行")

# 查看列名
print(f"\nK-SNP_Indel_CandidateGene 列名: {list(df1.columns)}")
print(f"K-SV(自己结果) 列名: {list(df2.columns)}")

# 提取基因名（假设基因名在某一列，需要查看实际数据）
print("\n" + "=" * 80)
print("查看前几行数据")
print("=" * 80)

print("\nK-SNP_Indel_CandidateGene 前5行:")
print(df1.head())

print("\nK-SV(自己结果) 前5行:")
print(df2.head())

# 尝试找到基因名列
gene_col1 = None
gene_col2 = None

for col in df1.columns:
    if 'gene' in col.lower() or 'Gene' in col or '基因' in col:
        gene_col1 = col
        break

for col in df2.columns:
    if 'gene' in col.lower() or 'Gene' in col or '基因' in col:
        gene_col2 = col
        break

if gene_col1:
    genes1 = df1[gene_col1].dropna().tolist()
    print(f"\nK-SNP_Indel_CandidateGene 基因列: {gene_col1}")
    print(f"基因数量: {len(genes1)}")
    print(f"基因列表: {genes1}")
else:
    print("\n未在 K-SNP_Indel_CandidateGene 中找到基因列")

if gene_col2:
    genes2 = df2[gene_col2].dropna().tolist()
    print(f"\nK-SV(自己结果) 基因列: {gene_col2}")
    print(f"基因数量: {len(genes2)}")
    print(f"基因列表: {genes2}")
else:
    print("\n未在 K-SV(自己结果) 中找到基因列")

# 合并两个 sheet 的基因
all_genes = []
if gene_col1:
    all_genes.extend(genes1)
if gene_col2:
    all_genes.extend(genes2)

all_genes = list(set(all_genes))  # 去重
print(f"\n" + "=" * 80)
print(f"老师所需的基因总数: {len(all_genes)}")
print(f"基因列表: {sorted(all_genes)}")
print("=" * 80)

# 检查数据库目录
db_dir = r'd:\Desktop\project1\database'
print(f"\n检查数据库目录: {db_dir}")

if os.path.exists(db_dir):
    files = os.listdir(db_dir)
    print(f"数据库目录中的文件: {files}")
    
    # 查找数据库文件（可能是 .db, .sqlite, .sqlite3 等）
    db_files = [f for f in files if f.endswith(('.db', '.sqlite', '.sqlite3'))]
    
    if db_files:
        print(f"\n找到数据库文件: {db_files}")
        
        import sqlite3
        
        for db_file in db_files:
            db_path = os.path.join(db_dir, db_file)
            print(f"\n{'=' * 80}")
            print(f"检查数据库: {db_file}")
            print(f"{'=' * 80}")
            
            try:
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                
                # 获取所有表名
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = cursor.fetchall()
                print(f"\n数据库中的表: {[t[0] for t in tables]}")
                
                # 检查每个表
                for table in tables:
                    table_name = table[0]
                    print(f"\n--- 表: {table_name} ---")
                    
                    # 获取列名
                    cursor.execute(f"PRAGMA table_info({table_name})")
                    columns = cursor.fetchall()
                    print(f"列名: {[col[1] for col in columns]}")
                    
                    # 查找可能的基因列
                    gene_column = None
                    for col in columns:
                        col_name = col[1]
                        if 'gene' in col_name.lower() or 'Gene' in col_name or '基因' in col_name:
                            gene_column = col_name
                            break
                    
                    if gene_column:
                        # 获取该列的所有值
                        cursor.execute(f"SELECT DISTINCT {gene_column} FROM {table_name}")
                        db_genes = [row[0] for row in cursor.fetchall() if row[0]]
                        print(f"\n数据库中的基因数量: {len(db_genes)}")
                        
                        # 检查老师所需的基因是否在数据库中
                        found_genes = []
                        missing_genes = []
                        
                        for gene in all_genes:
                            if gene in db_genes:
                                found_genes.append(gene)
                            else:
                                missing_genes.append(gene)
                        
                        print(f"\n找到的基因 ({len(found_genes)}): {found_genes}")
                        print(f"\n缺失的基因 ({len(missing_genes)}): {missing_genes}")
                        
                        if len(found_genes) == len(all_genes):
                            print("\n✅ 所有基因都在数据库中！")
                        else:
                            print(f"\n⚠️  缺失 {len(missing_genes)} 个基因")
                    else:
                        print("未找到基因列")
                
                conn.close()
                
            except Exception as e:
                print(f"错误: {e}")
    else:
        print("\n未找到数据库文件 (.db, .sqlite, .sqlite3)")
else:
    print(f"目录不存在: {db_dir}")
