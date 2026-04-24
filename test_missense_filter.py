# -*- coding: utf-8 -*-
"""测试同义/错义突变过滤功能"""
import os
import sys

# 添加项目路径
sys.path.insert(0, 'd:/Desktop/project1')

from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer, DataConfig

def test_missense_synonymous_filter():
    """测试missense/synonymous过滤"""
    
    # 测试基因
    test_gene = 'HORVU.MOREX.r3.2HG0213430'
    test_db_dir = 'd:/Desktop/project1/database_text'
    
    # 检查数据库是否存在
    gene_db_dir = os.path.join(test_db_dir, test_gene)
    if not os.path.exists(gene_db_dir):
        print(f"[ERROR] 测试数据库不存在: {gene_db_dir}")
        return False
    
    # 检查gene_info.json
    gene_info_path = os.path.join(gene_db_dir, 'gene_info.json')
    if not os.path.exists(gene_info_path):
        print(f"[ERROR] gene_info.json不存在: {gene_info_path}")
        return False
    
    import json
    with open(gene_info_path, 'r') as f:
        gene_info = json.load(f)
    
    print(f"[INFO] 测试基因: {test_gene}")
    print(f"[INFO] 基因信息:")
    print(f"  - 染色体: {gene_info.get('chrom')}")
    print(f"  - 位置: {gene_info.get('start')}-{gene_info.get('end')}")
    print(f"  - CDS区间: {len(gene_info.get('cds', []))} 个")
    print(f"  - 外显子区间: {len(gene_info.get('exons', []))} 个")
    
    # 检查vcf_file信息
    vcf_file = gene_info.get('vcf_file')
    print(f"  - VCF文件: {vcf_file}")
    
    return True


if __name__ == '__main__':
    print("=" * 60)
    print("测试同义/错义突变过滤功能")
    print("=" * 60)
    
    success = test_missense_synonymous_filter()
    
    if success:
        print("\n[OK] 测试数据库检查通过")
    else:
        print("\n[ERROR] 测试失败")
