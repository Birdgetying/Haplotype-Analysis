# -*- coding: utf-8 -*-
"""给haplotype_phenotype_analysis.py添加SV VCF识别"""
import re

target_file = 'd:/Desktop/project1/haplotype_phenotype_analysis.py'

with open(target_file, 'r', encoding='utf-8') as f:
    content = f.read()

# 在 VCF 加载块之后（第9380行后）添加 SV VCF 识别
# 找到这段代码的结束位置
old_text = """                # 5. **新增**: VCF文件（优先使用，与原始数据完全一致）
                vcf_path = os.path.join(gene_db_dir, 'variants.vcf.gz')
                if os.path.exists(vcf_path):
                    try:
                        preloaded_data['vcf_file'] = vcf_path
                        logger.info(f"[数据库] 已找到 VCF 文件: {vcf_path}")
                    except Exception as e:
                        logger.warning(f"[数据库] 加载 VCF 失败: {e}")"""

new_text = """                # 5. **新增**: VCF文件（优先使用，与原始数据完全一致）
                vcf_path = os.path.join(gene_db_dir, 'variants.vcf.gz')
                if os.path.exists(vcf_path):
                    try:
                        preloaded_data['vcf_file'] = vcf_path
                        logger.info(f"[数据库] 已找到 VCF 文件: {vcf_path}")
                    except Exception as e:
                        logger.warning(f"[数据库] 加载 VCF 失败: {e}")

                # 6. **新增**: SV VCF文件（结构变异）
                sv_vcf_path = os.path.join(gene_db_dir, 'sv_variants.vcf.gz')
                if os.path.exists(sv_vcf_path):
                    try:
                        preloaded_data['sv_vcf_file'] = sv_vcf_path
                        logger.info(f"[数据库] 已找到 SV VCF 文件: {sv_vcf_path}")
                    except Exception as e:
                        logger.warning(f"[数据库] 加载 SV VCF 失败: {e}")"""

if old_text in content:
    content = content.replace(old_text, new_text)
    with open(target_file, 'w', encoding='utf-8') as f:
        f.write(content)
    print("[OK] 已添加 SV VCF 识别逻辑到 haplotype_phenotype_analysis.py")
else:
    print("[ERROR] 未找到目标文本")
    # 尝试找到附近的代码
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if 'variants.vcf.gz' in line and 'gene_db_dir' in line:
            print(f"行{i+1}: {line}")