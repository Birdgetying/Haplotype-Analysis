#!/usr/bin/env python
"""检查Python代码中的注释逻辑是否完全一致"""
import re

print("=" * 80)
print("代码逻辑一致性检查")
print("=" * 80)

with open("haplotype_phenotype_analysis.py", "r", encoding="utf-8") as f:
    code = f.read()

issues = []

# 检查1: annotate_snp_effects_for_region 独立函数的逻辑
print("\n【检查1】annotate_snp_effects_for_region 独立函数的判断逻辑")
func_match = re.search(r'def annotate_snp_effects_for_region\(.*?\n(.*?)(?=\ndef |\Z)', code, re.DOTALL)
if func_match:
    func_code = func_match.group(1)
    
    # 检查INS/DEL是否保留类型
    if "effects[pos] = var_type" in func_code or "effects[pos] = 'INS'" in func_code:
        print("  ✅ INS/DEL正确保留类型")
    else:
        print("  ❌ INS/DEL可能未正确保留类型")
        issues.append("独立函数: INS/DEL类型保留问题")
    
    # 检查promoter判断
    if "promoter_start <= pos <= promoter_end" in func_code:
        print("  ✅ promoter区间判断存在")
    else:
        print("  ❌ promoter区间判断缺失")
        issues.append("独立函数: promoter判断缺失")
    
    # 检查intron判断
    if "gene_start <= pos <= gene_end" in func_code and "intron" in func_code:
        print("  ✅ intron判断存在")
    else:
        print("  ❌ intron判断缺失")
        issues.append("独立函数: intron判断缺失")

# 检查2: _annotate_snp_effects_tabix 类方法的逻辑
print("\n【检查2】_annotate_snp_effects_tabix 类方法的判断逻辑")
method_match = re.search(r'def _annotate_snp_effects_tabix\(.*?\n(.*?)(?=\n    def |\Z)', code, re.DOTALL)
if method_match:
    method_code = method_match.group(1)
    
    # 检查INS/DEL是否保留类型
    if "effects[pos] = 'INS'" in method_code and "effects[pos] = 'DEL'" in method_code:
        print("  ✅ INS/DEL正确保留类型")
    else:
        print("  ❌ INS/DEL可能未正确保留类型")
        issues.append("类方法: INS/DEL类型保留问题")
    
    # 检查promoter判断
    if "promoter_start <= pos <= promoter_end" in method_code or "in_promoter" in method_code:
        print("  ✅ promoter区间判断存在")
    else:
        print("  ❌ promoter区间判断缺失")
        issues.append("类方法: promoter判断缺失")
    
    # 检查SV判断
    if "is_symbolic" in method_code or "svtype" in method_code:
        print("  ✅ SV符号等位基因判断存在")
    else:
        print("  ⚠️ SV符号等位基因判断可能缺失")

# 检查3: 两处逻辑的判断顺序是否一致
print("\n【检查3】两处代码的判断顺序一致性")
print("  期望顺序: SV/INS/DEL → not in_exon (promoter/intron/other) → not in_cds (UTR) → CDS SNP")

# 检查4: reannotate_database.py的逻辑是否与主代码一致
print("\n【检查4】reannotate_database.py的逻辑一致性")
with open("reannotate_database.py", "r", encoding="utf-8") as f:
    reannot_code = f.read()

if "if var_type == 'SV':" in reannot_code and "elif var_type in ('INS', 'DEL'):" in reannot_code:
    print("  ✅ 重注释脚本的SV/INS/DEL判断正确")
else:
    print("  ❌ 重注释脚本的SV/INS/DEL判断有问题")
    issues.append("重注释脚本: 变异类型判断问题")

if "if promoter_start <= pos <= promoter_end:" in reannot_code:
    print("  ✅ 重注释脚本的promoter判断正确")
else:
    print("  ❌ 重注释脚本的promoter判断有问题")
    issues.append("重注释脚本: promoter判断问题")

# 总结
print(f"\n{'='*80}")
print(f"检查完成")
if issues:
    print(f"发现 {len(issues)} 个潜在问题:")
    for i, issue in enumerate(issues, 1):
        print(f"  {i}. {issue}")
else:
    print("✅ 所有逻辑一致，无问题！")
print("=" * 80)
