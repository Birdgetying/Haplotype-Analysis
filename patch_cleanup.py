"""
清理调试代码：移除临时添加的DEBUG print语句
"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 1. 移除 DEBUG-LD-CALL print
old1 = '            print(f"[DEBUG-LD-CALL] merged_df.shape={assoc_module.merged_df.shape}, cols={list(assoc_module.merged_df.columns)[:6]}")\n            self.reporter.generate_integrated_html('
new1 = '            self.reporter.generate_integrated_html('
if old1 in content:
    content = content.replace(old1, new1, 1)
    print("移除 DEBUG-LD-CALL")
else:
    print("未找到 DEBUG-LD-CALL")

# 2. 移除 DEBUG-LD-START print
old2 = '        print(f"[DEBUG-LD-START] n_dp={n_dp}, has_Haplotype_Seq={\'Haplotype_Seq\' in hap_sample_df.columns}, shape={hap_sample_df.shape}")\n        if n_dp >= 2'
new2 = '        if n_dp >= 2'
if old2 in content:
    content = content.replace(old2, new2, 1)
    print("移除 DEBUG-LD-START")
else:
    print("未找到 DEBUG-LD-START")

with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("文件保存完成")

# 验证
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    verify = f.read()
print("DEBUG-LD-CALL 已删除:", 'DEBUG-LD-CALL' not in verify)
print("DEBUG-LD-START 已删除:", 'DEBUG-LD-START' not in verify)
print("all_bases 仍在:", 'all_bases = []' in verify)

