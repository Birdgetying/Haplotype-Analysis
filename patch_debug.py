"""
直接用Python脚本修改haplotype_phenotype_analysis.py
在关键位置添加debug print
"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 1. 在generate_integrated_html调用前添加debug print
old = '            self.reporter.generate_integrated_html(\n                hap_sample_df=assoc_module.merged_df,  # 使用 merged_df，包含表型数据'
new = '            print(f"[DEBUG-LD-CALL] merged_df.shape={assoc_module.merged_df.shape}, cols={list(assoc_module.merged_df.columns)[:6]}")\n            self.reporter.generate_integrated_html(\n                hap_sample_df=assoc_module.merged_df,  # 使用 merged_df，包含表型数据'
if old in content:
    content = content.replace(old, new, 1)
    print("成功添加DEBUG-LD-CALL")
else:
    print("未找到目标字符串1，尝试搜索...")
    idx = content.find('self.reporter.generate_integrated_html(')
    print(f"找到 generate_integrated_html 在字符 {idx}")
    print(repr(content[idx-100:idx+100]))

# 2. 在LD计算前添加debug print  
old2 = '        ld_r2_matrix = []\n        ld_positions_list = list(display_positions)  # 与display_positions完全一致\n        n_dp = len(ld_positions_list)\n        if n_dp >= 2'
new2 = '        ld_r2_matrix = []\n        ld_positions_list = list(display_positions)  # 与display_positions完全一致\n        n_dp = len(ld_positions_list)\n        print(f"[DEBUG-LD-START] n_dp={n_dp}, has_Haplotype_Seq={\'Haplotype_Seq\' in hap_sample_df.columns}, shape={hap_sample_df.shape}")\n        if n_dp >= 2'
if old2 in content:
    content = content.replace(old2, new2, 1)
    print("成功添加DEBUG-LD-START")
else:
    print("未找到目标字符串2")
    # 尝试找到这段
    idx2 = content.find('ld_r2_matrix = []\n        ld_positions_list')
    print(f"找到 ld_r2_matrix 在字符 {idx2}")
    if idx2 >= 0:
        print(repr(content[idx2:idx2+200]))

with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("文件已保存")
