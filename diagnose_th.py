"""诊断文件中序列列th生成代码的实际内容"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

# 找到for pos in display_positions的位置
occurrences = []
idx = 0
while True:
    idx = content.find('for pos in display_positions', idx)
    if idx == -1:
        break
    occurrences.append(idx)
    idx += 1

print(f"'for pos in display_positions' 出现{len(occurrences)}次")
for occ_idx in occurrences:
    snippet = content[occ_idx:occ_idx+500]
    print(f"\n=== at {occ_idx} ===")
    print(repr(snippet))
