"""FIX6: Replace snp_effects rebuild to use annotation column"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find the exact block: "for _, row in variant_info_df.iterrows():" inside the database rebuild section
target_found = False
for i in range(len(lines)):
    if ('for _, row in variant_info_df.iterrows():' in lines[i] and 
        i > 9700 and i < 9800):
        # Verify this is inside the snp_effects rebuild
        context = ''.join(lines[max(0,i-10):i])
        if 'snp_effects' in context and 'database_dir' in context:
            # Find block end: logger.info line
            block_end = i + 1
            for k in range(i+1, min(i+30, len(lines))):
                if 'logger.info' in lines[k] and '\u91cd\u5efaSNP\u6ce8\u91ca' in lines[k]:
                    block_end = k
                    break
            
            indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
            inner = indent + '    '
            inner2 = inner + '    '
            inner3 = inner2 + '    '
            
            print(f"Found snp_effects rebuild at L{i+1}-L{block_end+1}")
            print("Old code:")
            for j in range(i, block_end+1):
                print(f"  {lines[j].rstrip()}")
            
            # Build new block
            new_block = []
            new_block.append(f"{indent}has_annotation_col = 'annotation' in variant_info_df.columns\n")
            new_block.append(f"{indent}for _, row in variant_info_df.iterrows():\n")
            new_block.append(f"{inner}pos = row['position']\n")
            new_block.append(f"{inner}ann = row.get('annotation', 'other') if has_annotation_col else 'other'\n")
            new_block.append(f"{inner}len_diff = row.get('len_diff', 0)\n")
            new_block.append(f"{inner}is_sv = row.get('is_sv', False)\n")
            new_block.append(f"{inner}\n")
            new_block.append(f"{inner}# \u5173\u952e\u4fee\u590d: \u4f18\u5148\u4f7f\u7528\u6570\u636e\u5e93\u4e2d\u5df2\u8ba1\u7b97\u7684 annotation\n")
            new_block.append(f"{inner}if has_annotation_col and ann and ann != 'other':\n")
            new_block.append(f"{inner2}snp_effects[pos] = ann\n")
            new_block.append(f"{inner}else:\n")
            new_block.append(f"{inner2}# annotation \u4e3a 'other' \u6216\u7f3a\u5931\u65f6\uff0c\u7528 len_diff/is_sv \u56de\u9000\u5224\u65ad\n")
            new_block.append(f"{inner2}if is_sv or abs(len_diff) >= 50:\n")
            new_block.append(f"{inner2}    snp_effects[pos] = 'SV'\n")
            new_block.append(f"{inner2}elif len_diff > 0:\n")
            new_block.append(f"{inner2}    snp_effects[pos] = 'INS'\n")
            new_block.append(f"{inner2}elif len_diff < 0:\n")
            new_block.append(f"{inner2}    snp_effects[pos] = 'DEL'\n")
            new_block.append(f"{inner2}else:\n")
            new_block.append(f"{inner2}    snp_effects[pos] = 'other'\n")
            # Add annotation distribution logging
            new_block.append(f"{indent}_ann_dist = {{}}\n")
            new_block.append(f"{indent}for v in snp_effects.values():\n")
            new_block.append(f"{indent}    _ann_dist[v] = _ann_dist.get(v, 0) + 1\n")
            new_block.append(f'{indent}logger.info(f"  - \u6570\u636e\u5e93\u6ce8\u91ca\u5206\u5e03: {{_ann_dist}}")\n')
            
            lines[i:block_end+1] = new_block
            target_found = True
            print(f"\nFIX6 applied: replaced L{i+1}-L{block_end+1} with annotation-aware rebuild")
            break

if not target_found:
    print("ERROR: Could not find snp_effects rebuild block!")

with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print("Done.")
