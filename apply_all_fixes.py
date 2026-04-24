"""Apply all 7 fixes to haplotype_phenotype_analysis.py"""
import re

with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

total_fixes = 0

# ============================================================
# FIX 1: annotate_snp_effects_for_region - add .csi support
# Line 1354-1355: replace single tbi check with tbi+csi check
# ============================================================
for i in range(len(lines)):
    if (lines[i].strip() == 'tabix_ok = False' and 
        i+1 < len(lines) and 
        "vcf_file.endswith('.gz') and os.path.exists(vcf_file + '.tbi') and PYSAM_AVAILABLE:" in lines[i+1] and
        'def ' not in lines[i].strip()):
        indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
        lines[i] = f"{indent}tabix_ok = False\n"
        lines.insert(i+1, f"{indent}tbi_exists = vcf_file and os.path.exists(vcf_file + '.tbi')\n")
        lines.insert(i+2, f"{indent}csi_exists = vcf_file and os.path.exists(vcf_file + '.csi')\n")
        lines[i+3] = f"{indent}if vcf_file and vcf_file.endswith('.gz') and (tbi_exists or csi_exists) and PYSAM_AVAILABLE:\n"
        total_fixes += 1
        print(f"FIX1 applied at line {i+1}: .csi index support added")
        break

# ============================================================
# FIX 2: _annotate_snp_effects_tabix - add gene_body_start/end params
# ============================================================
for i in range(len(lines)):
    if ('def _annotate_snp_effects_tabix' in lines[i] and 
        'promoter_end: int = None) -> dict:' in lines[i+4 if i+4 < len(lines) else i]):
        # Check if gene_body_start already exists
        sig_text = ''.join(lines[i:i+6])
        if 'gene_body_start' not in sig_text:
            # Replace the closing of the signature
            for j in range(i, min(i+8, len(lines))):
                if 'promoter_end: int = None) -> dict:' in lines[j]:
                    indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                    lines[j] = f"{indent}promoter_end: int = None,\n"
                    lines.insert(j+1, f"{indent}gene_body_start: int = None,\n")
                    lines.insert(j+2, f"{indent}gene_body_end: int = None) -> dict:\n")
                    total_fixes += 1
                    print(f"FIX2 applied at line {j+1}: gene_body_start/end params added")
                    break
        break

# ============================================================
# FIX 3: _annotate_snp_effects_tabix - fix non-exon classification
# Replace 'indel' with INS/DEL, add intron detection
# ============================================================
for i in range(len(lines)):
    if 'def _annotate_snp_effects_tabix' in lines[i]:
        method_start = i
        # Find the non-exon classification block
        for j in range(i, min(i+300, len(lines))):
            if j+1 < len(lines) and 'not in_exon' in lines[j].strip() and lines[j].strip().startswith('if not in_exon'):
                # Find the block from 'if not in_exon' to the elif/else
                block_start = j
                # Determine indentation
                indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                inner = indent + '    '
                inner2 = inner + '    '
                
                # Find block end - look for the next elif at same or lesser indent
                block_end = j + 1
                for k in range(j+1, min(j+30, len(lines))):
                    stripped = lines[k].strip()
                    line_indent = lines[k][:len(lines[k]) - len(lines[k].lstrip())]
                    if stripped.startswith('elif') and len(line_indent) <= len(indent):
                        block_end = k
                        break
                    elif stripped.startswith('else:') and len(line_indent) <= len(indent):
                        block_end = k
                        break
                
                # Read existing block to understand structure
                old_block = lines[block_start:block_end]
                print(f"  Old non-exon block (L{block_start+1}-L{block_end}):")
                for line in old_block:
                    print(f"    {line.rstrip()}")
                
                # Build new block with INS/DEL distinction and intron detection
                new_block = []
                new_block.append(f"{indent}if not in_exon:\n")
                new_block.append(f"{inner}# \u975e\u5916\u663e\u5b50\u533a\u57df: \u533a\u5206promoter/intron/SV/INS/DEL\n")
                new_block.append(f"{inner}in_promoter = (promoter_start is not None and\n")
                new_block.append(f"{inner}              promoter_end is not None and\n")
                new_block.append(f"{inner}              promoter_start <= pos <= promoter_end)\n")
                new_block.append(f"{inner}in_gene_body = (gene_body_start is not None and\n")
                new_block.append(f"{inner}                gene_body_end is not None and\n")
                new_block.append(f"{inner}                gene_body_start <= pos <= gene_body_end)\n")
                new_block.append(f"{inner}if len_diff >= 50:\n")
                new_block.append(f"{inner2}effects[pos] = 'SV'\n")
                new_block.append(f"{inner}elif in_promoter:\n")
                new_block.append(f"{inner2}effects[pos] = 'promoter'\n")
                new_block.append(f"{inner}elif in_gene_body:\n")
                new_block.append(f"{inner2}effects[pos] = 'intron'\n")
                new_block.append(f"{inner}elif len(alt_allele) > len(ref):\n")
                new_block.append(f"{inner2}effects[pos] = 'INS'\n")
                new_block.append(f"{inner}elif len(ref) > len(alt_allele):\n")
                new_block.append(f"{inner2}effects[pos] = 'DEL'\n")
                new_block.append(f"{inner}else:\n")
                new_block.append(f"{inner2}effects[pos] = 'other'\n")
                
                lines[block_start:block_end] = new_block
                total_fixes += 1
                print(f"FIX3 applied at line {block_start+1}: INS/DEL/intron classification")
                break
        break

# ============================================================
# FIX 4: VCF fallback - add promoter/intron detection
# ============================================================
for i in range(len(lines)):
    if "VCF\u6253\u5f00\u5931\u8d25" in lines[i]:
        # Find the simple fallback block
        for j in range(i+1, min(i+10, len(lines))):
            if "effects[pos] = 'UTR' if (in_exon and not in_cds) else 'other'" in lines[j]:
                indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                inner = indent + '    '
                # Replace with proper classification
                new_lines = []
                new_lines.append(f"{indent}in_exon = _pos_in_any_interval(pos, exon_intervals)\n")
                new_lines.append(f"{indent}in_cds  = _pos_in_any_interval(pos, cds_intervals)\n")
                new_lines.append(f"{indent}if in_exon and not in_cds:\n")
                new_lines.append(f"{inner}effects[pos] = 'UTR'\n")
                new_lines.append(f"{indent}elif promoter_start and promoter_end and promoter_start <= pos <= promoter_end:\n")
                new_lines.append(f"{inner}effects[pos] = 'promoter'\n")
                new_lines.append(f"{indent}elif gene_body_start and gene_body_end and gene_body_start <= pos <= gene_body_end and not in_exon:\n")
                new_lines.append(f"{inner}effects[pos] = 'intron'\n")
                new_lines.append(f"{indent}else:\n")
                new_lines.append(f"{inner}effects[pos] = 'other'\n")
                # Replace the 3 lines (in_exon, in_cds, effects line) with new block
                lines[j-2:j+1] = new_lines
                total_fixes += 1
                print(f"FIX4 applied at line {j+1}: VCF fallback promoter/intron")
                break
        break

# ============================================================
# FIX 5: preloaded_data - add annotation field
# ============================================================
for i in range(len(lines)):
    if "'missing_rate': row.get('missing_rate', 0.0)" in lines[i]:
        # Check it's inside preloaded_data block
        context = ''.join(lines[max(0,i-10):i+1])
        if 'preloaded_data' in context:
            indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
            # Check if annotation already exists
            if 'annotation' not in lines[i+1]:
                old_line = lines[i].rstrip()
                lines[i] = old_line + ',\n'  # Add comma
                lines.insert(i+1, f"{indent}'annotation': row.get('annotation', 'other')\n")
                total_fixes += 1
                print(f"FIX5 applied at line {i+1}: annotation field in preloaded_data")
        break

# ============================================================
# FIX 6: analyze_gene snp_effects rebuild - use annotation column
# ============================================================
for i in range(len(lines)):
    if "# **\u4f18\u5148\u4ece\u6570\u636e\u5e93\u91cd\u5efa**" in lines[i] or ("\u4f18\u5148\u4ece\u6570\u636e\u5e93\u91cd\u5efa" in lines[i] and "snp_effects" not in lines[i]):
        # Found the database rebuild section, now find the for loop
        for j in range(i, min(i+20, len(lines))):
            if 'for _, row in variant_info_df.iterrows():' in lines[j]:
                # Find the block end (logger.info line)
                block_end = j + 1
                for k in range(j+1, min(j+40, len(lines))):
                    if 'logger.info' in lines[k] and '\u4ece\u6570\u636e\u5e93\u91cd\u5efaSNP\u6ce8\u91ca' in lines[k]:
                        block_end = k
                        break
                
                indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                inner = indent + '    '
                inner2 = inner + '    '
                
                # Build new block
                new_block = []
                new_block.append(f"{indent}has_annotation_col = 'annotation' in variant_info_df.columns\n")
                new_block.append(f"{indent}for _, row in variant_info_df.iterrows():\n")
                new_block.append(f"{inner}pos = row['position']\n")
                new_block.append(f"{inner}ann = row.get('annotation', 'other') if has_annotation_col else 'other'\n")
                new_block.append(f"{inner}len_diff = row.get('len_diff', 0)\n")
                new_block.append(f"{inner}is_sv = row.get('is_sv', False)\n")
                new_block.append(f"{inner}\n")
                new_block.append(f"{inner}# **\u5173\u952e\u4fee\u590d**: \u4f18\u5148\u4f7f\u7528\u6570\u636e\u5e93\u4e2d\u5df2\u8ba1\u7b97\u7684 annotation\n")
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
                
                lines[j:block_end] = new_block
                total_fixes += 1
                print(f"FIX6 applied at line {j+1}: snp_effects rebuild uses annotation column")
                
                # Add annotation distribution logging after the rebuild
                for m in range(j + len(new_block), min(j + len(new_block) + 5, len(lines))):
                    if '\u4ece\u6570\u636e\u5e93\u91cd\u5efaSNP\u6ce8\u91ca' in lines[m]:
                        old_log = lines[m]
                        ann_dist_log = f"{indent}_ann_dist = {{}}\n{indent}for v in snp_effects.values():\n{indent}    _ann_dist[v] = _ann_dist.get(v, 0) + 1\n"
                        lines.insert(m, ann_dist_log)
                        # Update the logger line
                        lines[m+1] = lines[m+1].replace(
                            '\u4ece\u6570\u636e\u5e93\u91cd\u5efaSNP\u6ce8\u91ca: {len(snp_effects)} \u4e2a\u4f4d\u70b9',
                            '\u6570\u636e\u5e93\u6ce8\u91ca\u5206\u5e03: {_ann_dist}'
                        )
                        break
                break
        break

# ============================================================
# FIX 7: _annotate_snp_effects_tabix call site - add gene_body params
# ============================================================
for i in range(len(lines)):
    if '_annotate_snp_effects_tabix(' in lines[i] and 'def ' not in lines[i]:
        # Find the closing parenthesis
        for j in range(i, min(i+20, len(lines))):
            if 'promoter_end=' in lines[j] and ')' in lines[j+1 if j+1 < len(lines) else j]:
                indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                # Check if gene_body_start already there
                if 'gene_body_start' not in lines[j+1]:
                    # Remove the closing paren from next line and add params
                    old_promoter_line = lines[j].rstrip()
                    if old_promoter_line.endswith(','):
                        pass  # already has comma
                    else:
                        lines[j] = old_promoter_line + ',\n'
                    # Insert gene_body params before closing paren
                    close_line = lines[j+1]
                    close_indent = close_line[:len(close_line) - len(close_line.lstrip())]
                    lines.insert(j+1, f"{indent}gene_body_start=start,\n")
                    lines.insert(j+2, f"{indent}gene_body_end=end\n")
                    # The closing paren line is now at j+3
                    total_fixes += 1
                    print(f"FIX7 applied at line {j+1}: gene_body params in call site")
                break
        break

# Write fixed file
with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print(f"\n{'='*60}")
print(f"Total fixes applied: {total_fixes}")
print(f"{'='*60}")
