"""扫描所有需要修复的位置，检查磁盘上的真实代码状态"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

print("=== 1. annotate_snp_effects_for_region: .csi check (around L1355) ===")
for i in range(1352, 1368):
    print(f"  L{i+1}: {lines[i].rstrip()}")

print("\n=== 2. _annotate_snp_effects_tabix signature (around L9041) ===")
for i in range(9040, 9055):
    print(f"  L{i+1}: {lines[i].rstrip()}")

print("\n=== 3. _annotate_snp_effects_tabix INS/DEL/intron (search) ===")
in_tabix_method = False
for i, line in enumerate(lines):
    if '_annotate_snp_effects_tabix' in line and 'def ' in line:
        in_tabix_method = True
        print(f"  Method starts at L{i+1}")
    if in_tabix_method:
        if "'INS'" in line or "'DEL'" in line or "'intron'" in line:
            print(f"  L{i+1}: {line.rstrip()}")
        if i > 9040 and 'def ' in line and '_annotate_snp_effects_tabix' not in line:
            in_tabix_method = False
            break

print("\n=== 4. analyze_gene: has_annotation_col check ===")
found = False
for i, line in enumerate(lines):
    if 'has_annotation_col' in line:
        found = True
        print(f"  L{i+1}: {line.rstrip()}")
if not found:
    print("  NOT FOUND!")

print("\n=== 5. preloaded_data annotation field ===")
found = False
for i, line in enumerate(lines):
    if 'preloaded_data' in line and 'annotation' in line:
        found = True
        print(f"  L{i+1}: {line.rstrip()}")
if not found:
    print("  NOT FOUND - checking variant_info loading in preloaded_data...")
    for i, line in enumerate(lines):
        if "preloaded_data['variant_info']" in line:
            print(f"  L{i+1}: {line.rstrip()}")
            for j in range(i+1, min(i+15, len(lines))):
                print(f"  L{j+1}: {lines[j].rstrip()}")
            break

print("\n=== 6. analyze_gene snp_effects rebuild (search 'ann_dist') ===")
found = False
for i, line in enumerate(lines):
    if 'ann_dist' in line or ("ann" in line and "!= 'other'" in line):
        found = True
        print(f"  L{i+1}: {line.rstrip()}")
if not found:
    print("  NOT FOUND - searching snp_effects rebuild in analyze_gene...")
    for i, line in enumerate(lines):
        if 'snp_effects' in line and 'variant_info' in line and 'annotation' in line:
            print(f"  L{i+1}: {line.rstrip()}")

print("\n=== 7. VCF fallback intron detection ===")
for i, line in enumerate(lines):
    if 'VCF打开失败' in line or 'VCF open fail' in line.lower():
        print(f"  L{i+1}: {line.rstrip()}")
        for j in range(i+1, min(i+15, len(lines))):
            print(f"  L{j+1}: {lines[j].rstrip()}")
        break
