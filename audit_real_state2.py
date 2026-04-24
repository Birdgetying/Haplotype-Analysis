"""详细检查_annotate_snp_effects_tabix和analyze_gene"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# 检查 _annotate_snp_effects_tabix 方法签名
print("=== _annotate_snp_effects_tabix signature ===")
for i, line in enumerate(lines):
    if 'def _annotate_snp_effects_tabix' in line:
        for j in range(i, min(i+10, len(lines))):
            print(f"  L{j+1}: {lines[j].rstrip()}")
        break

# 检查 _annotate_snp_effects_tabix 中的变异分类逻辑
print("\n=== _annotate_snp_effects_tabix: non-exon classification ===")
for i, line in enumerate(lines):
    if 'def _annotate_snp_effects_tabix' in line:
        start = i
        for j in range(i+1, min(i+300, len(lines))):
            if lines[j].strip().startswith('def ') and j > i+5:
                break
            if 'not in_exon' in lines[j] or 'intron' in lines[j] or "'INS'" in lines[j] or "'DEL'" in lines[j] or 'indel' in lines[j]:
                print(f"  L{j+1}: {lines[j].rstrip()}")
        break

# 检查 analyze_gene 中 snp_effects 重建逻辑
print("\n=== analyze_gene: snp_effects rebuild ===")
in_analyze = False
for i, line in enumerate(lines):
    if 'def analyze_gene' in line:
        in_analyze = True
        continue
    if in_analyze and line.strip().startswith('def ') and 'analyze_gene' not in line:
        break
    if in_analyze and 'snp_effects' in line:
        print(f"  L{i+1}: {line.rstrip()}")

# 检查 _annotate_snp_effects_tabix 调用处
print("\n=== _annotate_snp_effects_tabix call sites ===")
for i, line in enumerate(lines):
    if '_annotate_snp_effects_tabix(' in line and 'def ' not in line:
        for j in range(max(0, i-2), min(i+15, len(lines))):
            print(f"  L{j+1}: {lines[j].rstrip()}")
        print("  ---")
