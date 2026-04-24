"""提取所有需要修复区域的精确代码"""
with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()
    lines = content.split('\n')

def show(start, end, label):
    print(f"\n{'='*60}\n{label} (L{start+1}-L{end})\n{'='*60}")
    for i in range(start, min(end, len(lines))):
        print(f"L{i+1}|{lines[i]}")

# 1. annotate_snp_effects_for_region .csi support
show(1352, 1358, "FIX1: .csi support")

# 2. _annotate_snp_effects_tabix signature 
show(9013, 9019, "FIX2: _annotate_snp_effects_tabix sig")

# 3. _annotate_snp_effects_tabix non-exon classification
show(9080, 9115, "FIX3: non-exon classification")

# 4. VCF fallback path
show(9046, 9054, "FIX4: VCF fallback")

# 5. preloaded_data 
show(9320, 9334, "FIX5: preloaded_data")

# 6. analyze_gene snp_effects rebuild  
show(9733, 9770, "FIX6: snp_effects rebuild")

# 7. _annotate_snp_effects_tabix call site
show(9774, 9790, "FIX7: tabix call site")
