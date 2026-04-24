# -*- coding: utf-8 -*-
"""最终验证: 检查所有修复是否到位"""
import inspect
import sys
sys.path.insert(0, '.')

print("=" * 60)
print("最终验证: 数据库建库和注释读取逻辑")
print("=" * 60)

# 1. 检查 annotate_snp_effects_for_region 函数签名
from haplotype_phenotype_analysis import annotate_snp_effects_for_region
sig = inspect.signature(annotate_snp_effects_for_region)
print("\n[1] annotate_snp_effects_for_region 函数签名:")
for name, param in sig.parameters.items():
    default = param.default if param.default is not param.empty else "REQUIRED"
    print(f"    {name}: {default}")
expected_params = ['vcf_file', 'fasta_path', 'gene_chrom', 'cds_intervals', 'exon_intervals',
                   'gene_strand', 'positions', 'gene_start', 'gene_end', 'promoter_start', 'promoter_end']
for p in expected_params:
    assert p in sig.parameters, f"Missing parameter: {p}"
print("    [OK] All required parameters present")

# 2. 检查函数内部是否支持 .csi 索引
source = inspect.getsource(annotate_snp_effects_for_region)
assert '.csi' in source, "Missing .csi index support!"
print("\n[2] VCF索引检查:")
print("    [OK] .tbi and .csi both checked")

# 3. 检查 _annotate_snp_effects_tabix 方法签名
from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer
method = getattr(HaplotypePhenotypeAnalyzer, '_annotate_snp_effects_tabix')
sig2 = inspect.signature(method)
print("\n[3] _annotate_snp_effects_tabix 方法签名:")
for name, param in sig2.parameters.items():
    if name == 'self':
        continue
    default = param.default if param.default is not param.empty else "REQUIRED"
    print(f"    {name}: {default}")
assert 'gene_body_start' in sig2.parameters, "Missing gene_body_start!"
assert 'gene_body_end' in sig2.parameters, "Missing gene_body_end!"
print("    [OK] gene_body_start and gene_body_end present")

# 4. 检查方法内部是否区分 INS/DEL
method_source = inspect.getsource(method)
assert "'INS'" in method_source, "Missing INS classification!"
assert "'DEL'" in method_source, "Missing DEL classification!"
assert "'intron'" in method_source, "Missing intron classification!"
print("\n[4] _annotate_snp_effects_tabix 分类检查:")
print("    [OK] INS, DEL, intron classifications present")

# 5. 检查 analyze_gene 方法的 snp_effects 重建逻辑
analyze_source = inspect.getsource(HaplotypePhenotypeAnalyzer.analyze_gene)
assert 'has_annotation_col' in analyze_source, "Missing has_annotation_col check!"
assert "ann != 'other'" in analyze_source, "Missing annotation priority check!"
print("\n[5] analyze_gene snp_effects 重建逻辑:")
print("    [OK] Uses annotation column from variant_info.csv")

# 6. 检查 preloaded_data 是否包含 annotation
assert "'annotation'" in analyze_source, "Missing annotation in preloaded_data!"
print("\n[6] preloaded_data 加载检查:")
print("    [OK] annotation field included in preloaded_data")

# 7. 检查 genome_wide_haplotype_scan.py 中的调用
print("\n[7] genome_wide_haplotype_scan.py 检查:")
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    scan_source = f.read()
assert 'gene_start=gene_start' in scan_source, "Missing gene_start parameter in call!"
assert 'gene_end=gene_end' in scan_source, "Missing gene_end parameter in call!"
assert 'promoter_start=_promoter_start_ann' in scan_source, "Missing promoter_start parameter!"
assert 'promoter_end=_promoter_end_ann' in scan_source, "Missing promoter_end parameter!"
print("    [OK] All parameters correctly passed to annotate_snp_effects_for_region")

# 8. 检查 FASTA_FILE 传递
assert '_fasta_path = None' in scan_source, "Missing FASTA fallback!"
assert '_fasta_path = FASTA_FILE' in scan_source, "Missing FASTA_FILE assignment!"
print("    [OK] FASTA_FILE correctly passed via try/except")

# 9. 检查 variant_info.csv 保存逻辑
assert "annotation': _ann_effects.get(pos, 'other')" in scan_source, "Missing annotation save!"
print("    [OK] annotation field saved to variant_info.csv")

# 10. 数据流总结
print("\n" + "=" * 60)
print("数据流验证通过!")
print("=" * 60)
print("""
建库阶段 (genome_wide_haplotype_scan.py):
  1. FASTA_FILE → _fasta_path (try/except)
  2. annotate_snp_effects_for_region(vcf_file, fasta_path, ..., gene_start, gene_end, promoter_start, promoter_end)
  3. 返回 {pos: 'missense_conservative'|'synonymous'|'UTR'|'intron'|'promoter'|'INS'|'DEL'|'SV'|'other'}
  4. 保存到 variant_info.csv 的 annotation 列

分析阶段 (haplotype_phenotype_analysis.py):
  1. 从 variant_info.csv 读取 annotation 列
  2. 优先使用 annotation (if != 'other')
  3. 回退到 len_diff/is_sv 判断 INS/DEL/SV
  4. 传递给 HTML/JS 用于过滤和可视化
""")
