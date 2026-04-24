with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

checks = [
    "--sv-vcf",
    "'sv_variants.vcf.gz'",
    "sv_vcf_file: str = None",
    "sv_extractor = HaplotypeExtractor(sv_vcf_file)",
    "'sv_vcf_file': sv_vcf_file if sv_vcf_file else None",
    "sv_vcf_file=args.sv_vcf",
    "sv_subset_path",
    "create_subset_vcf(sv_vcf_file",
]
for kw in checks:
    idx = c.find(kw)
    status = "OK" if idx >= 0 else "MISS"
    print(f"[{status}] {kw}")

# 检查语法
import ast
try:
    compile(c, 'genome_wide_haplotype_scan.py', 'exec')
    print("\n[OK] Python syntax check PASSED")
except SyntaxError as e:
    print(f"\n[FAIL] Syntax error: {e}")
