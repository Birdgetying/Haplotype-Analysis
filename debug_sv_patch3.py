with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

# 检查已知的修改
checks = [
    ("sv_vcf_file in gene_info", "'sv_vcf_file': sv_vcf_file if sv_vcf_file else None"),
    ("sv_vcf_file param in sig", "sv_vcf_file: str = None) -> dict:"),
    ("sv_vcf_file in call", "sv_vcf_file=args.sv_vcf"),
    ("--sv-vcf arg", "--sv-vcf"),
    ("sv_variants.vcf.gz in req", "'sv_variants.vcf.gz'"),
    ("SV VCF extractor", "sv_extractor = HaplotypeExtractor(sv_vcf_file)"),
]

for name, test in checks:
    idx = c.find(test)
    print(f'{name}: {"FOUND at " + str(idx) if idx >= 0 else "NOT FOUND"}')

print()

# 检查 required_files 区域
idx = c.find("required_files = [")
print("=== required_files at", idx, "===")
print(repr(c[idx:idx+400]))
