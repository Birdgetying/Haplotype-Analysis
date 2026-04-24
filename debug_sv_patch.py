with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    c = f.read()

tests = [
    ("vcf_mtime", "'vcf_file': vcf_file,\n            'vcf_mtime':"),
    ("process_single_gene sig", "def process_single_gene(gene_info: dict"),
    ("process_single_gene call", "result = process_single_gene(gene_info, vcf_file"),
    ("extractor extract_region", "positions, hap_df, hap_sample_df = extractor.extract_region"),
    ("variants.vcf.gz", "'variants.vcf.gz',  # 关键"),
]

for name, test in tests:
    idx = c.find(test)
    print(f'{name}: {"FOUND at " + str(idx) if idx >= 0 else "NOT FOUND"}')
    if idx < 0:
        # Find nearby text
        short_test = test[:30].replace('\n', '\\n')
        idx2 = c.find(short_test.replace('\\n', ''))
        if idx2 >= 0:
            print(f'  Hint: found "{short_test}" is different from actual at {idx2}')
            print(f'  Actual: {repr(c[idx2:idx2+100])}')
