with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()
    for i in range(1352, 1370):
        print(f"L{i+1}: {repr(lines[i])}")
