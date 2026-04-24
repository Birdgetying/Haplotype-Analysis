with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    c = f.read()
keys = ['all_bases = []', 'seq-col-th', 'r2ToColor', 'ldR2Matrix', 'drawLDTriangle', 'ld-triangle-wrapper']
for k in keys:
    print('OK' if k in c else 'MISS', repr(k))
