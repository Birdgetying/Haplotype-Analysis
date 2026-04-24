# -*- coding: utf-8 -*-
"""验证数据库完整性和annotation字段"""
import pandas as pd
import json, os

test_gene = 'HORVU.MOREX.r3.2HG0213430'
db_dir = 'database_text'
variant_info_path = os.path.join(db_dir, test_gene, 'variant_info.csv')
gene_info_path = os.path.join(db_dir, test_gene, 'gene_info.json')

# Load variant_info
df = pd.read_csv(variant_info_path)
print('=== variant_info.csv ===')
print(f'Columns: {list(df.columns)}')
has_ann = 'annotation' in df.columns
print(f'Has annotation column: {has_ann}')
if has_ann:
    print(f'Annotation distribution: {dict(df["annotation"].value_counts())}')

# Load gene_info
with open(gene_info_path) as f:
    gi = json.load(f)
print(f'\n=== gene_info.json ===')
print(f'gene_start={gi["gene_start"]}, gene_end={gi["gene_end"]}, strand={gi["strand"]}')
print(f'CDS intervals: {gi["cds"]}')
print(f'Exon intervals: {gi["exons"]}')

cds = gi['cds']
exons = gi['exons']
gene_start = gi['gene_start']
gene_end = gi['gene_end']

# Count by region
in_cds_count = 0
in_exon_not_cds = 0
in_gene_not_exon = 0
upstream = 0
for _, row in df.iterrows():
    pos = row['position']
    in_cds = any(s <= pos <= e for s, e in cds)
    in_exon = any(s <= pos <= e for s, e in exons)
    if in_cds:
        in_cds_count += 1
    elif in_exon:
        in_exon_not_cds += 1
    elif gene_start <= pos <= gene_end:
        in_gene_not_exon += 1
    else:
        upstream += 1

print(f'\nPosition analysis:')
print(f'  In CDS: {in_cds_count} (should get missense/synonymous with FASTA)')
print(f'  In exon not CDS (UTR): {in_exon_not_cds}')
print(f'  In gene not exon (intron): {in_gene_not_exon}')
print(f'  Upstream/promoter: {upstream}')
print(f'  Total: {len(df)}')

snp_count = sum(1 for _, r in df.iterrows() if r['len_diff'] == 0)
indel_count = sum(1 for _, r in df.iterrows() if r['len_diff'] != 0)
print(f'\nSNP: {snp_count}, Indel: {indel_count}')

# CDS SNPs
cds_snps = [(r['position'], r['ref'], r['alt']) for _, r in df.iterrows() 
            if r['len_diff'] == 0 and any(s <= r['position'] <= e for s, e in cds)]
print(f'\nCDS SNPs (candidates for missense/synonymous): {len(cds_snps)}')
for pos, ref, alt in cds_snps:
    print(f'  pos={pos}, ref={ref}, alt={alt}')

# Test the snp_effects rebuild logic
print('\n=== Test snp_effects rebuild (new logic) ===')
snp_effects = {}
for _, row in df.iterrows():
    pos = row['position']
    ann = row.get('annotation', 'other') if has_ann else 'other'
    len_diff = row.get('len_diff', 0)
    is_sv = row.get('is_sv', False)
    
    if has_ann and ann and ann != 'other':
        snp_effects[pos] = ann
    else:
        if is_sv or abs(len_diff) >= 50:
            snp_effects[pos] = 'SV'
        elif len_diff > 0:
            snp_effects[pos] = 'INS'
        elif len_diff < 0:
            snp_effects[pos] = 'DEL'
        else:
            snp_effects[pos] = 'other'

# Annotation distribution
ann_dist = {}
for v in snp_effects.values():
    ann_dist[v] = ann_dist.get(v, 0) + 1
print(f'Rebuilt annotation distribution: {ann_dist}')
print(f'\nExpected after FASTA rebuild on server:')
print(f'  - {in_cds_count} CDS SNPs should become missense_*/synonymous')
print(f'  - {in_exon_not_cds} UTR positions')
print(f'  - {in_gene_not_exon} intron positions')
print(f'  - {upstream} promoter positions')
