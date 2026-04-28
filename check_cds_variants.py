#!/usr/bin/env python
"""检查特定变异的ref/alt"""
import pandas as pd

gene_id = "HORVU.MOREX.r3.1HG0080220"
df = pd.read_csv(f"database_text/{gene_id}/variant_info.csv")

# 检查在CDS内的变异
cds_positions = [485703621, 485704064, 485704073, 485704076, 485706228, 485707213, 485707234, 485707318, 485707334]

print("CDS内变异的ref/alt:")
for pos in cds_positions:
    row = df[df['position'] == pos]
    if len(row) > 0:
        r = row.iloc[0]
        print(f"  Pos {pos}: ref={r['ref']}, alt={r['alt']}, len_diff={r['len_diff']}, is_sv={r['is_sv']}, annotation={r['annotation']}")
