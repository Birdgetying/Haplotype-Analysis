"""
修复 haplotype_phenotype_analysis.py 中的两个注释Bug:
1. 无FASTA路径: 缺少 promoter/intron 检测
2. 有FASTA路径 (tabix+scan): 缺少 promoter 检测

同时利用现有 database_text 中的 gene_info.json + variant_info.csv
重新注释所有基因（无需重建数据库）
"""
import os, re

# ===================== PART 1: 修复代码 =====================

with open('haplotype_phenotype_analysis.py', 'r', encoding='utf-8') as f:
    content = f.read()

original_len = len(content)

# --- FIX A: 无FASTA路径补充 promoter/intron 检测 ---
old_nofasta = '''    if not PYSAM_AVAILABLE or not fasta_path or not os.path.exists(fasta_path):
        # 没有 FASTA，只能判断 UTR vs other
        for pos in positions:
            in_cds = _pos_in_any_interval(pos, cds_intervals)
            in_exon = _pos_in_any_interval(pos, exon_intervals)
            if in_exon and not in_cds:
                effects[pos] = 'UTR'
            elif in_cds:
                effects[pos] = 'other'  # 无法判断 missense/synonymous
        return effects'''

new_nofasta = '''    if not PYSAM_AVAILABLE or not fasta_path or not os.path.exists(fasta_path):
        # 没有 FASTA，无法判断 missense/synonymous，但仍可判断 promoter/intron/UTR
        for pos in positions:
            in_cds = _pos_in_any_interval(pos, cds_intervals)
            in_exon = _pos_in_any_interval(pos, exon_intervals)
            if in_exon and not in_cds:
                effects[pos] = 'UTR'
            elif in_cds:
                effects[pos] = 'other'  # 无法判断 missense/synonymous
            elif promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:
                effects[pos] = 'promoter'
            elif gene_start is not None and gene_end is not None and gene_start <= pos <= gene_end:
                effects[pos] = 'intron'
            # else: keep 'other' (intergenic)
        return effects'''

if old_nofasta in content:
    content = content.replace(old_nofasta, new_nofasta, 1)
    print('[FIX A] 无FASTA路径 promoter/intron 检测已修复')
else:
    print('[ERROR] 无FASTA路径未找到匹配文本，请手动检查')

# --- FIX B: 有FASTA路径(tabix)补充 promoter 检测 ---
# 当前tabix路径: elif not in_exon: check gene_start/gene_end → 'intron' or 'other'
# 需要在 intron 判断前先检查 promoter

old_tabix_nonexon = '''                        elif not in_exon:
                            # 区分内含子和基因间区
                            if gene_start and gene_end and gene_start <= pos <= gene_end:
                                effects[pos] = 'intron'  # 在基因边界内但不在外显子中 = 内含子
                            else:
                                effects[pos] = 'other'  # 基因间区
                        elif not in_cds:
// This is the omitted part'''

# Let's search more carefully
pattern_tabix = (
    '                        elif not in_exon:\n'
    '                            # 区分内含子和基因间区\n'
    '                            if gene_start and gene_end and gene_start <= pos <= gene_end:\n'
    '                                effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
    '                            else:\n'
    '                                effects[pos] = \'other\'  # 基因间区\n'
    '                        elif not in_cds:'
)

replace_tabix = (
    '                        elif not in_exon:\n'
    '                            # 区分启动子/内含子/基因间区\n'
    '                            if promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:\n'
    '                                effects[pos] = \'promoter\'\n'
    '                            elif gene_start and gene_end and gene_start <= pos <= gene_end:\n'
    '                                effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
    '                            else:\n'
    '                                effects[pos] = \'other\'  # 基因间区\n'
    '                        elif not in_cds:'
)

if pattern_tabix in content:
    content = content.replace(pattern_tabix, replace_tabix, 1)
    print('[FIX B] tabix路径 promoter 检测已修复')
else:
    print('[WARN B] tabix路径未找到匹配，搜索变体...')
    # Try searching without comment
    p2 = (
        '                        elif not in_exon:\n'
        '                            if gene_start and gene_end and gene_start <= pos <= gene_end:\n'
        '                                effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
        '                            else:\n'
        '                                effects[pos] = \'other\'  # 基因间区\n'
        '                        elif not in_cds:'
    )
    if p2 in content:
        r2 = (
            '                        elif not in_exon:\n'
            '                            # 区分启动子/内含子/基因间区\n'
            '                            if promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:\n'
            '                                effects[pos] = \'promoter\'\n'
            '                            elif gene_start and gene_end and gene_start <= pos <= gene_end:\n'
            '                                effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
            '                            else:\n'
            '                                effects[pos] = \'other\'  # 基因间区\n'
            '                        elif not in_cds:'
        )
        content = content.replace(p2, r2, 1)
        print('[FIX B] tabix路径 promoter 检测已修复(变体)')
    else:
        print('[ERROR B] tabix路径所有尝试均失败，手动检查行1391附近')

# --- FIX C: 有FASTA路径(scan)补充 promoter 检测 ---
old_scan_nonexon = (
    '                    elif not in_exon:\n'
    '                        # 区分内含子和基因间区\n'
    '                        if gene_start and gene_end and gene_start <= pos <= gene_end:\n'
    '                            effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
    '                        else:\n'
    '                            effects[pos] = \'other\'  # 基因间区\n'
    '                    elif not in_cds:'
)

new_scan_nonexon = (
    '                    elif not in_exon:\n'
    '                        # 区分启动子/内含子/基因间区\n'
    '                        if promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:\n'
    '                            effects[pos] = \'promoter\'\n'
    '                        elif gene_start and gene_end and gene_start <= pos <= gene_end:\n'
    '                            effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
    '                        else:\n'
    '                            effects[pos] = \'other\'  # 基因间区\n'
    '                    elif not in_cds:'
)

if old_scan_nonexon in content:
    content = content.replace(old_scan_nonexon, new_scan_nonexon, 1)
    print('[FIX C] scan路径 promoter 检测已修复')
else:
    # Try without comment
    p3 = (
        '                    elif not in_exon:\n'
        '                        if gene_start and gene_end and gene_start <= pos <= gene_end:\n'
        '                            effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
        '                        else:\n'
        '                            effects[pos] = \'other\'  # 基因间区\n'
        '                    elif not in_cds:'
    )
    if p3 in content:
        r3 = (
            '                    elif not in_exon:\n'
            '                        # 区分启动子/内含子/基因间区\n'
            '                        if promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:\n'
            '                            effects[pos] = \'promoter\'\n'
            '                        elif gene_start and gene_end and gene_start <= pos <= gene_end:\n'
            '                            effects[pos] = \'intron\'  # 在基因边界内但不在外显子中 = 内含子\n'
            '                        else:\n'
            '                            effects[pos] = \'other\'  # 基因间区\n'
            '                    elif not in_cds:'
        )
        content = content.replace(p3, r3, 1)
        print('[FIX C] scan路径 promoter 检测已修复(变体)')
    else:
        print('[ERROR C] scan路径未找到，手动检查行1458附近')

print(f'\n代码长度: {original_len} -> {len(content)} (差值: {len(content)-original_len})')

# 写回文件
with open('haplotype_phenotype_analysis.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('代码修复完成，已写入文件\n')

# ===================== PART 2: 语法检查 =====================
import py_compile, sys
try:
    py_compile.compile('haplotype_phenotype_analysis.py', doraise=True)
    print('[OK] 语法检查通过')
except py_compile.PyCompileError as e:
    print(f'[ERROR] 语法错误: {e}')
    sys.exit(1)

# ===================== PART 3: 重新注释 database_text =====================
print('\n开始重新注释 database_text...')
import json, csv

def pos_in_intervals(pos, intervals):
    for s, e in intervals:
        if s <= pos <= e:
            return True
    return False

db_dir = 'database_text'
genes = [d for d in os.listdir(db_dir) if os.path.isdir(os.path.join(db_dir, d))]
total_updated = 0

for g in sorted(genes):
    gdir = os.path.join(db_dir, g)
    gene_info_path = os.path.join(gdir, 'gene_info.json')
    variant_info_path = os.path.join(gdir, 'variant_info.csv')
    
    if not os.path.exists(gene_info_path) or not os.path.exists(variant_info_path):
        print(f'  {g}: 跳过（文件缺失）')
        continue
    
    with open(gene_info_path, encoding='utf-8') as f:
        gi = json.load(f)
    
    cds_list = gi.get('cds', [])
    exon_list = gi.get('exons', [])
    gene_start = gi.get('gene_start')
    gene_end = gi.get('gene_end')
    promoter_start = gi.get('promoter_start')
    promoter_end = gi.get('promoter_end')
    
    # 转换格式
    cds_intervals = [(int(c['start']), int(c['end'])) for c in cds_list if 'start' in c and 'end' in c]
    exon_intervals = [(int(e['start']), int(e['end'])) for e in exon_list if 'start' in e and 'end' in e]
    
    # 读取 variant_info.csv
    with open(variant_info_path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    
    if not rows:
        continue
    
    updated = 0
    for row in rows:
        try:
            pos = int(row.get('position', 0))
        except:
            continue
        
        try:
            len_diff = int(row.get('len_diff', 0))
        except:
            len_diff = 0
        
        is_sv = str(row.get('is_sv', '')).lower() in ('true', '1', 'yes')
        
        in_cds = pos_in_intervals(pos, cds_intervals)
        in_exon = pos_in_intervals(pos, exon_intervals)
        
        # 先按位置分类
        if in_exon and not in_cds:
            loc_ann = 'UTR'
        elif in_cds:
            loc_ann = 'cds'  # 内部标记，后续再按变异类型区分
        elif promoter_start is not None and promoter_end is not None and promoter_start <= pos <= promoter_end:
            loc_ann = 'promoter'
        elif gene_start is not None and gene_end is not None and gene_start <= pos <= gene_end:
            loc_ann = 'intron'
        else:
            loc_ann = 'other'
        
        # 按变异类型和位置综合分类
        if is_sv or abs(len_diff) >= 50:
            new_ann = 'SV'
        elif len_diff > 0:
            new_ann = 'INS'
        elif len_diff < 0:
            new_ann = 'DEL'
        else:
            # SNP (len_diff == 0)
            if loc_ann == 'UTR':
                new_ann = 'UTR'
            elif loc_ann == 'cds':
                # 没有FASTA，无法判断 missense/synonymous
                new_ann = 'other'
            elif loc_ann in ('promoter', 'intron'):
                new_ann = loc_ann
            else:
                new_ann = 'other'
        
        old_ann = row.get('annotation', 'other')
        if old_ann != new_ann:
            row['annotation'] = new_ann
            updated += 1
        else:
            row['annotation'] = new_ann  # 确保字段存在
    
    # 确保 annotation 在 fieldnames 中
    if 'annotation' not in fieldnames:
        fieldnames = list(fieldnames) + ['annotation']
    
    # 写回
    with open(variant_info_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    # 统计注释分布
    ann_dist = {}
    for row in rows:
        a = row.get('annotation', 'other')
        ann_dist[a] = ann_dist.get(a, 0) + 1
    print(f'  {g}: {len(rows)}个位点, 更新{updated}个, 分布={ann_dist}')
    total_updated += updated

print(f'\n重新注释完成，共更新 {total_updated} 个位点')
