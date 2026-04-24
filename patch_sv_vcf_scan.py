"""
修改 genome_wide_haplotype_scan.py: 添加 SV VCF 支持
使用简单的行搜索和修改
"""
with open('genome_wide_haplotype_scan.py', 'r', encoding='utf-8') as f:
    content = f.read()

results = []

# ========== 1. 添加 --sv-vcf 命令行参数 ==========
# 在 args = parser.parse_args() 之前插入
TARGET = '    args = parser.parse_args()'
NEW_ARG = '''    parser.add_argument("--sv-vcf", type=str, default=None,
                        help="结构变异VCF文件路径（如Bubble检测结果）")

'''
if TARGET in content:
    content = content.replace(TARGET, NEW_ARG + TARGET, 1)
    results.append("OK: 添加 --sv-vcf 命令行参数")
else:
    results.append("FAIL: 找不到 'args = parser.parse_args()'")

# ========== 2. 添加 sv_variants.vcf.gz 到 required_files ==========
# 在 'variants.vcf.gz' 行的前面添加新条目
TARGET2 = "'variants.vcf.gz',  # 关键：VCF子集文件"
NEW_ITEM = "            'sv_variants.vcf.gz',  # 结构变异VCF子集\n"
if TARGET2 in content:
    content = content.replace(TARGET2, NEW_ITEM + TARGET2, 1)
    results.append("OK: 添加 sv_variants.vcf.gz 到 required_files")
else:
    # 尝试不带空格的版本
    TARGET2b = "'variants.vcf.gz'"
    if TARGET2b in content:
        idx = content.find(TARGET2b)
        print(f"variants.vcf.gz found at {idx}")
        print(repr(content[idx-30:idx+80]))
    results.append("FAIL: 找不到 'variants.vcf.gz' 关键注释")

# ========== 3. 修改 process_single_gene 签名 ==========
TARGET3 = "                        cophe_files: list = None) -> dict:"
NEW3 = "                        cophe_files: list = None,\n                        sv_vcf_file: str = None) -> dict:"
if TARGET3 in content:
    content = content.replace(TARGET3, NEW3, 1)
    results.append("OK: 修改 process_single_gene 签名")
else:
    # 查找实际内容
    idx = content.find("cophe_files: list = None)")
    if idx >= 0:
        print(f"cophe_files 找到 at {idx}")
        print(repr(content[idx:idx+80]))
    results.append("FAIL: 找不到 process_single_gene 签名 cophe_files")

# ========== 4. 在 extractor.extract_region 后添加 SV VCF 处理 ==========
# 找 extractor 调用后的空行，然后找"# 关键新增"注释
ANCHOR = "        # 关键新增：如果没有变异"
SV_CODE = '''
        # ---- SV VCF 处理：同时从 SV VCF 提取结构变异并合并 ----
        if sv_vcf_file and os.path.exists(sv_vcf_file):
            try:
                print(f"[INFO] {gene_id}: 同时从 SV VCF 提取结构变异: {sv_vcf_file}")
                sv_extractor = HaplotypeExtractor(sv_vcf_file)
                sv_positions, sv_hap_df, sv_hap_sample_df = sv_extractor.extract_region(
                    chrom, start, end, min_samples=min_samples, snp_only=False
                )

                if sv_positions and len(sv_positions) > 0:
                    print(f"[INFO] {gene_id}: 从 SV VCF 提取到 {len(sv_positions)} 个结构变异")

                    # 合并 SV 位置到主位置列表（去重）
                    if positions:
                        pos_set = set(positions)
                        new_sv_pos = [p for p in sv_positions if p not in pos_set]
                        positions = sorted(positions + new_sv_pos)
                    else:
                        positions = list(sv_positions)
                        hap_df = sv_hap_df
                        hap_sample_df = sv_hap_sample_df

                    # 合并 SV 序列到主单倍型序列
                    if sv_hap_sample_df is not None and len(sv_hap_sample_df) > 0:
                        sv_seq_map = dict(zip(sv_hap_sample_df['SampleID'], sv_hap_sample_df['Haplotype_Seq']))
                        new_seqs = []
                        for _, row in hap_sample_df.iterrows():
                            sid = row['SampleID']
                            orig_seq = row['Haplotype_Seq']
                            sv_seq = sv_seq_map.get(sid, '')
                            sv_seq_clean = sv_seq.replace('|', '') if sv_seq else ''
                            new_seqs.append(orig_seq + sv_seq_clean)
                        hap_sample_df = hap_sample_df.copy()
                        hap_sample_df['Haplotype_Seq'] = new_seqs

                        # 更新 hap_df 的 Haplotype_Seq（按 Hap_Name 关联）
                        if hap_df is not None and len(hap_df) > 0:
                            for idx in hap_df.index:
                                hname = hap_df.at[idx, 'Hap_Name']
                                first_sample = hap_sample_df[hap_sample_df['Hap_Name'] == hname]
                                if len(first_sample) > 0:
                                    hap_df.at[idx, 'Haplotype_Seq'] = first_sample.iloc[0]['Haplotype_Seq']

                    # 保存 SV VCF 子集
                    try:
                        sv_subset_path = os.path.join(gene_data_dir, 'sv_variants.vcf.gz')
                        vcf_sample_ids = (hap_sample_df['SampleID'].tolist()
                                         if hap_sample_df is not None else [])
                        create_subset_vcf(sv_vcf_file, chrom, start, end, sv_subset_path,
                                        sample_ids=vcf_sample_ids)
                        sv_size = os.path.getsize(sv_subset_path)
                        if sv_size > 1024:
                            print(f"[INFO] SV VCF子集已保存: {sv_subset_path} ({sv_size/1024:.1f} KB)")
                        else:
                            print(f"[WARNING] SV VCF子集文件过小 ({sv_size} bytes)")
                    except Exception as sv_e:
                        print(f"[WARNING] 保存SV VCF子集失败: {sv_e}")
                else:
                    print(f"[INFO] {gene_id}: SV VCF 在该区域无变异")
            except Exception as sv_ext_e:
                print(f"[WARNING] SV VCF 提取失败: {sv_ext_e}")

'''

if ANCHOR in content:
    content = content.replace(ANCHOR, SV_CODE + ANCHOR, 1)
    results.append("OK: 添加 SV VCF 处理逻辑")
else:
    results.append("FAIL: 找不到 '# 关键新增' 锚点")

# ========== 5. 修改 gene_info_dict ==========
TARGET5 = "            'vcf_mtime': os.path.getmtime(vcf_file),"
NEW5 = "            'sv_vcf_file': sv_vcf_file if sv_vcf_file else None,\n            'vcf_mtime': os.path.getmtime(vcf_file),"
if TARGET5 in content:
    content = content.replace(TARGET5, NEW5, 1)
    results.append("OK: 添加 sv_vcf_file 到 gene_info_dict")
else:
    # 尝试只找 vcf_mtime 前面
    idx = content.find("'vcf_mtime':")
    if idx >= 0:
        print(f"vcf_mtime found at {idx}")
        print(repr(content[idx-50:idx+80]))
    results.append("FAIL: 找不到 gene_info_dict vcf_mtime")

# ========== 6. 修改 process_single_gene 调用 ==========
TARGET6 = "                                     cophe_files=cophe_files)"
NEW6 = "                                     cophe_files=cophe_files,\n                                     sv_vcf_file=args.sv_vcf)"
if TARGET6 in content:
    content = content.replace(TARGET6, NEW6, 1)
    results.append("OK: 修改 process_single_gene 调用，传递 sv_vcf_file")
else:
    # 找调用
    idx = content.find("result = process_single_gene(gene_info, vcf_file, pheno_df,")
    if idx >= 0:
        end_idx = content.find(")", idx + 100)
        print(f"Call found {idx}-{end_idx}")
        print(repr(content[idx:end_idx+10]))
    results.append("FAIL: 找不到 process_single_gene 调用 cophe_files 参数")

# ========== 写回文件 ==========
with open('genome_wide_haplotype_scan.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("=" * 60)
for r in results:
    print(r)
print("=" * 60)
