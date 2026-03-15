import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import gzip
import pysam

def parse_gtf(gtf_file, target_gene_id):
    """从 GTF 解析目标基因的外显子和 CDS 位置"""
    exons = []
    cds = []
    gene_strand = '+'
    gene_chrom = None
    
    # 支持 .gz 或纯文本格式
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    with open_func(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            attr = parts[8]
            if target_gene_id not in attr:
                continue

            if gene_chrom is None:
                gene_chrom = parts[0]
            gene_strand = parts[6]

            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])

            if feature == 'exon':
                exons.append((start, end))
            elif feature == 'CDS':
                cds.append((start, end))

    exons = sorted(exons)
    cds = sorted(cds)
    return exons, cds, gene_strand, gene_chrom

def _gt_to_allele(rec, sample_name, ref, alt0):
    """从 pysam 记录的 GT 得到该样本在该位点的等位基因串（REF/ALT/N）。"""
    try:
        gt = rec.samples[sample_name].get("GT")
    except Exception:
        return "N"
    if gt is None or None in gt:
        return "N"
    # 纯合 REF
    if gt == (0, 0) or gt == (0,) or gt == (0, 0, 0):
        return ref if ref else "N"
    # 纯合 ALT (第一个 ALT)
    if gt == (1, 1) or gt == (1,) or gt == (1, 1, 1):
        return alt0 if alt0 else "N"
    # 杂合按 REF 算
    if 0 in gt and 1 in gt:
        return ref if ref else "N"
    return "N"


def parse_vcf_and_get_haps(vcf_file, min_samples=5):
    """
    使用 pysam 流式解析 VCF，不将整个文件读入内存。
    返回:
      - positions: 所有变异的物理位置列表
      - hap_df:    单倍型汇总表 (Haplotype_Seq, Count, Hap_Name, Alleles)
      - hap_sample_df: 单倍型-样本明细表 (SampleID, Hap_Name, Haplotype_Seq)
    """
    vcf = pysam.VariantFile(vcf_file)
    sample_cols = list(vcf.header.samples)

    positions = []
    sample_alleles = {s: [] for s in sample_cols}

    for rec in vcf:
        pos = rec.pos
        ref = rec.ref or ""
        alt0 = (rec.alts[0] if rec.alts else "") or ""
        positions.append(pos)
        for s in sample_cols:
            allele = _gt_to_allele(rec, s, ref, alt0)
            if not allele:
                allele = "N"
            sample_alleles[s].append(allele)

    vcf.close()

    hap_dict = {}
    sample_to_hap = {}

    for sample in sample_cols:
        allele_seq = sample_alleles[sample]
        if "N" in allele_seq:
            continue
        hap_key = "|".join(allele_seq)
        hap_dict[hap_key] = hap_dict.get(hap_key, 0) + 1
        sample_to_hap[sample] = hap_key

    hap_df = pd.DataFrame(list(hap_dict.items()), columns=["Haplotype_Seq", "Count"])
    hap_df = hap_df[hap_df["Count"] >= min_samples].sort_values(by="Count", ascending=False).reset_index(drop=True)
    hap_df["Hap_Name"] = [f"Hap{i+1}" for i in range(len(hap_df))]
    hap_df["Alleles"] = hap_df["Haplotype_Seq"].str.split("|")

    hap_seq_to_name = dict(zip(hap_df["Haplotype_Seq"], hap_df["Hap_Name"]))
    rows = []
    for sample, hap_seq in sample_to_hap.items():
        hap_name = hap_seq_to_name.get(hap_seq, "Other")
        rows.append({"SampleID": sample, "Hap_Name": hap_name, "Haplotype_Seq": hap_seq})
    hap_sample_df = pd.DataFrame(rows, columns=["SampleID", "Hap_Name", "Haplotype_Seq"])

    return positions, hap_df, hap_sample_df


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_codon(codon: str) -> str:
    codon = codon.upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return GENETIC_CODE.get(codon, "X")


def build_coding_context(chrom: str, cds_intervals, strand: str, fasta_path: str):
    """
    根据 CDS 区间构建:
      - cds_seq: 整个 CDS 的参考序列 (按转录本方向拼接)
      - cds_pos_to_idx: 基因组坐标 -> CDS 中的碱基索引 (0-based)
    """
    if not cds_intervals:
        return "", {}

    fasta = pysam.FastaFile(fasta_path)
    cds_intervals = sorted(cds_intervals)
    cds_pos_to_idx = {}
    seq_frags = []
    idx = 0

    if strand == "+":
        for start, end in cds_intervals:
            frag = fasta.fetch(chrom, start - 1, end)
            seq_frags.append(frag)
            for pos in range(start, end + 1):
                cds_pos_to_idx[pos] = idx
                idx += 1
    else:
        # 负链按从高到低坐标遍历，并取反向互补
        for start, end in reversed(cds_intervals):
            frag = fasta.fetch(chrom, start - 1, end)
            frag_rc = revcomp(frag)
            seq_frags.append(frag_rc)
            for pos in range(end, start - 1, -1):
                cds_pos_to_idx[pos] = idx
                idx += 1

    fasta.close()
    cds_seq = "".join(seq_frags)
    return cds_seq, cds_pos_to_idx


def pos_in_any_interval(pos: int, intervals) -> bool:
    for s, e in intervals:
        if s <= pos <= e:
            return True
    return False


def annotate_snp_effects(vcf_file, fasta_path, gene_chrom, cds_intervals, exon_intervals, gene_strand):
    """
    粗略注释 SNP 功能:
      - missense: CDS 内，改变氨基酸
      - synonymous: CDS 内，不改变氨基酸
      - UTR: exon 内但不在 CDS 中
      - other: 其它情况 (intron / intergenic / indel 等)
    返回: {pos -> effect_str}
    """
    effects = {}

    # 若既没有 CDS 也没有 exon，直接返回空
    if not cds_intervals and not exon_intervals:
        return effects

    cds_seq, cds_pos_to_idx = build_coding_context(gene_chrom, cds_intervals, gene_strand, fasta_path)

    fasta = pysam.FastaFile(fasta_path)
    open_func = gzip.open if vcf_file.endswith(".gz") else open
    with open_func(vcf_file, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos_str, _id, ref, alt = parts[:5]
            if chrom != gene_chrom:
                continue
            try:
                pos = int(pos_str)
            except Exception:
                continue

            alts = alt.split(",")
            alt_allele = alts[0]

            # 只考虑简单 SNP
            if len(ref) != 1 or len(alt_allele) != 1:
                effects[pos] = "other"
                continue

            # 先判断是否在 CDS / exon / UTR
            in_cds = pos_in_any_interval(pos, cds_intervals)
            in_exon = pos_in_any_interval(pos, exon_intervals)

            if not in_exon:
                effects[pos] = "other"
                continue

            if in_exon and not in_cds:
                effects[pos] = "UTR"
                continue

            # coding 区，尝试判断是否同义/错义
            if pos not in cds_pos_to_idx or not cds_seq:
                effects[pos] = "other"
                continue

            idx = cds_pos_to_idx[pos]
            codon_start = (idx // 3) * 3
            if codon_start + 3 > len(cds_seq):
                effects[pos] = "other"
                continue

            codon_ref = cds_seq[codon_start : codon_start + 3]

            # alt 碱基按转录本方向修正
            if gene_strand == "+":
                alt_base = alt_allele
            else:
                alt_base = revcomp(alt_allele)[0]

            offset = idx % 3
            codon_alt = codon_ref[:offset] + alt_base + codon_ref[offset + 1 : ]

            aa_ref = translate_codon(codon_ref)
            aa_alt = translate_codon(codon_alt)

            if aa_ref == "X" or aa_alt == "X":
                effects[pos] = "other"
            elif aa_ref == aa_alt:
                effects[pos] = "synonymous"
            else:
                effects[pos] = "missense"

    fasta.close()
    return effects

def _allele_display(allele, indel_style):
    """
    indel_style: 'seq' -> 显示完整序列; 'bp' -> 单碱基照常，插入显示 +Nbp，缺失显示 -Nbp。
    单碱基不变；len>1 视为插入显示 +Nbp（插入长度 = len-1，因首碱基多为 anchor）；len==1 且为 ref 时缺失用 -Nbp 需上游传 ref_len，此处仅做 +Nbp。
    """
    if indel_style == "seq":
        return allele
    if not allele:
        return "N"
    if len(allele) == 1:
        return allele
    return f"+{len(allele)-1}bp"


def plot_haplotype_map(positions, hap_df, exons, cds, region_start, region_end, snp_effects=None,
                       gene_id=None, gene_chrom=None, strand=None, indel_style="seq",
                       output_file="haplotype_plot.pdf"):
    """
    改进版：优化坐标轴尺度，解决基因结构压缩问题
    gene_id, gene_chrom, strand 若提供则显示在基因结构图左上角；SNP 颜色图例显示在右上角。
    indel_style: 'seq' 显示真实序列，'bp' 显示 +Nbp 形式。
    """
    # PDF/PS 中文字为可编辑字符（Adobe 中可选中编辑），非轮廓
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    num_snps = len(positions)
    num_haps = len(hap_df)
    
    # --- 坐标转换逻辑 ---
    # 回到“相对起点”的写法：以 region_start 为 0，单位 kb
    # 然后根据横轴跨度自适应基因结构的高度，保证在任何区间长度下都能看清楚 exon。
    def to_kb(pos):
        return (pos - region_start) / 1000.0
    
    rel_positions_kb = [to_kb(p) for p in positions]
    rel_exons_kb = [(to_kb(s), to_kb(e)) for s, e in exons]

    # X 轴左端固定为 0，右端取 SNP/exon 中的最大相对位置
    rel_start_kb = 0.0
    x_candidates = rel_positions_kb + [e_kb for _, e_kb in rel_exons_kb]
    if x_candidates:
        rel_end_kb = max(x_candidates)
    else:
        rel_end_kb = max((region_end - region_start) / 1000.0, 1.0)

    # 根据横轴长度调整基因结构子图的纵轴范围和 exon 高度（变细、缩小 Y 跨度）
    span = max(rel_end_kb - rel_start_kb, 1e-3)
    gene_half_height = 0.5
    exon_height = 0.28
    cds_height = 0.28
    gene_y_center = -0.22   # 基因结构和变异线整体下移

    # --- 画布设置 ---
    fig = plt.figure(figsize=(max(12, num_snps * 0.5), max(8, num_haps * 0.6 + 4)))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 1.0, num_haps * 0.6], hspace=0.08)
    ax_gene = fig.add_subplot(gs[0])
    ax_line = fig.add_subplot(gs[1])
    ax_table = fig.add_subplot(gs[2])

    effect_colors = {
        "missense": "red",
        "synonymous": "orange",
        "UTR": "purple",
        "other": "gray",
    }
    effect_labels = {
        "missense": "Missense",
        "synonymous": "Synonymous",
        "UTR": "UTR",
        "other": "Other (intron/intergenic)",
    }

    # === 1. 绘制基因结构（变细）===
    ax_gene.set_xlim(rel_start_kb, rel_end_kb)
    ax_gene.set_ylim(-gene_half_height, gene_half_height)
    ax_gene.spines["bottom"].set_visible(False)
    ax_gene.spines["right"].set_visible(False)
    ax_gene.spines["left"].set_visible(False)
    ax_gene.get_yaxis().set_visible(False)
    # 把 Relative Position 坐标放到图上面
    ax_gene.xaxis.tick_top()
    ax_gene.set_xlabel("Relative Position (kb)", fontsize=10, labelpad=6)
    ax_gene.xaxis.set_label_position("top")

    ax_gene.hlines(gene_y_center, rel_start_kb, rel_end_kb, color="gray", linewidth=1.2, zorder=1)
    if rel_exons_kb:
        for s_kb, e_kb in rel_exons_kb:
            w = e_kb - s_kb
            rect = patches.Rectangle(
                (s_kb, gene_y_center - exon_height / 2), w, exon_height,
                linewidth=0.8, edgecolor="navy", facecolor="white", zorder=2,
            )
            ax_gene.add_patch(rect)
    if cds:
        rel_cds_kb = [(to_kb(s), to_kb(e)) for s, e in cds]
        for s_kb, e_kb in rel_cds_kb:
            w = e_kb - s_kb
            rect = patches.Rectangle(
                (s_kb, gene_y_center - cds_height / 2), w, cds_height,
                linewidth=0, edgecolor="navy", facecolor="#5DADE2", zorder=3,
            )
            ax_gene.add_patch(rect)

    # SNP 线直接画在基因结构上（穿过 exon 高度）
    for pos, p_kb in zip(positions, rel_positions_kb):
        effect = snp_effects.get(pos, "other") if snp_effects else "other"
        color = effect_colors.get(effect, "gray")
        ax_gene.vlines(p_kb, gene_y_center - exon_height / 2, gene_y_center + exon_height / 2, color=color, linewidth=1.2, alpha=0.9, zorder=4)

    # 基因信息：左上角
    if gene_id is not None or gene_chrom is not None or strand is not None:
        info_lines = []
        if gene_id is not None:
            info_lines.append(f"Gene: {gene_id}")
        if gene_chrom is not None:
            info_lines.append(f"Chr: {gene_chrom}")
        if strand is not None:
            info_lines.append(f"Strand: {'+' if strand == '+' else '-'}")
        if info_lines:
            ax_gene.text(
                0.02, 0.98, "\n".join(info_lines),
                transform=ax_gene.transAxes,
                fontsize=10,
                verticalalignment="top",
                horizontalalignment="left",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9),
            )

    # SNP 变异类型颜色图例：右上角
    legend_handles = [
        patches.Patch(facecolor=effect_colors[k], edgecolor="black", label=effect_labels[k])
        for k in ["missense", "synonymous", "UTR", "other"]
    ]
    ax_gene.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=8,
        framealpha=0.9,
    )

    # === 2. 绘制连接线，按变异类型上色 ===
    ax_line.set_xlim(rel_start_kb, rel_end_kb)
    ax_line.set_ylim(0, 1)
    ax_line.axis("off")
    table_x_coords = np.linspace(rel_start_kb + rel_end_kb * 0.05, rel_end_kb * 0.95, num_snps)
    for i, p_kb in enumerate(rel_positions_kb):
        tx = table_x_coords[i]
        effect = snp_effects.get(positions[i], "other") if snp_effects else "other"
        color = effect_colors.get(effect, "gray")
        ax_line.plot([p_kb, tx], [1, 0], color=color, linestyle="-", linewidth=0.9, alpha=0.7)

    # === 3. 绘制单倍型序列矩阵 ===
    ax_table.set_xlim(rel_start_kb, rel_end_kb)
    ax_table.set_ylim(-num_haps - 1, 1)
    ax_table.axis('off')
    
    # 颜色映射
    base_colors = {'A': '#2ECC71', 'T': '#E74C3C', 'C': '#3498DB', 'G': '#F1C40F', 'N': '#BDC3C7'}
    
    # 表头：SNP 位置 (这里显示原始 bp 位置，但排布是均匀的)
    for i, pos_bp in enumerate(positions):
        ax_table.text(table_x_coords[i], 0.5, str(pos_bp), ha='center', va='bottom', 
                      rotation=90, fontsize=8, family='monospace')

    for row_idx, row in hap_df.iterrows():
        y_pos = -row_idx - 1
        # 行背景色
        if row_idx % 2 == 0:
            ax_table.axhspan(y_pos - 0.5, y_pos + 0.5, facecolor='#F8F9F9', zorder=-1)
            
        # 单倍型名称和样本数
        ax_table.text(rel_start_kb, y_pos, row['Hap_Name'], ha='left', va='center', fontweight='bold')
        ax_table.text(rel_end_kb, y_pos, f"n={row['Count']}", ha='right', va='center', fontsize=9)
        
        # 序列碱基 / indel：可选真实序列或 +Nbp 显示
        alleles = row.get("Alleles", list(row["Haplotype_Seq"].split("|")))
        for col_idx, allele in enumerate(alleles):
            if col_idx >= len(table_x_coords):
                break
            display_text = _allele_display(allele, indel_style)
            base_key = allele[0] if allele else "N"
            ax_table.text(
                table_x_coords[col_idx],
                y_pos,
                display_text,
                ha="center",
                va="center",
                color=base_colors.get(base_key, "black"),
                fontweight="bold",
                family="monospace",
            )

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ 优化后的单倍型图已生成: {output_file}")

def _parse_args():
    parser = argparse.ArgumentParser(
        description="单倍型区块图：从 VCF 提取单倍型，绘制基因结构 + SNP 功能注释 + 单倍型矩阵。支持 pysam 流式读 VCF，可选 indel 显示方式。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--vcf",
        default="./tmp.vcf",
        help="区间 VCF 文件路径（.vcf 或 .vcf.gz / .bcf）",
    )
    parser.add_argument(
        "--gtf",
        default="/home/qinz/data/genomes/IWGSC_v1/IWGSC_v1.1_HC_LC.gtf",
        help="GTF 注释文件（用于基因结构 exon/CDS）",
    )
    parser.add_argument(
        "--fasta",
        default="/home/qinz/data/genomes/IWGSC_v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta",
        help="参考基因组 FASTA（用于 SNP 功能注释）",
    )
    parser.add_argument(
        "--gene-id",
        default="TraesCS6B02G268300",
        help="目标基因 ID（与 GTF 中 gene_id 匹配）",
    )
    parser.add_argument(
        "--region-start",
        type=int,
        default=482468518,
        help="区间起始位置（bp）",
    )
    parser.add_argument(
        "--region-end",
        type=int,
        default=482472579,
        help="区间终止位置（bp）",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="图中仅展示样本数 >= 该值的单倍型",
    )
    parser.add_argument(
        "--indel-style",
        choices=["seq", "bp"],
        default="seq",
        help="indel 显示方式: seq=真实序列, bp=插入显示 +Nbp（单碱基不变）",
    )
    parser.add_argument(
        "--output",
        default="Gene_Haplotypes.pdf",
        help="输出 PDF 路径",
    )
    parser.add_argument(
        "--summary-tsv",
        default="Haplotypes_summary.tsv",
        help="单倍型汇总表输出路径",
    )
    parser.add_argument(
        "--samples-tsv",
        default="Haplotype_samples.tsv",
        help="单倍型-样本表输出路径",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    VCF_FILE = args.vcf
    GTF_FILE = args.gtf
    FASTA_FILE = args.fasta
    GENE_ID = args.gene_id
    REGION_START = args.region_start
    REGION_END = args.region_end

    print("正在解析 GTF...")
    exons, cds, strand, gene_chrom = parse_gtf(GTF_FILE, GENE_ID)
    if not exons and not cds:
        print(f"⚠️ 注意: 在 GTF 中没有找到基因 {GENE_ID} 的 exon/CDS 注释，基因结构子图将为空或只有一条基线。")
    else:
        print(f"✅ 在 GTF 中找到 {len(exons)} 个 exon, {len(cds)} 个 CDS; 染色体 {gene_chrom}")
        if exons:
            print(f"   exon 坐标范围约为 {exons[0][0]} - {exons[-1][1]}")

    print("正在解析 VCF 鉴定单倍型 (pysam 流式读取)...")
    positions, hap_df, hap_sample_df = parse_vcf_and_get_haps(VCF_FILE, min_samples=args.min_samples)

    hap_df.to_csv(args.summary_tsv, sep="\t", index=False)
    hap_sample_df.to_csv(args.samples_tsv, sep="\t", index=False)
    print(f"✅ 已生成单倍型汇总表: {args.summary_tsv}")
    print(f"✅ 已生成单倍型-样本表: {args.samples_tsv}")

    print("正在注释 SNP 功能 (missense / synonymous / UTR / other)...")
    snp_effects = annotate_snp_effects(
        VCF_FILE,
        FASTA_FILE,
        gene_chrom,
        cds,
        exons,
        strand,
    )
    if snp_effects:
        summary = {k: sum(1 for v in snp_effects.values() if v == k) for k in set(snp_effects.values())}
        print("   SNP 功能注释统计：", summary)
    else:
        print("   警告：没有可注释的 SNP（可能没有 CDS/exon 或坐标不匹配）")

    print("正在生成可视化图表...")
    plot_haplotype_map(
        positions,
        hap_df,
        exons,
        cds,
        REGION_START,
        REGION_END,
        snp_effects=snp_effects,
        gene_id=GENE_ID,
        gene_chrom=gene_chrom,
        strand=strand,
        indel_style=args.indel_style,
        output_file=args.output,
    )

