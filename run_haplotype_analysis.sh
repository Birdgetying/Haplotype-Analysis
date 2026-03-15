#!/bin/bash

# ============================================================================
# 单倍型-表型关联分析平台 (单机版本)
# 基于haplotype_phenotype_analysis.py
# ============================================================================

echo "========================================="
echo "  Haplotype-Phenotype Association Analysis"
echo "========================================="
echo "Start time: $(date)"
echo "Host: $HOSTNAME"

# ============================================================================
# 环境配置
# ============================================================================

echo ""
echo "[Step 1] Checking Python environment..."

# 检查Python
if ! command -v python &> /dev/null; then
    echo "ERROR: Python not found"
    exit 1
fi

# 安装缺失的依赖
echo "[Installing required packages...]"

# 检查并安装必要依赖
for pkg in numpy pandas scipy matplotlib pysam; do
    if python -c "import $pkg" 2>/dev/null; then
        echo "  ✓ $pkg already installed"
    else
        echo "  ✗ $pkg not found, please install manually"
        echo "    conda install $pkg"
        exit 1
    fi
done

# sklearn (scikit-learn) 特殊处理
if python -c "import sklearn" 2>/dev/null; then
    echo "  ✓ scikit-learn already installed"
else
    echo "  ✗ scikit-learn not found, please install manually"
    echo "    conda install scikit-learn"
    exit 1
fi

echo "Python: $(python --version)"


# ============================================================================
# 工作目录
# ============================================================================

echo ""
echo "[Step 2] Setting up working directory..."

cd ~/project1 2>/dev/null || {
    echo "ERROR: Cannot find project directory"
    exit 1
}

echo "Working directory: $(pwd)"


# ============================================================================
# 创建输出目录
# ============================================================================

mkdir -p results/haplotype_analysis
mkdir -p logs

echo "Output directories created"


# ============================================================================
# 数据检查
# ============================================================================

echo ""
echo "[Step 3] Checking data files..."

# 数据路径配置（CS-IAAS T2T v1.1）
VCF_FILE="/storage/public/home/2024110093/data/Variation/CSIAAS/Core819Samples_ALL.vcf.gz"
PHENO_FILE="/storage/public/home/2024110093/data/Variation/CSIAAS/Phe.txt"
GFF_FILE="/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1_HC.gff3"
FASTA_FILE="/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1.fasta"

if [ -f "${VCF_FILE}" ]; then
    echo "  ✓ VCF file found: ${VCF_FILE}"
else
    echo "  ✗ VCF NOT FOUND: ${VCF_FILE}"
    exit 1
fi

if [ -f "${PHENO_FILE}" ]; then
    echo "  ✓ Phenotype file found: ${PHENO_FILE}"
else
    echo "  ✗ Phenotype NOT FOUND: ${PHENO_FILE}"
    exit 1
fi

if [ -f "${GFF_FILE}" ]; then
    echo "  ✓ GFF3 annotation found: ${GFF_FILE}"
else
    echo "  ✗ GFF3 NOT FOUND: ${GFF_FILE}"
    exit 1
fi

if [ -f "${FASTA_FILE}" ]; then
    echo "  ✓ FASTA file found: ${FASTA_FILE}"
else
    echo "  ✗ FASTA NOT FOUND: ${FASTA_FILE}"
    exit 1
fi

# 检查分析脚本
if [ -f "haplotype_phenotype_analysis.py" ]; then
    echo "  ✓ haplotype_phenotype_analysis.py found"
else
    echo "  ✗ haplotype_phenotype_analysis.py - NOT FOUND"
    exit 1
fi


# ============================================================================
# 函数：从 GTF 查找基因坐标
# ============================================================================

get_gene_coords() {
    local GENE_ID=$1
    local GFF=$2
    
    # GFF3格式：匹配 ID=基因ID 或 Parent=基因ID
    # 先尝试匹配 gene 行
    local GENE_ROW=$(grep -E "ID=${GENE_ID}[;\t]" "${GFF}" | awk '$3=="gene"' | head -1)
    
    if [ -n "${GENE_ROW}" ]; then
        echo "[DEBUG] Found gene row" >&2
        echo "${GENE_ROW}" | awk '{print $1"\t"$4"\t"$5"\t"$7}'
        return
    fi
    
    # 如果没有 gene 行，从 mRNA/exon/CDS 计算边界
    local GREP_RESULT=$(grep -E "(ID|Parent)=${GENE_ID}" "${GFF}" | head -5)
    
    if [ -z "${GREP_RESULT}" ]; then
        echo "[DEBUG] No match for ${GENE_ID}" >&2
        return
    fi
    
    echo "[DEBUG] Calculating from mRNA/exon/CDS" >&2
    grep -E "(ID|Parent)=${GENE_ID}" "${GFF}" | awk '
        BEGIN { chr=""; strand=""; min_start=999999999; max_end=0 }
        {
            if (chr=="") { chr=$1; strand=$7 }
            if ($4 < min_start) min_start=$4
            if ($5 > max_end) max_end=$5
        }
        END { if (chr!="") print chr"\t"min_start"\t"max_end"\t"strand }
    '
}


# ============================================================================
# 运行单倍型-表型关联分析
# ============================================================================

echo ""
echo "========================================="
echo "  Running Haplotype-Phenotype Association Analysis"
echo "========================================="
echo "Start: $(date '+%Y-%m-%d %H:%M:%S')"


# ==== 分析测试基因 ====
# di19: CSIAAS1BG1157200HC
# hox: CSIAAS4BG0701800HC

for GENE_ID in "CSIAAS1BG1157200HC" "CSIAAS4BG0701800HC"; do

    echo ""
    echo "-----------------------------------------"
    echo "  Gene: ${GENE_ID}"
    echo "-----------------------------------------"

    # 从 GFF3 获取坐标
    echo "  [DEBUG] Looking up ${GENE_ID} in GFF3..."
    COORDS=$(get_gene_coords "${GENE_ID}" "${GFF_FILE}")
    echo "  [DEBUG] COORDS='${COORDS}'"
    
    if [ -z "${COORDS}" ]; then
        echo "  ERROR: Cannot find ${GENE_ID} in GTF, skipping..."
        continue
    fi

    CHROM=$(echo "${COORDS}" | cut -f1)
    START=$(echo "${COORDS}" | cut -f2)
    END=$(echo "${COORDS}" | cut -f3)
    STRAND=$(echo "${COORDS}" | cut -f4)

    echo "  Coordinates: ${CHROM}:${START}-${END} (${STRAND})"

    # 输出目录按基因分开
    OUT_DIR="./results/${GENE_ID}"
    mkdir -p "${OUT_DIR}"

    LOG_FILE="logs/${GENE_ID}_$(date +%s).log"

    echo "  Output: ${OUT_DIR}"
    echo "  Log: ${LOG_FILE}"
    echo "  Running..."

    python haplotype_phenotype_analysis.py \
        --vcf ${VCF_FILE} \
        --phenotype ${PHENO_FILE} \
        --gtf ${GFF_FILE} \
        --fasta ${FASTA_FILE} \
        --chrom ${CHROM} \
        --start ${START} \
        --end ${END} \
        --gene-id ${GENE_ID} \
        --strand ${STRAND} \
        --output-dir ${OUT_DIR} \
        2>&1 | tee "${LOG_FILE}"

    EXIT_CODE=${PIPESTATUS[0]}

    if [ $EXIT_CODE -eq 0 ]; then
        echo "  ✓ ${GENE_ID} completed successfully!"
        echo "  Results: ls -lh ${OUT_DIR}/"
        ls -lh "${OUT_DIR}/" 2>/dev/null
    else
        echo "  ✗ ${GENE_ID} failed (exit code: $EXIT_CODE)"
        echo "  Last 20 lines of log:"
        tail -20 "${LOG_FILE}"
    fi

done


# ============================================================================
# 完成
# ============================================================================

echo ""
echo "========================================="
echo "  All Analyses Complete"
echo "========================================="
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Results:"
echo "  results/CSIAAS1BG1157200HC/ (di19)"
echo "  results/CSIAAS4BG0701800HC/ (hox)"
