#!/bin/bash

# ============================================================================
# 全基因组单倍型扫描 - 生成单倍型数据库
# ============================================================================

echo "========================================="
echo "  Genome-Wide Haplotype Database Builder"
echo "========================================="
echo "Start time: $(date)"
echo "Host: $HOSTNAME"

# ============================================================================
# 配置
# ============================================================================

# 数据路径 (CS-IAAS T2T v1.1)
VCF_FILE="/storage/public/home/2024110093/data/Variation/CSIAAS/Core819Samples_ALL.vcf.gz"
GFF_FILE="/storage/public/home/2024110093/data/genomes/CS_T2T_v1.1/CS-IAAS_v1.1_HC.gff3"
PHENO_FILE="/storage/public/home/2024110093/data/Variation/CSIAAS/Phe.txt"

# 输出目录
OUTPUT_DIR="./results/genome_scan"

# 指定基因列表（可选，不指定则使用脚本默认值）
# GENES="CSIAAS1BG1157200HC CSIAAS4BG0701800HC"

# ============================================================================
# 环境检查
# ============================================================================

echo ""
echo "[Step 1] Checking environment..."

# 检查Python
if ! command -v python &> /dev/null; then
    echo "ERROR: Python not found"
    exit 1
fi
echo "  ✓ Python: $(python --version)"

# 检查必要模块
for pkg in numpy pandas pysam; do
    if python -c "import $pkg" 2>/dev/null; then
        echo "  ✓ $pkg available"
    else
        echo "  ✗ $pkg not found"
        exit 1
    fi
done

# ============================================================================
# 数据检查
# ============================================================================

echo ""
echo "[Step 2] Checking data files..."

if [ -f "${VCF_FILE}" ]; then
    echo "  ✓ VCF: ${VCF_FILE}"
else
    echo "  ✗ VCF NOT FOUND"
    exit 1
fi

if [ -f "${GFF_FILE}" ]; then
    echo "  ✓ GFF: ${GFF_FILE}"
else
    echo "  ✗ GFF NOT FOUND"
    exit 1
fi

if [ -f "${PHENO_FILE}" ]; then
    echo "  ✓ Phenotype: ${PHENO_FILE}"
else
    echo "  ✗ Phenotype NOT FOUND"
    exit 1
fi

# ============================================================================
# 运行扫描
# ============================================================================

echo ""
echo "[Step 3] Running genome-wide scan..."

cd ~/project1 2>/dev/null || {
    echo "ERROR: Cannot find project directory"
    exit 1
}

mkdir -p ${OUTPUT_DIR}
mkdir -p logs

LOG_FILE="logs/genome_scan_$(date +%s).log"

echo "  Output: ${OUTPUT_DIR}"
echo "  Log: ${LOG_FILE}"
echo ""

# 运行扫描脚本
python genome_wide_haplotype_scan.py \
    --vcf ${VCF_FILE} \
    --gff ${GFF_FILE} \
    --phenotype ${PHENO_FILE} \
    --output-dir ${OUTPUT_DIR} \
    --min-samples 5 \
    2>&1 | tee "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

# ============================================================================
# 完成
# ============================================================================

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "  Scan completed successfully!"
    echo "  Results: ${OUTPUT_DIR}/"
    ls -lh ${OUTPUT_DIR}/ 2>/dev/null
else
    echo "  Scan failed (exit code: $EXIT_CODE)"
    echo "  Check log: ${LOG_FILE}"
fi
echo "========================================="
echo "End time: $(date)"
