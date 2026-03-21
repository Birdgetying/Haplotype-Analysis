#!/bin/bash

# ============================================================================
# 单倍型数据集构建 - 生成单倍型数据库
# ============================================================================
# 
# 目录结构:
#   database/                           # 数据库文件夹
#   ├── {gene_id}/
#   │   ├── gene_info.json              # 基因基本信息
#   │   ├── haplotype_data.csv          # 单倍型数据
#   │   ├── haplotype_samples.csv       # 样本-单倍型对应
#   │   ├── haplotype_stats.csv         # 单倍型统计
#   │   ├── phenotype_data.csv          # 表型数据
#   │   └── association_result.csv      # 关联分析结果
#   └── summary.csv                      # 所有基因汇总
#
#   results/                            # 结果文件夹
#   └── {gene_id}/
#       └── integrated_analysis.html    # 综HTML图
#
# 默认处理基因:
#   - CSIAAS1BG1157200HC (di19)
#   - CSIAAS4BG0701800HC (hox)
# ============================================================================

echo "========================================="
echo "  Haplotype Database Builder"
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
DATABASE_DIR="./database"
RESULTS_DIR="./results"

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

mkdir -p ${DATABASE_DIR}
mkdir -p ${RESULTS_DIR}
mkdir -p logs

LOG_FILE="logs/genome_scan_$(date +%s).log"

echo "  Database: ${DATABASE_DIR}"
echo "  Results:  ${RESULTS_DIR}"
echo "  Log: ${LOG_FILE}"
echo ""

# 运行扫描脚本（完整基因区间）
python genome_wide_haplotype_scan.py \
    --vcf ${VCF_FILE} \
    --gff ${GFF_FILE} \
    --phenotype ${PHENO_FILE} \
    --database-dir ${DATABASE_DIR} \
    --results-dir ${RESULTS_DIR} \
    --min-samples 5 \
    --test-region 0 \
    2>&1 | tee "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

# ============================================================================
# 完成
# ============================================================================

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "  Build completed successfully!"
    echo "  Database: ${DATABASE_DIR}/"
    ls -lh ${DATABASE_DIR}/ 2>/dev/null
    echo ""
    echo "  Results: ${RESULTS_DIR}/"
    ls -lh ${RESULTS_DIR}/ 2>/dev/null
else
    echo "  Build failed (exit code: $EXIT_CODE)"
    echo "  Check log: ${LOG_FILE}"
fi
echo "========================================="
echo "End time: $(date)"
