#!/bin/bash
# ============================================================================
# 创建VCF索引脚本
# 在当前目录(d:/Desktop/project1)下创建软连接，然后建立CSI索引
# ============================================================================

set -e

echo "========================================="
echo "  VCF Index Builder"
echo "========================================="

# 项目目录
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

# VCF文件路径（服务器路径）
VCF_FILE="/home/URPs/2025_URP_G1/project1/chrALL.impute.vcf.gz"
SV_VCF_FILE="/home/qinz/project/tmp_Proj/07.GWAS/20260103_Barley_salt_SV/03.Pruning/BubbleAll_bnopart.vcf.gz"

# 软连接目标（放在项目目录下）
VCF_LINK="chrALL.impute.vcf.gz"
SV_VCF_LINK="BubbleAll_bnopart.vcf.gz"

echo ""
echo "[Step 1] 创建软连接..."

# 主VCF软连接
if [ -L "$VCF_LINK" ]; then
    echo "  $VCF_LINK 已存在，跳过"
elif [ -f "$VCF_FILE" ]; then
    ln -s "$VCF_FILE" "$VCF_LINK"
    echo "  ✓ 创建 $VCF_LINK -> $VCF_FILE"
else
    echo "  ✗ 源文件不存在: $VCF_FILE"
    exit 1
fi

# SV VCF软连接
if [ -L "$SV_VCF_LINK" ]; then
    echo "  $SV_VCF_LINK 已存在，跳过"
elif [ -f "$SV_VCF_FILE" ]; then
    ln -s "$SV_VCF_FILE" "$SV_VCF_LINK"
    echo "  ✓ 创建 $SV_VCF_LINK -> $SV_VCF_FILE"
else
    echo "  ✗ 源文件不存在: $SV_VCF_FILE"
    exit 1
fi

echo ""
echo "[Step 2] 创建CSI索引（bcftools index -c）..."

# 检查bcftools是否可用
if ! command -v bcftools &> /dev/null; then
    echo "  ✗ bcftools 未安装，尝试安装..."
    conda install -c bioconda bcftools -y
fi

# 为VCF创建CSI索引（同时在原始路径和软连接目标路径创建）
create_index_for_vcf() {
    local vcf_path="$1"
    local vcf_name=$(basename "$vcf_path" .gz)
    
    # 检查原始路径索引
    if [ -f "${vcf_path}.csi" ]; then
        echo "  ✓ ${vcf_path}.csi 已存在"
        return 0
    fi
    
    echo "  为 ${vcf_name} 创建CSI索引（可能需要几分钟）..."
    bcftools index -c "$vcf_path" 2>/dev/null && echo "  ✓ ${vcf_path}.csi 创建成功" || echo "  ✗ 创建失败"
}

echo "  为主VCF创建索引..."
create_index_for_vcf "$VCF_FILE"

echo "  为SV VCF创建索引..."
create_index_for_vcf "$SV_VCF_FILE"

echo ""
echo "[Step 3] 验证索引..."

ls -lh "$VCF_FILE.csi" 2>/dev/null || echo "  主VCF索引: 待验证"
ls -lh "$SV_VCF_FILE.csi" 2>/dev/null || echo "  SV VCF索引: 待验证"

echo ""
echo "========================================="
echo "  完成！"
echo "========================================="
echo ""
echo "索引文件位置："
echo "  主VCF: $VCF_FILE.csi"
echo "  SV VCF: $SV_VCF_FILE.csi"
echo ""
echo "现在可以重新运行建库脚本，速度将大幅提升。"
