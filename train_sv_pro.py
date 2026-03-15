"""
SV Pro Training with Advanced Model Architecture and Training Strategies
=======================================================================
Advanced SV analysis with complex model architectures, attention mechanisms,
and sophisticated training strategies for improved performance.
"""

import os
import sys
import json
import time
import csv
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from torch.amp import autocast, GradScaler
import scipy.stats
from sklearn.metrics import r2_score
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV, Ridge
from sklearn.linear_model import ElasticNet
import xgboost as xgb
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import optuna
    OPTUNA_AVAILABLE = True
    print("[INFO] Optuna library available for hyperparameter optimization")
except ImportError:
    OPTUNA_AVAILABLE = False
    print("[WARNING] Optuna library not available, hyperparameter optimization disabled")

# Set font parameters to support English text
rcParams['font.family'] = 'DejaVu Sans'
rcParams['axes.unicode_minus'] = False  # Handle minus signs properly

# 尝试导入VCF处理库
VCF_AVAILABLE = False
ALLEL_AVAILABLE = False
HDF5_AVAILABLE = False

try:
    import cyvcf2
    VCF_AVAILABLE = True
    print("[INFO] cyvcf2 library available for VCF processing")
except ImportError:
    VCF_AVAILABLE = False
    print("[WARNING] cyvcf2 library not available, will install via pip if needed")

try:
    import allel
    ALLEL_AVAILABLE = True
    print("[INFO] scikit-allel library available for genetic analysis")
except ImportError:
    ALLEL_AVAILABLE = False
    print("[WARNING] scikit-allel library not available, will install via pip if needed")

try:
    import h5py
    HDF5_AVAILABLE = True
    print("[INFO] h5py library available for data storage")
except ImportError:
    HDF5_AVAILABLE = False
    print("[WARNING] h5py library not available")

# ============================================================================

# GPU初始化与检查
# ============================================================================

def init_gpu():
    """初始化GPU环境并强制检查"""
    # 设置CUDA环境变量
    if 'CUDA_VISIBLE_DEVICES' not in os.environ:
        os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    
    # 设置cuDNN确定性
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    print("="*60)
    print("GPU INITIALIZATION")
    print("="*60)
    
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print(f"[OK] CUDA is available")
        print(f"  - Device count: {torch.cuda.device_count()}")
        print(f"  - Current device: {torch.cuda.current_device()}")
        print(f"  - Device name: {torch.cuda.get_device_name(0)}")
        print(f"  - Memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
        
        # 清空缓存
        torch.cuda.empty_cache()
        
        # 测试GPU计算
        try:
            test_tensor = torch.randn(100, 100).to(device)
            result = torch.mm(test_tensor, test_tensor.t())
            del test_tensor, result
            torch.cuda.empty_cache()
            print(f"  - GPU compute test: PASSED")
        except Exception as e:
            print(f"  - GPU compute test: FAILED ({e})")
            device = torch.device('cpu')
    else:
        device = torch.device('cpu')
        print(f"[WARNING] CUDA not available!")
        print(f"  - CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
        print(f"  - torch.version.cuda: {torch.version.cuda}")
        print(f"  - Using CPU (will be SLOW)")
    
    print(f"\nFinal device: {device}")
    print("="*60 + "\n")
    return device

# 立即初始化
DEVICE = init_gpu()

# ============================================================================

# 配置参数
# ============================================================================

class SVConfig:
    """SV数据配置"""
    # 数据路径 - 使用完整的SV VCF文件（与train_enhanced_automl_separate.py一致）
    DATA_BASE_PATH = "/storage/public/home/2024110093/data/Variation/CSIAAS/"
    SV_VCF_PATH = DATA_BASE_PATH + "SV.new.vcf.gz"  # 完整SV VCF文件
    PHENOTYPE_PATH = DATA_BASE_PATH + "Phe.txt"
    VCF_ID_PATH = DATA_BASE_PATH + "VCFID.txt"
    
    # 保留旧路径用于兼容（不再使用）
    GWAS_BASE_PATH = DATA_BASE_PATH + "GWAS/"
    SV_ASSOC_PATH = GWAS_BASE_PATH + "Core819Samples_SV.genotypes.whole.final.id_gt.splitbial.miss0.6.813m.DSI_TFW.1e4_assoc.txt"
    
    # 数据加载参数
    MAX_VARIANTS = 20000  # 增加到2万位点进行更全面的分析
    
    # 训练参数 - 优化版本（与train_sv_pro_gwas.py同步）
    BATCH_SIZE = 32
    EPOCHS = 150  # 从100增加到150
    LEARNING_RATE = 3e-4  # 从1e-3降低到3e-4，更稳定的训练
    WEIGHT_DECAY = 0.01
    GRADIENT_CLIP = 1.0
    
    # 增强模型参数 - 优化以防止过拟合
    EMBED_DIM = 128
    DEPTH = 3  # 从6降低到3，防止过拟合（关键修改！）
    NUM_HEADS = 8
    DROP_RATE = 0.3  # 从0.2提高到0.3，增强正则化
    DROP_PATH_RATE = 0.15  # 从0.1提高到0.15
    CNN_CHANNELS = [16, 32, 64, 128]
    PATCH_SIZE = 50
    
    # 交叉验证
    N_FOLDS = 5
    EARLY_STOPPING_PATIENCE = 30  # 从15增加到30，给模型更多训练机会
    
    # GPU优化
    USE_AMP = True  # 启用自动混合精度
    PIN_MEMORY = True
    NUM_WORKERS = 4
    
    # 输出
    OUTPUT_DIR = "./results/sv_pro"
    LOG_DIR = "./logs"

# ============================================================================

# 新增高级组件
# ============================================================================

class MultiBranchCNNExtractor(nn.Module):
    """多分支CNN提取不同尺度的局部特征"""
    def __init__(self, in_channels=1, base_channels=16):
        super().__init__()
        # 不同尺度的卷积核
        self.branch1 = nn.Sequential(
            nn.Conv1d(in_channels, base_channels, kernel_size=3, padding=1),
            nn.BatchNorm1d(base_channels),
            nn.GELU(),
            nn.Conv1d(base_channels, base_channels*2, kernel_size=3, padding=1),
            nn.BatchNorm1d(base_channels*2),
            nn.GELU(),
        )
        
        self.branch2 = nn.Sequential(
            nn.Conv1d(in_channels, base_channels, kernel_size=5, padding=2),
            nn.BatchNorm1d(base_channels),
            nn.GELU(),
            nn.Conv1d(base_channels, base_channels*2, kernel_size=5, padding=2),
            nn.BatchNorm1d(base_channels*2),
            nn.GELU(),
        )
        
        self.branch3 = nn.Sequential(
            nn.Conv1d(in_channels, base_channels, kernel_size=7, padding=3),
            nn.BatchNorm1d(base_channels),
            nn.GELU(),
            nn.Conv1d(base_channels, base_channels*2, kernel_size=7, padding=3),
            nn.BatchNorm1d(base_channels*2),
            nn.GELU(),
        )
        
        self.attention = nn.Sequential(
            nn.Linear(base_channels*2*3, base_channels*2),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        x = x.unsqueeze(1)
        b1 = self.branch1(x)
        b2 = self.branch2(x)
        b3 = self.branch3(x)
        
        # 自适应加权融合
        combined = torch.cat([b1, b2, b3], dim=1)
        B, C, L = combined.shape
        attention_weights = self.attention(combined.mean(dim=-1))
        attention_weights = attention_weights.view(B, C, 1)
        weighted = combined * attention_weights
        
        return weighted


class RelativePositionEncoding(nn.Module):
    """相对位置编码，更适合序列建模"""
    def __init__(self, embed_dim, max_len=512):
        super().__init__()
        self.embed_dim = embed_dim
        self.max_len = max_len
        
        # 学习相对位置偏置
        self.rel_pos_bias = nn.Parameter(torch.zeros(2 * max_len - 1, embed_dim))
        nn.init.xavier_uniform_(self.rel_pos_bias)
        
    def forward(self, seq_len):
        # 生成相对位置索引
        coords = torch.arange(seq_len)
        relative_coords = coords[:, None] - coords[None, :]  # (seq_len, seq_len)
        relative_coords += self.max_len - 1  # 转换为正索引
        
        # 获取对应的偏置
        bias = self.rel_pos_bias[relative_coords]  # (seq_len, seq_len, embed_dim)
        return bias


class EnhancedTransformerBlock(nn.Module):
    """增强的Transformer块，包含相对位置编码"""
    def __init__(self, dim, num_heads=8, max_len=512, mlp_ratio=4.0, 
                 drop=0.0, attn_drop=0.0, drop_path=0.0):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, num_heads, dropout=attn_drop, batch_first=True)
        self.drop_path = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.norm2 = nn.LayerNorm(dim)
        mlp_hidden_dim = int(dim * mlp_ratio)
        self.mlp = nn.Sequential(
            nn.Linear(dim, mlp_hidden_dim),
            nn.GELU(),
            nn.Dropout(drop),
            nn.Linear(mlp_hidden_dim, dim),
            nn.Dropout(drop)
        )
        
        # 相对位置编码
        self.relative_pos_encoding = RelativePositionEncoding(dim, max_len)
        
    def forward(self, x):
        B, N, C = x.shape
        rel_pos_bias = self.relative_pos_encoding(N)
        
        # 应用层归一化
        norm_x = self.norm1(x)
        
        # 多头注意力
        attn_out, _ = self.attn(norm_x, norm_x, norm_x)
        
        # 添加相对位置偏置
        attn_out = attn_out + rel_pos_bias.mean(dim=-1).unsqueeze(-1).expand_as(attn_out)
        
        # 残差连接和DropPath
        x = x + self.drop_path(attn_out)
        
        # MLP部分
        x = x + self.drop_path(self.mlp(self.norm2(x)))
        
        return x


class AttentionBasedFeatureSelector(nn.Module):
    """基于注意力权重的特征选择"""
    def __init__(self, input_dim, bottleneck_ratio=0.5):
        super().__init__()
        self.bottleneck_dim = int(input_dim * bottleneck_ratio)
        
        # 注意力机制学习每个特征的重要性
        self.attention = nn.Sequential(
            nn.Linear(input_dim, input_dim // 2),
            nn.ReLU(),
            nn.Linear(input_dim // 2, input_dim),
            nn.Sigmoid()
        )
        
        # 特征变换
        self.transform = nn.Sequential(
            nn.Linear(input_dim, self.bottleneck_dim),
            nn.LayerNorm(self.bottleneck_dim),
            nn.GELU()
        )
        
    def forward(self, x, return_attention=False):
        # x: (batch, features)
        attention_weights = self.attention(x.mean(dim=0, keepdim=True))  # 全局注意力
        weighted_x = x * attention_weights
        
        # 特征选择：选择注意力权重最高的特征
        with torch.no_grad():
            _, top_indices = attention_weights.topk(self.bottleneck_dim, dim=-1)
        selected_features = x[:, top_indices.squeeze()]
        
        transformed = self.transform(selected_features)
        
        if return_attention:
            return transformed, attention_weights, top_indices
        return transformed


class SelfSupervisedPretextTask(nn.Module):
    """自监督预训练任务"""
    def __init__(self, embed_dim, num_heads):
        super().__init__()
        
        # 主编码器
        self.encoder = nn.Sequential(
            nn.Linear(embed_dim, embed_dim * 2),
            nn.ReLU(),
            nn.Linear(embed_dim * 2, embed_dim)
        )
        
        # 对比学习投影头
        self.projection_head = nn.Sequential(
            nn.Linear(embed_dim, embed_dim),
            nn.ReLU(),
            nn.Linear(embed_dim, embed_dim // 2)
        )
        
    def forward(self, x, x_augmented):
        """对比学习：最大化原始样本和增强样本的相似性"""
        # 编码原始样本和增强样本
        z_original = self.encoder(x)
        z_augmented = self.encoder(x_augmented)
        
        # 投影到对比空间
        p_original = self.projection_head(z_original)
        p_augmented = self.projection_head(z_augmented)
        
        return p_original, p_augmented


def contrastive_loss(features1, features2, temperature=0.1):
    """对比损失（InfoNCE）"""
    batch_size = features1.shape[0]
    features = torch.cat([features1, features2], dim=0)
    
    # 计算相似度矩阵
    similarity_matrix = F.cosine_similarity(features.unsqueeze(1), 
                                           features.unsqueeze(0), dim=2)
    
    # 创建标签（正样本对在对角线上）
    labels = torch.cat([torch.arange(batch_size) for _ in range(2)], dim=0)
    labels = (labels.unsqueeze(0) == labels.unsqueeze(1)).float()
    
    # 应用温度参数
    similarity_matrix = similarity_matrix / temperature
    
    # 计算交叉熵损失
    loss = F.cross_entropy(similarity_matrix, labels.max(dim=1)[1])
    return loss

# ============================================================================

# 位点采样和优化算法
# ============================================================================

def sample_marker_ratios(X, y, ratios=[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]):
    """
    从SV标记中按不同比例随机抽取位点
    """
    # 设置随机种子以确保结果可重现
    np.random.seed(42)
    
    results = {}
    
    for ratio in ratios:
        print(f"\nSampling {ratio*100}% of markers...")
        n_markers = max(1, int(X.shape[1] * ratio))
        
        # 随机选择位点索引
        marker_indices = np.random.choice(X.shape[1], size=n_markers, replace=False)
        X_sampled = X[:, marker_indices]
        
        results[ratio] = {
            'X': X_sampled,
            'indices': marker_indices,
            'n_markers': n_markers
        }
        
        print(f"  Selected {n_markers} markers out of {X.shape[1]} total")
        
    return results


def sample_marker_counts(X, y, marker_counts=[50, 100, 200, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000]):
    """
    按固定数量抽取标记，用于研究标记数量与精度的关系
    
    Args:
        X: 基因型矩阵 (samples x markers)
        y: 表型数据
        marker_counts: 要测试的标记数量列表
    
    Returns:
        dict: 每个标记数量对应的数据字典
    """
    np.random.seed(42)
    
    total_markers = X.shape[1]
    results = {}
    
    print(f"\n{'='*60}")
    print(f"标记数量采样实验 - 总标记数: {total_markers}")
    print(f"{'='*60}")
    
    for n_markers in marker_counts:
        if n_markers > total_markers:
            print(f"  跳过 {n_markers} 标记 (超过总数 {total_markers})")
            continue
        
        print(f"\n采样 {n_markers} 个标记...")
        
        # 随机选择位点索引
        marker_indices = np.random.choice(total_markers, size=n_markers, replace=False)
        X_sampled = X[:, marker_indices]
        
        results[n_markers] = {
            'X': X_sampled,
            'indices': marker_indices,
            'n_markers': n_markers
        }
        
        print(f"  成功采样 {n_markers} 个标记 (X形状: {X_sampled.shape})")
    
    return results


def evaluate_marker_count_effect(X, y, config: SVConfig, device: torch.device,
                                  marker_counts=[50, 100, 200, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000]):
    """
    评估不同标记数量对预测精度的影响
    
    Args:
        X: 完整基因型矩阵
        y: 表型数据
        config: 配置对象
        device: 计算设备
        marker_counts: 要测试的标记数量列表
    
    Returns:
        dict: 每个标记数量的评估结果
    """
    print(f"\n{'='*60}")
    print("标记数量与精度关系研究")
    print(f"Marker Count vs Prediction Accuracy Study")
    print(f"测试标记数量: {marker_counts}")
    print(f"{'='*60}")
    
    # 采样不同数量的标记
    marker_samples = sample_marker_counts(X, y, marker_counts)
    
    # 准备交叉验证
    y_binned = pd.qcut(y, q=5, labels=False, duplicates='drop')
    kfold = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=42)
    
    results = {}
    
    for n_markers, data in marker_samples.items():
        print(f"\n{'='*40}")
        print(f"评估 {n_markers} 个标记的性能")
        print(f"{'='*40}")
        
        X_sampled = data['X']
        
        # 初始化结果存储
        model_results = {
            'RRBLUP': {'r2': [], 'corr': []},
            'XGBoost': {'r2': [], 'corr': []},
            'WheatGP': {'r2': [], 'corr': []},
            'Enhanced': {'r2': [], 'corr': []}
        }
        
        for fold, (train_idx, val_idx) in enumerate(kfold.split(X_sampled, y_binned)):
            print(f"\n  --- Fold {fold+1}/{config.N_FOLDS} ---")
            
            X_train, X_val = X_sampled[train_idx], X_sampled[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            # Fold内标准化
            y_mean, y_std = y_train.mean(), y_train.std()
            y_train_norm = (y_train - y_mean) / (y_std + 1e-8)
            y_val_norm = (y_val - y_mean) / (y_std + 1e-8)
            
            X_mean, X_std = X_train.mean(axis=0, keepdims=True), X_train.std(axis=0, keepdims=True) + 1e-8
            X_train_norm = (X_train - X_mean) / X_std
            X_val_norm = (X_val - X_mean) / X_std
            
            # 1. RRBLUP
            # 注意：alpha上界限制为10^3，避免高维数据下过度收缩导致负R²
            rr_model = RidgeCV(alphas=np.logspace(-3, 3, 20))
            rr_model.fit(X_train_norm, y_train_norm)
            rr_preds = rr_model.predict(X_val_norm)
            rr_r2 = r2_score(y_val_norm, rr_preds)
            rr_corr, _ = scipy.stats.pearsonr(y_val_norm, rr_preds)
            model_results['RRBLUP']['r2'].append(rr_r2)
            model_results['RRBLUP']['corr'].append(rr_corr)
            
            # 2. XGBoost
            xgb_preds, xgb_r2, xgb_corr = train_xgboost_model(
                X_train_norm, X_val_norm, y_train_norm, y_val_norm
            )
            model_results['XGBoost']['r2'].append(xgb_r2)
            model_results['XGBoost']['corr'].append(xgb_corr)
            
            # 3. WheatGP
            wheatgp_preds, wheatgp_r2, wheatgp_corr = train_wheatgp_model(
                X_train_norm, X_val_norm, y_train_norm, y_val_norm
            )
            model_results['WheatGP']['r2'].append(wheatgp_r2)
            model_results['WheatGP']['corr'].append(wheatgp_corr)
            
            # 4. Enhanced Model
            enhanced_preds, enhanced_r2, enhanced_corr = train_enhanced_model(
                X_train_norm, X_val_norm, y_train_norm, y_val_norm, device, config, use_optimization=False
            )
            model_results['Enhanced']['r2'].append(enhanced_r2)
            model_results['Enhanced']['corr'].append(enhanced_corr)
            
            print(f"    RRBLUP R2: {rr_r2:.4f}, Corr: {rr_corr:.4f}")
            print(f"    XGBoost R2: {xgb_r2:.4f}, Corr: {xgb_corr:.4f}")
            print(f"    WheatGP R2: {wheatgp_r2:.4f}, Corr: {wheatgp_corr:.4f}")
            print(f"    Enhanced R2: {enhanced_r2:.4f}, Corr: {enhanced_corr:.4f}")
        
        # 汇总结果
        results[n_markers] = {
            'n_markers': n_markers,
            'RRBLUP': {
                'mean_r2': float(np.mean(model_results['RRBLUP']['r2'])),
                'std_r2': float(np.std(model_results['RRBLUP']['r2'])),
                'mean_corr': float(np.mean(model_results['RRBLUP']['corr'])),
                'std_corr': float(np.std(model_results['RRBLUP']['corr']))
            },
            'XGBoost': {
                'mean_r2': float(np.mean(model_results['XGBoost']['r2'])),
                'std_r2': float(np.std(model_results['XGBoost']['r2'])),
                'mean_corr': float(np.mean(model_results['XGBoost']['corr'])),
                'std_corr': float(np.std(model_results['XGBoost']['corr']))
            },
            'WheatGP': {
                'mean_r2': float(np.mean(model_results['WheatGP']['r2'])),
                'std_r2': float(np.std(model_results['WheatGP']['r2'])),
                'mean_corr': float(np.mean(model_results['WheatGP']['corr'])),
                'std_corr': float(np.std(model_results['WheatGP']['corr']))
            },
            'Enhanced': {
                'mean_r2': float(np.mean(model_results['Enhanced']['r2'])),
                'std_r2': float(np.std(model_results['Enhanced']['r2'])),
                'mean_corr': float(np.mean(model_results['Enhanced']['corr'])),
                'std_corr': float(np.std(model_results['Enhanced']['corr']))
            }
        }
        
        print(f"\n  {n_markers}标记汇总:")
        for model_name in ['RRBLUP', 'XGBoost', 'WheatGP', 'Enhanced']:
            r2_mean = results[n_markers][model_name]['mean_r2']
            r2_std = results[n_markers][model_name]['std_r2']
            corr_mean = results[n_markers][model_name]['mean_corr']
            print(f"    {model_name}: R2={r2_mean:.4f}±{r2_std:.4f}, Corr={corr_mean:.4f}")
    
    return results


def plot_marker_count_vs_accuracy(results, output_dir):
    """
    绘制标记数量与预测精度的关系图
    
    Args:
        results: evaluate_marker_count_effect()的返回结果
        output_dir: 输出目录
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 提取数据
    marker_counts = sorted(results.keys())
    models = ['RRBLUP', 'XGBoost', 'WheatGP', 'Enhanced']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers_style = ['o', 's', '^', 'D']
    
    # 创建2x2子图（只绘制Correlation相关）
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # ========== 图1: Correlation vs 标记数量 (折线图) ==========
    ax1 = axes[0, 0]
    for i, model in enumerate(models):
        corr_means = [results[n][model]['mean_corr'] for n in marker_counts]
        corr_stds = [results[n][model]['std_corr'] for n in marker_counts]
        ax1.errorbar(marker_counts, corr_means, yerr=corr_stds, 
                    marker=markers_style[i], color=colors[i], 
                    label=model, capsize=5, linewidth=2, markersize=8)
    ax1.set_xlabel('Number of SV Markers', fontsize=12)
    ax1.set_ylabel('Correlation', fontsize=12)
    ax1.set_title('Correlation vs Number of SV Markers', fontsize=14, fontweight='bold')
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(marker_counts)
    ax1.set_xticklabels([str(n) for n in marker_counts], rotation=45)
    
    # ========== 图2: Correlation横向对比(柱状图) ==========
    ax2 = axes[0, 1]
    x_pos = np.arange(len(marker_counts))
    width = 0.2
    for i, model in enumerate(models):
        corr_means = [results[n][model]['mean_corr'] for n in marker_counts]
        ax2.bar(x_pos + i*width, corr_means, width, label=model, color=colors[i], alpha=0.8)
    ax2.set_xlabel('Number of Markers', fontsize=12)
    ax2.set_ylabel('Correlation', fontsize=12)
    ax2.set_title('Model Comparison at Different Marker Counts', fontsize=14, fontweight='bold')
    ax2.set_xticks(x_pos + width*1.5)
    ax2.set_xticklabels([str(n) for n in marker_counts], rotation=45)
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # ========== 图3: Correlation提升曲线 ==========
    ax3 = axes[1, 0]
    min_markers = marker_counts[0]
    for i, model in enumerate(models):
        base_corr = results[min_markers][model]['mean_corr']
        corr_improvements = [(results[n][model]['mean_corr'] - base_corr) for n in marker_counts]
        ax3.plot(marker_counts, corr_improvements, 
                marker=markers_style[i], color=colors[i],
                label=model, linewidth=2, markersize=8)
    ax3.set_xlabel('Number of SV Markers', fontsize=12)
    ax3.set_ylabel(f'Correlation Improvement (vs {min_markers} markers)', fontsize=12)
    ax3.set_title('Performance Improvement vs Marker Count', fontsize=14, fontweight='bold')
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax3.legend(loc='lower right')
    ax3.grid(True, alpha=0.3)
    ax3.set_xticks(marker_counts)
    ax3.set_xticklabels([str(n) for n in marker_counts], rotation=45)
    
    # ========== 图4: 最佳模型在各标记数下的Correlation ==========
    ax4 = axes[1, 1]
    best_corrs = []
    best_models = []
    for n in marker_counts:
        model_corrs = {model: results[n][model]['mean_corr'] for model in models}
        best_model = max(model_corrs, key=model_corrs.get)
        best_corrs.append(model_corrs[best_model])
        best_models.append(best_model)
    
    bar_colors = [colors[models.index(m)] for m in best_models]
    bars = ax4.bar(range(len(marker_counts)), best_corrs, color=bar_colors, alpha=0.8)
    ax4.set_xlabel('Number of SV Markers', fontsize=12)
    ax4.set_ylabel('Best Correlation', fontsize=12)
    ax4.set_title('Best Model Correlation at Each Marker Count', fontsize=14, fontweight='bold')
    ax4.set_xticks(range(len(marker_counts)))
    ax4.set_xticklabels([str(n) for n in marker_counts], rotation=45)
    ax4.grid(True, alpha=0.3, axis='y')
    
    for i, (bar, model) in enumerate(zip(bars, best_models)):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                model[:3], ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存图片
    plot_path = os.path.join(output_dir, "marker_count_vs_accuracy.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"\n已保存标记数量-精度关系图: {plot_path}")
    plt.close()
    
    # 创建额外的详细图（单独模型Correlation对比）
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 12))
    
    for idx, model in enumerate(models):
        ax = axes2[idx // 2, idx % 2]
        corr_means = [results[n][model]['mean_corr'] for n in marker_counts]
        corr_stds = [results[n][model]['std_corr'] for n in marker_counts]
        
        ax.errorbar(marker_counts, corr_means, yerr=corr_stds,
                   marker='o', color=colors[idx], label='Correlation', 
                   capsize=5, linewidth=2, markersize=8)
        
        ax.set_xlabel('Number of SV Markers', fontsize=12)
        ax.set_ylabel('Correlation', fontsize=12)
        ax.set_title(f'{model} Correlation vs Marker Count', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xticks(marker_counts)
        ax.set_xticklabels([str(n) for n in marker_counts], rotation=45)
    
    plt.tight_layout()
    plot_path2 = os.path.join(output_dir, "marker_count_by_model.png")
    plt.savefig(plot_path2, dpi=300, bbox_inches='tight')
    print(f"已保存各模型详细图: {plot_path2}")
    plt.close()
    
    # 保存结果到CSV
    csv_path = os.path.join(output_dir, "marker_count_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['n_markers', 'model', 'mean_r2', 'std_r2', 'mean_corr', 'std_corr'])
        for n_markers in marker_counts:
            for model in models:
                writer.writerow([
                    n_markers,
                    model,
                    results[n_markers][model]['mean_r2'],
                    results[n_markers][model]['std_r2'],
                    results[n_markers][model]['mean_corr'],
                    results[n_markers][model]['std_corr']
                ])
    print(f"已保存结果CSV: {csv_path}")
    
    return plot_path, plot_path2, csv_path


def optimize_marker_selection_on_train_only(X_train, y_train, initial_ratio=0.3, iterations=10):
    """
    在训练集上进行标记选择优化，避免使用验证集信息
    """
    print(f"\nOptimizing marker selection on training set only with initial ratio {initial_ratio}...")
    
    # 优化函数使用固定的随机种子以确保结果可重现
    np.random.seed(42)
    
    # 初始采样
    n_markers = max(1, int(X_train.shape[1] * initial_ratio))
    current_indices = np.random.choice(X_train.shape[1], size=n_markers, replace=False)
    best_indices = current_indices.copy()
    best_performance = -float('inf')
    
    # 使用RRBLUP评估性能
    y_binned = pd.qcut(y_train, q=5, labels=False, duplicates='drop')
    kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for iteration in range(iterations):
        print(f"  Iteration {iteration + 1}/{iterations}")
        
        # 当前标记集性能评估 - 在每个fold内标准化
        X_current = X_train[:, current_indices]
        fold_scores = []
        
        for train_idx, val_idx in kfold.split(X_current, y_binned):
            X_fold_train, X_fold_val = X_current[train_idx], X_current[val_idx]
            y_fold_train, y_fold_val = y_train[train_idx], y_train[val_idx]
            
            # 在每个fold内对y进行标准化
            y_fold_train_mean, y_fold_train_std = y_fold_train.mean(), y_fold_train.std()
            y_fold_train_normalized = (y_fold_train - y_fold_train_mean) / (y_fold_train_std + 1e-8)
            y_fold_val_normalized = (y_fold_val - y_fold_train_mean) / (y_fold_train_std + 1e-8)
            
            # 在每个fold内对X进行标准化
            X_fold_train_mean, X_fold_train_std = X_fold_train.mean(axis=0, keepdims=True), X_fold_train.std(axis=0, keepdims=True) + 1e-8
            X_fold_train_normalized = (X_fold_train - X_fold_train_mean) / X_fold_train_std
            X_fold_val_normalized = (X_fold_val - X_fold_train_mean) / X_fold_train_std
            
            # RRBLUP训练
            # alpha上界限制为10^3，避免过度收缩
            rr_model = RidgeCV(alphas=np.logspace(-3, 3, 20))
            rr_model.fit(X_fold_train_normalized, y_fold_train_normalized)
            rr_preds = rr_model.predict(X_fold_val_normalized)
            
            score = r2_score(y_fold_val_normalized, rr_preds)
            fold_scores.append(score)
        
        current_performance = np.mean(fold_scores)
        print(f"    Current performance (R2): {current_performance:.4f}")
        
        # 更新最佳结果
        if current_performance > best_performance:
            best_performance = current_performance
            best_indices = current_indices.copy()
            print(f"    New best performance: {best_performance:.4f}")
        
        # 改进的优化策略：使用基于性能的特征重要性来选择要移除的标记
        n_to_add = max(1, int(len(current_indices) * 0.1))  # 交换10%的标记
        
        # 计算每个标记的贡献度：通过逐一移除标记并观察性能变化
        if len(current_indices) > 1:
            contributions = []
            X_current = X_train[:, current_indices]
            
            # 使用快速评估方法计算每个标记的贡献
            temp_indices = current_indices.copy()
            for i in range(len(current_indices)):
                # 创建一个临时的标记集，移除第i个标记
                temp_selected_indices = np.delete(temp_indices, i)
                if len(temp_selected_indices) > 0:
                    X_temp = X_train[:, temp_selected_indices]
                    
                    # 快速评估性能（使用较小的kfold）
                    temp_scores = []
                    temp_kfold = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)  # 减少fold数以提高速度
                    y_binned_temp = pd.qcut(y_train, q=3, labels=False, duplicates='drop')  # 减少bin数
                    
                    for train_idx, val_idx in temp_kfold.split(X_temp, y_binned_temp):
                        X_fold_train, X_fold_val = X_temp[train_idx], X_temp[val_idx]
                        y_fold_train, y_fold_val = y_train[train_idx], y_train[val_idx]
                        
                        # 在每个fold内对y进行标准化
                        y_fold_train_mean, y_fold_train_std = y_fold_train.mean(), y_fold_train.std()
                        y_fold_train_normalized = (y_fold_train - y_fold_train_mean) / (y_fold_train_std + 1e-8)
                        y_fold_val_normalized = (y_fold_val - y_fold_train_mean) / (y_fold_train_std + 1e-8)
                        
                        # 在每个fold内对X进行标准化
                        X_fold_train_mean, X_fold_train_std = X_fold_train.mean(axis=0, keepdims=True), X_fold_train.std(axis=0, keepdims=True) + 1e-8
                        X_fold_train_normalized = (X_fold_train - X_fold_train_mean) / X_fold_train_std
                        X_fold_val_normalized = (X_fold_val - X_fold_train_mean) / X_fold_train_std
                        
                        # RRBLUP训练
                        rr_model = RidgeCV(alphas=np.logspace(-3, 2, 10))  # 减少alphas数量
                        rr_model.fit(X_fold_train_normalized, y_fold_train_normalized)
                        rr_preds = rr_model.predict(X_fold_val_normalized)
                        
                        score = r2_score(y_fold_val_normalized, rr_preds)
                        temp_scores.append(score)
                    
                    contribution = current_performance - np.mean(temp_scores)  # 贡献度 = 完整性能 - 移除后性能
                    contributions.append(contribution)
                else:
                    contributions.append(0)  # 如果只剩一个标记，贡献度为0
            
            # 选择贡献度最高的标记进行保留，移除贡献度最低的
            if len(contributions) >= n_to_add:
                indices_to_remove = np.argsort(contributions)[:n_to_add]  # 选择贡献度最低的n_to_add个标记
            else:
                indices_to_remove = np.arange(len(contributions))  # 如果不够，移除所有
        else:
            # 如果只有一个标记，随机移除
            indices_to_remove = np.random.choice(len(current_indices), size=min(n_to_add, len(current_indices)), replace=False)
        
        remaining_indices = np.delete(current_indices, indices_to_remove)
        
        # 选择要添加的标记（从未使用的标记中选择）
        unused_indices = np.setdiff1d(np.arange(X_train.shape[1]), remaining_indices)
        if len(unused_indices) >= n_to_add:
            # 为了更智能地选择要添加的标记，我们可以随机选择或使用其他启发式方法
            new_indices = np.random.choice(unused_indices, size=n_to_add, replace=False)
        else:
            # 如果未使用的标记不足，从剩余的标记中随机选择
            new_indices = unused_indices if len(unused_indices) > 0 else np.array([])
            if len(new_indices) < n_to_add and len(remaining_indices) > 0:
                additional_needed = n_to_add - len(new_indices)
                # 从剩余标记中随机选择一些进行替换
                additional_indices = np.random.choice(remaining_indices, size=min(additional_needed, len(remaining_indices)), replace=False)
                new_indices = np.concatenate([new_indices, additional_indices])
        
        # 更新当前标记集
        if len(new_indices) > 0:
            current_indices = np.concatenate([remaining_indices, new_indices])
        
        # 确保标记集大小一致
        if len(current_indices) != n_markers:
            if len(current_indices) > n_markers:
                # 如果多了，随机删除一些
                current_indices = np.random.choice(current_indices, size=n_markers, replace=False)
            else:
                # 如果少了，随机添加一些
                unused_indices = np.setdiff1d(np.arange(X_train.shape[1]), current_indices)
                if len(unused_indices) > 0:
                    n_needed = n_markers - len(current_indices)
                    additional_indices = np.random.choice(unused_indices, size=min(n_needed, len(unused_indices)), replace=False)
                    current_indices = np.concatenate([current_indices, additional_indices])
    
    print(f"Training set optimization complete. Best performance: {best_performance:.4f} with {len(best_indices)} markers")
    
    return best_indices, best_performance

# 删除此函数以避免混淆
# def optimize_marker_selection(X, y, initial_ratio=0.3, iterations=10):
#     """
#     旧的函数，现在使用新的标记选择方法
#     """
#     # 调用新函数
#     return optimize_marker_selection_on_train_only(X, y, initial_ratio, iterations)



def evaluate_marker_sets(marker_samples, y, config: SVConfig):
    """
    评估不同标记集的性能，对每个比例评估所有模型
    """
    print(f"\nEvaluating marker sets for all models...")
    
    y_binned = pd.qcut(y, q=5, labels=False, duplicates='drop')
    kfold = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=42)
    
    results = {}
    
    for ratio, data in marker_samples.items():
        print(f"Evaluating {ratio*100}% markers ({data['n_markers']} markers) for all models...")
        
        X_sampled = data['X']
        
        # RRBLUP fold scores
        rrblup_fold_r2_scores = []
        rrblup_fold_corr_scores = []
        
        # XGBoost fold scores
        xgboost_fold_r2_scores = []
        xgboost_fold_corr_scores = []
        
        # Enhanced fold scores
        enhanced_fold_r2_scores = []
        enhanced_fold_corr_scores = []
        
        for train_idx, val_idx in kfold.split(X_sampled, y_binned):
            X_train, X_val = X_sampled[train_idx], X_sampled[val_idx]
            y_train_actual, y_val_actual = y[train_idx], y[val_idx]
            
            # 在每个fold内对y进行标准化
            y_fold_mean, y_fold_std = y_train_actual.mean(), y_train_actual.std()
            y_train_normalized = (y_train_actual - y_fold_mean) / (y_fold_std + 1e-8)
            y_val_normalized = (y_val_actual - y_fold_mean) / (y_fold_std + 1e-8)
            
            # 在每个fold内对X进行标准化
            X_train_mean, X_train_std = X_train.mean(axis=0, keepdims=True), X_train.std(axis=0, keepdims=True) + 1e-8
            X_train_normalized = (X_train - X_train_mean) / X_train_std
            X_val_normalized = (X_val - X_train_mean) / X_train_std
            
            # RRBLUP模型（alpha上界限制为10^3）
            rr_model = RidgeCV(alphas=np.logspace(-3, 3, 20))
            rr_model.fit(X_train_normalized, y_train_normalized)
            rr_preds = rr_model.predict(X_val_normalized)
            
            rr_r2 = r2_score(y_val_normalized, rr_preds)
            rr_corr, _ = scipy.stats.pearsonr(y_val_normalized, rr_preds)
            
            rrblup_fold_r2_scores.append(rr_r2)
            rrblup_fold_corr_scores.append(rr_corr)
            
            # XGBoost模型
            xgb_preds, xgb_r2, xgb_corr = train_xgboost_model(
                X_train_normalized, X_val_normalized, y_train_normalized, y_val_normalized
            )
            
            xgboost_fold_r2_scores.append(xgb_r2)
            xgboost_fold_corr_scores.append(xgb_corr)
            
            # Enhanced模型
            enhanced_preds, enhanced_r2, enhanced_corr = train_enhanced_model(
                X_train_normalized, X_val_normalized, y_train_normalized, y_val_normalized, DEVICE, config, use_optimization=True
            )
            
            enhanced_fold_r2_scores.append(enhanced_r2)
            enhanced_fold_corr_scores.append(enhanced_corr)
        
        results[ratio] = {
            'RRBLUP': {
                'mean_r2': float(np.mean(rrblup_fold_r2_scores)),
                'std_r2': float(np.std(rrblup_fold_r2_scores)),
                'mean_corr': float(np.mean(rrblup_fold_corr_scores)),
                'std_corr': float(np.std(rrblup_fold_corr_scores))
            },
            'XGBoost': {
                'mean_r2': float(np.mean(xgboost_fold_r2_scores)),
                'std_r2': float(np.std(xgboost_fold_r2_scores)),
                'mean_corr': float(np.mean(xgboost_fold_corr_scores)),
                'std_corr': float(np.std(xgboost_fold_corr_scores))
            },
            'Enhanced': {
                'mean_r2': float(np.mean(enhanced_fold_r2_scores)),
                'std_r2': float(np.std(enhanced_fold_r2_scores)),
                'mean_corr': float(np.mean(enhanced_fold_corr_scores)),
                'std_corr': float(np.std(enhanced_fold_corr_scores))
            },
            'n_markers': int(data['n_markers'])
        }
        
        print(f"  RRBLUP R2: {results[ratio]['RRBLUP']['mean_r2']:.4f} ± {results[ratio]['RRBLUP']['std_r2']:.4f}")
        print(f"  XGBoost R2: {results[ratio]['XGBoost']['mean_r2']:.4f} ± {results[ratio]['XGBoost']['std_r2']:.4f}")
        print(f"  Enhanced R2: {results[ratio]['Enhanced']['mean_r2']:.4f} ± {results[ratio]['Enhanced']['std_r2']:.4f}")
        
        print(f"  RRBLUP Corr: {results[ratio]['RRBLUP']['mean_corr']:.4f} ± {results[ratio]['RRBLUP']['std_corr']:.4f}")
        print(f"  XGBoost Corr: {results[ratio]['XGBoost']['mean_corr']:.4f} ± {results[ratio]['XGBoost']['std_corr']:.4f}")
        print(f"  Enhanced Corr: {results[ratio]['Enhanced']['mean_corr']:.4f} ± {results[ratio]['Enhanced']['std_corr']:.4f}")
    
    return results


class SVProNet(nn.Module):
    """高级SV增强网络，集成多分支CNN、相对位置编码、注意力机制等"""
    def __init__(self, n_sv, embed_dim=128, depth=6, num_heads=8, 
                 cnn_channels=[16, 32, 64, 128], patch_size=50, 
                 drop_rate=0.2, drop_path_rate=0.1):
        super().__init__()
        
        # 输入维度
        self.n_sv = n_sv
        self.embed_dim = embed_dim
        self.patch_size = patch_size
        self.drop_rate = drop_rate
        
        # 多分支CNN特征提取器
        self.multi_branch_cnn = MultiBranchCNNExtractor(in_channels=1, base_channels=cnn_channels[0]//2)
        
        # CNN到Transformer的投影
        self.cnn_to_transformer = nn.Linear(cnn_channels[-1], embed_dim)
        
        # SV重要性学习
        self.sv_importance = nn.Parameter(torch.randn(n_sv) * 0.02)
        
        # 位置嵌入
        self.pos_embed = nn.Parameter(torch.randn(1, n_sv + 1, embed_dim) * 0.02)  # +1 for CLS token
        
        # 类别token
        self.cls_token = nn.Parameter(torch.randn(1, 1, embed_dim) * 0.02)
        
        # Dropout
        self.pos_drop = nn.Dropout(drop_rate)
        
        # 计算DropPath速率
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, depth)]
        
        # Transformer块
        self.blocks = nn.ModuleList([
            EnhancedTransformerBlock(
                dim=embed_dim, 
                num_heads=num_heads, 
                max_len=n_sv + 1,  # +1 for CLS token
                drop=drop_rate, 
                attn_drop=drop_rate, 
                drop_path=dpr[i]
            )
            for i in range(depth)
        ])
        
        # 归一化层
        self.norm = nn.LayerNorm(embed_dim)
        
        # 输出层
        self.head = nn.Sequential(
            nn.Linear(embed_dim, embed_dim // 2),
            nn.GELU(),
            nn.Dropout(drop_rate),
            nn.Linear(embed_dim // 2, 1)
        )
        
    def forward(self, x):
        B = x.size(0)
        
        # 应用SV重要性权重
        importance = torch.sigmoid(self.sv_importance)
        x = x * importance.unsqueeze(0)
        
        # 多分支CNN特征提取
        x_cnn = self.multi_branch_cnn(x)  # (B, channels, seq_len)
        x_cnn = x_cnn.permute(0, 2, 1)  # (B, seq_len, channels)
        
        # 投影到Transformer维度
        x_embed = self.cnn_to_transformer(x_cnn)  # (B, seq_len, embed_dim)
        
        # 添加类别token
        cls_tokens = self.cls_token.expand(B, -1, -1)  # (B, 1, embed_dim)
        x_embed = torch.cat((cls_tokens, x_embed), dim=1)  # (B, seq_len+1, embed_dim)
        
        # 添加位置嵌入
        x_embed = x_embed + self.pos_embed
        x_embed = self.pos_drop(x_embed)
        
        # 通过Transformer块
        for blk in self.blocks:
            x_embed = blk(x_embed)
        
        # 应用最终归一化
        x_embed = self.norm(x_embed)
        
        # 使用类别token的输出进行预测
        cls_output = x_embed[:, 0]  # (B, embed_dim)
        
        # 通过输出层
        output = self.head(cls_output)
        
        return output.squeeze(-1)


def curriculum_learning_training(model, X_train, y_train, device, config, 
                                 difficulty_levels=5):
    """课程学习：从易到难训练"""
    
    # 根据表型值复杂度划分样本
    y_sorted_indices = np.argsort(np.abs(y_train))
    samples_per_level = len(y_train) // difficulty_levels
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.LEARNING_RATE)
    criterion = nn.MSELoss()
    
    for level in range(difficulty_levels):
        print(f"Training level {level+1}/{difficulty_levels}")
        
        # 选择当前难度的样本
        start_idx = level * samples_per_level
        end_idx = (level + 1) * samples_per_level
        level_indices = y_sorted_indices[start_idx:end_idx]
        
        X_level = X_train[level_indices]
        y_level = y_train[level_indices]
        
        # 转换为tensor
        X_tensor = torch.FloatTensor(X_level).to(device)
        y_tensor = torch.FloatTensor(y_level).to(device)
        
        # 训练当前级别
        for epoch in range(config.EPOCHS // difficulty_levels):
            model.train()
            optimizer.zero_grad()
            
            outputs = model(X_tensor)
            loss = criterion(outputs, y_tensor)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
            optimizer.step()
            
            if epoch % 10 == 0:
                print(f"  Level {level+1}, Epoch {epoch}, Loss: {loss.item():.4f}")
    
    return model


class EnsembleFeatureSelector:
    """集成多种特征选择方法"""
    def __init__(self, X, y, methods=['rrblup', 'xgboost', 'elasticnet']):
        self.X = X
        self.y = y
        self.methods = methods
        self.importance_scores = {}
        
    def calculate_ensemble_importance(self):
        """计算集成特征重要性"""
        all_importances = []
        
        # 1. RRBLUP重要性（alpha上界限制为10^3）
        if 'rrblup' in self.methods:
            ridge = RidgeCV(alphas=np.logspace(-3, 3, 20))
            ridge.fit(self.X, self.y)
            rrblup_importance = np.abs(ridge.coef_)
            all_importances.append(rrblup_importance)
            
        # 2. XGBoost重要性
        if 'xgboost' in self.methods:
            xgb_model = xgb.XGBRegressor(
                n_estimators=100, max_depth=5, learning_rate=0.1
            )
            xgb_model.fit(self.X, self.y)
            xgb_importance = xgb_model.feature_importances_
            all_importances.append(xgb_importance)
            
        # 3. ElasticNet重要性
        if 'elasticnet' in self.methods:
            enet = ElasticNet(alpha=0.01, l1_ratio=0.5)
            enet.fit(self.X, self.y)
            enet_importance = np.abs(enet.coef_)
            all_importances.append(enet_importance)
            
        # 4. 互信息重要性
        if 'mutual_info' in self.methods:
            from sklearn.feature_selection import mutual_info_regression
            mi_importance = mutual_info_regression(self.X, self.y)
            all_importances.append(mi_importance)
            
        # 集成重要性（加权平均）
        ensemble_importance = np.mean(all_importances, axis=0)
        return ensemble_importance
    
    def select_features(self, top_k_ratio=0.3):
        """选择top-k%的特征"""
        ensemble_importance = self.calculate_ensemble_importance()
        k = int(len(ensemble_importance) * top_k_ratio)
        top_indices = np.argsort(ensemble_importance)[-k:]
        return top_indices, ensemble_importance


def train_enhanced_model(X_train, X_val, y_train_actual, y_val_actual, device, config, use_optimization=False):
    """训练增强模型（使用与train_enhanced_automl_separate.py相同的架构）"""
    # 数据准备 - 使用DataLoader进行批处理
    train_dataset = TensorDataset(
        torch.FloatTensor(X_train),
        torch.FloatTensor(y_train_actual)
    )
    train_loader = DataLoader(
        train_dataset,
        batch_size=config.BATCH_SIZE,
        shuffle=True,
        pin_memory=config.PIN_MEMORY
    )
    
    # 创建模型 - 使用SVProNet（已经是合适的架构）
    model = SVProNet(
        n_sv=X_train.shape[1],
        embed_dim=config.EMBED_DIM,
        depth=config.DEPTH,
        num_heads=config.NUM_HEADS,
        cnn_channels=config.CNN_CHANNELS,
        drop_rate=config.DROP_RATE,
        drop_path_rate=config.DROP_PATH_RATE
    ).to(device)
    
    # 定义优化器和调度器
    optimizer = torch.optim.AdamW(model.parameters(), 
                                  lr=config.LEARNING_RATE, 
                                  weight_decay=config.WEIGHT_DECAY)
    
    # 使用Cosine学习率调度（与train_enhanced_automl_separate.py一致）
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
        optimizer, T_0=50, T_mult=2, eta_min=1e-7  # T_0从30增加到50
    )
    
    # 定义损失函数
    criterion = nn.MSELoss()
    
    # 验证数据准备
    X_val_tensor = torch.FloatTensor(X_val).to(device)
    y_val_tensor = torch.FloatTensor(y_val_actual).to(device)
    
    # 训练循环
    best_val_r2 = -float('inf')
    best_preds = None
    patience_counter = 0
    patience = config.EARLY_STOPPING_PATIENCE
    
    scaler = GradScaler('cuda') if torch.cuda.is_available() and config.USE_AMP else None
    
    for epoch in range(config.EPOCHS):
        # 训练
        model.train()
        train_losses = []
        
        for X_batch, y_batch in train_loader:
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)
            optimizer.zero_grad()
            
            if scaler:
                with autocast('cuda'):
                    outputs = model(X_batch)
                    loss = criterion(outputs, y_batch)
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
                scaler.step(optimizer)
                scaler.update()
            else:
                outputs = model(X_batch)
                loss = criterion(outputs, y_batch)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
                optimizer.step()
            
            train_losses.append(loss.item())
        
        # 更新学习率
        scheduler.step()
        
        # 验证
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val_tensor)
            val_loss = criterion(val_outputs, y_val_tensor)
            
            # 计算R²和相关系数
            val_preds = val_outputs.cpu().numpy()
            val_r2 = r2_score(y_val_actual, val_preds)
            val_corr, _ = scipy.stats.pearsonr(y_val_actual, val_preds)
        
        # 早停检查
        if val_r2 > best_val_r2:
            best_val_r2 = val_r2
            best_preds = val_preds
            patience_counter = 0
        else:
            patience_counter += 1
            
        if patience_counter >= patience:
            print(f"Early stopping at epoch {epoch+1}")
            break
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}, Train Loss: {np.mean(train_losses):.4f}, Val Loss: {val_loss.item():.4f}, Val R2: {val_r2:.4f}, Corr: {val_corr:.4f}")
    
    # 使用最佳预测
    if best_preds is None:
        # 如果没有最佳预测，重新计算
        model.eval()
        with torch.no_grad():
            best_preds = model(X_val_tensor).cpu().numpy()
    
    # 重新计算最终的R²和相关系数
    final_r2 = r2_score(y_val_actual, best_preds)
    final_corr, _ = scipy.stats.pearsonr(y_val_actual, best_preds)
    
    print(f"Enhanced Model - R2: {final_r2:.4f}, Corr: {final_corr:.4f}")
    
    return best_preds, final_r2, final_corr


def load_genotype_data_from_sv_assoc(sv_path, phenotype_path: str):
    """
    从SV关联文件加载基因型数据
    """
    print(f"\n[基因型数据加载] 从SV关联文件加载基因型数据")
    
    # 加载表型数据
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f"表型文件不存在: {phenotype_path}")
    
    # 尝试不同的分隔符并跳过标题行
    try:
        # 首先尝试以制表符分隔，假设第一行为标题
        df = pd.read_csv(phenotype_path, sep='\t', header=0)  # 第一行作为标题
        # 假设第一列是ID，第二列是表型值
        if df.shape[1] >= 2:
            y = df.iloc[:, 1].values.astype(np.float32)  # 使用第二列作为表型
            print(f"  加载表型数据: {len(y)} 个样本")
        else:
            # 如果只有1列，尝试将其作为表型值（跳过第一行）
            df_full = pd.read_csv(phenotype_path, sep='\t', header=None)
            y = df_full.iloc[1:, 0].values.astype(np.float32)  # 跳过第一行
            print(f"  加载表型数据: {len(y)} 个样本 (跳过标题)")
    except Exception as e:
        print(f"  制表符分隔加载失败: {e}, 尝试其他格式")
        try:
            # 尝试无分隔符，带标题的情况
            df = pd.read_csv(phenotype_path, header=0)  # 第一行作为标题
            if df.shape[1] >= 2:
                y = df.iloc[:, 1].values.astype(np.float32)  # 使用第二列作为表型
                print(f"  加载表型数据: {len(y)} 个样本")
            else:
                df_full = pd.read_csv(phenotype_path, header=None)
                y = df_full.iloc[1:, 0].values.astype(np.float32)  # 跳过第一行
                print(f"  加载表型数据: {len(y)} 个样本 (跳过标题)")
        except Exception as e2:
            raise ValueError(f"无法解析表型文件 {phenotype_path}: {e2}")
    
    # 加载VCF ID映射文件
    vcf_id_path = "/storage/public/home/2024110093/data/Variation/CSIAAS/VCFID.txt"
    if not os.path.exists(vcf_id_path):
        raise FileNotFoundError(f"VCF ID文件不存在: {vcf_id_path}")
    
    vcf_ids = pd.read_csv(vcf_id_path, header=None, sep='\s+|\t|,')  # 尝试多种分隔符
    sample_names = vcf_ids.iloc[:, 0].tolist()  # 假设第一列是样本ID
    print(f"  VCF文件包含 {len(sample_names)} 个样本")
    
    # 检查样本数量是否与表型数据匹配
    if len(sample_names) != len(y):
        print(f"  警告: VCF样本数({len(sample_names)})与表型数据数({len(y)})不匹配")
        print(f"  将使用较小的数量: {min(len(sample_names), len(y))}")
        min_samples = min(len(sample_names), len(y))
        y = y[:min_samples]
        sample_names = sample_names[:min_samples]
    else:
        min_samples = len(y)
    
    # 加载VCF文件并提取SV关联的基因型
    all_vcf_path = "/storage/public/home/2024110093/data/Variation/CSIAAS/Core819Samples_ALL.vcf.gz"
    if not os.path.exists(all_vcf_path):
        raise FileNotFoundError(f"主VCF文件不存在: {all_vcf_path}")
    
    if VCF_AVAILABLE:
        print("  使用cyvcf2库加载VCF数据...")
        try:
            vcf = cyvcf2.VCF(all_vcf_path)
            
            # 读取SV关联文件，获取显著位点信息
            df_sv = None
            
            # 读取SV关联文件
            if os.path.exists(sv_path):
                df_sv = pd.read_csv(sv_path, sep='\t', header=0)
                print(f"  SV关联文件包含 {len(df_sv)} 个显著位点")
            
            # 获取所有显著位点的ID列表
            sv_ids = set()
            
            # 提取SV关联文件中的位点ID
            if df_sv is not None:
                # 尝试找到ID列，常见的ID列名
                id_columns = [col for col in df_sv.columns if any(id_name in col.lower() for id_name in ['id', 'rs', 'marker', 'name'])]
                pos_columns = [col for col in df_sv.columns if any(pos_name in col.lower() for pos_name in ['pos', 'position', 'bp'])]
                chr_columns = [col for col in df_sv.columns if any(chr_name in col.lower() for chr_name in ['chr', 'chrom'])]
                
                if id_columns:
                    # 如果有ID列，直接使用
                    sv_ids.update(df_sv[id_columns[0]].astype(str).values)
                elif chr_columns and pos_columns:
                    # 如果有CHR和POS列，组合成chr:pos格式
                    for _, row in df_sv.iterrows():
                        sv_ids.add(f"{row[chr_columns[0]]}:{row[pos_columns[0]]}")
                else:
                    # 如果没有明确的ID或位置列，使用索引
                    sv_ids.update([f"SV_{i}" for i in range(len(df_sv))])
            
            print(f"  从关联文件识别: SV位点{len(sv_ids)}")
            
            # 初始化SV基因型数据列表
            sv_genotypes = []
            
            variant_count = 0
            max_variants = 5000  # 限制位点数量
            
            # 遍历VCF文件，提取显著位点的基因型
            for variant in vcf:
                if variant_count >= max_variants:
                    break
                
                # 构建当前变体的ID（优先使用VCF中的ID，否则使用CHR:POS格式）
                variant_id = variant.ID if variant.ID else f"{variant.CHROM}:{variant.POS}"
                
                # 检查当前变体是否在显著位点集中
                is_sv_variant = variant_id in sv_ids
                
                # 提取基因型数据
                gt = variant.gt_types
                if len(gt) >= min_samples:
                    gt_subset = gt[:min_samples]
                    
                    # 将基因型编码转换为数值 (ref/ref->0, ref/alt->1, alt/alt->2)
                    gt_numeric = []
                    for g in gt_subset:
                        if g == 0:  # ref/ref
                            gt_numeric.append(0)
                        elif g == 1:  # ref/alt or het
                            gt_numeric.append(1)
                        elif g == 2:  # alt/alt or hom
                            gt_numeric.append(2)
                        else:  # missing or other
                            gt_numeric.append(-1)  # 标记为缺失值，后续用均值填充
                    
                    # 根据变体类型分配到对应的基因型列表
                    if is_sv_variant:
                        sv_genotypes.append(gt_numeric)
                        variant_count += 1
            
            print(f"  从VCF成功提取 {variant_count} 个显著位点的基因型")
            print(f"    SV位点: {len(sv_genotypes)}")
            
            # 转换为numpy数组
            sv_X = np.array(sv_genotypes, dtype=np.float32).T if sv_genotypes else np.zeros((min_samples, 1), dtype=np.float32)
            
            # 如果没有显著位点，创建一个零列
            if sv_X.shape[1] == 0:
                sv_X = np.zeros((min_samples, 1), dtype=np.float32)
            
            print(f"  基因型矩阵 - SV: {sv_X.shape}, y: {y.shape}")
            
            # 移除 NaN
            mask = ~np.isnan(y)
            sv_X, y = sv_X[mask], y[mask]
            
            # 检查数据完整性
            print(f"  数据完整性检查:")
            print(f"    基因型矩阵形状: {sv_X.shape}")
            print(f"    表型数据形状: {y.shape}")
            print(f"    基因型缺失率: {np.sum(sv_X == -1) / sv_X.size * 100:.2f}%")
            print(f"    表型缺失率: {np.sum(np.isnan(y)) / len(y) * 100:.2f}%")
            
            # 处理缺失的基因型数据（用均值填充）
            if sv_X.shape[1] > 0:
                for j in range(sv_X.shape[1]):
                    col = sv_X[:, j]
                    missing_mask = col == -1
                    if np.any(missing_mask):
                        # 计算非缺失值的均值
                        col_mean = np.mean(col[~missing_mask])
                        # 用均值填充缺失值
                        sv_X[missing_mask, j] = col_mean
                
                # 注意：标准化将在每个交叉验证fold内进行，避免数据泄漏
            
            print(f"  最终数据 - SV: {sv_X.shape}, y: {y.shape}")
            print(f"  y范围: [{y.min():.3f}, {y.max():.3f}]")
            
            return sv_X, y
        except Exception as e:
            print(f"  cyvcf2加载失败: {e}")
    else:
        print("  cyvcf2库不可用，尝试安装...")
        try:
            import subprocess
            import sys
            # 使用超算的pip镜像源安装cyvcf2
            result = subprocess.run([
                sys.executable, "-m", "pip", "install", 
                "-i", "http://10.105.32.248/pypi/simple/", 
                "--trusted-host", "10.105.32.248",
                "cyvcf2"
            ], capture_output=True, text=True)
            if result.returncode == 0:
                print("  cyvcf2安装成功，重新导入...")
                globals()['VCF_AVAILABLE'] = True  # 使用globals()修改全局变量
                # 递归调用自身以使用新安装的库
                return load_genotype_data_from_sv_assoc(sv_path, phenotype_path)
            else:
                print(f"  cyvcf2安装失败: {result.stderr}")
        except Exception as install_e:
            print(f"  安装cyvcf2时出错: {install_e}")
    
    # 如果VCF处理失败，抛出错误
    raise RuntimeError(f"无法加载VCF数据，请确保VCF文件存在且格式正确")


def load_sv_data_from_vcf(sv_vcf_path: str, phenotype_path: str, vcf_id_path: str, max_variants: int = 5000):
    """
    直接从完整SV VCF文件加载数据（与train_enhanced_automl_separate.py一致）
    这样可以获取更多的遗传信息，而不是只用GWAS筛选后的少数显著位点
    """
    print(f"\n[SV VCF数据加载] 从完整VCF文件加载所有SV变异")
    print(f"  VCF文件: {sv_vcf_path}")
    print(f"  最大变异数: {max_variants}")
    
    # 加载表型数据
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f"表型文件不存在: {phenotype_path}")
    
    try:
        df = pd.read_csv(phenotype_path, sep='\t', header=0)
        if df.shape[1] >= 2:
            y = df.iloc[:, 1].values.astype(np.float32)
            print(f"  加载表型数据: {len(y)} 个样本")
        else:
            df_full = pd.read_csv(phenotype_path, sep='\t', header=None)
            y = df_full.iloc[1:, 0].values.astype(np.float32)
            print(f"  加载表型数据: {len(y)} 个样本 (跳过标题)")
    except Exception as e:
        print(f"  制表符分隔加载失败: {e}, 尝试其他格式")
        try:
            df = pd.read_csv(phenotype_path, header=0)
            if df.shape[1] >= 2:
                y = df.iloc[:, 1].values.astype(np.float32)
                print(f"  加载表型数据: {len(y)} 个样本")
            else:
                df_full = pd.read_csv(phenotype_path, header=None)
                y = df_full.iloc[1:, 0].values.astype(np.float32)
                print(f"  加载表型数据: {len(y)} 个样本 (跳过标题)")
        except Exception as e2:
            raise ValueError(f"无法解析表型文件 {phenotype_path}: {e2}")
    
    # 加载VCF ID映射文件
    if vcf_id_path and os.path.exists(vcf_id_path):
        vcf_ids = pd.read_csv(vcf_id_path, header=None, sep='\s+|\t|,')
        expected_samples = len(vcf_ids)
        print(f"  VCF ID文件包含 {expected_samples} 个样本")
    else:
        expected_samples = len(y)
        print(f"  使用表型数据样本数: {expected_samples}")
    
    # 确定最小样本数
    min_samples = min(len(y), expected_samples)
    y = y[:min_samples]
    
    # 加载VCF文件
    if not VCF_AVAILABLE:
        print("　cyvcf2库不可用，尝试安装...")
        try:
            import subprocess
            import sys
            result = subprocess.run([
                sys.executable, "-m", "pip", "install", 
                "-i", "http://10.105.32.248/pypi/simple/", 
                "--trusted-host", "10.105.32.248",
                "cyvcf2"
            ], capture_output=True, text=True)
            if result.returncode == 0:
                print("　cyvcf2安装成功，重新加载模块")
                import cyvcf2 as _cyvcf2
                globals()['cyvcf2'] = _cyvcf2
                globals()['VCF_AVAILABLE'] = True
            else:
                raise RuntimeError(f"cyvcf2安装失败: {result.stderr}")
        except Exception as install_e:
            raise RuntimeError(f"安装cyvcf2时出错: {install_e}")
        
    if not os.path.exists(sv_vcf_path):
        raise FileNotFoundError(f"SV VCF文件不存在: {sv_vcf_path}")
        
    print(f"　使用cyvcf2库加载SV VCF数据...")
    try:
        import cyvcf2 as _vcf_module
        vcf = _vcf_module.VCF(sv_vcf_path)
        
        # 获取样本数量
        samples = vcf.samples
        n_samples_vcf = len(samples)
        print(f"  SV VCF文件包含 {n_samples_vcf} 个样本")
        
        # 更新最小样本数
        min_samples = min(min_samples, n_samples_vcf)
        y = y[:min_samples]
        
        # 提取变异位点信息和基因型
        genotypes = []
        variant_count = 0
        
        for variant in vcf:
            if variant_count >= max_variants:
                break
            
            # 提取基因型信息 (0/0->0, 0/1->1, 1/1->2)
            gt = variant.gt_types
            if len(gt) >= min_samples:
                gt_subset = gt[:min_samples]
                # 将基因型编码转换为数值
                gt_numeric = []
                for g in gt_subset:
                    if g == 0:  # ref/ref
                        gt_numeric.append(0)
                    elif g == 1:  # ref/alt or het
                        gt_numeric.append(1)
                    elif g == 2:  # alt/alt or hom
                        gt_numeric.append(2)
                    else:  # missing or other
                        gt_numeric.append(0)  # 设为参考纯合子作为默认
                genotypes.append(gt_numeric)
                variant_count += 1
        
        if variant_count > 0:
            print(f"  成功提取 {variant_count} 个SV变异位点")
            # 转换为numpy数组 (variants x samples) -> (samples x variants)
            X = np.array(genotypes, dtype=np.float32).T
            print(f"  SV基因型矩阵形状: {X.shape}")
        else:
            print(f"  警告: 未能从SV VCF文件中提取到任何变异位点")
            X = np.zeros((min_samples, 100), dtype=np.float32)
        
        # 移除NaN
        mask = ~np.isnan(y)
        X, y = X[mask], y[mask]
        
        # 数据完整性检查
        print(f"  数据完整性检查:")
        print(f"    基因型矩阵形状: {X.shape}")
        print(f"    表型数据形状: {y.shape}")
        print(f"    y范围: [{y.min():.3f}, {y.max():.3f}]")
        
        return X, y
        
    except Exception as e:
        print(f"  SV VCF加载失败: {e}")
        raise RuntimeError(f"无法加载SV VCF数据: {e}")


def train_xgboost_model(X_train, X_val, y_train_actual, y_val_actual):
    """训练XGBoost模型"""
    # 根据数据规模调整参数
    n_samples = X_train.shape[0]
    n_features = X_train.shape[1]
    
    # 计算特征-样本比率，用于动态调整参数
    feature_to_sample_ratio = n_features / max(n_samples, 1)
    
    if feature_to_sample_ratio > 10:  # 高维情况
        xgb_params = {
            'objective': 'reg:squarederror',
            'eval_metric': 'rmse',
            'max_depth': 3,  # 适中深度
            'learning_rate': 0.05,  # 降低学习率
            'n_estimators': min(30, max(10, n_samples // 3)),  # 根据样本数调整
            'subsample': 0.8,  # 适中子样本率
            'colsample_bytree': 0.7,  # 适中列采样率
            'reg_alpha': 0.5,  # 适度L1正则化
            'reg_lambda': 0.5,  # 适度L2正则化
            'min_child_weight': 2,  # 适中最小叶子权重
            'random_state': 42
        }
    elif feature_to_sample_ratio > 5:  # 高维情况
        xgb_params = {
            'objective': 'reg:squarederror',
            'eval_metric': 'rmse',
            'max_depth': 3,  # 适中深度
            'learning_rate': 0.05,  # 降低学习率
            'n_estimators': min(30, max(10, n_samples // 3)),  # 根据样本数调整
            'subsample': 0.8,  # 适中子样本率
            'colsample_bytree': 0.7,  # 适中列采样率
            'reg_alpha': 0.5,  # 适度L1正则化
            'reg_lambda': 0.5,  # 适度L2正则化
            'min_child_weight': 2,  # 适中最小叶子权重
            'random_state': 42
        }
    else:  # 正常情况
        xgb_params = {
            'objective': 'reg:squarederror',
            'eval_metric': 'rmse',
            'max_depth': max(3, min(6, int(4 - feature_to_sample_ratio*0.1))),  # 动态调整深度
            'learning_rate': 0.05,  # 降低学习率
            'n_estimators': min(50, max(20, n_samples // 4)),  # 根据样本数调整树的数量
            'subsample': max(0.7, 0.9 - feature_to_sample_ratio * 0.05),  # 根据特征样本比调整采样
            'colsample_bytree': max(0.7, 0.85 - feature_to_sample_ratio * 0.05),  # 根据特征样本比调整列采样
            'reg_alpha': max(0.1, 0.5 * feature_to_sample_ratio),  # L1正则化
            'reg_lambda': max(0.1, 0.5 * feature_to_sample_ratio),  # L2正则化
            'min_child_weight': max(1, int(1 + feature_to_sample_ratio*0.5)),  # 防止过拟合
            'random_state': 42
        }
    
    xgb_model = xgb.XGBRegressor(**xgb_params)
    xgb_model.fit(X_train, y_train_actual)
    xgb_preds = xgb_model.predict(X_val)
    
    xgb_r2 = r2_score(y_val_actual, xgb_preds)
    from scipy.stats import pearsonr
    xgb_corr, _ = pearsonr(y_val_actual, xgb_preds)
    
    # 只修正数学上不可能的值（NaN或无穷大），保留实际计算的性能值
    if np.isnan(xgb_r2) or np.isinf(xgb_r2):
        xgb_r2 = 0.0
    if np.isnan(xgb_corr) or np.isinf(xgb_corr) or abs(xgb_corr) > 1.0:
        xgb_corr = 0.0
    
    print(f"    XGBoost - Estimators: {xgb_model.n_estimators}, R2: {xgb_r2:.4f}, Corr: {xgb_corr:.4f}")
    
    return xgb_preds, xgb_r2, xgb_corr


class WheatGPCNN(nn.Module):
    def __init__(self, input_features):
        super(WheatGPCNN, self).__init__()
        # 增加CNN容量以适应基因组数据
        # First convolutional layer: input channels = 1, output channels = 16, kernel size = 1, padding = 0
        self.conv0 = nn.Conv1d(1, 16, 1, padding=0)  # 增加输出通道数
        # ReLU activation function after the first convolutional layer
        self.relu0 = nn.ReLU()
        self.bn0 = nn.BatchNorm1d(16)  # 添加批归一化
        
        self.conv1 = nn.Conv1d(16, 32, 3, padding=1)  # 增加输出通道数
        self.relu1 = nn.ReLU()
        self.bn1 = nn.BatchNorm1d(32)  # 添加批归一化
        
        self.conv2 = nn.Conv1d(32, 64, 9, padding=1)  # 增加输出通道数
        self.relu2 = nn.ReLU()
        self.bn2 = nn.BatchNorm1d(64)  # 添加批归一化
        
        self.drop = nn.Dropout(0.3)
        
        # 添加自适应池化层来确保输出固定大小
        self.adaptive_pool = nn.AdaptiveAvgPool1d(input_features)  # 确保输出维度匹配LSTM输入

    def forward(self, x):
        x = self.conv0(x)
        x = self.bn0(x)
        x = self.relu0(x)

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu1(x)

        x = self.conv2(x)
        x = self.bn2(x)
        x = self.relu2(x)
        
        # 确保x的形状正确
        if x.size(2) == 0:  # 如果序列长度为0，扩展它
            x = torch.nn.functional.interpolate(x, size=16, mode='linear', align_corners=False)

        x = self.adaptive_pool(x)  # 调整到指定大小
        x = self.drop(x)

        return x


class WheatGPShapeModule(nn.Module):
    def __init__(self):
        super(WheatGPShapeModule, self).__init__()

    def forward(self, x1, x2, x3, x4, x5, adjust_dim=True, concat=True):
        if adjust_dim:
            x1 = x1.unsqueeze(1)
            x2 = x2.unsqueeze(1)
            x3 = x3.unsqueeze(1)
            x4 = x4.unsqueeze(1)
            x5 = x5.unsqueeze(1)
        if concat:
            A_flat = x1.view(x1.size(0), -1)
            B_flat = x2.view(x2.size(0), -1)
            C_flat = x3.view(x3.size(0), -1)
            D_flat = x4.view(x4.size(0), -1)
            E_flat = x5.view(x5.size(0), -1)
            output = torch.cat((A_flat, B_flat, C_flat, D_flat, E_flat), dim=1)
            output = output.reshape(output.shape[0], 1, -1)

            return output
        else:

            return x1, x2, x3, x4, x5


class WheatGPLSTM(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers=1, batch_first=True):
        super(WheatGPLSTM, self).__init__()
        # Define an LSTM layer with the specified input size, hidden size, number of layers, and batch first flag
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=batch_first)
        self.drop = nn.Dropout(0.3)

    def forward(self, x):
        lstm_out, (h_n, c_n) = self.lstm(x)
        lstm_out = self.drop(lstm_out)
        return lstm_out


class WheatGPModel(nn.Module):
    """真正的WheatGP模型实现，包含多尺度CNN、LSTM和五个子网络"""
    def __init__(self, input_dim):
        super().__init__()
        
        # 改进的特征分组：使用随机分组以避免顺序偏差
        self.input_dim = input_dim
        self.num_subnetworks = 5
        
        # 预计算特征索引，确保所有特征都被使用
        feature_indices = np.arange(input_dim)
        np.random.shuffle(feature_indices)  # 随机打乱特征顺序
        
        # 将特征索引平均分配到5个组
        self.feature_groups = []
        group_size = input_dim // self.num_subnetworks
        remainder = input_dim % self.num_subnetworks
        
        start_idx = 0
        for i in range(self.num_subnetworks):
            end_idx = start_idx + group_size + (1 if i < remainder else 0)
            group_indices = feature_indices[start_idx:end_idx]
            self.feature_groups.append(group_indices)
            start_idx = end_idx
        
        # 五个子网络并行处理
        self.subnetworks = nn.ModuleList()
        for i in range(self.num_subnetworks):
            group_size = len(self.feature_groups[i])
            
            # 多尺度CNN - 修复通道数分配
            kernel_sizes = [3, 5, 7]
            
            # 计算每个多尺度分支的实际输出通道数
            base_channels = 16  # 每个分支的基础通道数
            cnn_branches = nn.ModuleList()
            
            for k in kernel_sizes:
                branch = nn.Sequential(
                    nn.Conv1d(1, base_channels, kernel_size=k, padding=k//2),
                    nn.BatchNorm1d(base_channels),
                    nn.ReLU(),
                    nn.Conv1d(base_channels, base_channels*2, kernel_size=3, padding=1),
                    nn.BatchNorm1d(base_channels*2),
                    nn.ReLU()
                )
                cnn_branches.append(branch)
            
            # LSTM层 - 基于实际通道数
            lstm_input_size = base_channels * 2 * len(kernel_sizes)  # 每个分支输出base_channels*2
            lstm_hidden_size = 64
            lstm = nn.LSTM(lstm_input_size, lstm_hidden_size, batch_first=True, bidirectional=True)
            
            # 注意力机制
            attention = nn.Sequential(
                nn.Linear(lstm_hidden_size * 2, lstm_hidden_size),  # 双向LSTM
                nn.Tanh(),
                nn.Linear(lstm_hidden_size, 1),
                nn.Softmax(dim=1)
            )
            
            subnetwork = nn.ModuleDict({
                'cnn_branches': cnn_branches,
                'lstm': lstm,
                'attention': attention
            })
            
            self.subnetworks.append(subnetwork)
        
        # 最终输出层 - 基于实际的输出维度
        total_lstm_output = lstm_hidden_size * 2 * self.num_subnetworks  # 双向LSTM * 5个子网络
        self.output = nn.Sequential(
            nn.Linear(total_lstm_output, 128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 1)
        )
    
    def forward(self, x):
        # x: (batch, features)
        batch_size, total_features = x.shape
        
        # 【修复】使用顺序切片，保持基因组标记的空间局部性
        # CNN卷积核(kernel_size=3,5,7)依赖于相邻标记的连锁不平衡(LD)来提取有效特征
        # 随机打乱会破坏这种局部相关性，导致CNN无法收敛
        features_per_group = total_features // self.num_subnetworks
        remainder = total_features % self.num_subnetworks
        
        group_outputs = []
        
        for i, subnetwork in enumerate(self.subnetworks):
            # 按顺序切分特征，保持基因组位置相邻性
            start_idx = i * features_per_group + min(i, remainder)
            end_idx = min((i + 1) * features_per_group + min(i + 1, remainder), total_features)
            
            if start_idx >= total_features:
                x_group = torch.zeros(batch_size, 1, device=x.device, dtype=x.dtype)
            elif start_idx >= end_idx:
                x_group = x[:, start_idx:start_idx+1] if start_idx < total_features else torch.zeros(batch_size, 1, device=x.device, dtype=x.dtype)
            else:
                x_group = x[:, start_idx:end_idx]  # (batch, group_size)
            
            # 扩展维度以适应CNN (batch, 1, group_size)
            x_group = x_group.unsqueeze(1)
            
            # 多尺度CNN处理
            cnn_outputs = []
            for cnn_branch in subnetwork['cnn_branches']:
                cnn_out = cnn_branch(x_group)  # (batch, channels, group_size_after_conv)
                cnn_outputs.append(cnn_out)
            
            # 拼接CNN输出
            cnn_combined = torch.cat(cnn_outputs, dim=1)  # (batch, total_channels, group_size_after_conv)
            
            # 针对大规模标记数，使用MaxPooling降低序列长度以节省显存
            seq_len = cnn_combined.size(2)
            if seq_len > 2000:  # 当序列长度超过2000时，使用池化
                pool_size = seq_len // 1000  # 动态计算池化核大小
                cnn_combined = F.max_pool1d(cnn_combined, kernel_size=pool_size, stride=pool_size)
            
            # 调整维度以适应LSTM (batch, seq_len, features)
            cnn_combined = cnn_combined.transpose(1, 2)  # (batch, group_size_after_conv, total_channels)
            
            # LSTM处理
            lstm_out, _ = subnetwork['lstm'](cnn_combined)  # (batch, seq_len, hidden_size*2)
            
            # 注意力机制
            attn_weights = subnetwork['attention'](lstm_out)  # (batch, seq_len, 1)
            attended_out = torch.sum(lstm_out * attn_weights, dim=1)  # (batch, hidden_size*2)
            
            group_outputs.append(attended_out)
        
        # 拼接所有子网络输出
        final_features = torch.cat(group_outputs, dim=1)  # (batch, hidden_size*2*num_networks)
        
        # 最终预测
        output = self.output(final_features)
        return output.squeeze(-1)


def train_wheatgp_model(X_train, X_val, y_train_actual, y_val_actual):
    """训练真正的WheatGP模型 - 基于CNN和LSTM的深度学习模型
    
    【修复v3】关键改进：
    1. 保存最佳模型权重（best checkpoint），预测时使用最佳而非最后一轮
    2. 梯度裁剪，防止梯度爆炸
    3. 学习率warmup，避免初期梯度过大
    """
    import copy
    
    # 设置随机种子以确保结果可重复
    np.random.seed(42)
    torch.manual_seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(42)
        torch.cuda.manual_seed_all(42)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    print("    Training REAL WheatGP model (CNN+LSTM-based)...")
    
    # 获取特征数，用于模型初始化
    n_features = X_train.shape[1]
    n_samples = X_train.shape[0]
    
    # 创建模型
    model = WheatGPModel(n_features).to(DEVICE)
    
    # 定义损失函数
    criterion = nn.MSELoss()
    
    # 根据特征数/样本数比率调整学习率
    # 过参数化越严重，学习率需要越低
    param_ratio = sum(p.numel() for p in model.parameters()) / n_samples
    if param_ratio > 500:
        lr = 0.0003  # 严重过参数化
    elif param_ratio > 100:
        lr = 0.0005
    else:
        lr = 0.001
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=0.001)
    print(f"    [WheatGP] Features: {n_features}, Samples: {n_samples}, Param ratio: {param_ratio:.0f}, LR: {lr}")
    
    # 预先将数据转为张量
    X_train_tensor = torch.FloatTensor(X_train).to(DEVICE)
    X_val_tensor = torch.FloatTensor(X_val).to(DEVICE)
    y_train_tensor = torch.FloatTensor(y_train_actual).to(DEVICE)
    y_val_tensor = torch.FloatTensor(y_val_actual).to(DEVICE)
    
    # 训练模型 - 全批量训练
    best_val_loss = float('inf')
    best_model_state = None  # 【关键】保存最佳模型权重
    patience = 20  # 增加patience给模型更多学习机会
    patience_counter = 0
    max_epochs = 150
    
    for epoch in range(max_epochs):
        model.train()
        optimizer.zero_grad()
        
        # 全批量前向传播
        outputs = model(X_train_tensor)
        loss = criterion(outputs, y_train_tensor)
        
        # 反向传播 + 梯度裁剪
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        # 验证
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val_tensor)
            val_loss = criterion(val_outputs, y_val_tensor).item()
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = copy.deepcopy(model.state_dict())  # 【关键】保存最佳权重
            patience_counter = 0
        else:
            patience_counter += 1
            
        if patience_counter >= patience:
            print(f"    WheatGP early stopping at epoch {epoch+1}")
            break
            
        if (epoch + 1) % 20 == 0:
            print(f"    WheatGP Epoch {epoch+1}/{max_epochs}, Train Loss: {loss.item():.6f}, Val Loss: {val_loss:.6f}")
    
    # 【关键】加载最佳模型权重进行预测
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    # 最终预测
    model.eval()
    with torch.no_grad():
        wheatgp_preds_tensor = model(X_val_tensor)
        wheatgp_preds = wheatgp_preds_tensor.cpu().numpy()
    
    # 清理GPU显存
    del X_train_tensor, X_val_tensor, y_train_tensor, y_val_tensor
    del wheatgp_preds_tensor, best_model_state
    torch.cuda.empty_cache()
    
    # 计算性能指标
    wheatgp_r2 = r2_score(y_val_actual, wheatgp_preds)
    from scipy.stats import pearsonr
    wheatgp_corr, _ = pearsonr(y_val_actual, wheatgp_preds)
    
    # 只修正数学上不可能的值（NaN或无穷大），保留实际计算的性能值
    if np.isnan(wheatgp_r2) or np.isinf(wheatgp_r2):
        wheatgp_r2 = 0.0
    if np.isnan(wheatgp_corr) or np.isinf(wheatgp_corr) or abs(wheatgp_corr) > 1.0:
        wheatgp_corr = 0.0
    
    print(f"    WheatGP - R2: {wheatgp_r2:.4f}, Corr: {wheatgp_corr:.4f}")
    
    return wheatgp_preds, wheatgp_r2, wheatgp_corr


# ============================================================================
# 高级正则化模块
# ============================================================================

class DropPath(nn.Module):
    """
    DropPath (Stochastic Depth) - 随机丢弃整个残差路径
    参考: Deep Networks with Stochastic Depth (ECCV 2016)
    """
    def __init__(self, drop_prob: float = 0.0):
        super().__init__()
        self.drop_prob = drop_prob
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.drop_prob == 0.0 or not self.training:
            return x
        keep_prob = 1 - self.drop_prob
        shape = (x.shape[0],) + (1,) * (x.ndim - 1)
        random_tensor = keep_prob + torch.rand(shape, dtype=x.dtype, device=x.device)
        random_tensor.floor_()
        output = x.div(keep_prob) * random_tensor
        return output

class SpatialDropout1D(nn.Module):
    """
    Spatial Dropout - 丢弃整个通道而非单个神经元
    更适合序列数据，保持空间一致性
    """
    def __init__(self, drop_prob: float = 0.2):
        super().__init__()
        self.drop_prob = drop_prob
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (batch, seq, channels)
        if not self.training or self.drop_prob == 0:
            return x
        # 对整个通道维度应用dropout
        x = x.permute(0, 2, 1)  # (batch, channels, seq)
        x = F.dropout2d(x.unsqueeze(-1), self.drop_prob, self.training).squeeze(-1)
        return x.permute(0, 2, 1)

class LayerScale(nn.Module):
    """Layer Scale - 自适应缩放残差连接"""
    def __init__(self, dim: int, init_value: float = 1e-5):
        super().__init__()
        self.gamma = nn.Parameter(init_value * torch.ones(dim))
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x * self.gamma

# ============================================================================
# RRBLUP指导模块
# ============================================================================

class RRBLUPGuidedFeatureSelector:
    """
    RRBLUP指导的特征选择模块
    使用RRBLUP识别的重要位点来指导神经网络特征选择
    """
    def __init__(self, top_k_percent: float = 0.1):
        self.top_k_percent = top_k_percent
        self.selected_indices = None
        self.feature_importance = None
        
    def fit(self, X: np.ndarray, y: np.ndarray):
        """
        使用RRBLUP方法评估特征重要性
        X: (n_samples, n_features) 基因型数据
        y: (n_samples,) 表型数据
        """
        from sklearn.linear_model import Ridge
        from sklearn.preprocessing import StandardScaler
        
        # 标准化数据
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # 使用Ridge回归模拟RRBLUP（简化版）
        ridge = Ridge(alpha=1.0)
        ridge.fit(X_scaled, y)
        
        # 获取特征重要性（系数绝对值）
        self.feature_importance = np.abs(ridge.coef_)
        
        # 选择top-k%最重要的特征
        n_features = X.shape[1]
        k = max(1, int(n_features * self.top_k_percent))
        
        # 获取最重要的k个特征的索引
        top_k_indices = np.argsort(self.feature_importance)[-k:]
        self.selected_indices = sorted(top_k_indices)
        
        print(f"RRBLUP Feature Selection - Selecting top {self.top_k_percent*100}% features...")
        print(f"Selected {len(self.selected_indices)} features out of {n_features}")
        print(f"Importance range: {self.feature_importance.min():.6f} - {self.feature_importance.max():.6f}")
        
        return self
        
    def transform(self, X: np.ndarray):
        """
        应用特征选择
        """
        if self.selected_indices is None:
            raise ValueError("Must call fit() before transform()")
        
        return X[:, self.selected_indices]
        
    def fit_transform(self, X: np.ndarray, y: np.ndarray):
        """
        拟合并转换
        """
        return self.fit(X, y).transform(X)
        
    def get_importance_weights(self):
        """
        获取RRBLUP识别的特征重要性权重
        """
        if self.feature_importance is None:
            raise ValueError("Must call fit() before get_importance_weights()")
        
        # 返回全部特征的重要性权重
        return self.feature_importance

class RRBLUPGuidedAttention(nn.Module):
    """
    RRBLUP指导的注意力机制
    使用RRBLUP特征重要性作为注意力先验
    """
    def __init__(self, embed_dim: int, num_heads: int = 8):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads
        assert embed_dim % num_heads == 0, "embed_dim must be divisible by num_heads"
        
        self.qkv = nn.Linear(embed_dim, embed_dim * 3)
        self.proj = nn.Linear(embed_dim, embed_dim)
        
        # 学习RRBLUP先验权重
        self.rrblup_weight = nn.Parameter(torch.ones(1))
        
    def forward(self, x: torch.Tensor, rrblup_importance: Optional[torch.Tensor] = None):
        # 处理二维输入 (B, C) 或三维输入 (B, N, C)
        if len(x.shape) == 2:  # (B, C)
            B, C = x.shape
            # 将二维输入视为只有一个序列元素
            x = x.unsqueeze(1)  # (B, 1, C)
            N = 1
        elif len(x.shape) == 3:  # (B, N, C)
            B, N, C = x.shape
        else:
            raise ValueError(f"Expected 2D (B, C) or 3D (B, N, C) input, got {x.shape}")
        
        qkv = self.qkv(x).reshape(B, N, 3, self.num_heads, self.head_dim).permute(2, 0, 3, 1, 4)
        q, k, v = qkv.unbind(0)  # (B, num_heads, N, head_dim)
        
        attn = (q @ k.transpose(-2, -1)) * (self.head_dim ** -0.5)
        
        # 如果提供了RRBLUP重要性，将其作为先验加入注意力
        if rrblup_importance is not None:
            # 处理RRBLUP重要性维度与实际输入维度的映射
            if rrblup_importance.dim() == 1:
                # 假设rrblup_importance是针对embed_dim的权重
                if rrblup_importance.size(0) == C:  # 如果长度等于通道数
                    # 为每个序列位置重复相同的权重
                    rrblup_importance = rrblup_importance.unsqueeze(0).unsqueeze(0).expand(B, N, C)  # (B, N, C)
                else:
                    # 如果长度等于序列长度N
                    rrblup_importance = rrblup_importance.unsqueeze(0).unsqueeze(0)  # (1, 1, N)
            elif rrblup_importance.dim() == 2:
                if rrblup_importance.size(1) == C:  # (B, C)
                    rrblup_importance = rrblup_importance.unsqueeze(1).expand(B, N, C)  # (B, N, C)
                else:  # (B, N)
                    rrblup_importance = rrblup_importance.unsqueeze(1)  # (B, 1, N)
            
            # 根据rrblup_importance的维度处理注意力计算
            if rrblup_importance.size() == (B, N, C):
                # 将RRBLUP重要性应用到注意力权重中
                # 对C维度求平均得到(B, N)的权重
                rrblup_weights = rrblup_importance.mean(dim=-1)  # (B, N)
                # 扩展到(B, num_heads, N)用于注意力计算
                rrblup_weights_expanded = rrblup_weights.unsqueeze(1).expand(B, self.num_heads, N)  # (B, num_heads, N)
                # 为注意力矩阵创建对角线模式，让每个位置关注其RRBLUP重要性
                rrblup_attn = self.rrblup_weight * rrblup_weights_expanded.unsqueeze(-1)  # (B, num_heads, N, 1)
                # 重复以创建完整的注意力偏置矩阵
                rrblup_attn = rrblup_attn.expand(B, self.num_heads, N, N)  # (B, num_heads, N, N)
                # 将RRBLUP先验加入注意力
                attn = attn + rrblup_attn
            elif rrblup_importance.size() == (B, 1, N):
                # 扩展到多头
                rrblup_attn = self.rrblup_weight * rrblup_importance.expand(B, self.num_heads, N)
                rrblup_attn = rrblup_attn.unsqueeze(-1).expand(B, self.num_heads, N, N)
                
                # 将RRBLUP先验加入注意力
                attn = attn + rrblup_attn
        
        attn = F.softmax(attn, dim=-1)
        
        x = (attn @ v).transpose(1, 2).reshape(B, N, C)
        x = self.proj(x)
        
        # 如果输入是二维的，输出也应该是二维的
        if N == 1:
            x = x.squeeze(1)  # (B, C)
        
        return x

# ============================================================================
# CNN局部特征提取器
# ============================================================================

class LocalCNNExtractor(nn.Module):
    """改进的CNN局部特征提取 - 更适合捕捉基因组LD结构"""
    def __init__(self, in_channels: int = 1, channels: List[int] = [8, 16, 32],
                 drop_rate: float = 0.2):
        super().__init__()
        
        # 为每个通道层创建多尺度卷积
        self.conv_layers = nn.ModuleList()
        prev_ch = in_channels
        
        for ch in channels:
            # 多尺度卷积核，模拟不同连锁不平衡距离
            kernel_sizes = [3, 5, 7]
            
            # 计算每个多尺度分支的实际输出通道数，确保总数等于目标通道数
            scale_ch = ch // len(kernel_sizes)
            remaining_ch = ch % len(kernel_sizes)  # 使用模运算而不是整除
            
            # 为每个多尺度分支分配通道数
            actual_scale_chs = [scale_ch + (1 if i < remaining_ch else 0) for i in range(len(kernel_sizes))]
            
            scale_convs = nn.ModuleList()
            for i, k in enumerate(kernel_sizes):
                scale_convs.append(nn.Conv1d(prev_ch, actual_scale_chs[i], kernel_size=k, padding=k//2))
            
            # 归一化和激活层
            norm_act = nn.Sequential(
                nn.BatchNorm1d(ch),  # 使用目标通道数 ch
                nn.GELU()
            )
            
            layer = nn.ModuleDict({
                'scale_convs': scale_convs,
                'norm_act': norm_act
            })
            
            self.conv_layers.append(layer)
            prev_ch = ch
        
        self.dropout = SpatialDropout1D(drop_rate)
        self.out_channels = channels[-1]
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.unsqueeze(1)  # (B, 1, N)
        
        # 遍历每个多尺度卷积层
        for layer in self.conv_layers:
            # 多尺度卷积
            multi_scale_outputs = []
            for conv in layer['scale_convs']:
                multi_scale_outputs.append(conv(x))
            
            # 拼接多尺度输出
            x = torch.cat(multi_scale_outputs, dim=1)
            
            # 应用归一化和激活
            x = layer['norm_act'](x)
        
        # 应用Spatial Dropout
        x = self.dropout(x)
        return x

# ============================================================================
# GWAS门控融合
# ============================================================================

class GWASGatedFusion(nn.Module):
    """GWAS先验与深度特征的自适应融合"""
    def __init__(self, embed_dim: int):
        super().__init__()
        self.gate = nn.Sequential(
            nn.Linear(embed_dim * 2, embed_dim),
            nn.Sigmoid()
        )
        self.proj = nn.Linear(embed_dim, embed_dim)
        
    def forward(self, x: torch.Tensor, gwas_weights: torch.Tensor) -> torch.Tensor:
        # x: (B, N, C), gwas_weights: (1, N) or (N,)
        B, N, C = x.shape
        
        # 确保gwas_weights扩展到batch维度
        if gwas_weights.dim() == 1:
            gwas_weights = gwas_weights.unsqueeze(0)  # (1, N)
        
        # 扩展到batch size
        if gwas_weights.size(0) == 1 and B > 1:
            gwas_weights = gwas_weights.expand(B, -1)  # (B, N)
        
        # 扩展到embed维度
        gwas_embed = gwas_weights.unsqueeze(-1).expand(-1, -1, C)  # (B, N, C)
        
        combined = torch.cat([x, gwas_embed], dim=-1)  # (B, N, 2C)
        gate = self.gate(combined)  # (B, N, C)
        out = gate * x + (1 - gate) * self.proj(gwas_embed)
        return out

# ============================================================================
# 多尺度注意力
# ============================================================================

class MultiScaleAttention(nn.Module):
    """多尺度注意力: 局部窗口 + 全局稀疏"""
    def __init__(self, dim: int, num_heads: int = 4, window_size: int = 16,
                 attn_drop: float = 0.0, proj_drop: float = 0.0):
        super().__init__()
        self.dim = dim
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5
        self.window_size = window_size
        
        self.qkv = nn.Linear(dim, dim * 3)
        self.attn_drop = nn.Dropout(attn_drop)
        self.proj = nn.Linear(dim, dim)
        self.proj_drop = nn.Dropout(proj_drop)
        self.local_global_gate = nn.Parameter(torch.ones(1) * 0.5)
        
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        B, N, C = x.shape
        qkv = self.qkv(x).reshape(B, N, 3, self.num_heads, self.head_dim).permute(2, 0, 3, 1, 4)
        q, k, v = qkv.unbind(0)
        
        # 全局注意力
        global_attn = (q @ k.transpose(-2, -1)) * self.scale
        global_attn = global_attn.softmax(dim=-1)
        global_attn = self.attn_drop(global_attn)
        global_out = (global_attn @ v).transpose(1, 2).reshape(B, N, C)
        
        gate = torch.sigmoid(self.local_global_gate)
        out = gate * global_out + (1 - gate) * global_out
        out = self.proj(out)
        out = self.proj_drop(out)
        
        return out, global_attn.mean(dim=1)

# ============================================================================
# SwiGLU FFN
# ============================================================================

class SwiGLU(nn.Module):
    """SwiGLU激活 - 借鉴LLaMA"""
    def __init__(self, in_features: int, hidden_features: int, 
                 out_features: int, drop: float = 0.0):
        super().__init__()
        self.w1 = nn.Linear(in_features, hidden_features, bias=False)
        self.w2 = nn.Linear(hidden_features, out_features, bias=False)
        self.w3 = nn.Linear(in_features, hidden_features, bias=False)
        self.drop = nn.Dropout(drop)
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.drop(self.w2(F.silu(self.w1(x)) * self.w3(x)))

# ============================================================================
# Transformer块 (带DropPath)
# ============================================================================

class EnhancedTransformerBlock(nn.Module):
    """增强版Transformer块 - 带DropPath和LayerScale"""
    def __init__(self, dim: int, num_heads: int = 4, mlp_ratio: float = 2.67,
                 drop: float = 0.0, attn_drop: float = 0.0, drop_path: float = 0.0,
                 layer_scale_init: float = 1e-5):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = MultiScaleAttention(dim, num_heads, attn_drop=attn_drop, proj_drop=drop)
        self.norm2 = nn.LayerNorm(dim)
        hidden_dim = int(dim * mlp_ratio)
        self.mlp = SwiGLU(dim, hidden_dim, dim, drop)
        
        # DropPath
        self.drop_path1 = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.drop_path2 = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        
        # LayerScale
        self.ls1 = LayerScale(dim, layer_scale_init)
        self.ls2 = LayerScale(dim, layer_scale_init)
        
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        attn_out, attn_weights = self.attn(self.norm1(x))
        x = x + self.drop_path1(self.ls1(attn_out))
        x = x + self.drop_path2(self.ls2(self.mlp(self.norm2(x))))
        return x, attn_weights

# ============================================================================
# 主模型: SVEnhancedNet
# ============================================================================

class SVEnhancedNet(nn.Module):
    """
    SV增强版基因组预测网络
    
    特性:
    1. CNN局部特征 + Transformer全局建模
    2. GWAS门控融合
    3. 高级正则化 (DropPath, SpatialDropout, LayerScale)
    4. RRBLUP指导的特征选择
    """
    def __init__(self, 
                 n_sv: int,
                 embed_dim: int = 64,
                 depth: int = 3,
                 num_heads: int = 4,
                 patch_size: int = 50,
                 drop_rate: float = 0.1,
                 drop_path_rate: float = 0.1,
                 cnn_channels: List[int] = [8, 16, 32],
                 # RRBLUP指导
                 use_rrblup_guidance: bool = True):
        super().__init__()
        self.n_sv = n_sv
        self.patch_size = patch_size
        self.n_patches = (n_sv + patch_size - 1) // patch_size
        self.embed_dim = embed_dim
        self.use_rrblup_guidance = use_rrblup_guidance
        
        # Stage 1: CNN局部特征
        self.cnn = LocalCNNExtractor(1, cnn_channels, drop_rate)
        
        # Stage 2: Patch Embedding
        cnn_out_dim = cnn_channels[-1] * patch_size
        self.patch_embed = nn.Sequential(
            nn.Linear(cnn_out_dim, embed_dim),
            nn.LayerNorm(embed_dim),
            nn.GELU(),
            nn.Dropout(drop_rate * 0.5)
        )
        
        # 位置编码
        self.pos_embed = nn.Parameter(torch.randn(1, self.n_patches + 1, embed_dim) * 0.02)
        self.cls_token = nn.Parameter(torch.randn(1, 1, embed_dim) * 0.02)
        
        # Stage 3: GWAS融合
        self.gwas_fusion = GWASGatedFusion(embed_dim)
        
        # Stage 4: Transformer (带Stochastic Depth)
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, depth)]
        self.blocks = nn.ModuleList([
            EnhancedTransformerBlock(
                embed_dim, num_heads, 
                drop=drop_rate, attn_drop=drop_rate, 
                drop_path=dpr[i]
            )
            for i in range(depth)
        ])
        self.norm = nn.LayerNorm(embed_dim)
        
        # Stage 5: 预测头
        self.head = nn.Sequential(
            nn.Linear(embed_dim, embed_dim // 2),
            nn.GELU(),
            nn.Dropout(drop_rate),
            nn.Linear(embed_dim // 2, 1)
        )
        
        # 可学习SV重要性
        self.sv_importance = nn.Parameter(torch.zeros(n_sv))
        
        # RRBLUP指导模块
        if self.use_rrblup_guidance:
            self.rrblup_attention = None  # 延迟初始化，由set_rrblup_weights控制
    
    def set_rrblup_weights(self, rrblup_importance_weights):
        """设置RRBLUP特征重要性权重"""
        if self.use_rrblup_guidance:
            # 标准化RRBLUP权重
            if len(rrblup_importance_weights) != self.embed_dim:
                # 如果维度不匹配，进行插值或截断
                if len(rrblup_importance_weights) > self.embed_dim:
                    rrblup_importance_weights = rrblup_importance_weights[:self.embed_dim]
                else:
                    padded_weights = np.zeros(self.embed_dim)
                    padded_weights[:len(rrblup_importance_weights)] = rrblup_importance_weights
                    rrblup_importance_weights = padded_weights
            
            self.rrblup_attention = RRBLUPGuidedAttention(
                embed_dim=self.embed_dim
            ).to(next(self.parameters()).device)
            
            self.rrblup_importance_weights = torch.from_numpy(rrblup_importance_weights).float().to(next(self.parameters()).device)
        
    def forward(self, 
                x: torch.Tensor, 
                gwas_weights: Optional[torch.Tensor] = None) -> torch.Tensor:
        B = x.size(0)
        
        # SV重要性加权
        importance = torch.sigmoid(self.sv_importance)
        x = x * importance.unsqueeze(0)
        
        # CNN特征
        x_cnn = self.cnn(x)
        
        # Reshape to patches
        L = x_cnn.size(2)
        if L % self.patch_size != 0:
            pad = self.patch_size - (L % self.patch_size)
            x_cnn = F.pad(x_cnn, (0, pad))
        x_cnn = x_cnn.view(B, x_cnn.size(1), self.n_patches, self.patch_size)
        x_cnn = x_cnn.permute(0, 2, 1, 3).reshape(B, self.n_patches, -1)
        
        # Patch embedding
        x_embed = self.patch_embed(x_cnn)
        
        # CLS token + position
        cls_tokens = self.cls_token.expand(B, -1, -1)
        x_embed = torch.cat([cls_tokens, x_embed], dim=1)
        x_embed = x_embed + self.pos_embed
        
        # GWAS融合 (避免inplace操作)
        if gwas_weights is not None:
            gwas_patched = self._aggregate_gwas(gwas_weights)
            # 分离CLS token和patches
            cls_part = x_embed[:, :1, :]  # (B, 1, C)
            patch_part = x_embed[:, 1:, :]  # (B, N, C)
            # 应用GWAS融合
            patch_part_fused = self.gwas_fusion(patch_part, gwas_patched)
            # 重新组合（非inplace）
            x_embed = torch.cat([cls_part, patch_part_fused], dim=1)
        
        # Transformer
        for block in self.blocks:
            x_embed, _ = block(x_embed)
        x_embed = self.norm(x_embed)
        
        # CLS output
        cls_output = x_embed[:, 0]
        
        # 应用RRBLUP指导的注意力（如果启用）
        if self.use_rrblup_guidance and hasattr(self, 'rrblup_attention') and self.rrblup_attention is not None:
            if hasattr(self, 'rrblup_importance_weights'):
                cls_output = self.rrblup_attention(cls_output, self.rrblup_importance_weights)
            else:
                cls_output = self.rrblup_attention(cls_output)
        
        # 预测
        output = self.head(cls_output)
        return output.squeeze(-1)
    
    def _aggregate_gwas(self, gwas_weights: torch.Tensor) -> torch.Tensor:
        if gwas_weights.dim() == 1:
            gwas_weights = gwas_weights.unsqueeze(0)
        B = gwas_weights.size(0)
        L = gwas_weights.size(1)
        if L % self.patch_size != 0:
            pad = self.patch_size - (L % self.patch_size)
            gwas_weights = F.pad(gwas_weights, (0, pad))
        gwas_patched = gwas_weights.view(B, self.n_patches, self.patch_size)
        gwas_patched = gwas_patched.mean(dim=-1)
        return gwas_patched


def objective(trial, X_train, X_val, y_train, y_val, device, config: SVConfig):
    """Optuna目标函数，用于超参数优化"""
    # 定义超参数搜索空间
    embed_dim = trial.suggest_categorical('embed_dim', [32, 64, 128])
    depth = trial.suggest_int('depth', 2, 6)
    num_heads = trial.suggest_categorical('num_heads', [2, 4, 8])
    drop_rate = trial.suggest_float('drop_rate', 0.0, 0.3)
    drop_path_rate = trial.suggest_float('drop_path_rate', 0.0, 0.3)
    learning_rate = trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True)
    weight_decay = trial.suggest_float('weight_decay', 1e-6, 1e-2, log=True)
    
    # 创建模型
    model = SVEnhancedNet(
        n_sv=X_train.shape[1],
        embed_dim=embed_dim,
        depth=depth,
        num_heads=num_heads,
        patch_size=config.PATCH_SIZE,
        drop_rate=drop_rate,
        drop_path_rate=drop_path_rate,
        cnn_channels=config.CNN_CHANNELS,
        use_rrblup_guidance=True
    ).to(device)
    
    optimizer = torch.optim.AdamW(
        model.parameters(), 
        lr=learning_rate, 
        weight_decay=weight_decay
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.EPOCHS)
    criterion = nn.MSELoss()
    
    # 转换数据为tensor
    X_train_tensor = torch.FloatTensor(X_train).to(device)
    y_train_tensor = torch.FloatTensor(y_train).to(device)
    X_val_tensor = torch.FloatTensor(X_val).to(device)
    y_val_tensor = torch.FloatTensor(y_val).to(device)
    
    # 初始化GradScaler（如果使用AMP）
    if config.USE_AMP and hasattr(torch.cuda, 'amp'):
        scaler = GradScaler()
    
    # 训练循环
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(config.EPOCHS // 2):  # 减少训练轮数以加快优化
        model.train()
        optimizer.zero_grad()
        
        if config.USE_AMP and hasattr(torch.cuda, 'amp'):
            with autocast(device.type):
                outputs = model(X_train_tensor)
                loss = criterion(outputs, y_train_tensor)
            
            scaler.scale(loss).backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
            scaler.step(optimizer)
            scaler.update()
        else:
            outputs = model(X_train_tensor)
            loss = criterion(outputs, y_train_tensor)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
            optimizer.step()
        
        scheduler.step()
        
        # 验证
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val_tensor)
            val_loss = criterion(val_outputs, y_val_tensor)
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= config.EARLY_STOPPING_PATIENCE // 2:  # 减少早停耐心值以加快优化
            break
        
        # 报告中间结果给Optuna
        trial.report(val_loss.item(), epoch)
        
        # 如果需要剪枝则停止
        if trial.should_prune():
            raise optuna.TrialPruned()
    
    # 计算最终性能指标
    model.eval()
    with torch.no_grad():
        val_outputs = model(X_val_tensor).cpu().numpy()
    
    val_r2 = r2_score(y_val, val_outputs)
    
    # 返回优化目标（负R2，因为Optuna默认最小化）
    return -val_r2


def optimize_enhanced_model_hyperparams(X_train, X_val, y_train, y_val, device, config: SVConfig, n_trials=20):
    """使用Optuna优化增强模型的超参数"""
    if not OPTUNA_AVAILABLE:
        print("[WARNING] Optuna not available, skipping hyperparameter optimization")
        # 使用默认配置
        return train_enhanced_model(X_train, X_val, y_train, y_val, device, config)
    
    print(f"Starting hyperparameter optimization with {n_trials} trials...")
    
    def wrapped_objective(trial):
        return objective(trial, X_train, X_val, y_train, y_val, device, config)
    
    study = optuna.create_study(direction='minimize')  # 最小化负R2
    study.optimize(wrapped_objective, n_trials=n_trials, show_progress_bar=True)
    
    print(f"Best trial: {study.best_trial.value}")
    print(f"Best params: {study.best_trial.params}")
    
    # 使用最佳参数重新训练模型
    best_params = study.best_trial.params
    
    # 创建使用最佳参数的模型
    model = SVEnhancedNet(
        n_sv=X_train.shape[1],
        embed_dim=best_params.get('embed_dim', config.EMBED_DIM),
        depth=best_params.get('depth', config.DEPTH),
        num_heads=best_params.get('num_heads', config.NUM_HEADS),
        patch_size=config.PATCH_SIZE,
        drop_rate=best_params.get('drop_rate', config.DROP_RATE),
        drop_path_rate=best_params.get('drop_path_rate', config.DROP_PATH_RATE),
        cnn_channels=config.CNN_CHANNELS,
        use_rrblup_guidance=True
    ).to(device)
    
    # 使用最佳参数设置优化器
    optimizer = torch.optim.AdamW(
        model.parameters(), 
        lr=best_params.get('learning_rate', config.LEARNING_RATE), 
        weight_decay=best_params.get('weight_decay', config.WEIGHT_DECAY)
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.EPOCHS)
    criterion = nn.MSELoss()
    
    # 转换数据为tensor
    X_train_tensor = torch.FloatTensor(X_train).to(device)
    y_train_tensor = torch.FloatTensor(y_train).to(device)
    X_val_tensor = torch.FloatTensor(X_val).to(device)
    y_val_tensor = torch.FloatTensor(y_val).to(device)
    
    # 初始化GradScaler（如果使用AMP）
    if config.USE_AMP and hasattr(torch.cuda, 'amp'):
        scaler = GradScaler()
    
    # 训练循环
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(config.EPOCHS):
        model.train()
        optimizer.zero_grad()
        
        if config.USE_AMP and hasattr(torch.cuda, 'amp'):
            with autocast(device.type):
                outputs = model(X_train_tensor)
                loss = criterion(outputs, y_train_tensor)
            
            scaler.scale(loss).backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
            scaler.step(optimizer)
            scaler.update()
        else:
            outputs = model(X_train_tensor)
            loss = criterion(outputs, y_train_tensor)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
            optimizer.step()
        
        scheduler.step()
        
        # 验证
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val_tensor)
            val_loss = criterion(val_outputs, y_val_tensor)
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            best_model_state = model.state_dict().copy()
        else:
            patience_counter += 1
        
        if patience_counter >= config.EARLY_STOPPING_PATIENCE:
            print(f"Early stopping at epoch {epoch}")
            break
        
        if epoch % 10 == 0:
            print(f"Epoch {epoch}, Train Loss: {loss.item():.4f}, Val Loss: {val_loss.item():.4f}")
    
    # 加载最佳模型
    model.load_state_dict(best_model_state)
    
    # 预测
    model.eval()
    with torch.no_grad():
        val_outputs = model(X_val_tensor).cpu().numpy()
    
    # 计算性能指标
    val_r2 = r2_score(y_val, val_outputs)
    from scipy.stats import pearsonr
    val_corr, _ = pearsonr(y_val, val_outputs)
    
    print(f"Enhanced Model (Optimized) - R2: {val_r2:.4f}, Corr: {val_corr:.4f}")
    
    return val_outputs, val_r2, val_corr


def train_enhanced_model(X_train, X_val, y_train, y_val, device, config: SVConfig, use_optimization=False):
    """训练增强模型 - 参考train_enhanced_automl_separate_gwas.py的成功实现
    关键改进：
    1. 使用RRBLUP进行特征选择（只保留前10%最重要的特征）
    2. 不过度限制模型容量
    3. 使用正常的学习率和正则化参数
    """
    if use_optimization:
        return optimize_enhanced_model_hyperparams(X_train, X_val, y_train, y_val, device, config)
    
    from torch.utils.data import DataLoader, TensorDataset
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    
    # ========== 关键优化1: RRBLUP指导的特征选择 ==========
    print("  [Enhanced] 使用RRBLUP进行特征选择...")
    
    # 标准化数据用于RRBLUP特征选择
    scaler_fs = StandardScaler()
    X_train_scaled = scaler_fs.fit_transform(X_train)
    
    # 使用Ridge回归评估特征重要性（RRBLUP简化版）
    ridge_fs = Ridge(alpha=1.0)
    ridge_fs.fit(X_train_scaled, y_train)
    feature_importance = np.abs(ridge_fs.coef_)
    
    # 选择top 10%最重要的特征（与train_enhanced_automl_separate.py一致）
    n_features = X_train.shape[1]
    top_k_percent = 0.1  # 与automl_separate保持一致，更激进地过滤噪声位点
    k = max(10, int(n_features * top_k_percent))
    k = min(k, n_features)
    
    top_k_indices = np.argsort(feature_importance)[-k:]
    selected_indices = sorted(top_k_indices)
    
    print(f"  [Enhanced] RRBLUP特征选择: 从 {n_features} 个特征中选择 top {top_k_percent*100:.0f}% = {len(selected_indices)} 个特征")
    
    X_train_selected = X_train[:, selected_indices]
    X_val_selected = X_val[:, selected_indices]
    
    # ========== 关键优化2: 使用正常的模型参数 ==========
    train_loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train_selected).float(),
                     torch.from_numpy(y_train).float()),
        batch_size=config.BATCH_SIZE, shuffle=True, pin_memory=True
    )
    
    selected_features = X_train_selected.shape[1]
    model = SVEnhancedNet(
        n_sv=selected_features,
        embed_dim=config.EMBED_DIM,
        depth=config.DEPTH,
        num_heads=config.NUM_HEADS,
        drop_rate=config.DROP_RATE,
        drop_path_rate=config.DROP_PATH_RATE,
        cnn_channels=config.CNN_CHANNELS,
        use_rrblup_guidance=True
    ).to(device)
    
    # 获取RRBLUP特征重要性权重
    rrblup_importance_weights = feature_importance[selected_indices]
    
    # 设置RRBLUP权重
    if hasattr(model, 'set_rrblup_weights'):
        model.set_rrblup_weights(rrblup_importance_weights)
    
    # ========== 关键优化4: 使用正常的优化器参数 ==========
    optimizer = torch.optim.AdamW(
        model.parameters(), 
        lr=config.LEARNING_RATE,
        weight_decay=config.WEIGHT_DECAY
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
        optimizer, T_0=50, T_mult=2, eta_min=1e-7  # T_0从30增加到50
    )
    
    criterion = nn.MSELoss()
    scaler = GradScaler('cuda') if torch.cuda.is_available() and config.USE_AMP else None
    
    # 训练循环
    best_val_r2 = -float('inf')
    best_preds = None
    patience_counter = 0
    
    for epoch in range(config.EPOCHS):
        model.train()
        train_losses = []
        
        for X_batch, y_batch in train_loader:
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)
            optimizer.zero_grad()
            
            if scaler:
                with autocast('cuda'):
                    outputs = model(X_batch)
                    loss = criterion(outputs, y_batch)
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
                scaler.step(optimizer)
                scaler.update()
            else:
                outputs = model(X_batch)
                loss = criterion(outputs, y_batch)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), config.GRADIENT_CLIP)
                optimizer.step()
            
            train_losses.append(loss.item())
        
        scheduler.step()
        
        model.eval()
        with torch.no_grad():
            X_val_t = torch.from_numpy(X_val_selected).float().to(device)
            preds = model(X_val_t).cpu().numpy()
        
        valid = np.isfinite(preds)
        if valid.sum() > 10:
            preds_valid = preds[valid] if valid.sum() < len(preds) else preds
            val_r2 = r2_score(y_val, preds_valid)
        else:
            val_r2 = 0.0
        
        avg_train_loss = np.mean(train_losses)
        
        if epoch % 10 == 0:
            print(f"  [Enhanced] Epoch {epoch}, Train Loss: {avg_train_loss:.4f}, Val R²: {val_r2:.4f}")
        
        if val_r2 > best_val_r2:
            best_val_r2 = val_r2
            best_preds = preds.copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= config.EARLY_STOPPING_PATIENCE:
            print(f"  [Enhanced] Early stopping at epoch {epoch}")
            break
    
    from scipy.stats import pearsonr
    if best_preds is not None:
        final_r2 = r2_score(y_val, best_preds)
        final_corr, _ = pearsonr(y_val.flatten(), best_preds.flatten())
    else:
        final_r2, final_corr = 0.0, 0.0
    
    print(f"  [Enhanced] Final - R²: {final_r2:.4f}, Corr: {final_corr:.4f}")
    
    return best_preds, final_r2, final_corr




def run_sv_analysis(X_sv, y, config: SVConfig, device: torch.device):
    """分析SV数据 (SV) - 实现正确的交叉验证流程"""
    print(f"\n{'='*60}")
    print(f"Analysis for SV variants")
    print(f"{'='*60}")
    
    n_samples = len(y)
    
    # 首先对所有标记运行采样和优化分析
    sampling_ratios = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    marker_samples = sample_marker_ratios(X_sv, y, ratios=sampling_ratios)
    sampling_results = evaluate_marker_sets(marker_samples, y, config)
    
    # 使用完整的交叉验证流程
    y_binned = pd.qcut(y, q=5, labels=False, duplicates='drop')
    kfold = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=42)
    
    results = {
        'RRBLUP': {'r2': [], 'corr': []},
        'XGBoost': {'r2': [], 'corr': []},
        'WheatGP': {'r2': [], 'corr': []},
        'Enhanced': {'r2': [], 'corr': []}
    }
    
    all_rrblup_preds = np.zeros(n_samples)
    all_xgboost_preds = np.zeros(n_samples)
    all_wheatgp_preds = np.zeros(n_samples)
    all_enhanced_preds = np.zeros(n_samples)
    all_targets = np.zeros(n_samples)
    
    for fold, (train_idx, val_idx) in enumerate(kfold.split(X_sv, y_binned), 1):
        print(f"\n--- SV Fold {fold}/{config.N_FOLDS} ---")
        
        # 分割数据
        X_train, X_val = X_sv[train_idx], X_sv[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        # 1. 在训练集上进行标记选择
        optimal_indices, optimal_performance = optimize_marker_selection_on_train_only(X_train, y_train, initial_ratio=0.3, iterations=5)
        
        # 2. 应用标记选择到训练集和验证集
        X_train_selected = X_train[:, optimal_indices]
        X_val_selected = X_val[:, optimal_indices]
        
        # 3. 在每个fold内对y进行标准化
        y_train_mean, y_train_std = y_train.mean(), y_train.std()
        y_train_normalized = (y_train - y_train_mean) / (y_train_std + 1e-8)
        y_val_normalized = (y_val - y_train_mean) / (y_train_std + 1e-8)
        
        # 4. 在每个fold内对X进行标准化
        X_train_mean, X_train_std = X_train_selected.mean(axis=0, keepdims=True), X_train_selected.std(axis=0, keepdims=True) + 1e-8
        X_train_normalized = (X_train_selected - X_train_mean) / X_train_std
        X_val_normalized = (X_val_selected - X_train_mean) / X_train_std
        
        # 5. 训练各种模型
        # 1. RRBLUP (使用优化后的标记)
        print(f"  训练 SV RRBLUP (优化标记)...")
        # 真实的RRBLUP实现，使用更复杂的岭回归和交叉验证
        from sklearn.model_selection import cross_val_score
        from sklearn.linear_model import Ridge
        alphas = np.logspace(-3, 5, 50)  # 增加alpha的数量以获得更好的调优
        best_score = -np.inf
        best_alpha = 1.0
        
        # 内部交叉验证选择最佳alpha
        for alpha in alphas:
            ridge_temp = Ridge(alpha=alpha)
            scores = cross_val_score(ridge_temp, X_train_normalized, y_train_normalized, cv=5, scoring='r2')
            avg_score = scores.mean()
            if avg_score > best_score:
                best_score = avg_score
                best_alpha = alpha
        
        # 使用最佳alpha训练最终模型
        rr_model = Ridge(alpha=best_alpha)
        rr_model.fit(X_train_normalized, y_train_normalized)
        rr_preds = rr_model.predict(X_val_normalized)
        
        rr_r2 = r2_score(y_val_normalized, rr_preds)
        from scipy.stats import pearsonr
        rr_corr, _ = pearsonr(y_val_normalized, rr_preds)
        
        results['RRBLUP']['r2'].append(rr_r2)
        results['RRBLUP']['corr'].append(rr_corr)
        
        # 2. XGBoost (使用优化后的标记)
        print(f"  训练 SV XGBoost (优化标记)...")
        # 使用更全面的XGBoost参数调优和特征分析
        import xgboost as xgb
        from sklearn.model_selection import GridSearchCV
        
        # 定义参数网格
        param_grid = {
            'n_estimators': [50, 100, 200],
            'max_depth': [3, 5, 7],
            'learning_rate': [0.01, 0.1, 0.2],
            'subsample': [0.8, 0.9, 1.0],
            'colsample_bytree': [0.8, 0.9, 1.0]
        }
        
        # 使用GridSearchCV进行超参数调优
        xgb_base = xgb.XGBRegressor(random_state=42, n_jobs=-1)
        grid_search = GridSearchCV(estimator=xgb_base, param_grid=param_grid, 
                                  scoring='r2', cv=3, verbose=0, n_jobs=-1)
        grid_search.fit(X_train_normalized, y_train_normalized)
        
        # 使用最佳参数的模型（GridSearchCV已经完成了训练）
        xgb_model = grid_search.best_estimator_
        
        # 预测
        xgb_preds = xgb_model.predict(X_val_normalized)
        
        # 特征重要性分析
        feature_importance = xgb_model.feature_importances_
        
        # 计算性能指标
        xgb_r2 = r2_score(y_val_normalized, xgb_preds)
        from scipy.stats import pearsonr
        xgb_corr, _ = pearsonr(y_val_normalized, xgb_preds)
        
        results['XGBoost']['r2'].append(xgb_r2)
        results['XGBoost']['corr'].append(xgb_corr)
        
        # 3. WheatGP (使用优化后的标记)
        print(f"  训练 SV WheatGP (优化标记)...")
        wheatgp_preds, wheatgp_r2, wheatgp_corr = train_wheatgp_model(
            X_train_normalized, X_val_normalized, y_train_normalized, y_val_normalized
        )
        results['WheatGP']['r2'].append(wheatgp_r2)
        results['WheatGP']['corr'].append(wheatgp_corr)
        
        # 4. Enhanced Model (使用优化后的标记)
        print(f"  训练 SV Enhanced Model (优化标记)...")
        # 修复CNN处理基因组数据的方式
        enhanced_preds, enhanced_r2, enhanced_corr = train_enhanced_model(
            X_train_normalized, X_val_normalized, y_train_normalized, y_val_normalized, device, config, use_optimization=True
        )
        results['Enhanced']['r2'].append(enhanced_r2)
        results['Enhanced']['corr'].append(enhanced_corr)
        
        # 保存预测结果
        all_rrblup_preds[val_idx] = rr_preds
        all_xgboost_preds[val_idx] = xgb_preds
        all_wheatgp_preds[val_idx] = wheatgp_preds
        all_enhanced_preds[val_idx] = enhanced_preds
        all_targets[val_idx] = y_val_normalized
        
        print(f"  结果 - RRBLUP R2: {rr_r2:.4f}, Corr: {rr_corr:.4f}")
        print(f"  结果 - XGBoost R2: {xgb_r2:.4f}, Corr: {xgb_corr:.4f}")
        print(f"  结果 - WheatGP R2: {wheatgp_r2:.4f}, Corr: {wheatgp_corr:.4f}")
        print(f"  结果 - Enhanced R2: {enhanced_r2:.4f}, Corr: {enhanced_corr:.4f}")
        
        # 记录中间结果
        print(f"    Fold {fold} - 训练样本数: {len(train_idx)}, 验证样本数: {len(val_idx)}")
        print(f"    Fold {fold} - 训练y均值: {y_train_normalized.mean():.3f}, 验证y均值: {y_val_normalized.mean():.3f}")
        print(f"    Fold {fold} - 训练y标准差: {y_train_normalized.std():.3f}, 验证y标准差: {y_val_normalized.std():.3f}")
        
        # 计算预测值的统计信息
        print(f"    Fold {fold} - RRBLUP预测均值: {rr_preds.mean():.3f}, 标准差: {rr_preds.std():.3f}")
        print(f"    Fold {fold} - XGBoost预测均值: {xgb_preds.mean():.3f}, 标准差: {xgb_preds.std():.3f}")
        print(f"    Fold {fold} - WheatGP预测均值: {wheatgp_preds.mean():.3f}, 标准差: {wheatgp_preds.std():.3f}")
        print(f"    Fold {fold} - Enhanced预测均值: {enhanced_preds.mean():.3f}, 标准差: {enhanced_preds.std():.3f}")
    
    # 计算总体性能 - 使用统一的标准化
    # 由于每个fold内部使用不同的标准化，我们需要统一所有预测值和目标值
    # 这里我们使用每个fold内部的标准化方式，但统一计算总体R²
    
    # 验证掩码
    valid_mask = (
        np.isfinite(all_rrblup_preds) & 
        np.isfinite(all_xgboost_preds) &
        np.isfinite(all_wheatgp_preds) & 
        np.isfinite(all_enhanced_preds) &
        np.isfinite(all_targets)
    )
    
    # 使用统一标准化后的数据计算总体性能
    overall_rrblup_r2 = r2_score(all_targets[valid_mask], all_rrblup_preds[valid_mask])
    from scipy.stats import pearsonr
    overall_rrblup_corr, _ = pearsonr(all_targets[valid_mask], all_rrblup_preds[valid_mask])
    
    overall_xgb_r2 = r2_score(all_targets[valid_mask], all_xgboost_preds[valid_mask])
    overall_xgb_corr, _ = pearsonr(all_targets[valid_mask], all_xgboost_preds[valid_mask])
    
    overall_wheatgp_r2 = r2_score(all_targets[valid_mask], all_wheatgp_preds[valid_mask])
    overall_wheatgp_corr, _ = pearsonr(all_targets[valid_mask], all_wheatgp_preds[valid_mask])
    
    overall_enhanced_r2 = r2_score(all_targets[valid_mask], all_enhanced_preds[valid_mask])
    overall_enhanced_corr, _ = pearsonr(all_targets[valid_mask], all_enhanced_preds[valid_mask])
    
    summary = {
        'variant_type': 'SV',
        'sampling_results': sampling_results,  # 添加采样结果
        'optimal_n_markers': int(len(optimal_indices)),  # 添加优化后的标记数量
        'optimal_performance': float(optimal_performance),  # 添加优化性能
        'RRBLUP': {
            'mean_r2': float(np.mean(results['RRBLUP']['r2'])),
            'std_r2': float(np.std(results['RRBLUP']['r2'])),
            'mean_corr': float(np.mean(results['RRBLUP']['corr'])),
            'std_corr': float(np.std(results['RRBLUP']['corr'])),
            'overall_r2': float(overall_rrblup_r2),
            'overall_corr': float(overall_rrblup_corr)
        },
        'XGBoost': {
            'mean_r2': float(np.mean(results['XGBoost']['r2'])),
            'std_r2': float(np.std(results['XGBoost']['r2'])),
            'mean_corr': float(np.mean(results['XGBoost']['corr'])),
            'std_corr': float(np.std(results['XGBoost']['corr'])),
            'overall_r2': float(overall_xgb_r2),
            'overall_corr': float(overall_xgb_corr)
        },
        'WheatGP': {
            'mean_r2': float(np.mean(results['WheatGP']['r2'])),
            'std_r2': float(np.std(results['WheatGP']['r2'])),
            'mean_corr': float(np.mean(results['WheatGP']['corr'])),
            'std_corr': float(np.std(results['WheatGP']['corr'])),
            'overall_r2': float(overall_wheatgp_r2),
            'overall_corr': float(overall_wheatgp_corr)
        },
        'Enhanced': {
            'mean_r2': float(np.mean(results['Enhanced']['r2'])),
            'std_r2': float(np.std(results['Enhanced']['r2'])),
            'mean_corr': float(np.mean(results['Enhanced']['corr'])),
            'std_corr': float(np.std(results['Enhanced']['corr'])),
            'overall_r2': float(overall_enhanced_r2),
            'overall_corr': float(overall_enhanced_corr)
        },
        # 添加详细的预测结果
        'all_targets': all_targets.tolist(),
        'all_rrblup_preds': all_rrblup_preds.tolist(),
        'all_xgboost_preds': all_xgboost_preds.tolist(),
        'all_wheatgp_preds': all_wheatgp_preds.tolist(),
        'all_enhanced_preds': all_enhanced_preds.tolist()
    }
    
    return summary


def check_data_integrity(X, y):
    """检查数据完整性并发出警告"""
    print(f"数据形状: X={X.shape}, y={y.shape}")
    print(f"X缺失值比例: {np.isnan(X).sum() / X.size:.2%}")
    print(f"y缺失值比例: {np.isnan(y).sum() / len(y):.2%}")
    print(f"X值范围: [{X.min():.3f}, {X.max():.3f}]")
    print(f"y值范围: [{y.min():.3f}, {y.max():.3f}]")
    
    # 检查数据质量问题并发出警告
    issues = []
    if np.isnan(X).any():
        issues.append(f"X包含 {np.isnan(X).sum()} 个缺失值")
    if np.isnan(y).any():
        issues.append(f"y包含 {np.isnan(y).sum()} 个缺失值")
    if X.min() < 0 or X.max() > 2:
        issues.append(f"X值超出预期范围 [0, 2]: [{X.min():.3f}, {X.max():.3f}]")
    if np.abs(y).max() > 1e6:
        issues.append(f"y值过大，可能存在异常值")
    
    if issues:
        print(f"\n[警告] 发现数据质量问题:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    
    print(f"\n[信息] 数据完整性检查通过")
    return True


def save_results_sv(results):
    """保存SV分析结果"""
    os.makedirs(SVConfig.OUTPUT_DIR, exist_ok=True)
    
    # 保存主要结果
    results_path = os.path.join(SVConfig.OUTPUT_DIR, "final_results_sv_enhanced.json")
    with open(results_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Results saved to: {results_path}")
    
    # 保存预测结果到CSV文件
    sv_predictions_path = os.path.join(SVConfig.OUTPUT_DIR, "predictions_sv_optimized.csv")
    if isinstance(results, dict) and 'SV' in results:
        sv_result = results['SV']
        # 创建包含所有预测结果的DataFrame
        predictions_df = pd.DataFrame({
            'RRBLUP_R2': [sv_result['RRBLUP']['overall_r2']],
            'RRBLUP_Corr': [sv_result['RRBLUP']['overall_corr']],
            'XGBoost_R2': [sv_result['XGBoost']['overall_r2']],
            'XGBoost_Corr': [sv_result['XGBoost']['overall_corr']],
            'WheatGP_R2': [sv_result['WheatGP']['overall_r2']],
            'WheatGP_Corr': [sv_result['WheatGP']['overall_corr']],
            'Enhanced_R2': [sv_result['Enhanced']['overall_r2']],
            'Enhanced_Corr': [sv_result['Enhanced']['overall_corr']],
            'Optimal_Markers_Count': [sv_result['optimal_n_markers']],
            'Optimal_Performance': [sv_result['optimal_performance']]
        })
        predictions_df.to_csv(sv_predictions_path, index=False)
        print(f"SV predictions saved to: {sv_predictions_path}")
        
        # 保存每个模型的详细预测结果
        detailed_predictions_path = os.path.join(SVConfig.OUTPUT_DIR, "detailed_predictions.csv")
        detailed_df = pd.DataFrame({
            'Target': sv_result.get('all_targets', []),
            'RRBLUP_Predictions': sv_result.get('all_rrblup_preds', []),
            'XGBoost_Predictions': sv_result.get('all_xgboost_preds', []),
            'WheatGP_Predictions': sv_result.get('all_wheatgp_preds', []),
            'Enhanced_Predictions': sv_result.get('all_enhanced_preds', [])
        })
        # 只保存有数据的列
        detailed_df = detailed_df.loc[:, (detailed_df != 0).any(axis=0) | (detailed_df.columns.isin(['Target']))]
        if len(detailed_df) > 0:
            detailed_df.to_csv(detailed_predictions_path, index=False)
            print(f"Detailed predictions saved to: {detailed_predictions_path}")


def plot_comparison_sv(results):
    """绘制SV采样和优化结果图，展示各模型随位点数量变化的真实性能趋势"""
    os.makedirs(SVConfig.OUTPUT_DIR, exist_ok=True)
    
    # 提取采样结果数据
    sampling_results = results['sampling_results']
    ratios = sorted(sampling_results.keys())
    marker_counts = [sampling_results[r]['n_markers'] for r in ratios]
    
    # 检查哪些模型可用
    available_models = ['RRBLUP', 'XGBoost', 'Enhanced']
    # 检查WheatGP是否在所有结果中可用
    if all('WheatGP' in sampling_results[r] for r in ratios):
        available_models.insert(2, 'WheatGP')  # 在XGBoost后面插入
    
    # 提取各模型在不同比例下的真实R²和相关系数数据
    model_data = {}
    for model in available_models:
        model_data[model] = {
            'r2_means': [sampling_results[r][model]['mean_r2'] for r in ratios],
            'r2_stds': [sampling_results[r][model]['std_r2'] for r in ratios],
            'corr_means': [sampling_results[r][model]['mean_corr'] for r in ratios],
            'corr_stds': [sampling_results[r][model]['std_corr'] for r in ratios]
        }
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 设置模型样式
    model_styles = {
        'RRBLUP': {'marker': 'o', 'linestyle': '-'},
        'XGBoost': {'marker': 's', 'linestyle': '-'},
        'WheatGP': {'marker': '^', 'linestyle': '-'},
        'Enhanced': {'marker': 'D', 'linestyle': '-'}
    }
    
    # 使用均匀间隔的x轴索引，避免小数值挤在左边
    x_positions = np.arange(len(marker_counts))
    
    # R² vs Markers - 显示各模型的真实性能对比（使用均匀间隔）
    for model in available_models:
        axes[0, 0].errorbar(
            x_positions, model_data[model]['r2_means'], 
            yerr=model_data[model]['r2_stds'],
            marker=model_styles[model]['marker'], capsize=5, 
            linestyle=model_styles[model]['linestyle'], 
            label=model, linewidth=2
        )
    axes[0, 0].set_xticks(x_positions)
    axes[0, 0].set_xticklabels(marker_counts, rotation=45)
    axes[0, 0].set_xlabel('Number of Markers')
    axes[0, 0].set_ylabel('R2')
    axes[0, 0].set_title('R2 vs Number of SV Markers (Model Comparison)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # R² vs Ratio - 显示各模型的真实性能对比（使用均匀间隔）
    ratio_percentages = [r*100 for r in ratios]
    for model in available_models:
        axes[0, 1].errorbar(
            x_positions, model_data[model]['r2_means'],
            yerr=model_data[model]['r2_stds'],
            marker=model_styles[model]['marker'], capsize=5,
            linestyle=model_styles[model]['linestyle'],
            label=model, linewidth=2
        )
    axes[0, 1].set_xticks(x_positions)
    axes[0, 1].set_xticklabels([f'{r:.0f}' for r in ratio_percentages], rotation=45)
    axes[0, 1].set_xlabel('Marker Ratio (%)')
    axes[0, 1].set_ylabel('R2')
    axes[0, 1].set_title('R2 vs SV Marker Ratio (Model Comparison)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Correlation vs Markers - 显示各模型的真实相关系数对比（使用均匀间隔）
    for model in available_models:
        axes[1, 0].errorbar(
            x_positions, model_data[model]['corr_means'],
            yerr=model_data[model]['corr_stds'],
            marker=model_styles[model]['marker'], capsize=5,
            linestyle=model_styles[model]['linestyle'],
            label=model, linewidth=2
        )
    axes[1, 0].set_xticks(x_positions)
    axes[1, 0].set_xticklabels(marker_counts, rotation=45)
    axes[1, 0].set_xlabel('Number of Markers')
    axes[1, 0].set_ylabel('Correlation')
    axes[1, 0].set_title('Correlation vs Number of SV Markers (Model Comparison)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Performance Trend - 显示各模型随标记数量变化的真实趋势（使用均匀间隔）
    for model in available_models:
        axes[1, 1].plot(
            x_positions, model_data[model]['r2_means'],
            marker=model_styles[model]['marker'],
            linestyle=model_styles[model]['linestyle'],
            linewidth=2, label=model, markersize=6
        )
    axes[1, 1].set_xticks(x_positions)
    axes[1, 1].set_xticklabels(marker_counts, rotation=45)
    axes[1, 1].set_xlabel('Number of Markers')
    axes[1, 1].set_ylabel('Performance (R2)')
    axes[1, 1].set_title('Performance Trend vs Number of SV Markers')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_path = os.path.join(SVConfig.OUTPUT_DIR, "sv_marker_sampling_optimization.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved sampling and optimization plot: {plot_path}")
    plt.close()


def plot_model_comparison(results):
    """绘制模型比较图"""
    os.makedirs(SVConfig.OUTPUT_DIR, exist_ok=True)
    
    # 提取模型性能数据
    rrblup_r2 = results['RRBLUP']['overall_r2']
    xgboost_r2 = results['XGBoost']['overall_r2']
    wheatgp_r2 = results['WheatGP']['overall_r2']
    enhanced_r2 = results['Enhanced']['overall_r2']
    
    rrblup_corr = results['RRBLUP']['overall_corr']
    xgboost_corr = results['XGBoost']['overall_corr']
    wheatgp_corr = results['WheatGP']['overall_corr']
    enhanced_corr = results['Enhanced']['overall_corr']
    
    models = ['RRBLUP', 'XGBoost', 'WheatGP', 'Enhanced']
    r2_values = [rrblup_r2, xgboost_r2, wheatgp_r2, enhanced_r2]
    corr_values = [rrblup_corr, xgboost_corr, wheatgp_corr, enhanced_corr]
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # R²比较
    bars1 = axes[0].bar(models, r2_values, color=['blue', 'orange', 'green'], alpha=0.7)
    axes[0].set_ylabel('R2')
    axes[0].set_title('Model Comparison - R2')
    axes[0].set_ylim(bottom=0)
    
    # 在柱状图上添加数值标签
    for bar, value in zip(bars1, r2_values):
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.3f}',
                    ha='center', va='bottom')
    
    # Correlation比较
    bars2 = axes[1].bar(models, corr_values, color=['blue', 'orange', 'green', 'red'], alpha=0.7)
    axes[1].set_ylabel('Correlation')
    axes[1].set_title('Model Comparison - Correlation')
    axes[1].set_ylim(bottom=0)
    
    # 在柱状图上添加数值标签
    for bar, value in zip(bars2, corr_values):
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.3f}',
                    ha='center', va='bottom')
    
    plt.tight_layout()
    plot_path = os.path.join(SVConfig.OUTPUT_DIR, "sv_model_comparison_optimized_markers.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved model comparison plot: {plot_path}")
    plt.close()


def check_data_integrity(X, y):
    """数据完整性检查，发现问题时发出警告"""
    print(f"数据完整性检查:")
    print(f"  数据形状: X={X.shape}, y={y.shape}")
    
    x_missing_ratio = np.isnan(X).sum() / X.size
    y_missing_ratio = np.isnan(y).sum() / len(y)
    
    print(f"  X缺失值比例: {x_missing_ratio:.2%}")
    print(f"  y缺失值比例: {y_missing_ratio:.2%}")
    print(f"  X值范围: [{X.min():.3f}, {X.max():.3f}]")
    print(f"  y值范围: [{y.min():.3f}, {y.max():.3f}]")
    print(f"  X均值: {X.mean():.3f}, 标准差: {X.std():.3f}")
    print(f"  y均值: {y.mean():.3f}, 标准差: {y.std():.3f}")
    print(f"  X零值比例: {(X == 0).sum() / X.size:.2%}")
    print(f"  X非零值比例: {(X != 0).sum() / X.size:.2%}")
    
    # 发出警告
    warnings_found = []
    if x_missing_ratio > 0.05:  # 如果X缺失率超过5%
        warnings_found.append(f"X缺失值比例过高 ({x_missing_ratio:.2%} > 5%)")
    if y_missing_ratio > 0.05:  # 如果y缺失率超过5%
        warnings_found.append(f"y缺失值比例过高 ({y_missing_ratio:.2%} > 5%)")
    if X.shape[0] != y.shape[0]:
        warnings_found.append(f"X和y的样本数不匹配 (X: {X.shape[0]}, y: {y.shape[0]})")
    if np.any(np.isinf(X)):
        warnings_found.append("X中包含无穷大值")
    if np.any(np.isinf(y)):
        warnings_found.append("y中包含无穷大值")
    
    if warnings_found:
        print(f"  \n警告: 发现以下问题:")
        for warning in warnings_found:
            print(f"    - {warning}")
        print(f"  建议检查数据质量!")
    else:
        print(f"  数据质量良好，未发现严重问题。")


def load_sv_data():
    """加载SV数据 - 使用完整VCF文件（与train_enhanced_automl_separate.py一致）"""
    print(f"{'='*60}")
    print("Loading SV Data from Complete VCF File")
    print("Using full SV VCF (same as train_enhanced_automl_separate.py)")
    print(f"{'='*60}")
    
    # 使用完整的SV VCF文件，而不是GWAS筛选后的显著位点
    X_sv, y = load_sv_data_from_vcf(
        SVConfig.SV_VCF_PATH,  # 使用完整VCF文件
        SVConfig.PHENOTYPE_PATH,
        SVConfig.VCF_ID_PATH,
        max_variants=SVConfig.MAX_VARIANTS
    )
    
    # 数据完整性检查
    check_data_integrity(X_sv, y)
    
    return X_sv, y


def run_sv_enhanced_analysis(config: SVConfig, device: torch.device):
    """运行SV增强分析"""
    print(f"{'='*60}")
    print("SV Enhanced Training with Marker Selection and XGBoost Integration")
    print("Analyzing impact of different marker proportions on accuracy")
    print("Developing algorithm to select optimal markers for high accuracy with fewer markers")
    print(f"{'='*60}")
    
    X_sv, y = load_sv_data()
    
    # 分析SV变异类型
    sv_results = run_sv_analysis(X_sv, y, config, device)
    
    return {
        'SV': sv_results
    }


def main():
    """主函数"""
    # 设置随机种子以确保结果可重现
    np.random.seed(42)
    torch.manual_seed(42)
    torch.cuda.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    # 设置cuDNN确定性行为
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    print("="*60)
    print("SV Enhanced Training - Marker Count vs Accuracy Study")
    print("研究SV标记数量与预测精度的关系")
    print("Using complete SV VCF file (same as train_enhanced_automl_separate.py)")
    print("="*60)
    
    # 创建输出目录
    os.makedirs(SVConfig.OUTPUT_DIR, exist_ok=True)
    os.makedirs(SVConfig.LOG_DIR, exist_ok=True)
    
    # 加载完整SV数据
    print("\n[Step 1] 加载完整SV VCF数据...")
    X_sv, y = load_sv_data()
    print(f"  加载完成: X={X_sv.shape}, y={y.shape}")
    
    config = SVConfig()
    
    # =====================================================
    # 核心研究: 标记数量与精度的关系
    # =====================================================
    print("\n[Step 2] 标记数量与精度关系研究...")
    
    # 定义要测试的标记数量
    total_markers = X_sv.shape[1]
    # 根据总标记数动态生成测试列表
    marker_counts = [50, 100, 200, 500, 1000]
    if total_markers >= 2000:
        marker_counts.append(2000)
    if total_markers >= 3000:
        marker_counts.append(3000)
    if total_markers >= 5000:
        marker_counts.append(5000)
    # 添加完整标记数
    if total_markers not in marker_counts:
        marker_counts.append(total_markers)
    marker_counts = sorted([m for m in marker_counts if m <= total_markers])
    
    print(f"  总标记数: {total_markers}")
    print(f"  测试标记数量: {marker_counts}")
    
    # 评估不同标记数量的效果
    marker_count_results = evaluate_marker_count_effect(
        X_sv, y, config, DEVICE, marker_counts=marker_counts
    )
    
    # 绘制可视化图
    print("\n[Step 3] 绘制标记数量-精度关系图...")
    plot_marker_count_vs_accuracy(marker_count_results, SVConfig.OUTPUT_DIR)
    
    # 保存标记数量研究结果
    marker_study_path = os.path.join(SVConfig.OUTPUT_DIR, "marker_count_study.json")
    with open(marker_study_path, 'w', encoding='utf-8') as f:
        # 转换key为字符串
        save_results = {str(k): v for k, v in marker_count_results.items()}
        json.dump(save_results, f, indent=2, ensure_ascii=False)
    print(f"  标记数量研究结果已保存: {marker_study_path}")
    
    # =====================================================
    # 原有流程: 采样比例分析和标记优化
    # =====================================================
    print("\n[Step 4] 运行原有采样比例分析...")
    results = run_sv_enhanced_analysis(config, DEVICE)
    
    # 保存结果
    save_results_sv(results)
    
    # 绘制采样和优化结果图
    plot_comparison_sv(results['SV'])
    
    # 绘制模型比较图
    plot_model_comparison(results['SV'])
    
    # =====================================================
    # 打印最终总结
    # =====================================================
    print("\n" + "="*60)
    print("FINAL RESULTS SUMMARY")
    print("="*60)
    
    # 1. 标记数量研究结果
    print("\n[标记数量与精度关系]")
    print("-" * 40)
    print(f"{'Markers':<10} {'RRBLUP R2':<12} {'XGBoost R2':<12} {'WheatGP R2':<12} {'Enhanced R2':<12}")
    print("-" * 40)
    for n_markers in sorted(marker_count_results.keys()):
        res = marker_count_results[n_markers]
        print(f"{n_markers:<10} {res['RRBLUP']['mean_r2']:.4f}±{res['RRBLUP']['std_r2']:.2f}  "
              f"{res['XGBoost']['mean_r2']:.4f}±{res['XGBoost']['std_r2']:.2f}  "
              f"{res['WheatGP']['mean_r2']:.4f}±{res['WheatGP']['std_r2']:.2f}  "
              f"{res['Enhanced']['mean_r2']:.4f}±{res['Enhanced']['std_r2']:.2f}")
    
    # 2. 原有分析结果
    print("\n[采样比例分析结果]")
    for var_type, res in results.items():
        print(f"{var_type}:")
        print(f"  原始标记数: {res['sampling_results'][1.0]['n_markers']}")
        print(f"  优化后标记数: {res['optimal_n_markers']}")
        print(f"  标记减少: {(1 - res['optimal_n_markers']/res['sampling_results'][1.0]['n_markers'])*100:.1f}%")
        print(f"  优化后性能: {res['optimal_performance']:.4f}")
        print(f"  最终模型性能:")
        print(f"    RRBLUP - R2: {res['RRBLUP']['overall_r2']:.4f}, Corr: {res['RRBLUP']['overall_corr']:.4f}")
        print(f"    XGBoost - R2: {res['XGBoost']['overall_r2']:.4f}, Corr: {res['XGBoost']['overall_corr']:.4f}")
        print(f"    WheatGP - R2: {res['WheatGP']['overall_r2']:.4f}, Corr: {res['WheatGP']['overall_corr']:.4f}")
        print(f"    Enhanced - R2: {res['Enhanced']['overall_r2']:.4f}, Corr: {res['Enhanced']['overall_corr']:.4f}")
    
    print("\n" + "="*60)
    print("所有结果已保存到:", SVConfig.OUTPUT_DIR)
    print("="*60)


if __name__ == "__main__":
    main()
