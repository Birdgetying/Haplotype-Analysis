# 可见内容裁剪与初始视角设计

## 问题与根因

VRN-B1 综合报告在初始范围 `12,379-12,466` 下实际只显示 25 个变异位点，
但页面打开后首屏可能完全空白，并且横向滚动区仍包含约 2,000 px 的无内容区域。

浏览器实测显示这是两个独立问题：

1. `scrollToAppliedRange()` 把基因图内部中心坐标除以 SVG 宽度后，再乘整个
   `.main-data-section` 的 `scrollWidth`。基因图宽 1170 px、容器滚动宽 3052 px
   时，该换算把初始 `scrollLeft` 设置为约 1215 px，导致基因/GWAS 面板完全
   滚出视口。
2. 范围外序列列的 `th/td` 已隐藏，`col` 也被写成零宽，但 Chrome 的固定表格
   布局仍给这些列分配固有宽度，因此 `.data-table` 仍有 3042 px。仅使用
   `display:none` 或 `visibility:collapse` 不能同时保证宽度裁剪和可见列 20 px
   对齐。

## 已确认方案

采用“精确可见表格宽度 + 真实像素滚动定位”，不删除隐藏列，也不重建单倍型
数据节点。

每次 `applyFilters()` 得到最终可见位点后：

```text
table_width = 90 + 180 + 180 + visible_variant_count * 20 + 60
```

- 90 px：Haplotype 列；
- 两个 180 px：Effect 与 Phenotype；
- 20 px：每个实际可见序列位点；
- 60 px：样本数 n。

该宽度同时写入表格的 `width`、`min-width` 和 `max-width`。隐藏列继续留在 DOM
中，以便 Apply 新范围时重新显示；它们不再通过表格自动宽度制造空白滚动区。

## 滚动规则

`scrollToAppliedRange()` 不再使用“SVG 内部比例 × 整个容器 scrollWidth”。它应：

1. 读取基因/GWAS 面板相对于主数据容器的真实水平位置；
2. 计算面板可见基因区域的像素中心；
3. 仅在内容宽于容器时，将该真实中心滚到视口中心；
4. 若裁剪后的内容不宽于容器，将 `scrollLeft` 归零。

初次加载、Apply、Reset 和外部 filter 消息共用同一规则。修改数字输入框或拖动
范围手柄仍只更新 pending 状态，点击 Apply 前不得调整表格宽度、滚动位置或重绘。

## 视图同步

现有动态宽度规则保持不变：

```text
gene_area_width = max(320 px, visible_variant_count * 20 px)
svg_total_width = 450 px + gene_area_width + 220 px
```

GWAS、基因结构、connector、序列表格和 LD 图仍以同一组 applied visible positions
为准。表格精确宽度更新必须发生在 connector 和 LD 的延迟重绘之前。

## 边界情况

- 0 个可见位点：表格宽 510 px；基因区域仍使用 320 px 最小宽度。
- 1 或 2 个位点：序列列保持每列 20 px，视图不产生无意义横向滚动。
- 23/25 个位点：内容窄于常见桌面视口时首屏从左侧直接显示，不自动滚入空白。
- 大量位点：表格和 SVG 可继续产生必要的横向滚动，真实像素中心定位仍有效。
- 浏览器缩放、打印和 SVG 导出继续使用现有动态 SVG 尺寸，不改变科学数据。

## 验证

1. 先增加失败回归测试，证明旧代码仍使用比例滚动且未设置精确表格宽度。
2. 覆盖 0、2、23、25 和大量可见位点的宽度公式。
3. 验证 pending 数字/滑块修改不会调用宽度更新；Apply/Reset 各更新一次。
4. 运行 `py_compile`、完整 `python -m unittest test_star_gene_data` 和项目快速测试。
5. 重跑真实 VRN-B1 robust-discovery HTML，并按项目规则更新
   `star_gene_validation_record.md`；本次仅改变布局，不改变负验证结论。
6. 在 Chromium 中确认初始 `scrollLeft=0`、无多千像素空白滚动区，并检查
   2/23/25 位点下 GWAS、基因结构、序列列、connector 和 LD 对齐以及控制台错误。

## 范围外事项

- 不改变 robust discovery 打分、候选排序或文献位点验证边界。
- 不改变初始范围选择算法和每列 20 px 的序列可读性。
- 不删除范围外位点数据，不重建整个表格，不修改配色或侧边栏结构。
