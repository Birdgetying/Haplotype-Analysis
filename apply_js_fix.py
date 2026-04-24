#!/usr/bin/env python3
"""修复 JS annNorm 函数和添加 INS/DEL 复选框"""
import re

filepath = 'haplotype_phenotype_analysis.py'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

changes = 0

# FIX 1: annNorm 函数中 return 'other' 改为 return a (保留 INS/DEL/SV 原始类型)
# 搜索: 在 if (d.functional_ann && (a === 'indel' ... 块末尾的 return 'other';
# 具体位置: "return 'other';" 在 "return 'promoter';" 之后的那一个
old_js = """        return 'promoter';
        }
        return 'other';
    }"""
new_js = """        return 'promoter';
        }
        // 保留原始类型 (INS/DEL/SV/indel)
        return a;
    }"""

if old_js in content:
    content = content.replace(old_js, new_js, 1)
    changes += 1
    print(f"[FIX1] annNorm: return 'other' -> return a (保留 INS/DEL/SV)")
else:
    print(f"[SKIP] FIX1: pattern not found")

# FIX 2: 集成页添加 INS 和 DEL 复选框 (在 SV 复选框之前)
old_cb = """<label><input type="checkbox" class="ann-cb" value="SV" checked onchange="applyFilters()"> SV</label>"""
new_cb = """<label><input type="checkbox" class="ann-cb" value="INS" checked onchange="applyFilters()"> INS</label>
                <label><input type="checkbox" class="ann-cb" value="DEL" checked onchange="applyFilters()"> DEL</label>
                <label><input type="checkbox" class="ann-cb" value="SV" checked onchange="applyFilters()"> SV</label>"""

if old_cb in content:
    # Only replace the first occurrence (integrated page)
    content = content.replace(old_cb, new_cb, 1)
    changes += 1
    print(f"[FIX2] 集成页添加 INS/DEL 复选框")
else:
    print(f"[SKIP] FIX2: pattern not found")

# FIX 3: 多面板页 - 将 indel 复选框改为 INS + DEL + indel (兼容旧数据)
old_mp_cb = """<label><input type="checkbox" class="ann-cb-mp" value="indel" checked onchange="applyFilters()"> Indel</label>"""
new_mp_cb = """<label><input type="checkbox" class="ann-cb-mp" value="INS" checked onchange="applyFilters()"> INS</label>
                            <label><input type="checkbox" class="ann-cb-mp" value="DEL" checked onchange="applyFilters()"> DEL</label>
                            <label><input type="checkbox" class="ann-cb-mp" value="indel" checked onchange="applyFilters()"> Indel</label>"""

if old_mp_cb in content:
    content = content.replace(old_mp_cb, new_mp_cb, 1)
    changes += 1
    print(f"[FIX3] 多面板页添加 INS/DEL 复选框")
else:
    print(f"[SKIP] FIX3: pattern not found")

if changes > 0:
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"\n已保存 {changes} 个修复到 {filepath}")
else:
    print("\n没有需要修复的内容")
