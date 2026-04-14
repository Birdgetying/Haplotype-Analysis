#!/usr/bin/env python3
"""彻底修复表格排版问题"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

file_path = 'haplotype_phenotype_analysis.py'

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# 问题1：表头中 Copy All 按钮可能被遮挡
# 解决：增加表头行高度，确保按钮完全可见
old_thead = '''        html += f\'\'\'<table class="data-table" style="width:{table_width}px;">
<thead><tr>
    <th style="width:90px;min-width:90px;max-width:90px;text-align:left;padding-left:10px;vertical-align:middle;">Haplotype</th>
    <th class="effect-cell" style="vertical-align:middle;">Effect (vs Grand Mean)</th>
    <th class="box-cell" style="vertical-align:middle;">Phenotype</th>
    <th style="width:70px;min-width:70px;max-width:70px;text-align:center;vertical-align:middle;">
        <button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:5px 10px;border-radius:3px;cursor:pointer;font-size:11px;font-weight:600;">Copy All</button>
    </th>\\n\'\'\''''

new_thead = '''        html += f\'\'\'<table class="data-table" style="width:{table_width}px;">
<thead><tr style="height:45px;">
    <th style="width:90px;min-width:90px;max-width:90px;text-align:left;padding-left:10px;vertical-align:middle;">Haplotype</th>
    <th class="effect-cell" style="vertical-align:middle;">Effect (vs Grand Mean)</th>
    <th class="box-cell" style="vertical-align:middle;">Phenotype</th>
    <th style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;">
        <button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:12px;font-weight:700;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy All</button>
    </th>\\n\'\'\''''

content = content.replace(old_thead, new_thead)

# 问题2：Copy 按钮列和 n 列重叠
# 解决：增加 Copy 列宽度到 80px，确保与 n 列有明显间隔
old_copy_btn = '''            html += f'<td style="width:70px;min-width:70px;max-width:70px;text-align:center;vertical-align:middle;"><button class="copy-samples-btn" onclick="copySamples(\\'{hap}\\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:4px 10px;border-radius:3px;cursor:pointer;font-size:10px;font-weight:500;white-space:nowrap;">Copy</button></td>\\n'
            
            html += f'<td class="n-cell">{cnt}</td></tr>\\n\''''

new_copy_btn = '''            html += f'<td style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;padding:5px;"><button class="copy-samples-btn" onclick="copySamples(\\'{hap}\\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:5px 12px;border-radius:4px;cursor:pointer;font-size:11px;font-weight:600;white-space:nowrap;transition:all 0.2s;box-shadow:0 2px 4px rgba(0,0,0,0.1);">Copy</button></td>\\n'
            
            html += f'<td class="n-cell" style="padding-left:15px;">{cnt}</td></tr>\\n\''''

content = content.replace(old_copy_btn, new_copy_btn)

# 问题3：底部坐标轴行的列对齐问题
old_bottom_row = '''        # 第4列：复制样本列（空）
        html += '<td style="border:none;"></td>\\n'
        
        # 剩余列（序列+n）
        html += f'<td colspan="{len(display_positions)+1}" style="border:none;"></td>\\n\''''

new_bottom_row = '''        # 第4列：复制样本列（空）
        html += '<td style="width:80px;min-width:80px;max-width:80px;border:none;"></td>\\n'
        
        # 第5列：n 列（空）
        html += '<td style="width:40px;min-width:40px;max-width:40px;border:none;"></td>\\n'
        
        # 剩余列（序列）
        html += f'<td colspan="{len(display_positions)}" style="border:none;"></td>\\n\''''

content = content.replace(old_bottom_row, new_bottom_row)

# 问题4：添加点击按钮的视觉反馈
old_js_copy = '''// ==================== 复制样本功能 ====================
function copySamples(hapName) {
    // 找到对应的按钮
    var btn = document.querySelector(`button.copy-samples-btn[onclick="copySamples('${hapName}')"]`);
    if (!btn) return;
    
    var samples = btn.getAttribute('data-samples');
    if (!samples) {
        alert('No samples for ' + hapName);
        return;
    }
    
    // 复制到剪贴板
    navigator.clipboard.writeText(samples).then(function() {
        // 静默成功，不显示提示
    }).catch(function(err) {
        // 降级方案：使用传统方法
        var textArea = document.createElement('textarea');
        textArea.value = samples;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        // 静默成功
    });
}'''

new_js_copy = '''// ==================== 复制样本功能 ====================
function copySamples(hapName) {
    // 找到对应的按钮
    var btn = document.querySelector(`button.copy-samples-btn[onclick="copySamples('${hapName}')"]`);
    if (!btn) return;
    
    // 添加点击反馈
    btn.style.transform = 'scale(0.95)';
    btn.style.boxShadow = '0 1px 2px rgba(0,0,0,0.2)';
    
    var samples = btn.getAttribute('data-samples');
    if (!samples) {
        alert('No samples for ' + hapName);
        btn.style.transform = '';
        btn.style.boxShadow = '';
        return;
    }
    
    // 复制到剪贴板
    navigator.clipboard.writeText(samples).then(function() {
        // 恢复按钮样式
        setTimeout(function() {
            btn.style.transform = '';
            btn.style.boxShadow = '';
        }, 150);
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = samples;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        setTimeout(function() {
            btn.style.transform = '';
            btn.style.boxShadow = '';
        }, 150);
    });
}'''

content = content.replace(old_js_copy, new_js_copy)

# 问题5：Copy All 按钮也需要点击反馈
old_js_copyall = '''    navigator.clipboard.writeText(allText).then(function() {
        // 静默成功
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = allText;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
    });
}'''

new_js_copyall = '''    var btn = document.getElementById('copyAllBtn');
    
    // 添加点击反馈
    if (btn) {
        btn.style.transform = 'scale(0.95)';
        btn.style.boxShadow = '0 1px 2px rgba(0,0,0,0.2)';
    }
    
    navigator.clipboard.writeText(allText).then(function() {
        setTimeout(function() {
            if (btn) {
                btn.style.transform = '';
                btn.style.boxShadow = '';
            }
        }, 150);
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = allText;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        setTimeout(function() {
            if (btn) {
                btn.style.transform = '';
                btn.style.boxShadow = '';
            }
        }, 150);
    });
}'''

content = content.replace(old_js_copyall, new_js_copyall)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)

print('✓ 修复完成')
print('  1. 表头高度45px，Copy All 按钮完全可见（80px宽，12px字体，加粗700）')
print('  2. Copy 按钮80px宽，与 n 列间隔15px，不重叠')
print('  3. 底部坐标轴行列对齐正确')
print('  4. 点击按钮有缩放反馈效果（scale 0.95）')
print('  5. Copy All 按钮也有点击反馈')
