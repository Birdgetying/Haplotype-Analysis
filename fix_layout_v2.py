#!/usr/bin/env python3
import sys
sys.stdout.reconfigure(encoding='utf-8')

file_path = 'haplotype_phenotype_analysis.py'

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# 1. 修改表头 - 删除 height:60px，调整 Copy All 按钮
content = content.replace(
    '<th style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;height:60px;">\n        <button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:4px 8px;border-radius:3px;cursor:pointer;font-size:10px;">Copy All</button>',
    '<th style="width:70px;min-width:70px;max-width:70px;text-align:center;vertical-align:middle;">\n        <button id="copyAllBtn" onclick="copyAllSamples()" style="background:#3498db;color:white;border:none;padding:5px 10px;border-radius:3px;cursor:pointer;font-size:11px;font-weight:600;">Copy All</button>'
)

# 2. 删除其他表头的 height:60px
content = content.replace('vertical-align:middle;height:60px;', 'vertical-align:middle;')

# 3. 修改 Copy 按钮
content = content.replace(
    '<td style="width:80px;min-width:80px;max-width:80px;text-align:center;vertical-align:middle;"><button class="copy-samples-btn" onclick="copySamples(\'{hap}\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:4px 8px;border-radius:3px;cursor:pointer;font-size:9px;white-space:nowrap;">Copy</button></td>',
    '<td style="width:70px;min-width:70px;max-width:70px;text-align:center;vertical-align:middle;"><button class="copy-samples-btn" onclick="copySamples(\'{hap}\')" data-samples="{samples_str}" style="background:#3498db;color:white;border:none;padding:4px 10px;border-radius:3px;cursor:pointer;font-size:10px;font-weight:500;white-space:nowrap;">Copy</button></td>'
)

# 4. 简化 JavaScript - 移除"已复制"提示
# copySamples 函数
old_copy1 = """    navigator.clipboard.writeText(samples).then(function() {
        // 成功提示
        var originalText = btn.innerText;
        btn.innerText = 'Copied!';
        btn.style.background = '#27ae60';
        setTimeout(function() {
            btn.innerText = originalText;
            btn.style.background = '#3498db';
        }, 1500);
    }).catch(function(err) {
        // 降级方案：使用传统方法
        var textArea = document.createElement('textarea');
        textArea.value = samples;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        
        var originalText = btn.innerText;
        btn.innerText = 'Copied!';
        btn.style.background = '#27ae60';
        setTimeout(function() {
            btn.innerText = originalText;
            btn.style.background = '#3498db';
        }, 1500);
    });"""

new_copy1 = """    navigator.clipboard.writeText(samples).then(function() {
        // 静默成功
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = samples;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
    });"""

content = content.replace(old_copy1, new_copy1)

# copyAllSamples 函数  
old_copy2 = """    navigator.clipboard.writeText(allText).then(function() {
        var btn = document.getElementById('copyAllBtn');
        var originalText = btn.innerText;
        btn.innerText = 'Copied All!';
        btn.style.background = '#27ae60';
        setTimeout(function() {
            btn.innerText = originalText;
            btn.style.background = '#3498db';
        }, 1500);
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = allText;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        
        var btn = document.getElementById('copyAllBtn');
        var originalText = btn.innerText;
        btn.innerText = 'Copied All!';
        btn.style.background = '#27ae60';
        setTimeout(function() {
            btn.innerText = originalText;
            btn.style.background = '#3498db';
        }, 1500);
    });"""

new_copy2 = """    navigator.clipboard.writeText(allText).then(function() {
        // 静默成功
    }).catch(function(err) {
        var textArea = document.createElement('textarea');
        textArea.value = allText;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
    });"""

content = content.replace(old_copy2, new_copy2)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)

print('✓ 修改完成')
print('  - 表头: 删除 height:60px，自适应高度')
print('  - Copy All 按钮: 70px 宽，11px 字体，加粗')
print('  - Copy 按钮: 70px 宽，10px 字体，加粗')
print('  - JavaScript: 移除"已复制"提示，静默复制')
