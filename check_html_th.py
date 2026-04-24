import re
with open('test_ld_output/[TestLD]_HORVU.MOREX.r3.1HG0081480/integrated_analysis.html', 'r', encoding='utf-8') as f:
    content = f.read()

# 检查thead中的th
thead_match = re.search(r'<thead>(.*?)</thead>', content, re.DOTALL)
if thead_match:
    thead = thead_match.group(1)
    ths = re.findall(r'<th[^>]*>', thead)
    print(f'thead中共{len(ths)}个th:')
    for i, th in enumerate(ths[:8]):
        print(f'  [{i}] {th[:150]}')
    has_seq = sum(1 for th in ths if 'seq-col-th' in th)
    print(f'含seq-col-th的th数量: {has_seq}')
else:
    print('未找到thead')

# 也检查全文中
all_counts = content.count('seq-col-th')
print(f'\n全文seq-col-th出现次数: {all_counts}')
# 查找前几次出现位置
idx = 0
for i in range(min(5, all_counts)):
    idx = content.find('seq-col-th', idx)
    print(f'  出现{i+1}: pos={idx}, context={content[max(0,idx-30):idx+50]!r}')
    idx += 1
