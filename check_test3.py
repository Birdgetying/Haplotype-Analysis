with open('result_test3/HORVU.MOREX.r3.1HG0080220.html', 'r', encoding='utf-8') as f:
    content = f.read()
print('文件大小:', len(content), 'bytes')

checks = [
    ('annGetTags函数', 'function annGetTags'),
    ('多标签注释', 'annGetTags(d)'),
    ('INS多标签', 'INS'),
    ('DEL多标签', 'DEL'),
    ('promoter过滤', 'promoter'),
    ('intron过滤', 'intron'),
    ('annAllowed多标签循环', 'for (var i = 0; i < tags.length'),
    ('SV勾选框', 'value="SV"'),
    ('DEL勾选框', 'value="DEL"'),
    ('INS勾选框', 'value="INS"'),
]
for name, key in checks:
    status = 'OK' if key in content else 'MISS'
    print(f'  [{status}] {name}')
