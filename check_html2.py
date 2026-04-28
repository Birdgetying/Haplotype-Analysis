import re, json

with open('result_test/[NoExtend]_HORVU.MOREX.r3.1HG0080220/integrated_analysis.html', 'r', encoding='utf-8') as f:
    content = f.read()

# 找 variantsData
idx = content.find('var variantsData')
if idx < 0:
    print('未找到variantsData')
else:
    # 找到数组结束
    start = content.index('[', idx)
    depth = 0
    end = start
    for i, c in enumerate(content[start:], start):
        if c == '[': depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    arr_str = content[start:end]
    try:
        data = json.loads(arr_str)
        print(f"variantsData 共 {len(data)} 条")
        print()
        # 统计 annotation 和 functional_ann
        from collections import Counter
        ann_cnt = Counter(d.get('annotation','') for d in data)
        fan_cnt = Counter(d.get('functional_ann','') for d in data)
        print("annotation 分布:", dict(ann_cnt))
        print("functional_ann 分布:", dict(fan_cnt))
        print()
        # 打印前5条
        for d in data[:5]:
            print(f"  pos={d.get('pos')}, ann={d.get('annotation')}, func={d.get('functional_ann')}, maf={d.get('maf')}")
    except Exception as e:
        print(f"JSON解析失败: {e}")
        print(arr_str[:500])

# 检查 annNorm 函数内容
idx2 = content.find('function annNorm')
if idx2 > 0:
    print()
    print("=== annNorm 函数 ===")
    print(content[idx2:idx2+600])
