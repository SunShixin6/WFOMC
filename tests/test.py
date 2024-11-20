import re
from collections import defaultdict

# 这些代码是将输出的log文件，排序打印
def sort_by_depth(file_path):
    # 用于存储每个 ki_div_Mi 的数据
    results = defaultdict(list)

    # 打开并逐行读取文件
    with open(file_path, 'r') as file:
        current_key = None  # 当前的 ki_div_Mi 值
        for line in file:
            # 匹配 ki_div_Mi 行
            ki_match = re.match(r"\s*ki_div_Mi:\s*(.+)", line)
            if ki_match:
                current_key = ki_match.group(1).strip()
                results[current_key] = []  # 初始化该 ki_div_Mi 的数据列表
            # 匹配 depth 行
            elif current_key and re.match(r"\s*depth:\s*(\d+)\s+value:\s+(.+)", line):
                depth_match = re.match(r"\s*depth:\s*(\d+)\s+value:\s+(.+)", line)
                if depth_match:
                    depth = int(depth_match.group(1))
                    value = depth_match.group(2).strip()
                    results[current_key].append((depth, value))

    # 对每个 ki_div_Mi 的数据按 depth 排序
    for key in results:
        results[key].sort(key=lambda x: x[0])  # 按 depth 排序

    return results

def write_sorted_results(results, output_path):
    # 写入排序后的数据到文件
    with open(output_path, 'w') as file:
        for ki_div_Mi, entries in results.items():
            file.write(f"ki_div_Mi:  {ki_div_Mi}\n")
            for depth, value in entries:
                file.write(f"  depth:  {depth}   value:  {value}\n")
            file.write("\n")  # 分隔不同的 ki_div_Mi

# 文件路径
input_file = "../wfomc/log.txt"  # 替换为您的文件路径
# input_file = "./log.txt"  # 替换为您的文件路径
output_file = "output.txt"

# 排序并写入结果
sorted_results = sort_by_depth(input_file)
write_sorted_results(sorted_results, output_file)

print(f"已完成排序，结果保存在 {output_file}")
