import pickle
import pandas as pd

with open("../data.pkl", "rb") as f:
    n_decomposed = pickle.load(f)

domain_size = 5

decomposed = n_decomposed[domain_size - 1]  # 获取domain_size - 1对应的分解结果

# 将数据转换为适合处理的格式
processed_data = []
for item in decomposed:
    ccs = item[0][0]
    ivec = item[1]
    coeff = item[2]
    processed_data.append((ccs, ivec, coeff))

# 创建DataFrame
df = pd.DataFrame(processed_data, columns=["ccs", "ivec", "coeff"])

# 按照first_element和second_element分组，并将coeff求和
grouped_df = df.groupby(['ccs', 'ivec'], as_index=False)['coeff'].sum()

# 按照ccs和ivec降序排序
sorted_df = grouped_df.sort_values(by=['ccs', 'ivec'], ascending=[False, False])

# 打印结果
print(sorted_df.to_latex(index=False))



# all_decomposed = dict()  # 初始化一个空字典用于存储所有分解结果
# for item in decomposed:  # 遍历分解结果
#     if item[-1] not in all_decomposed:  # 如果item的最后一个元素不在all_decomposed中
#         tmp = dict()  # 初始化一个空字典
#     else:  # 否则
#         tmp = all_decomposed[item[-1]]  # 获取对应的字典
#     if item[1] not in tmp:  # 如果item的第二个元素不在tmp中
#         tmp[item[1]] = item[0]  # 将item的第一个元素赋值给tmp的对应键
#     else:  # 否则
#         tmp[item[1]] = tmp[item[1]] + item[0]  # 累加item的第一个元素到tmp的对应键
#     all_decomposed[item[-1]] = tmp  # 更新all_decomposed字典
# df = []  # 初始化一个空列表用于存储数据
# for k, v in all_decomposed.items():  # 遍历all_decomposed字典
#     for kk, vv in v.items():  # 遍历内部字典
#         if vv != 0:  # 如果值不为0
#             df.append([k, kk, vv])  # 将键和值添加到df列表中
# df = pd.DataFrame(df, columns=["ccs", "ivec", "coeff"])  # 将df列表转换为DataFrame并指定列名
# print(df.to_latex())  # 打印DataFrame的LaTeX表示
# print(df)
