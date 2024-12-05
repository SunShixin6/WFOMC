from collections import defaultdict
import pandas as pd
from typing import Callable
from functools import reduce
from wfomc.cell_graph import build_cell_graphs
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from wfomc.utils import multinomial
from wfomc.utils.multinomial import MultinomialCoefficients
from wfomc.utils.polynomial import coeff_dict, expand

def incremental_wfomc_new2(context: WFOMCContext, leq_pred) -> RingElement:
    new_ccs = True

    formula = context.formula
    domain = context.domain
    get_weight = context.get_weight
    res = Rational(0, 1)
    domain_size = len(domain)

    # 准备基数约束
    ccs: list[int] = [] # ccs 用于存储约束条件k
    gen_vars: list[RingElement] = [] # gen_vars 用于存储生成变量symbolic weight
    if context.contain_cardinality_constraint() and new_ccs: # 检查上下文是否包含基数约束，并且 new_ccs 为 True。如果是，则继续以下操作。
        pred2var = {pred: var for var, pred in context.cardinality_constraint.var2pred.items()} # 根据基数约束中的 var2pred 字典创建一个新的字典 pred2var，将谓词映射到变量。
        constraints = context.cardinality_constraint.constraints # 获取基数约束的具体约束列表。
        for constraint in constraints: # 遍历每个基数约束。
            coeffs, comp, param = constraint # 解包约束中的系数、比较符号和参数。
            assert len(coeffs) == 1 and comp == '=' # 确保约束中的系数只有一个，且比较符号为 =。
            param = int(param) # 将参数转换为整数。
            pred, coef = next(iter(coeffs.items())) # 获取系数中的第一个键值对，pred 是谓词，coef 是系数。
            assert coef == 1 # 确保系数为 1。
            ccs.append(param) # 将约束的参数添加到 ccs 中
            gen_vars.append(pred2var[pred]) # 将谓词对应的变量添加到 gen_vars 中。

    def helper(poly, gen_vars, ccs): # 定义一个辅助函数 helper，用于处理多项式和生成变量。
        for degrees, coeff in coeff_dict(poly, gen_vars): # 遍历通过 coeff_dict 函数获取的多项式的系数字典，degrees 是多项式的各项的指数，coeff 是系数。
            if all(degree <= cc for degree, cc in zip(degrees, ccs)):
                yield degrees, coeff # 如果多项式的各项指数小于等于约束条件 ccs 中的相应值，则返回该项的指数和系数。

    cell_graph, _ = next(build_cell_graphs(formula, get_weight, leq_pred=leq_pred))
    cells = cell_graph.get_cells()
    print("单元格信息:", cells)

    n_decomposed = {n: [] for n in range(domain_size + 1)} # # 初始化一个字典 n_decomposed，用于记录中间结果。键是从 0 到 domain_size 的整数，值是空列表。



    # 初始化DP表
    dp_table = defaultdict(lambda: Rational(0, 1))
    initial_ivec = tuple([0] * len(cells))
    initial_ccs = tuple(ccs)
    dp_table[(initial_ivec, initial_ccs)] = Rational(1, 1)

    # 生成所有可能的ivec配置
    all_ivecs = list(multinomial(len(cells), domain_size))

    # 按ivec元素之和排序，确保依赖关系已解决
    all_ivecs.sort(key=lambda x: sum(x))

    for ivec in all_ivecs:
        for current_ccs in [tuple(ccs)]:
            state = (ivec, current_ccs)
            current_value = dp_table.get(state, Rational(0, 1))
            if current_value == Rational(0, 1):
                continue  # 跳过贡献为零的状态

            # 遍历每个单元格以进行状态转移
            for i, cell in enumerate(cells):
                if ivec[i] == 0:
                    continue  # 如果当前单元格计数为零，则跳过

                # 创建新的ivec，减少第i个单元格的计数
                new_ivec = list(ivec)
                new_ivec[i] -= 1
                new_ivec = tuple(new_ivec)

                # 计算转移的权重
                w = cell_graph.get_cell_weight(cell)
                mul = w
                for k, other_cell in enumerate(cells):
                    mul *= cell_graph.get_two_table_weight((cell, other_cell)) ** new_ivec[k]
                mul = expand(mul)

                # 处理多项式以获取度数和系数
                for deg, coeff in helper(mul, gen_vars, current_ccs):
                    # 计算新的基数约束
                    new_ccs = tuple(n - k for n, k in zip(current_ccs, deg))

                    # 确保基数约束不被违反
                    if any(n < 0 for n in new_ccs):
                        continue

                    # 更新DP表
                    new_state = (new_ivec, new_ccs)
                    dp_table[new_state] += current_value * coeff

    # 汇总结果
    for state, value in dp_table.items():
        ivec, ccs_state = state
        if all(ccs_state):
            res += value

    print("结果:", res)

    # 可选：保存分解结果（如果需要）
    # with open('./tmp.pkl', 'wb') as f:
    #     import pickle
    #     pickle.dump(n_decomposed, f)

    # 可选：处理和显示分解结果
    # decomposed = n_decomposed[domain_size - 1]
    # all_decomposed = {}
    # for item in decomposed:
    #     if item[-1] not in all_decomposed:
    #         tmp = {}
    #     else:
    #         tmp = all_decomposed[item[-1]]
    #     if item[1] not in tmp:
    #         tmp[item[1]] = item[0]
    #     else:
    #         tmp[item[1]] += item[0]
    #     all_decomposed[item[-1]] = tmp
    # df = []
    # for k, v in all_decomposed.items():
    #     for kk, vv in v.items():
    #         if vv != 0:
    #             df.append([k, kk, vv])
    # df = pd.DataFrame(df, columns=["ccs", "ivec", "coeff"])
    # print(df.to_latex())

    return res
