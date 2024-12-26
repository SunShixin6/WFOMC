from collections import defaultdict
import functools
import pandas as pd
from typing import Callable
from functools import reduce
# from rcviz import callgraph, viz

from wfomc.cell_graph import build_cell_graphs
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from wfomc.utils import multinomial
from wfomc.utils.multinomial import MultinomialCoefficients
from wfomc.utils.polynomial import coeff_dict, expand


def incremental_wfomc_teacher(context: WFOMCContext, leq_pred) -> RingElement:  # 接收一个 WFOMCContext 对象 context，该对象包含了所有上下文信息。函数返回类型为 RingElement。
    new_ccs = True

    formula = context.formula  # 从 context 中提取公式 formula、域 domain、获取权重的函数 get_weight 和可选的谓词 leq_pred。
    domain = context.domain
    get_weight = context.get_weight
    res = Rational(0, 1)  # 初始化结果 res 为有理数 0。
    domain_size = len(domain)  # 获取 domain 的大小，表示域的大小（元素个数）

    # prepare for cardinality constraint
    ccs: list[int] = list()  # k, 约束是多少 # ccs 用于存储约束条件，
    gen_vars: list[RingElement] = list()  # 这里symbolic weight只是为了表示出现次数有多少个 # gen_vars 用于存储生成变量。
    if context.contain_cardinality_constraint() and new_ccs:  # 检查上下文是否包含基数约束，并且 new_ccs 为 True。如果是，则继续以下操作。
        pred2var = dict((pred, var) for var, pred in context.cardinality_constraint.var2pred.items())  # 根据基数约束中的 var2pred 字典创建一个新的字典 pred2var，将谓词映射到变量。
        constraints = context.cardinality_constraint.constraints  # 获取基数约束的具体约束列表。
        for constraint in constraints:  # 遍历每个基数约束。
            coeffs, comp, param = constraint  # 解包约束中的系数、比较符号和参数。
            assert len(coeffs) == 1 and comp == '='  # 确保约束中的系数只有一个，且比较符号为 =。
            param = int(param)  # 将参数转换为整数。
            pred, coef = next(iter(coeffs.items()))  # 获取系数中的第一个键值对，pred 是谓词，coef 是系数。
            assert coef == 1  # 确保系数为 1。
            ccs.append(param)  # 将约束的参数添加到 ccs 中
            gen_vars.append(pred2var[pred])  # 将谓词对应的变量添加到 gen_vars 中。

    def helper(poly, gen_vars, ccs):  # 定义一个辅助函数 helper，用于处理多项式和生成变量。
        for degrees, coeff in coeff_dict(poly, gen_vars):  # 遍历通过 coeff_dict 函数获取的多项式的系数字典，degrees 是多项式的各项的指数，coeff 是系数。
            if all(
                    degree <= cc for degree, cc in zip(degrees, ccs)
            ):
                yield degrees, coeff  # 如果多项式的各项指数小于等于约束条件 ccs 中的相应值，则返回该项的指数和系数。

    cell_graph, _ = next(
        build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
        )
    )  # 通过 build_cell_graphs 构建单元格图 cell_graph，并提取图形的权重。
    cells = cell_graph.get_cells()  # 获取图中的所有单元格。
    # print(cells)  # 打印单元格信息。

    n_decomposed = dict(
        (n, list()) for n in range(domain_size + 1)
    )  # 初始化一个字典 n_decomposed，用于记录中间结果。键是从 0 到 domain_size 的整数，值是空列表。

    # 记录了中间结果
    # def dp(old_ivec: tuple[int], old_ccs: tuple[int],
    #        old_coeff: int = Rational(1, 1),
    #        save_decomposed: bool = True) -> RingElement:
    #     decomposed = n_decomposed[sum(old_ivec)]
    #     if all(i == 0 for i in old_ivec):
    #         if all(i == 0 for i in old_ccs):
    #             return Rational(1, 1)
    #         else:
    #             return Rational(0, 1)
    #     ret = Rational(0, 1)
    #     for i, cell in enumerate(cells):
    #         if old_ivec[i] == 0:
    #             continue
    #         new_ivec = list(old_ivec)
    #         new_ivec[i] -= 1
    #         new_ivec = tuple(new_ivec)
    #         w = cell_graph.get_cell_weight(cell)
    #         mul = w * reduce(
    #             lambda x, y: x * y,
    #             (
    #                 cell_graph.get_two_table_weight(
    #                     (cell, cells[k])
    #                 )
    #                 ** int(new_ivec[k]) for k in range(len(cells))
    #             ),
    #             Rational(1, 1)
    #         )
    #         mul = expand(mul)
    #         for deg, coeff in helper(mul, gen_vars, old_ccs):
    #             new_ccs = tuple(
    #                 n - k for n, k in zip(old_ccs, deg)
    #             )
    #             new_coeff = coeff * old_coeff
    #             sub = dp(
    #                 new_ivec, new_ccs, new_coeff,
    #                 sum(new_ccs) != sum(new_ivec)
    #             )
    #             if save_decomposed:
    #                 decomposed.append(
    #                     (new_coeff, new_ivec, new_ccs)
    #                 )
    #             ret = ret + sub * coeff
    #     return ret

    # 没有记录中间结果，其实就是一个递归，
    @functools.lru_cache(maxsize=None)  # 这个dp求的是在这个1type configuration条件下，满足这个cc的weight是什么
    def dp(old_ivec: tuple[int], old_ccs: tuple[int]) -> RingElement:  # 定义一个递归函数 dp，用于计算递归值。old_ivec 是当前的配置，old_ccs 是当前的基数约束。

        ret = Rational(0, 1)  # 初始化结果为 0。

        if all(i == 0 for i in old_ivec):  # 如果所有配置和基数约束都为零，返回 1 或 0，表示递归的终止条件。
            if all(i == 0 for i in old_ccs):
                return Rational(1, 1)
            else:
                return Rational(0, 1)


        for i, cell in enumerate(cells):  # 遍历这个element可以取的cell type，即遍历所有单元格。

            if old_ivec[i] == 0:  # 如果当前配置中该单元格的数量为零，则跳过该单元格。
                continue

            new_ivec = list(old_ivec)
            new_ivec[i] -= 1  # 更新配置 new_ivec，将当前单元格的数量减 1。
            new_ivec = tuple(new_ivec)

            w = cell_graph.get_cell_weight(cell)  # 获取当前单元格的权重。
            mul = w * reduce(  # 计算当前单元格的权重与其他单元格权重的乘积。# 计算和其他已有的element之间的weight，原始的inc WFOMC就是一个值，但是这里要考虑，每一个predicate的cardinality
                lambda x, y: x * y, (cell_graph.get_two_table_weight((cell, cells[k])) ** int(new_ivec[k]) for k in range(len(cells))),
                Rational(1, 1)
            )
            mul = expand(mul)  # 展开计算结果。

            for deg, coeff in helper(mul, gen_vars, old_ccs):  # 遍历cardinality的取值情况，deg是cardinality的值，coeff是此时的weight是什么，就是model counting的值是什么， # 使用 helper 函数遍历多项式的每一项，deg 是指数，coeff 是系数。
                new_ccs = tuple(n - k for n, k in zip(old_ccs, deg))  # 现在有了CC（old_ccs现在要满足的）和已经满足有的边deg，减一下，得到新的CC（new_ccs） # 根据指数 deg 和基数约束 old_ccs 计算新的约束 new_ccs。
                sub = dp(new_ivec, new_ccs)  # 递归调用 dp 函数，计算下一个子问题的结果。
                ret = ret + sub * coeff  # 将子问题的结果加权后累加到 ret 中。
        return ret  # 返回计算结果。


    for ivec in multinomial(len(cells), domain_size):  # 遍历所有可能的配置 ivec，通过 multinomial 函数生成配置。
        sub = dp(ivec, tuple(ccs))  # 计算当前配置下的结果。
        res = res + sub  # 将子问题的结果累加到总结果 res 中。


    print(res)
    return res
