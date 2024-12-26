from collections import defaultdict
import functools
import pandas as pd
from typing import Callable, Tuple, Dict
from functools import reduce
# from rcviz import callgraph, viz

from wfomc.cell_graph import build_cell_graphs
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from wfomc.utils import multinomial
from wfomc.utils.multinomial import MultinomialCoefficients
from wfomc.utils.polynomial import coeff_dict, expand


def incremental_wfomc_teacher_dp(context: WFOMCContext, leq_pred) -> RingElement:
    new_ccs = True

    formula = context.formula
    domain = context.domain
    get_weight = context.get_weight
    res = Rational(0, 1)
    domain_size = len(domain)

    ccs: list[int] = list()
    gen_vars: list[RingElement] = list()
    if context.contain_cardinality_constraint() and new_ccs:
        pred2var = dict(
            (pred, var) for var, pred in
            context.cardinality_constraint.var2pred.items()
        )
        constraints = context.cardinality_constraint.constraints
        for constraint in constraints:
            coeffs, comp, param = constraint
            assert len(coeffs) == 1 and comp == '='
            param = int(param)
            pred, coef = next(iter(coeffs.items()))
            assert coef == 1
            ccs.append(param)
            gen_vars.append(pred2var[pred])

    def helper(poly, gen_vars, ccs):
        for degrees, coeff in coeff_dict(poly, gen_vars):
            if all(
                    degree <= cc for degree, cc in zip(degrees, ccs)
            ):
                yield degrees, coeff

    cell_graph, _ = next(
        build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
        )
    )
    cells = cell_graph.get_cells()
    print(cells)

    n_decomposed = dict(
        (n, list()) for n in range(domain_size + 1)
    )

    # 以下为迭代版本的dp函数实现, 使用一个堆栈来模拟递归过程。当我们 dp(ivec, ccs) 时，如果该状态还未计算，就将其入栈处理。在出栈时可以利用已经计算好的子问题结果来求解当前状态的值。
    dp_memo = {}  # 用于存储dp结果: dp_memo[(ivec, ccs)] = RingElement值

    def dp_iter(ivec: tuple[int], ccs: tuple[int]) -> RingElement:
        """
        定义一个函数用于返回dp(ivec, ccs)的值。如果值尚未计算，则通过堆栈方法迭代计算。
        如果已计算过该状态，直接返回已缓存结果。如果未计算过，入栈该状态，并根据转移条件寻找子状态。如果子状态也未计算过，则递归入栈子状态。如此往复直到遇到最基本的状态（全0向量时判断ccs是否全0）。
        最终从栈顶往回计算上层状态的值。
        """
        if (ivec, ccs) in dp_memo: # 若已有结果则直接返回
            return dp_memo[(ivec, ccs)]

        # 堆栈元素存放的结构：(
        #   state, # (ivec, ccs)
        #   "processing" or "resolved",
        #   用于中间存储转移所需信息的临时变量
        # )
        # 当标记为"processing"表示该状态尚在处理子问题阶段，
        # 当转移子问题全部解决后标记"resolved"以存结果。
        stack = []
        stack.append(((ivec, ccs), "processing", None)) # 将当前状态标记为 "processing"（处理中），并入栈

        while stack: # 开始模拟递归的迭代过程
            (cur_ivec, cur_ccs), status, temp_data = stack.pop() # 从栈中弹出一个状态进行处理

            if status == "processing": # 若当前状态尚未计算，需要根据状态公式进行求解
                if all(i == 0 for i in cur_ivec):
                    # 基本情况：ivec全为0时，如果ccs也全为0，则dp = 1，否则dp = 0
                    val = Rational(1, 1) if all(x == 0 for x in cur_ccs) else Rational(0, 1)
                    dp_memo[(cur_ivec, cur_ccs)] = val # 将结果存储到缓存中
                    continue

                # 对当前状态计算需要的子状态列表和中间数据准备
                ret = Rational(0, 1)
                child_states = []  # 存储当前状态的所有子状态，用于后续计算
                for i, cell in enumerate(cells): # 遍历所有单元格，尝试减少每个单元格的使用次数。 ivec[i] -> new_ivec
                    if cur_ivec[i] == 0:
                        continue
                    new_ivec = list(cur_ivec)
                    new_ivec[i] -= 1 # 减少当前单元格的使用次数
                    new_ivec = tuple(new_ivec)
                    w = cell_graph.get_cell_weight(cell) # 计算当前单元格的权重
                    mul = w * reduce( # 计算权重的乘积，包括当前单元格和其他单元格之间的权重关系
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight(
                                (cell, cells[k])
                            ) ** int(new_ivec[k]) for k in range(len(cells))
                        ),
                        Rational(1, 1)
                    )
                    mul = expand(mul)
                    for deg, coeff in helper(mul, gen_vars, cur_ccs): # 遍历多项式中满足约束的所有项
                        new_ccs = tuple(
                            n - k for n, k in zip(cur_ccs, deg)
                        )
                        # dp(new_ivec, new_ccs)是子状态，需要确保其已计算  若子状态未计算，则先处理子状态
                        if (new_ivec, new_ccs) not in dp_memo: # 如果子状态未计算，则将其添加到待处理的子状态列表中
                            child_states.append((new_ivec, new_ccs))

                if not child_states:
                    # 没有需要先计算的子状态，直接计算ret
                    # 再次遍历与上面相同的转移，以求和最终结果
                    ret = Rational(0, 1)
                    for i, cell in enumerate(cells):
                        if cur_ivec[i] == 0:
                            continue
                        new_ivec = list(cur_ivec)
                        new_ivec[i] -= 1
                        new_ivec = tuple(new_ivec)
                        w = cell_graph.get_cell_weight(cell)
                        mul = w * reduce(
                            lambda x, y: x * y,
                            (
                                cell_graph.get_two_table_weight(
                                    (cell, cells[k])
                                ) ** int(new_ivec[k]) for k in range(len(cells))
                            ),
                            Rational(1, 1)
                        )
                        mul = expand(mul)
                        for deg, coeff in helper(mul, gen_vars, cur_ccs):
                            new_ccs = tuple(
                                n - k for n, k in zip(cur_ccs, deg)
                            )
                            sub = dp_memo[(new_ivec, new_ccs)]
                            ret = ret + sub * coeff
                    dp_memo[(cur_ivec, cur_ccs)] = ret
                else:
                    # 若还有子状态未计算，先把当前状态和子状态全部压栈处理，先把当前状态再次入栈标记为"resolved"等待子状态结果
                    stack.append(((cur_ivec, cur_ccs), "resolved", None))
                    # 子状态处理
                    for s in child_states:
                        # 若子状态已存在结果，则无需入栈
                        if s not in dp_memo:
                            stack.append((s, "processing", None))

            else:
                # status == "resolved" 表示所有子状态都已计算完毕，可以计算当前状态的dp值
                ret = Rational(0, 1)
                for i, cell in enumerate(cells):
                    if cur_ivec[i] == 0:
                        continue
                    new_ivec = list(cur_ivec)
                    new_ivec[i] -= 1
                    new_ivec = tuple(new_ivec)
                    w = cell_graph.get_cell_weight(cell)
                    mul = w * reduce(
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight(
                                (cell, cells[k])
                            ) ** int(new_ivec[k]) for k in range(len(cells))
                        ),
                        Rational(1, 1)
                    )
                    mul = expand(mul)
                    for deg, coeff in helper(mul, gen_vars, cur_ccs):
                        new_ccs = tuple(
                            n - k for n, k in zip(cur_ccs, deg)
                        )
                        sub = dp_memo[(new_ivec, new_ccs)]
                        ret = ret + sub * coeff
                dp_memo[(cur_ivec, cur_ccs)] = ret

        return dp_memo[(ivec, ccs)]

    # 计算总和
    for ivec in multinomial(len(cells), domain_size):
        sub = dp_iter(ivec, tuple(ccs))
        res = res + sub
    print("--->",res)

    # with open('./tmp.pkl', 'wb') as f:
    #     import pickle
    #     pickle.dump(n_decomposed, f)
    # decomposed = n_decomposed[domain_size - 1]
    # all_decomposed = dict()
    # for item in decomposed:
    #     if item[-1] not in all_decomposed:
    #         tmp = dict()
    #     else:
    #         tmp = all_decomposed[item[-1]]
    #     if item[1] not in tmp:
    #         tmp[item[1]] = item[0]
    #     else:
    #         tmp[item[1]] = tmp[item[1]] + item[0]
    #     all_decomposed[item[-1]] = tmp
    # df = []
    # for k, v in all_decomposed.items():
    #     for kk, vv in v.items():
    #         if vv != 0:
    #             df.append([k, kk, vv])
    # df = pd.DataFrame(df, columns=["ccs", "ivec", "coeff"])
    # print(df.to_latex())

    return res

