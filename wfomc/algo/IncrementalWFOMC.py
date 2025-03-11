from typing import Callable, List, Tuple
from functools import reduce
import pickle
import pydot
from wfomc.utils.polynomial import coeff_dict, expand
from wfomc.cell_graph import build_cell_graphs
from wfomc.context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from collections import defaultdict



def incremental_wfomc(context: WFOMCContext,
                      formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred], Tuple[RingElement, RingElement]],
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)
    domain_size = len(domain)


    ccs: List[int] = list()
    gen_vars: List[RingElement] = list()
    if context.contain_cardinality_constraint():
        pred2var = dict((pred, var) for var, pred in context.cardinality_constraint.var2pred.items())
        constraints = context.cardinality_constraint.constraints
        for constraint in constraints:
            coeffs, comp, param = constraint
            assert len(coeffs) == 1 and comp == '='
            param = int(param)
            pred, coef = next(iter(coeffs.items()))
            assert coef == 1
            ccs.append(param)
            gen_vars.append(pred2var[pred])

    def GetCoef(poly, gen_vars):  # 定义一个辅助函数 GetCoef，用于处理多项式和生成变量。
        for degrees, coeff in coeff_dict(poly, gen_vars):  # 遍历通过 coeff_dict 函数获取的多项式的系数字典，degrees 是多项式的各项的指数，coeff 是系数。
            yield degrees, coeff  # 如果多项式的各项指数小于等于约束条件 ccs 中的相应值，则返回该项的指数和系数。

    for cell_graph, weight in build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
    ):  # 从这里开始相当于对算法进行一次或者多次求和（若有多个cell_graph） NOTE: 下面是Algo2, ccs 就是算法2里面的输入的|P|=c
        cells = cell_graph.get_cells()  # 获取cell_graph中的所有cell。返回一个列表
        n_cells = len(cells)  # 计算单元格的数量
        domain_size = len(domain)

        T = defaultdict(lambda:Rational(0,1))
        for j, cell in enumerate(cells):  # NOTE: line 2
            w_j = expand(cell_graph.get_cell_weight(cell))  # 多项式，要expand一下
            delta_j = tuple(int(i == j) for i in range(n_cells))  # 构造 δ_j（在 j 位置为1，其余为0的元组）。
            for degrees, coeff in GetCoef(w_j, gen_vars):  # 提取每一项的指数和系数，degrees是一个元组（对文章算法里面的h），如（0,0）（0,1）（1,1）。1x+y + x*y deg:(0,1) # NOTE: line 3
                T[(delta_j, degrees)] = coeff  # 这里面w_jh就是约束为css的时候对应的权重 # NOTE line 4

        for i in range(2, domain_size + 1):  # domain中剩下的元素 NOTE: line 7
            old_T = T  # 把上一轮的结果作为 T_{i-1}。
            T = defaultdict(lambda:Rational(0,1))  # 新建一个空字典，用于存储 T_i。

            for j, cell in enumerate(cells):  # 遍历每个可以放元素的cell type C_j。NOTE: line 8
                w_j = cell_graph.get_cell_weight(cell)  # 获得一个选定的单元 C_j 的 w_j权重。

                for (k_old, d_old), W_old in old_T.items():  # NOTE: line 9

                    # if W_old == 0.0:  # w_old = 0 可以跳过， # 这里在存储数据的时候，才注释掉，让w=0的也存储
                    #     continue

                    k_new = tuple(k_old[m] + (1 if m == j else 0) for m in range(n_cells))  # 新的cell config NOTE: line 10

                    symbolic_weight = expand(w_j * reduce(
                        lambda x, y: x * y, (
                            cell_graph.get_two_table_weight((cell, cells[k]))
                            ** int(k_old[k]) for k in range(n_cells)
                        ), Rational(1, 1)
                    )) # 这是选中的cell_j和所有cell的r权重成绩，表示新引入这个元素带来了哪些新权重

                    for degrees, W_delta in GetCoef(symbolic_weight, gen_vars):  # 系数coeff也就是W_k_old_j_t
                        d_new = tuple(d_old + degree for d_old, degree in zip(d_old, degrees))  # 更新d_old
                        if not all(d <= c for d, c in zip(d_new, ccs)):  # d_new不断增加，防止得到的d_new 超过给的ccs
                            continue
                        W_new = W_old * W_delta
                        T[(k_new, d_new)] = T.get((k_new, d_new), Rational(0, 1)) + W_new  # 更新T


        res_sum = sum(
            W for (k, d), W in T.items() if d == tuple(ccs)  # 对T中的元素遍历，选取满足domain_size和k的 W 求和
        )

        res += weight * res_sum

    # with open('./data.pkl', 'wb') as f:
    #     pickle.dump(data, f)
    #     print("保存成功")

    return res
