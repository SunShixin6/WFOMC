from typing import Callable, List, Tuple
from functools import reduce
from wfomc.utils.polynomial import coeff_dict, expand
from wfomc.cell_graph import build_cell_graphs
from wfomc.context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from collections import defaultdict
from data_analysis.dp_tree import *
import data_analysis.dp_tree as dp_tree

def incremental_wfomc(context: WFOMCContext,
                      formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred], Tuple[RingElement, RingElement]],
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)
    domain_size = len(domain)
    dp_tree.domain_size = domain_size

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

    def get_coef(poly, gen_vars):  # 定义一个辅助函数 get_coef，用于处理多项式和生成变量。
        for degrees, coeff in coeff_dict(poly, gen_vars):  # 遍历通过 coeff_dict 函数获取的多项式的系数字典，degrees 是多项式的各项的指数，coeff 是系数。
            yield degrees, coeff  # 如果多项式的各项指数小于等于约束条件 ccs 中的相应值，则返回该项的指数和系数。

    for cell_graph, weight in build_cell_graphs(formula, get_weight, leq_pred=leq_pred):  # 从这里开始相当于对算法进行一次或者多次求和（若有多个cell_graph）
        cells = cell_graph.get_cells()
        n_cells = len(cells)

        T = defaultdict(lambda: Rational(0, 1))
        for j, cell in enumerate(cells):
            w_j = expand(cell_graph.get_cell_weight(cell))
            k_new = tuple(int(i == j) for i in range(n_cells))  # 构造 δ_j（在 j 位置为1，其余为0的元组）。
            for d_new, W_new in get_coef(w_j, gen_vars):
                T[(d_new, k_new)] = W_new
                new_node = TreeNode(d_new, k_new, W_new, sum(k_new))
                dp_tree.node_dict[(d_new, k_new)] = new_node # 临时缓存

        for i in range(2, domain_size + 1):  # domain中剩下的元素
            old_T = T  # 把上一轮的结果作为 T_{i-1}。
            T = defaultdict(lambda: Rational(0, 1))  # 新建一个空字典，用于存储 T_i。

            for j, cell in enumerate(cells):  # 遍历每个可以放元素的cell type C_j。
                w_j = cell_graph.get_cell_weight(cell)  # 获得一个选定的单元 C_j 的 w_j权重。
                for (d_old, k_old), W_old in old_T.items():

                    # if W_old == 0.0:  # w_old = 0 可以跳过， # 这里在存储数据的时候，才注释掉，让w=0的也存储
                    #     continue

                    k_new = tuple(k_old[m] + (1 if m == j else 0) for m in range(n_cells))  # 新的cell config

                    symbolic_weight = expand(w_j * reduce(
                        lambda x, y: x * y, (
                            cell_graph.get_two_table_weight((cell, cells[k]))
                            ** int(k_old[k]) for k in range(n_cells)
                        ), Rational(1, 1)
                    ))  # 这是选中的cell_j和所有cell的r权重成绩，表示新引入这个元素带来了哪些新权重

                    for d_delta, W_delta in get_coef(symbolic_weight, gen_vars):
                        d_new = tuple(d_old + d for d_old, d in zip(d_old, d_delta))  # 更新d_old
                        if not all(d <= c for d, c in zip(d_new, ccs)):  # d_new不断增加，防止得到的d_new 超过给的ccs
                            continue
                        W_new = W_old * W_delta
                        T[(d_new, k_new)] = T.get((d_new, k_new), Rational(0, 1)) + W_new  # 更新T

                        new_node = TreeNode(d_new, k_new, W_new, sum(k_new))  # 创建新节点（父节点）
                        old_node = dp_tree.node_dict.get((d_old, k_old), None)  # 查找旧节点（子节点）
                        if old_node:
                            new_node.add_child(old_node)  # 将子节点添加到父节点的 children 列表中
                        dp_tree.node_dict[(d_new, k_new)] = new_node  # 将新节点添加到字典中

        res_sum = sum(W for (d, k), W in T.items() if d == tuple(ccs))  # 对T中的元素遍历，选取满足domain_size和k的 W 求和

        res += weight * res_sum

    # save_to_file()
    # save_to_svg()

    return res
