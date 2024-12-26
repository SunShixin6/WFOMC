from typing import Callable
from functools import reduce

from wfomc.cell_graph import build_cell_graphs
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula


def incremental_wfomc(formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred],
                      tuple[RingElement, RingElement]],  # 一个回调函数，它接收一个谓词 Pred 并返回一个元组，元组包含两个 RingElement 元素，表示权重。
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)  # 初始化一个 Rational(0, 1)，即 0 作为最终结果的初值。
    domain_size = len(domain)
    for cell_graph, weight in build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
    ):  # 通过 build_cell_graphs 函数生成一组 cell_graph 和对应的 weight # 从这里开始相当于对算法进行一次或者多次求和（若有多个cell_graph）
        # cell_graph.show()
        cells = cell_graph.get_cells()  # 获取 cell_graph 中的所有单元格（cells）。
        n_cells = len(cells)  # 计算单元格的数量
        domain_size = len(domain)

        table = dict( # 初始化T1，与算法中 line 1-3 对应 # T1[δ_j] = w_j
            (tuple(int(k == i) for k in range(n_cells)), # 构造 δ_j（在 j 位置为1，其余为0的元组）。
             cell_graph.get_cell_weight(cell)) # 给出 w_j，从而完成 T1[δ_j] = w_j 的初始化。
            for i, cell in enumerate(cells)
        )

        for _ in range(domain_size - 1):  # 对应从 i=2 到 n，循环 n-1 次，每次相当于 i 的增加。
            old_table = table  # 把上一轮的结果作为 T_{i-1}。
            table = dict()  # 新建一个空字典，用于存储 T_i。
            for j, cell in enumerate(cells):  # 遍历每个单元 C_j。
                w = cell_graph.get_cell_weight(cell)  # 获得当前单元 C_j 的 w_j。
                for ivec, w_old in old_table.items():  # 遍历 T_{i-1} 中的所有条目 (k_old, W_old)。
                    w_new = w_old * w * reduce( # 对应行7中 W_new 的计算。
                        lambda x, y: x * y, # 计算新的权重 w_new。首先，w_old 与当前单元格的权重 w 相乘，然后使用 reduce 函数对所有其他相关单元格的权重进行累积计算。
                        (
                            cell_graph.get_two_table_weight((cell, cells[k]))  # cell_graph.get_two_table_weight((cell, cells[k])) 获取每对单元格的权重，
                            ** int(ivec[k]) for k in range(n_cells)  # ** int(ivec[k]) 用于对权重进行幂运算，累乘得到 w_new。
                        ),
                        Rational(1, 1)
                    )

                    ivec = list(ivec)
                    ivec[j] += 1  # 对应行8中 k_new 的构造(k_new = k_old + δ_j)。
                    ivec = tuple(ivec)

                    w_new = w_new + table.get(ivec, Rational(0, 1)) # 将新的权重 w_new 与 table 中对应状态的权重相加，如果 ivec 不在 table 中，则使用默认值 Rational(0, 1)。
                    table[tuple(ivec)] = w_new # 将更新后的权重 w_new 存储到 table 中，键是 ivec。


        res = res + weight * sum(table.values()) # 将 table 中所有权重的和乘以 weight 后，加到 res 中，逐步累加计算结果。

    # if leq_pred is not None: # 这行被注释掉的代码表示，如果 leq_pred 不为 None，则将 res 乘以 domain_size 的阶乘，可能是为了调整结果，但此处被注释掉。
    #     res *= Rational(math.factorial(domain_size), 1)
    return res # 返回最终的结果 res，这是一个 RingElement 类型的值，表示经过递增计算的结果。
