from typing import Callable
from functools import reduce
from wfomc.utils.polynomial import coeff_dict, expand
from wfomc.cell_graph import build_cell_graphs
from wfomc.context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula


def incremental_wfomc(context: WFOMCContext,  # 注意这里要和原始incremental_wfomc，多一个参数
                      formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred],
                      tuple[RingElement, RingElement]],
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)  # 初始化一个 Rational(0, 1)，即 0 作为最终结果的初值。
    domain_size = len(domain)

    # NOTE: 为基数约束做准备
    ccs: list[int] = list()  # 列表ccs用于存储约束是k多少
    gen_vars: list[RingElement] = list()  # 列表gen_vars 用于存储生成变量symbolic weight。
    if context.contain_cardinality_constraint():  # 检查上下文是否包含基数约束。如果是，则继续以下操作。
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

    def helper(poly, gen_vars):  # 定义一个辅助函数 helper，用于处理多项式和生成变量。
        for degrees, coeff in coeff_dict(poly, gen_vars):  # 遍历通过 coeff_dict 函数获取的多项式的系数字典，degrees 是多项式的各项的指数，coeff 是系数。
            yield degrees, coeff  # 如果多项式的各项指数小于等于约束条件 ccs 中的相应值，则返回该项的指数和系数。

    n_decomposed = dict(
        (n, list()) for n in range(domain_size + 1)
    )  # 初始化一个字典 n_decomposed，用于记录中间结果。键是从 0 到 domain_size 的整数，值是空列表。

    for cell_graph, weight in build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
    ):  # 从这里开始相当于对算法进行一次或者多次求和（若有多个cell_graph） NOTE: 下面是Algo2, ccs 就是算法2里面的输入的 k
        cells = cell_graph.get_cells()  # 获取cell_graph中的所有cell。返回一个列表
        n_cells = len(cells)  # 计算单元格的数量
        domain_size = len(domain)

        T = dict()  # 初始化T1， 对应 # T1[δ_j,h] = w_jh
        for j, cell in enumerate(cells):  # NOTE: line 2
            w_j = expand(cell_graph.get_cell_weight(cell))  # 多项式，要expand一下
            delta_j = tuple(int(i == j) for i in range(n_cells))  # 构造 δ_j（在 j 位置为1，其余为0的元组）。
            for degrees, coeff in helper(w_j, gen_vars):  # 提取每一项的指数和系数，degrees是一个元组（对文章算法里面的h），如（0,0）（0,1）（1,1）。1x+y + x*y deg:(0,1) # NOTE: line 3
                T[(delta_j, degrees)] = coeff  # 这里面w_jh就是约束为css的时候对应的权重 # NOTE line 4

        for _ in range(domain_size - 1):  # 对应从 i=2 到 n，循环 n-1 次，每次相当于 i 的增加。 NOTE: line 7
            old_T = T  # 把上一轮的结果作为 T_{i-1}。
            T = dict()  # 新建一个空字典，用于存储 T_i。

            for j, cell in enumerate(cells):  # 遍历每个单元 C_j。NOTE: line 8
                w_j = cell_graph.get_cell_weight(cell)  # 获得当前单元 C_j 的 w_j。

                for (k_old, h_old), W_old in old_T.items():  # NOTE: line 9
                    if W_old == 0.0: # w_old = 0 可以跳过，
                        continue

                    k_new = tuple(k_old[m] + (1 if m == j else 0) for m in range(n_cells))  # NOTE: line 10     k_new

                    mul = w_j * reduce(
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight((cell, cells[k]))  # 获取每对单元格的权重，
                            ** int(k_old[k]) for k in range(n_cells)  # ** int(k_old[k]) 用于对权重进行幂运算
                        ),
                        Rational(1, 1)
                    )
                    mul = expand(mul)  # 展开计算结果。

                    for degrees, W_k_old_j_t in helper(mul, gen_vars):  # 系数coeff也就是W_k_old_j_t  NOTE: line 11
                        h_new = tuple(h_old + degree for h_old, degree in zip(h_old, degrees))  # 更新h_old
                        if not all(h_new <= h_new for h_new, h_new in zip(h_new, h_new)):  # h_new不断增加，防止得到的h_new 超过给的ccs
                            continue

                        T[(k_new, h_new)] = T.get((k_new, h_new), Rational(0, 1)) + W_old * W_k_old_j_t  # NOTE: line 12

        res_sum = sum(
            W for (k_vec, h_val), W in T.items() if h_val == tuple(ccs)  # 对T中的元素遍历，选取满足domain_size和k的 W 求和
        )  # NOTE: line 17

        res += weight * res_sum
        print("res----->", res)

    return res
