from typing import Callable
from functools import reduce

from wfomc.utils.polynomial import coeff_dict, expand
from wfomc.cell_graph import build_cell_graphs
from wfomc.context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from collections import defaultdict
from wfomc.utils import multinomial
from wfomc.utils.multinomial import MultinomialCoefficients

def incremental_wfomc(context: WFOMCContext,
                      formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred],
                      tuple[RingElement, RingElement]],  # 一个回调函数，它接收一个谓词 Pred 并返回一个元组，元组包含两个 RingElement 元素，表示权重。
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)  # 初始化一个 Rational(0, 1)，即 0 作为最终结果的初值。
    domain_size = len(domain)

    # 准备基数约束
    new_ccs = True
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

    n_decomposed = {n: [] for n in range(domain_size + 1)}  # # 初始化一个字典 n_decomposed，用于记录中间结果。键是从 0 到 domain_size 的整数，值是空列表。

    for cell_graph, weight in build_cell_graphs(formula, get_weight, leq_pred=leq_pred):  # 通过 build_cell_graphs 函数生成一组 cell_graph 和对应的 weight，该函数的参数包括 formula、get_weight 和可选的 leq_pred。
        # cell_graph.show()
        cells = cell_graph.get_cells()  # 通过 cell_graph.get_cells() 获取 cell_graph 中的所有单元格（cells）。
        n_cells = len(cells)  # 计算单元格的数量，并存储在 n_cells 变量中。
        domain_size = len(domain)

        # 初始化DP表
        dp_table = defaultdict(lambda: Rational(0, 1))
        initial_ivec = tuple([0] * len(cells))
        initial_ccs = tuple(ccs)
        dp_table[(initial_ivec, initial_ccs)] = Rational(1, 1)

        # 生成所有可能的ivec配置
        all_ivecs = list(multinomial(len(cells), domain_size))

        # 按ivec元素之和排序，确保依赖关系已解决
        all_ivecs.sort(key=lambda x: sum(x))


        # 下面是原始代码
        table = dict(
            (tuple(int(k == i) for k in range(n_cells)),
             cell_graph.get_cell_weight(cell))
            for i, cell in enumerate(cells)
        )  # 初始化 table 字典，其中键是元组，表示单元格在图中的状态。每个元组的元素是 int(k == i)，用于表示每个单元格的状态（例如是否为目标单元格）。字典的值是该单元格的权重，使用 cell_graph.get_cell_weight(cell) 获取。
        for _ in range(domain_size - 1):  # 进行 domain_size - 1 次循环。每次循环表示对域大小的逐步递增操作。
            old_table = table  # 将当前的 table 保存为 old_table，用于后续计算中比较和更新。
            table = dict()  # 清空并重新初始化 table 字典，用于存储更新后的表。
            for j, cell in enumerate(cells):  # 遍历每个单元格 cell，通过 enumerate 获取单元格的索引 j。
                w = cell_graph.get_cell_weight(cell)  # 获取当前单元格的权重 w，通过 cell_graph.get_cell_weight(cell) 函数获得。
                for ivec, w_old in old_table.items():  # 遍历 old_table 中的每一项，ivec 是单元格状态的索引元组，w_old 是对应的权重。
                    w_new = w_old * w * reduce(  # 计算新的权重 w_new。首先，w_old 与当前单元格的权重 w 相乘，然后使用 reduce 函数对所有其他相关单元格的权重进行累积计算。
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight((cell, cells[k]))  # cell_graph.get_two_table_weight((cell, cells[k])) 获取每对单元格的权重，
                            ** int(ivec[k]) for k in range(n_cells)  # ** int(ivec[k]) 用于对权重进行幂运算，累乘得到 w_new。
                        ),
                        Rational(1, 1)
                    )
                    ivec = list(ivec)  # 将 ivec 转换为列表，
                    ivec[j] += 1  # 修改第 j 个元素（即当前单元格的状态）
                    ivec = tuple(ivec)  # 然后将其转换回元组。

                    w_new = w_new + table.get(ivec, Rational(0, 1))  # 将新的权重 w_new 与 table 中对应状态的权重相加，如果 ivec 不在 table 中，则使用默认值 Rational(0, 1)。
                    table[tuple(ivec)] = w_new  # 将更新后的权重 w_new 存储到 table 中，键是 ivec。
        res = res + weight * sum(table.values())  # 将 table 中所有权重的和乘以 weight 后，加到 res 中，逐步累加计算结果。

    return res  # 返回最终的结果 res，这是一个 RingElement 类型的值，表示经过递增计算的结果。
