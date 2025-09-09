from typing import Callable
from functools import reduce

from wfomc.cell_graph import build_cell_graphs
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.network.constraint import UnaryEvidenceEncoding
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
from wfomc.utils.multinomial import MultinomialCoefficients


def incremental_wfomc(context: WFOMCContext) -> RingElement:
    formula = context.formula
    domain = context.domain
    get_weight = context.get_weight
    leq_pred = context.leq_pred
    res = Rational(0, 1)
    domain_size = len(domain)
    for cell_graph, weight in build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
    ):
        # cell_graph.show()
        cells = cell_graph.get_cells()
        n_cells = len(cells)

        def helper(cell, pc_pred, pc_ccs):
            for i, p in enumerate(pc_pred):
                if cell.is_positive(p) and pc_ccs[i] > 0:
                    return i
            return None

        if context.unary_evidence_encoding == \
                UnaryEvidenceEncoding.PC:
            pc_pred, pc_ccs = zip(*context.partition_constraint.partition)
            table = dict()
            for i, cell in enumerate(cells):
                j = helper(cell, pc_pred, pc_ccs)
                if j is None:
                    continue
                table[
                    tuple(int(k == i) for k in range(n_cells))
                ] = (
                    cell_graph.get_cell_weight(cell),
                    tuple(cc - 1 if k == j else cc
                          for k, cc in enumerate(pc_ccs))
                )
        else:
            table = dict(
                (
                    tuple(int(k == i) for k in range(n_cells)), # <-- 这是 𝐤   (ivec)
                    (
                        cell_graph.get_cell_weight(cell), # <-- 第一个分量对应 T₁(𝐤) 的权重
                        None
                    )
                )
                for i, cell in enumerate(cells)
            )

        for _ in range(domain_size - 1): # 在 “加一个新元素” 的循环里，对每个目标 cell 重新计算权重并累加到新表 table 中：
            old_table = table # 取旧 table # 保存上一轮的哈希表，键就是 1-type 配置 𝐤（长度 p 的元组），值是 (已累积权重 T_h(𝐤), 约束状态)
            table = dict() # 为下一轮新建空表 T_{h+1}
            for j, cell in enumerate(cells): # 尝试把 本轮新元素 放进第 j 个 1-type cell
                w = cell_graph.get_cell_weight(cell)
                for ivec, (w_old, old_ccs) in old_table.items(): # 尝试把 本轮新元素 放进第 j 个 1-type cell
                    if old_ccs is not None: # 如果当前 cell 不能再放元素，则 continue。否则更新剩余配额 → new_ccs
                        idx = helper(cell, pc_pred, old_ccs)
                        if idx is None:
                            continue
                        new_ccs = tuple(
                            cc - 1 if k == idx else cc
                            for k, cc in enumerate(old_ccs)
                        )
                    else:
                        new_ccs = None

                    w_new = w_old * w * reduce(
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight((cell, cells[k]))
                            ** int(ivec[k]) for k in range(n_cells)
                        ),
                        Rational(1, 1)
                    )
                    ivec = list(ivec)
                    ivec[j] += 1
                    ivec = tuple(ivec)
                    w_new = w_new + table.get(ivec, (Rational(0, 1), ()))[0]
                    table[tuple(ivec)] = (
                        w_new, new_ccs
                    )
        res = res + weight * sum(w for w, _ in table.values())

    if context.unary_evidence_encoding == \
            UnaryEvidenceEncoding.PC:
        res = res / MultinomialCoefficients.coef(
            tuple(
                i for _, i in context.partition_constraint.partition
            )
        )
    return res
