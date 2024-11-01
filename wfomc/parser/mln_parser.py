# 这段代码实现了一个用于解析 Markov 逻辑网络（MLN）的解析器，它通过 lark 库定义语法树的转换，并将其转换为 MLN 数据结构。以下是逐行解释：
# Lark 是一个 Python 的解析库，用于定义和解析上下文无关文法（CFG）。在这段代码中，Lark 用于解析 MLN 文法。
from lark import Lark
# MLN：Markov 逻辑网络的基础类
from wfomc.network.mln import MLN
# Const 和 Pred：分别表示常量和谓词。
from wfomc.fol.syntax import Const, Pred
# CardinalityConstraint：表示基数约束。
from wfomc.network.constraint import CardinalityConstraint
# CCTransfomer：用于处理基数约束语法的转换器。
from wfomc.parser.cardinality_constraints_parser import CCTransfomer
# grammar：定义了 MLN 的语法规则。FOLTransformer：用于将谓词逻辑语法树转换为相应的数据结构。
from wfomc.parser.mln_grammar import grammar
# # Rational：表示有理数的工具类。
from wfomc.utils import Rational


from wfomc.parser.fol_parser import FOLTransformer
# # MLNProblem：用于封装一个完整的 MLN 问题。
from wfomc.problems import MLNProblem

from wfomc.fol.syntax import *

class MLNTransformer(FOLTransformer, CCTransfomer): # MLNTransformer 类继承了 FOLTransformer 和 CCTransfomer，表示它可以处理谓词逻辑和基数约束的语法转换。

    def domain_elements(self, args): # 处理领域中的元素，将解析到的领域元素转换为列表形式。
        return list(args)

    def int_domain(self, args): # 处理整数类型的领域定义，将解析到的字符串转换为整数。
        return int(args[0])

    def element(self, args): # 处理领域中的单个元素，返回其值。
        return args[0].value

    def set_domain(self, args): # 处理集合类型的领域定义，将元素列表转换为集合。
        return set(args[0])

    def domain_name(self, args): # 返回领域名称。
        return args[0].value

    def domain(self, args): # 处理领域定义。如果领域是一个整数，转换为一组以领域名称加索引的字符串集合；否则直接返回领域名称和领域集合。
        domain_name, domain_spec = args
        if isinstance(domain_spec, int):
            domain_spec = set(f'{domain_name}{i}' for i in range(domain_spec))
        return (domain_name, domain_spec)

    def weighting(self, args):  # 处理软规则的权重，转换为浮点数。
        return float(args[0])

    def rules(self, args): # 处理多个规则，分别提取权重和逻辑公式。
        rules = args
        weightings = []
        formulas = []
        for w, r in rules:
            weightings.append(w)
            formulas.append(r)
        return weightings, formulas

    def rule(self, args): # 处理单个规则，返回权重和逻辑公式。
        w, r = args[0]
        return w, r

    def hard_rule(self, args): # 处理硬规则，硬规则的权重为无穷大，表示它必须严格满足。
        return float('inf'), args[0]

    def soft_rule(self, args): # 处理软规则，返回权重和逻辑公式。
        return args[0], args[1]

    def mln(self, args): # 处理整个 MLN，提取规则、领域和基数约束。只支持一个领域定义。
        rules = args[0]
        domain = args[1][1] # Only one definition domain is supported
        cardinality_constraints = args[2]

        ccs: list[tuple[dict[Pred, float], str, float]] = list() # 处理基数约束，将谓词名转换为谓词对象，并构建基数约束对象。如果没有基数约束，则返回 None。
        if len(cardinality_constraints) > 0:
            for cc in cardinality_constraints:
                new_expr = dict()
                expr, comparator, param = cc
                for pred_name, coef in expr.items():
                    pred = self.name2pred.get(pred_name, None)
                    if not pred:
                        raise ValueError(f'Predicate {pred_name} not found')
                    new_expr[pred] = coef
                ccs.append((new_expr, comparator, param))
            cardinality_constraint = CardinalityConstraint(ccs)
        else:
            cardinality_constraint = None
        # 返回解析得到的规则、领域和基数约束。
        return rules, domain, cardinality_constraint

def parse(text: str) -> MLNProblem:# parse 函数解析输入的文本，使用 Lark 库根据 mln 语法解析文本，生成语法树，并通过 MLNTransformer 进行转换。
    mln_parser = Lark(grammar,
                        start='mln')
    tree = mln_parser.parse(text)
    (rules, domain, cardinality_constraint) = MLNTransformer().transform(tree)

    return MLNProblem( # 最终返回一个 MLNProblem 对象，封装解析得到的规则、领域和基数约束。
        rules,
        domain,
        cardinality_constraint
    )
# 总结：
# 这段代码实现了一个用于解析 MLN 的解析器，能够将文本形式的 MLN 表达式解析为相应的 Markov 逻辑网络数据结构，支持软规则、硬规则、领域定义和基数约束。 ​