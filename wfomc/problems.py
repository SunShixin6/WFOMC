from __future__ import annotations # 该行代码引入了未来特性，允许在 Python 3.7+ 中使用前向引用，即可以在类型注解中使用尚未定义的类或类型。

from wfomc.fol.sc2 import SC2, to_sc2 # SC2 和 to_sc2：用于处理 SC2 公式（约束二阶逻辑公式）。
from wfomc.fol.syntax import Const, Pred, top, AUXILIARY_PRED_NAME, \
    Formula, QuantifiedFormula, Universal, Equivalence
from wfomc.fol.utils import new_predicate
from wfomc.network.constraint import CardinalityConstraint
from wfomc.utils.polynomial import Rational # Rational：用于处理有理数。
from fractions import Fraction # Fraction 和 math：用于处理数学运算。
import math


class WFOMCProblem(object): # 该类定义了加权一阶模型计数/采样问题。它是用来描述具体的 WFOMC 问题的类。
    """
    A weighted first-order model counting problem.
    """

    def __init__(self, sentence: SC2,
                 domain: set[Const],
                 weights: dict[Pred, tuple[Rational, Rational]],
                 cardinality_constraint: CardinalityConstraint = None):
        self.domain: set[Const] = domain
        self.sentence: SC2 = sentence # sentence：SC2 公式，代表问题的逻辑句子。这个就是代表FOL某种类型问题，不需要理解了解。详情见论On exact sa...的论文
        self.weights: dict[Pred, tuple[Rational, Rational]] = weights # weights：字典，键为谓词，值为两个有理数的元组，表示加权。
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint # cardinality_constraint：可选的基数约束对象。

    def __str__(self) -> str: # __str__ 方法定义了如何将对象转换为字符串。返回领域、句子、权重以及基数约束的字符串表示。
        s = ''
        s += 'Domain: \n'
        s += '\t' + str(self.domain) + '\n'
        s += 'Sentence: \n'
        s += '\t' + str(self.sentence) + '\n'
        s += 'Weights: \n'
        s += '\t' + str(self.weights) + '\n'
        if self.cardinality_constraint is not None:
            s += 'Cardinality Constraint: \n'
            s += '\t' + str(self.cardinality_constraint) + '\n'
        return s

    def __repr__(self) -> str: # __repr__ 方法返回对象的字符串表示，直接调用了 __str__ 方法。
        return str(self)


class MLNProblem(object):# 该类定义了一个Markov 逻辑网络问题。
    """
    A Markov Logic Network problem.
    """

    def __init__(self, rules: tuple[list[tuple[Rational, Rational]], list[Formula]],
                 domain: set[Const],
                 cardinality_constraint: CardinalityConstraint):
        self.rules = rules # rules：包含规则的元组，规则由公式和权重组成。
        # self.formulas: rules[1]
        # self.formula_weights: = dict(zip(rules[1], rules[0]))
        self.domain: set[Const] = domain # domain：常量的集合，表示问题的领域
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint # cardinality_constraint：基数约束对象。


def MLN_to_WFOMC(mln: MLNProblem): # 定义了一个函数 MLN_to_WFOMC，该函数将接收一个 MLNProblem 类型的参数 mln。这个函数可能用于将一个马尔可夫逻辑网络（MLN）转换为一个加权一阶模型计数（WFOMC）问题。
    sentence = top # 初始化 sentence 为 top，这可能表示一个布尔逻辑的恒真公式，即所有条件均为真
    weightings: dict[Pred, tuple[Rational, Rational]] = dict() # 初始化了一个名为 weightings 的空字典，用来存储每个谓词的权重信息。字典的键为谓词 Pred，值为一对有理数 tuple[Rational, Rational]，这对值可能表示谓词的正反两种情况的权重。
    for weighting, formula in zip(*mln.rules): # 对于每条规则，提取其权重和公式。对于非硬规则（权重不是无穷大），引入辅助谓词，并将公式转换为其与辅助谓词的等价形式。同时，将权重添加到 weightings 字典中。
        free_vars = formula.free_vars() # 提取公式 formula 中的自由变量并存储到 free_vars 变量中。自由变量是指未绑定的变量，可能在公式中需要进一步处理。
        if weighting != float('inf'):# 检查当前规则的权重 weighting 是否为无穷大。如果权重不是无穷大，则表示该规则不是一个硬约束（硬规则），意味着其可能可以被违反。
            aux_pred = new_predicate(len(free_vars), AUXILIARY_PRED_NAME)  # 定义一个辅助谓词 aux_pred，其名称为 AUXILIARY_PRED_NAME，并且参数个数等于自由变量的数量 len(free_vars)。该辅助谓词将与原公式等价地替换，作为帮助转化规则的一部分。
            formula = Equivalence(formula, aux_pred(*free_vars)) # 将原公式 formula 替换为一个新的公式，表示公式与辅助谓词 aux_pred 之间的等价性。这是公式转换的一部分，用于将规则表示成辅助谓词的等价形式。
            # 将权重存储在 weightings 字典中。对于辅助谓词 aux_pred，权重被转换为一个有理数对 (Rational(...), Rational(1, 1))。
            # 这里使用了 math.exp(weighting) 来处理规则的权重，并将其转换为分数形式 Fraction，再取出分子和分母分别作为有理数 Rational。这代表了正反两种情况的权重，Rational(1, 1) 可能表示谓词为假的权重是 1。
            weightings[aux_pred] = (Rational(Fraction(math.exp(weighting)).numerator,
                                             Fraction(math.exp(weighting)).denominator), Rational(1, 1))
        for free_var in free_vars:# 遍历每个自由变量 free_var，为它们添加全称量化符号（Universal Quantifier），表示公式中所有的自由变量都要满足该公式。这将公式封装为量化公式 QuantifiedFormula，进一步增强逻辑约束。
            formula = QuantifiedFormula(Universal(free_var), formula)
        sentence = sentence & formula # 使用逻辑与（&）操作符将新的公式与 sentence 结合，累积构建一个包含所有规则的逻辑句子。

    try: # 尝试将句子转换为 SC2 格式。如果失败，抛出异常。
        sentence = to_sc2(sentence)
    except:
        raise ValueError('Sentence must be a valid SC2 formula.')
    return WFOMCProblem(sentence, mln.domain, weightings, mln.cardinality_constraint) # 返回一个 WFOMCSProblem 对象，包含转换后的句子、领域、权重和基数约束。
# 总结：
# 这段代码定义了两个类，WFOMCSProblem 用于描述加权一阶模型计数问题，MLNProblem 用于描述 Markov 逻辑网络问题。MLN_to_WFOMC 函数则实现了将 MLN 问题转换为 WFOMC 问题的逻辑。