from __future__ import annotations
from logzero import logger
from wfomc.fol.sc2 import SC2
# new_predicate: 用于创建新的谓词。convert_counting_formula: 可能用于将某种公式转换为计数公式。
from wfomc.fol.utils import new_predicate, convert_counting_formula
# 从 sampling_fo2.network.constraint 模块中导入 CardinalityConstraint，这个类可能用于定义某种基于基数约束的逻辑条件。
from wfomc.network.constraint import CardinalityConstraint
# 从 sampling_fo2.fol.syntax 导入所有定义，可能包含一阶逻辑的语法元素，如谓词、函数符号、量词等。
from wfomc.fol.syntax import *
# 导入 WFOMCSProblem，这很可能是一个用于处理加权一阶模型计数问题的类。
from wfomc.problems import WFOMCProblem
from wfomc.fol.syntax import AUXILIARY_PRED_NAME, SKOLEM_PRED_NAME
from wfomc.utils.third_typing import RingElement, Rational

# 定义了 WFOMCContext 类，这是一个用于加权一阶模型计数（WFOMC）算法的上下文类。文档字符串简单说明了该类的作用——为 WFOMC 算法提供上下文环境。
class WFOMCContext(object):
    """
    Context for WFOMC algorithm
    """

    def __init__(self, problem: WFOMCProblem):
        self.domain: set[Const] = problem.domain # 这是一个 set 集合，包含问题中的常量 Const。这些常量定义了域的范围。{'domain0', 'domain1', 'domain2', 'domain3', 'domain4', 'domain5', 'domain6', 'domain7', 'domain8', 'domain9'}
        self.sentence: SC2 = problem.sentence # 这是 SC2 类型的句子，表示问题中涉及的逻辑公式。
        self.weights: dict[Pred, tuple[Rational, Rational]] = problem.weights # self.weights 是一个字典，用于存储谓词（Pred）及其对应的权重对（tuple[Rational, Rational]）。这个权重对通常表示谓词为真的权重和为假的权重。这里将从 problem 对象中获取权重并赋值给 self.weights。
        self.cardinality_constraint: CardinalityConstraint = problem.cardinality_constraint # 将问题中的基数约束 cardinality_constraint 存储在上下文对象中。这个约束通常用于限制某些谓词出现的频率或数量，表示对某个领域中的元素计数的约束条件。
        self.repeat_factor = 1 # 初始化 self.repeat_factor 为 1，表示一个重复因子，可能用于某种重复处理或扩展操作，在后续代码中可能会用到。

        logger.info('sentence: \n%s', self.sentence) # 记录日志，输出 self.sentence，这是存储在上下文中的逻辑公式。
        logger.info('domain: \n%s', self.domain) # 记录日志，输出 self.domain，这是问题的领域集合。
        logger.info('weights:') # 记录日志，输出谓词及其对应的权重信息。通过遍历 self.weights 字典，逐一输出每个谓词 pred 及其对应的权重 w。
        for pred, w in self.weights.items():
            logger.info('%s: %s', pred, w)
        logger.info('cardinality constraint: %s', self.cardinality_constraint) # 记录日志，输出 self.cardinality_constraint，即基数约束信息。

        self.formula: QFFormula # 这行代码定义了 self.formula，其类型为 QFFormula，表明类中使用了该公式对象。
        self._build() # 用类的内部方法 _build()，该方法应该是用来构建公式或模型的内部结构。
        logger.info('Skolemized formula for WFOMC: \n%s', self.formula)
        logger.info('weights for WFOMC: \n%s', self.weights)

    def contain_cardinality_constraint(self) -> bool: # 该方法检查是否存在基数约束（cardinality constraint）。如果 self.cardinality_constraint 不为 None，则表明存在基数约束，返回 True，否则返回 False。
        return self.cardinality_constraint is not None and \
            not self.cardinality_constraint.empty()

    def contain_existential_quantifier(self) -> bool: # 该方法检查公式中是否包含存在量词。通过调用 self.sentence.contain_existential_quantifier() 来确定。
        return self.sentence.contain_existential_quantifier()

    # 该方法用于获取给定谓词 pred 的权重。权重由一个二元组（RingElement, RingElement）表示。
    # 如果权重在 self.weights 中找不到，则返回默认权重 (Rational(1, 1), Rational(1, 1))。
    def get_weight(self, pred: Pred) -> tuple[RingElement, RingElement]:
        default = Rational(1, 1)
        if pred in self.weights: # 这里通过检查 self.weights 来判断是否已经存储了这个谓词 pred 的权重。self.weights 是一个字典或类似的结构，用于保存不同谓词的权重。
            return self.weights[pred] # 如果 pred 存在于 self.weights 中，那么直接返回存储的权重值 self.weights[pred]。该权重值应为一个 tuple[RingElement, RingElement]，代表两个数值，可能表示正反权重或不同条件下的权重。
        return (default, default) # 如果 pred 不存在于 self.weights 中，函数就会返回默认值 (default, default)，也就是 (Rational(1, 1), Rational(1, 1))。# 这表示如果没有特定的权重信息，谓词将被赋予 1:1 的默认权重。

    def decode_result(self, res: RingElement): # 该方法对结果进行解码。如果没有基数约束，则结果 res 除以重复因子 self.repeat_factor。
        if not self.contain_cardinality_constraint():
            return res / self.repeat_factor
        res = self.cardinality_constraint.decode_poly(res)
        return res / self.repeat_factor

    def _skolemize_one_formula(self, formula: QuantifiedFormula) -> QFFormula:
        """
        Only need to deal with \forall X \exists Y: f(X,Y) or \exists X: f(X,Y)
        """
        quantified_formula = formula.quantified_formula
        quantifier_num = 1
        while(not isinstance(quantified_formula, QFFormula)):
            quantified_formula = quantified_formula.quantified_formula
            quantifier_num += 1

        formula: QFFormula = top
        ext_formula = quantified_formula
        if not isinstance(ext_formula, AtomicFormula):
            aux_pred = new_predicate(quantifier_num, AUXILIARY_PRED_NAME)
            aux_atom = aux_pred(X, Y) if quantifier_num == 2 else aux_pred(X)
            formula = formula & (ext_formula.equivalent(aux_atom))
            ext_formula = aux_atom

        if quantifier_num == 2:
            skolem_pred = new_predicate(1, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred(X)
        elif quantifier_num == 1:
            skolem_pred = new_predicate(0, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred()
        formula = formula & (skolem_atom | ~ext_formula)
        self.weights[skolem_pred] = (Rational(1, 1), Rational(-1, 1))
        return formula

    def _skolemize(self) -> QFFormula:
        """
        Skolemize the sentence
        """
        formula = self.sentence.uni_formula
        while(not isinstance(formula, QFFormula)):
            formula = formula.quantified_formula
        for ext_formula in self.sentence.ext_formulas:
            formula = formula & self._skolemize_one_formula(ext_formula)
        return formula

    def _build(self): # 定义了一个 _build 方法，用于构建当前对象的公式结构以及相关的约束和权重。
        self.formula = self.sentence.uni_formula # 将 self.sentence 中的全称量化公式 uni_formula 赋值给 self.formula，用于开始构建逻辑公式。
        while(not isinstance(self.formula, QFFormula)): # 进入一个 while 循环，检查 self.formula 是否是 QFFormula 类型。如果不是，则继续解包公式。
            self.formula = self.formula.quantified_formula # 在循环中，将公式替换为其内部的量化公式部分，直到找到一个具体的量化公式 QFFormula。

        self.ext_formulas = self.sentence.ext_formulas # 将句子中的存在量化公式列表 ext_formulas 赋值给 self.ext_formulas，这些公式将在后续构建过程中使用。
        if self.sentence.contain_counting_quantifier(): # 检查句子中是否包含计数量化公式。如果存在，将进入 if 语句。
            logger.info('translate SC2 to SNF') # 记录日志，说明当前正在将 SC2 公式转换为 Skolem Normal Form（SNF），这是一种用于逻辑公式转换的标准形式。
            if not self.contain_cardinality_constraint():  # 检查当前对象是否包含基数约束条件。如果不存在，则进入 if 语句。
                self.cardinality_constraint = CardinalityConstraint() #  如果基数约束条件不存在，创建一个新的 CardinalityConstraint 对象并赋值给 self.cardinality_constraint。
            for cnt_formula in self.sentence.cnt_formulas: # 遍历句子中的所有计数量化公式 cnt_formulas。
                # 调用 convert_counting_formula 函数，将计数量化公式 cnt_formula 转换为标准形式。该函数返回全称公式 uni_formula，存在公式 ext_formulas，基数约束 cardinality_constraint，以及重复因子 repeat_factor。
                uni_formula, ext_formulas, cardinality_constraint, repeat_factor = \
                    convert_counting_formula(cnt_formula, self.domain)
                self.formula = self.formula & uni_formula # 将转换后的全称量化公式 uni_formula 与现有的 self.formula 进行逻辑与操作，更新公式。
                self.ext_formulas = self.ext_formulas + ext_formulas # 将转换后的存在量化公式 ext_formulas 添加到现有的 self.ext_formulas 列表中。
                self.cardinality_constraint.add_simple_constraint(*cardinality_constraint)
                self.repeat_factor *= repeat_factor

        if self.contain_cardinality_constraint(): # 检查是否包含基数约束（cardinality constraint）
            self.cardinality_constraint.build()# 如果包含基数约束，则调用 self.cardinality_constraint.build() 来构建该约束。# 基数约束通常用于限制集合的大小或满足某些条件的元素数量，可能与数据集合或约束求解相关。

        for ext_formula in self.ext_formulas:# 遍历 self.ext_formulas 列表中的每个扩展公式（ext_formula）
            self.formula = self.formula & self._skolemize_one_formula(ext_formula)# 然后将每个扩展公式经过 _skolemize_one_formula 方法处理后，与当前的 self.formula 进行逻辑与运算（&）。# Skolem化通常用于消除逻辑公式中的存在量词，使公式更便于处理或简化，常见于自动定理证明领域。

        # self.formula = self.formula.simplify()

        if self.contain_cardinality_constraint():# 首先检查是否包含基数约束（cardinality constraint），如果包含
            self.weights.update(# 然后将返回的结果用于更新当前对象的 self.weights 属性
                self.cardinality_constraint.transform_weighting(
                    self.get_weight,
                )
            )
