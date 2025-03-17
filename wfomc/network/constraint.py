from __future__ import annotations

from abc import ABC
from typing import Callable
from logzero import logger
from dataclasses import dataclass

from wfomc.fol.syntax import Pred
from wfomc.utils import Rational
from wfomc.utils.polynomial import coeff_dict, create_vars, Symbol, expand
from wfomc.utils.third_typing import RingElement


class Constraint(ABC):
    pass


@dataclass(frozen=True)
class TreeConstraint(Constraint):
    pred: Pred

    def __str__(self):
        return "Tree({})".format(self.pred)

    def __repr__(self):
        return str(self)


class CardinalityConstraint(Constraint):
    def __init__(self, constraints: list[tuple[dict[Pred, float], str, float]] = None): # 参数 constraints 是可选的，默认值为 None，其类型是一个列表，列表中包含若干个元组。
        self.constraints: list[tuple[dict[Pred, float], str, float]] = constraints # 每个元组包含三部分：dict[Pred, float] 表示谓词 Pred 到浮点数的映射，str 是一个字符串，float 是一个浮点数。
        if self.constraints is None:
            self.constraints = list() # 如果 constraints 参数为 None，则将 self.constraints 初始化为空列表，确保不为空值。

        self.preds: set[Pred] = set() # 初始化 self.preds 为一个空集合，用于存储谓词 Pred。
        if self.constraints is not None: # 如果 self.constraints 不是 None，则遍历每个约束 constraint。
            for constraint in self.constraints:
                self.preds.update(constraint[0].keys()) # constraint[0] 是一个 dict[Pred, float]，获取它的所有键（这些键是谓词 Pred）。使用 self.preds.update() 将这些键加入到 self.preds 集合中，以收集所有涉及的谓词。

        self.gen_vars: list[Symbol] # 声明 self.gen_vars，类型是 list[Symbol]，但未初始化。这意味着稍后应该赋值一个包含 Symbol 类型元素的列表给 self.gen_vars。
        self.var2pred: dict[Symbol, Pred] = dict() # 初始化 self.var2pred 为一个空字典，用于将 Symbol 映射到 Pred。self.var2pred 的类型是 dict[Symbol, Pred]。
        self.validator: str = "" # 初始化 self.validator 为一个空字符串，用于存储验证信息。类型为 str。

    def empty(self) -> bool: # 定义一个名为 empty 的方法，用于检查约束列表是否为空。
        return len(self.constraints) == 0 # 返回布尔值 True 或 False：如果 self.constraints 长度为 0，则返回 True，表示约束为空。

    def transform_weighting(self, get_weight: Callable[[Pred], tuple[Rational, Rational]]) \
            -> dict[Pred, tuple[Rational, Rational]]: # 定义一个名为 transform_weighting 的方法，接受一个参数 get_weight，它是一个调用签名为 Callable[[Pred], tuple[Rational, Rational]] 的函数。这个函数接收一个 Pred 类型的参数并返回一个包含两个 Rational 类型的元组。
        new_weights: dict[Pred, tuple[RingElement, RingElement]] = {} # 初始化一个空字典 new_weights，用于存储新的权重。键是 Pred 类型，值是包含两个 RingElement 类型的元组（假设 Rational 是 RingElement 的子类或相关类型）。
        self.gen_vars = create_vars('x0:{}'.format( # 调用 create_vars 函数，生成一组变量。变量名格式为 'x0:0' 到 'x0:n'，其中 n 是 self.preds 的长度。生成的变量将被赋值给 self.gen_vars。
            len(self.preds))
        )
        for sym, pred in zip(self.gen_vars, self.preds): # 使用 zip 函数将 self.gen_vars 和 self.preds 配对，遍历每一对（sym 是生成的变量，pred 是预测变量）。
            weight = get_weight(pred) # 调用 get_weight 函数，传入当前的 pred，获取相应的权重。返回的 weight 是一个包含两个 Rational 类型的元组。
            #将 pred 作为键，新的权重作为值存入 new_weights 字典。新的权重由两个部分组成：
            # 第一个部分是 weight[0] * sym，将权重的第一个值与生成的变量相乘。
            # 第二个部分是 weight[1]，直接使用权重的第二个值。
            new_weights[pred] = (weight[0] * sym, weight[1])
            self.var2pred[sym] = pred # 将生成的变量 sym 与对应的预测变量 pred 存储在 self.var2pred 字典中，形成变量与预测之间的映射关系。
        return new_weights # 返回构建好的 new_weights 字典，其中包含了新的权重信息。
    # 这个方法的主要功能是根据给定的权重函数和一组预测变量，生成新的权重，并建立变量与预测之间的映射。返回的字典 new_weights 结构化了权重信息，以便在后续处理或计算中使用。

    def decode_poly(self, poly: RingElement) -> RingElement: # 定义一个名为 decode_poly 的方法，接受一个参数 poly，类型为 RingElement。该方法返回一个 RingElement 类型的值。
        poly = expand(poly) # 调用 expand 函数，展开多项式 poly。这一步通常用于将多项式形式转换为标准形式，确保所有项都明确展开。
        coeffs = coeff_dict(poly, self.gen_vars) # 调用 coeff_dict 函数，传入展开后的多项式 poly 和生成的变量 self.gen_vars。该函数返回一个字典 coeffs，其中包含每个变量的幂次和相应的系数。
        # logger.debug('coeffs: %s', list(coeffs))
        res = Rational(0, 1) # 创建一个 Rational 类型的对象 res，初始值为 0。这个对象用于累加有效的多项式系数。
        for degrees, coeff in coeffs: # 遍历 coeffs 字典中的每一对（degrees 和 coeff）。degrees 是一个列表，表示变量的幂次，coeff 是对应的系数。
            if self.valid(degrees): # 调用 self.valid 方法，检查当前的 degrees 是否有效。有效性根据预设的约束条件判断。
                res += coeff # 如果 degrees 有效，将当前的系数 coeff 加到结果 res 上，累积有效的多项式系数。
        return res # 返回累积后的结果 res，即有效的多项式系数的总和。

    def valid(self, degrees: list[int]) -> bool: # 定义一个名为 valid 的方法，接受一个参数 degrees，类型为整数列表。该方法返回一个布尔值，表示给定的 degrees 是否有效。
        # kwargs = zip((self.var2pred[sym].name for sym in self.gen_vars), degrees) # 使用 zip 函数将生成的变量名与其对应的幂次 degrees 配对。这里通过生成器表达式 self.var2pred[sym].name for sym in self.gen_vars 获取每个变量对应的预测名称，将它们与 degrees 中的值组合成元组。
        # return eval(self.validator.format(**dict(kwargs))) # 使用 eval 函数动态执行构建的验证器表达式。首先用 kwargs 字典替换 self.validator 中的占位符，生成完整的表达式，然后评估该表达式的真假并返回结果。这一步会根据之前 build 方法构建的逻辑表达式进行验证。

        # 第一步：创建一个生成器，用于从 `self.gen_vars` 中获取每个符号变量的名称
        # 并与对应的 `degrees` 值配对
        predicate_names = (self.var2pred[sym].name for sym in self.gen_vars) # predicate_names 是一个生成器，它从 self.gen_vars 中获取每个符号变量，并使用 self.var2pred 将符号变量映射到相应的谓词名称。

        # 第二步：将 `predicate_names` 与 `degrees` 配对生成键值对，形成一个字典 `kwargs`
        kwargs = dict(zip(predicate_names, degrees)) # kwargs 是一个字典，用 zip 函数将 predicate_names 和 degrees 组合成键值对。每个键是谓词的名称，值是相应的指数。

        # 第三步：将 `self.validator` 中的占位符替换为实际的 `kwargs` 值
        # `self.validator` 是预定义的字符串模板，包含约束表达式
        formatted_validator = self.validator.format(**kwargs) # formatted_validator 是格式化后的字符串。通过 self.validator.format(**kwargs) 将 self.validator 中的每个占位符替换为 kwargs 中的对应值。

        # 第四步：使用 `eval` 计算替换后的表达式，返回布尔值
        result = eval(formatted_validator) # eval(formatted_validator) 计算并返回 formatted_validator 表达式的布尔值，即是否满足约束条件。

        # 返回最终的验证结果
        return result

    def extend_simple_constraints(self, ccs: list[tuple[Pred, str, int]]):
        for pred, comp, card in ccs:
            self.add_simple_constraint(pred, comp, card)

    def add_simple_constraint(self, pred: Pred, comp: str, card: int):
        """
        Add a constraint of the form |pred| comp card
        """
        self.constraints.append(({pred: 1}, comp, card))
        self.preds.add(pred)

    def add(self, expr: dict[Pred, float], comp: str, param: float):
        self.constraints.append((expr, comp, param))
        self.preds.update(expr.keys())

    def build(self): # 这段代码是一个类中的方法，名为 build，其主要功能是构建一个用于验证约束条件的表达式。
        validator_list: list[str] = [] # 创建一个空列表 validator_list，用于存储生成的验证表达式。类型注解 list[str] 指定这个列表将只包含字符串。
        for expr, comp, param in self.constraints:# expr、comp 和 param 分别代表约束条件的表达式、比较符号和参数值。
            single_validator = [] # 初始化一个空列表 single_validator，用于存储当前约束的单个验证表达式。
            for pred, coef in expr.items(): # 遍历 expr 中的每个项，pred 是约束的一个预测变量（可能是一个对象），coef 是与该预测变量相乘的系数。
                single_validator.append(f'{coef} * {{{pred.name}}}') # 将格式化字符串 '{coef} * {{{pred.name}}}' 添加到 single_validator 列表中，表示用系数乘以预测变量的名称（用大括号 {} 包裹变量名可能是为了在某些模板系统中进行替换）。
            single_validator = ' + '.join(single_validator) # 将 single_validator 列表中的所有字符串用 ' + ' 连接成一个完整的表达式，例如 2 * {x} + 3 * {y}。
            if comp == '=': # 如果比较符号 comp 是 '='，则将其替换为 '=='，因为在 Python 中，== 用于相等比较。
                comp = '=='
            validator_list.append(f'{single_validator} {comp} {param}') # 将当前的验证表达式（例如 2 * {x} + 3 * {y} == 5）添加到 validator_list 列表中。
        self.validator = ' and '.join(validator_list) #将所有的验证表达式用 ' and ' 连接起来，形成一个完整的验证条件，赋值给 self.validator。
        logger.info('cardinality validator: \n%s', self.validator)
        # 这段代码的作用是从约束条件中生成一个综合的验证表达式，方便后续的验证操作。生成的表达式是由多个条件组成的逻辑与（and）表达式，每个条件由预测变量的线性组合与参数进行比较。
    def __str__(self):
        s = ''
        for expr, comp, param in self.constraints:
            s += ' + '.join(f'{coef} * |{pred.name}|' for pred, coef in expr.items())
            s += ' {} {}'.format(comp, param)
            s += '\n'
        return s

    def __repr__(self):
        return str(self)
