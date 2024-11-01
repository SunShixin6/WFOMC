from __future__ import annotations
# 导入 reduce 函数，它用于对序列中的所有元素进行累计操作，可以将函数连续应用于序列的项，得到单一结果。
from functools import reduce
from logzero import logger
# 解释: 从 typing 模块中导入 Callable 和 Union 类型，用于类型提示。Callable 表示可调用对象（如函数），Union 表示可以是多种类型之一。
from typing import Callable, Union
# 解释: 从 sampling_fo2.fol.utils 模块中导入 new_scott_predicate 函数，可能是用于创建 Scott 范畴谓词的工具函数。
from wfomc.fol.utils import new_predicate, new_scott_predicate
from .syntax import * #  从当前模块的 syntax 文件中导入所有内容，包含了用于公式操作的类和函数。
from .syntax import FOLSyntaxError # 仅从 syntax 文件中导入 FOLSyntaxError，这是一个表示一阶逻辑语法错误的自定义异常类。


class SC2(Formula):# 定义一个名为 SC2 的类，继承自 Formula 类，表示一种特定类型的公式结构。
    def __init__(self, uni_formula: QuantifiedFormula = None,# 定义类的构造函数，初始化 SC2 对象，接收三个参数：uni_formula (全称量化公式)、ext_formulas (存在量化公式的列表) 和 cnt_formulas (计数量化公式的列表)，都可以为空。
                       ext_formulas: list[QuantifiedFormula] = None,
                       cnt_formulas: list[QuantifiedFormula] = None):
        self.uni_formula: QuantifiedFormula = uni_formula # 将 uni_formula 参数赋值给类的实例变量 self.uni_formula，用于存储全称量化公式。
        self.ext_formulas: list[QuantifiedFormula] = ext_formulas or [] # 如果传入的 ext_formulas 为空，则将 ext_formulas 设置为一个空列表，否则使用传入的值。
        self.cnt_formulas: list[QuantifiedFormula] = cnt_formulas or [] # 同样的逻辑应用于 cnt_formulas，如果为空，则初始化为空列表。
        self.index = 0 # 初始化 index 变量为 0，用于在遍历公式或操作时记录索引。

    def contain_existential_quantifier(self) -> bool:# 定义一个方法 contain_existential_quantifier，返回值为布尔类型，用于检查是否包含存在量化公式。
        return len(self.ext_formulas) > 0 # 如果 ext_formulas 列表的长度大于 0，则返回 True，表示包含存在量化公式；否则返回 False。

    def contain_counting_quantifier(self) -> bool: # 定义一个方法 contain_counting_quantifier，返回值为布尔类型，用于检查是否包含计数量化公式。
        return len(self.cnt_formulas) > 0 # 如果 cnt_formulas 列表的长度大于 0，则返回 True，表示包含计数量化公式；否则返回 False。

    def append_ext(self, formula: QuantifiedFormula): # 定义方法 append_ext，用于向 ext_formulas 列表中添加一个存在量化公式。
        self.ext_formulas.append(formula) # 将传入的 formula 添加到 ext_formulas 列表中。

    def append_cnt(self, formula: QuantifiedFormula): # 定义方法 append_cnt，用于向 cnt_formulas 列表中添加一个计数量化公式。
        self.cnt_formulas.append(formula) # 将传入的 formula 添加到 cnt_formulas 列表中。

    def preds(self): # 定义方法 preds，用于返回公式中包含的所有谓词。
        p = self.uni_formula.preds() if self.uni_formula is not None else set() # 如果 uni_formula 不为空，则调用其 preds 方法获取所有谓词，否则返回一个空集合 set()。
        return reduce(lambda x, y: x.union(y), map(lambda f: f.preds(), self.ext_formulas), p) # 使用 reduce 函数合并 ext_formulas 中所有公式的谓词集合，最终返回包含所有谓词的集合。

    def pred_by_name(self, name: str) -> Union[Pred, None]: # 定义方法 pred_by_name，根据谓词名 name 查找并返回与该名字匹配的谓词。
        for pred in self.preds(): # 遍历所有谓词。
            if pred.name == name: # 如果谓词的名称与传入的 name 匹配。
                return pred # 返回匹配的谓词对象。
        return None

    def __str__(self) -> str:
        s = ''
        if self.uni_formula is not None:
            s += f'Universally quantified formula: {self.uni_formula}'
        if self.ext_formulas:
            s += '\nExistentially quantified formulas:\n'
            s += '\n'.join(map(str, self.ext_formulas))
        if self.cnt_formulas:
            s += '\nCounting quantified formulas:\n'
            s += '\n'.join(map(str, self.cnt_formulas))
        return s

    def __repr__(self) -> str:
        return str(self)


def dfs(formula: Formula, f: Callable[[Formula], Formula]) -> Formula:
    if isinstance(formula, QFFormula):
        formula, additional_formulas = f(formula)
        return formula, additional_formulas
    elif isinstance(formula, QuantifiedFormula):
        quantified_formula, additional_formulas = dfs(formula.quantified_formula, f)
        formula.quantified_formula = quantified_formula
        formula, a_formulas = f(formula)
        additional_formulas.extend(a_formulas)
        return formula, additional_formulas
    else:
        if isinstance(formula, Negation):
            sub_formula, additional_formulas = dfs(formula.sub_formula, f)
            formula = ~sub_formula
            formula, a_formulas = f(formula)
            additional_formulas.extend(a_formulas)
            return formula, additional_formulas
        else:
            left_formula, additional_formulas = dfs(formula.left_formula, f)
            right_formula, a_formulas = dfs(formula.right_formula, f)
            formula = formula.op(left_formula, right_formula)
            additional_formulas.extend(a_formulas)
            formula, a_formulas = f(formula)
            additional_formulas.extend(a_formulas)
            return formula, additional_formulas


def bfs(formula: Formula, f: Callable[[Formula], Formula]) -> Formula:
    formula, additional_formulas = f(formula)
    if isinstance(formula, QFFormula):
        return formula, additional_formulas
    elif isinstance(formula, QuantifiedFormula):
        formula.quantified_formula, a_formulas = bfs(formula.quantified_formula, f)
        additional_formulas.extend(a_formulas)
        return formula, additional_formulas
    else:
        if isinstance(formula, Negation):
            sub_formula, a_formulas = bfs(formula.sub_formula, f)
            formula = ~sub_formula
            additional_formulas.extend(a_formulas)
            return formula, additional_formulas
        else:
            left_formula, a_formulas = bfs(formula.left_formula, f)
            additional_formulas.extend(a_formulas)
            right_formula, a_formulas = bfs(formula.right_formula, f)
            formula = formula.op(left_formula, right_formula)
            additional_formulas.extend(a_formulas)
            return formula, additional_formulas


def transformer(name: str, formula_clses: tuple[type[Formula]]):
    def decorator(f):
        def wrapper(formula: Formula):
            if isinstance(formula, formula_clses):
                logger.debug("%s: %s", name, formula)
                formula = f(formula)
                if isinstance(formula, tuple):
                    return formula
                return formula, []
            return formula, []
        return wrapper
    return decorator


@transformer("Convert implication and equivalence", (Implication, Equivalence))
def convert_implies_equiv(formula: Formula) -> Formula:
    """
    After this transformation, the formula should not contain any implication or equivalence
    """
    if isinstance(formula, Implication):
        return ~formula.left_formula | formula.right_formula
    elif isinstance(formula, Equivalence):
        return (~formula.left_formula | formula.right_formula) & \
            (~formula.right_formula | formula.left_formula)


@transformer("Push negation", (Negation,))
def push_negation(formula: Formula) -> Formula:
    """
    After this transformation, the negation only appears in quantifer-free formulas
    """
    sub_formula = formula.sub_formula
    if isinstance(sub_formula, Negation):
        # NOTE: in case of multiple negations, e.g., ~~~A(x), we remove them in one go
        return push_negation(~sub_formula.sub_formula)
    elif isinstance(sub_formula, BinaryFormula):
        if isinstance(sub_formula, Conjunction):
            return ~sub_formula.left_formula | ~sub_formula.right_formula
        elif isinstance(sub_formula, Disjunction):
            return ~sub_formula.left_formula & ~sub_formula.right_formula
    elif isinstance(sub_formula, QuantifiedFormula):
        return QuantifiedFormula(
            sub_formula.quantifier_scope.complement(),
            ~sub_formula.quantified_formula
        )


@transformer("Push quantifier-free formula", (BinaryFormula,))
def push_qfformula(formula: Formula) -> Formula:
    """
    After this transformation, the quantifier-free formula only appears in quantified formulas
    That is, the resulting formula only contains quantified formulas and compound formulas
    """
    left = formula.left_formula
    right = formula.right_formula
    if isinstance(right, QFFormula):
        left, right = right, left
    if not isinstance(left, QFFormula):
        return formula
    if isinstance(right, BinaryFormula):
        _left = right.left_formula
        _right = right.right_formula
        if isinstance(right, type(formula)):
            #       &             |
            #      / \    or     / \
            #     F   &         F   |
            if isinstance(_left, BinaryFormula):
                return formula.op(_left, formula.op(left, _right))
            else:
                return formula.op(_right, formula.op(left, _left))
        else:
            # distribute quantifier over conjunction/disjunction
            #       &             |
            #      / \    or     / \
            #     F   |         F   &
            return right.op(
                formula.op(left, _left),
                formula.op(left, _right)
            )
    elif isinstance(right, QuantifiedFormula):
        #       &             |
        #      / \    or     / \
        #     F  Qx:        F  Qx:
        if right.quantified_var not in left.vars():
            return QuantifiedFormula(
                right.quantifier_scope,
                formula.op(right.quantified_formula, left)
            )
        else:
            raise FOLSyntaxError(
                'Not support variable renaming yet, '
                'please rename the variable {} in {}'.format(
                    right.quantified_var, left)
            )
    else:
        raise FOLSyntaxError(
            'Should not reach here'
        )


@transformer("Pop quantifier", (Conjunction, Disjunction))
def pop_quantifier_once(formula: Formula) -> Formula:
    left = formula.left_formula
    right = formula.right_formula
    if isinstance(right, QuantifiedFormula):
        left, right = right, left
    # assert isinstance(left, QuantifiedFormula)
    if not isinstance(left, QuantifiedFormula):
        return formula
    if isinstance(right, QuantifiedFormula):
        if left.quantifier_scope == right.quantifier_scope and \
                isinstance(right.quantifier_scope, Universal):
            return QuantifiedFormula(
                left.quantifier_scope,
                formula.op(left.quantified_formula, right.quantified_formula)
            )
    elif isinstance(right, QFFormula) and not isinstance(left.quantifier_scope, Counting):
        if left.quantified_var not in right.vars():
            return QuantifiedFormula(
                left.quantifier_scope,
                formula.op(left.quantified_formula, right)
            )
        else:
            raise FOLSyntaxError(
                'Not support variable renaming yet, '
                'please rename the variable {} in {}'.format(
                    left.quantified_var, right)
            )

    return formula


@transformer("Distribute quantifier", (QuantifiedFormula, ))
def distribute_quantifier(formula: Formula) -> Formula:
    """
    After this transformation, the quantified formula should not contain any compound formula
    """
    quantified_formula = formula.quantified_formula
    if isinstance(quantified_formula, CompoundFormula):
        if (
            isinstance(quantified_formula, Conjunction) and
            isinstance(formula.quantifier_scope, Universal)
        ) or (
            isinstance(quantified_formula, Disjunction) and
            isinstance(formula.quantifier_scope, Existential)
        ):
            return quantified_formula.op(
                QuantifiedFormula(
                    formula.quantifier_scope,
                    quantified_formula.left_formula
                ),
                QuantifiedFormula(
                    formula.quantifier_scope,
                    quantified_formula.right_formula
                )
            )
        else:
            raise FOLSyntaxError(
                'Not support quantifier distribution for {}'.format(
                    quantified_formula)
            )
    return formula


@transformer("Replace disjunction", (Disjunction, ))
def replace_disjunction(formula: Formula) -> Formula:
    left = formula.left_formula
    right = formula.right_formula
    additional_formulas = []
    if isinstance(left, QuantifiedFormula):
        vars = left.free_vars()
        aux_pred = new_predicate(len(vars), AUXILIARY_PRED_NAME)
        aux_atom = aux_pred(*vars)
        additional_formulas.append(left.equivalent(aux_atom))
        left = aux_atom
    if isinstance(right, QuantifiedFormula):
        vars = right.free_vars()
        aux_pred = new_predicate(len(vars), AUXILIARY_PRED_NAME)
        aux_atom = aux_pred(*vars)
        additional_formulas.append(right.equivalent(aux_atom))
        right = aux_atom
    return left | right, additional_formulas


@transformer('Remove existential quantifier', (QuantifiedFormula,))
def remove_existential_quantifier(formula: Formula) -> Formula:
    # NOTE: only remove existential quantifier, we don't want
    # to complicate the resulting formulas
    assert not isinstance(formula.quantified_formula, CompoundFormula), \
        'Compound formula should be distributed already'
    # NOTE: here we only support FO2, i.e., quantified formulas
    # with two recursive structures
    additional_formulas = []
    if isinstance(formula.quantifier_scope, Existential):
        quantified_formula = formula.quantified_formula
        if isinstance(quantified_formula, QFFormula):
            return formula
            # free_vars = formula.free_vars()
            # aux_pred = new_scott_predicate(len(free_vars))
            # aux_atom = aux_pred(*free_vars)
            # additional_formula = formula.equivalent(aux_atom)
            # for var in free_vars:
            #     additional_formula = QuantifiedFormula(
            #         Universal(var), additional_formula
            #     )
            # additional_formulas.append(additional_formula)
            # return aux_atom, additional_formulas
        elif isinstance(quantified_formula.quantifier_scope, (Universal, Existential)):
            # \exists X \forall Y: f(X,Y) ==>
            # \exists X: S(X) & \forall X\forall Y: (S(X) <-> f(X,Y))
            quantified_var = formula.quantified_var
            aux_pred = new_scott_predicate(1)
            aux_atom = aux_pred(quantified_var)
            additional_formula = quantified_formula.equivalent(aux_atom)
            additional_formula = QuantifiedFormula(
                Universal(quantified_var), additional_formula
            )
            additional_formulas.append(additional_formula)
            formula.quantified_formula = aux_atom
            # formula, a_formulas = remove_quantifier(formula)
            # additional_formulas.extend(a_formulas)
            return formula, additional_formulas
        else:
            raise FOLSyntaxError(
                'Not support quantified formula with more than two quantifiers'
            )
    return formula, additional_formulas


def rename_variables(formula: Formula, depth: int,
                     default_vars: list[Var],
                     substitution: dict[Var, Term]) -> Formula:
    if isinstance(formula, QFFormula):
        return formula.substitute(substitution)
    elif isinstance(formula, QuantifiedFormula):
        quantified_var = formula.quantified_var
        if isinstance(formula.quantifier_scope, Counting):
            quantifier_scope = type(formula.quantifier_scope)(
                default_vars[depth],
                formula.quantifier_scope.comparator,
                formula.quantifier_scope.count_param
            )
        else:
            quantifier_scope = type(formula.quantifier_scope)(
                default_vars[depth]
            )
        substitution[quantified_var] = default_vars[depth]
        quantified_formula = rename_variables(
            formula.quantified_formula,
            depth+1,
            default_vars,
            substitution
        )
        return QuantifiedFormula(quantifier_scope, quantified_formula)
    else:
        if isinstance(formula, Negation):
            return ~rename_variables(formula.sub_formula, depth, default_vars, substitution)
        elif isinstance(formula, BinaryFormula):
            return formula.op(
                rename_variables(formula.left_formula, depth, default_vars, substitution),
                rename_variables(formula.right_formula, depth, default_vars, substitution)
            )


def pop_quantifier(formula: Formula) -> Formula:
    # poping quantifier twice is enough for FO2
    formula, _ = dfs(formula, pop_quantifier_once)
    formula, _ = dfs(formula, pop_quantifier_once)
    return formula


@transformer("Check all conjunction", (BinaryFormula, ))
def check_all_conjunction(formula: Formula) -> bool:
    if not isinstance(formula, Conjunction):
        raise FOLSyntaxError(
            f'Found {type(formula).__name__} formula: {formula}'
        )
    return formula


def standardize(formula: Formula) -> Formula:
    """
    Standardize the given formula to a compound of quantified formulas
    """
    formula, _ = dfs(formula, convert_implies_equiv)
    logger.debug("After convert implication and equivalence: %s", formula)
    formula, _ = bfs(formula, push_negation)
    logger.debug("After push negation: %s", formula)
    formula, _ = bfs(formula, distribute_quantifier)
    logger.debug("After distribute quantifier: %s", formula)
    # formula, _ = bfs(formula, push_qfformula)
    # logger.debug("After push quantified-free formula: %s", formula.pretty())
    # formula, _ = bfs(formula, pop_quantifier)
    # logger.debug("After pop quantifier: %s", formula.pretty())
    return formula


def to_sc2(formula: Formula) -> SC2:
    """
    The formula must satisify that a compound formula has at most one quantifier-free subformula
    """
    if isinstance(formula, QFFormula):
        raise FOLSyntaxError(
            f'Found quantified-free formula: {formula}'
        )
    # elif isinstance(formula, QuantifiedFormula):
    #     if isinstance(formula.quantifier_scope, Universal):
    #         return SNF(uni_formula=formula)
    #     elif isinstance(formula.quantifier_scope, Existential):
    #         return SNF(ext_formulas=[formula])
    #     else:
    #         raise FOLSyntaxError(
    #             f'Found counting quantifier {formula.quantifier_scope}'
    #         )

    logger.debug("Before standardize: %s", formula)
    formula = standardize(formula)
    logger.debug("After standardize: \n%s", pretty_print(formula))
    formula, additional_formulas = dfs(formula, replace_disjunction)
    for additional_formula in additional_formulas:
        formula = formula & additional_formula
    formula = standardize(formula)
    logger.debug("After replace disjunction: \n%s", pretty_print(formula))
    formula, additional_formulas = dfs(formula, remove_existential_quantifier)
    scott_formula = formula # top
    for additional_formula in additional_formulas:
        scott_formula = scott_formula & additional_formula
    logger.debug("After remove existential quantifier: \n%s", pretty_print(scott_formula))
    formula = standardize(scott_formula)
    logger.debug("After standardize: \n%s", pretty_print(formula))
    formula = rename_variables(
        formula, 0, [U, V, W], {}
    )
    formula = rename_variables(
        formula, 0, [X, Y, Z], {}
    )
    logger.debug("After rename variables: \n%s", pretty_print(formula))
    # TODO: disable due to https://github.com/lucienwang1009/lifted_sampling_fo2/issues/8
    # formula = pop_quantifier(formula)
    # logger.debug("After pop quantifier: \n%s", pretty_print(formula))
    # here, the formula must only contain conjunctions
    bfs(formula, check_all_conjunction)

    sc2 = SC2()
    uni_formulas: list[QFFormula] = []
    uni_quantifier_scopes: list[Universal] = []
    def collect_formula(formula: Formula,
                        quantifier_scopes: tuple[Quantifier] = ()):
        if isinstance(formula, QuantifiedFormula):
            quantifier_scopes += (formula.quantifier_scope, )
            collect_formula(formula.quantified_formula, quantifier_scopes)
        elif isinstance(formula, QFFormula):
            if all(isinstance(quantifier_scope, Universal) for quantifier_scope in quantifier_scopes):
                uni_formulas.append(formula)
                uni_quantifier_scopes.append(quantifier_scopes)
            elif any(isinstance(quantifier_scope, Counting) for quantifier_scope in quantifier_scopes):
                collected_formula = reduce(
                    lambda x, y: QuantifiedFormula(y, x),
                    quantifier_scopes[::-1],
                    formula
                )
                if len(quantifier_scopes) == 2 and \
                    isinstance(quantifier_scopes[0], Universal) and \
                        isinstance(quantifier_scopes[1], Counting):
                            sc2.append_cnt(collected_formula)
                else:
                    raise FOLSyntaxError(f"Not support fomula \"{collected_formula}\"")
            else:
                collected_formula = reduce(
                    lambda x, y: QuantifiedFormula(y, x),
                    quantifier_scopes[::-1],
                    formula
                )
                sc2.append_ext(collected_formula)
        elif isinstance(formula, Conjunction):
            collect_formula(formula.left_formula, quantifier_scopes)
            collect_formula(formula.right_formula, quantifier_scopes)
        else:
            raise FOLSyntaxError(
                f'Found {type(formula).__name__} formula: {formula}'
            )
    collect_formula(formula)
    uni_formula = top
    for formula in uni_formulas:
        uni_formula &= formula
    if uni_formula == top:
        sc2.uni_formula = top
    else:
        max_uni_quantifier_scope = max(uni_quantifier_scopes, key=len)
        sc2.uni_formula = reduce(
            lambda x, y: QuantifiedFormula(y, x),
            max_uni_quantifier_scope[::-1],
            uni_formula
        )
    return sc2
