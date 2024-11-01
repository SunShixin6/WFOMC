from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import Callable, Iterable
from collections import OrderedDict
from PrettyPrint import PrettyPrintTree

from . import boolean_algebra as backend


__all__ = [
    'Pred',
    'Term',
    'Var',
    'Const',
    'Formula',
    'QFFormula',
    'AtomicFormula',
    'Quantifier',
    'Universal',
    'Existential',
    'Counting',
    'QuantifiedFormula',
    'CompoundFormula',
    'Conjunction',
    'Disjunction',
    'Implication',
    'Equivalence',
    'Negation',
    'BinaryFormula',
    'SCOTT_PREDICATE_PREFIX',
    'AUXILIARY_PRED_NAME',
    'TSEITIN_PRED_NAME',
    'SKOLEM_PRED_NAME',
    'EVIDOM_PRED_NAME',
    'PREDS_FOR_EXISTENTIAL',
    'pretty_print',
    'X', 'Y', 'Z',
    'U', 'V', 'W',
    'top', 'bot',
]


class FOLSyntaxError(Exception):
    pass


class Term(object):
    """
    First-order logic terms, including constants and variables
    """
    name: str


SCOTT_PREDICATE_PREFIX = '@scott'
AUXILIARY_PRED_NAME = '@aux'
TSEITIN_PRED_NAME = '@tseitin'
SKOLEM_PRED_NAME = '@skolem'
EVIDOM_PRED_NAME = '@evidom'
PREDS_FOR_EXISTENTIAL = [
    TSEITIN_PRED_NAME, SKOLEM_PRED_NAME, EVIDOM_PRED_NAME
]


RESERVED_PRED_NAMES: tuple[str] = (
    'true',
    'false',
    SCOTT_PREDICATE_PREFIX,
    AUXILIARY_PRED_NAME,
    TSEITIN_PRED_NAME,
    SKOLEM_PRED_NAME,
    EVIDOM_PRED_NAME
)

RESERVED_VAR_NAMES: tuple[str] = (
)

@dataclass(frozen=True)
class Pred:
    """
    Predicate
    """
    name: str
    arity: int

    def __post_init__(self):
        if self.name.split('_')[0].lower() in RESERVED_PRED_NAMES:
            raise FOLSyntaxError("Predicate name cannot be %s" % self.name)
        if self.arity < 0:
            raise FOLSyntaxError("Arity must be a natural number")

    def __call__(self, *args: Term):
        # NOTE(hack): the callable obj cannot be the column of dataframe
        # if len(args) == 0 or not isinstance(args[0], (Var, Const)):
        #     return self
        if self.arity != len(args):
            raise FOLSyntaxError(
                "Mismatching number of arguments and predicate %s: %s != %s", str(self), self.arity, len(args))
        return AtomicFormula(pred=self, args=tuple(args), positive=True)

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


@dataclass(frozen=True)
class Var(Term):
    """
    Variable
    """
    name: str

    def __post_init__(self):
        if self.name in RESERVED_VAR_NAMES:
            raise FOLSyntaxError("Variable name cannot be %s" % self.name)

    def substitute(self, substitution: dict[Var, Term]) -> Term:
        if self in substitution:
            return substitution[self]
        return self

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)


@dataclass(frozen=True)
class Const(Term):
    """
    Constant
    """
    name: str

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)


class Formula(object):
    """
    Base class for first-order logic formulas
    """
    ...

# 这是一个 Python 的装饰器，来自 dataclasses 模块。@dataclass 用于自动生成类的初始化方法、比较方法等。frozen=True 参数表示这个数据类是不可变的，即类实例一旦创建，它的属性就不能再被修改。这在数学表达式的场景中非常有用，因为数学公式往往不应该在创建后被修改。
@dataclass(frozen=True)
class QFFormula(Formula):
    """
    Quantifier-free formula 定义了一个名为 QFFormula 的类，并且它继承自 Formula 类。QFFormula 用于表示量词自由（quantifier-free）的公式。量词自由公式是逻辑中不含有量词（如存在量词 ∃ 或全称量词 ∀）的表达式。
    """
    expr: backend.Expr # 定义了 QFFormula 类的一个属性 expr，它的类型是 backend.Expr。这里的 backend.Expr 很可能是一个来自于某个外部模块或库的类，表示一个逻辑或数学表达式。

    def __invert__(self) -> QFFormula: # 定义了一个魔术方法 __invert__，它允许对 QFFormula 对象使用逻辑“非”操作符（~）。这个方法返回一个新的 QFFormula 对象。
        return QFFormula(backend.Not(self.expr))

    def __or__(self, other: Formula) -> Formula: # 定义了另一个魔术方法 __or__，用于实现逻辑公式的“或”操作（| 运算符）。它接受另一个 Formula 对象，并返回一个新的逻辑公式。
        if isinstance(other, (Top, Bot)): # 检查 other 是否是 Top 或 Bot 类的实例，Top 和 Bot 很可能是逻辑中的特殊值，分别表示恒真和恒假。
            return other | self # 如果 other 是 Top 或 Bot，则调用 other 的 | 运算符方法，执行逻辑“或”操作。
        if isinstance(other, QFFormula): # 检查 other 是否是 QFFormula 的实例，表示如果 other 也是一个量词自由公式。
            return QFFormula(backend.Or(self.expr, other.expr)) # 对两个量词自由公式的表达式 self.expr 和 other.expr 进行逻辑“或”操作，使用 backend.Or，然后返回一个新的 QFFormula 对象。
        return other | self # 如果 other 既不是 Top 也不是 QFFormula，则调用 other 的 | 运算符方法，将 self 作为参数执行“或”操作。

    def __and__(self, other: Formula) -> Formula: # 定义了 __and__ 方法，用于实现逻辑公式的“与”操作（& 运算符）。它接受另一个 Formula 对象，返回一个新的逻辑公式。
        if isinstance(other, (Top, Bot)):
            return other & self
        if isinstance(other, QFFormula):
            return QFFormula(backend.And(self.expr, other.expr))
        return other & self

    def implies(self, other: Formula) -> Formula:
        if isinstance(other, (Top, Bot)):
            return ~self | other
        if isinstance(other, QFFormula):
            return QFFormula(backend.Implies(self.expr, other.expr))
        return ~self | other # 根据逻辑蕴含的定义：A implies B 等价于 ~A or B。这里通过对 self 取反并与 other 执行逻辑“或”操作来实现。

    def equivalent(self, other: Formula) -> Formula: # 定义了 equivalent 方法，用于实现逻辑公式的“等价”关系，接受一个 Formula 对象 other，返回一个新的公式。
        if isinstance(other, (Top, Bot)):
            return (~self | other) & (self | ~other) # 根据逻辑等价的定义：A 等价于 B 等价于 (~A or B) and (A or ~B)。这里直接返回两个公式的逻辑组合。
        if isinstance(other, QFFormula):
            return QFFormula(backend.Equivalent(self.expr, other.expr)) # 如果 other 是量词自由公式，则通过 backend.Equivalent 对两个表达式 self.expr 和 other.expr 执行逻辑等价操作，并返回一个新的 QFFormula 对象。
        else:
            return other.equivalent(self)

    def __str__(self) -> str:
        return str(self.expr)

    def __repr__(self) -> str:
        return str(self)

    def atoms(self) -> frozenset[AtomicFormula]: # 定义了 atoms 方法，用于返回公式中的原子公式集合，返回类型是不可变的集合 frozenset，其中包含 AtomicFormula 对象。
        return backend.get_atoms(self.expr) # 调用 backend.get_atoms(self.expr)，从当前表达式 self.expr 中提取所有原子公式，并将其作为不可变集合返回。

    def terms(self) -> Iterable[Term]: # 定义了 terms 方法，用于生成公式中的所有术语（Term）。返回类型是可迭代对象 Iterable，其中包含 Term 对象。
        for atom in self.atoms(): # 遍历当前公式中的每个原子公式，通过调用之前定义的 atoms() 方法获得所有原子公式。
            for term in atom.args: # 对每个原子公式的参数 args 进行遍历，args 通常包含原子公式中的所有术语。
                yield term # 逐个返回原子公式中的术语，这使得 terms() 方法成为一个生成器，能够懒加载公式中的所有术语。

    def vars(self) -> frozenset[Var]: # 定义了 vars 方法，用于提取公式中的所有变量，返回类型是一个包含变量对象 Var 的不可变集合 frozenset。
        return frozenset(filter(lambda x: isinstance(x, Var), self.terms())) # 通过调用 terms() 方法提取公式中的所有术语，并使用 filter() 过滤出所有变量 Var，最后将这些变量转换为不可变集合 frozenset 返回。

    def free_vars(self) -> frozenset[Var]: # 定义了 free_vars 方法，可能用于提取公式中的所有自由变量，返回类型也是不可变集合 frozenset，其中包含 Var 对象。
        return self.vars() # 在 free_vars() 方法中，直接调用 vars() 方法返回所有变量。这意味着当前公式中的所有变量都被视为自由变量。

    def consts(self) -> frozenset[Const]: # 定义了 consts 方法，用于提取公式中的所有常量，返回一个包含常量 Const 对象的不可变集合 frozenset。
        return frozenset(filter(lambda x: isinstance(x, Const), self.terms())) # 通过调用 terms() 方法提取公式中的所有术语，并使用 filter() 过滤出所有常量 Const，然后将这些常量转换为不可变集合 frozenset 返回。

    def preds(self) -> frozenset[Pred]: # 定义了 preds 方法，用于提取公式中的所有谓词，返回一个包含谓词 Pred 对象的不可变集合 frozenset。
        return frozenset(atom.pred for atom in self.atoms()) # 通过遍历当前公式中的每个原子公式 atom，提取其中的谓词 pred，并将这些谓词转换为不可变集合 frozenset 返回。

    def satisfiable(self) -> bool:
        return backend.satisfiable(self.expr)

    def models(self) -> Iterable[frozenset[AtomicFormula]]: # 定义了 models 方法，用于生成当前公式的所有模型，返回类型是可迭代对象 Iterable，其中包含不可变集合 frozenset，每个集合中包含若干 AtomicFormula。
        """
        Yield all models of the formula
# 该方法会生成公式的所有模型，并且返回类型是一个可迭代对象 Iterable，其中每个元素是一个包含 Lit 对象（即文字）的不可变集合 frozenset。
        :rtype Iterable[frozenset[Lit]]: models
        """
        if not self.satisfiable():
            raise RuntimeError("Formula is not satisfiable")

        for model in backend.get_models(self.expr): # 调用 backend.get_models(self.expr)，从当前公式的表达式 self.expr 中获取所有模型。model 是一个字典，表示每个符号的布尔值分配。
            yield frozenset( # 将模型的每个符号和其对应的值转换为不可变集合 frozenset 并返回。
                backend.get_atom(symbol) if value else ~backend.get_atom(symbol) # 对于模型中的每个符号 symbol，如果其布尔值 value 为真，调用 backend.get_atom(symbol) 返回该符号对应的原子公式；如果为假，则取反 ~ 操作，返回该符号的否定形式。
                for symbol, value in model.items() # 遍历模型中的每个符号 symbol 和其布尔值 value，并根据其值生成对应的原子公式或其否定形式。
            )

    def substitute(self, substitution: dict[Term, Term]) -> QFFormula:# 定义了 substitute 方法，用于对当前公式中的术语进行替换，接受一个字典 substitution，其中键和值都是 Term 类型，表示术语的替换规则。返回一个新的 QFFormula 对象。
        atom_substitutions = OrderedDict() # 初始化一个 OrderedDict 对象，用于存储每个原子公式及其替换后的表达式，确保替换操作按顺序进行。
        for atom in self.atoms(): # 遍历当前公式中的每个原子公式，调用 atoms() 方法获取所有的原子公式。
            atom_substitutions[atom.expr] = atom.substitute(substitution).expr # 对每个原子公式，调用其 substitute() 方法，根据提供的替换规则 substitution 对原子公式进行替换，并将替换后的表达式存储在 atom_substitutions 字典中。
        return QFFormula(backend.substitute(self.expr, atom_substitutions)) # 使用 backend.substitute() 方法，对整个公式的表达式 self.expr 应用替换操作，并返回一个新的 QFFormula 对象，包含替换后的表达式。

    def sub_nullary_atoms(self, substitution: dict[AtomicFormula, bool]) -> QFFormula:
        substitution = dict((atom.expr, value) for atom, value in substitution.items())
        return QFFormula(backend.substitute(self.expr, substitution))

    def simplify(self) -> QFFormula: # 定义了 simplify 方法，用于简化当前的公式。返回类型是 QFFormula，表示简化后的量词自由公式。
        return QFFormula(backend.simplify(self.expr))


@dataclass(frozen=True)
class AtomicFormula(QFFormula):
    """
    Atomic formula, i.e. a predicate applied to a tuple of terms.
    It is a subclass of QFFormula.
    It actually acts as a literal.
    """
    pred: Pred
    args: tuple[Term]
    positive: bool
    expr: backend.Expr = field(init=False, default=None,
                               hash=False, compare=False)

    def __post_init__(self):
        if len(self.args) != self.pred.arity:
            raise FOLSyntaxError(
                "Number of terms does not match the predicate's arity")
        atom = self
        if not self.positive:
            atom = ~self
        expr = backend.get_symbol(atom)
        expr = expr if self.positive else backend.Not(expr)
        object.__setattr__(self, 'expr', expr)

    @functools.lru_cache(maxsize=None)
    def __invert__(self):
        return AtomicFormula(self.pred, self.args, not self.positive)

    @functools.lru_cache(maxsize=None)
    def make_positive(self):
        if self.positive:
            return self
        return AtomicFormula(self.pred, self.args, True)

    def __str__(self):
        s = '{}({})'.format(self.pred,
                               ','.join([str(arg) for arg in self.args]))
        return s if self.positive else '~' + s

    def __repr__(self):
        return str(self)

    def vars(self) -> frozenset[Var]:
        return frozenset(filter(lambda x: isinstance(x, Var), self.args))

    def consts(self) -> frozenset[Const]:
        return frozenset(filter(lambda x: isinstance(x, Const), self.args))

    def substitute(self, substitution: dict[Term, Term]) -> AtomicFormula:
        substituted_args = []
        for arg in self.args:
            if arg in substitution:
                substituted_args.append(substitution[arg])
            else:
                substituted_args.append(arg)
        return AtomicFormula(self.pred, tuple(substituted_args), self.positive)

    def simplify(self) -> QFFormula:
        return self


@dataclass(frozen=True)
class Bot(QFFormula):
    expr: backend.Expr = field(init=False, default=None)

    def __and__(self, other: Formula) -> Formula:
        return self

    def __or__(self, other: Formula) -> Formula:
        return other

    def __invert__(self) -> QFFormula:
        return top

    def implies(self, other: Formula) -> Formula:
        return top

    def equivalent(self, other: Formula) -> Formula:
        return ~other

    def __str__(self) -> str:
        return '⊥'

    def preds(self) -> frozenset[Pred]:
        return frozenset()


@dataclass(frozen=True)
class Top(QFFormula):
    expr: backend.Expr = field(init=False, default=None)

    def __and__(self, other: Formula) -> Formula:
        return other

    def __or__(self, other: Formula) -> Formula:
        return self

    def __invert__(self) -> QFFormula:
        return bot

    def implies(self, other: Formula) -> Formula:
        return other

    def equivalent(self, other: Formula) -> Formula:
        return other

    def __str__(self) -> str:
        return '⊤'

    def preds(self) -> frozenset[Pred]:
        return frozenset()


@dataclass(frozen=True)
class Quantifier(object):
    quantifier: str = field(init=False)
    quantified_var: Var

    def __str__(self):
        return '{} {}'.format(self.quantifier, self.quantified_var)

    def __repr__(self):
        return str(self)


@dataclass(frozen=True)
class Universal(Quantifier):
    def __post_init__(self):
        object.__setattr__(self, 'quantifier', '\\forall')

    def complement(self) -> Existential:
        return Existential(self.quantified_var)


@dataclass(frozen=True)
class Existential(Quantifier):
    def __post_init__(self):
        object.__setattr__(self, 'quantifier', '\\exists')

    def complement(self) -> Universal:
        return Universal(self.quantified_var)


@dataclass(frozen=True)
class Counting(Quantifier):
    comparator: str
    count_param: int

    def __post_init__(self):
        assert self.comparator in ['='], 'Only equality is supported'
        object.__setattr__(self, 'quantifier', '\\exists')

    def complement(self) -> Counting:
        raise FOLSyntaxError('Complement of counting quantifier is not supported')

    def __str__(self):
        return '{}_{{{}{}}} {}'.format(
            self.quantifier, self.comparator,
            self.count_param, self.quantified_var
        )


class QuantifiedFormula(Formula):
    """
    Quantified formula, e.g. \\forall x P(x),
    \\exists x P(x) and \\exists_{=2} x P(x)
    """
    def __init__(self, quantifier_scope: Quantifier, quantified_formula: Formula):
        self.quantifier_scope = quantifier_scope
        self.quantified_formula = quantified_formula

    def vars(self) -> frozenset[Var]:
        return self.quantified_formula.vars()

    @property
    def quantified_var(self) -> Var:
        return self.quantifier_scope.quantified_var

    def free_vars(self) -> frozenset[Var]:
        return self.quantified_formula.free_vars() - frozenset([self.quantified_var])

    def rename(self, substitution: dict[Term, Term]) -> QuantifiedFormula:
        # filter out the quantified variable
        substitution = {k: v for k, v in substitution.items() if k in self.free_vars()}
        inverse_substitution = {v: k for k, v in substitution.items()}
        if self.quantified_var in inverse_substitution:
            raise FOLSyntaxError('Subsituting variable {} with {} will cause collision in the formula: {}'.format(
                inverse_substitution[self.quantified_var],
                self.quantified_var, self
            )
        )
        quantifier_scope = self.quantifier_scope.rename_quantified_var(
            substitution.get(self.quantified_var, self.quantified_var)
        )
        if isinstance(self.quantified_formula, QFFormula):
            return QuantifiedFormula(quantifier_scope, self.quantified_formula.substitute(substitution))
        elif isinstance(self.quantified_formula, QuantifiedFormula):
            return QuantifiedFormula(quantifier_scope, self.quantified_formula.rename(substitution))
        else:
            raise FOLSyntaxError('Compound quantified formula is not supported')

    def consts(self) -> frozenset[Const]:
        return self.quantified_formula.consts()

    def atoms(self) -> frozenset[AtomicFormula]:
        return self.quantified_formula.atoms()

    def __invert__(self) -> QuantifiedFormula:
        # return QuantifiedFormula(self.quantifier_scope.complement(), ~self.formula)
        if isinstance(self.quantified_formula, QFFormula):
            return QuantifiedFormula(self.quantifier_scope.complement(), ~self.quantified_formula)
        return Negation(self)

    def __or__(self, other: Formula) -> Formula:
        if isinstance(other, QFFormula) and \
                self.quantified_var not in other.vars() and \
                not isinstance(self.quantifier_scope, Counting):
            return QuantifiedFormula(self.quantifier_scope, self.quantified_formula | other)
        return Disjunction(self, other)

    def __and__(self, other: Formula) -> Formula:
        if isinstance(other, QFFormula) and \
                self.quantified_var not in other.vars() and \
                not isinstance(self.quantifier_scope, Counting):
            return QuantifiedFormula(self.quantifier_scope, self.quantified_formula & other)
        return Conjunction(self, other)

    def implies(self, other: Formula) -> Formula:
        return Implication(self, other)

    def equivalent(self, other: Formula) -> Formula:
        return Equivalence(self, other)

    def ext_uni_vars(self) -> tuple[frozenset[Var], frozenset[Var]]:
        all_vars = self.vars()
        if self.exist is None:
            ext_vars = None
        else:
            ext_vars = self.exist.quantified_vars
        return (ext_vars, all_vars - ext_vars)

    def preds(self) -> frozenset[Pred]:
        return self.quantified_formula.preds()

    def is_exist(self) -> bool:
        return self.exist is not None

    def __str__(self):
        return '{}: {}'.format(self.quantifier_scope, str(self.quantified_formula))

    def __repr__(self):
        return str(self)


class CompoundFormula(Formula):
    def __init__(self, op_name: str, op: Callable[..., Formula]):
        self.op_name: str = op_name
        self.op: Callable[..., Formula] = op

    def __invert__(self):
        return Negation(self)

    def __or__(self, other):
        return Disjunction(self, other)

    def __and__(self, other):
        return Conjunction(self, other)

    def implies(self, other):
        return Implication(self, other)

    def equivalent(self, other):
        return Equivalence(self, other)


class Negation(CompoundFormula):
    def __init__(self, formula: Formula) -> None:
        super().__init__('~', lambda x: ~x)
        self.sub_formula: Formula = formula

    def __str__(self):
        return f'~{self.sub_formula}'

    def __repr__(self):
        return str(self)

    def vars(self) -> frozenset[Var]:
        return self.sub_formula.vars()


class BinaryFormula(CompoundFormula):
    def __init__(self, op_name: str, op: Callable[[Formula, Formula], Formula],
                 left: Formula, right: Formula) -> None:
        super().__init__(op_name, op)
        self.left_formula: Formula = left
        self.right_formula: Formula = right

    def __str__(self):
        return f'({self.left_formula}) {self.op_name} ({self.right_formula})'

    def __repr__(self):
        return str(self)

    def vars(self) -> frozenset[Var]:
        return self.left_formula.vars() | self.right_formula.vars()


class Conjunction(BinaryFormula):
    def __init__(self, left: Formula, right: Formula) -> None:
        super().__init__('&', lambda x, y: x & y,
                         left, right)


class Disjunction(BinaryFormula):
    def __init__(self, left: Formula, right: Formula) -> None:
        super().__init__('|', lambda x, y: x | y,
                         left, right)


class Implication(BinaryFormula):
    def __init__(self, left: Formula, right: Formula) -> None:
        super().__init__('->', lambda x, y: x.implies(y),
                         left, right)


class Equivalence(BinaryFormula):
    def __init__(self, left: Formula, right: Formula) -> None:
        super().__init__('<->', lambda x, y: x.equivalent(y),
                         left, right)


def pretty_print(formula: Formula) -> None:
    def get_children(formula: Formula) -> list[Formula]:
        if isinstance(formula, QFFormula):
            return []
        elif isinstance(formula, QuantifiedFormula):
            return [formula.quantified_formula]
        elif isinstance(formula, Negation):
            return [formula.sub_formula]
        elif isinstance(formula, BinaryFormula):
            return [formula.left_formula, formula.right_formula]

    def get_value(formula: Formula) -> str:
        if isinstance(formula, QFFormula):
            return str(formula)
        elif isinstance(formula, QuantifiedFormula):
            return f'{formula.quantifier_scope}:'
        elif isinstance(formula, Negation):
            return '~'
        elif isinstance(formula, BinaryFormula):
            return formula.op_name

    pt = PrettyPrintTree(get_children, get_value, return_instead_of_print=True)
    return pt(formula)


X, Y, Z = Var('X'), Var('Y'), Var('Z')
U, V, W = Var('U'), Var('V'), Var('W')
a, b, c = Const('a'), Const('b'), Const('c')
top, bot = Top(), Bot()
