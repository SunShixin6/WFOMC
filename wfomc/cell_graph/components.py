from __future__ import annotations

import functools

from functools import reduce
from typing import FrozenSet, List, Tuple
from dataclasses import dataclass, field
from logzero import logger
from sympy import Poly
from wfomc.cell_graph.utils import conditional_on

from wfomc.fol.syntax import AtomicFormula, Pred, Term, a, b, X
from wfomc.fol.utils import get_predicates
from wfomc.utils import Rational
from wfomc.utils.third_typing import RingElement


@dataclass(frozen=True)
class Cell(object):
    """
    In other words, the Unary types
    """
    code: Tuple[bool] = field(hash=False, compare=False)
    preds: Tuple[Pred] = field(hash=False, compare=False)
    # for hashing
    _identifier: FrozenSet[Tuple[Pred, bool]] = field(
        default=None, repr=False, init=False, hash=True, compare=True)

    def __post_init__(self):
        object.__setattr__(self, '_identifier',
                           frozenset(zip(self.preds, self.code)))

    @functools.lru_cache(maxsize=None)
    def get_evidences(self, term: Term) -> FrozenSet[AtomicFormula]:
        evidences: set[AtomicFormula] = set()
        for i, p in enumerate(self.preds):
            atom = p(*([term] * p.arity))
            if (self.code[i]):
                evidences.add(atom)
            else:
                evidences.add(~atom)
        return frozenset(evidences)

    @functools.lru_cache(maxsize=None)
    def is_positive(self, pred: Pred) -> bool:
        return self.code[self.preds.index(pred)]

    def negate(self, pred: Pred) -> Cell:
        idx = self.preds.index(pred)
        new_code = list(self.code)
        new_code[idx] = not new_code[idx]
        return Cell(tuple(new_code), self.preds)

    def drop_preds(self, preds: List[Pred] = None, prefixes: List[str] = None) -> Cell:
        if not preds and not prefixes:
            raise RuntimeError(
                'Dropped pred is not assigned'
            )
        if preds is not None:
            all_preds = [pred for pred in preds]
        else:
            all_preds = []

        if prefixes is not None:
            for prefix in prefixes:
                all_preds.extend(
                    get_predicates(prefix)
                )
        new_code, new_preds = zip(
            *[(c, p) for c, p in zip(self.code, self.preds)
              if p not in all_preds]
        )
        return Cell(tuple(new_code), tuple(new_preds))

    def __str__(self):
        evidences: frozenset[AtomicFormula] = self.get_evidences(X)
        # evidences = filter(lambda x: x.pred.name.startswith('skolem') or x.pred.name.startswith('aux') or x.pred.name == 'E', evidences)
        lits = [str(lit) for lit in evidences]
        lits.sort()
        return '^'.join(lits)

    def __repr__(self):
        return self.__str__()


class TwoTable(object):
    def __init__(self, models: dict[frozenset[AtomicFormula], RingElement], # models：一个字典，键为 frozenset[AtomicFormula] 类型，值为 RingElement 类型，表示一些模型及其相关权重。
                 gnd_lits: frozenset[AtomicFormula]): # gnd_lits：一个 frozenset（不可变集合），其中的元素为 AtomicFormula，表示一组基础的文字（ground literals）。
        self.models = models # 将传入的 models 参数赋值给实例变量 self.models，用于存储模型数据。
        self.gnd_lits = gnd_lits # 将传入的 gnd_lits 参数赋值给实例变量 self.gnd_lits，用于存储基础文字。

    # 定义方法 get_weight，用于计算给定证据下的模型权重。此方法接受一个可选参数 evidence，类型为 FrozenSet[AtomicFormula]，默认为 None。返回值类型是 Poly，表示一个多项式。
    def get_weight(self, evidence: FrozenSet[AtomicFormula] = None) -> Poly:
        if not self.satisfiable(evidence):# 检查当前 evidence 是否使模型可满足，调用 satisfiable 方法来判断。如果 evidence 使模型不可满足，则进入条件语句内部。
            return Rational(0, 1) # 如果 evidence 使模型不可满足，则返回分数值 Rational(0, 1)（即 0）。这意味着在给定证据下，模型的权重为零。
        conditional_models = conditional_on(self.models, self.gnd_lits, evidence) # 调用 conditional_on 函数（假设已定义），使用 self.models、self.gnd_lits 和 evidence 作为参数，得到一个在给定证据下的条件模型 conditional_models。
        ret = reduce( # 使用 reduce 函数将 conditional_models 中所有权重相加，初始值为 Rational(0, 1)。ret 表示所有符合条件的模型的总权重。
            lambda a, b: a + b,
            conditional_models.values(),
            Rational(0, 1)
        )
        return ret # 返回 ret 作为 get_weight 方法的结果，即给定证据下的总权重。

    # 定义方法 get_two_tables，用于获取符合证据的模型数据表。该方法接受一个可选参数 evidence，类型为 FrozenSet[AtomicFormula]，默认为 None。返回类型为 Tuple[FrozenSet[AtomicFormula], Poly]，表示满足条件的模型及其对应权重的元组。
    def get_two_tables(self, evidence: FrozenSet[AtomicFormula] = None) \
            -> Tuple[FrozenSet[AtomicFormula], Poly]:
        if not self.satisfiable(evidence): # 检查给定证据 evidence 是否使模型可满足，调用 satisfiable 方法。如果不可满足，则进入条件语句内部。
            return tuple() # 若模型在给定证据下不可满足，则返回一个空的元组 tuple()，表示无满足条件的模型。
        conditional_models = conditional_on(self.models, self.gnd_lits, evidence) # 调用 conditional_on 函数，根据 self.models、self.gnd_lits 和 evidence 生成符合条件的模型 conditional_models。
        return tuple(conditional_models.items()) # 将 conditional_models 中的项目（模型及其权重）转换为元组形式，并作为 get_two_tables 的返回值。

    # 定义 satisfiable 方法，用于判断在给定证据下是否有满足条件的模型。接受一个可选参数 evidence，类型为 FrozenSet[AtomicFormula]，默认为 None，返回一个布尔值。
    def satisfiable(self, evidence: FrozenSet[AtomicFormula] = None) -> bool:
        # 调用 conditional_on 函数，根据 self.models、self.gnd_lits 和 evidence 生成符合条件的模型 conditional_models。
        conditional_models = conditional_on(self.models, self.gnd_lits, evidence)
        if len(conditional_models) == 0: # 检查 conditional_models 的长度。如果没有符合条件的模型（即长度为 0），则进入条件语句内部。
            return False # 如果 conditional_models 为空，则返回 False，表示模型在给定证据下不可满足。
        return True # 若 conditional_models 非空，则返回 True，表示模型在给定证据下是可满足的。
