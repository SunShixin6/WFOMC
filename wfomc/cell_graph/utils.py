from __future__ import annotations

from wfomc.fol.syntax import AtomicFormula
from wfomc.utils.third_typing import RingElement


def conditional_on(models: dict[frozenset[AtomicFormula], RingElement],
                   gnd_lits: frozenset[AtomicFormula],
                   evidence: frozenset[AtomicFormula] = None) \
        -> dict[frozenset[AtomicFormula], RingElement]: #
    """
    :param models: a dictionary of models and their weights # models: 一个字典，键是不可变集合（frozenset[AtomicFormula]），代表一组公式（即模型）；值是该模型的权重（RingElement 类型）。
    :param gnd_lits: all possible ground literals in the formula # gnd_lits: 一个不可变集合，包含所有可能的基本文字（ground literals）。
    :param evidence: a set of ground literals # evidence: 一个不可变集合，表示证据。如果没有提供证据，默认为 None。
    :return: a dictionary of models and their weights conditioned on the evidence
    """
    if evidence is None: # 如果 evidence 为 None，直接返回原始的 models 字典，因为没有证据可用于过滤。
        return models
    # 如果 evidence 不为空，则过滤 models 中的键，使得这些键包含证据中所有的文字
    filtered_models = dict(filter( # filter 函数对 models 的每个键值对进行检查。
        lambda x: evidence.intersection(gnd_lits).issubset(x[0]), # evidence.intersection(gnd_lits).issubset(x[0]) 用于检查证据中的基本文字是否是当前模型（x[0]）的子集。如果是，则保留该模型。
        models.items()
    ))
    return filtered_models
