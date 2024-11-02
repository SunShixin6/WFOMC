from wfomc.algo import recursive_wfomc
from wfomc.fol.syntax import Const, Pred, QFFormula
from typing import Callable
from wfomc.utils import RingElement, Rational


# DFT:求g（k），对于任意k属于D
def get_dft(
        formula: QFFormula,  # # 要求解的逻辑公式。skolem之后的
        domain: set[Const],
        get_weight: Callable[[Pred], tuple[RingElement, RingElement]],  # 一个正 一个负
        leq_pred: Pred,  # leq_pred: 谓词，用于指定“≤”关系。线性阶，这个本质上是一个binary predicate
        real_version: bool = True) -> RingElement:
    # 调用WFOMC，只是让里面的权重增加一项而已。相当于在原始WFOMC外面包了一层
    res = recursive_wfomc(formula, domain, get_weight, leq_pred, True)
    return res
    pass


# 傅里叶反变换
def get_reverse_DFT(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    # 先调用上面的变换，然后在这个函数里面完成反变换
    res = get_dft(formula, domain, get_weight, leq_pred, real_version)
    # 反变换



    return res
    pass


# 主函数
def dft(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred, real_version: bool = True):
    res =get_reverse_DFT(formula, domain, get_weight, leq_pred, real_version)
    return  res
    pass



