from wfomc.algo import recursive_wfomc
from wfomc.context import WFOMCContext
from wfomc.fol.syntax import Const, Pred, QFFormula
from typing import Callable

from wfomc.network import CardinalityConstraint
from wfomc.utils import RingElement, Rational
from itertools import product
import numpy as np
from sympy import Rational, I, exp, symbols, pi  # 新导入的
from wfomc.utils.simplify import my_simplify
from functools import partial
from wfomc.utils.polynomial import expand


def dft(context: WFOMCContext,
        formula: QFFormula,
        domain: set[Const],
        get_weight: Callable[[Pred], tuple[RingElement, RingElement]],
        leq_pred: Pred):
    final_res = Rational(0, 1)  # 总的求和，是对q 求和
    res_cache = dict()  # 用于缓存相同的ki_div_Mi 的CCG结果

    n, D, M = context.get_coef()  # 获取三个参数
    for k in D:
        dot_res = context.get_dot(n, k, M)  # n 和 k/M 做点乘
        for index in range(len(k)):
            ki_div_Mi = k[index] / M[index]  # 获得每一个ki/Mi, 传入CCG中
            exp_coef = my_simplify(exp(I * 2 * pi * dot_res))  # 因为e指数有周期性，可以进行化简到0到2pi之间
            new_get_weight = partial(get_weight, ki_div_Mi)  # 这里get_weight 传入第一个参数ki_div_Mi

            # with open("./log.txt", "a") as f:  # 这里记录一下
            #     print("\n", "ki_div_Mi: ", ki_div_Mi, ",exp_coef: ", exp_coef, file=f)

            if ki_div_Mi not in res_cache:  # 没有计算过ki_div_Mi的CCG, 就在下面算一下
                tmp_res = recursive_wfomc(formula, domain, new_get_weight, leq_pred, exp_coef)  # 每个WFOMC乘以指数项
                res_cache[ki_div_Mi] = tmp_res  # 保存这次的WFOMC结果
            else:  # 如果计算过ki_div_Mi的CCG，从之前保存的WFOMC结果中获得
                tmp_res = res_cache[ki_div_Mi]
            final_res += tmp_res  # 累加到sum里面
    return expand(final_res / sum(M))  # 除以一个分母
