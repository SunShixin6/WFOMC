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


def dft(context: WFOMCContext,
        formula: QFFormula,
        domain: set[Const],
        get_weight: Callable[[Pred], tuple[RingElement, RingElement]],
        leq_pred: Pred):
    tmp = Rational(0, 1) # 记录每次k 得到的结果
    sum_q = Rational(0, 1) # 总的求和，是对q 求和
    tmp_cache = dict()  # 用于缓存相同的ki_div_Mi 的CCG结果

    n, D, M = context.get_coef() # 获取三个参数
    for k in D:
        dot_res = context.get_dot(n, k, M)  # n 和 k/M 做点乘
        for index in range(len(k)):
            ki_div_Mi = k[index] / M[index]  # 获得每一个ki/Mi, 传入CCG中
            with open("./log.txt", "a") as f:  # 这里记录一下
                print("\n", "ki_div_Mi: ", ki_div_Mi, file=f)
            exp_coef = my_simplify(exp(I * 2 * pi * dot_res))  # 因为e指数有周期性，可以进行化简到0到2pi之间
            new_get_weight = partial(get_weight, ki_div_Mi)  # 这里get_weight 传入第一个参数ki_div_Mi
            if ki_div_Mi not in tmp_cache:  # 没有计算过ki_div_Mi的CCG
                tmp = recursive_wfomc(formula, domain, new_get_weight, leq_pred) * exp_coef  # 每个WFOMC乘以指数项
                tmp_cache[ki_div_Mi] = tmp # 保存这次的WFOMC结果
            else:
                tmp = tmp_cache[ki_div_Mi] # 获取之前保存的WFOMC结果
            sum_q += tmp # 累加到sum里面
    return sum_q / sum(M) # 除以一个分母


def dft_vector():
    pass
