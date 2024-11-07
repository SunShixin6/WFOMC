from wfomc.algo import recursive_wfomc
from wfomc.context import WFOMCContext
from wfomc.fol.syntax import Const, Pred, QFFormula
from typing import Callable

from wfomc.network import CardinalityConstraint
from wfomc.utils import RingElement, Rational
from itertools import product
import numpy as np
from sympy import Rational, I, exp, symbols, pi # 新导入的
from wfomc.utils.simplify import my_simplify

def g(formula: QFFormula, domain: set[Const], get_weight: Callable[[Pred], tuple[RingElement, RingElement]], leq_pred: Pred, ki_div_Mi, real_version: bool = True) -> RingElement:
    # 调用WFOMC，只是让里面的权重增加一项而已。相当于在原始WFOMC外面包了一层
    res = recursive_wfomc(formula, domain, get_weight, leq_pred, real_version, ki_div_Mi, ues_dft = True ) # ki_div_Mi, ues_dft = True 这两个参数是新加的
    return res


def get_dot(n: list[Rational], k: tuple[Rational, Rational], M: list[Rational]):
    tmp = Rational(0, 1) # 确保每次计算都使用 Rational 类型
    for i in range(len(n)):
        tmp += n[i] * k[i] / M[i]
    return tmp


# 傅里叶反变换得到q
def q(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    tmp = Rational(0,1)
    sum_q = Rational(0,1 )
    tmp_cache = dict() # 用于缓存相同的ki_div_Mi 的CCG结果
    for k in D:
        dot_res = get_dot(n, k, M)  # n 和 k/M 做点乘 # 这里不能是numpy类型

        for index in range(length):
            ki_div_Mi = k[index] / M[index] # 获得每一个ki/Mi, 传入CCG中
            exp_ = my_simplify(exp(I * 2 * pi * dot_res)) # 因为e指数有周期性，可以进行化简到0到2pi之间
            if ki_div_Mi not in tmp_cache: # 没有计算过ki_div_Mi的CCG
                tmp = g(formula, domain, get_weight, leq_pred, ki_div_Mi, real_version)  * exp_# 先调用上面的变换，然后在这个函数里面完成反变换 #
                tmp_cache[ki_div_Mi] = tmp
            else:
                tmp = tmp_cache[ki_div_Mi] # 计算过ki_div_Mi的CCG
            sum_q += tmp
    return sum_q / sum(M)


def generate_D(domain_size, var_counts): # var_counts 是一个列表，其中包含每个公式的变量数量
    """
    生成集合 D，其中 D 是一个多维整数向量的集合。
    """
    ranges = [range(int(Rational(domain_size) ** Rational(var_count) + 1)) for var_count in var_counts]  # ranges 中的每个 range 表示每个公式可以生成的计数范围。
    D = [tuple(Rational(x) for x in combo) for combo in product(*ranges)]  # 将所有可能的组合转换为包含 Rational 类型元素的元组
    return D


def count_variables_in_formulas(formula: QFFormula): #  从formula里面找每个公式的变量数量
    tmp = []
    for clause in formula.expr.args:
        tmp.append(len(clause.free_symbols))
    return tmp

def generate_M(Delta, vars_list):
    return [Rational(Delta) ** Rational(vars_count) + Rational(1, 1) for vars_count in vars_list]


# 主函数
def dft(cons: CardinalityConstraint, formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred):
    pred_cons_dict = {list(item[0].keys())[0]: item[2] for item in cons.constraints} # 构建字典，键是谓词，值是基数约束
    pred_arity_dict = {pred: pred.arity for pred in cons.preds} # 构建字典，键是变量，值是它们的 arity 属性 # 只需要关注被约束的谓词的变量的个数
    global n
    n = [] # 基数约束，每个谓词对应的约束
    var_counts = [] # 存储 每个谓词 对应的变量数目
    for _ , item in enumerate(pred_arity_dict):
        n.append(Rational(pred_cons_dict[item]))
        var_counts.append(pred_arity_dict[item])
    # 确定 k n M 约束谓词 的维度length
    global length
    length = len(var_counts)
    # 确定D
    global D #
    D = generate_D(len(domain), var_counts) # 获取D = {0, 1, . . . , |∆||vars(α1)|} × · · · × {0, 1, . . . , |∆||vars(αm)|}.
    D = D[27:]
    # 确定M
    global M
    M = generate_M(len(domain), var_counts) # M = [|∆||vars(α1)| + 1, . . . , |∆||vars(αm)| + 1, 1],
    return q(formula, domain, get_weight, leq_pred)



# get_weight 中的X0是什么