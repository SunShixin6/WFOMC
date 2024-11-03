from wfomc.algo import recursive_wfomc
from wfomc.context import WFOMCContext
from wfomc.fol.syntax import Const, Pred, QFFormula
from typing import Callable
from wfomc.utils import RingElement, Rational, k
from itertools import product
import numpy as np
from symengine import Rational, I, exp, symbols, pi # 新导入的
# # DFT:求g（k），对于任意k属于D
# def get_dft(
#         formula: QFFormula,  # # 要求解的逻辑公式。skolem之后的
#         domain: set[Const],
#         get_weight: Callable[[Pred], tuple[RingElement, RingElement]],  # 一个正 一个负
#         leq_pred: Pred,  # leq_pred: 谓词，用于指定“≤”关系。线性阶，这个本质上是一个binary predicate
#         real_version: bool = True) -> RingElement:
#     # 调用WFOMC，只是让里面的权重增加一项而已。相当于在原始WFOMC外面包了一层
#     res = recursive_wfomc(formula, domain, get_weight, leq_pred, True)
#     return res
#     pass
#
#
# # 傅里叶反变换
# def get_reverse_DFT(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
#     # 先调用上面的变换，然后在这个函数里面完成反变换
#     res = get_dft(formula, domain, get_weight, leq_pred, real_version)
#     # 反变换
#
#
#
#     return res
#     pass
#
#
# # 主函数
# def dft(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred, real_version: bool = True):
#     res =get_reverse_DFT(formula, domain, get_weight, leq_pred, real_version)
#     return  res
#     pass

# DFT:求g（k），对于任意k属于D
def get_dft(
        formula: QFFormula,  # # 要求解的逻辑公式。skolem之后的
        domain: set[Const],
        get_weight: Callable[[Pred], tuple[RingElement, RingElement]],  # 一个正 一个负
        leq_pred: Pred,  # leq_pred: 谓词，用于指定“≤”关系。线性阶，这个本质上是一个binary predicate
        k_div_M,
        real_version: bool = True) -> RingElement:
    # 调用WFOMC，只是让里面的权重增加一项而已。相当于在原始WFOMC外面包了一层
    res = recursive_wfomc(formula, domain, get_weight, leq_pred, real_version, k_div_M )
    return res
    pass


# 傅里叶反变换
def get_reverse_DFT(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    tmp_res = Rational(0,1)
    for k in D:
        k = np.array(k)
        k_div_M = k / M
        dot_res = np.dot(n, k / M )
        # 先调用上面的变换，然后在这个函数里面完成反变换
        tmp_res += get_dft(formula, domain, get_weight, leq_pred, k_div_M, real_version)  * np.exp(1j * 2 * np.pi * dot_res) # fixme  这里好像可以直接传入coef 应该传入什么呢，向量还是元素
    return tmp_res
    pass

def generate_D(domain_size, var_counts): # var_counts 是一个列表，其中包含每个公式的变量数量
    """
    生成集合 D，其中 D 是一个多维整数向量的集合。

    参数:
    - domain_size (int): 常量集合 Delta 的大小 |Delta|
    - var_counts (list of int): 每个公式中变量的数量 |vars(α_i)|
    # # 示例用法
    # domain_size = 3  # 假设 |Delta| = 3
    # var_counts = [2, 1, 2]  # 假设公式 α1 有 2 个变量, α2 有 1 个变量, α3 有 2 个变量
    # D = generate_D(domain_size, var_counts)
    # print(D)
    返回:
    - D (list of tuples): 集合 D 的所有可能组合
    """
    ranges = [range(domain_size ** var_count + 1) for var_count in var_counts] # ranges 中的每个 range(domain_size ** var_count + 1) 表示每个公式可以生成的计数范围。
    D = list(product(*ranges)) #  product(*ranges) 生成了所有可能的组合，即集合 𝐷
    return D

def count_variables_in_formulas(formula: QFFormula): #  从formula里面找每个公式的变量数量
    tmp = []
    for clause in formula.expr.args:
        tmp.append(len(clause.free_symbols))
    return tmp

def generate_M(Delta, vars_list):
    return [Delta ** vars_count + 1 for vars_count in vars_list]


# 主函数
def dft(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    # 这里n是基数约束 # TODO 看看代码中哪里处理基数约束？ 某个谓词没有基数约束，代码如何表示？ 去看看decode部分 ,感觉好像要在decode部分处理
    global n
    n = np.array([1,1,1,1,1])

    var_counts = count_variables_in_formulas(formula) #  从formula里面找每个公式的变量数量
    # 确定D
    global D
    D = generate_D(len(domain), var_counts) # 获取D = {0, 1, . . . , |∆||vars(α1)|} × · · · × {0, 1, . . . , |∆||vars(αm)|}.
    # 确定M
    global M
    M = np.array(generate_M(len(domain), var_counts)) # M = [|∆||vars(α1)| + 1, . . . , |∆||vars(αm)| + 1, 1],
    tmp = get_reverse_DFT(formula, domain, get_weight, leq_pred)
    return  tmp / sum(M)


# fixme 代码中什么时候要保留公式，什么时候要带入值呢？ 同构检测的时候，使用公式还是使用具体的值呢？
# fixme 代码都有哪些地方需要加复数，也就是引入正权重的时候，也就是构建图的时候
# get_weight 中的X0是什么