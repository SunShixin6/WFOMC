from wfomc.algo import recursive_wfomc
from wfomc.context import WFOMCContext
from wfomc.fol.syntax import Const, Pred, QFFormula
from typing import Callable
from wfomc.utils import RingElement, Rational
from itertools import product
import numpy as np
from sympy import Rational, I, exp, symbols, pi # 新导入的
from wfomc.utils.simplify import my_simplify

def g(formula: QFFormula, domain: set[Const], get_weight: Callable[[Pred], tuple[RingElement, RingElement]], leq_pred: Pred, ki_div_Mi, real_version: bool = True) -> RingElement:
    # 调用WFOMC，只是让里面的权重增加一项而已。相当于在原始WFOMC外面包了一层
    res = recursive_wfomc(formula, domain, get_weight, leq_pred, real_version, ki_div_Mi, ues_dft = True ) # ki_div_Mi, ues_dft = True 这两个参数是新加的
    return res


# 傅里叶反变换
def q(formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    tmp_list = [Rational(0, 1) for _ in range(length)] # 生成一个列表，长度为length，也就是k的长度，每个元素都是Rational(0, 1)
    D1 = D[1000000:] # 从100000开始
    tmp_cache = dict()
    for k in D1: # k:(2, 4, 23, 7, 14)
        dot_res = np.dot(n, np.array(k) / np.array(M))  # n 和 k/M 做点乘
        for index in range(length):
            ki_div_Mi = k[index] / M[index] # 获得每一个ki/Mi, 传入CCG中
            e_coef = exp(I * 2 * pi * dot_res)
            e_coef_simple = my_simplify(e_coef)[0]
            if ki_div_Mi not in tmp_cache: # 没有计算过ki_div_Mi的CCG
                tmp_res = g(formula, domain, get_weight, leq_pred, ki_div_Mi, real_version)  * e_coef_simple# 先调用上面的变换，然后在这个函数里面完成反变换 # fixme  这里好像可以直接传入coef 应该传入什么呢，向量还是元素
                tmp_cache[ki_div_Mi] = tmp_res
            else:
                tmp_res = tmp_cache[ki_div_Mi] # 计算过ki_div_Mi的CCG
            tmp_list[index] += tmp_res #累加
    return tmp_list / sum(M)


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



# 获取基数约束数组
def get_n(context):
    # fixme 定义一个字典，里面对应的是每个谓词和约束
    if not context.contain_cardinality_constraint(): #没有约束，返回None
        return None

    preds = context.cardinality_constraint.preds #
    constraints = context.cardinality_constraint.constraints
    result_dict = {list(item[0].keys())[0]: item[2] for item in constraints}
    pass



# 主函数
def dft(context: WFOMCContext, formula: QFFormula,domain: set[Const],get_weight: Callable[[Pred], tuple[RingElement, RingElement]],leq_pred: Pred,real_version: bool = True):
    # 这里n是基数约束
    get_n(context) # fixme 看看代码中哪里处理基数约束？ 某个谓词没有基数约束，代码如何表示？ 去看看decode部分 ,感觉好像要在decode部分处理
    global n
    n = np.array([1,1,1,1,1]) # fixme  没有约束的话，这里是None吗

    var_counts = count_variables_in_formulas(formula) # fixme  是应该从formula里面找每个公式的变量数量吗，这里面formula有五个，但是原始文件中有4个
    global length
    length = len(var_counts)
    # 确定D
    global D # fixme D 和J 搞不清楚
    D = generate_D(len(domain), var_counts) # 获取D = {0, 1, . . . , |∆||vars(α1)|} × · · · × {0, 1, . . . , |∆||vars(αm)|}.
    # 确定M
    global M
    M = generate_M(len(domain), var_counts) # M = [|∆||vars(α1)| + 1, . . . , |∆||vars(αm)| + 1, 1],

    return q(formula, domain, get_weight, leq_pred)



# get_weight 中的X0是什么