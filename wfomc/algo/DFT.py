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

# TODO 每一个g是一个CCG，实际上求一个q 就是求了D的D个CCG，也就是对每一个k算一个CCG，就是对每一个k构建了一个cell_graph,而这些cell graph的weight是不一样，
#       在每个cell graph上做CCG时候计算的同构信息，是在这个cell graph上产生的同构的图。我们猜想虽然产生了size个cell graph，也就是产生了size个nauty的无向图，但是他们之间会有很多同构信息，
#       所以可以把每次创建的nauty无向图保存下来，但是有个问题是nauty graph的边权重不同。也就是说，在两个图上做CCG，如果这两个图一开始边权重就是不同的，再怎么做CCG都不可能有同构的。
#       我们对边染色的时候，边权重不同，颜色也不同，我们可以用一个全局变量，存储边权重和对应的颜色，这是一个方法。
#  或者可以只调用一次CCG，也就是对原始的算法改进.也就是把size个cell graph聚合成一个，这就要改进CCG的定义，就是论文4.2中，式子5，就是有size个CCG加在一起，每个前面都有个e权重，然后把连乘号提到式子最左边，这样图就会很大，并且和domain size相关。不过只需要递归n次，
#  所以现在可以设计一个新的CCG的算法，还是只算一次比较好，因为检查同构也是在n是相等的情况下检测，就算对于不同的CCG，也是在n是相同的情况下去检测，这就需要更改CCG的调用过程。有可能不是一个整的图，可能是一组cell graph，
#  因为求同构的时候对整个大图求同构也没有意义，其实还是对每一个cell 求同构，本质上还是在一个cell graph，而不是一个聚合的cell graph上求同构。我们接下来要构建很多个cell graph合在一起，构建一个大的，求nauty会很慢。
#
#  TODO 另一种方法是调用D次CCG，实现统计一下边颜色，每次取一个k，构建cell grpah，计算label，然后再取下一个k，构建cell graph，然后label。我们用一个全局的变量，储存这个边的权重的映射，比如多少是颜色几、上面两种方法速度不会相差太大，只不过实现起来难度不同。
#   总结来说，对每一个k调用一次CCG，在调用CCG的时候，因为CCG递归的时候，边的颜色是不会变的。调用CCG的一开始的时候，会把边的颜色求出来，求颜色的时候 不仅仅考虑当前的cell graph，还要考虑别的cell gralh。
#   也就是对所有cellgraph构建一个映射。比如一维k=0，对应的边是000，k=1,111，k=2,101，要把1 0 先转化为color，然后再看cell graph的边用什么颜色。
#   现在的代码中只需要处理一个cell graph，就不许保存权重到颜色的映射，就需要把权重换成颜色就可以，但是如果现在有多个cell graph，边的权重不一样了，就要考虑边不同颜色的影响了，就是不同颜色构建一个字典。
#   现在要做一个实验，看一下边的权重有多少种，如果有很多种，那就没戏了， 得想另外的方法，那就是没有同构。我们不期望完全相同，如果第一个图有一条红边两条蓝边，第二个图也有一个红边两个蓝边，那么就是有可能同构的，先算一下边权重种类个数有多少个，


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
    # 确定M
    global M
    M = generate_M(len(domain), var_counts) # M = [|∆||vars(α1)| + 1, . . . , |∆||vars(αm)| + 1, 1],
    return q(formula, domain, get_weight, leq_pred)



# get_weight 中的X0是什么