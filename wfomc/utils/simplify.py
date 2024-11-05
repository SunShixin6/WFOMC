from sympy import symbols, exp, I, pi, simplify, Mod, re, im, Add, Mul, Expr
from symengine import sympify as symengine_sympify  # 从 symengine 导入 sympify
import sympy as sp

def my_simplify(weight):
    # 检查是否为单个表达式，如果是，划为一个元素的列表便于后续处理
    if not isinstance(weight, list):
        weight = [weight]

    # 将 symengine 表达式转换为 sympy 表达式
    weight = [sp.sympify(str(expr)) for expr in weight]

    # 定义一个函数来归一化指数部分
    def normalize_expression(expr):
        # print(type(expr))
        # 如果表达式是加法项，要分解并遍历每一项
        if isinstance(expr, Add):
            return Add(*[normalize_expression(arg) for arg in expr.args])
        # 如果表达式是乘法项，需要进一步检查指数
        elif isinstance(expr, Mul):
            new_args = []
            for arg in expr.args:
                if isinstance(arg, Expr) and arg.func == exp:
                    # 获取指数部分并进行规范化
                    exponent = arg.args[0]
                    new_exponent = Mod(re(exponent), 2) * pi + I * (im(exponent) % (2 * pi))
                    new_args.append(exp(new_exponent))
                else:
                    new_args.append(arg)
            return Mul(*new_args)
        # 处理单独的指数项
        elif isinstance(expr, Expr) and expr.func == exp:
            exponent = expr.args[0]
            new_exponent = Mod(re(exponent), 2) * pi + I * (im(exponent) % (2 * pi))
            return exp(new_exponent)
        # 对于其他类型的表达式，直接返回
        else:
            return expr

    # 处理每个表达式
    normalized_expressions = [normalize_expression(expr) for expr in weight]

    # 输出归一化后的表达式
    # for expr in normalized_expressions:
    #     print(expr)
    return normalized_expressions


if __name__ == '__main__':
    x0, x1 = symbols("x0 x1")
    # 定义你的表达式列表
    expressions = [
        exp((-0.0 - 5.33333333333333 * I) * pi) * x0 ** 2 * x1 ** 2 + 2 * exp((-0.0 - 3.33333333333333 * I) * pi) * x0 * x1 + exp((-0.0 - 1.33333333333333 * I) * pi),
        -exp((-0.0 - 2.66666666666667 * I) * pi) * x0 * x1 - exp((-0.0 - 0.666666666666667 * I) * pi),
        -exp((-0.0 - 2.66666666666667 * I) * pi) * x0 * x1 - exp((-0.0 - 0.666666666666667 * I) * pi),
        1  # 这是一个整数
    ]
    print(my_simplify(expressions))
    """
    [x0**2*x1**2*exp(0.66666666666667*I*pi) + 2*x0*x1*exp(0.66666666666667*I*pi) + exp(0.66666666666667*I*pi), 
    -x0*x1*exp(1.33333333333333*I*pi) - exp(1.33333333333333*I*pi), 
    -x0*x1*exp(1.33333333333333*I*pi) - exp(1.33333333333333*I*pi), 
    1]
    """
