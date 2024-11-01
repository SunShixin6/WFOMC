from __future__ import annotations # 引入未来特性annotations，使得函数返回值和参数类型的注解可以使用 Python 3.9+ 的新特性。即使当前 Python 版本低于3.9，也可以使用新语法。
# 作用：导入 functools 模块，该模块包含一些用于函数操作的工具，这里会用到它的 lru_cache 装饰器。
import functools

# 生成一个由 length 个数字组成的元组，其元素之和等于 total_sum。
def multinomial(length: int, total_sum: int) -> tuple[int]:
    """
    Generate a list of numbers, whose size is `length` and sum is `total_sum`

    :param length int: length of the generated list
    :param total_sum int: the summation over the list
    :rtype tuple[int]:
    """
    if length == 1: # 如果长度为 1，则直接返回一个元组，元素是 total_sum。
        yield (total_sum, )
    else: # 对于长度大于 1 的情况，遍历从 0 到 total_sum 的值，并递归调用 multinomial 函数生成其他部分的排列组合。
        for value in range(total_sum + 1):
            for permutation in multinomial(length - 1, total_sum - value):
                yield (value, ) + permutation

# 生成一个由 length 个数字组成的元组，其元素之和小于 total_sum。
def multinomial_less_than(length: int, total_sum: int) -> tuple[int]:
    """
    Generate a list of numbers, whose size is `length` and sum is less than `total_sum`

    :param length int: length of the generated list
    :param total_sum int: the summation over the list
    :rtype tuple[int]:
    """
    if length == 0: # 如果长度为 0，返回一个空元组。
        yield ()
        return
    if length == 1: # 如果长度为 1，返回从 0 到 total_sum 的单个数字元组。
        for i in range(total_sum + 1):
            yield (i, )
    else: # 与 multinomial 类似，递归生成元组，其和小于 total_sum。
        for value in range(total_sum + 1):
            for permutation in multinomial_less_than(length - 1, total_sum - value):
                yield (value, ) + permutation


class MultinomialCoefficients(object):
    """
    Multinomial coefficients

    Usage:
    ```
    MultinomialCoefficients.precompute_pascal(n)
    ...
    MultinomialCoefficients.coef(list)
    ```


    """
    pt: list[list[int]] = None
    n: int = 0

    @staticmethod
    # @jit
    def setup(n: int):
        """
        Pre-compute the pascal triangle.

        :param n int: the maximal total sum
        """
        if n <= MultinomialCoefficients.n:
            return
        pt: list[list[int]] = []
        lst: list[int] = [1]
        for i in range(n + 1):
            pt.append(lst)
            newlist = []
            newlist.append(lst[0])
            for i in range(len(lst) - 1):
                newlist.append(lst[i] + lst[i + 1])
            newlist.append(lst[-1])
            lst = newlist
        MultinomialCoefficients.pt = pt
        MultinomialCoefficients.n = n

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def coef(lst: tuple[int]) -> int:
        """
        Compute the multinomial coefficient of `lst`
        """
        if MultinomialCoefficients.pt is None:
            raise RuntimeError(
                'Please initialize MultinomialCoefficients first by `MultinomialCoefficients.setup(n)`'
            )
        if sum(lst) > MultinomialCoefficients.n:
            raise RuntimeError(
                'The sum %d of input is larger than precomputed maximal sum %d, '
                'please re-initialized MultinomialCoefficients using bigger n',
                sum(lst), MultinomialCoefficients.n
            )
        ret = 1
        tmplist = lst
        while len(tmplist) > 1:
            ret *= MultinomialCoefficients.comb(sum(tmplist), tmplist[-1])
            tmplist = tmplist[:-1]
        return ret

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def comb(a, b):
        if a < b:
            return 0
        elif b == 0:
            return 1
        else:
            return MultinomialCoefficients.pt[a][b]
