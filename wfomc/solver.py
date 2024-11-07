from __future__ import annotations

import os
import argparse
import logging
import logzero

from logzero import logger
from contexttimer import Timer

from wfomc.algo.DFT import dft
from wfomc.problems import WFOMCProblem
from wfomc.algo import Algo, standard_wfomc, fast_wfomc, incremental_wfomc, recursive_wfomc

from wfomc.utils import MultinomialCoefficients, Rational, round_rational
from wfomc.context import WFOMCContext
from wfomc.parser import parse_input
from wfomc.fol.syntax import Pred

# 这段代码实现了加权一阶模型计数（Weighted First-Order Model Counting, WFOMC）的算法，使用不同的计数方法处理给定的输入
# 定义一个名为 wfomc 的函数，它接收两个参数：problem：一个 WFOMCSProblem 实例，表示需要解决的具体问题。algo：算法选择，默认为 STANDARD。返回值类型为 Rational，即有理数结果。
def wfomc(problem: WFOMCProblem, algo: Algo = Algo.STANDARD) -> Rational: # 这里不仅仅是有理数，还可能是多项式
    # both standard and fast WFOMCs need precomputation
    if algo == Algo.STANDARD or algo == Algo.FAST or \
            algo == algo.FASTv2:
        MultinomialCoefficients.setup(len(problem.domain))

    context = WFOMCContext(problem) # 使用 problem 创建一个 WFOMCContext 上下文对象，用于存储公式、域、权重等问题的相关信息。
    leq_pred = Pred('LEQ', 2) # 创建一个谓词 LEQ，它带有两个参数，表示一个“≤”比较的关系。
    if leq_pred in context.formula.preds(): # 检查当前公式中是否包含 LEQ 谓词
        logger.info('Linear order axiom with the predicate LEQ is found')
        if algo == Algo.RECURSIVE:
            logger.info('Invoke recursive WFOMC')
        else:
            algo = Algo.INCREMENTAL
            logger.info('Invoke incremental WFOMC')
    else:
        leq_pred = None

    with Timer() as t:
        if algo == Algo.STANDARD:
            res = standard_wfomc(
                context.formula, context.domain, context.get_weight
            )
        elif algo == Algo.FAST:
            res = fast_wfomc(
                context.formula, context.domain, context.get_weight
            )
        elif algo == Algo.FASTv2:
            res = fast_wfomc(
                context.formula, context.domain, context.get_weight, True
            )
        elif algo == Algo.INCREMENTAL: # 增量 WFOMC，适合处理包含 LEQ 谓词的情景。
            res = incremental_wfomc(
                context.formula, context.domain,
                context.get_weight, leq_pred
            )
        elif algo == Algo.RECURSIVE:
            res = recursive_wfomc(
                context.formula, context.domain,
                context.get_weight, leq_pred,
            )
        elif algo == Algo.DFT: # TODO 表示使用dft
            res = dft(
                context.cardinality_constraint,
                context.formula, context.domain,
                context.get_weight, leq_pred,
            )
    # res = context.decode_result(res) # 将结果通过上下文解码。
    # logger.info('WFOMC time: %s', t.elapsed) # 记录计算所花费的时间，并返回结果。
    # return res


def parse_args():# 使用 argparse 模块定义了命令行参数：
    parser = argparse.ArgumentParser(
        description='WFOMC for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file') # --inputfile：MLN 输入文件，必需参数。
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points') # --output_dir：输出目录，默认是 ./check-points。可能会输出一些日志之类的
    parser.add_argument('--algo', '-a', type=Algo,
                        choices=list(Algo), default=Algo.FASTv2) # --algo：选择使用的算法，默认是 FASTv2。
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # import sys
    # sys.setrecursionlimit(int(1e6)) # 设置 Python 递归限制为 1,000,000。
    args = parse_args() # 调用 parse_args 获取命令行参数
    if not os.path.exists(args.output_dir): # 检查输出目录是否存在。如果不存在则创建。
        os.makedirs(args.output_dir)
    if args.debug: # 根据 debug 参数设置日志级别，调试模式下记录详细信息，否则记录一般或严重错误信息。
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    logzero.logfile('{}/log.txt'.format(args.output_dir), mode='w') # 设置日志文件路径，并将日志信息写入到指定目录中的 log.txt 文件。

    with Timer() as t: # 使用 parse_input 函数解析输入文件（MLN 文件），并记录解析所需的时间。
        problem = parse_input(args.input) # 把文件转成一个类
    logger.info('Parse input: %ss', t)

    res = wfomc( # 调用 wfomc 函数执行计算，记录结果（有理数精度）。
        problem, algo=args.algo
    )
    # logger.info('WFOMC (arbitrary precision): %s', res)
    # round_val = round_rational(res) # 然后使用 round_rational 对结果进行四舍五入，并记录该四舍五入后的结果。
    # logger.info('WFOMC (round): %s (exp(%s))', round_val, round_val.ln())
