# 从 sampling_fo2.problems 模块中导入了 WFOMCSProblem 类和 MLN_to_WFOMC 函数。这意味着代码使用了外部模块 sampling_fo2.problems，其中定义了一个问题类 WFOMCSProblem 和一个将 MLN 转换为 WFOMC 的函数 MLN_to_WFOMC。
from wfomc.problems import WFOMCProblem, MLN_to_WFOMC
# # 从当前目录下的 wfomcs_parser 模块中导入了 parse 函数，并将其重命名为 wfomcs_parse。这个函数可能是解析 .wfomcs 文件内容的工具。
from .wfomcs_parser import parse as wfomcs_parse
# # 同样从当前目录下的 mln_parser 模块中导入了 parse 函数，并重命名为 mln_parse。这个函数大概率是用于解析 .mln 文件的。
from .mln_parser import parse as mln_parse

# # 定义了一个名为 parse_input 的函数，接受一个 input_file 字符串参数，并返回一个 WFOMCSProblem 对象。该函数用于根据文件类型解析输入文件。
def  parse_input(input_file: str) -> WFOMCProblem:
    if input_file.endswith('.mln'):# 判断输入文件的文件名是否以 .mln 结尾，即是否为 MLN 文件格式。
        with open(input_file, 'r') as f:
            input_content = f.read()
        mln_problem = mln_parse(input_content)
        wfomcs_problem = MLN_to_WFOMC(mln_problem)
        return wfomcs_problem
    elif input_file.endswith('.wfomcs'):
        with open(input_file, 'r') as f:
            input_content = f.read()
        return wfomcs_parse(input_content)
    else:
        raise RuntimeError(f'Unknown input file type: {input_file}')
