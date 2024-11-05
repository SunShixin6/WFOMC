from __future__ import annotations

import numpy as np
import pandas as pd
import functools
import networkx as nx
from itertools import product

from typing import Callable, Dict, FrozenSet, Generator, List, Tuple
from logzero import logger
from symengine import symbols
from sympy import Poly
from copy import deepcopy
from wfomc.cell_graph.utils import conditional_on

from wfomc.fol.syntax import AtomicFormula, Const, Pred, QFFormula, a, b, c
from wfomc.utils import Rational, RingElement
from wfomc.utils.multinomial import MultinomialCoefficients
from .components import Cell, TwoTable
from sympy import Rational, I, exp, symbols, pi # 新导入的
from wfomc.utils.simplify import my_simplify

class CellGraph(object):
    # 定义了一个 CellGraph 类，用于表示单元格图（CellGraph），主要用于处理单元格和它们之间的 WMC 计算。
    def __init__(self, formula: QFFormula, # 传入的 formula 是一个量词消去的公式（QFFormula）。
                 get_weight: Callable[[Pred], Tuple[RingElement, RingElement]], # 用于获取给定谓词（Pred）的权重，返回一个二元组 (RingElement, RingElement)。这个函数用于确定加权逻辑。
                 leq_pred: Pred = None): # 可选的谓词，用于比较两个值的大小
        """
        Cell graph that handles cells (1-types) and the WMC between them
        :param sentence QFFormula: the sentence in the form of quantifier-free formula
        :param get_weight Callable[[Pred], Tuple[RingElement, RingElement]]: the weighting function
        :param conditional_formulas List[CNF]: the optional conditional formula appended in WMC computing
        """
        self.formula: QFFormula = formula
        self.get_weight: Callable[[Pred],
                                  Tuple[RingElement, RingElement]] = get_weight # 这里将传入的 get_weight 回调函数赋值给 self.get_weight，使其成为类的一个属性。 get_weight 是一个函数，用于接收一个谓词并返回与其关联的两个 RingElement 作为权重。
        self.leq_pred: Pred = leq_pred # 这个谓词 leq_pred 可能用于在逻辑推理过程中进行大小比较，若为 None，则不启用比较功能。
        self.preds: Tuple[Pred] = tuple(self.formula.preds())  # self.preds 存储了从公式 self.formula 中提取的所有谓词，格式为一个 Tuple[Pred]。
        logger.debug('prednames: %s', self.preds) # 使用 logger.debug 记录一个调试日志，显示谓词的名称。

        # 这里调用了类的内部方法 _ground_on_tuple，该方法的功能可能是对公式进行“grounding”，即将公式中的逻辑量词替换为特定的变量或常量。
        gnd_formula_ab1: QFFormula = self._ground_on_tuple(
            self.formula, a, b
        )
        gnd_formula_ab2: QFFormula = self._ground_on_tuple( # 与之前的公式替换方向相反，这里将 b 和 a 互换。
            self.formula, b, a
        )
        self.gnd_formula_ab: QFFormula = gnd_formula_ab1 & gnd_formula_ab2 # 逻辑与表示只有两个公式都为真时，合并后的公式才为真。这个步骤可能是为了确保在 a 和 b 上的某些条件都满足时才继续推理。
        self.gnd_formula_cc: QFFormula = self._ground_on_tuple(
            self.formula, c
        )

        if self.leq_pred is not None: # 这段代码开始检查类属性 self.leq_pred 是否不为空（即初始化时是否传入了用于比较的谓词）。
            self.gnd_formula_cc = self.gnd_formula_cc & self.leq_pred(c, c) # 首先将 leq_pred(c, c) 添加到 self.gnd_formula_cc 中，将谓词 leq_pred(c, c) 以逻辑与的方式添加到公式 self.gnd_formula_cc 中。
            self.gnd_formula_ab = self.gnd_formula_ab & \
                self.leq_pred(b, a) & \
                (~self.leq_pred(a, b)) # 然后将 leq_pred(b, a) 和 leq_pred(a, b) 添加到 self.gnd_formula_ab 中，这个过程是在给公式加入额外的约束条件。# self.leq_pred(b, a) 表示 b <= a，将其添加到公式中。# (~self.leq_pred(a, b)) 表示 a > b，通过逻辑非运算符 ~ 对 leq_pred(a, b) 取反，将其也加入约束条件中。# 最终的公式 self.gnd_formula_ab 是在前面的公式基础上增加了这些比较条件。
        logger.info('ground a b: %s', self.gnd_formula_ab) # 通过 logger.info 输出 self.gnd_formula_ab 和 self.gnd_formula_cc 公式的内容，便于查看 grounding 操作后的公式。
        logger.info('ground c: %s', self.gnd_formula_cc)

        # build cells # 调用类的内部方法 _build_cells()，生成单元 cells，并将其存储在 self.cells 中，类型为 List[Cell]。
        self.cells: List[Cell] = self._build_cells()

        # filter cells
        logger.info('the number of valid cells: %s',
                    len(self.cells)) # 使用 logger.info 记录生成的有效单元（cells）的数量，便于检查结果。

        logger.info('computing cell weights')
        # 调用 _compute_cell_weights() 方法，为每个单元计算权重并存储在 self.cell_weights 中。
        # self.cell_weights 是一个字典，其键是 Cell 对象，值是多项式（Poly），可能表示与该单元关联的权重。
        self.cell_weights: Dict[Cell, Poly] = self._compute_cell_weights()

        logger.info('computing two table weights')
        # 调用 _build_two_tables() 方法，生成 two_tables，这是一个以 (Cell, Cell) 元组为键，TwoTable 为值的字典。
        # 这些表格可能用于表示两个单元之间的某种关系，可能是用于加权模型计数的中间计算表格。
        self.two_tables: Dict[Tuple[Cell, Cell], TwoTable] = self._build_two_tables()


    # 这是 _ground_on_tuple 方法，用于对传入的公式 formula 进行变量替换（grounding），将公式中的逻辑量词变量替换为传入的常量 args。
    # 该方法通过调用公式的 ground() 方法来实现，返回替换后的公式。
    def _ground_on_tuple(self, formula: QFFormula,
                         c1: Const, c2: Const = None) -> QFFormula:
        """
        将公式 formula 进行“grounding”（对其逻辑量词中的变量进行替换），使其只作用于给定的常量元组（*args）。
        formula: 类型为 QFFormula，表示一个量化自由公式 (Quantifier-Free Formula)。
        c1: 类型为 Const，表示一个常量。c2: 类型为 Const，默认为 None，表示可选的第二个常量。返回值类型为 QFFormula。
        """
        variables = formula.vars() # 通过调用 formula.vars()，获取该公式中涉及的所有变量，并将其存储在 variables 中。
        if len(variables) > 2: # 检查变量数量：如果 variables 中的变量数量超过2个，抛出运行时错误（RuntimeError），并显示错误消息 "Can only ground out FO2"。该方法仅处理一阶逻辑（FO2），最多只能处理2个变量的情况。
            raise RuntimeError(
                "Can only ground out FO2"
            )

        if len(variables) == 1: # 单变量情况：如果 variables 中只有一个变量，则将 c1 作为唯一的常量，并将其存储在列表 constants 中。
            constants = [c1]
        else: # 双变量情况
            if c2 is not None: # 如果变量数量为2个，并且传入了第二个常量 c2，则将 c1 和 c2 分别存入 constants 列表中。
                constants = [c1, c2]
            else: # 如果 c2 为 None，则将 c1 作为两个变量的常量，存储在 constants 列表中。
                constants = [c1, c1]
        substitution = dict(zip(variables, constants)) # 创建替换映射：使用 zip() 将 variables 和 constants 逐一配对，创建一个替换字典 substitution。这个字典表示要将哪些变量替换为哪些常量。
        gnd_formula = formula.substitute(substitution) # 执行替换：通过调用 formula.substitute(substitution)，将 formula 中的变量按照 substitution 中定义的规则替换为常量，生成一个新的量化自由公式 gnd_formula。
        return gnd_formula

    # show 方法：调用该对象的 __str__ 方法，将其字符串表示记录到日志中，方便调试和查看对象的状态。
    def show(self):
        logger.info(str(self))


    def __str__(self):
        s = 'CellGraph:\n'
        s += 'predicates: {}\n'.format(self.preds)
        cell_weight_df = []
        twotable_weight_df = []
        for _, cell1 in enumerate(self.cells):
            cell_weight_df.append(
                [str(cell1), self.get_cell_weight(cell1)]
            )
            twotable_weight = []
            for _, cell2 in enumerate(self.cells):
                # if idx1 < idx2:
                #     twotable_weight.append(0)
                #     continue
                twotable_weight.append(
                    self.get_two_table_weight(
                        (cell1, cell2))
                )
            twotable_weight_df.append(twotable_weight)
        cell_str = [str(cell) for cell in self.cells]
        cell_weight_df = pd.DataFrame(cell_weight_df, index=None,
                                      columns=['Cell', 'Weight'])
        twotable_weight_df = pd.DataFrame(twotable_weight_df, index=cell_str,
                                          columns=cell_str)
        s += 'cell weights: \n'
        s += cell_weight_df.to_markdown() + '\n'
        s += '2table weights: \n'
        s += twotable_weight_df.to_markdown()
        return s

    def __repr__(self):
        return str(self)

    def get_cells(self, cell_filter: Callable[[Cell], bool] = None) -> List[Cell]:
        if cell_filter is None:
            return self.cells
        return list(filter(cell_filter, self.cells))

    @functools.lru_cache(maxsize=None, typed=True)
    def get_cell_weight(self, cell: Cell) -> Poly:# 定义一个名为 get_cell_weight 的方法，接受一个参数 cell，类型为 Cell。该方法返回一个多项式 Poly，表示指定细胞的权重。
        if cell not in self.cell_weights: # 检查 cell 是否存在于 self.cell_weights 字典中，以确定该细胞是否已经计算过权重。
            logger.warning( # 如果细胞不存在，使用 logger 记录一个警告信息，指出该细胞未找到。
                "Cell %s not found", cell
            )
            return 0 # 返回 0，表示该细胞的权重不可用。
        return self.cell_weights.get(cell) # 返回 self.cell_weights 字典中与 cell 相关联的权重值。如果细胞存在，则返回对应的权重。

    def _check_existence(self, cells: Tuple[Cell, Cell]):
        if cells not in self.two_tables:
            raise ValueError(
                f"Cells {cells} not found, note that the order of cells matters!"
            )

    @functools.lru_cache(maxsize=None, typed=True)
    def get_two_table_weight(self, cells: Tuple[Cell, Cell],
                             evidences: FrozenSet[AtomicFormula] = None) -> RingElement:
        # 定义一个名为 get_two_table_weight 的方法，接受两个参数：cells 是一个元组，包含两个 Cell 对象；
        # evidences 是一个可选的 FrozenSet，包含原子公式。该方法返回一个 RingElement，表示这两个细胞之间的权重。
        self._check_existence(cells) # 调用私有方法 _check_existence 来检查提供的 cells 是否在 two_tables 字典中存在。如果不存在，会抛出异常。
        return self.two_tables.get(cells).get_weight(evidences) # 从 self.two_tables 字典中获取与 cells 相关联的 TwoTable 对象，然后调用该对象的 get_weight 方法，传入 evidences 参数，返回这两个细胞之间的权重。

    def get_all_weights(self) -> Tuple[List[RingElement], List[RingElement]]: # 定义一个名为 get_all_weights 的方法，该方法没有参数，返回一个包含两个列表的元组。每个列表的元素类型为 RingElement，表示细胞图中所有细胞的权重。
        cell_weights = [] # 初始化一个空列表 cell_weights，用于存储每个细胞的权重。
        twotable_weights = [] # 初始化另一个空列表 twotable_weights，用于存储细胞对之间的权重。
        for cell_i in self.cells: # 遍历 self.cells 列表中的每个细胞 cell_i。
            cell_weights.append(self.get_cell_weight(cell_i)) # 对于当前细胞 cell_i，调用 get_cell_weight 方法获取其权重，并将其添加到 cell_weights 列表中。
            twotable_weight = [] # 初始化一个空列表 twotable_weight，用于存储当前细胞 cell_i 与其他细胞之间的权重。
            for cell_j in self.cells: # 再次遍历 self.cells 列表中的每个细胞 cell_j，以便计算 cell_i 与 cell_j 之间的权重。
                twotable_weight.append(self.get_two_table_weight( # 调用 get_two_table_weight 方法，获取当前细胞 cell_i 和 cell_j 之间的权重，并将其添加到 twotable_weight 列表中。
                    (cell_i, cell_j)
                ))
            twotable_weights.append(twotable_weight) # 将 twotable_weight 列表添加到 twotable_weights 列表中，记录当前细胞 cell_i 的所有细胞对权重。
        return cell_weights, twotable_weights # 返回一个元组，包含 cell_weights 和 twotable_weights 列表，表示所有细胞的权重和细胞对之间的权重。

    @functools.lru_cache(maxsize=None, typed=True)
    def satisfiable(self, cells: Tuple[Cell, Cell],
                    evidences: FrozenSet[AtomicFormula] = None) -> bool:
        self._check_existence(cells)
        return self.two_tables.get(cells).satisfiable(evidences)

    @functools.lru_cache(maxsize=None)
    def get_two_tables(self, cells: Tuple[Cell, Cell],
                       evidences: FrozenSet[AtomicFormula] = None) \
            -> Tuple[FrozenSet[AtomicFormula], RingElement]:
        self._check_existence(cells)
        return self.two_tables.get(cells).get_two_tables(evidences)

    # 这个函数的目的是构建一个单元（Cell）列表。
    def _build_cells(self):
        cells = [] # 首先，cells 被初始化为空列表，用于存储构建好的 Cell 对象。
        code = {} # code 是一个字典，用来保存每个谓词的布尔值，表示某个谓词在模型中的真假状态。
        for model in self.gnd_formula_cc.models(): # 遍历 self.gnd_formula_cc.models()，其中每个 model 都是公式的一个模型（解释），代表公式在某些特定条件下为真的解释。
            for lit in model: # 对于每个模型，再遍历其包含的文字（lit），每个文字表示谓词及其正负性。
                code[lit.pred] = lit.positive # lit.pred 是谓词，lit.positive 表示该谓词在该模型中的真假状态。将谓词和它的真假状态保存在 code 字典中，键是谓词，值是布尔值。
            cells.append(Cell(tuple(code[p] for p in self.preds), self.preds)) # 将当前模型中所有谓词的真假状态按顺序保存为一个元组 tuple(code[p] for p in self.preds)，并创建一个 Cell 对象，包含这个元组和相应的谓词列表 self.preds。
        return cells # 最终，函数返回构建好的单元（cells）列表。

    def _compute_cell_weights(self): # 该函数用于计算每个单元（cell）的权重
        weights = dict() # 首先初始化一个空字典 weights，用于存储每个单元的计算结果。
        for cell in self.cells: # 遍历类中的所有单元 self.cells。
            weight = Rational(1, 1) # 为每个单元初始化一个默认的权重 Rational(1, 1)，即权重为 1
            for i, pred in zip(cell.code, cell.preds): # 遍历单元的 code 和 preds（谓词），zip 将单元的布尔值与相应的谓词配对。
                assert pred.arity > 0, "Nullary predicates should have been removed" # 如果谓词的元数（arity）大于 0，即谓词涉及至少一个参数：
                if i: # 如果 i（单元的代码中的某个值）为真，
                    if USE_DFT:
                        weight = weight * self.get_weight(pred)[0] * exp( - I * 2 * pi * coef) # TODO coef是一个列表，那么也会生成一个列表
                    else:
                        weight = weight * self.get_weight(pred)[0]  # 使用 self.get_weight(pred)[0]，即获取该谓词的第一个权重，并将其乘以当前的 weight。
                else: # 如果 i 为假，
                    weight = weight * self.get_weight(pred)[1] # 使用 self.get_weight(pred)[1]，即获取第二个权重并相乘。
            weights[cell] = weight # 将计算出的权重存储到 weights 字典中，键是单元，值是其对应的权重。
        return weights # 最后返回计算得到的权重字典。

    @functools.lru_cache(maxsize=None) # 这是一个装饰器，使用了 Python 的 lru_cache（Least Recently Used Cache），用于缓存函数的返回值，从而避免多次重复计算，提高性能。maxsize=None 表示缓存的数量没有限制。
    def get_nullary_weight(self, cell: Cell) -> RingElement: # 参数 cell 是一个类型为 Cell 的对象，返回类型为 RingElement。该方法的作用是计算给定单元的零元权重（nullary weight）。
        weight = Rational(1, 1)
        for i, pred in zip(cell.code, cell.preds): # cell.code 是一个布尔值列表，cell.preds 是谓词（predicate）列表。zip 会将 cell.code 和 cell.preds 一一配对，分别赋值给 i 和 pred。
            if pred.arity == 0: # 如果当前谓词的元数（arity）为 0（即没有参数），继续进行计算
                if i:  # 如果 i 为 True，意味着当前的谓词处于“真”的状态。
                    weight = weight * self.get_weight(pred)[0] # 获取当前谓词的第一个权重，并将其与当前权重相乘。
                else: # 如果 i 为 False，意味着当前谓词处于“假”的状态。
                    weight = weight * self.get_weight(pred)[1] # 获取当前谓词的第二个权重（表示“假”状态时的权重），并将其与当前权重相乘。
        return weight

    def _build_two_tables(self): # 这个函数的作用是构建 "two_tables"。
        # 构建一个包含所有模型和权重的pd.DataFrame
        models = dict() # 首先初始化一个空字典 models，用于存储每个模型及其对应的权重。
        gnd_lits = self.gnd_formula_ab.atoms() # gnd_lits 存储了公式 self.gnd_formula_ab 中的所有原子命题（谓词及其变量）。
        gnd_lits = gnd_lits.union(
            frozenset(map(lambda x: ~x, gnd_lits))
        ) # 通过使用 map 函数和取反运算符 ~，将每个原子命题的正负版本都加入到 gnd_lits 中，确保包括每个原子命题的正负两种形式。
        for model in self.gnd_formula_ab.models(): # 遍历公式 self.gnd_formula_ab 的所有模型（解释）。 # models 该方法会生成公式的所有模型
            weight = Rational(1, 1)
            for lit in model: # 对每个模型中的每个文字（谓词及其真假性）：
                # ignore the weight appearing in cell weight # 忽略那些在单元权重中已经考虑的权重（通过检查谓词的参数）。
                if (not (len(lit.args) == 1 or all(arg == lit.args[0] for arg in lit.args))):
                     # 根据文字的正负性，选择对应的权重进行相乘
                    if lit.positive: # 如果谓词为正
                        if USE_DFT:
                            weight *= self.get_weight(lit.pred)[0] * exp( - I * 2 * pi * coef) # 使用第一个权重
                        else:
                            weight *= self.get_weight(lit.pred)[0]  # 使用第一个权重
                    else: # ；如果谓词为负
                        weight *= self.get_weight(lit.pred)[1] # 使用第二个权重
            models[frozenset(model)] = weight # 将模型及其计算的权重添加到 models 字典中，键为模型的 frozenset（使模型顺序无关），值为对应的权重。

        # build twotable tables
        tables = dict() # 初始化一个空字典 tables，用于存储生成的 "two_tables"。
        for i, cell in enumerate(self.cells): # 遍历所有单元 self.cells，
            models_1 = conditional_on(models, gnd_lits, cell.get_evidences(a)) # 对每个单元 cell 调用 conditional_on 函数，根据模型和单元的证据（get_evidences(a)），计算条件模型 models_1
            for j, other_cell in enumerate(self.cells): # 对于每个单元 cell，再次遍历 self.cells 中的另一个单元 other_cell，检查它们的相对位置：
                # NOTE: leq is sensitive to the order of cells
                if i > j and self.leq_pred is None: # 如果 i > j（即单元 cell 在单元 other_cell 之后），则使用已经计算过的结果，避免重复计算。
                    tables[(cell, other_cell)] = tables[(other_cell, cell)]
                models_2 = conditional_on(models_1, gnd_lits, # 否则，调用 conditional_on 函数，进一步基于 other_cell 的证据计算新的条件模型 models_2。
                                          other_cell.get_evidences(b))
                tables[(cell, other_cell)] = TwoTable(models_2, gnd_lits) # 对于每对单元 (cell, other_cell)，创建一个 TwoTable 对象，并将其存储到 tables 字典中，键为单元对，值为计算出的 TwoTable 对象
        return tables # 返回构建好的 "two_tables" 字典。


class OptimizedCellGraph(CellGraph):
    def __init__(self, formula: QFFormula,
                 get_weight: Callable[[Pred], Tuple[RingElement, RingElement]],
                 domain_size: int,
                 modified_cell_symmetry: bool = False):
        """
        Optimized cell graph for FastWFOMC
        :param formula: the formula to be grounded
        :param get_weight: a function that returns the weight of a predicate
        :param domain_size: the domain size
        """
        super().__init__(formula, get_weight)
        self.modified_cell_symmetry = modified_cell_symmetry
        self.domain_size: int = domain_size
        MultinomialCoefficients.setup(self.domain_size)

        if self.modified_cell_symmetry:
            i1_ind_set, i2_ind_set, nonind_set = self.find_independent_sets()
            self.cliques, [self.i1_ind, self.i2_ind, self.nonind] = \
                self.build_symmetric_cliques_in_ind([i1_ind_set, i2_ind_set, nonind_set])
            self.nonind_map: dict[int, int] = dict(zip(self.nonind, range(len(self.nonind))))
        else:
            self.cliques: list[list[Cell]] = self.build_symmetric_cliques()
            self.i1_ind, self.i2_ind, self.ind, self.nonind \
                = self.find_independent_cliques()
            self.nonind_map: dict[int, int] = dict(
                zip(self.nonind, range(len(self.nonind))))

        logger.info("Found i1 independent cliques: %s", self.i1_ind)
        logger.info("Found i2 independent cliques: %s", self.i2_ind)
        logger.info("Found non-independent cliques: %s", self.nonind)

        self.term_cache = dict()

    def build_symmetric_cliques(self) -> List[List[Cell]]:
        cliques: list[list[Cell]] = []
        cells = deepcopy(self.get_cells())
        while len(cells) > 0:
            cell = cells.pop()
            clique = [cell]
            for other_cell in cells:
                if self._matches(clique, other_cell):
                    clique.append(other_cell)
            for other_cell in clique[1:]:
                cells.remove(other_cell)
            cliques.append(clique)
        cliques.sort(key=len)
        logger.info("Built %s symmetric cliques: %s", len(cliques), cliques)
        return cliques

    def build_symmetric_cliques_in_ind(self, cell_indices_list) -> \
            tuple[list[list[Cell]], list[list[int]]]:
        i1_ind_set = deepcopy(cell_indices_list[0])
        cliques: list[list[Cell]] = []
        ind_idx: list[list[int]] = []
        for cell_indices in cell_indices_list:
            idx_list = []
            while len(cell_indices) > 0:
                cell_idx = cell_indices.pop()
                clique = [self.cells[cell_idx]]
                # for cell in I1 independent set, we dont need to built sysmmetric cliques
                if cell_idx not in i1_ind_set:
                    for other_cell_idx in cell_indices:
                        other_cell = self.cells[other_cell_idx]
                        if self._matches(clique, other_cell):
                            clique.append(other_cell)
                    for other_cell in clique[1:]:
                        cell_indices.remove(self.cells.index(other_cell))
                cliques.append(clique)
                idx_list.append(len(cliques) - 1)
            ind_idx.append(idx_list)
        logger.info("Built %s symmetric cliques: %s", len(cliques), cliques)
        return cliques, ind_idx

    def find_independent_sets(self) -> tuple[list[int], list[int], list[int], list[int]]:
        g = nx.Graph()
        g.add_nodes_from(range(len(self.cells)))
        for i in range(len(self.cells)):
            for j in range(i + 1, len(self.cells)):
                if self.get_two_table_weight(
                        (self.cells[i], self.cells[j])
                ) != Rational(1, 1):
                    g.add_edge(i, j)

        self_loop = set()
        for i in range(len(self.cells)):
            if self.get_two_table_weight((self.cells[i], self.cells[i])) != Rational(1, 1):
                self_loop.add(i)

        non_self_loop = g.nodes - self_loop
        if len(non_self_loop) == 0:
            i1_ind = set()
        else:
            i1_ind = set(nx.maximal_independent_set(g.subgraph(non_self_loop)))
        g_ind = set(nx.maximal_independent_set(g, nodes=i1_ind))
        i2_ind = g_ind.difference(i1_ind)
        non_ind = g.nodes - i1_ind - i2_ind
        logger.info("Found i1 independent set: %s", i1_ind)
        logger.info("Found i2 independent set: %s", i2_ind)
        logger.info("Found non-independent set: %s", non_ind)
        return list(i1_ind), list(i2_ind), list(non_ind)

    def find_independent_cliques(self) -> tuple[list[int], list[int], list[int], list[int]]:
        g = nx.Graph()
        g.add_nodes_from(range(len(self.cliques)))
        for i in range(len(self.cliques)):
            for j in range(i + 1, len(self.cliques)):
                if self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[j][0])
                ) != Rational(1, 1):
                    g.add_edge(i, j)

        self_loop = set()
        for i in range(len(self.cliques)):
            for j in range(self.domain_size):
                if self.get_J_term(i, j) != Rational(1, 1):
                    self_loop.add(i)
                    break

        non_self_loop = g.nodes - self_loop
        if len(non_self_loop) == 0:
            g_ind = set()
        else:
            g_ind = set(nx.maximal_independent_set(g.subgraph(non_self_loop)))
        i2_ind = g_ind.intersection(self_loop)
        i1_ind = g_ind.difference(i2_ind)
        non_ind = g.nodes - i1_ind - i2_ind
        return list(i1_ind), list(i2_ind), list(g_ind), list(non_ind)

    def _matches(self, clique, other_cell) -> bool:
        cell = clique[0]
        if not self.modified_cell_symmetry:
            if self.get_cell_weight(cell) != self.get_cell_weight(other_cell) or \
                    self.get_two_table_weight((cell, cell)) != self.get_two_table_weight((other_cell, other_cell)):
                return False

        if len(clique) > 1:
            third_cell = clique[1]
            r = self.get_two_table_weight((cell, third_cell))
            for third_cell in clique:
                if r != self.get_two_table_weight((other_cell, third_cell)):
                    return False

        for third_cell in self.get_cells():
            if other_cell == third_cell or third_cell in clique:
                continue
            r = self.get_two_table_weight((cell, third_cell))
            if r != self.get_two_table_weight((other_cell, third_cell)):
                return False
        return True

    def setup_term_cache(self):
        self.term_cache = dict()

    def get_term(self, iv: int, bign: int, partition: tuple[int]) -> RingElement:
        if (iv, bign) in self.term_cache:
            return self.term_cache[(iv, bign)]

        if iv == 0:
            accum = Rational(0, 1)
            for j in self.i1_ind:
                tmp = self.get_cell_weight(self.cliques[j][0])
                for i in self.nonind:
                    tmp = tmp * self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[j][0])) ** partition[self.nonind_map[i]]
                accum = accum + tmp
            accum = accum ** (self.domain_size - sum(partition) - bign)
            self.term_cache[(iv, bign)] = accum
            return accum
        else:
            sumtoadd = 0
            s = self.i2_ind[len(self.i2_ind) - iv]
            for nval in range(self.domain_size - sum(partition) - bign + 1):
                smul = MultinomialCoefficients.comb(
                    self.domain_size - sum(partition) - bign, nval
                )
                smul = smul * self.get_J_term(s, nval)
                if not self.modified_cell_symmetry:
                    smul = smul * self.get_cell_weight(self.cliques[s][0]) ** nval

                for i in self.nonind:
                    smul = smul * self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[s][0])
                    ) ** (partition[self.nonind_map[i]] * nval)
                smul = smul * self.get_term(
                    iv - 1, bign + nval, partition
                )
                sumtoadd = sumtoadd + smul
            self.term_cache[(iv, bign)] = sumtoadd
            return sumtoadd

    @functools.lru_cache(maxsize=None)
    def get_J_term(self, l: int, nhat: int) -> RingElement:
        if len(self.cliques[l]) == 1:
            thesum = self.get_two_table_weight(
                (self.cliques[l][0], self.cliques[l][0])
            ) ** (int(nhat * (nhat - 1) / 2))
            if self.modified_cell_symmetry:
                thesum = thesum * self.get_cell_weight(self.cliques[l][0]) ** nhat
        else:
            thesum = self.get_d_term(l, nhat)
        return thesum

    @functools.lru_cache(maxsize=None)
    def get_d_term(self, l: int, n: int, cur: int = 0) -> RingElement:
        clique_size = len(self.cliques[l])
        r = self.get_two_table_weight((self.cliques[l][0], self.cliques[l][1]))
        s = self.get_two_table_weight((self.cliques[l][0], self.cliques[l][0]))
        if cur == clique_size - 1:
            if self.modified_cell_symmetry:
                w = self.get_cell_weight(self.cliques[l][cur]) ** n
                s = self.get_two_table_weight((self.cliques[l][cur], self.cliques[l][cur]))
                ret = w * s ** MultinomialCoefficients.comb(n, 2)
            else:
                ret = s ** MultinomialCoefficients.comb(n, 2)
        else:
            ret = 0
            for ni in range(n + 1):
                mult = MultinomialCoefficients.comb(n, ni)
                if self.modified_cell_symmetry:
                    w = self.get_cell_weight(self.cliques[l][cur]) ** ni
                    s = self.get_two_table_weight((self.cliques[l][cur], self.cliques[l][cur]))
                    mult = mult * w
                mult = mult * (s ** MultinomialCoefficients.comb(ni, 2))
                mult = mult * r ** (ni * (n - ni))
                mult = mult * self.get_d_term(l, n - ni, cur + 1)
                ret = ret + mult
        return ret


def build_cell_graphs(formula: QFFormula, # formula: 量化自由公式 (QFFormula)。
                      get_weight: Callable[[Pred],
                                           Tuple[RingElement, RingElement]], # get_weight: 一个函数，接受一个谓词并返回其权重。
                      leq_pred: Pred = None,
                      use_dft = False, #表示是否使用dft构建
                      k_div_M = None,
                      optimized: bool = False, # 当 optimized 为 True 时，构建优化的 cell graph；否则构建标准的 cell graph。
                      domain_size: int = 0, # domain_size: 整数，表示领域大小。
                      # 定义一个布尔参数 modified_cell_symmetry，默认为 False。用于控制优化模式下的 cell symmetry（单元对称性）修改行为。
                      modified_cell_symmetry: bool = False) \
        -> Generator[tuple[CellGraph, RingElement]]:  # 方法的返回类型是一个生成器，生成 CellGraph 和 RingElement 类型的元组。

    global USE_DFT # 全局声明是否用DFT，就不用来回传参了
    USE_DFT = use_dft # 全局声明是否用DFT，就不用来回传参了
    if use_dft:
        global coef
        coef = k_div_M

    nullary_atoms = [atom for atom in formula.atoms() if atom.pred.arity == 0] # 从 formula 中获取所有空元谓词（arity 为 0 的谓词），并将它们存储在 nullary_atoms 列表中。
    if len(nullary_atoms) == 0: # 判断是否存在空元谓词。如果没有空元谓词，则进入该条件分支。
        logger.info('No nullary atoms found, building a single cell graph') # 记录日志信息，说明没有找到空元谓词，将构建一个单一的 cell graph。
        if not optimized: # 检查 optimized 参数。如果 optimized 为 False，则按照标准模式构建 cell graph。
            yield CellGraph(formula, get_weight, leq_pred), Rational(1, 1) # 生成一个 CellGraph 实例及其对应的权重 Rational(1, 1)（表示权重为 1）。
        else: # 如果 optimized 为 True，进入此分支，构建一个优化后的 cell graph。
            yield OptimizedCellGraph( #
                formula, get_weight, domain_size, modified_cell_symmetry
            ), Rational(1, 1) # 生成一个 OptimizedCellGraph 实例及其对应的权重 Rational(1, 1)。
    else: # 如果存在空元谓词，则进入此分支。
        logger.info('Found nullary atoms %s', nullary_atoms) # 记录日志，说明找到空元谓词，并输出这些谓词的详细信息。
        for values in product(*([[True, False]] * len(nullary_atoms))): # 遍历所有空元谓词的布尔值组合，使用笛卡尔积生成不同的布尔组合（每个谓词取 True 或 False）。
            substitution = dict(zip(nullary_atoms, values)) # 将布尔组合与空元谓词绑定生成字典 substitution，表示谓词的特定布尔值替换方案。
            logger.info('Building cell graph with values %s', substitution) # 记录日志信息，说明当前正在基于布尔值替换方案 substitution 构建 cell graph。
            subs_formula = formula.sub_nullary_atoms(substitution).simplify() # 使用 substitution 替换 formula 中的空元谓词，并简化得到新的公式 subs_formula。
            if not subs_formula.satisfiable(): # 检查 subs_formula 是否可满足。如果不可满足，则进入此分支。
                logger.info('Formula is unsatisfiable, skipping') # 记录日志，说明该替换方案下的公式不可满足，因此跳过该方案。
                continue # 跳过不可满足的替换方案，继续下一个替换组合。
            if not optimized: # 检查 optimized 参数，如果为 False，则按照标准模式构建 cell graph。
                cell_graph = CellGraph(subs_formula, get_weight, leq_pred) # 创建一个标准的 CellGraph 实例 cell_graph，基于替换后的公式 subs_formula。
            else: # 如果 optimized 为 True，进入此分支，构建优化后的 cell graph。
                cell_graph = OptimizedCellGraph(
                    subs_formula, get_weight, domain_size,
                    modified_cell_symmetry
                ) # 创建一个 OptimizedCellGraph 实例 cell_graph，基于替换后的公式 subs_formula。
            weight = Rational(1, 1) # 初始化 weight 为 Rational(1, 1)，表示初始权重为 1。
            for atom, val in zip(nullary_atoms, values): # 遍历每个空元谓词及其布尔值组合，计算当前组合下的权重。
                weight = weight * (get_weight(atom.pred)[0] if val else get_weight(atom.pred)[1]) # 根据谓词布尔值 val，从 get_weight 中获取对应的权重并累乘至 weight。
            yield cell_graph, weight # 生成一个包含 cell_graph 和计算得出的权重 weight 的元组。
