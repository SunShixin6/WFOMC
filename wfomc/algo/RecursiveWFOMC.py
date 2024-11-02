import math
from collections import Counter # collections.Counter: Counter 是 collections 模块中的一个类，用来统计可哈希对象的数量，比如可以用来统计列表中元素出现的次数。例子：Counter([1, 2, 2, 3]) 会返回 Counter({2: 2, 1: 1, 3: 1})，表示 2 出现了两次，1 和 3 各出现一次。
from typing import Callable # 是 Python 的类型提示系统，表示一个可以调用的对象（比如函数）。它可以帮助我们在编写函数时，明确指定哪些参数需要是函数类型。例子：Callable[[int, int], str] 表示一个接受两个 int 参数并返回 str 类型的函数。
import pynauty # 导入用于图同构检查的库，该库可以帮助判断两个图是否是同构的。

from sympy import factorint, factor_list # sympy 是一个符号数学库，提供了factorint和factor_list用于因式分解。factorint: 将整数进行因式分解，返回分解结果。例如，factorint(12) 会返回 {2: 2, 3: 1}，表示 12 = 2² × 3。factor_list: 对多项式进行因式分解。它返回一个分解后的多项式列表。
from wfomc.cell_graph import CellGraph, build_cell_graphs
from wfomc.utils import RingElement, Rational
from wfomc.utils.polynomial import expand # expand: 用于将多项式展开为更加易于计算的形式。
from wfomc.fol.syntax import Const, Pred, QFFormula # Const、Pred 和 QFFormula: 这些是用于逻辑运算的一阶逻辑对象：Const: 表示逻辑中的常量。Pred: 表示谓词，它定义了一些逻辑关系。QFFormula: 表示量化逻辑公式，可能包含存在量化符号或者通用量化符号（比如 ∃、∀）。

# TreeNode类表示树中的一个节点，用于递归结构中的表示。它有三个主要属性：
class TreeNode(object):
    def __init__(self, cell_weights, depth):
        self.cell_weights = cell_weights # cell_weights: 节点的权重信息（与图的“单元格”权重有关）。
        self.depth = depth # depth: 节点的深度，表示节点在树中的层级，用来追踪递归中的层级。
        self.cell_to_children = dict[int, TreeNode]() # cell_to_children: 一个字典，键是整数（通常表示某种编号或 ID），值是TreeNode，表示该节点的子节点。

# print_tree 是递归函数，用于以缩进格式打印树结构：
def print_tree(node: TreeNode):
    print("  " * node.depth, node.depth, " ", node.cell_weights) # 打印当前节点的深度和它的权重（cell_weights）。深度越大，缩进越多。
    for k,v in node.cell_to_children.items():  # 遍历该节点的子节点，并递归调用print_tree来打印每个子树。
        print_tree(v)

class IsomorphicGraphCache(object): # 这个类用于实现同构图的缓存，避免对相同的图重复计算。
    def init(self, domain_size: int): # 定义一个名为 init 的方法，接收一个整数参数 domain_size，表示缓存的大小。。
        self.cache = [] # 列表中的每个元素是一个字典，字典用于存储不同层级的图的同构信息。
        self.cache_hit_count = []  # 用于记录每个层级缓存命中的次数。
        for _ in range(domain_size): # 使用 for 循环遍历 domain_size 范围，创建 domain_size 个层级的缓存。
            self.cache.append({}) # 为每个层级添加一个空字典到 cache 列表中，字典将用于存储不同颜色和标签的缓存结果。
            self.cache_hit_count.append(0) # 为每个层级的缓存命中计数器初始化为0，表示该层级的缓存尚未命中。

    # 定义一个名为 get 的方法，接收参数 level（层级）、color_kind（颜色类型元组）、color_count（颜色计数元组）和 can_label（是否可以标记）。该方法用于获取缓存中的值。
    def get(self, level: int, color_kind: tuple[int], color_count: tuple[int], can_label):
        if color_kind not in self.cache[level]: # 检查给定的 color_kind 是否存在于指定层级的缓存中。
            self.cache[level][color_kind] = {} # 如果 color_kind 不在缓存中，则在当前层级的缓存中为其创建一个空字典。
            return None
        if color_count not in self.cache[level][color_kind]: # 检查给定的 color_count 是否存在于 color_kind 的缓存中。
            self.cache[level][color_kind][color_count] = {}
            return None
        if can_label not in self.cache[level][color_kind][color_count]: # 检查 can_label 是否存在于 color_count 的缓存中。
            return None
        self.cache_hit_count[level] += 1 # 如果所有检查都通过，表示命中缓存，则将当前层级的命中计数器加1。
        return self.cache[level][color_kind][color_count][can_label] # 返回对应的缓存值。

    # 定义一个名为 set 的方法，接收参数 level（层级）、color_kind（颜色类型元组）、color_count（颜色计数元组）、can_label（是否可以标记）和 value（要存储的值）。该方法用于将值存储到缓存中。
    def set(self, level: int, color_kind: tuple[int], color_count: tuple[int], can_label, value):
        if color_kind not in self.cache[level]: # 检查给定的 color_kind 是否存在于指定层级的缓存中。
            self.cache[level][color_kind] = {}  # 如果 color_kind 不在缓存中，则在当前层级的缓存中为其创建一个空字典。
        if color_count not in self.cache[level][color_kind]: # 检查给定的 color_count 是否存在于 color_kind 的缓存中。
            self.cache[level][color_kind][color_count] = {}  # 如果 color_count 不在缓存中，则在 color_kind 的字典中为其创建一个空字典。
        self.cache[level][color_kind][color_count][can_label] = value # 将给定的 value 存储到对应层级的缓存中，键为 can_label。
    # 综上所述，IsomorphicGraphCache 类用于缓存图的同构信息，通过 get 和 set 方法来获取和存储缓存数据，以提高同构图检测的性能。


# only for Alog 'dfs_wfomc_real'
PRINT_TREE = False
ROOT = TreeNode([], 0)
# 缓存同构图
IG_CACHE = IsomorphicGraphCache()
# 原细胞图的权值邻接矩阵
ORI_WEIGHT_ADJ_MAT = []
# 转换边色图时需要多少层
LAYERS_NUM_FOR_CONVERT = 0
# 用于记录边的颜色数量，初始值为 0。
EDGE_COLOR_NUM = 0
# 是一个字典，用于将边的权重映射到边的颜色上。初始映射为 {1: 0}，即权重为 1 的边映射到颜色 0。
EDGE_WEIGHT2COLOR_MAP = {1:0}
# 是原始单元图的边颜色矩阵。这个矩阵是通过将 ORI_WEIGHT_ADJ_MAT（权重邻接矩阵）和 EDGE_WEIGHT2COLOR_MAP（边权重到颜色的映射）组合得到的。# edge color matrix of original cell graph (ORI_WEIGHT_ADJ_MAT + EDGE_WEIGHT2COLOR_MAP = COLOR_ADJ_MAT)
COLOR_ADJ_MAT = []
# 表示顶点的颜色数量，初始值为 0。
VERTEX_COLOR_NO = 0
# 表示原始单元图中的顶点数量，初始值为 0。
CELLS_NUM = 0
# 表示扩展单元图中的顶点数量，初始值为 0
EXT_CELLS_NUM = 0
# 是一个字典，用于将顶点的权重映射到顶点的颜色上。初始值为空字典 {}。
VERTEX_WEIGHT2COLOR_MAP: dict[any, int] = {}
# 是邻接字典，用于存储图的邻接关系，初始值为空字典 {}。
ADJACENCY_DICT = {}
# 用于缓存 pynauty.certificate 的调用结果，以减少重复调用，初始值为空字典 {}。
CACHE_FOR_NAUTY = {}
# 是一个布尔变量，表示是否启用同构图的缓存，初始值为 True。
ENABLE_ISOMORPHISM = True
#  是一个字典，用于将因子映射到其索引，初始值为空字典 {}。
FACTOR2INDEX_MAP = {}
#  表示因子 0 的索引，初始值为 -1，用于表示尚未定义。
ZERO_FACTOR_INDEX = -1
# 是边权重的因子矩阵，用于表示边的因子信息，初始值为空列表 []。
FACTOR_ADJ_MAT = []

def update_factor_dict(factor):
    global FACTOR2INDEX_MAP, ZERO_FACTOR_INDEX
    if FACTOR2INDEX_MAP.get(factor) is None:
        FACTOR2INDEX_MAP[factor] = len(FACTOR2INDEX_MAP)
        if factor == 0:
            ZERO_FACTOR_INDEX = FACTOR2INDEX_MAP[factor]
def prime_init_factors(cell_weights, edge_weights):
    '''
    prime init factors for the cell weights and edge weights (including expression with symbols)
    all factors are stored in FACTOR_DICT
    '''
    for w in cell_weights:
        factored_list = factor_list(w)
        coef = factored_list[0]
        syms = factored_list[1]
        for k,_ in factorint(coef).items():
            update_factor_dict(k)
        for sym in syms:
            update_factor_dict(sym[0])
    for rs in edge_weights:
        for r in rs:
            factored_list = factor_list(r)
            coef = factored_list[0]
            syms = factored_list[1]
            for k,_ in factorint(coef).items():
                update_factor_dict(k)
            for sym in syms:
                update_factor_dict(sym[0])

def get_init_factor_set(cell_weights, edge_weights):
    cell_factor_set = []
    for w in cell_weights:
        vector = [0] * len(FACTOR2INDEX_MAP)
        factored_list = factor_list(w)
        coef = factored_list[0]
        syms = factored_list[1]
        for k,v in factorint(coef).items():
            vector[FACTOR2INDEX_MAP[k]] = v
        for sym in syms:
            vector[FACTOR2INDEX_MAP[sym[0]]] = int(sym[1])
        cell_factor_set.append(tuple(vector))
    global FACTOR_ADJ_MAT
    for i in range(len(edge_weights)):
        rs = edge_weights[i]
        vecs = []
        for j in range(len(rs)):
            r = rs[j]
            vector = [0] * len(FACTOR2INDEX_MAP)
            factored_list = factor_list(r)
            coef = factored_list[0]
            syms = factored_list[1]
            for k,v in factorint(coef).items():
                vector[FACTOR2INDEX_MAP[k]] = v
            for sym in syms:
                vector[FACTOR2INDEX_MAP[sym[0]]] = int(sym[1])
            vecs.append(tuple(vector))
        FACTOR_ADJ_MAT.append(vecs)
    return cell_factor_set

# 定义一个名为 get_vertex_color 的函数，接受一个参数 weight，用于处理顶点的颜色分配。
def get_vertex_color(weight):
    global VERTEX_COLOR_NO # 使用 global 关键字声明全局变量 VERTEX_COLOR_NO，该变量用于记录当前分配的顶点颜色编号。
    if weight not in VERTEX_WEIGHT2COLOR_MAP: # 检查 weight 是否不在 VERTEX_WEIGHT2COLOR_MAP 字典中，VERTEX_WEIGHT2COLOR_MAP 是一个将顶点权重映射到颜色编号的字典。
        VERTEX_WEIGHT2COLOR_MAP[weight] = VERTEX_COLOR_NO # 如果权重 weight 不在 VERTEX_WEIGHT2COLOR_MAP 中，将当前的颜色编号 VERTEX_COLOR_NO 分配给该权重，并存储在字典中。
        VERTEX_COLOR_NO += 1 # 将 VERTEX_COLOR_NO 的值加一，为下一个未分配的权重准备新的颜色编号。
    return VERTEX_WEIGHT2COLOR_MAP[weight] # 返回 VERTEX_WEIGHT2COLOR_MAP 字典中与 weight 对应的颜色编号。


def edgeWeight_To_edgeColor(): # 定义一个名为 edgeWeight_To_edgeColor 的函数，用于处理边权重到颜色的映射。
    global EDGE_COLOR_NUM, EDGE_WEIGHT2COLOR_MAP, COLOR_ADJ_MAT # 这些变量用于边的颜色编号、权重到颜色的映射，以及颜色邻接矩阵。
    EDGE_COLOR_NUM = 0 # 初始化边的颜色数量为 0。
    EDGE_WEIGHT2COLOR_MAP = {1:0}  # 初始化边权重到颜色的映射，默认为权重 1 映射到颜色 0。
    COLOR_ADJ_MAT = [] # 用于存储转换后的颜色邻接矩阵。
    # 开始遍历 ORI_WEIGHT_ADJ_MAT 中的每个列表（lst），ORI_WEIGHT_ADJ_MAT 是原始的权重邻接矩阵，用于表示图中各边的权重。
    for lst in ORI_WEIGHT_ADJ_MAT:
        tmp_list = [] # 初始化一个名为 tmp_list 的空列表，临时存储当前行中的边颜色信息。
        for w in lst: # 遍历当前列表 lst 中的每个元素 w，每个元素 w 表示一条边的权重。
            if w not in EDGE_WEIGHT2COLOR_MAP: # 检查当前权重 w 是否未在 EDGE_WEIGHT2COLOR_MAP 字典中，即是否尚未分配颜色编号。
                EDGE_COLOR_NUM += 1 # 如果当前权重 w 未被分配颜色编号，将 EDGE_COLOR_NUM 的值加一，为该权重分配一个新的颜色编号。
                EDGE_WEIGHT2COLOR_MAP[w] = EDGE_COLOR_NUM # 将新增的颜色编号 EDGE_COLOR_NUM 赋给当前权重 w，并将该映射存储到 EDGE_WEIGHT2COLOR_MAP 字典中。
            tmp_list.append(EDGE_WEIGHT2COLOR_MAP[w]) # 将 EDGE_WEIGHT2COLOR_MAP 中与权重 w 对应的颜色编号添加到 tmp_list 列表中。
        COLOR_ADJ_MAT.append(tmp_list) # 将 tmp_list 列表（表示当前行的边颜色信息）添加到 COLOR_ADJ_MAT 中，构建最终的颜色邻接矩阵。
    
    global LAYERS_NUM_FOR_CONVERT, EXT_CELLS_NUM # 声明全局变量，LAYERS_NUM_FOR_CONVERT 用于存储转换时所需的层数，EXT_CELLS_NUM 用于存储扩展单元的数量。
    # 通过计算以2为底的对数，确定将边着色图转换为特定层数的所需层数。LAYERS_NUM_FOR_CONVERT 使用 EDGE_COLOR_NUM+1 来计算出最低的层数，使其适应所有可能的边颜色组合。
    # LAYERS_NUM_FOR_CONVERT 通过确定层数，将多色边转换为每层内单一色的多层图。这种转换能使复杂的边着色图简化为每层内部只有单色边的图，从而方便进行进一步的图算法计算或分析。
    LAYERS_NUM_FOR_CONVERT = math.ceil(math.log2(EDGE_COLOR_NUM+1))
    EXT_CELLS_NUM = CELLS_NUM * LAYERS_NUM_FOR_CONVERT # EXT_CELLS_NUM 是扩展单元图中的顶点数。它通过将原始图的顶点数 CELLS_NUM 乘以所需的层数 LAYERS_NUM_FOR_CONVERT 计算得出。

def cellWeight_To_vertexColor(cell_weights):  # cell_weights 是一个列表，其中每个元素代表一个单元的权重。
    """
    # 它的主要功能是将单元的权重转换为对应的顶点颜色。这个过程可以帮助我们根据不同的权重值区分顶点的颜色，并统计各个颜色的数量。这个函数的主要目的是：
    # 将单元的权重转换为顶点颜色。
    # 统计每种颜色的数量。
    # 返回包含颜色信息的三个数据结构。
    """
    vertex_colors = [] # 创建一个空列表 vertex_colors，用于存储每个单元权重对应的顶点颜色。
    for w in cell_weights: # 遍历输入的 cell_weights 列表中的每个权重 w。
        vertex_colors.append(get_vertex_color(w)) # 对每个权重，调用 get_vertex_color(w) 函数，获取该权重对应的顶点颜色，并将其添加到 vertex_colors 列表中。
    color_dict = Counter(vertex_colors)  # 使用 Counter 统计 vertex_colors 列表中每种颜色的出现次数，结果存储在 color_dict 中，color_dict 是一个字典，键为颜色，值为该颜色的计数。
    color_kind = tuple(sorted(color_dict)) # color_kind 是一个元组，包含所有不同的颜色（按排序后）。
    color_count = tuple(color_dict[num] for num in color_kind) # color_count 是一个元组，包含对应于 color_kind 中每种颜色的计数。
    return vertex_colors, color_kind, color_count # vertex_colors：表示每个单元权重对应的顶点颜色的列表。# color_kind：不同颜色的元组。# color_count：每种颜色出现次数的元组。

def calculate_adjacency_dict():
    """
    calculate_adjacency_dict 函数通过生成邻接字典，将边着色图转换为多层图，并确保原始图的结构在扩展图中得到保留。
    这一过程在扩展图的每层中建立了单色边连接，并通过层内完全子图的方式保留顶点之间的“垂直”连接关系。
    这种多层邻接结构为进一步的图分析或同构性判断奠定了基础。
    """
    # Generate new edges
    adjacency_dict = {} # 初始化一个空字典 adjacency_dict，用于存储图的邻接关系。
    for i in range(EXT_CELLS_NUM): # 为每个扩展单元 i 初始化一个空列表，表示该单元的邻接节点。
        adjacency_dict[i] = []
    
    c2layers = {} # 初始化一个空字典 c2layers，用于存储颜色到层的映射。
    for k in range(EDGE_COLOR_NUM+1): # 遍历所有边的颜色（包括颜色数量的一个额外值，通常是 0），以构建颜色到层的映射。
        bi = bin(k) # 将当前颜色 k 转换为二进制字符串表示。
        bi = bi[2:] # 去掉二进制字符串的前缀 0b。
        bi = bi[::-1]  # 反转二进制字符串，以便从低位到高位进行遍历。
        layers = [i for i in range(len(bi)) if bi[i] == '1'] # 生成一个列表 layers，包含在反转二进制字符串中所有为 1 的位对应的索引。这些索引表示颜色 k 对应的层。
        c2layers[k] = layers # 将颜色 k 及其对应的层列表存入字典 c2layers 中。
    
    for i in range(CELLS_NUM): # 遍历原始单元的每一个节点 i
        for j in range(CELLS_NUM): # 对于每个节点 i，遍历其他每一个节点 j
            layers = c2layers[COLOR_ADJ_MAT[i][j]] # 获取节点 i 和节点 j 之间的颜色对应的层列表。
            for l in layers: # 遍历这些层。
                #  将扩展单元 l * CELLS_NUM + i（表示节点 i 在层 l 中的扩展单元）与 l * CELLS_NUM + j（表示节点 j 在层 l 中的扩展单元）之间的邻接关系添加到邻接字典中。
                adjacency_dict[l*CELLS_NUM+i].append(l*CELLS_NUM+j)
    
    # The vertical threads (each corresponding to one vertex of the original graph) 
    # can be connected using either paths or cliques.将每个原始图的顶点连接到对应的垂直层。
    for i in range(CELLS_NUM): # 再次遍历原始单元的每一个节点 i。
        clique = [i + j*CELLS_NUM for j in range(LAYERS_NUM_FOR_CONVERT)] # 为节点 i 生成一个完整的子图（clique），其中每个节点 j 的扩展单元是 i + j * CELLS_NUM，j 在层数范围内。
        for ii in clique: # 遍历当前的 clique 中的每个节点 ii。
            for jj in clique:  # 对于 clique 中的每个节点 jj，建立它们之间的邻接关系。
                if ii == jj:  # 如果 ii 和 jj 是同一个节点，则跳过，不连接自身。
                    continue
                adjacency_dict[ii].append(jj) # 将节点 jj 添加到节点 ii 的邻接列表中，表示它们之间有直接连接。
    global ADJACENCY_DICT # 声明全局变量 ADJACENCY_DICT，用于存储计算得到的邻接字典。
    ADJACENCY_DICT = adjacency_dict # 将计算得到的邻接字典赋值给全局变量 ADJACENCY_DICT。
# create_graph 函数的功能是创建一个无向图对象，并将其存储在全局变量 GRAPH 中，图的结构由扩展单元格的数量和邻接字典定义。这个图可以在后续的操作中用于图的算法或分析。
def create_graph():
    global GRAPH
    GRAPH = pynauty.Graph(EXT_CELLS_NUM,  # 这是图中的顶点数量
                          directed=False, # 表示创建一个无向图。
                          adjacency_dict=ADJACENCY_DICT) # 使用已定义的邻接字典来表示图的连接关系。

def adjust_vertex_coloring(colored_vertices):# 接受一个参数 colored_vertices，它是一个列表，表示每个顶点的颜色编号。
    '''
    文档字符串，说明该函数的功能是调整顶点的颜色编号，使其从 0 开始并且连续。
    
    Args:
        colored_vertices: list[int]
            The color no. of vertices.
            
        Returns:
            new_colored_vertices: list[int]
                The adjusted color no. of vertices.
            num_color: int
                The number of colors. 
    
    Example:
        colored_vertices = [7, 5, 7, 3, 5, 7]
        new_colored_vertices, num_color = adjust_vertex_coloring(colored_vertices)
        print(new_colored_vertices)  # [0, 1, 0, 2, 1, 0]
        print(num_color)  # 3
    '''
    color_map = {} # 初始化一个空字典 color_map，用于映射原始颜色编号到新的连续颜色编号。
    new_colored_vertices = [] # 初始化一个空列表 new_colored_vertices，用于存储调整后的颜色编号。
    num_color = 0 # 初始化一个变量 num_color 为 0，用于记录新的颜色数量。
    for c in colored_vertices: # 开始遍历输入的 colored_vertices 列表中的每个颜色编号 c。
        if c not in color_map: # 检查当前颜色 c 是否已经存在于 color_map 字典中。
            color_map[c] = num_color # 如果当前颜色 c 不在字典中，将其映射到当前的 num_color 值。
            num_color += 1 # 将 num_color 增加 1，表示发现了一种新颜色。
        new_colored_vertices.append(color_map[c]) # 将 color_map 中当前颜色 c 对应的新颜色编号添加到 new_colored_vertices 列表中。
    return new_colored_vertices, num_color # 返回调整后的颜色列表和新的颜色数量

def extend_vertex_coloring(colored_vertices, no_color):
    '''
    文档字符串，说明该函数的功能是扩展顶点集以转换带颜色的边。
    # 接受两个参数：colored_vertices（顶点的颜色编号列表）和 no_color（颜色的数量）。
    Args:
        colored_vertices: list[int]
            The color no. of vertices.
        no_color: int
            The number of colors.
    
    Returns:
        vertex_coloring: list[set[int]]
            The color set of vertices.
    
    Example:
        colored_vertices = [0, 1, 0, 2, 1, 0]
        no_color = 3
        vertex_coloring = extend_vertex_coloring(colored_vertices, no_color)
        print(vertex_coloring)  # [{0, 2, 5}, {1, 4}, {3}]
    '''
    # Extend the vertex set to convert colored edge
    ext_colored_vertices = [] # 初始化一个空列表 ext_colored_vertices，用于存储扩展后的颜色编号。
    for i in range(LAYERS_NUM_FOR_CONVERT):  # 开始一个循环，遍历层数 LAYERS_NUM_FOR_CONVERT，为每一层扩展颜色。
        ext_colored_vertices += [x + no_color * i for x in colored_vertices] # 将每个顶点的颜色编号 x 加上 no_color * i，以为每一层的顶点生成不同的颜色编号，并将结果添加到 ext_colored_vertices 列表中。
    
    # Get color set of vertices
    no_color *= LAYERS_NUM_FOR_CONVERT # 将 no_color 乘以 LAYERS_NUM_FOR_CONVERT，以更新颜色数量，反映扩展后的颜色数。
    vertex_coloring = [ set() for _ in range(no_color)] # 创建一个包含 no_color 个空集合的列表 vertex_coloring，用于存储每个颜色对应的顶点集合。
    for i in range(len(ext_colored_vertices)):  # 遍历扩展后的颜色列表 ext_colored_vertices 的索引。
        c = ext_colored_vertices[i] # 获取扩展颜色列表中索引为 i 的颜色 c。
        vertex_coloring[c].add(i)# 将索引 i 添加到对应颜色 c 的集合中，以记录该颜色的顶点。
    
    return vertex_coloring # 返回包含每种颜色对应顶点集合的列表 vertex_coloring。

def update_graph(colored_vertices):
    # for speed up, we have modified the function 'set_vertex_coloring' in graph.py of pynauty
    GRAPH.set_vertex_coloring(colored_vertices)
    return GRAPH

# not used
def get_gcd(cell_weights, cell_factor_tuple_list):
    gcd = Rational(1,1)
    gcd_vector = [0 for _ in range(len(FACTOR2INDEX_MAP))]
    for i in range(len(FACTOR2INDEX_MAP)):
        gcd_vector[i] = min([cell_factor_tuple_list[j][i] for j in range(CELLS_NUM)])
    if sum(gcd_vector) > 0:
        for i in range(len(cell_factor_tuple_list)):
            l_new = [x-y for x, y in zip(cell_factor_tuple_list[i], gcd_vector)]
            cell_factor_tuple_list[i] = tuple(l_new)
        for k,v in FACTOR2INDEX_MAP.items():
            if gcd_vector[v] > 0:
                gcd = expand(gcd * k**gcd_vector[v])
        for j in range(CELLS_NUM):
            cell_weights[j] = expand(cell_weights[j] / gcd)
    return gcd, cell_weights, cell_factor_tuple_list

# if all weights are integers, we can use this function to speed up by using cardinalities of fasctors
def dfs_wfomc(cell_weights, domain_size, cell_factor_tuple_list):  
    res = 0
    for l in range(CELLS_NUM):
        w_l = cell_weights[l]
        new_cell_weights = [cell_weights[i] * ORI_WEIGHT_ADJ_MAT[l][i] for i in range(CELLS_NUM)]
        if domain_size - 1 == 1:
            value = sum(new_cell_weights)   
        else:
            new_cell_factor_tuple_list = []
            for i in range(CELLS_NUM):
                new_cell_factor_tuple_list.append([x+y for x, y in zip(cell_factor_tuple_list[i], FACTOR_ADJ_MAT[l][i])])
                if ZERO_FACTOR_INDEX >= 0 and new_cell_factor_tuple_list[-1][ZERO_FACTOR_INDEX] > 0:
                    new_cell_factor_tuple_list[-1][ZERO_FACTOR_INDEX] = 1
                new_cell_factor_tuple_list[-1] = tuple(new_cell_factor_tuple_list[-1])
            # gcd, new_cell_weights, new_cell_factor_tuple_list = get_gcd(new_cell_weights, new_cell_factor_tuple_list)
            original_vertex_colors, vertex_color_kind, vertex_color_count = cellWeight_To_vertexColor(new_cell_factor_tuple_list) # convert cell weights to vertex colors
            if ENABLE_ISOMORPHISM:
                adjust_vertex_colors, no_color = adjust_vertex_coloring(original_vertex_colors) # adjust the color no. of vertices to make them start from 0 and be continuous
                if tuple(adjust_vertex_colors) not in CACHE_FOR_NAUTY:
                    can_label = pynauty.certificate(update_graph(extend_vertex_coloring(adjust_vertex_colors, no_color)))
                    CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)] = can_label
                else:
                    can_label = CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)]
            else:
                can_label = tuple(original_vertex_colors)
            
            value = IG_CACHE.get(domain_size-1, vertex_color_kind, vertex_color_count, can_label)
            if value is None:
                value = dfs_wfomc(new_cell_weights, domain_size - 1, new_cell_factor_tuple_list)
                value = expand(value)
                IG_CACHE.set(domain_size-1, vertex_color_kind, vertex_color_count, can_label, value)
        res += w_l * value # * expand(gcd**(domain_size - 1))
    return res


# dfs_wfomc_real 是一个递归算法，表示在加权模型计数 (WFOMC) 问题上进行深度优先搜索的实际递归函数。
# ComputeCCG(G, h)
# 它遍历领域中的每个元素，对其权重进行修改，并通过图同构来优化计算。
# 该函数是加权一阶模型计数（WFOMC）递归求解的真实版本。它接受三个参数：
# cell_weights: 当前图的单元格权重。[1, -1, -1, 1, x0*x1]
# domain_size: 当前领域的大小，表示常量集合中的元素个数。10
# treeNode: 递归过程中使用的树节点，用于存储中间结果，类型为 TreeNode。
def dfs_wfomc_real(cell_weights, domain_size, node: TreeNode = None):
    res = 0
    for l in range(CELLS_NUM): # 代码遍历一个范围为 CELLS_NUM 的循环，表示程序正处理多个“细胞”。
        w_l = cell_weights[l] # 获取当前细胞的权重。
        # 计算新的细胞权重 new_cell_weights。具体方法是遍历 CELLS_NUM 内的每个 i，将 cell_weights[i] 与 ORI_WEIGHT_ADJ_MAT[l][i] 相乘后
        # 调用 expand 函数来调整和扩展权重。这生成了一个新的细胞权重列表。
        new_cell_weights = [expand(cell_weights[i] * ORI_WEIGHT_ADJ_MAT[l][i]) for i in range(CELLS_NUM)]
        if PRINT_TREE: # 检查是否启用了 PRINT_TREE 标志。如果启用，程序将创建一个新的 TreeNode 对象，传递 new_cell_weights 和 node.depth+1 作为参数。这表示树的深度递增一层，并将其添加到 node.cell_to_children[l] 中。
            node.cell_to_children[l] = TreeNode(new_cell_weights, node.depth+1)
        if domain_size - 1 == 1:
            value = sum(new_cell_weights)   
        else: # 也就是构建nauty的graph # 调用 cellWeight_To_vertexColor 函数将 new_cell_weights 转换为图结构中的顶点颜色。cellWeight_To_vertexColor 返回的值包含原始顶点颜色、颜色种类和颜色计数。
            original_vertex_colors, vertex_color_kind, vertex_color_count = cellWeight_To_vertexColor(new_cell_weights) # convert cell weights to vertex colors
            if ENABLE_ISOMORPHISM: # 如果启用了 ENABLE_ISOMORPHISM（同构处理），代码会调整顶点颜色，使颜色编号从 0 开始并连续。
                adjust_vertex_colors, no_color = adjust_vertex_coloring(original_vertex_colors) #adjust_vertex_coloring() 函数用于执行颜色调整。返回颜色列表和新的颜色数量
                if tuple(adjust_vertex_colors) not in CACHE_FOR_NAUTY: # 判断 adjust_vertex_colors 的元组形式是否已经存在于 CACHE_FOR_NAUTY 中。# CACHE_FOR_NAUTY 是一个缓存，用来保存之前已经生成的图的特征，以避免重复计算 # 这个判断是为了确保对于一个已经处理过的颜色分配，我们不需要重复调用 pynauty 进行图同构检查。
                    # 使用certificate()函数生成图的同构证书，传入的参数是调用 update_graph() 和 extend_vertex_coloring() 函数处理过的顶点颜色数据。
                    # 首先，调用 extend_vertex_coloring(adjust_vertex_colors, no_color) 扩展顶点颜色，以便图转换过程中可以准确反映每个顶点的颜色信息。
                    # extend_vertex_coloring 会根据层数将顶点颜色扩展为多层结构，使每层的顶点颜色信息可以反映在多层图中。
                    # update_graph 将更新后的顶点颜色信息设置到图中（GRAPH），使其可以用于图同构检查。
                    can_label = pynauty.certificate(update_graph(extend_vertex_coloring(adjust_vertex_colors, no_color)))
                    CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)] = can_label # 将计算出的 can_label 存储到缓存 CACHE_FOR_NAUTY 中，并以 adjust_vertex_colors 的元组形式为键，以避免未来对相同的顶点颜色组合进行重复计算。
                else: # 如果调整后的顶点颜色已经在缓存中，则直接从缓存 CACHE_FOR_NAUTY 中读取对应的证书，避免重新计算。
                    can_label = CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)]
            else: # 如果没有启用同构处理（即 ENABLE_ISOMORPHISM 为 False）
                can_label = tuple(original_vertex_colors) # 直接将原始的顶点颜色作为标签 can_label，跳过图同构相关的处理。
            # 这里通过 IG_CACHE.get 从缓存中获取当前 domain_size-1 和顶点颜色（vertex_color_kind, vertex_color_count）对应的 can_label 值。
            value = IG_CACHE.get(domain_size-1, vertex_color_kind, vertex_color_count, can_label) # 如果缓存中已有该组合的计算结果，则会直接返回值 value。
            if value is None: # 如果 value 为 None（即缓存中没有找到对应的值），程序会递归调用 dfs_wfomc_real 函数，
                value = dfs_wfomc_real(new_cell_weights, domain_size - 1, node.cell_to_children[l] if PRINT_TREE else None) # 传入新的细胞权重 new_cell_weights 和 domain_size - 1，同时检查是否启用了 PRINT_TREE 来决定是否传递 node.cell_to_children[l]。
                value = expand(value) # 调用 expand() 函数扩展 value，然后将扩展后的值存储到 IG_CACHE 缓存中，以便将来快速查找和使用。
                IG_CACHE.set(domain_size-1, vertex_color_kind, vertex_color_count, can_label, value) # 将 value 存入 IG_CACHE 缓存中，索引由 domain_size-1、vertex_color_kind、vertex_color_count 和 can_label 组合而成。
        # 将当前细胞的权重 w_l 与计算得到的 value 相乘，并将结果累加到 res 中。
        res += w_l * value # * expand(gcd**(domain_size - 1)) # 注释中的部分可能表明程序曾考虑过对结果进行更复杂的扩展处理（例如乘以某个最大公约数的幂），但目前这一部分被注释掉了
    return res # 最后，函数返回累加的结果 res，这是深度优先搜索和权重乘积计算的结果。

def get_cache_size():
    total_size = 0
    for n_level in IG_CACHE.cache: # k0 is domain size
        for k1,v1 in n_level.items(): # k1 is color kind
            for k2,v2 in v1.items(): # k2 is color count
                for k3,v3 in v2.items(): # k3 is can_label
                    total_size += 1
    return total_size

def clean_global_variables():
    
    global PRINT_TREE, ROOT, IG_CACHE, ORI_WEIGHT_ADJ_MAT, \
        LAYERS_NUM_FOR_CONVERT, EDGE_COLOR_NUM, EDGE_WEIGHT2COLOR_MAP, \
            COLOR_ADJ_MAT, VERTEX_COLOR_NO, CELLS_NUM, EXT_CELLS_NUM, \
                VERTEX_WEIGHT2COLOR_MAP, ADJACENCY_DICT, CACHE_FOR_NAUTY, \
                    ENABLE_ISOMORPHISM, FACTOR2INDEX_MAP, ZERO_FACTOR_INDEX, FACTOR_ADJ_MAT
    
    PRINT_TREE = False
    ROOT = TreeNode([], 0)
    IG_CACHE = IsomorphicGraphCache()
    ORI_WEIGHT_ADJ_MAT = []
    LAYERS_NUM_FOR_CONVERT = 0
    EDGE_COLOR_NUM = 0
    EDGE_WEIGHT2COLOR_MAP = {1:0}
    COLOR_ADJ_MAT = []
    VERTEX_COLOR_NO = 0
    CELLS_NUM = 0
    EXT_CELLS_NUM = 0
    VERTEX_WEIGHT2COLOR_MAP = {}
    ADJACENCY_DICT = {}
    CACHE_FOR_NAUTY = {}
    ENABLE_ISOMORPHISM = True
    FACTOR2INDEX_MAP = {}
    ZERO_FACTOR_INDEX = -1
    FACTOR_ADJ_MAT = []


# recursive_wfomc 是这个算法的核心函数，用于根据给定的公式、领域、权重和谓词计算加权模型计数（WFOMC）。
def recursive_wfomc(formula: QFFormula, # # 要求解的逻辑公式。skolem之后的
                  domain: set[Const],
                  get_weight: Callable[[Pred], tuple[RingElement, RingElement]], #一个正 一个负
                  leq_pred: Pred, # leq_pred: 谓词，用于指定“≤”关系。线性阶，这个本质上是一个binary predicate
                  ues_dft: bool = False, # TODO 用一个参数接收
                  real_version: bool = True) -> RingElement: # real_version: 是否使用真实版本的递归。权重为整数的时候，可以做优化，把它因式分解。为True，表示不做因式分解，只有特殊情况为False，才因式分解。因式分解没有写在文章里面
    domain_size = len(domain)
    res = Rational(0, 1)
    # 这行代码使用 for 循环迭代 build_cell_graphs 函数的返回值。每次迭代时，从返回的元组中解包出 cell_graph 和 weight。
    # build_cell_graphs 函数接收三个参数：formula、get_weight 和 leq_pred，并返回一系列的单元图及其对应的权重。
    for cell_graph, weight in build_cell_graphs( # 这里假设就build出一个cell，因为有的句子里面零元谓词，要先处理零元谓词 ，比如\forall: X:(P(X)|Q()), if Q() = T, T else \forall X:P(X)
        formula, get_weight, leq_pred=leq_pred, use_dft = ues_dft # 所以在这个函数里面分情况讨论，遍历所有的为T还是F，得到化简后的句子，然后再cell_graph
    ): # 从这个下面就是文章里面的主要构成 # TODO
        cell_weights = cell_graph.get_all_weights()[0] # 获取顶点权重a
        edge_weights = cell_graph.get_all_weights()[1] # 获取边权重b
        
        clean_global_variables() # 文章算法里面有很多全局变量，比如cache等等，
        
        IG_CACHE.init(domain_size) # 然后初始化同构图缓存IG_CACHE。
        global ORI_WEIGHT_ADJ_MAT, CELLS_NUM

        # not change in the same problem
        CELLS_NUM = len(cell_weights)  # 细胞数量
        ORI_WEIGHT_ADJ_MAT = edge_weights # 原始权重邻接矩阵，用于存储图的边的权重信息。
        edgeWeight_To_edgeColor() # 用于处理边权重到颜色的映射 # 需要了解一下nauty大概的思路，第四个有一个演示
        calculate_adjacency_dict() # 构建多层图的邻接字典 ADJACENCY_DICT，该字典存储了扩展图中每个顶点的邻接关系。
        create_graph() # 创建一个无向图对象
        
        # disable isomorphism
        # global ENABLE_ISOMORPHISM
        # ENABLE_ISOMORPHISM = False
        
        if not real_version: # 如果 real_version 是 False，则调用非真实版本的递归 dfs_wfomc。这通常是一个优化版本，避免一些计算。
            prime_init_factors(cell_weights, edge_weights) # prime_init_factors：初始化因式分解。
            cell_factor_tuple_list = get_init_factor_set(cell_weights, edge_weights) # get_init_factor_set：获取因式分解后的初始因子集，并将其传递给 dfs_wfomc 进行递归计算。
            res_ = dfs_wfomc(cell_weights, domain_size, cell_factor_tuple_list)
        else: # 看这里，上面的不要看 # 如果 real_version 是 True，则调用真实版本的递归 dfs_wfomc_real，并将结果返回。
            global ROOT
            ROOT.cell_weights = cell_weights
            res_ = dfs_wfomc_real(cell_weights, domain_size, ROOT) # 这里就是算法的第五行 类似于一个树的深度优先搜索，其实就是递归，
            if PRINT_TREE:
                print_tree(ROOT) 
        res = res + weight * res_ # TODO 这里有bug
        print(weight * res_)
    return res