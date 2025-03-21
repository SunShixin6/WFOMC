import json
import os.path
import sys
from pprint import pprint
import logging as logger

sys.path.append("/root/pycharm_workspace/WFOMC")  # 确保 wfomc 模块可以被找到

import pickle
import pydot  # 用于绘制树结构的图形。
import graphviz
from collections import deque  # 双端队列，用于层次遍历树。
import pandas as pd
import shutil

pd.set_option('display.max_rows', None)  # 显示所有行
# 全局变量
node_dict = {}  # 定义一个字典临时缓存之前变量过的节点
domain_size = None
cwd = None  # 当前脚本所在的目录。


def init_script_dir():  # 初始化 cwd，获取当前脚本所在的目录。
    global cwd
    cwd = os.path.dirname(os.path.abspath(__file__))  # 获取当前脚本所在的目录
    ensure_clean_directory(os.path.join(cwd, "after"))  # 确保目标目录是空的，如果存在文件则删除。
    ensure_clean_directory(os.path.join(cwd, "before"))
    logger.info("旧数据清除成功")


def init_domain_size():  # 初始化 domain_size
    global domain_size
    with open(os.path.join(cwd, "domain_size.txt"), "r") as f:
        domain_size = int(f.read())


def save_domain_size(domain_size):
    global cwd
    with open(os.path.join(cwd, "domain_size.txt"), "w") as f:
        f.write(str(domain_size))


def ensure_clean_directory(directory):
    """
    确保目标目录是空的，如果存在文件则删除。
    :param directory: 目标目录路径
    """
    if os.path.exists(directory):
        # 删除目录及其所有内容
        shutil.rmtree(directory)
    # 重新创建目录
    os.makedirs(directory)
    return directory


class TreeNode:  # NOTE: 创建类
    def __init__(self, d_new, k_new, w_new, level, parent=None):
        self.d_new = d_new  # 基数约束
        self.k_new = k_new  # cell_config
        self.w_new = w_new  # 权重
        self.level = level  # 也就是sum(cell_config) 表示有几个元素
        self.parent = parent
        self.children = []

    def add_child(self, child):
        self.children.append(child)  # 添加子节点到当前节点的children列表中

    def __str__(self):
        return f"d_new={self.d_new}, k_new={self.k_new}, w_new={self.w_new}, level={self.level}"

    __repr__ = __str__


def tree_to_file():
    ## 将树保存到文件中
    global cwd
    virtual_root = TreeNode(d_new="root", k_new="root", w_new="root", level="root")  ## 创建一个虚拟的根节点
    last_layer_nodes = [node for node in node_dict.values() if node.level == domain_size]  ## 收集所有最后一层的节点
    for node in last_layer_nodes:  # 将所有最后一层的节点添加到虚拟根节点的子节点中
        virtual_root.add_child(node)  # 添加子节点
    # 将虚拟根节点作为最终的根节点
    root = virtual_root  # 设置根节点
    if root is None:
        print("错误：未找到根节点，请检查动态规划的逻辑")
    with open(os.path.join(cwd, "tree.pkl"), 'wb') as f:
        pickle.dump(root, f)
        logger.info("树已保存到文件 data_analysis/tree.pkl中")


def tree_to_svg():
    ## 读取保存的树并绘制成svg
    global cwd
    with open(os.path.join(cwd, "tree.pkl"), "rb") as f:
        root = pickle.load(f)
        if not root:
            print("root为空")

    def visualize_tree(root, output_file=os.path.join(cwd, "tree"), dpi=1000):
        ## 使用 pydot 绘制树的可视化图形。
        graph = pydot.Dot(graph_type='digraph')  # 创建一个有向图

        def add_nodes_edges(node, parent_id=None):
            node_id = str(id(node))  # 生成节点ID
            label = f"Config: {node.k_new}\nCCS: {node.d_new[0]}\nW: {node.w_new}"  # 生成节点标签
            graph.add_node(pydot.Node(node_id, label=label))  # 添加节点到图中

            if parent_id is not None:
                graph.add_edge(pydot.Edge(parent_id, node_id))  # 添加边到图中

            for child in node.children:
                add_nodes_edges(child, node_id)  # 递归添加子节点

        add_nodes_edges(root)  # 调用函数添加节点和边
        dot_data = graph.to_string()  # 获取图的字符串表示
        dot = graphviz.Digraph(format='svg')  # 创建一个SVG格式的图
        dot.body.extend([line for line in dot_data.split(' ') if line.strip() and not line.startswith('strict')])  # 添加图的主体内容
        dot.render(output_file, cleanup=True)  # 渲染并保存图
        logger.info(f"树的可视化图形已保存到文件: data_analysis/{output_file}.svg")  # 打印保存成功信息

    visualize_tree(root)


def tree_to_dict():
    ## 将树转为字典
    global cwd  # 查找路径
    with open(os.path.join(cwd, "tree.pkl"), "rb") as f:  # 加载文件
        root = pickle.load(f)
        if not root:
            print("root为空")
            return None
    all_dict = dict((n, list()) for n in range(domain_size + 1))  ## 存储所有从根节点到叶子节点的dfs遍历结果，domain_size作为键

    def dfs(node):  # 深度优先遍历，将节点信息存储到字典中。
        if node.level != "root":  # 忽略虚拟根节点
            all_dict[node.level].append((node.d_new, node.k_new, node.w_new))
        for child in node.children:  # 递归遍历子节点
            dfs(child)

    dfs(root)  # 从根节点开始遍历

    ## 将all_dict字典保存为csv文件
    for level in range(domain_size + 1):  ## 根据domain_size的值，提取列表，转为df对象，存储到csv文件中
        df = pd.DataFrame(all_dict[level], columns=["ccs", "ivec(cell config)", "weight"])
        df.to_csv(os.path.join(cwd, "before", f"before_{level + 1}.csv"), index=False)  # 保存为csv文件
    ## 将整个all_dict字典存储
    with open(os.path.join(cwd, "tree_to_dict.pkl"), "wb") as f:
        pickle.dump(all_dict, f)
    print("树转为字典，存储在文件tree_to_dict.pkl中")
    # pprint(all_dict)  ## 格式化打印字典


def decompose_data():
    ## 字典数据处理
    global cwd
    with open(os.path.join(cwd, "tree_to_dict.pkl"), "rb") as f:  # 读取存储的字典tree_to_dict.pkl文件
        all_dict = pickle.load(f)

    for subproblem_size in range(0, domain_size + 1):
        subproblem_list = all_dict[subproblem_size]
        d_k_w = dict()
        for item in subproblem_list:  # ((0,), (1, 0, 0, 0, 0), 1),
            if item[0] not in d_k_w:  # d:(0,)
                k_w = dict()
            else:
                k_w = d_k_w[item[0]]

            if item[1] not in k_w:  # k:(1, 0, 0, 0, 0)
                k_w[item[1]] = item[2]  # w:1
            else:
                k_w[item[1]] += item[2]

            d_k_w[item[0]] = k_w

        df = []
        for d, k_w in d_k_w.items():
            for k, w in k_w.items():
                # if w != 0:
                df.append([d[0], k, w])

        df = pd.DataFrame(df, columns=["ccs", "ivec(cell config)", "weight"])
        df.to_csv(os.path.join(cwd, "after", f"after_{subproblem_size + 1}.csv"), index=False)
    print("数据处理成功成功，结果存放在after文件夹中")


if __name__ == '__main__':
    init_script_dir()  # 初始化 cwd，获取当前脚本所在的目录。
    init_domain_size()  # 从文件中读取 domain_size 是多少
    tree_to_dict()
    decompose_data()
