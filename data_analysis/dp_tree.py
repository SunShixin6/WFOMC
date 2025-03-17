import csv
import sys

sys.path.append("/root/pycharm_workspace/WFOMC")  # 确保 wfomc 模块可以被找到

import pickle
import pydot
import graphviz
from collections import deque
import pandas as pd

node_dict = {}  # 临时缓存之前的节点
domain_size = 5


class TreeNode:

    def __init__(self, d_new, k_new, w_new, level, parent=None):
        self.d_new = d_new
        self.k_new = k_new
        self.w_new = w_new
        self.level = level
        self.parent = parent
        self.children = []

    def add_child(self, child):
        # 添加子节点到当前节点的children列表中
        self.children.append(child)  # 添加子节点

    def __str__(self):
        return f"d_new={self.d_new}, k_new={self.k_new}, w_new={self.w_new}, level={self.level}"

    __repr__ = __str__


def save_to_file():
    """
    将树保存到文件中
    """
    # 创建一个虚拟的根节点
    virtual_root = TreeNode(d_new="root", k_new="root", w_new="root", level="root")  # 创建虚拟根节点

    # 收集所有最后一层的节点
    last_layer_nodes = [node for node in node_dict.values() if node.level == domain_size]  # 收集最后一层节点
    # 将所有最后一层的节点添加到虚拟根节点的子节点中
    for node in last_layer_nodes:
        virtual_root.add_child(node)  # 添加子节点

    # 将虚拟根节点作为最终的根节点
    root = virtual_root  # 设置根节点
    if root is None:
        print("错误：未找到根节点，请检查动态规划的逻辑")
    # else:
    # print(f"根节点: d_new={root.d_new}, k_new={root.k_new}, w_new={root.w_new}")
    with open("/root/pycharm_workspace/WFOMC/data_analysis/tree.pkl", 'wb') as f:
        pickle.dump(root, f)
        print("树已保存到文件中")


def save_to_svg():
    """
    读取保存的树并绘制成svg
    """
    with open("/root/pycharm_workspace/WFOMC/data_analysis/tree.pkl", "rb") as f:
        root = pickle.load(f)
        if not root:
            print("root为空")

    def visualize_tree(root, output_file='/root/pycharm_workspace/WFOMC/data_analysis/tree', dpi=1000):
        """
        使用 pydot 绘制树的可视化图形。
        :param root: 树的根节点
        :param output_file: 输出的图片文件名
        """
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

        # 使用 graphviz 渲染图形并设置 dpi 参数
        dot_data = graph.to_string()  # 获取图的字符串表示
        dot = graphviz.Digraph(format='svg')  # 创建一个SVG格式的图
        dot.body.extend([line for line in dot_data.split(' ') if line.strip() and not line.startswith('strict')])  # 添加图的主体内容
        dot.render(output_file, cleanup=True)  # 渲染并保存图
        print(f"树的可视化图形已保存到文件: {output_file}")  # 打印保存成功信息

    # 调用函数绘制树
    visualize_tree(root)


def decompose_tree():
    """
    将树转为字典
    """
    with open("/root/pycharm_workspace/WFOMC/data_analysis/data.pkl", "rb") as f:
        root = pickle.load(f)
        if not root:
            print("root为空")
            # break
        # else:
        #     print(root)
    all_dict = dict((n, list()) for n in range(domain_size + 1))
    # 层次遍历树
    queue = deque([root])
    # print("queue:", queue)
    while queue:
        node = queue.popleft()
        if node.level != "root":
            all_dict[node.level].append((node.d_new, node.k_new, node.w_new))
            # if node.d_new[0] > sum(node.k_new):
            #     print("error")
        # print("node:", node)
        # print("node.children:", node.children)

        for child in node.children:
            # print("child",child)
            queue.append(child)
    # print("all_dict",all_dict)
    pd.set_option('display.max_rows', None)  # 显示所有行

    # with open("/root/pycharm_workspace/WFOMC/data_analysis/tree_to_dict.csv", "w", newline="") as file:
    #     writer = csv.writer(file)
    #     for key, value in all_dict.items():
    #         writer.writerow([key, value])
    #
    for i in range(domain_size):
        df = pd.DataFrame(all_dict[i], columns=["ccs", "ivec(cell config)", "weight"])
        df.to_csv(f"/root/pycharm_workspace/WFOMC/data_analysis/before/before_{i}.csv", index=False)

    with open("/root/pycharm_workspace/WFOMC/data_analysis/tree_to_dict.pkl", "wb") as f:
        pickle.dump(all_dict, f)

    print("tree_to_dict树转为字典")


def data_process():
    """
    字典数据处理
    """
    with open("/root/pycharm_workspace/WFOMC/data_analysis/tree_to_dict.pkl", "rb") as f:
        all_dict = pickle.load(f)

        for subproblem_size in range(0, 5):
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
            df.to_csv(f"/root/pycharm_workspace/WFOMC/data_analysis/after/after_{subproblem_size}.csv", index=False)
        print("data_process成功")


if __name__ == '__main__':
    decompose_tree()
    data_process()
