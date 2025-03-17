class TreeNode:
    def __init__(self, k_new, d_new, w_new, parent=None):
        self.k_new = k_new
        self.d_new = d_new
        self.w_new = w_new  # 修改为w_new
        self.parent = parent
        self.level = 0 if parent is None else parent.level + 1

    def __str__(self):
        return f"d_new={self.d_new}, k_new={self.k_new}, w_new={self.w_new}, level={self.level}"  # 修改为w_new