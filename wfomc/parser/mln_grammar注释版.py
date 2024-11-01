# 这段代码定义了Markov逻辑网络（MLN）的语法，通过导入其他模块中的逻辑语法和约束条件，并结合领域定义和规则定义，描述了MLN的结构。下面逐行解释代码的含义：

# 这两行代码从其他模块中导入了函数自由逻辑语法（function_free_logic_grammar）和基数约束语法（cc_grammar），后者可能用于定义一些关于数量的约束条件。
from .fol_grammar import function_free_logic_grammar
from .cardinality_constraints_grammar import cc_grammar


# 这段代码定义了领域语法：
domain_grammar = r"""
    # domain 表示领域定义，格式为domain_name = domain_spec， 
    domain: domain_name "=" domain_spec
    
    # domain_name 是一个常量名称（CNAME），表示领域的名称。
    domain_name: CNAME
    
    # domain_spec 定义了领域的具体类型，可以是整数域（int_domain）或集合域（set_domain），其中集合域由一组元素构成。
    ?domain_spec: INT               -> int_domain
        | ("{" domain_elements "}") -> set_domain

    # domain_elements 是元素的列表，用逗号分隔。        
    domain_elements: element ("," element)*
    
    # element 是集合中的单个元素，也用常量名称表示。
    element: CNAME
"""

# 这段代码定义了规则语法：
rule_grammar = r"""
    # rules 表示多个规则的集合，规则可以是硬规则（hard_rule）或软规则（soft_rule）。
    rules: rule*
    rule: hard_rule | soft_rule

    # hard_rule 定义了硬规则，由函数自由逻辑表达式（ffl）和一个句号表示。这种规则必须严格满足    
    hard_rule: ffl "."

    # soft_rule 定义了软规则，包含一个**权重值（weighting）**和逻辑表达式。权重值是有符号的数字，表示规则的强度。
    soft_rule: weighting ffl
    weighting: SIGNED_NUMBER

    # function_free_logic_grammar 是之前导入的函数自由逻辑语法，扩展了规则的定义。    
""" + function_free_logic_grammar

# 这段代码定义了最终的MLN语法：
grammar = r"""

    # mln 表示一个完整的 Markov 逻辑网络，由规则（rules）、领域定义（domain）和基数约束（cardinality_constraints）组成。
    ?mln: rules domain cardinality_constraints
    
    # rule_grammar、domain_grammar 和 cc_grammar 分别是规则语法、领域语法和基数约束语法，它们组合在一起构成了 Markov 逻辑网络的完整语法定义。
""" + rule_grammar + domain_grammar + cc_grammar
