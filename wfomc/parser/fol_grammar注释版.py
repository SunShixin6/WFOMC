function_free_logic_grammar = r"""
    # ffl（函数自由逻辑表达式）可以是三种形式之一：atomic_ffl（原子逻辑表达式）、compound_ffl（复合逻辑表达式）、exactlyone（ExactlyOne结构）。
    ?ffl: atomic_ffl | compound_ffl | exactlyone
     
    # atomic_ffl 定义了一个原子逻辑表达式，即一个谓词（predicate），后面可选地跟着括号包裹的一组项（terms），这是最基本的逻辑表达式形式。
    atomic_ffl: predicate [left_parenthesis terms right_parenthesis]

    # exactlyone 定义了一个特殊结构，表示谓词集合中“恰好有一个”谓词为真。这个表达式的形式是 "ExactlyOne[predicate1, predicate2, ...]"    
    exactlyone: "ExactlyOne" left_square_bracket predicates right_square_bracket

    # predicates 是一个谓词列表，包含一个或多个谓词，用逗号分隔。这意味着可以定义多个谓词来表达逻辑结构。    
    predicates: predicate ("," predicate)*

    # terms 是一个项列表，包含一个或多个项（常量或变量），也是用逗号分隔。    
    terms: term ("," term)*

    # negation 定义了逻辑的否定运算（即“非”），由符号~表示。它作用于一个函数自由逻辑表达式ffl。    
    negation: not ffl

    # conjunction 表示逻辑的合取（即“与”运算），由符号&表示，表示两个逻辑表达式ffl同时为真。    
    conjunction: ffl and ffl

    # disjunction 表示逻辑的析取（即“或”运算），由符号|表示，表示两个逻辑表达式ffl中至少一个为真。    
    disjunction: ffl or ffl

    # implication 表示逻辑的蕴涵（即“如果……那么……”），由符号->表示。如果第一个逻辑表达式ffl为真，则第二个逻辑表达式ffl也为真。    
    implication: ffl implies ffl

    # equivalence 表示逻辑的等价（即“当且仅当”），由符号<->表示。如果两个逻辑表达式ffl同时为真或同时为假，则它们是等价的。    
    equivalence: ffl iff ffl

    # compound_ffl 定义了复合逻辑表达式的几种形式：
    # 使用括号括起的逻辑表达式(ffl)。
    # 含有量词的表达式（如∀x (ffl)），表示逻辑量化。
    # 其他复合逻辑，如等价、蕴涵、析取、合取、否定。    
    ?compound_ffl: left_parenthesis ffl right_parenthesis -> parenthesis
       | quantifier_variable ":" left_parenthesis ffl right_parenthesis -> quantification
       | equivalence
       | implication
       | disjunction
       | conjunction
       | negation

    # term 表示逻辑中的项，可以是常量（constant）或变量（variable）。
    ?term: constant
        | variable

    # 定义了方括号和圆括号，表示逻辑表达式中使用的符号。
    left_square_bracket: "["
    right_square_bracket: "]"
    left_parenthesis: "("
    right_parenthesis: ")"

    # quantifier_variable 表示量词加变量的组合，例如 ∀x 或 ∃y。    
    quantifier_variable: quantifier variable

    # quantifier 定义了三种量词：
    # universal_quantifier：全称量词，∀，表示“对于所有”。
    # existential_quantifier：存在量词，∃，表示“存在某个”。
    # counting_quantifier：计数量词，表示“至少/至多n个”。
    ?quantifier: universal_quantifier | existential_quantifier | counting_quantifier

    # 分别定义全称量词∀和存在量词∃的符号表示。
    universal_quantifier: "\\forall"
    existential_quantifier: "\\exists"

    # 计数量词的表示形式，类似“存在至少n个”的结构，其中包含比较运算符（comparator）和计数参数（count_parameter）。    
    counting_quantifier: "\\exists_{" comparator count_parameter "}"

    # constant 是小写字母（LCASE_CNAME），表示常量。    
    constant: LCASE_CNAME

    # variable 是大写字母（UCASE_LETTER），表示变量。
    variable: UCASE_LETTER

    # predicate 是谓词，用于表达逻辑中的命题。    
    predicate: CNAME

    # 定义了逻辑运算符的符号：
    # ~：非。
    # &：与。
    # |：或。
    # ->：蕴涵。
    # <->：等价。    
    not: "~"
    and: "&"
    or: "|"
    implies: "->"
    iff: "<->"

    # count_parameter 是整数类型（INT），用于计数量词。
    count_parameter: INT

    # comparator 是比较运算符，定义了等于、不等于、大于等符号。    
    ?comparator: equality | le | ge | lt | gt | nequality
    equality: "="
    nequality: "!="
    le: "<="
    ge: ">="
    lt: "<"
    gt: ">"

    # 定义了小写字母开头的名称（LCASE_CNAME），可以包含下划线、小写字母或数字。    
    LCASE_CNAME: LCASE_LETTER ("_"|LCASE_LETTER|DIGIT)*

    # 导入了一些通用的定义，如小写字母、大写字母、常量名称、数字、浮点数、带符号的数字、注释等。
    %import common.LCASE_LETTER
    %import common.UCASE_LETTER
    %import common.CNAME
    %import common.DIGIT
    %import common.FLOAT
    %import common.INT
    %import common.SIGNED_NUMBER
    %import common.NUMBER
    %import common.WS
    %import common.SH_COMMENT

    # 忽略空白字符（WS）和 shell 风格的注释（SH_COMMENT）。
    %ignore WS
    %ignore SH_COMMENT
"""
