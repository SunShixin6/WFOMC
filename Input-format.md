## Input Format

Input files use the `.wfomcs` extension. 

Comments beginning with `#` are supported and are ignored by the parser.

### First-Order Sentence

#### Variables

A variable is a **single uppercase letter**: `X`, `Y`, `Z`, etc. Only two distinct variables may appear in any sentence (two-variable logic).

#### Constants

Constants must **start with a lowercase letter**, optionally followed by letters, digits, or underscores: `alice`, `bob_1`, `node2`. Constants appear only in ground literals (e.g., in unary evidence).

#### Predicates

Predicate names are identifiers that may start with an uppercase or lowercase letter: `E`, `R`, `Odd`, `blue`, `neighbor`. Predicate names are case-sensitive.

```text
E(X, Y)     # binary predicate
Odd(X)      # unary predicate
R(X)        # unary predicate
```

#### Linear Order Predicate `LEQ`

In this repository, `LEQ(X, Y)` is the conventional name for the binary predicate intended to denote a linear order on the domain. If the input uses `LEQ`, you should solve it with `incremental`, `recursive`, or `incremental3`. Example:

```text
\forall X: (\forall Y: ((H(Y) & LEQ(X,Y)) -> H(X)))
```

#### Logical Connectives

| Syntax | Meaning     | Example                                  |
| ------ | ----------- | ---------------------------------------- |
| `~`    | negation    | `~E(X, X)`                               |
| `&`    | conjunction | `R(X) & ~B(X)`                           |
| `|`    | disjunction | `R(X) | B(X)`                            |
| `->`   | implication | `E(X,Y) -> E(Y,X)`                       |
| `<->`  | iff         | `Odd(X) <-> \exists_{1mod2} Y: (E(X,Y))` |

To avoid ambiguity, we recommend using parentheses explicitly when combining `~`, `&`, `|`, `->`, and `<->` in the same formula.

Parentheses `(...)` are used for grouping. The quantifier body **must** be enclosed in parentheses:

```text
\forall X: (...)     # correct
\forall X: ...       # parse error — parentheses are required
```

Multiple top-level conjuncts are commonly written one per line, joined by `&`:

```text
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (\exists_{=2} Y: (E(X,Y)))
```

#### Quantifiers

| Syntax                     | Meaning                             |
| -------------------------- | ----------------------------------- |
| `\forall X: (...)`         | Universal quantification over `X`   |
| `\exists X: (...)`         | Existential quantification over `X` |
| `\exists_{=k} X: (...)`    | Exactly *k* witnesses               |
| `\exists_{<=k} X: (...)`   | At most *k* witnesses               |
| `\exists_{rmodk} X: (...)` | Number of witnesses ≡ *r* (mod *k*) |

Here, `r` and `k` must be non-negative integers in the input syntax. Semantically, `\exists_{rmodk}` means that the number of witnesses is congruent to `r` modulo `k`. In standard mathematical notation, one usually assumes `k > 0` and typically `0 <= r < k`.

**Examples:**

```text
# Every vertex has exactly 2 neighbors
\forall X: (\exists_{=2} Y: (E(X,Y)))

# Every vertex has at most 3 neighbors
\forall X: (\exists_{<=3} Y: (E(X,Y)))

# Every vertex has an even number of neighbors
\forall X: (\exists_{0mod2} Y: (E(X,Y)))

# Exactly 0 vertices satisfy Odd(X)
\exists_{=0} X: (Odd(X))
```

### Domain

Declare the domain on a single line as `name = value`. The name can be any identifier (`V`, `n`, `domain`, `person`, …).

**Integer size** — the domain is implicitly `{1, 2, ..., n}`:

```text
V = 10
n = 5
domain = 100
```

**Explicit constant set** — a comma-separated list of identifiers inside `{...}`:

```text
domain = {alice, bob, charlie}
```

For consistency, we recommend using lowercase names for explicit constants.

### Weights (optional)

Each weight line assigns a positive weight and a negative weight to one predicate:

```text
<positive_weight>  <negative_weight>  <PredicateName>
```

Weights are signed integers or floating-point numbers. Predicates not listed default to weight `1` for both true and false groundings.

```text
2.7 1 E       # true groundings of E get weight 2.7; false get weight 1
1 -1 C        # true groundings of C get weight 1; false get weight -1
0.5 0.5 R     # uniform weight 0.5
```

Negative weights are valid and arise naturally in inclusion-exclusion encodings (e.g., `m-odd-degree` uses `1 -1 C`).

### Cardinality Constraints (optional)

These constraints restrict the total number of true groundings of a predicate in the model.

```text
|P| = k
|P| >= k
|P| <= k
|P| > k
|P| < k
```

