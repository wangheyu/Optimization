## 可行集的几何特性

现在让我们来考虑一下在标准形式下，线性规划可行集的几何特性，并引出正确的求解思路。首先我们约定(13.1)中的$A$是行满秩的。
$$
\begin{equation}
\min_{x \in \mathbb{R}^n} f(x), \quad \mbox{s.t.}~ Ax = b, x\geq0.
\tag{13.12}
\end{equation}
$$
由于$A$是$m \times n$矩阵，如果我们在$A$中确定$m$列，在指标$\{1, 2, \cdots, n\}$中选择$m$个指标构成子集$\mathcal{B}$，那么我们相当于从$A$中切割出一个方阵：
$$
B_{m \times m} = [A_i]_{i\in\mathcal{B}},
$$
$A_i$是$A$的第$i$列。若$B$非奇异，则显然可以唯一确定一个点
$$
x_B = B^{-1}b,
$$
这里$x_B$相当于只对$i \in \mathcal{B}$中的分量给了赋值，对于那些没有被$\mathcal{B}$选中的分量，我们记为
$$
\mathcal{N} = \{1, 2, \cdots, n\} \backslash \mathcal{B},
$$
它同样对应$A$的一些列，并可以组成一个矩阵
$$
N_{m \times (n - m)} = [A_i]_{i \in \mathcal{N}},
$$
以及对应的$x$分量组成的向量$x_N$。（$i \in \mathcal{N}$的$x$的分量组成，相当于现在$x$都可以分成$x_B$和$x_N$两部分。）对于$x_N$，我们令它的全部分量都为零，从而满足约束$x \geq 0$。但对于$x_B$，我们并不清楚它是否可行（是否各分量都大于等于零）。

（以上记号对$A$和$x$都按指标$\mathcal{B}$和$\mathcal{N}$做了一个划分，$\mathcal{B}\cup\mathcal{N} = \{1, 2, \cdots, n\}$。）

现在我们规定：

**定义** 称(13.1)的可行点$x$是一个基础可行点(basic feasible point)，若存在指标集$\{1, 2, \cdots, n\}$的子集$\mathcal{B}$满足：

1. $\mathcal{B}$包含$m$的指标；
2. $i \notin \mathcal{B}$，$x_i = 0$；
3. $B_{m \times m} = [A_i]_{i \in \mathcal{B}}$非奇异。

而$\mathcal{B}$称为(13.1)的一个基(basis)，$B$称为基矩阵。

（$\mathcal{N}$有时称为“非基”，$N$有时称为“非基矩阵”。）