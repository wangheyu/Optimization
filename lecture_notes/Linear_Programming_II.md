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

这里有一个严重的问题是，这样的可行点的存在性。下面这个定理回答了这个问题：

**定理 13.2**

1. 若(13.1)的可行域非空，则它至少有一个基础可行点；
2. 若(13.1)有解，则至少一个解是基础可行点；
3. 若(13.1)可行且有界，则它必有解。

这里我们不再重复讲义上的证明，而是直接讨论一下，从基础可行点的几何意义究竟是什么。实际上，在
$$
A x = b
$$
中对$A$选取$m$列的操作我们是熟悉的，通常对于一个欠定线性系统，我们都会这么去做。因为系统欠定，我们会有多个解，事实上全部的解集构成一个线性空间，它的维数就是$n - m$。现在结合不等值约束$x \geq 0$，于是可行域$\Omega$就是这个解空间落在$\mathbb{R}^n$中的“第一象限”的部分。大家可以想象一下，一个空间被一个$n$维的第一象限切割下的部分，实际上是$n$维空间的一个“凸多面体”（在$x$的正方向可以开放，但如果封闭，那么整个可行域就是一个凸多面体）。现在对于全部$i \in \mathcal{N}$，我们令$x_i = 0$，从而构成了基础可行点，这时我们发现：

1. 基础可行点是确定的单点，而不是空间。因为它的每一个分量都是唯一确定的；
2. 它是$x_B$空间和$x_N$空间的交点，而且是一个单点，如果我们将$\Omega$视为$\mathbb{R}^n$空间中落在第一象限内的一个凸多面体，那么基础可行点实际上就是这个凸多面体和$x_B$空间相交的顶点。不同的基础可行点代表了不同的顶点；
3. 从上述几何意义出发，定理13.2的结论是显然的。同时，提示我们可以通过某种方式遍历全部的基础可行点（顶点）来寻找问题(13.1)的全局最优解（如果有，如果$\Omega$在$x$正方向是开放的，同时$f$沿$x$的正方向是下降的，那么显然最优值是$-\infty$）。

## 单纯形方法(the simplex method)

现在我们来具体实现上面第3点讨论的方法。我们已经根据下标$\mathcal{B}$和$\mathcal{N}$划分了$B$和$N$，以及$x_B$和$x_N$。而判定一个点是否为全局最优解的条件是(13.4)，为方便讨论起见，我们进一步将(13.4)中的$s$和$c$也做对应的划分：
$$
\begin{array}{ll}
s_B = [s_i]_{i \in \mathcal{B}}, & s_N = [s_i]_{i \in \mathcal{N}},\\
c_B = [c_i]_{i \in \mathcal{B}}, & c_N = [c_i]_{i \in \mathcal{N}}.
\end{array}
$$
结合(13.4)，我们有
$$
A x = Bx_B + Nx_N = b.
$$
而作为一个算法的起点，我们需要第一个基础可行点做初值，或者说我们需要选择一个$\mathcal{B}$，然后我们就有
$$
\begin{equation}
x_B= B^{-1}b, x_N = 0.
\tag{13.18}
\end{equation}
$$
这里显然我们应该考虑那些使$B$在求逆时尽可能简单的$\mathcal{B}$，这个也可以通过算法手段得到，我们后面再讨论，目前我们只在所有实际可能中选取一个“看上去最好”的$\mathcal{B}$。

**例 13.1** 
$$
\begin{array}{ll}
\min & -4 x_1 - 2x_2, \\
\mbox{s. t.} & x_1 + x_2 + x_3 = 5,\\
& 2x_1 + \frac{1}{2}x_2 + x_4 = 8,\\
& x \geq0.
\end{array}
$$
首先我们整理一下，这里
$$
A = \left[
\begin{array}{cccc}
1&1&1&0\\
2&\frac{1}{2}&0&1
\end{array}
\right],\quad 
b = \left[
\begin{array}{c}
5\\8
\end{array}
\right],\quad 
c = \left(-4, -2, 0, 0\right)^T.
$$
这里我们选择一个“尽可能好的”$\mathcal{B} = \{3, 4\}$，于是$\mathcal{N} = \{1, 2\}$，
$$
\begin{array}{rr}
B = \left[
\begin{array}{cc}
1 & 0\\
0 & 1
\end{array}
\right], & 
N = \left[
\begin{array}{cc}
1 & 1\\ 2 & \frac12
\end{array}
\right],\\
c_B = \left[
\begin{array}{c}
0\\0
\end{array}
\right], &
c_N = \left[
\begin{array}{c}
-4 \\ -2
\end{array}
\right].
\end{array}
$$
而我们的计算将从第一个基本可行点开始：
$$
x_B = B^{-1}\cdot b = \left[
\begin{array}{c}
x_3\\x_4
\end{array}
\right] = \left[
\begin{array}{c}
5\\8
\end{array}
\right], \quad 
x_N = \left[
\begin{array}{c}
x_1 \\ x_2
\end{array}
\right] = \left[
\begin{array}{c}
0 \\ 0
\end{array}
\right].
$$
即$x = (0, 0, 5, 8)^T$，注意此时
$$
f(x) = c^Tx = 0
$$
即是当前目标值。接下去迭代应该设法找到新的基础可行点来降低这个值，如果它还不是最优的话。

对一个基础可行点$x$判定是否最优，我们直接采用(13.4)。