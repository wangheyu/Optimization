## 评价函数和筛选(Merit Functions and Filters)

另一种通用的策略是将目标函数和约束函数看作一个整体，通过迭代法来求解。我们的目的是寻求满足约束并使目标函数最小化的点。如果一个迭代步，它目标函数有显著下降，但是约束却有微不足道的违背，我们是否应该接受这一步？

这是一个没有一般性答案的问题。我们一般通过算法技术来平衡这两种需求。其中一种技术称为评价函数，相当于综合下降和违背的效应来考虑；另一种称为筛选，相当于在一定的原则之下做妥协。

### 评价函数

提出评价函数本质上是试图将约束优化问题转成一个无约束优化问题。将目标下降和约束违背整合成一个评价目标。比如针对(15.1)提出
$$
\begin{equation}
\phi_1(x; \mu) = f(x) + \mu \sum_{i \in \mathcal{E}}|c_i(x)| + \mu\sum_{i \in \mathcal{I}}[c_i(x)]^-, \tag{15.24}
\end{equation}
$$
这里记号$[z]^- = \max\{0, -z\}$，正数$\mu$称为惩罚因子，(15.24)称为$l1$惩罚函数。尽管这种惩罚函数是不光滑的，但它具有精确性(exact)。

**定义 15.1（精确评价函数）**评价函数$\phi(x; \mu)$称为是精确的，若存在正数$\mu^*$使任意$\mu > \mu^*$，(15.1)的任何局部最优解都是$\phi(x; \mu)$的无约束局部最优解。（这里主要强调$\mu^* < +\infty$。）

对$l1$评价函数$\phi_1(x; \mu)$，令
$$
\mu^* = \max\{|\lambda_i^*|, i \in \mathcal{E} \cup \mathcal{I}\},
$$
这里$\lambda_i^*$是对应局部最优解$x^*$的Lagrange乘子，则可以证明它是精确的（定理17.3）。由于$\lambda^*$不是预知的，故在实际算法中，我们通过自适应调节使惩罚因子充分大。

另一个有用的评价函数是精确$l2$函数，它也是针对等值约束的，
$$
\begin{equation}
\phi_2(x; \mu) = f(x) + \mu \|c(x)\|_2. \tag{15.25}
\end{equation}
$$
注意这里$2-$范数项并没有平方，这使得它成为精确但不光滑评价函数。它在$c(x) = 0$点的导数不存在。

既光滑又精确的评价函数也是存在的。对等值约束，有Fletcher增广Lagrange函数
$$
\begin{equation}
\phi_F(x; \mu) = f(x) - \lambda(x)^Tc(x) + \frac{1}{2} \mu \sum_{i \in \mathcal{E}} c_i(x)^2,
\tag{15.26}
\end{equation}
$$
其中$\mu > 0$是惩罚因子，
$$
\begin{equation}
\lambda(x) = \left[A(x)A(x)^T\right]^{-1}A(x)\nabla f(x). 
\tag{15.27}
\end{equation}
$$
这里$A(x)$表示$c(x)$的Jacobi矩阵，和之前的定义是相关的。尽管这个函数性质很好，但是涉及到$\lambda(x)$的计算代价过大，限制了它的应用。在实际计算中，标准的增广Lagrange函数和相应的自适应迭代格式更受欢迎。

（书上这里有点啰嗦，重复了一些已经讨论过的内容，被我忽略了。）

**筛选(Filters)** （这个翻译可以讨论） 

过滤的思想来源是多目标优化(multiobjective optimization)。我们做约束优化其实并行地有两个目标：最小化目标函数，同时满足约束。如果我们定义一个不可行测度(measure of infeasibility)如下：
$$
\begin{equation}
h(x) = \sum_{i \in \mathcal{E}}|c_i(x)| + \sum_{i \in \mathcal{I}}[c_i(x)]^-, \tag{15.31}
\end{equation}
$$
则可以把我们的两个目标写成
$$
\begin{equation}
\left\{
\begin{array}{ll}
\displaystyle \min_x & f(x), \\
\displaystyle \min_x & h(x).
\end{array}
\right. \tag{15.32}
\end{equation}
$$
和评价函数的思路不同，筛选并不将两个优化目标统一成一个优化问题，而是分别处理。这里我们引入一个被控制(dominated)的概念：

**定义 15.2**

1. 第$k$步的目标对$(f_k, h_k)$称为是控制另一个目标对$(f_l, h_l)$的，如果$f_k \leq f_l$且$h_k \leq h_l$；（即两个目标都更好。）
2. 一个筛选有一系列目标对$(f_l, h_l)$构成，其中任两个目标对互相之间都不存在控制；
3. 一个迭代值$x_k$被称为是被筛选接受的，如果$(f_k, h_k)$没有被筛选中的任何一个目标对所控制。

当一个迭代值$x_k$是被筛选接受的，我们将$(f_k, h_k)$加入该筛选，同时将所有被$(f_k, h_k)$所控制的目标对移除筛选。如此我们构建了一个迭代法的更新策略，我们通过调整算法的其它参数，扩充筛选序列，直到达到最优解。在实际计算中，筛选策略有时具有非常好的效果（然而理论上仍然在研究中）。这一策略会在第18章得到应用和进一步讨论。

## Maratos效应

不论是评价函数还是筛选策略，在实际计算中都可能不收敛，因为它们有时会排除能取得最大进展的迭代步。这一现象被称为Maratos效应(Maratos effect)，它由Maratos第一次发现并发表于[199]。我们用以下例子来具体描述：

**例 15.4**（Powell [255]）
$$
\begin{equation}
\begin{array}{ll}
\min & f(x_1, x_2)  = 2(x_1^2 + x_2^2 - 1) - x_1,\\
\mbox{s.t.} & x_1^2 + x_2^2 = 1.

\end{array} \tag{15.34}
\end{equation}
$$
显然最优解是$x^* = (1, 0)^T$，对应的Lagrange乘子是$\lambda^* = \frac{3}{2}$，并且
$$
\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*) = I.
$$
考虑迭代序列$x_k = (\cos\theta, \sin\theta)^T$，对任何$\theta$，该点都是可行的。对一个迭代步长
$$
\begin{equation}
p_k = \left(
\begin{array}{c}
\sin^2\theta \\
-\sin\theta\cos\theta
\end{array}
\right), \tag{15.35}
\end{equation}
$$
我们有迭代改进
$$
x_k + p_k = \left(
\begin{array}{c}
\cos\theta + \sin^2\theta \\
\sin\theta(1 - \cos\theta)
\end{array}
\right).
$$
由基本三角恒等式(elementary trigonometric identities)，有
$$
\|x_k + p_k - x^*\|_2 = 2 \sin^2\left(\frac{\theta}{2}\right), \quad \|x_k - x^*\|_2 = 2\left|\sin\left(\frac{\theta}{2}\right)\right|,
$$
于是
$$
\frac{\|x_k + p_k - x^*\|_2}{\|x_k - x^*\|_2^2} = \frac{1}{2}.
$$
也即这个迭代是$Q-$二阶收敛的。然而，
$$
\begin{array}{rcl}
f(x_k + p_k) &=& \sin^2\theta - \cos\theta > -\cos\theta = f(x_k) \\
c(x_k + p_k) &=& \sin^2\theta > c(x_k) = 0.
\end{array}
$$
也即同时有目标函数上升，且约束违背增强。而且对任何$\theta$非零，这一现象都存在，不论$x_k$有多靠近$x^*$。

我们看到，只要我们采取形如
$$
\phi(x; \mu) = f(x) + \mu h(c(x))
$$
的评价函数，其中$h(0) = 0$，那么它就会拒绝(15.35)那样的高效率改进。显然，筛选策略也会拒绝这一改进。在实际计算中，Maratos效应会造成收敛变慢，因为它会拒绝本质上造成全部二阶收敛的步长。解决的办法如下：

1. 使用没有Maratos效应的评价函数，比如Fletcher增广Lagrange函数(15.26)；
2. 在评价和筛选过程中加入二阶信息；
3. 允许目标函数或约束违背有适当的上升。

