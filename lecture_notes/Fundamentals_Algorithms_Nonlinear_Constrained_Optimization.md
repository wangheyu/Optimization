# 约束非线性优化算法基础

本章讨论如何在第12章的理论基础上求解一般形式的约束优化问题，也即
$$
\begin{equation}
\begin{array}{ll}
\displaystyle \min_{x\in\mathbb{R}^n} & f(x) \\
\mbox{s.t.} & c_i(x) = 0, i \in \mathcal{E},\\
& c_i(x) \geq 0, i \in \mathcal{I},
\end{array}
\tag{15.1}
\end{equation}
$$
其中目标函数$f$和约束$c_i$都是光滑、实值，定义域是$\mathbb{R}^n$的子集，此外 $\mathcal{E}$和$\mathcal{I}$均为有限指标集。

## 优化算法分类

优化算法并没有标准分类，我们这里的分类只是基于本书的剩余部分进行的。

1. 第16章介绍求解二次规划(quadratic)的算法。和线性规划一样，这是一类在建模和计算中广泛应用的模型和算法，因此单独分成一章讨论。二次规划几乎是所有（经典）高效率算法的基础，很多一般问题可以用序列化的二次规划子问题得到求解。在这一章，我们会介绍活跃集方法(active set method)，内点法(interior-point methods)和梯度投影方法(gradient projection method)。

2. 第17章介绍罚函数方法(penalty methods)和增广Lagrange方法(augmented Lagrangian methods)。罚函数的基本思想是联合目标和约束函数构建一个罚函数(penalty function)，然后将一个约束优化问题转为一个无约束优化问题。比如，如果只考虑(15.1)中的等值约束，那么我们可以构建二次罚函数
   $$
   \begin{equation}
   f(x) + \frac{\mu}{2} \sum_{i \in \mathcal{E}}c_i^2(x), \tag{15.2}
   \end{equation}
   $$
   其中$\mu > 0$称为惩罚因子(penalty parameter)。我们对一个递减的$\mu$求解(15.2)的无约束极小值，直到解充分接近约束问题(15.1)的解为止。罚函数(15.2)的问题是只有$\mu \to +\infty$时，它的解才收敛到(15.1)的解，而在对任何一个具体的$\mu > 0$，它的最优解和(15.1)总是有差异的。另一个思路是采用一个“精确”(exact)的罚函数
   $$
   f(x) + \mu \sum_{i \in \mathcal{E}}|c_i(x)|.
   $$
   当$\mu$充分大时，该函数提供的精确的最优解位置。然而该函数是不光滑的，因此给数值求解带来额外困难（但是是可以解决的）。所以这也是一种可以选择的方案。

   而增广Lagrange法则是针对Lagrange函数构建罚函数并将其转为无约束优化问题，称为增广Lagrange函数
   $$
   \mathcal{L}_A(x, \lambda;\mu) = f(x) - \sum_{i \in \mathcal{E}}\lambda_i c_i(x) + \frac{\mu}{2}\sum_{i \in \mathcal{E}}c_i^2(x).
   $$
   接下去我们先固定$(\mu, \lambda)$的值，寻找$\mathcal{L}_A$关于$x$的无约束最小值，一旦找到，就可以更新$\lambda$值，并且增大$\mu$。然后在新的$(\mu, \lambda)$值下寻找新的$x$。这个过程可以不断重复，直到找到正确的$(x, \lambda)$。这个方法相比(15.2)的优点在于它的惩罚因子$\mu$不用趋于$+\infty$，而是有一个确定的上界。

3. 第18章介绍序列化二阶规划(sequential quadratic programming, SQP)方法。它的思路是将(15.1)转为一系列二阶规划子问题，从而得到一系列搜索方向，最终收敛到正确的解。最基本的SQP方法，关于$(x_k, \lambda_k)$的搜索方向$p_k$如下定义：
   $$
   \begin{eqnarray}
   &\min& \frac{1}{2}p^T\nabla_{xx}^2\mathcal{L}(x_k, \lambda_k)p + \nabla f(x_k)^Tp, \tag{15.3a}\\
   &\mbox{s. t.}& \nabla c_i(x_k)^Tp + c_i(x_k) = 0, i \in \mathcal{E}, \tag{15.3b}\\
   && \nabla c_i(x_k)^Tp + c_i(x_k) \geq 0, i \in \mathcal{I}, \tag{15.3c}
   \end{eqnarray}
   $$
   其中$\mathcal{L}$是(12.33)所定义的Lagrange函数。这是一个二次规划问题，它的解$p$用来将$x_k$更新到$x_k + p$。对它的进一步改进可以考虑施加一个信任域半径，从而确保收敛的稳定性。也可以考虑对$\nabla_{xx}^2 \mathcal{L}$采用拟Newton逼近以提高计算效率。

   这种方法的一种变体被称为序列化线性-二次规划(sequential linear-quadratic programming)，思路是将$p_k$的计算分成两个阶段。在第一阶段，我们忽略(15.3a)的第一项（二次项），从而将问题退化为一个线性规划，然后在一个信任域的约束下求解。第二阶段，我们将第一阶段解的全部活跃约束作为等值约束，并忽略其余约束，求解这个二次规划来得到真正的$p_k$。

4. 第19章我们研究非线性规划的内点法(interior-point methods for nonlinear programming)。这个方法可以看作是第14章中介绍的针对线性规划的主-对偶内点法的一个扩展。我们也可以将它看作是障碍法(barrier methods)，障碍函数为
   $$
   \begin{eqnarray}
   &\min_{x, s}& f(x) - \mu \sum_{i = 1}^m \log s_i, \tag{15.4a} \\
   &\mbox{s.t.}& c_i(x) = 0, i \in \mathcal{E}, \tag{15.4b} \\
   && c_i(x) - s_i = 0, i \in \mathcal{I}, \tag{15.4c}
   \end{eqnarray}
   $$
   这里$\mu$是障碍因子，$s_i$是松弛变量。障碍函数法本质上也是一种罚函数方法，它要求迭代点在可行域内部，当迭代路径接近可行域边界时，由于障碍项的存在，导致目标值快速趋于$+\infty$，从而阻挡迭代突破边界。由于它的迭代序列整体保持在可行域内部，故称内点法。

本章的剩余部分，我们讨论一下用于以上各方法的一些基础算法和建模技巧。

## 不等值约束的组合困难

我们可能以及注意到，约束优化问题的主要困难来自如何处理不等值约束，特别是如何确定哪些不等值约束是活跃的。而活跃集方法的本质就是先猜测一个活跃集$\mathcal{A}(x^*)$，称为工作集(working set)，记作$\mathcal{W}$。然后我们就可以把在$\mathcal{W}$内的约束都设成等值。基于猜测$\mathcal{W}$，我们计算等值约束问题得到一个$x^*$，然后我们在完整的问题中检查它是否满足KKT条件。如果满足，则确定它为(15.1)的一个局部最优解（其实只是KKT点，但作为程序这里可以退出了）；如果不满足，则我们调整$\mathcal{W}$，重复操作直到找到最优点。这个问题就是基于这样一种认知：求解等值约束问题相比不等值约束问题要简单的多。然而工作集可能的组合至多会有$2^{|\mathcal{I}|}$个，因为对每一个$i \in \mathcal{I}$，都可能有$i \in \mathcal{W}$和$i \notin \mathcal{W}$两种情况。如此巨大的基数，使得构建一个算法来遍历全部的工作集组合实际不可能。这一点我们称为不等值约束的组合困难。

下面的例子表明，即便是很少个数的不等值约束，活跃集的判定仍然不是一个简单的任务。

**例 15.1** 考虑
$$
\begin{equation}
\begin{array}{ll}
\displaystyle \min_{x,y} & f(x, y) \overset{\mathrm{def}}{=} \frac{1}{2} (x - 2)^2 + \frac{1}{2}(y - \frac{1}{2})^2, \\
\mbox{s. t.} & (x + 1)^{-1} - y - \frac{1}{4} \geq 0, \\
& x \geq 0,\\
& y \geq 0.
\end{array} \tag{15.5}
\end{equation}
$$
从草图我们可以看到，在解$(x^*, y^*)^T = (1.953, 0.089)^T$，活跃的约束只有$c_1$。但是这个问题的全部可能的活跃集有$2^3 =8$种。首先我们可以判定，由于目标函数的无约束全局最优点在可行域之外，所以活跃集不可能是空集；此外，三个约束不可能同时活跃。因此剩下需要逐一判定的可能是3种单一约束活跃，和3种两个约束一起活跃的情形，共$6$种。这里我们选择三种情况稍加分析：

+ $\mathcal{W} = \{2\}$，于是$x = 0$。我们只保留这个等值约束求解问题，即
  $$
  \begin{equation}
  \begin{array}{ll}
  \displaystyle \min_{x,y} & f(x, y) \overset{\mathrm{def}}{=} \frac{1}{2} (x - 2)^2 + \frac{1}{2}(y - \frac{1}{2})^2, \\
  \mbox{s. t.} & x = 0.\\
  \end{array}
  \end{equation}
  $$
  显然有$(x^*, y^*)^T = (0, \frac{1}{2})^T$。检查KKT条件，$\lambda_1 = \lambda_3 = 0$，于是
  $$
  \nabla f(x^*, y^*) - \lambda_2 \nabla c_2(x^*, y^*)= \left[
  \begin{array}{c}
  x^* - 2\\
  y^* - \frac{1}{2}
  \end{array}
  \right] - \lambda_2 \left[
  \begin{array}{c}
  1 \\ 0
  \end{array}
  \right]=\left[
  \begin{array}{c}
   - 2\\
  0
  \end{array}
  \right] - \left[
  \begin{array}{c}
   \lambda_2\\
  0
  \end{array}
  \right] = 0 \Rightarrow \lambda_2 = -2.
  $$
  不满足KKT条件。
+ $\mathcal{W} = \{1, 3\}$，于是
  $$
  \left\{\begin{array}{rcl}
  (x + 1)^{-1} - y - \frac{1}{4} &=& 0\\
  y &=& 0
  \end{array}\right. \Rightarrow \left\{
  \begin{array}{rcl}
  x &=& 3 \\
  y &=& 0
  \end{array}
  \right..
  $$
  即$(x^*, y^*)^T = (3, 0)^T$。检查KKT条件，有$\lambda_1 = -16$，$\lambda_3 = -16.5$。不满足KKT条件。

+ $\mathcal{W} = \{1\}$，只有等值约束$c_1(x) = 0$。可以验证它满足KKT条件。（后两种情况自己验证补全一下。）

即便是这么简单的例子，我们要完全遍历全部可能情况也是困难的。不过我们注意到在实际求解的时候，我们先通过函数的局部信息，排除了一些实际不可能的情况。第16章介绍的活跃集方法，就会充分利用这种技巧（本质上单纯形法也属于这一类算法）。而19章介绍的内点法，则从根本上试图避开这一组合困难，而是直接从非线性问题的迭代行为入手去求解。

## 变量消去

对于线性问题或问题的线性部分利用消元手段降低问题的维数是一个很自然的想法。但由于我们的问题本质上是非线性的，消元一定要小心，因为对问题的不同代数形式表示可能会导致病态的计算行为。

我们先看一个有好结局的例子。
$$
\begin{array}{ll}
\min & f(x) = f(x_1, x_2, x_3, x_4), \\
\mbox{s.t.} & x_1 + x_3^2 - x_4x_3 = 0,\\
& -x_2 + x_4 + x_3^2 = 0,
\end{array}
$$
这里令
$$
x_1 = x_4x_3 - x_3^2, \quad x_2 = x_4 + x_3^2,
$$
则原问题化为等价的无约束问题
$$
\min h(x_3, x_4) = f(x_4x_3 - x_3^2, x_4 + x_3^2, x_3, x_4).
$$
这似乎是一个很不错的办法，但下一个例子表明了这种做法的危险性。

**例 15.2(Flectcher[101])** 
$$
\begin{array}{ll}
\min & x^2 + y^2 \\
\mbox{s.t.} & (x - 1)^3 = y^2.
\end{array}
$$
从草图中不难发现，最优解为$(x^*, y^*) = (1, 0)^T$。这个问题似乎也可以用消元法变为无约束问题
$$
h(x) = x^2 + (x - 1)^3.
$$
然而，因为$x \to -\infty$，$h(x) \to -\infty$，问题是无解的。毛病出在哪里？

其实还是表现形式的问题。原问题中的
$$
(x - 1)^3 = y^2 \geq 0,
$$
但在消元之后，这个隐含的信息$x - 1 \geq 0$被忽略了。这就导致了问题解的不一致。这个例子表明，对非线性方程组进行消元是非常危险的，因为可能会导致一些隐含的信息被不正确地抛弃，而且这种意外很难在推导和计算过程中被发现和回溯。由于这个原因，一般不主张对非线性约束直接消元，而是要先把约束做线性化，然后再对线性问题消元。

### 对线性约束的简单消元

现在考虑目标函数是非线性，而约束是等值线性的特例，
$$
\begin{equation}
\begin{array}{ll}
\min & f(x) \\
\mbox{s.t.} & Ax = b.
\end{array}
\tag{15.6}
\end{equation}
$$
这里$A$是$m \times n$矩阵，且$m \leq n$，且不妨设$A$行满秩。注意我们目前的情况和之前线性规划的情况很接近。事实上，这种思路可以认为是单纯形法的一种泛化。我们选择$A$的$m$列构成一个非奇异的$m \times m$方阵$B$，为了将这$m$列交换到$A$的前$m$列，我们需要一个$n \times n$的交换阵$P$，于是有
$$
\begin{equation}
A P = [B | N], \tag{15.7}
\end{equation}
$$
其中$N$是$A$中剩下的$n - m$列。对应这两个矩阵的列标，我们有$x_B \in \mathbb{R}^m$和$x_N \in \mathbb{R}^{n - m}$：
$$
\begin{equation}
\left[
\begin{array}{c}
x_B \\
x_N
\end{array}
\right] = P^Tx, \tag{15.8}
\end{equation}
$$
这里我们称$x_B$为基础向量(basic variables)，$B$为基础矩阵(basic matrix)。注意$P^TP = I$，于是有
$$
b = Ax = AP(P^Tx) = Bx_B + Nx_N.
$$
这里基础向量直接可以消去
$$
\begin{equation}
x_B = B^{-1}b - B^{-1}Nx_N. \tag{15.9}
\end{equation}
$$
也就是任何可行点，本质上只需要确定$x_N$就行了，而$x_B$可以通过(15.9)确定。直接代入(15.6)，得等价的无约束优化问题
$$
\begin{equation}
\min_{x_N} h(x_N) \overset{\mathrm{def}}{=} f\left(P\left[
\begin{array}{c}
B^{-1}b - B^{-1}Nx_N\\
x_N
\end{array}
\right]\right). \tag{15.10}
\end{equation}
$$
这种方法就是所谓的简单变量消去法(simple elimination of variables)。以上讨论显示了只含等值线性约束的优化问题，和一个无约束优化问题的等价性。

**例 15.3**
$$
\begin{eqnarray}
\min && \sin(x_1 + x_2) + x_3^2 + \frac{1}{3}(x_4 + x_5^4 + \frac{x_6}{2}), \tag{15.11a}\\
\mbox{s.t.}&&8x_1 - 6x_2 + x_3 + 9x_4 + 4x_5 = 6, \tag{15.11b}\\
&& 3x_1 + 2 x_2 - x_4 + 6x_5 + 4x_6 = -4.
\end{eqnarray}
$$


我们通过列交换阵$P$，使得$x^T = (x_3, x_6, x_1, x_2, x_4, x_5)^T$，也即
$$
AP = \left[
\begin{array}{cc|cccc}
1 & 0 & 8 & -6 & 9&4\\
0&4 &3 &2&-1 &6
\end{array}
\right].
$$
基础矩阵$B$是一个对角阵，求逆简单。根据(15.9)有
$$
\begin{equation}
\left[
\begin{array}{c}
x_3 \\
x_6
\end{array}
\right] = -\left[
\begin{array}{cccc}
8&-6&9&4\\
\frac{3}{4}&\frac{1}{2}&-\frac{1}{4}&\frac{3}{2}
\end{array}
\right]
\left[
\begin{array}{c}
x_1\\x_2\\x_4\\x_5
\end{array}
\right] + \left[
\begin{array}{c}6\\-1\end{array}
\right]. \tag{15.12}
\end{equation}
$$
将$x_3$，$x_6$代入(15.11a)，原问题变为
$$
\begin{equation}
\begin{array}{l}
\displaystyle \min_{x_1, x_2, x_4, x_5} & \sin(x_1 + x_2) + (8 x_1 - 6x_2 + 9x_4 + 4x_5 -6)^2 \\
&\displaystyle + \frac{1}{3}\left(x_4 + x_5^4 -\left[\frac{1}{2} + \frac{3}{8}x_1 + \frac{1}{4}x_2 -\frac{1}{8}x_4 + \frac{3}{4}x_5\right]\right) 
\end{array} \tag{15.13}
\end{equation}
$$
当然，$B$的选择是不唯一的，但如果另选$A$的两列作为$B$，那么求逆矩阵就不会那么简单。

以上的消去过程实际上就是Gauss消去法。在理想情况下，我们可以选择一个即容易分解，条件数又不大的矩阵$B$。有些技术甚至能将这一算法应用在大规模稀疏矩阵上，并在消去过程中保持稀疏性。比如著名的MA48算法和相应的计算库HSL[96]（我居然没听说过，恰饭嫌疑）。但是如接下去的讨论，Gauss消去策略并不能保证总是找到最好的基础矩阵。

以下，我们忽略选列和交换过程已简化讨论。也就是我们总是不妨设$P = I$。由(15.8)和(15.9)，我们知道对任何可行点$x$，以及线性约束(15.6)，可以有如下表示
$$
\begin{equation}
\left[
\begin{array}{c}
x_B \\ x_N
\end{array}
\right] = x = Yb + Zx_N, 
\end{equation}\tag{15.14}
$$
其中
$$
\begin{equation}
Y = \left[
\begin{array}{c}
B^{-1} \\ 0
\end{array}
\right], \quad Z = \left[
\begin{array}{c}-B^{-1}N \\ I\end{array}
\right]. \tag{15.15}
\end{equation}
$$
注意$Z$有$n - m$个线性无关的列（注意它下方的单位块），且$AZ = 0$，因此$Z$的列事实上是$A$的零空间的基。并且
$$
Yb = \left[
\begin{array}{c}
B^{-1} \\ 0
\end{array}
\right]b = \left[
\begin{array}{c}B^{-1}b \\ 0 \end{array}\right]
$$
事实上是欠定方程组$A x = b$的一个特解。因此Gauss消去法的本质，是将一个可行点分解成一个$Ax = b$的特解和零空间的一个平移之和（见(15.14)）。如果把$Ax = b$想象成一条直线，$x_B = x_1$，$x_N = x_2$，那么这种分解本质上就是取了直线在$x_1$轴上的截距。但是如果$Ax = b$本身是一条和$x_1$轴几乎平行的直线，那么这个截距的计算必然是很不稳定的（表现为数值扰动极大，或者矩阵病态）。这种情况下，我们更应该选择$x_2$而不是$x_1$作为$B$。但是究竟选择哪个$B$更好，这一点在消去过程中无法直接体现。为此，我们有以下更一般的策略讨论。

### 对线性约束的一般性约减讨论

为了使(15.14)和(15.15)更加一般化，我们选择矩阵$Y \in \mathbb{R}^{n \times m}$和$Z \in R^{n \times (n - m)}$，且
$$
\begin{equation}
[Y | Z] \in \mathbb{R}^{n \times n}, \quad AZ = 0. \tag{15.16}
\end{equation}
$$
这里$[Y|Z]$非奇异。由(15.15)，$Z$的列仍然是$A$的零空间的基。由$A$列满秩，
$$
A[Y | Z] = [AY | 0],
$$
于是$m \times m$矩阵$AY$是非奇异的（否则$[Y|Z]$奇异）。我们现在将$Ax = b$的解展开成
$$
\begin{equation}
x = Y x_Y + Zx_Z, \tag{15.17}
\end{equation}
$$
其中$x_Y \in \mathbb{R}^m$，$x_Z \in \mathbb{R}^{n - m}$，将其代入$Ax = b$，有
$$
Ax = (AY)x_Y = b.
$$
因为$AY$非奇异，所以有
$$
\begin{equation}
x_Y = (AY)^{-1}b. \tag{15.18}
\end{equation}
$$
再代回(15.17)，有
$$
\begin{equation}
x = Y(AY)^{-1}b + Zx_Z \tag{15.19}
\end{equation}
$$
满足约束$Ax = b$，不论$x_Z$如何取值。于是，问题(15.6)可以写成无约束问题
$$
\begin{equation}
\min_{x_Z} f(Y(AY)^{-1}b + Zx_Z). \tag{15.20}
\end{equation}
$$
这里对$Y$和$Z$最理想的选择是通过对$A^T$的QR分解：
$$
\begin{equation}
A^T\Pi = [Q_1\quad Q_2]\left[
\begin{array}{c}R \\ 0\end{array}
\right], \tag{15.21}
\end{equation}
$$
其中$[Q_1 \quad Q_2]$是正交矩阵，$Q_1$和$Q_2$分别是$n \times m$和$n \times (n - m)$的正交块，而$R$是$m \times m$的上三角块。而$n\times  n$阵$\Pi$是交换矩阵（具体的QR分解算法见任何一本数值代数教材，或者本书附录(A.24)）。现在，可以令$Y = Q_1$，$Z = Q_2$。这样$Y$和$Z$的列向量一起构成了$\mathbb{R}^n$空间的一组正交基。由(15.21)，有
$$
AY = \Pi R^T, \quad AZ = 0.
$$
也即$Y$和$Z$就是我们寻找的分解，并且$AY$具有和$R$和$A$相同的条件数。由(15.19)，我们看到$Ax = b$可以表示为
$$
x = Q_1R^{-T}\Pi^Tb + Q_2x_Z,
$$
$x_Z$可任取，而$R^{-T}\Pi^Tb$的计算是简单的。（也即这一过程的实际计算量在QR分解。）

通过简单的计算，我们可以看出$Q_1R^{-T}\Pi^Tb$可以被写成$A^T(AA^T)^{-1}b$。而后者是问题
$$
\min\|x\|^2, \quad \mbox{s.t.} Ax = b;
$$
的解，即零元在$Ax = b$空间上的正交投影，或者说，$Ax = b$的最小二乘解。

尽管这种基于QR分解的消去算法具有极好的数值稳定性，但对于大规模稀疏矩阵，QR分解的代价要高于Gauss消去。因此Gauss消去仍然是一种活跃的算法，或者，我们有时寻求二者的一种平衡。（见习题15.7）

### 不等值约束的效果

和等值约束一起出现的时候，消去策略对不等值约束未必有效。假设(15.11)同时有不等值约束$x \geq 0$，那么对等值约束做了之前的消元之后，剩下问题的约束变为
$$
\begin{array}{rcl}
(x_1, x_2, x_4, x_5) &\geq &0, \\
8x_1 - 6x_2 + 9 x_4 + 4 x_5 &\leq& 6, \\
\frac{3}{4}x_1 + \frac{1}{2}x_2 - \frac{1}{4}x_4 + \frac{3}{2}x_5 &\leq&-1.
\end{array}
$$
变得更加复杂，这样的消元不如不消。

除非，比如(15.11)增加一个比较一般的不等值约束
$$
3x_1 + 2x_3 \geq 1,
$$
那么消元之后，新的问题附带一个不等值约束
$$
-13x_1 + 12x_2 - 18x_4-8x_5 \geq -11. \tag{15.23}
$$
这种情况至少约束没有变得更加复杂，也许可以接受。