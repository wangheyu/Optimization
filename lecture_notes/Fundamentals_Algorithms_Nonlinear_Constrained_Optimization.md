# 约束非线性优化算法基础

本章讨论如何在第12章的理论基础上求解一般形式的约束优化问题，也即
$$
\begin{equation}
\begin{array}{ll}
\min_{x\in\mathbb{R}^n} & f(x) \\
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
   其中$\mu > 0$称为惩罚因子(penalty parameter)。我们对一个递减的$\mu$求解(15.2)的无约束极小值，直到解充分接近约束问题(15.1)的解为止。罚函数(15.2)的问题是只有$\mu \to 0$时，它的解才收敛到(15.1)的解，而在对任何一个具体的$\mu > 0$，它的最优解和(15.1)总是有差异的。另一个思路是采用一个“精确”(exact)的罚函数
   $$
   f(x) + \mu \sum_{i \in \mathcal{E}}|c_i(x)|.
   $$
   当$\mu$充分大时，该函数提供的精确的最优解位置。然而该函数是不光滑的，因此给数值求解带来额外困难（但是是可以解决的）。所以这也是一种可以选择的方案。

   而增广Lagrange法则是针对Lagrange函数构建罚函数并将其转为无约束优化问题，称为增广Lagrange函数
   $$
   \mathcal{L}_A(x, \lambda;\mu) = f(x) - \sum_{i \in \mathcal{E}}\lambda_i c_i(x) + \frac{\mu}{2}\sum_{i \in \mathcal{E}}c_i^2(x).
   $$
   接下去我们先固定$(\mu, \lambda)$的值，寻找$\mathcal{L}_A$关于$x$的无约束最小值，一旦找到，就可以更新$\lambda$值，并且缩小$\mu$。然后在新的$(\mu, \lambda)$值下寻找新的$x$。这个过程可以不断重复，直到找到正确的$(x, \lambda)$。这个方法相比(15.2)的优点在于它的惩罚因子$\mu$不用趋于$0$，而是有一个确定的下界。

3. 第18章介绍序列化二阶规划(sequential quadratic programming, SQP)方法。它的思路是将(15.1)转为一系列二阶规划子问题，从而得到一系列搜索方向，最终收敛到正确的解。最基本的SQP方法，关于$(x_k, \lambda_k)$的搜索方向$p_k$如下定义：
   $$
   \begin{eqnarray}
   &\min& \frac{1}{2}p^T\nabla_{xx}^2\mathcal{L}(x_k, \lambda_k)p + \nabla f(x_k)^Tp, \tag{15.3a}\\
   &\mbox{s. t.}& \nabla c_i(x_k)^Tp + c_i(x_k) = 0, i \in \mathcal{E}, \tag{15.3b}\\
   && \nabla c_i(x_k)^Tp + c_i(x_k) = 0, i \in \mathcal{I}, \tag{15.3c}
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

   我们可能以及注意到，约束优化问题的主要困难来自如何处理不等值约束，特别是如何确定哪些不等值约束是活跃的。而活跃集方法的本质就是先猜测一个活跃集$\mathcal{A}(x^*)$，称为工作集(working set)，记作$\mathcal{W}$。然后我们就可以把在$\mathcal{W}$内的约束都设成等值。基于猜测$\mathcal{W}$，我们计算等值约束问题得到一个$x^*$，然后我们在完整的问题中检查它是否满足KKT条件。如果满足，则确定它为(15.1)的一个局部最优解（其实只是KKT点，但作为程序这里可以退出了）；如果不满足，则我们调整$\mathcal{W}$，重复操作直到找到最优点。这个问题就是基于这样一种认知：求解等值约束问题相比不等值约束问题要简单的多。