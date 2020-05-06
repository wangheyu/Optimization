## 二阶优化条件

在LICQ（或其他适当的约束规范）前提下，我们可以通过局部最优点$x^*$处的$\mathcal{F}(x^*)$来分析全部一阶局部信息应该具有的性质，也就是KKT条件。本节我们继续讨论二阶局部信息应该具有的性质，并最终给出二阶必要和充分条件。

首先，我们将焦点放在一阶条件中无法判定升降的方向，也即$w$满足
$$
w^T\nabla f(x^*) = 0.
$$
我们在KKT条件的框架下，将其定义为一个锥的形式，称为关键锥(critical cone)：
$$
\mathcal{C}(x^*, \lambda^*) = \left\{w \in \mathcal{F}(x^*) \left|\nabla c_i(x^*)^Tw = 0, \forall i \in \mathcal{A}(x^*)\cap\mathcal{I}, \lambda_i > 0\right.\right\}.
$$
或等价地，
$$
\begin{equation}
w \in \mathcal{C}(x^*, \lambda^*) \Leftrightarrow \left\{
\begin{array}{ll}
\nabla c_i(x^*)^T w = 0, & \forall i \in \mathcal{E},\\
\nabla c_i(x^*)^T w = 0, & \forall i \in \mathcal{A}(x^*)\cap\mathcal{I}, \lambda_i > 0\\
\nabla c_i(x^*)^T w \geq 0, &\forall i \in \mathcal{A}(x^*)\cap\mathcal{I}, \lambda_i = 0.
\end{array}
\right. \tag{12.53}
\end{equation}
$$
也即$w$方向在等值约束梯度方向上一阶不变，在活跃的且强互补的不等值约束梯度方向上一阶不变，在活跃且不严格互补的不等值约束梯度方向上上升（继续满足约束，但因为不严格互补，此信息没用）。由(12.53)，我们马上有结论
$$
\begin{equation}
w \in \mathcal{C}(x^*, \lambda^*) \Rightarrow \lambda_i^*\nabla c_i(x^*)^Tw = 0, \quad \forall i \in \mathcal{E}\cup\mathcal{I}.
\tag{12.54}
\end{equation}
$$
于是结合(12.34a)和(12.33)，得
$$
\begin{equation}
w \in \mathcal{C}(x^*, \lambda^*) \Rightarrow w^T\nabla f(x^*) = \sum_{i \in \mathcal{E} \cup \mathcal{I}} \lambda_i^* w^T\nabla c_i(x^*) = 0. \tag{12.55}
\end{equation}
$$
也即沿关键锥中的方向，一阶信息已经退化了，和之前无约束优化问题一样，如果要知道$\mathcal{C}$中的方向是上升还是下降，我们必须进一步考虑其二阶信息。

**例 12.7**

考虑问题
$$
\begin{equation}
\min x_1, \quad \mbox{s. t.} x_2 \geq 0, 1 - (x_1 - 1)^2 - x_2^2 \geq 0, \tag{12.56}
\end{equation}
$$
（草图）显然全局最优解是$x^* = (0, 0)^T$，而$\mathcal{A}(x^*) = \{1, 2\}$，对应有唯一的Lagrange乘子$\lambda^* = (0, 0.5)^T$。同时，注意到在$x^*$点，活跃约束的梯度分别为
$$
\nabla c_1(x^*) = (0, 1)^T, \quad \nabla c_2(x^*) = (2, 0)^T,
$$
满足LICQ条件。该问题的线性化可行方向为
$$
\mathcal{F}(x^*) = \{d | d \geq 0\},
$$
以及关键锥为
$$
\mathcal{C}(x^*, \lambda^*) = \{(0, w_2)^T | w_2 \geq 0\}. 
$$
我们注意到$f$在$\mathcal{F}(x^*)$方向内，都是非减的，因此和$x^*$是最优解不矛盾。但要判定$x^*$就是最优解，我们缺少了$\mathcal{C}(x^*, \lambda^*)$方向上的信息，因为在这个方向上，$\nabla f(x^*) = 0$，$f$究竟是上升还是下降，由Taylor展开的下一项（二阶项）决定。至此，二阶必要条件和充分条件都已经很清晰了。

**定理 12.5（二阶必要条件）** 假设$x^*$是问题(12.1)的局部最优解，且有LICQ成立。令$\lambda^*$是对应的Lagrange乘子，则
$$
\begin{equation}
w^T \nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)w \geq 0, \quad \forall w \in \mathcal{C}(x^*, \lambda^*). \tag{12.57}
\end{equation}
$$
**证明** 由$x^*$是(12.1)的局部最优解，全部趋于$x^*$的可行序列$\{z_k\}$对充分大的$k$均有
$$
f(z_k) \geq f(x^*).
$$
现在我们构建可行序列，使得其极限方向就是$w$。首先由于$w \in \mathcal{C}(x^*, \lambda^*) \subset \mathcal{F}(x^*)$，因此我们可以用证明引理12.2时的技巧（隐函数定理），构建可行序列和对应正序列，使得
$$
\begin{equation}
\lim_{k \to \infty} \frac{z_k - x^*}{t_k} = w, \tag{12.58}
\end{equation}
$$
也即
$$
\begin{equation}
z_k - x^* = t_k w + o(t_k). \tag{12.59}
\end{equation}
$$
而在构建$\{z_k\}$的时候，我们由(12.42)，已知
$$
c_i(z_k) = t_k \nabla c_i(x^*)^Tw, \forall i \in \mathcal{A}(x^*), \tag{12.60}
$$
（以上推导表明在$w$方向上，一阶Taylor展开仍然是成立的。）结合(12.33)，(12.60)和(12.54)，我们有
$$
\begin{eqnarray}
\mathcal{L}(z_k, \lambda^*) &=& f(z_k) - \sum_{i \in \mathcal{E}\cup\mathcal{I}}\lambda^*_i c_i(z_k) \\
&=& f(z_k) - t_k\sum_{i \in \mathcal{A}(x^*)} \lambda_i^* \nabla c_i(x^*)^Tw\\
&=& f(z_k), \tag{12.61}
\end{eqnarray}
$$
（然而一阶项实际是零，消失了。）为此我们继续将Taylor展开做下去，
$$
\begin{eqnarray}
\mathcal{L}(z_k, \lambda^*) &=& \mathcal{L}(x^*, \lambda^*) + (z_k - x^*)^T\nabla_x\mathcal{L}(x^*, \lambda^*) + \frac12(z_k - x^*)^T\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)(z_k - x^*) + o(\|z_k - x^*\|). \tag{12.62}
\end{eqnarray}
$$
在我们的例子中，考虑到互补性条件(12.34e)，已有$\mathcal{L}(x^*, \lambda^*) = f(x^*)$。而在(12.34a)中（一阶必要条件），全部一阶导数项也全部为零（关键锥），再由(12.59)，我们可以将(12.62)改写为
$$
\begin{equation}
\mathcal{L}(z_k, \lambda^*) = f(x^*) + \frac12t_k^2w^T\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*) + o(t_k^2). \tag{12.63}
\end{equation}
$$
结合(12.61)，有
$$
\begin{equation}
f(z_k) = f(x^*) + \frac12 t_k^2w^T\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)w + o(t_k^2).\tag{12.64}
\end{equation}
$$
若$w^T\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)w < 0$，由(12.64)，存在$f(z) < f(x^*)$对$k$充分大成立，和$x^*$是局部最优解矛盾。因此(12.57)必然成立。$\Box$

**定理 12.6（二阶充分条件）** 假设对可行点$x^* \in \mathbb{R}^n$存在Lagrange乘子$\lambda^*$使得KKT条件成立，并且同时有
$$
\begin{equation}
w^T \nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)w > 0, \quad \forall w \in \mathcal{C}(x^*, \lambda^*), w \neq 0. \tag{12.65}
\end{equation}
$$
则$x^*$是问题(12.1)的严格局部最优解。

证明自己看书。实际上书上的证明也不完整。但这个定理的结论还是很显然的。

**例12.8** 我们检查例12.2中出现的问题(12.8)的二阶局部信息。注意到
$$
f(x) = x_1 + x_2, \quad c_1(x) = 2 - x_1^2 - x_2^2 \geq 0,
$$
因此
$$
\mathcal{L}(x, \lambda) = (x_1 + x_2) - \lambda_1(2 - x_1^2- x_2^2),
$$
而$x^* = (-1, -1)^T$，$\lambda^* = \frac12$。于是Lagrange函数的Hessian阵为：
$$
\nabla_{xx}^2 \mathcal{L}(x^*, \lambda^*) = \left[
\begin{array}{cc}
2\lambda_1^* & 0 \\
0 & 2\lambda_1
\end{array}
\right] = \left[
\begin{array}{cc}
1 & 0\\
0 & 1
\end{array}
\right].
$$
这是一个正定矩阵，所以定理12.6的条件必然满足。于是$x^* = (-1, -1)^T$必然是严格局部最优点。（这个例子太绝对了一些）

**例 12.9** 更复杂一些的例子：
$$
\begin{equation}
\min -0.1(x_1 - 4)^2 + x_2^2, \quad \mbox{s.t.} x_1^2 + x_2^2 - 1 \geq 0, \tag{12.72}
\end{equation}
$$
（草图）我们的可行域是单位圆和以外部分。现在这是一个无界问题，而且目标函数可以取到$-\infty$，因此也无解。但它在边界上是否存在严格局部最优解呢？我们应该无法从草图上清晰地了解这一点，但仍然可以通过KKT条件和二阶条件来判定：
$$
\begin{eqnarray}
\nabla_x\mathcal{L}(x, \lambda) &=& \left[
\begin{array}{c}
-0.2(x_1 -4) - 2\lambda_1 x_1\\
2x_2 - 2\lambda_1x_2
\end{array}
\right],\tag{12.73a}\\
\nabla_{xx}^2\mathcal{L}(x, \lambda) &=& \left[
\begin{array}{cc}
-0.2 - 2 \lambda_1 & 0\\
0 & 2 - 2\lambda_2
\end{array}
\right]\tag{12.73b},
\end{eqnarray}
$$
这里$x^* = (1, 0)^T$，$\lambda^* = 0.3$满足KKT条件，且对应$\mathcal{A}(x^*) = \{1\}$。而
$$
\nabla c_1(x^*) = \left[\begin{array}{c}
2 \\ 0
\end{array}\right],
$$
所以关键锥为
$$
\mathcal{C}(x^*, \lambda^*) = \left\{(0, w_2)^T \left| w_2 \in \mathbb{R}\right.\right\}.
$$
确实是可行域的边缘切线，我们现在要判定沿这条边缘切线，目标函数是上升还是下降，因为一阶信息在这里退化了，所以要进一步考虑二阶信息。
$$
w^T\nabla_{xx}^2\mathcal{L}(x^*, \lambda^*)w = 
[0 \quad w_2]\left[
\begin{array}{cc}
-0.4 & 0 \\
0 & 1.4
\end{array}
\right]
\left[
\begin{array}{c}
0 \\ w_2
\end{array}
\right] = 1.4 2_2^2 > 0.
$$
因此由定理12.6，$x^* = (1, 0)^T$是严格局部最优解。

