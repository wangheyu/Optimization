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
也即$w$方向在等值约束梯度方向上一阶不变，在活跃的且强互补的不等值约束梯度方向上一阶不变，在活跃且不严格互补的不等值约束梯度方向上上升（继续满足约束，但因为不严格互补，此信息没用）。

