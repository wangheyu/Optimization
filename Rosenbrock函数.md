# Rosenbrock函数

多元 Rosenbrock 函数的一种形式如下：
$$
f(x) = \sum_{i = 1}^{n - 1}\left[(1 - x_i)^2 + 100 (x_{i + 1} - x_i^2)^2\right], x \in \mathbb{R}^n.
$$
其梯度为：
$$
\begin{eqnarray}
\frac{\partial f}{\partial x_1} &=& \frac{\partial}{\partial x_1}\left[(1 - x_1)^2 + 100(x_2 - x_1^2)^2\right] \\
&=& 2(x_1 - 1) + 400 x_1 (x_1^2 - x_2);\\
\frac{\partial f}{\partial x_i} &=& \frac{\partial}{\partial x_i}\left[(1 - x_i)^2 + 100(x_i - x_{i - 1}^2)^2 + 100(x_{i + 1} - x_i^2)^2\right]\\
&=&2(x_i - 1) + 400 x_i (x_i^2 - x_{i+1}) + 200 (x_i - x_{i - 1}^2),\\
&& i = 2,3, \cdots, n - 1;\\
\frac{\partial f}{\partial x_n} &=& 200(x_n - x_{n - 1}^2). 
\end{eqnarray}
$$
其 Hessian 为：
$$
\begin{eqnarray}
\frac{\partial^2f}{\partial x_1^2} &=& 1200 x_1^2 - 400 x_2 + 2;\\
\frac{\partial^2f}{\partial x_1 \partial x_2} &=& -400 x_1;\\
\frac{\partial^2f}{\partial x_i \partial x_{i - 1}} &=& -400 x_{i - 1};\\
\frac{\partial^2f}{\partial x_i^2} &=& 1200 x_i^2 - 400 x_{i + 1} + 202;\\
\frac{\partial^2f}{\partial x_i \partial x_{i + 1}} &=& -400 x_i,\\
&& i = 2, 3, \cdots, n - 1;\\
\frac{\partial^2f}{\partial x_n \partial x_{n - 1}} &=& -400 x_{n - 1};\\
\frac{\partial^2f}{\partial x_n^2} &=& 200;\\
\end{eqnarray}
$$

