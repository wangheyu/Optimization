# 关于两种正定修正的补充说明

问题：已知$A$为对称不充分正定矩阵，对$A$做适当的修正使修正后的矩阵$\tilde{A}$充分正定。这里充分正定要求$A$的最小特征值$\lambda_\min$满足$\lambda_\min \geq \delta > 0$，$\delta > 0$是一个显著大于机器精度的正数。

## 基于特征分解的修正

对$A$做谱分解：
$$
A = Q \Lambda Q^T,
$$
其中$Q$是正交阵，$\Lambda$是对角阵，且对角元$\lambda_i, i = 1, 2, \cdots, n, $是$A$全部特征值。由于$A$不充分正定，故存在$\lambda_i < \delta$，现令
$$
\tau_i = \left\{\begin{array}{ll} 0, &\lambda_i \geq \delta, \\ \delta - \lambda_i, &\lambda_i < \delta. \end{array}\right.
$$


+ 若令

$$
\Delta A = Q \mathrm{diag}(\tau_i)Q^T,
$$

则$\Delta A$是使$\tilde{A} = A + \Delta A$满足$\lambda_\min \geq \delta$的Frobenius范数最小的修正项；

+ 若令$\tau = \max(0, \delta - \lambda_\min(A))$，则

$$
\Delta A = \tau I
$$
则$\Delta A$是使$\tilde{A} = A + \Delta A$满足$\lambda_\min \geq \delta$的Euclidean范数（欧氏范数、2范数）最小的修正项.

基于谱分解的修正的优点是$\tilde{A}$的特征分解（特征值、特征向量）和$A$很接近。但$\tilde{A}$和$A$形式上并不接近。

## 基于Cholesky分解的修正

对$A$做$LDL^T$型Cholesky分解，则对比$MM^T$型分解，有
$$
M = LD^{\frac{1}{2}}.
$$
现令$\delta, \beta>0$，则算法（Matlab伪）：

```matlab
for j = 1:n
	C(j,j) = A(j,j)-sum(d(1:j-1)'.*L(j,1:j-1).^2);
	theta = 0;
	for i = j+1:n
		C(i,j) = A(i,j)-sum(d(1:j-1)'.*L(i,1:j-1).*L(j,1:j-1));
		if theta < abs(C(i,j))
			theta = abs(C(i,j));
		end
    end
    d(j) = max([abs(C(j,j)), (theta/beta)^2, \delta]);
    for i = j+1:n
    	L(i,j) = C(i,j)/d(j);
    end
    L(j,j) = 1.0;
end
D=diag(d);
```

满足

$\tilde{A} = LDL^T$正定，且有
$$
d_j \geq \delta, |m_{ij}| \leq \beta, j = 1,2,\cdots,n, i=j+1,\cdots,n.
$$
这种修正的特点是$\tilde{A}$形式上和$A$接近，但二者特征分解相去甚远。