<h1>
    <div align='center'>
    	第四章 线性方程组数值解法（2）
        J迭代法，G-S迭代法的敛散性
    </div>
</h1>

## Jacobi 迭代法的敛散性判断

### Jacobi 迭代法的敛散性

**将 Jacobi 迭代法写成矩阵迭代形式**：

对线性方程组 $Ax = b, \ A \in \mathbb{R}^{n \times n}$，将 $A$ 按照对角、下三角、上三角分裂
$$
A = D - L - U
$$
其中

- $D$：只包含矩阵 $A$ 的对角元 $a_{ii}$
- $L$：只包含矩阵 $A$ 的严格下三角部分元素，并全部取 **负号**（$-a_{ij}, \ i > j$）

- $U$：只包含矩阵 $A$ 的严格上三角部分元素，并全部取 **负号**（$-a_{ij}. \ i < j$）

因此 $Ax = b$ 可以表示成：
$$
(D - L - U)x = b
$$
改写为迭代格式：
$$
Dx = (L + U)x + b
$$

$$
x^{(k + 1)} = D^{-1}(L + U)x^{(k)} + c
$$

其中 $c = D^{-1}b$（$c$ 对于敛散性判断并不重要）

可以写成标准迭代格式：
$$
x^{(k + 1)} = B x^{(k)} + c
$$
其中 $B = D^{-1}(L + U)$ ，可以进一步化简，得到：
$$
B = D^{-1}(L + U) = I - D^{-1} A
$$
要判断 Jacobi 迭代法是否收敛（Jacobi 法对任意初值收敛），即判断 **矩阵 $B$ 的谱半径是否 $\rho(B) < 1$** .

**方法总结**：

1. 把 $A$ 分解为 $D - L - U$
2. 写出迭代格式进一步得到 Jacobi 迭代矩阵 $B = D^{-1}(L + U) = I - D^{-1} A$
3. 求 $B$ 的谱半径 $\rho(B)$
   - 若 $\rho(B) < 1$，则 Jacobi 迭代收敛
   - 否则，不收敛



### 例题 1

> $A = \begin{pmatrix} 1 & 0 & 2 \\ 3 & 2 & 1 \\ 0 & 4 & 3 \end{pmatrix}$，判断 Jacobi 迭代法对该矩阵是否收敛。

解答：

- 将 $A$ 分解为 $A = D - L - U$
  $$
  D = \begin{pmatrix}
  1 & 0 & 0 \\
  0 & 2 & 0 \\
  0 & 0 & 3 
  \end{pmatrix}
  $$

- Jacobi 迭代矩阵 $B$
  $$
  \begin{align*}
  B &= I - D^{-1}A \\
  &= \begin{pmatrix}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1
  \end{pmatrix} 
  - \begin{pmatrix}
  1 & 0 & 0 \\
  0 & \frac{1}{2} & 0 \\
  0 & 0 & \frac{1}{3}
  \end{pmatrix}
  \begin{pmatrix}
  1 & 0 & 2 \\ 
  3 & 2 & 1 \\ 
  0 & 4 & 3
  \end{pmatrix} \\
  &= \begin{pmatrix}
  0 & 0 & -2 \\\
  -\frac{3}{2} & 0 & -\frac{1}{2} \\
  0 & -\frac{4}{3} & 0
  \end{pmatrix}
  \end{align*}
  $$

- 求 $B$ 的谱半径判断收敛性
  $$
  \lambda I - B = \begin{pmatrix}
  \lambda & 0 & 2 \\
  \frac{3}{2} & \lambda & \frac{1}{2} \\
  0 & \frac{4}{3} & \lambda
  \end{pmatrix}
  $$

  $$
  \det(\lambda I - B) = \lambda \cdot (\lambda ^ 2 - \frac{4}{3} \times \frac{1}{2}) - 0 + 2\cdot (\frac{3}{2} \times \frac{4}{3}) = \lambda^3 - \frac{2}{3}\lambda + 4
  $$

  特征方程
  $$
  \lambda^3 - \frac{2}{3}\lambda + 4 = 0
  $$
  可以尝试取 $\lambda = -2, -1$
  $$
  f(-1) = 11 > 0, f(-2) = -8 < 0
  $$
  由于 $f(-2) < 0 < f(-1)$，由介值定理可知在 $[-2, -1]$ 范围内存在 $\lambda_1$ 使得特征方程为 $0$，且
  $$
  |\lambda_1| > 1
  $$
  因此
  $$
  \rho(B) = \max_{i} |\lambda_i| \ge |\lambda_1| > 1
  $$

由此得出结论，**Jacobi 迭代法对该矩阵 $A$ 是发散的**（不收敛）.

---



## G-S 迭代法的敛散性判断

### G-S 迭代法的敛散性

和之前的 Jacobi 迭代法一样的$Ax = b$ 可以表示成：
$$
(D - L - U)x = b
$$
改写为迭代格式：
与前面 Jacobi 迭代法 不同的是，G-S 法会用到之前更新后的数值。

因此，在矩阵形式的迭代公式中，第 $k+1$ 次迭代时，用“最新”的分量，所以把左边看成矩阵 $D-L$ 乘以 $x^{(k+1)}$，右边只用旧迭代 $x^{(k)}$：
$$
(D - L)x = Ux + b
$$

$$
x^{(k + 1)} = (D-L)^{-1}Ux^{(k)} + c
$$

其中 $c = (D - L)^{-1}b$ .

可以写成标准的迭代格式：
$$
x^{(k + 1)} = B_{GS}Ux^{(k)} + c
$$
其中 $B_{GS} = (D - L)^{-1}U$ .

要判断 G-S 迭代法是否收敛（G-S 法对任意初值收敛），即判断 **矩阵 $B_{GS}$ 的谱半径是否 $\rho(B_{GS}) < 1$** .

**方法总结**：

1. 把 $A$ 分解为 $D - L - U$
2. 写出迭代格式进一步得到 G-S 迭代矩阵 $B_{GS} = (D - L)^{-1}U$
3. 求 $B_{GS}$ 的谱半径 $\rho(B_{GS})$
   - 若 $\rho(B_{GS}) < 1$，则 G-S 迭代收敛
   - 否则，不收敛


### 例题 2

> $A = \begin{pmatrix} 1 & -2 & 2 \\ -1 & 1 & -1 \\ -2 & -2 & 1 \end{pmatrix}$，判断 G-S 迭代法对该矩阵是否收敛。

- 将 $A$ 分解为 $A = D - L - U$
  $$
  D = \begin{pmatrix}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1 
  \end{pmatrix}, \quad 
  L = \begin{pmatrix}
  0 & 0 & 0 \\
  1 & 0 & 0 \\
  2 & 2 & 0 
  \end{pmatrix}, \quad 
  U = \begin{pmatrix}
  0 & 2 & -2 \\
  0 & 0 & 1 \\
  0 & 0 & 0
  \end{pmatrix}
  $$

- G-S 迭代矩阵 $B_{GS}$
  $$
  B = (D - L)^{-1}U 
  $$
  先求 $(D - L)^{-1}$：
  $$
  D - L = \begin{pmatrix}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1 
  \end{pmatrix}
  - \begin{pmatrix}
  0 & 0 & 0 \\
  1 & 0 & 0 \\
  2 & 2 & 0 
  \end{pmatrix} 
  = \begin{pmatrix}
  1 & 0 & 0 \\
  -1 & 1 & 0 \\
  -2 & -2 & 1 
  \end{pmatrix}
  $$

  $$
  (D - L)^{-1} = \begin{pmatrix}
  1 & 0 & 0 \\
  1 & 1 & 0 \\
  4 & 2 & 1
  \end{pmatrix}
  $$

  再求 $(D - L)^{-1} U$：
  $$
  B = (D - L)^{-1}U = \begin{pmatrix}
  1 & 0 & 0 \\
  1 & 1 & 0 \\
  4 & 2 & 1
  \end{pmatrix}
  \begin{pmatrix}
  0 & 2 & -2 \\
  0 & 0 & 1 \\
  0 & 0 & 0
  \end{pmatrix}
  = \begin{pmatrix}
  0 & 2 & -2 \\
  0 & 2 & -1 \\
  0 & 8 & -6
  \end{pmatrix}
  $$
  

- 求 $B$ 的谱半径判断收敛性
  $$
  \lambda I - B = \begin{pmatrix}
  \lambda & -2 & 2 \\
  0 & \lambda - 2 & 1 \\
  0 & -8 & \lambda + 6
  \end{pmatrix}
  $$

  $$
  \begin{align*}
  \det(\lambda I - B) &= \lambda \cdot [(\lambda - 2)(\lambda + 6) + 8] - 0 + 0 \\
  &= \lambda \cdot (\lambda^2 + 4 \lambda - 4)
  \end{align*}
  $$

  特征方程
  $$
  \lambda \cdot (\lambda^2 + 4 \lambda - 4) = 0
  $$
  求得特征值为：
  $$
  \lambda_1 = 0, \quad \lambda_2 = -2 + 2\sqrt 2, \quad \lambda_3 = -2 - 2\sqrt 2
  $$
  因此
  $$
  \rho(B_{GS}) = \max_{i}|\lambda_i| = |\lambda_3| = 2 + 2\sqrt 2 > 1
  $$

由此得出结论，**G-S 迭代法对该矩阵 $A$ 是发散的**（不收敛）.

---

