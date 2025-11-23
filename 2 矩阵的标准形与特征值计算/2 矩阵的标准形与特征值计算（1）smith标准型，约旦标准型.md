<h1>
    <div align='center'>
        第二章 矩阵的标准形与特征值计算（1）
        Smith标准型，Jordan 标准型
    </div>
</h1>

## Smith 标准型

### Smith 标准型定义

对于一个整数矩阵 $A$，进行一系列初等行列变换可以转换为一个对角矩阵（$P$ 代表行变换，$Q$ 代表列变换）：
$$
PAQ = \begin{pmatrix}
d_1 & 0   & \cdots & 0 \\
0   & d_2 & \cdots & 0 \\
\vdots & & \ddots & \vdots \\
0   & 0   & \cdots & d_r \\
&&&& 0 \\
&&&&& \ddots
\end{pmatrix}
$$
满足：

- 每个 $d_i \not = 0$ ；
- $d_1 \mid d_2,\ d_2 \mid d_3,\ \dots,\ d_{r-1} \mid d_r$ （$d_1 \mid d_2$ 表示 $d_1$ 整除 $d_2$）；
- $d_k$ 是矩阵 $A$ 的所有 $k$ 阶余子式的最大公因子。

这个对角矩阵就是 $A$ 的 **Smith 标准型**。

### 不变因子

Smith 标准型中**非零对角线上的数**：
$$
d_1,\ d_2,\ \dots,\ d_r,\quad (d_1 \mid d_2 \mid \cdots\mid d_r)
$$
称为矩阵 $A$ 的 **不变因子**。

### 初等因子

将每个不变因子 $d_k$ 进行分解，然后提取出所有一次幂。

例如：
$$
d_1 = \lambda, \quad d_2 = \lambda^3 + \lambda = \lambda(\lambda^2 + 1)
$$
得到所有一次幂：
$$
\lambda, \ \lambda, \ \lambda + i, \ \lambda - i
$$
这些一次幂就是矩阵 $A$ 的 **初等因子**。

### 行列式因子

设 $A$ 是一个 $m\times n$ 的整数矩阵。对每个 $k=1,2,\dots,r$（其中 $r=\text{rank}(A)$），定义：

$$
\Delta_k = \gcd{\text{所有 }k\times k\text{ 子式的行列式}}.
$$

* $\Delta_1$：所有元素（1×1 子式）的 gcd
* $\Delta_2$：所有 2×2 子矩阵的行列式的 gcd
* …
* $\Delta_r$：所有 $r\times r$ 子式的行列式的 gcd（满秩时就是 $|\det A|$）

这串 $\Delta_k$ 就叫做矩阵的 **行列式因子序列**。

**与不变因子的关系**：

如果矩阵的 Smith 标准型对角元（不变因子）是

$$
d_1, d_2, \dots, d_r \quad (d_1\mid d_2\mid \cdots\mid d_r),
$$
那么有非常重要的关系：

$$
\Delta_k = d_1 d_2 \cdots d_k \quad (k=1,\dots,r).
$$
于是可以反过来求不变因子：

$$
d_1 = \Delta_1,\quad
d_2 = \frac{\Delta_2}{\Delta_1},\quad
d_3 = \frac{\Delta_3}{\Delta_2},\ \dots
$$
这就是所谓的 **gcd 法**/行列式因子法。

---



## 不变因子和初等因子求法

### 初等变换法

**步骤**：

1. 对矩阵 $A$ 进行一系列初等变换，转换为对角阵（注意最后还要满足相关性质）；
2. 对角阵上的非 0 元素即 **不变因子** ；
3. 将不变因子进行质因数分解，所有质数幂即 **初等因子**。

**初等变换的要求**：

* 交换两行 / 两列
* 一行（列）加「任意多项式倍」的另一行（列），比如 $R_i \leftarrow R_i + f(\lambda)R_j$
* 一行（列）乘以 $-1$



**示例**：

矩阵 $A$ 为：
$$
A = \begin{pmatrix}
\lambda & 0 & 0 \\
0 & \lambda & 0 \\
1 & 1 & \lambda
\end{pmatrix}
$$
初等变换：

- 把 1 移到左上角

  交换第 1、3 行： $R_1 \leftrightarrow R_3$

  $$
  \begin{pmatrix}
  1 & 1 & \lambda\\
  0 & \lambda & 0\\
  \lambda & 0 & 0
  \end{pmatrix}
  $$

- 消去第一列下面的 $\lambda$

  $R_3 \leftarrow R_3 - \lambda R_1$（加上“$-\lambda$ 倍”的第 1 行）

  $$
  \begin{pmatrix}
  1 & 1 & \lambda\\
  0 & \lambda & 0\\
  0 & -\lambda & -\lambda^2
  \end{pmatrix}
  $$

- 把右边整理成上三角

  $C_2 \leftarrow C_2 - C_1$

  $$
  \begin{pmatrix}
  1 & 0 & \lambda\\
  0 & \lambda & 0\\
  0 & -\lambda & -\lambda^2
  \end{pmatrix}
  $$

  $C_3 \leftarrow C_3 - \lambda C_1$

  $$
  \begin{pmatrix}
  1 & 0 & 0\\
  0 & \lambda & 0\\
  0 & -\lambda & -\lambda^2
  \end{pmatrix}
  $$

- 处理右下 $2\times2$ 子块

  先 $R_3 \leftarrow -R_3$：

  $$
  \begin{pmatrix}
  1 & 0 & 0\\
  0 & \lambda & 0\\
  0 & \lambda & \lambda^2
  \end{pmatrix}
  $$

  再 $R_3 \leftarrow R_3 - R_2$：

  $$
  \begin{pmatrix}
  1 & 0 & 0\\
  0 & \lambda & 0\\
  0 & 0 & \lambda^2
  \end{pmatrix}
  $$

此时已经是对角，且满足 $1 \mid \lambda \mid \lambda^2$，由此得到 **Smith 标准型**。

**不变因子**：
$$
d_1 = 1,\quad d_2 = \lambda,\quad d_3 = \lambda^2
$$
**初等因子**：
$$
\boxed{\lambda,\ \lambda,\ \lambda} \quad(\lambda 出现三次，带重数).
$$


### gcd 法（行列式因子法）

**步骤**：

1. 先求出 $\Delta_1$，即所有 1 阶余子式的最大公因式。由此得到 $d_1 = \Delta_1$
2. 求出 $\Delta_2$，即所有 2 阶余子式的最大公因式。然后可以得到 $d_2 = \frac{\Delta_2}{\Delta_1}$
3. 以此类推，求得 $d_3, \dots, d_r$ .



**示例**：

- $\Delta_1$：所有 1×1 子式（元素本身）的 gcd

  矩阵元素有 $\lambda, 0, 1$ 等，显然
  $$
  \Delta_1 = \gcd(\lambda,0,1,\dots)=1.
  $$
  所以
  $$
  d_1 = \Delta_1 = 1.
  $$

- $\Delta_2$：所有 2×2 子式行列式的 gcd

  列举代表性的 2×2 子式（行列式）：

  * 取第 1、2 行与第 1、2 列：
    $$
    \begin{vmatrix}
    \lambda & 0\\
    0 & \lambda
    \end{vmatrix}
    = \lambda^2
    $$

  * 取第 1、3 行与第 1、2 列：
    $$
    \begin{vmatrix}
    \lambda & 0\\
    1 & 1
    \end{vmatrix}
    = \lambda\cdot 1 - 0\cdot 1 = \lambda
    $$

  * 取第 2、3 行与第 1、2 列：
    $$
    \begin{vmatrix}
    0 & \lambda\\
    1 & 1
    \end{vmatrix}
    = 0\cdot 1 - \lambda\cdot 1 = -\lambda
    $$

  由 ${\lambda^2, \lambda, -\lambda, 0,\dots}$ 的最大公因子可知
  $$
  \Delta_2 = \gcd(\lambda^2,\lambda,-\lambda) = \lambda.
  $$

  于是
  $$
  d_2 = \frac{\Delta_2}{\Delta_1} = \frac{\lambda}{1} = \lambda.
  $$

- $\Delta_3$：所有 3×3 子式 —— 就是 $\det A$

  沿第一行展开行列式：

  $$
  \det A =
  \lambda \cdot
  \begin{vmatrix}
  \lambda & 0\\
  1 & \lambda
  \end{vmatrix}
  - 0 + 0
  = \lambda (\lambda^2 - 0)
  = \lambda^3.
  $$

  所以
  $$
  \Delta_3 = \det A = \lambda^3,
  \quad d_3 = \frac{\Delta_3}{\Delta_2} = \frac{\lambda^3}{\lambda} = \lambda^2.
  $$

**不变因子**：
$$
d_1 = 1,\quad d_2 = \lambda,\quad d_3 = \lambda^2
$$
**初等因子**：
$$
\boxed{\lambda,\ \lambda,\ \lambda} \quad(\lambda 出现三次，带重数).
$$

---



## Jordan 标准型

### Jordan 标准型定义

Jordan 标准型由若干个 **Jordan 块** 沿对角线拼接而成。一个大小为 $k \times k$，对应特征值为 $\lambda$ 的 Jordan 块形如：
$$
J_k(\lambda) =
\begin{pmatrix}
\lambda & 1      &        &        & 0 \\
& \lambda & 1      &        &   \\
&        & \ddots & \ddots &   \\
&        &        & \lambda & 1 \\
0       &        &        &        & \lambda
\end{pmatrix}
$$
特点：

* 主对角线上全是同一个特征值 $\lambda$；
* 紧靠主对角线的上一条对角线是 1；
* 其它地方全是 0。



可以通过求 $(\lambda I - A)$ 的初等因子得到矩阵 $A$ 的 Jordan 标准型：

> 初等因子 → 分解出每个特征值对应的 Jordan 块大小 → 写出 Jordan 标准型。



### 初等因子与 Jordan 块的关系

在一个代数闭域（比如复数域 $\mathbb{C}$）中，一个矩阵 $A$ 的初等因子都形如：

$$
(\lambda - \lambda_i)^{k}
$$
每个这样的因子与 Jordan 块对应：

* 每一个初等因子 → 一个 Jordan 块；
* 指数 $k$ → 这个 Jordan 块的大小（阶数）。

举例：

| 初等因子        | Jordan 块                                      |
| --------------- | ---------------------------------------------- |
| $(\lambda-3)^1$ | $J_1(3) = [3]$                                 |
| $(\lambda-3)^2$ | $J_2(3)=\begin{pmatrix}3&1\\ 0&3\end{pmatrix}$ |
| $(\lambda-5)^4$ | $J_4(5)$（4×4的大块）                          |



### 例题 1

> $A = \begin{pmatrix} 3 & 1 & -3 \\ -7 & -2 & 9 \\ -2 & -1 & 4 \end{pmatrix}$，求这个矩阵 $A$ 的 Jordan 标准型 $J$。

解答：

- 求 $\lambda I - A$ 
  $$
  \lambda I - A = \begin{pmatrix}
  \lambda - 3 & -1 & 3 \\
  7 & \lambda + 2 & -9 \\
  2 & 1 & \lambda - 4
  \end{pmatrix}
  $$

- 此处考虑用 gcd 法求 $\lambda I - A$ 的初等因子
  $$
  d_1 = \Delta_1 = \gcd(\lambda - 3, -1, 3, \dots) = 1
  $$
  列举代表性的 2×2 子式（行列式）：

  - 第 1 行，第 3 列

  $$
  \begin{vmatrix}
  7 & \lambda + 2 \\
  2 & 1
  \end{vmatrix}
  = 7 - 2(\lambda + 2) = -2\lambda + 3
  $$

  - 第 3 行，第 1 列
    $$
    \begin{vmatrix}
    -1 & 3 \\
    \lambda + 2 & -9
    \end{vmatrix}
    = 9 - 3(\lambda + 2) = -3\lambda + 3 = -3(\lambda - 1)
    $$

  可以发现，$\gcd(-2\lambda + 3, -3(\lambda - 1)) = 1$。

  故：
  $$
  \Delta_2 = \gcd(-2\lambda + 3, -3(\lambda - 1), \dots) = 1
  $$

  $$
  d_2 = \frac{\Delta_2}{\Delta_1} = 1
  $$

  求 3×3 子式（行列式）：
  
  $$
  \begin{align*}
  \det(\lambda I - A) &= 
  (\lambda - 3) \cdot \begin{vmatrix} \lambda + 2 & -9 \\ 1 & \lambda - 4 \end{vmatrix}
  -
  (-1) \cdot \begin{vmatrix} 7 & -9 \\ 2 & \lambda - 4 \end{vmatrix}
  +
  3 \cdot \begin{vmatrix} 7 & \lambda + 2 \\ 2 & 1 \end{vmatrix} \\
  &= (\lambda - 3)(\lambda^2 - 2\lambda + 1) + (7\lambda - 10) + 3(3 - 2\lambda) \\
  &= \lambda^3 - 5\lambda^2 + 8\lambda - 4
  \end{align*}
  $$

  故：

  $$
  \Delta_3 = \det(\lambda I - A) = \lambda^3 - 5\lambda^2 + 8\lambda - 4 = (\lambda - 2)^2(\lambda - 1)
  $$

  $$
  d_3 = \frac{\Delta_3}{\Delta_2} = (\lambda - 2)^2(\lambda - 1)
  $$

  **不变因子**：
  $$
  d_1 = 1, \ d_2 = 1, \ d_3 = (\lambda - 2)^2(\lambda - 1)
  $$
  **初等因子**：
  $$
  (\lambda - 2)^2, \quad \lambda - 1
  $$

- 初等因子对应的 Jordan 块
  $$
  \begin{pmatrix}
  2 & 1 \\
  0 & 2
  \end{pmatrix}, \quad 
  \begin{pmatrix}
  1
  \end{pmatrix}
  $$

- 最终得到 Jordan 标准型
  $$
  \boxed{
  J = \begin{pmatrix}
  2 & 1 & 0 \\
  0 & 2 & 0 \\
  0 & 0 & 1
  \end{pmatrix}
  }
  $$