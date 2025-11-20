<h1>
    <div align='center'>
        第三章 矩阵分解与广义逆矩阵（1）
        LU分解，PLU分解，满秩分解
    </div>
</h1>

## LU 分解（Doolittle 型）

**Doolittle 型 LU 分解**：将矩阵 $A$ 分解为
$$
A = LU
$$
其中

- $L$：**下三角矩阵**（通常主对角线全为 $1$）
- $U$：**上三角矩阵**

通常应用于求解方程 $Ax = b$，经过 $A = LU$ 分解后，可以先求解 $LUx = Ly = b$，然后再求解 $Ux = y$（回代），能够更加容易求解方程。

**一般要求**：所有**前 k 阶主子式（左上角 k×k 行列式）非零**，简单说就是消元时每一步主元都不为 0，就能顺利做 LU 而不用交换行。



**LU 分解基本过程**：

- 经过一系列初等行变换（**不能交换行**）得到一个上三角矩阵 $U$；
- 可以看出下三角矩阵 $L$：
  - 首先主对角线元素全为 $1$；
  - 第 1 列元素比例（从对角线元素 $1$ 开始）对应 初始矩阵 $A$ 的第 1 列元素的比例；
  - 第 2 列元素比例（从对角线元素 $1$ 开始）对应 经过初等行变换后第一列主元下方均为 $0$ 时的矩阵的第 2 列元素的比例
  - 以此类推，得到下三角矩阵 $L$。

**示例**：

> 求矩阵 $A = \begin{pmatrix} 2 & 5 & -6 \\ 4 & 13 & -19 \\ -6 & -3 & -6 \end{pmatrix}$ 的 LU 分解

- 用初等行变换将第 1 列主元下方全部变为 $0$
  $$
  \begin{pmatrix}
  2 & 5 & -6 \\ 
  4 & 13 & -19 \\ 
  -6 & -3 & -6
  \end{pmatrix}
  \xRightarrow{r_2 - 2r_1, \ r_3 + 3r1}
  \begin{pmatrix}
  2 & 5 & 6 \\
  0 & 3 & -7 \\
  0 & 12 & -24
  \end{pmatrix}
  $$

- 用初等行变换将第 2 列主元下方全部变为 $0$
  $$
  \begin{pmatrix}
  2 & 5 & 6 \\
  0 & 3 & -7 \\
  0 & 12 & -24
  \end{pmatrix}
  \xRightarrow{r_3 - 4r2}
  \begin{pmatrix}
  2 & 5 & 6 \\
  0 & 3 & -7 \\
  0 & 0 & 4
  \end{pmatrix}
  $$

此时得到的矩阵是一个上三角矩阵，即 $U$：
$$
U = \begin{pmatrix}
2 & 5 & 6 \\
0 & 3 & -7 \\
0 & 0 & 4
\end{pmatrix}
$$
变换之前，第 1 列主元及下方元素比例为 $[1, 2, -3]$；第 1 次变换后，第 2 列主元及下方元素比例为 $[1, 4]$。故下三角矩阵 $L$ 为：
$$
L = \begin{pmatrix}
1 & 0 & 0 \\
2 & 1 & 0 \\
-3 & 4 & 1
\end{pmatrix}
$$
 最终得到 LU 分解：
$$
A = LU =  \begin{pmatrix}
1 & 0 & 0 \\
2 & 1 & 0 \\
-3 & 4 & 1
\end{pmatrix}
\begin{pmatrix}
2 & 5 & 6 \\
0 & 3 & -7 \\
0 & 0 & 4
\end{pmatrix}
$$

---



## PLU 分解（Doolittle 型）

PLU 分解和 LU 分解的关键区别就在于 PLU 分解是带有行交换的：
$$
PA = LU
$$
其中

- $P$：**置换矩阵**（从单位矩阵通过交换行得到，代表进行的行交换操作）
- $L$：**下三角矩阵**（通常主对角线全为 $1$）
- $U$：**上三角矩阵**

**为什么需要 PLU 分解**：在进行初等行变换时，可能遇到某一步主元为 $0$。解决方法就是，每一次初等行变换先“选择主元”，通常的做法是将当前列最大元素所在行交换到主元所在行。这些行交换操作最终形成一个置换矩阵 $P$。



**示例**：

> 求矩阵 $A = \begin{pmatrix} 0.5 & 1 & 0 \\ 2 & 1.5 & 1 \\ 0.2 & 1 & 2.5 \end{pmatrix}$ 的 PLU 分解

- 交换行，并将第 1 列主元下方元素变为 $0$
  $$
  \begin{pmatrix} 
  0.5 & 1 & 0 \\ 
  2 & 1.5 & 1 \\ 
  0.2 & 1 & 2.5 
  \end{pmatrix}
  \xRightarrow{r_1\leftrightarrow r_2}
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0.5 & 1 & 0 \\ 
  0.2 & 1 & 2.5 
  \end{pmatrix}
  $$

  $$
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0.5 & 1 & 0 \\ 
  0.2 & 1 & 2.5 
  \end{pmatrix}
  \xRightarrow{r_2 - \frac{1}{4}r_1, \ r_3 - \frac{1}{10}r_1}
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0 & 0.625 & -0.25 \\ 
  0 & 0.85 & 2.4 
  \end{pmatrix}
  $$

- 交换行，并将第 2 列主元下方元素变为 $0$
  $$
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0 & 0.625 & -0.25 \\ 
  0 & 0.85 & 2.4 
  \end{pmatrix}
  \xRightarrow{r_2 \leftrightarrow r_3}
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0 & 0.85 & 2.4 \\
  0 & 0.625 & -0.25
  \end{pmatrix}
  $$

  $$
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0 & 0.85 & 2.4 \\
  0 & 0.625 & -0.25
  \end{pmatrix}
  \xRightarrow{r_3 - \frac{0.625}{0.85}r_2}
  \begin{pmatrix} 
  2 & 1.5 & 1 \\
  0 & 0.85 & 2.4 \\
  0 & 0 & -2.0147
  \end{pmatrix}
  $$

此时得到的矩阵是一个上三角矩阵，即 $U$：
$$
U = \begin{pmatrix} 
2 & 1.5 & 1 \\
0 & 0.85 & 2.4 \\
0 & 0 & -2.0147
\end{pmatrix}
$$
可以得到下三角矩阵 $L$（注意此处应该遵从多次行交换后的最终的比例）：
$$
L = \begin{pmatrix}
1 & 0 & 0 \\
0.1 & 1 & 0 \\
0.25 & 0.7353 & 1 
\end{pmatrix}
$$
行变换矩阵 $P$：
$$
I = \begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{pmatrix}
\xRightarrow{r_1 \leftrightarrow r_2} 
\begin{pmatrix}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1
\end{pmatrix}
\xRightarrow{r_2 \leftrightarrow r_3} 
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{pmatrix} = P
$$
最终得到 PLU 分解：
$$
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{pmatrix} A 
=
\begin{pmatrix}
1 & 0 & 0 \\
0.1 & 1 & 0 \\
0.25 & 0.7353 & 1 
\end{pmatrix}
\begin{pmatrix} 
2 & 1.5 & 1 \\
0 & 0.85 & 2.4 \\
0 & 0 & -2.0147
\end{pmatrix}
$$

---



## 满秩分解

把一个矩阵 $A$ 写成“两个**满秩**矩阵相乘”的形式，用这两个小一点的矩阵来体现它的秩结构。

设 $A$ 是一个 $m\times n$ 矩阵，$\operatorname{rank}(A)=r$。
若存在矩阵
$$
A = F G
$$
其中

* $F$ 是 $m\times r$ 矩阵，且 $\operatorname{rank}(F)=r$（列满秩）
* $G$ 是 $r\times n$ 矩阵，且 $\operatorname{rank}(G)=r$（行满秩）

那么就称 $A = F G$ 是 $A$ 的**满秩分解**。



**满秩分解的基本过程**：

- 对矩阵 $A$ 进行初等行变换（不要进行交换行），得到矩阵的 **行最简阶梯形**；

  **行最简阶梯形** 的条件：

  - 零行（全是 $0$ 的行）排在最下面；
  - 每个非零行的第一个非零元是 1。这个 1 叫做该行的 **主元**；
  - 主元所在的那一列，其他行都是 $0$；
  - 从上到下看，每一行的主元都在上一行主元的右边（呈阶梯型）

- 取主元在 **原矩阵** $A$ 中对应的 **列** 构成 $F$ ；

- 取主元在 **行最简阶梯形** 中对应的 **行** 构成 $G$ ；

- 最终得到 $A = BC$.

**示例**：

> 求矩阵 $A = \begin{pmatrix} 1 & -1 & -1 & 4 \\ 0 & 0 & 2 & 2 \\ -1 & 1 & 5 & 0 \end{pmatrix}$ 的满秩分解

- 初等行变换得到 **行最简阶梯型**
  $$
  \begin{pmatrix} 
  1 & -1 & -1 & 4 \\ 
  0 & 0 & 2 & 2 \\ 
  -1 & 1 & 5 & 0 
  \end{pmatrix}
  \xRightarrow{r_3 + r_1} 
  \begin{pmatrix} 
  1 & -1 & -1 & 4 \\ 
  0 & 0 & 2 & 2 \\ 
  0 & 0 & 4 & 4 
  \end{pmatrix}
  \xRightarrow{r_3 - 2r_2} 
  \begin{pmatrix} 
  1 & -1 & -1 & 4 \\ 
  0 & 0 & 2 & 2 \\ 
  0 & 0 & 0 & 0 
  \end{pmatrix}
  $$

  $$
  \begin{pmatrix} 
  1 & -1 & -1 & 4 \\ 
  0 & 0 & 2 & 2 \\ 
  0 & 0 & 0 & 0 
  \end{pmatrix}
  \xRightarrow{r_2 \times \frac{1}{2}}
  \begin{pmatrix} 
  1 & -1 & -1 & 4 \\ 
  0 & 0 & 1 & 1 \\ 
  0 & 0 & 0 & 0 
  \end{pmatrix}
  \xRightarrow{r_1 + r_2}
  \begin{pmatrix} 
  \bf{1} & -1 & 0 & 5 \\ 
  0 & 0 & \bf{1} & 1 \\ 
  0 & 0 & 0 & 0 
  \end{pmatrix}
  $$

- 取主元在 **原矩阵** 中对应的列构成 $F$：
  $$
  F = \begin{pmatrix}
  1 & -1 \\
  0 & 2 \\
  -1 & 5
  \end{pmatrix}
  $$

- 取主元在 **行最简阶梯形** 中对应的行构成 $G$：
  $$
  G = \begin{pmatrix}
  1 & -1 & 0 & 5 \\
  0 & 0 & 1 & 1
  \end{pmatrix}
  $$

最终得到满秩分解：
$$
A = FG = \begin{pmatrix}
1 & -1 \\
0 & 2 \\
-1 & 5
\end{pmatrix}
\begin{pmatrix}
1 & -1 & 0 & 5 \\
0 & 0 & 1 & 1
\end{pmatrix}
$$

---

