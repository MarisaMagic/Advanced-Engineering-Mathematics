<h1>
    <div align='center'>
        第五章 最优化方法（4）
        KKT点问题，二阶充分条件 
    </div>
</h1>

## 带约束规划问题的最优性条件（求 KKT 点问题）

### 标准问题形式 & Lagrange 函数

该类型优化问题的标准形：
$$
\begin{aligned}
\min \quad & f(x) \\
\text{s.t.} \quad & h_i(x)=0,\ i\in \mathcal E \\
& g_j(x)\ge 0,\ j\in \mathcal I.
\end{aligned}
$$
其中 $h_i(x)=0,\ i\in E; \ g_j(x)\ge 0,\ j\in I.$ 是约束的可行域。通常考试只考带有 $g_j(x)\ge 0,\ j\in I$ 的情形。

**定义 Lagrange 函数**：
$$
L(x,w,v)=f(x)-\sum_{i\in E}w_i h_i(x)-\sum_{j\in I}v_j g_j(x)
$$
其中 $w_i$ 是等式约束的拉格朗日乘子，$v_j$ 是不等式约束的拉格朗日乘子。
### KKT 一阶必要条件

在满足一定约束资格条件（线性独立性条件）的前提下，若 $x^*$ 是局部极小点，则存在乘子向量 $w^*,v^*$，使得：

1. **梯度条件**
   $$
   \nabla_x L(x^*,w^*,v^*)=0
   $$
   即： Lagrange 函数对所有变量求偏导并令其为 0。

2. **可行性条件**
   $$
   h_i(x^*)=0,\ i\in \mathcal E \\
   g_j(x^*)\ge 0,\ j\in \mathcal I.
   $$
   也就是问题本身自带的可行域约束

3. **互补松弛条件**
   $$
   v_j^*\ge 0,\ j\in \mathcal I. \\
   v_j^* g_j(x^*)=0,\ j\in \mathcal I
   $$
   包括拉格朗日乘子 $v_j$ 非负，并且和对应约束 $g_j(x)$ 互补松弛。

   - 若约束 $g_j(x)$ 不活跃（$g_j(x) > 0$，严格大于 0），则对应乘子 $v_j$ 必须为 0；
   - 若约束 $g_j(x)$ 活跃（$g_j(x) = 0$，即在边界），则乘子 $v_j$ 可非零。

总结可得 KKT 条件如下（适用于本校考试）：
$$
\begin{cases}
\nabla_x L(x, w, v) = 0, \\
h_i(x) = 0, \ i \in \mathcal E, \\
g_j(x) \ge 0, \ j \in \mathcal I, \\
v_j \ge 0, \ j \in \mathcal I, \\
v_j g_j(x)=0,\ j\in \mathcal I
\end{cases}
$$

---



## 二阶充分条件

**二阶充分条件**：

设 $x^*$ 是问题 的 KKT 点，对应乘子 $w^*,v^*$，且 $v_j^*\ge 0$。
若对所有 **非零** 方向 $d\in G(x^*,w^*,v^*)$，都有
$$
d^T\nabla_{xx}^2L(x^*,w^*,v^*)d>0
$$
则 $x^*$ 是问题的 **严格局部极小点**。

其中的方向集合 $G(x^*,w^*,v^*)$ 定义为：

$$
G(x^*,w^*,v^*) = \left\{ d\ \middle|\
\begin{aligned}
&\nabla h_i(x^*)^T d = 0,\quad i\in \mathcal E,\\
&\nabla g_j(x^*)^T d = 0,\quad j\in \mathcal I(x^*),\ v_j^*>0,\\
&\nabla g_j(x^*)^T d \ge 0,\quad j\in \mathcal I(x^*),\ v_j^*=0
\end{aligned} 
\right\}
$$

* $\mathcal I(x^*)={j\in \mathcal I\mid g_j(x^*)=0}$：**活跃的不等式约束集合**；
* 第一行：沿等式约束切面移动；
* 第二、三行：对活跃不等式约束，根据乘子是否 >0 再细分。



结合 KKT 点的求解和上述二阶充分条件，完整的解题步骤如下：
**Step 0：写成标准形式**

把原题整理成：
$$
\begin{aligned}
\min\quad & f(x) \\
\text{s.t. }\quad & h_i(x)=0,\\
& g_j(x)\ge0.
\end{aligned}
$$
**Step 1：引入乘子，写 Lagrange 函数**
$$
L(x,w,v)=f(x)-\sum w_i h_i(x)-\sum v_j g_j(x)
$$
**Step 2：列 KKT 条件并求解**

$$
\begin{cases}
\nabla_x L(x,w,v)=0, &\text{(梯度条件)}\\
h_i(x)=0, &\text{(等式可行)}\\
g_j(x)\ge0, &\text{(不等式可行)}\\
v_j g_j(x)=0,\ v_j\ge0,&\text{(互补松弛)}
\end{cases}
$$

**分类讨论活跃约束**，解出所有 $(x^*,w^*,v^*)$。得到的每个 $x^*$ 都是 **候选最优点**，接下来用二阶条件筛选。

**Step 3：对每个候选点构造 $G(x^*,w^*,v^*)$**

1. 先确定活跃集合：
   $$
   \mathcal I(x^*)={j\mid g_j(x^*)=0}
   $$

2. 根据乘子 $v_j^*$ 再细分：

   * 如果 $j\in I(x^*)$ 且 $v_j^*>0$，要求
     $$
     \nabla g_j(x^*)^T d = 0;
     $$
   * 如果 $j\in I(x^*)$ 且 $v_j^*=0$，要求
     $$
     \nabla g_j(x^*)^T d \ge 0.
     $$

3. 再加上所有等式约束的方向条件（不过考试一般没有这个）：
   $$
   \nabla h_i(x^*)^T d=0,\ i\in E.
   $$

**Step 4：计算二阶导（Hessian）**
$$
H=\nabla_{xx}^2L(x^*,w^*,v^*)
$$

注意：若 $h_i,g_j$ 是线性的，则它们二阶导为 0，$H$ 就等于 $\nabla^2 f(x^*)$。

**Step 5：检查二阶充分条件**

在所有非零 $d\in G(x^*,w^*,v^*)$ 上检查：

$$
d^T H d>0,\quad \forall, d\in G(x^*,w^*,v^*),\ d\neq0.
$$

* 若成立，则 $x^*$ 是问题的**严格局部极小点**；
* 若不成立：这个 KKT 点不一定是极小点，可能是鞍点或别的。

---



## 例题

> 求下列约束优化问题的局部极小值点
> $$
> \begin{aligned}
> \min\quad & (x_1-1)^2 + x_2\\
> \text{s.t. }\quad & -x_1 - x_2 + 2 \ge 0,\\
> & x_2 \ge 0.
> \end{aligned}
> $$

解答：

**Step 0：写成标准形式**

设：

$$
f(x)= (x_1-1)^2 + x_2,
\quad
g_1(x)=-x_1-x_2+2\ge0,\quad
g_2(x)=x_2\ge0.
$$

本题无等式约束。


**Step 1：写 Lagrange 函数**
$$
L(x,v)=f(x)-v_1g_1(x)-v_2g_2(x)
$$

即：

$$
L=(x_1-1)^2+x_2 - v_1(-x_1-x_2+2)-v_2x_2
$$

**Step 2：列 KKT 条件并求解**

**1. 梯度条件：**

$$
\nabla_x L=0\Rightarrow
\begin{cases}
2(x_1-1)+v_1=0\\
1+v_1-v_2=0
\end{cases}
\tag{A}
$$

**2. 可行性条件：**

$$
-x_1-x_2+2\ge0,\quad x_2\ge0
\tag{B}
$$

**3. 互补松弛条件：**

$$
v_1(-x_1-x_2+2)=0,\quad v_2x_2=0,\quad v_1,v_2\ge0
\tag{C}
$$
（此处省略分类讨论，此处本应需要 4 中情况分类讨论，每个约束为活跃或非活跃）
代入 $x = (1, 0)^T$ 发现满足 KKT 条件，且其它情况不满足。
故 KKT 候选点 $x^{*} = (1, 0)^T$，对应乘子 $v_1^*=0,\ v_2^*=1$ .



**Step 3：构造 $G(x^*,v^*)$**

首先确定 **活跃不等式集合**：

$$
\mathcal I(x^*)=\{j \mid g_j(x^*)=0\}=\{2\}
$$

也就是第 2 个条件 $g_2(x)$ 是活跃的。
因为 $v_2^*=1>0$，根据规则：
$$
\nabla g_2(x^*)^Td = 0 \Rightarrow (0,1)^T \cdot (d_1,d_2)=d_2=0
$$
由于 $d$ 是 **非零** 方向向量，因此：
$$
G(x^*,v^*)= \{d\mid d=(d_1,0)\neq(0,0)\}
$$

**Step 4：计算 Hessian**
$$
H=\nabla^2_{xx}L=\nabla^2 f(x)=
\begin{pmatrix}
2 & 0\\
0 & 0
\end{pmatrix}
$$

因为 $g_1,g_2$ 是线性函数，其二阶导为 0，因此 Hessian 只由 $f$ 决定。


**Step 5：检查二阶充分条件**

对所有 $d\in G(x^*,v^*)$，$d=(d_1,0)\neq0$：

$$
d^T H d = (d_1,0)
\begin{pmatrix}
2 & 0\\
0 & 0
\end{pmatrix}
\begin{pmatrix}
d_1 \\ 0
\end{pmatrix}
=2d_1^2>0
$$

**符合二阶充分条件**。


**综上**：
$$
\boxed{x^*=(1,0)\text{ 是该问题的严格局部极小点（并且是全局最优点）。}}
$$

---

