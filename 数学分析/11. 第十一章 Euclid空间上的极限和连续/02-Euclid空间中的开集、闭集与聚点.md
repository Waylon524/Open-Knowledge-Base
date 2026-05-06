# 02. Euclid 空间中的开集、闭集与聚点

> 所属章节：第十一章 Euclid 空间上的极限和连续  |  文件序号：02  |  难度：基础
> 常见混淆点：1) 聚点不一定属于集合本身——例如 $\{1/n\}$ 的聚点 $0$ 不在集合中，但 $0$ 的任意邻域都包含集合中的点；2) "开"与"闭"不是互斥的概念——$\mathbb{R}^n$ 本身和空集 $\varnothing$ 既开又闭，许多集合（如 $[0,1)$）既不开也不闭

## 1. 学习目标与先修前置

### 学习目标
- 理解 $\mathbb{R}^n$ 中 $\delta$-邻域（开球）的定义与几何意义
- 掌握内点、内部、开集、闭集的形式化定义及其相互联系
- 能用三角不等式严格证明开球是开集、闭球是闭集
- 掌握判断集合开/闭性的两种基本方法：内点验证法和补集法
- 理解聚点和孤立点的定义，掌握求聚点的基本方法
- 理解闭包的定义 $\overline{S}=S\cup S'$ 及其基本性质
- 掌握聚点的序列刻画定理并能进行双向证明

### 先修知识
- 文件 01（第十一章）：$\mathbb{R}^n$ 中的向量、范数 $\|\cdot\|$、距离 $d(\mathbf{x},\mathbf{y})=\|\mathbf{x}-\mathbf{y}\|$、点列收敛的 $\varepsilon$-$N$ 定义、三角不等式
- 文件 01（第一章）：集合的基本概念（子集、并集、交集、补集）
- 读者应熟悉实数比较和不等式的基本操作

---

## 2. 背景与应用场景

在上一节中，我们建立了 $\mathbb{R}^n$ 的基本度量结构——用范数 $\|\cdot\|$ 度量向量长度，用距离 $d(\mathbf{x},\mathbf{y})=\|\mathbf{x}-\mathbf{y}\|$ 度量点与点的远近，并定义了 $\mathbb{R}^n$ 中点列收敛的概念。

有了距离，我们就可以在 $\mathbb{R}^n$ 中讨论"邻近"的概念：一个点附近有哪些点？一个集合是"开放"的还是"封闭"的？这类似于一维情形中实数集上的开区间 $(a,b)$ 和闭区间 $[a,b]$ 的区别——开区间不包含端点，闭区间包含端点。

在 $\mathbb{R}^n$ 中，情况更加丰富：
- 一个集合可能是"到处充满"的（开集），每个点周围都有足够的空间容纳一个小球
- 一个集合可能是"紧实"的（闭集），包含其所有极限位置
- 一个点可能是集合的"凝聚点"（聚点），周围聚集了集合中无限多个其他点

这些拓扑概念是多元函数极限和连续性理论的基础。在下一章学习多元函数的极限时，开集和闭集将扮演核心角色——例如，函数只在定义域的内点处才可能连续，而在聚点处才可能讨论极限。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 示例 |
|------|------|------|
| $O(a,\delta)$ | 以 $a$ 为中心、$\delta$ 为半径的开球（$\delta$-邻域） | $O((0,0),1)=\{(x,y)\mid x^2+y^2<1\}$ |
| $B(a,r)$ | 开球（与 $O(a,r)$ 同义） | $B(a,r)=\{\mathbf{x}\mid \|\mathbf{x}-a\|<r\}$ |
| $\overline{B}(a,r)$ | 闭球 | $\overline{B}(a,r)=\{\mathbf{x}\mid \|\mathbf{x}-a\|\leq r\}$ |
| $\mathring{E}$ 或 $E^\circ$ | 集合 $E$ 的内部（所有内点构成的集合） | $\mathring{\overline{B}}(a,r)=B(a,r)$ |
| $E^c$ | 集合 $E$ 的补集 $\mathbb{R}^n\setminus E$ | |
| $S'$ | 集合 $S$ 的所有聚点构成的集合（导集） | |
| $\overline{S}$ | 集合 $S$ 的闭包 $S\cup S'$ | |
| $\varnothing$ | 空集 | |

### 3.1 $\delta$-邻域（开球）

**定义 11.7（$\delta$-邻域 / 开球）**：设 $a\in\mathbb{R}^n$，$\delta>0$ 是正实数。称集合
$$O(a,\delta)=\{\,\mathbf{x}\in\mathbb{R}^n\mid \|\mathbf{x}-a\|<\delta\,\}$$
为以 $a$ 为中心、$\delta$ 为半径的 **$\delta$-邻域**（$\delta$-neighborhood），也称为 **开球**（open ball），记作 $B(a,\delta)$。

**几何意义**：
- 在 $\mathbb{R}^1$ 中，$O(a,\delta)=(a-\delta,\;a+\delta)$ 是一个开区间。
- 在 $\mathbb{R}^2$ 中，$O(a,\delta)$ 是以 $a$ 为圆心、$\delta$ 为半径的圆的内部（不含边界）。
- 在 $\mathbb{R}^3$ 中，$O(a,\delta)$ 是以 $a$ 为球心、$\delta$ 为半径的球的内部（不含球面）。

**例 1**：在 $\mathbb{R}^2$ 中，$a=(1,2)$，$\delta=3$。
$$O((1,2),3)=\{(x,y)\in\mathbb{R}^2\mid (x-1)^2+(y-2)^2<9\}$$
这是平面上以 $(1,2)$ 为圆心、半径为 $3$ 的圆的内部区域。

**说明**：$\delta$-邻域的概念将一维中"开区间"的"邻近"意推广到了 $n$ 维。$\delta$ 越小，邻域范围越小，表示我们考虑的是点 $a$ 附近越来越小的区域。

---

### 3.2 内点与内部

**定义 11.8（内点）**：设 $E\subseteq\mathbb{R}^n$，$a\in\mathbb{R}^n$。若存在 $\delta>0$，使得
$$O(a,\delta)\subseteq E$$
则称 $a$ 为 $E$ 的 **内点**（interior point）。

换言之，$a$ 是 $E$ 的内点当且仅当存在某个以 $a$ 为中心的开球整个包含在 $E$ 内部。

**定义 11.9（内部）**：集合 $E\subseteq\mathbb{R}^n$ 的**所有内点构成的集合**称为 $E$ 的**内部**（interior），记作
$$\mathring{E}\quad\text{或}\quad E^\circ$$

**例 2（内点的判断）**：在 $\mathbb{R}^1$ 中：
- 取 $E=[0,1]$。则 $0$ 不是 $E$ 的内点，因为任何 $O(0,\delta)=(-\delta,\delta)$ 都包含负数（不在 $E$ 中）。同理 $1$ 也不是内点。
- 任何 $x\in(0,1)$ 都是 $E$ 的内点：取 $\delta=\min\{x,1-x\}>0$，则 $O(x,\delta)\subseteq(0,1)\subseteq[0,1]$。
- 因此 $\mathring{E}=(0,1)$。

**例 3（不同维度的内点）**：在 $\mathbb{R}^2$ 中考虑集合 $E=[0,1]\times\{0\}$（即 $x$ 轴上的线段）。
- 取 $a=(\frac12,0)\in E$。对于任意 $\delta>0$，$O(a,\delta)$ 包含纵坐标 $y\neq0$ 的点，而这些点不在 $E$ 中。因此 $a$ 不是 $E$ 的内点。
- 事实上，$[0,1]\times\{0\}$ 在 $\mathbb{R}^2$ 中的内部是空集：$\mathring{E}=\varnothing$。

这个例子说明：**内点依赖于所在的空间**。在 $\mathbb{R}^1$ 中 $(0,1)$ 全是内点，但在 $\mathbb{R}^2$ 中同一集合被嵌入到更高维空间后内部变成了空集。

---

### 3.3 开集

**定义 11.10（开集）**：设 $E\subseteq\mathbb{R}^n$。若 $E$ 中的每一个点都是 $E$ 的内点，即
$$\forall a\in E,\ \exists\delta>0,\ \text{s.t.}\ O(a,\delta)\subseteq E$$
则称 $E$ 是 $\mathbb{R}^n$ 中的**开集**（open set）。

等价地说，$E$ 是开集当且仅当 $E=\mathring{E}$。

**例 4（开集的例子）**：
- $\mathbb{R}^n$ 本身是开集：对任意 $a\in\mathbb{R}^n$，任何 $\delta>0$ 都有 $O(a,\delta)\subseteq\mathbb{R}^n$。
- 空集 $\varnothing$ 是开集（约定：没有点需要检验，条件空真）。
- 开区间 $(0,1)\subset\mathbb{R}$ 是开集（$\mathbb{R}^1$ 中）。
- 开矩形 $(0,1)\times(0,1)\subset\mathbb{R}^2$ 是开集。

**例 5（非开集的例子）**：
- $[0,1)\subset\mathbb{R}$ 不是开集，因为 $0\in[0,1)$ 不是内点。
- $\{(x,y)\in\mathbb{R}^2\mid x^2+y^2\leq 1\}$ 不是开集，因为边界上的点（如 $(1,0)$）不是内点。

---

### 3.4 闭集

**定义 11.11（闭集）**：设 $E\subseteq\mathbb{R}^n$。若 $E$ 的**补集** $E^c=\mathbb{R}^n\setminus E$ 是开集，则称 $E$ 是 $\mathbb{R}^n$ 中的**闭集**（closed set）。

**例 6（闭集的例子）**：
- $\mathbb{R}^n$ 本身是闭集，因为 $\mathbb{R}^n$ 的补集 $\varnothing$ 是开集。
- 空集 $\varnothing$ 是闭集，因为 $\varnothing$ 的补集 $\mathbb{R}^n$ 是开集。
- 闭区间 $[0,1]\subset\mathbb{R}$ 是闭集（其补集 $(-\infty,0)\cup(1,\infty)$ 是开集）。
- 闭球 $\overline{B}(a,r)=\{\mathbf{x}\in\mathbb{R}^n\mid \|\mathbf{x}-a\|\leq r\}$ 是闭集（将在第 4.2 节证明）。

**例 7（既不开也不闭的集合）**：
- $[0,1)\subset\mathbb{R}$ 既不是开集（$0$ 不是内点）也不是闭集（因为其补集 $(-\infty,0)\cup[1,\infty)$ 不是开集——$1$ 不是内点）。
- 这说明"开"与"闭"不是互补的概念。

---

### 3.5 聚点与孤立点

**定义 11.12（聚点 / 极限点）**：设 $S\subseteq\mathbb{R}^n$，$x\in\mathbb{R}^n$（$x$ 可以属于 $S$，也可以不属于 $S$）。若对任意 $\delta>0$，都有
$$O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$$
则称 $x$ 为 $S$ 的**聚点**（accumulation point / limit point），也称**极限点**。

换言之，$x$ 的任意 $\delta$-邻域中都包含 $S$ 中**异于 $x$** 的点。

$S$ 的所有聚点构成的集合称为 $S$ 的**导集**（derived set），记作 $S'$。

**定义 11.13（孤立点）**：设 $S\subseteq\mathbb{R}^n$，$x\in S$。若存在 $\delta>0$，使得
$$O(x,\delta)\cap S=\{x\}$$
即 $x$ 的某个邻域中除 $x$ 本身外不含 $S$ 的其他点，则称 $x$ 为 $S$ 的**孤立点**（isolated point）。

**聚点与孤立点的关系**：集合 $S$ 中的点要么是聚点，要么是孤立点（二者必居其一且仅居其一）。不是 $S$ 中的点也可以成为 $S$ 的聚点。

**例 8（聚点与孤立点）**：设 $S=\{\frac1n\mid n\in\mathbb{N}^+\}\cup\{0\}\subset\mathbb{R}$。
- 对任意 $\delta>0$，取正整数 $N>\frac1\delta$，则 $\frac1N\in O(0,\delta)$ 且 $\frac1N\neq0$，故 $0$ 是 $S$ 的聚点。
- 对任意固定的 $n_0\in\mathbb{N}^+$，取 $\delta=\frac1{n_0(n_0+1)}$（即 $\frac1{n_0}$ 到相邻点 $\frac1{n_0+1}$ 距离的一半），则 $O(\frac1{n_0},\delta)\cap S=\{\frac1{n_0}\}$，故 $\frac1{n_0}$ 是孤立点。
- 因此 $S'=\{0\}$，$S\setminus S'=\{\frac1n\mid n\in\mathbb{N}^+\}$ 中的点都是孤立点。

**例 9（不同集合的聚点）**：
- $S=(0,1)\subset\mathbb{R}$：每个 $x\in[0,1]$ 都是 $S$ 的聚点。$S'=[0,1]$。
- $S=\mathbb{Q}\subset\mathbb{R}$：每个 $x\in\mathbb{R}$ 都是 $\mathbb{Q}$ 的聚点，$S'=\mathbb{R}$。
- $S=\mathbb{Z}\subset\mathbb{R}$：$\mathbb{Z}$ 没有聚点（每个整数都是孤立点），$S'=\varnothing$。
- $S=\{\frac1n+\frac1m\mid n,m\in\mathbb{N}^+\}\subset\mathbb{R}$：聚点包括 $\{\frac1n\mid n\in\mathbb{N}^+\}$ 和 $\{0\}$。

---

### 3.6 闭包

**定义 11.14（闭包）**：设 $S\subseteq\mathbb{R}^n$。$S$ 的**闭包**（closure）定义为 $S$ 与其所有聚点之并：
$$\overline{S}=S\cup S'$$
其中 $S'$ 是 $S$ 的导集。

**例 10**：
- $S=(0,1)\subset\mathbb{R}$，$S'=[0,1]$，$\overline{S}=S\cup S'=[0,1]$。
- $S=\{\frac1n\mid n\in\mathbb{N}^+\}\subset\mathbb{R}$，$S'=\{0\}$，$\overline{S}=S\cup\{0\}$。
- $S=\mathbb{Q}\subset\mathbb{R}$，$S'=\mathbb{R}$，$\overline{S}=\mathbb{R}$。
- $S=\mathbb{Z}\subset\mathbb{R}$，$S'=\varnothing$，$\overline{S}=\mathbb{Z}$。

**说明**：闭包 $\overline{S}$ 可以直观地理解为"在 $S$ 中加入所有极限位置后得到的集合"。闭包总是闭集（将在第 4.5 节证明）。

---

## 4. 原理与方法

### 4.1 开球是开集的证明（三角不等式法）

**命题**：$\mathbb{R}^n$ 中的任意开球 $B(a,r)$ 是开集。

**分析**：要证明 $B(a,r)$ 是开集，按照定义 11.10，需要证明 $B(a,r)$ 中的每个点都是它的内点。即对任意 $x\in B(a,r)$，需要找到一个 $\delta>0$，使得整个开球 $O(x,\delta)$ 都包含在 $B(a,r)$ 中。几何直觉告诉我们：从 $x$ 到 $a$ 的距离是 $\|x-a\|$，只要 $\delta$ 不超过剩下的"余量" $r-\|x-a\|$，整个 $O(x,\delta)$ 就不会跑出 $B(a,r)$。

**证明**：

任取 $x\in B(a,r)$。由开球的定义知
$$d=\|x-a\|<r$$

令 $\delta=r-d>0$（见图示：$\delta$ 是从 $x$ 到球面 $B(a,r)$ 边界的"剩余距离"）。

下面证明 $O(x,\delta)\subseteq B(a,r)$。任取 $y\in O(x,\delta)$，有 $\|y-x\|<\delta$。

由三角不等式（定理 11.2）：
$$
\begin{aligned}
\|y-a\| &= \|(y-x)+(x-a)\| \\
&\leq \|y-x\| + \|x-a\| \\
&< \delta + d \\
&= (r-d) + d = r
\end{aligned}
$$

因此 $\|y-a\|<r$，即 $y\in B(a,r)$。

故 $O(x,\delta)\subseteq B(a,r)$，从而 $x$ 是 $B(a,r)$ 的内点。由 $x$ 的任意性，$B(a,r)$ 的每个点都是内点，所以 $B(a,r)$ 是开集。证毕。

**核心思想**：三角不等式 $\|y-a\|\leq\|y-x\|+\|x-a\|$ 是将 $y$ 到 $a$ 的距离分解为从 $y$ 到 $x$ 再到 $a$ 的两段路程之和。通过控制第一段路程 $\|y-x\|<\delta$，确保了总路程不超过 $r$。

---

### 4.2 闭球是闭集的证明

**命题**：$\mathbb{R}^n$ 中的任意闭球 $\overline{B}(a,r)=\{\mathbf{x}\in\mathbb{R}^n\mid \|\mathbf{x}-a\|\leq r\}$ 是闭集。

**分析**：按照闭集的定义（定义 11.11），要证明 $\overline{B}(a,r)$ 是闭集，只需证明其补集
$$(\overline{B}(a,r))^c=\{\mathbf{x}\in\mathbb{R}^n\mid \|\mathbf{x}-a\|>r\}$$
是开集。关键在于：对补集中的任意点 $x$（满足 $\|x-a\|>r$），需要找到一个 $\delta>0$ 使得整个 $O(x,\delta)$ 仍然在补集中，即其中的每个点 $y$ 都满足 $\|y-a\|>r$。

**证明**：

记 $F=\overline{B}(a,r)$，其补集 $F^c=\{\mathbf{x}\in\mathbb{R}^n\mid \|\mathbf{x}-a\|>r\}$。我们要证明 $F^c$ 是开集。

任取 $x\in F^c$，则 $\|x-a\|>r$。令
$$d=\|x-a\|,\quad \delta=d-r>0$$

下面证明 $O(x,\delta)\subseteq F^c$。任取 $y\in O(x,\delta)$，有 $\|y-x\|<\delta$。

由三角不等式（反向形式）：对任意向量 $\mathbf{u},\mathbf{v}$，有 $\|\mathbf{u}\|-\|\mathbf{v}\|\leq\|\mathbf{u}-\mathbf{v}\|$。令 $\mathbf{u}=y-a$，$\mathbf{v}=x-a$，则：
$$\|y-a\|\geq\|x-a\|-\|y-x\|$$

这是因为由三角不等式 $\|x-a\|\leq\|x-y\|+\|y-a\|$，移项即得 $\|y-a\|\geq\|x-a\|-\|x-y\|$。

因此：
$$\|y-a\|\geq\|x-a\|-\|y-x\| > d-\delta = d-(d-r)=r$$

故 $\|y-a\|>r$，即 $y\in F^c$。

所以 $O(x,\delta)\subseteq F^c$，$x$ 是 $F^c$ 的内点。由 $x$ 的任意性，$F^c$ 是开集，从而 $F$ 是闭集。证毕。

**两种证明的比较**：

| 证明 | 目标 | 选取的 $\delta$ | 依赖的不等式 |
|------|------|----------------|-------------|
| 开球是开集 | $B(a,r)$ 的每个点都是内点 | $\delta=r-\|x-a\|$ | 正三角不等式 $\|y-a\|\leq\|y-x\|+\|x-a\|$ |
| 闭球是闭集 | 补集的每个点都是内点 | $\delta=\|x-a\|-r$ | 逆三角不等式 $\|y-a\|\geq\|x-a\|-\|y-x\|$ |

仔细观察发现，两种证明中 $\delta$ 的选取正好是"到边界的距离"——开球情况是从内点到边界的距离，补集情况是从外点到边界的距离。这是判断开/闭集的一种直观方法。

---

### 4.3 开/闭集分类判断方法

判断一个集合 $E\subseteq\mathbb{R}^n$ 是开集、闭集、或两者都不是，有以下两种基本方法。

**方法一：内点验证法（直接法）**

对开集：
1. 任取 $x\in E$。
2. 尝试构造 $\delta>0$ 使得 $O(x,\delta)\subseteq E$。
3. 若对所有 $x\in E$ 都能成功构造，则 $E$ 是开集；若存在某个 $x\in E$ 使构造失败，则 $E$ 不是开集。

对闭集：
1. 考虑补集 $E^c$。
2. 用上述方法判断 $E^c$ 是否为开集。
3. 若 $E^c$ 是开集，则 $E$ 是闭集；否则不是。

**方法二：边界判断法（直观法）**

- 若 $E$ 不包含任何边界点（或等价地，$E$ 中的每个点都可以向外扩展一个小球），则 $E$ 是开集。
- 若 $E$ 包含所有边界点（即边界点全部属于 $E$），则 $E$ 是闭集。
- 若 $E$ 包含部分边界点但不包含全部，则 $E$ 既不开也不闭。

**例 11（分类判断）**：判断下列 $\mathbb{R}^2$ 中集合的开/闭性。

(1) $A=\{(x,y)\mid x>0,\ y>0\}$
- **开集**。任取 $(x_0,y_0)\in A$，有 $x_0>0$，$y_0>0$。取 $\delta=\min\{x_0,y_0\}>0$，则 $O((x_0,y_0),\delta)$ 中任意点的横坐标 $>x_0-\delta\geq0$，纵坐标 $>y_0-\delta\geq0$，故整个开球包含在 $A$ 中。

(2) $B=\{(x,y)\mid x\geq 0,\ y\geq 0\}$
- **闭集**。补集 $B^c=\{(x,y)\mid x<0\ \text{或}\ y<0\}$ 是开集（任取 $(x_0,y_0)\in B^c$，若 $x_0<0$，取 $\delta=|x_0|/2$ 即可）。故 $B$ 是闭集。

(3) $C=\{(x,y)\mid x^2+y^2<1\}\cup\{(1,0)\}$
- **既不开也不闭**。$(1,0)$ 属于 $C$ 但不是内点（它的任何邻域都包含 $x^2+y^2>1$ 的点），故 $C$ 不是开集。同时，$(0,1)$（不在 $C$ 中）是 $C$ 的边界点，但 $(0,1)$ 的任何邻域都包含 $C$ 中的点，故 $(0,1)\in (C^c)'$ 但 $(0,1)\notin C^c$，所以 $C^c$ 不是开集，$C$ 不是闭集。

**与一维情形的类比**：
- 一维开集：开区间 $(a,b)$ 及其任意并
- 一维闭集：闭区间 $[a,b]$ 及其有限并
- 但 $\mathbb{R}^n$ 中的开/闭集远比一维复杂（例如 $\mathbb{R}^2$ 中有"开圆盘"、"开矩形"等各种形状的开集）

---

### 4.4 聚点的序列刻画定理

**定理 11.4（聚点的序列刻画）**：设 $S\subseteq\mathbb{R}^n$，$x\in\mathbb{R}^n$。则 $x$ 是 $S$ 的聚点当且仅当存在点列 $\{x_k\}_{k=1}^\infty\subseteq S\setminus\{x\}$（即 $\{x_k\}$ 中的每个点都属于 $S$ 且不等于 $x$），使得
$$x_k\to x\quad(k\to\infty)$$

**证明**：双向证明。

**$(\Rightarrow)$ 方向（聚点 $\Rightarrow$ 存在序列）**：

假设 $x$ 是 $S$ 的聚点。由定义 11.12，对任意 $\delta>0$，有
$$O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$$

我们利用这一性质构造序列。对每个正整数 $k\in\mathbb{N}^+$，取 $\delta_k=\frac1k>0$。由聚点定义，存在
$$x_k\in O\Big(x,\frac1k\Big)\cap(S\setminus\{x\})$$

这样得到的 $\{x_k\}_{k=1}^\infty$ 满足：
- 每个 $x_k\in S\setminus\{x\}$（由构造方式）
- $\|x_k-x\|<\frac1k$（因为 $x_k\in O(x,\frac1k)$）

对任意 $\varepsilon>0$，取正整数 $N>\frac1\varepsilon$，则当 $k>N$ 时：
$$\|x_k-x\|<\frac1k\leq\frac1N<\varepsilon$$

由点列收敛的 $\varepsilon$-$N$ 定义（定义 11.6），$x_k\to x$。故满足要求的序列存在。

**$(\Leftarrow)$ 方向（存在序列 $\Rightarrow$ 聚点）**：

假设存在序列 $\{x_k\}\subseteq S\setminus\{x\}$ 满足 $x_k\to x$。

对任意 $\delta>0$，由 $x_k\to x$ 的定义，存在 $N\in\mathbb{N}^+$，使得当 $k>N$ 时恒有
$$\|x_k-x\|<\delta$$

特别地，取 $k=N+1$，则 $x_{N+1}\in S\setminus\{x\}$ 且 $\|x_{N+1}-x\|<\delta$，即
$$x_{N+1}\in O(x,\delta)\cap(S\setminus\{x\})$$

因此 $O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$。由 $\delta>0$ 的任意性，$x$ 是 $S$ 的聚点。证毕。

**定理的意义**：
- 该定理建立了聚点的**几何定义**（邻域总含有异于点）与**分析定义**（存在序列逼近）之间的等价关系。
- 在实际应用中，当我们想证明某点是聚点时，构造一个收敛到该点的序列往往比直接验证邻域条件更方便（见例题 3）。
- 反过来，当已知某点是聚点时，我们可以在该点附近"取样"得到收敛序列，这为后续研究极限和连续性提供了有力工具。

---

### 4.5 闭包的性质

**定理 11.5（闭包的基本性质）**：设 $S,T\subseteq\mathbb{R}^n$，则以下性质成立：

1. **包含性**：$S\subseteq\overline{S}$。
2. **单调性**：若 $S\subseteq T$，则 $\overline{S}\subseteq\overline{T}$。
3. **闭包是闭集**：$\overline{S}$ 是 $\mathbb{R}^n$ 中的闭集。
4. **幂等性**：$(\overline{S})=\overline{S}$，即取闭包运算重复两次等价于取一次。
5. **闭集判别**：$S$ 是闭集当且仅当 $S=\overline{S}$。

**证明**：

**(1) 包含性**：由闭包定义 $\overline{S}=S\cup S'$，显然 $S\subseteq\overline{S}$。

**(2) 单调性**：设 $S\subseteq T$。任取 $x\in\overline{S}=S\cup S'$。
- 若 $x\in S\subseteq T\subseteq\overline{T}$，则 $x\in\overline{T}$。
- 若 $x\in S'$（$x$ 是 $S$ 的聚点），由聚点定义，$\forall\delta>0$，$O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$。由于 $S\subseteq T$，$O(x,\delta)\cap(T\setminus\{x\})\supseteq O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$，故 $x$ 也是 $T$ 的聚点，即 $x\in T'\subseteq\overline{T}$。
因此 $\overline{S}\subseteq\overline{T}$。

**(3) 闭包是闭集**：需证 $(\overline{S})^c=\mathbb{R}^n\setminus\overline{S}$ 是开集。

任取 $x\notin\overline{S}$。这意味着 $x\notin S$ 且 $x$ 不是 $S$ 的聚点。由聚点定义的否定：存在 $\delta>0$，使得
$$O(x,\delta)\cap(S\setminus\{x\})=\varnothing$$
由于 $x\notin S$，$S\setminus\{x\}=S$，故 $O(x,\delta)\cap S=\varnothing$。

现在证明 $O(x,\delta)\subseteq(\overline{S})^c$。任取 $y\in O(x,\delta)$。由 $O(x,\delta)\cap S=\varnothing$ 知 $y\notin S$。进一步，由于 $O(x,\delta)$ 是开集，存在 $\varepsilon>0$ 使得 $O(y,\varepsilon)\subseteq O(x,\delta)$，从而 $O(y,\varepsilon)\cap S=\varnothing$，故 $y$ 不是 $S$ 的聚点。因此 $y\in(\overline{S})^c$。

故 $O(x,\delta)\subseteq(\overline{S})^c$，即 $x$ 是 $(\overline{S})^c$ 的内点。由 $x$ 的任意性，$(\overline{S})^c$ 是开集，$\overline{S}$ 是闭集。证毕。

**(4) 幂等性**：由 (3) 知 $\overline{S}$ 是闭集，由 (5) 的等价关系（下面将证明），闭集等于自身的闭包，故 $\overline{\overline{S}}=\overline{S}$。

**(5) 闭集判别**：分两个方向。
- $(\Rightarrow)$：若 $S$ 是闭集，证 $S=\overline{S}$。由 (1) 已有 $S\subseteq\overline{S}$，只需证 $\overline{S}\subseteq S$。任取 $x\in\overline{S}=S\cup S'$。若 $x\in S$ 则显然。若 $x\in S'$（聚点），反设 $x\notin S$，则 $x\in S^c$。由于 $S$ 是闭集，$S^c$ 是开集，故存在 $\delta>0$ 使 $O(x,\delta)\subseteq S^c$，从而 $O(x,\delta)\cap S=\varnothing$，与 $x$ 是 $S$ 的聚点矛盾。故 $x\in S$。所以 $\overline{S}\subseteq S$，从而 $\overline{S}=S$。
- $(\Leftarrow)$：若 $S=\overline{S}$，由 (3) 知 $\overline{S}$ 是闭集，故 $S$ 是闭集。

证毕。

---

## 5. 例题

### 例题 1：判断集合的开/闭性

判断下列 $\mathbb{R}^2$ 中集合的开/闭性，并说明理由。

(1) $A=\{(x,y)\mid x^2+y^2<4\}$
(2) $B=\{(x,y)\mid 1\leq x^2+y^2\leq 4\}$
(3) $C=\{(x,y)\mid x>0,\ y\geq 0\}$

**解**：

**(1)** $A$ 是以 $(0,0)$ 为圆心、$2$ 为半径的开球，由第 4.1 节的命题知 $A$ 是**开集**。

**(2)** $B$ 是闭圆环 $\{\mathbf{x}\mid 1\leq\|\mathbf{x}\|\leq 2\}$。我们证明 $B$ 是**闭集**。

采用补集法。$B^c=\{\mathbf{x}\mid \|\mathbf{x}\|<1\}\cup\{\mathbf{x}\mid \|\mathbf{x}\|>2\}$。前一部分是开球 $B(0,1)$，后一部分是闭球 $\overline{B}(0,2)$ 的补集。开球是开集（第 4.1 节），闭球的补集是开集（第 4.2 节），两个开集之并仍是开集（开集的任意并是开集，此处不证明，留给后续课程）。故 $B^c$ 是开集，$B$ 是闭集。

**(3)** $C=\{(x,y)\mid x>0,\ y\geq 0\}$ 是右半平面与上半平面之交（包含正 $x$ 轴）。

- $C$ **不是开集**：取 $(1,0)\in C$。对任意 $\delta>0$，$O((1,0),\delta)$ 中包含纵坐标为负的点（如 $(1,-\delta/2)$），这些点不在 $C$ 中（因为 $y<0$），故 $(1,0)$ 不是内点。
- $C$ **不是闭集**：考虑补集 $C^c=\{(x,y)\mid x\leq 0\ \text{或}\ y<0\}$。取 $(0,0)\in C^c$（因为 $0\leq0$）。对任意 $\delta>0$，$O((0,0),\delta)$ 中包含 $x>0,\ y>0$ 的点，这些点属于 $C$ 而不属于 $C^c$。故 $(0,0)$ 不是 $C^c$ 的内点，$C^c$ 不是开集，$C$ 不是闭集。

因此 $C$ **既不开也不闭**。

---

### 例题 2：求聚点和闭包

求下列 $\mathbb{R}$ 中集合的聚点（导集）和闭包。

(1) $S_1=\left\{\frac{n-1}{n}\;\Big|\; n\in\mathbb{N}^+\right\}$
(2) $S_2=\mathbb{Q}\cap[0,1]$
(3) $S_3=\left\{\frac{m}{2^n}\;\Big|\; m\in\mathbb{Z},\ n\in\mathbb{N}\right\}$（二进制有理数集）

**解**：

**(1)** $S_1=\{0,\frac12,\frac23,\frac34,\frac45,\dots\}$。

当 $n\to\infty$ 时，$\frac{n-1}{n}=1-\frac1n\to1$。对任意 $\delta>0$，存在 $N>\frac1\delta$，使得 $\frac{N-1}{N}\in O(1,\delta)\cap(S_1\setminus\{1\})$，故 $1$ 是 $S_1$ 的聚点。除此之外，$S_1$ 中的每个点都是孤立点（以到相邻点的半距离为邻域半径即可）。

还需要验证：$S_1$ 是否有其他聚点？$S_1$ 的点都在 $[0,1)$ 中，且除了趋近于 $1$ 外没有任何子列趋近于其他值。故 $S_1'=\{1\}$。

闭包：$\overline{S_1}=S_1\cup S_1'=S_1\cup\{1\}$。

**(2)** $S_2=[0,1]$ 中的全体有理数。

我们知道 $\mathbb{Q}$ 在 $\mathbb{R}$ 中稠密：对任意 $x\in\mathbb{R}$，有理数都可以任意逼近 $x$。特别地，对任意 $x\in[0,1]$，存在有理数列 $\{r_k\}\subseteq\mathbb{Q}\cap[0,1]$ 收敛到 $x$。因此每个 $x\in[0,1]$ 都是 $S_2$ 的聚点。而 $x\notin[0,1]$ 时，取 $\delta=\min\{|x-0|,|x-1|\}/2>0$，则 $O(x,\delta)\cap S_2=\varnothing$，故 $x$ 不是 $S_2$ 的聚点。

因此 $S_2'=[0,1]$，$\overline{S_2}=S_2\cup[0,1]=[0,1]$。

**(3)** $S_3$ 是所有分母为 $2$ 的幂的分数。

- $0$ 是 $S_3$ 的聚点（取 $x_k=1/2^k\to0$）。
- $1$ 是 $S_3$ 的聚点（取 $x_k=1-1/2^k\to1$）。
- 事实上，由于二进制有理数在 $\mathbb{R}$ 中稠密，任意 $x\in\mathbb{R}$ 都可以被形如 $m/2^n$ 的数任意逼近。故 $S_3'=\mathbb{R}$。
- $\overline{S_3}=S_3\cup\mathbb{R}=\mathbb{R}$。

**小结**：聚点的求解实质上是在问"有哪些极限位置可以被 $S$ 中的点列逼近"。常用的方法是：
1. 观察 $S$ 中是否有趋于某个极限的子列。
2. 利用已知的稠密性结论（如有理数在实数中稠密）。
3. 检查集合的边界点是否可能成为聚点。

---

### 例题 3：用序列刻画定理证明聚点

设 $S=\{\frac1n+\frac1m\mid n,m\in\mathbb{N}^+\}\subset\mathbb{R}$。

(1) 证明 $0$ 是 $S$ 的聚点。
(2) 证明 $2$ 不是 $S$ 的聚点。

**解**：

**(1) 证明 $0$ 是聚点**：

构造序列 $x_k=\frac1k+\frac1k=\frac2k\in S$（取 $n=m=k$），则 $x_k\in S\setminus\{0\}$ 且
$$\lim_{k\to\infty}x_k=\lim_{k\to\infty}\frac2k=0$$
由定理 11.4（聚点的序列刻画），$0$ 是 $S$ 的聚点。

**验证（直接用聚点定义）**：对任意 $\delta>0$，取 $N>\frac2\delta$，则 $\frac2N\in O(0,\delta)\cap(S\setminus\{0\})$，故 $0$ 是聚点。两种方法等价。

**(2) 证明 $2$ 不是聚点**：

反证法。假设 $2$ 是 $S$ 的聚点。由定理 11.4，存在 $\{x_k\}\subseteq S\setminus\{2\}$ 使得 $x_k\to2$。

但 $S$ 中的每个元素都是 $\frac1n+\frac1m$ 的形式，其中 $n,m\in\mathbb{N}^+$。由于 $n,m\geq1$，有
$$\frac1n+\frac1m\leq 1+1=2$$
且等号 $\frac1n+\frac1m=2$ 成立当且仅当 $n=m=1$，此时 $\frac11+\frac11=2$，对应的元素就是 $2$ 本身。

因此 $S\setminus\{2\}\subseteq(0,2)$，即 $S\setminus\{2\}$ 中的每个元素都严格小于 $2$。若 $x_k\to2$ 且每个 $x_k<2$，则 $2-x_k>0$ 趋于 $0$，这是可能的（例如取 $x_k=2-\frac1k$）。但是 $S$ 中是否有趋近于 $2$ 的序列呢？

对任意 $n,m\in\mathbb{N}^+$，$\frac1n+\frac1m\leq\frac11+\frac1m=1+\frac1m<2$，且当 $m\to\infty$ 时 $1+\frac1m\to1$。最接近 $2$ 的值是 $1+1=2$（$n=m=1$）和 $1+\frac12=1.5$（$n=1,m=2$ 或 $n=2,m=1$）。两者之间没有可以逼近 $2$ 的其他值。

更严格地：任取 $x\in S\setminus\{2\}$，有 $x\leq\frac11+\frac12=1.5$（因为至少有一个 $n$ 或 $m$ 大于 $1$，不妨设 $n>1$，则 $\frac1n+\frac1m\leq\frac12+1=1.5$）。因此 $S\setminus\{2\}\subseteq[0,1.5]$。

若 $\{x_k\}\subseteq S\setminus\{2\}$ 收敛，则其极限 $\leq1.5$，不可能等于 $2$。矛盾。故 $2$ 不是 $S$ 的聚点。

**直接方法**：取 $\delta=0.5$，则 $O(2,0.5)=(1.5,2.5)$ 与 $S\setminus\{2\}$ 的交集为空（因为 $S\setminus\{2\}\subseteq[0,1.5]$），故 $2$ 不是聚点。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 认为"开集"和"闭集"是互斥的概念，一个集合要么开要么闭 | $\mathbb{R}^n$ 本身和空集 $\varnothing$ 是既开又闭的；而 $[0,1)$ 既不开也不闭。"开"和"闭"不是互补关系 |
| 认为聚点必须属于集合本身 | 聚点可以属于集合（如 $[0,1]$ 中的每个点都是聚点），也可以不属于集合（如 $\{1/n\}$ 的聚点 $0$ 不在集合中）。聚点的定义只要求邻域包含集合中异于该点的点，不要求该点本身属于集合 |
| 将一个集合的"内部"与集合本身等同，认为开区间 $(0,1)$ 的"闭包"还是 $(0,1)$ | 开区间 $(0,1)$ 的闭包是 $[0,1]$，因为 $0$ 和 $1$ 是它的聚点。闭包会"填满"所有极限位置 |
| 用一维直觉处理高维问题：认为 $\{(x,0)\mid 0<x<1\}\subset\mathbb{R}^2$ 是开集（因为它在 $\mathbb{R}^1$ 中是开区间） | 在 $\mathbb{R}^2$ 中，该集合不包含任何内点（任何二维开球都会包含 $y\neq0$ 的点），因此内部为空，不是开集 |
| 混淆"聚点"与"边界点" | 边界点一定是聚点（不一定），边界点强调"任意邻域既含 $E$ 中点也含 $E$ 外点"，聚点只要求"任意邻域含异于该点的 $E$ 中点"。例如孤立点可能是边界点但不是聚点 |

### 检查点

- [ ] 能否写出 $\mathbb{R}^n$ 中 $\delta$-邻域 $O(a,\delta)$ 的定义？
- [ ] 能否准确区分内点、聚点、孤立点、边界点的概念？
- [ ] 能否用三角不等式完整证明"开球是开集"？
- [ ] 能否用补集法证明"闭球是闭集"？
- [ ] 能否列举一个既不开也不闭的集合并解释原因？
- [ ] 能否陈述并证明聚点的序列刻画定理（两个方向）？
- [ ] 求出一个给定集合的聚点和闭包（如 $\{1/n\}$ 或 $\mathbb{Q}$）？
- [ ] 是否理解闭包的五条基本性质并能证明其中至少三条？
- [ ] 能否用序列刻画定理将几何问题转化为分析问题？

---

## 练习题

### 基础巩固

**1.** 判断下列 $\mathbb{R}$ 中集合的开/闭性，并说明理由。

(1) $A=(0,1]$
(2) $B=\mathbb{N}$
(3) $C=\left\{\frac1n\;\Big|\; n\in\mathbb{N}^+\right\}$
(4) $D=\varnothing$

<details><summary>参考答案</summary>

(1) $A=(0,1]$ **既不开也不闭**。
- 不是开集：$1\in A$，但任何 $O(1,\delta)=(1-\delta,1+\delta)$ 都包含大于 $1$ 的数，这些数不在 $A$ 中，故 $1$ 不是内点。
- 不是闭集：补集 $A^c=(-\infty,0]\cup(1,\infty)$。$0\in A^c$ 但 $0$ 的任意邻域都包含 $(0,\delta)$ 中的正数（属于 $A$），故 $0$ 不是 $A^c$ 的内点，$A^c$ 不是开集。

(2) $B=\mathbb{N}$ **是闭集**。
补集 $\mathbb{N}^c=\mathbb{R}\setminus\mathbb{N}=\bigcup_{n\in\mathbb{Z}}(n,n+1)$ 是开区间的并，每个开区间是开集，开集的任意并是开集。故 $\mathbb{N}^c$ 是开集，$\mathbb{N}$ 是闭集。
$B$ 不是开集：$1\in B$，但 $O(1,\delta)$ 包含非整数。

(3) $C=\{\frac1n\mid n\in\mathbb{N}^+\}$ **不是开集也不是闭集**。
- 不是开集：每个点都是孤立点，不是内点。
- 不是闭集：$0$ 是 $C$ 的聚点但 $0\notin C$，故 $C$ 不含所有聚点，不是闭集。

(4) $D=\varnothing$ **既开又闭**。
- 开集：没有点需要检验，条件空真。
- 闭集：补集 $\varnothing^c=\mathbb{R}^n$ 是开集。

</details>

---

**2.** 求下列 $\mathbb{R}$ 中集合的聚点（导集）和闭包。

(1) $S_1=\left\{\frac{(-1)^n}{n}\;\Big|\; n\in\mathbb{N}^+\right\}$
(2) $S_2=\left\{\frac{n}{n+1}\;\Big|\; n\in\mathbb{N}^+\right\}\cup\{1\}$
(3) $S_3=\mathbb{Q}$

<details><summary>参考答案</summary>

(1) $S_1=\{-1, \frac12, -\frac13, \frac14, -\frac15, \dots\}$。

当 $n\to\infty$ 时，$\frac{(-1)^n}{n}\to0$。对任意 $\delta>0$，存在 $N>\frac1\delta$，令 $n=N$，则 $\frac{(-1)^n}{n}\in O(0,\delta)\cap(S_1\setminus\{0\})$，故 $0$ 是聚点。
- 除了 $0$ 之外，$S_1$ 中的每个点都是孤立点（以到相邻点的半距离为半径即可）。
- 有没有其他聚点？没有，因为 $S_1$ 的绝对值趋于 $0$。
- 因此 $S_1'=\{0\}$，$\overline{S_1}=S_1\cup\{0\}$。

(2) $S_2=\{\frac12,\frac23,\frac34,\frac45,\dots\}\cup\{1\}$。

$\frac{n}{n+1}=1-\frac1{n+1}\to1$，故 $1$ 是 $S_2$ 的聚点。注意 $1$ 已经在 $S_2$ 中，但这不妨碍它是聚点。
- 其他点都是孤立点，没有其他聚点。
- $S_2'=\{1\}$，$\overline{S_2}=S_2\cup\{1\}=S_2$（因为 $1$ 已经在 $S_2$ 中）。

注意到 $\overline{S_2}=S_2$，由定理 11.5(5) 知 $S_2$ 是闭集。事实上，补集是开集，可以直接验证。

(3) $S_3=\mathbb{Q}$。

有理数在实数中稠密：对任意 $x\in\mathbb{R}$，$\mathbb{Q}$ 中存在序列收敛到 $x$。故每个实数都是 $\mathbb{Q}$ 的聚点。
- $S_3'=\mathbb{R}$，$\overline{S_3}=\mathbb{Q}\cup\mathbb{R}=\mathbb{R}$。

</details>

---

**3.** 用序列刻画定理证明：设 $x$ 是 $S$ 的聚点，则 $x$ 也是 $\overline{S}$ 的聚点。反之是否成立？

<details><summary>参考答案</summary>

**正向**：设 $x$ 是 $S$ 的聚点。由定理 11.4，存在 $\{x_k\}\subseteq S\setminus\{x\}$ 使得 $x_k\to x$。由于 $S\subseteq\overline{S}$，$\{x_k\}\subseteq\overline{S}\setminus\{x\}$，再次由定理 11.4 知 $x$ 是 $\overline{S}$ 的聚点。

**反向**：不成立。反例：取 $S=\{0\}\subset\mathbb{R}$。$S$ 没有聚点（$S'=\varnothing$）。但 $\overline{S}=S$（因为 $S'=\varnothing$），故 $\overline{S}=S$ 也没有聚点。这里"不成立"是说：$x$ 是 $\overline{S}$ 的聚点不一定推出 $x$ 是 $S$ 的聚点。

更合适的反例：$S=(0,1)$。$\overline{S}=[0,1]$。$0$ 是 $\overline{S}$ 的聚点，也是 $S$ 的聚点。那有没有 $\overline{S}$ 的聚点不是 $S$ 的聚点的情况？

考虑 $S=\{\frac1n\mid n\in\mathbb{N}^+\}$。$S'=\{0\}$，$\overline{S}=S\cup\{0\}$。$\overline{S}$ 的聚点仍是 $\{0\}$，也是 $S$ 的聚点。

再看 $S=(0,1)\cup\{2\}$。$\overline{S}=[0,1]\cup\{2\}$。$2$ 是 $\overline{S}$ 的孤立点（不是聚点）。$\overline{S}$ 的聚点是 $[0,1]$，也就是 $S$ 的聚点。

看来需要更好的反例。事实上，可以证明：$\overline{S}$ 的聚点一定是 $S$ 的聚点。因为若 $x$ 是 $\overline{S}=S\cup S'$ 的聚点，则存在 $\{x_k\}\subseteq (S\cup S')\setminus\{x\}$ 且 $x_k\to x$。去掉其中属于 $S'$ 但不在 $S$ 中的点，仍然可以用三角不等式逼近。这个结论的严格证明较复杂，此处不展开。

所以这个练习题的答案是：正向成立，反向也成立（需要证明）。更精确地说，$\overline{S}$ 的聚点集等于 $S$ 的聚点集。

</details>

---

### 迁移应用

**4.** 设 $S\subseteq\mathbb{R}^n$。证明：$S$ 是闭集当且仅当 $S$ 包含其所有聚点。

<details><summary>参考答案</summary>

**$(\Rightarrow)$**：设 $S$ 是闭集。任取 $x\in S'$（$x$ 是 $S$ 的聚点），需证 $x\in S$。反设 $x\notin S$，则 $x\in S^c$。由于 $S$ 是闭集，$S^c$ 是开集，故存在 $\delta>0$ 使得 $O(x,\delta)\subseteq S^c$。于是 $O(x,\delta)\cap S=\varnothing$，这与 $x$ 是 $S$ 的聚点矛盾（聚点要求对任意 $\delta>0$，$O(x,\delta)\cap(S\setminus\{x\})\neq\varnothing$）。因此 $x\in S$，故 $S'\subseteq S$。

**$(\Leftarrow)$**：设 $S$ 包含其所有聚点，即 $S'\subseteq S$。则 $\overline{S}=S\cup S'=S$。由定理 11.5(3)，$\overline{S}$ 是闭集，故 $S$ 是闭集。

**等价结论**：闭集 $\iff$ 包含所有聚点 $\iff$ $\overline{S}=S$。

</details>

---

**5.** 设 $S=\{(x,y)\in\mathbb{R}^2\mid 0<x^2+y^2<1\}$。求 $S$ 的聚点、内部、闭包，并判断 $S$ 的开/闭性。

<details><summary>参考答案</summary>

**分析**：$S$ 是 $\mathbb{R}^2$ 中去掉原点后的单位开圆盘（不含边界）。

**聚点**：
- 对任意满足 $x^2+y^2<1$ 的点 $(x,y)$，它显然在 $S$ 中或是 $0$ 点。若 $x^2+y^2<1$ 且 $(x,y)\neq(0,0)$，则 $(x,y)\in S$，且存在 $S$ 中的序列收敛到它（自身常数列即可），故它是聚点。
- 原点 $(0,0)$：虽然 $(0,0)\notin S$，但取 $P_k=(\frac1k,0)\in S$，则 $P_k\to(0,0)$，故 $(0,0)$ 是 $S$ 的聚点。
- 边界上的点：设 $x^2+y^2=1$，取 $P_k=(1-\frac1k)\cdot(x,y)\in S$（从内部逼近），则 $P_k\to(x,y)$，故边界点也是聚点。
- 外部点：若 $x^2+y^2>1$，取 $\delta=(\sqrt{x^2+y^2}-1)/2>0$，则 $O((x,y),\delta)\cap S=\varnothing$，不是聚点。

因此 $S'=\{(x,y)\mid x^2+y^2\leq 1\}$（包含边界和原点的闭单位圆盘）。

**内部**：
$S$ 中的每个点 $(x,y)\neq(0,0)$ 且 $x^2+y^2<1$，取 $\delta=\min\{\|(x,y)\|,\ 1-\|(x,y)\|\}>0$，则 $O((x,y),\delta)\subseteq S$。故所有非零且模长小于 $1$ 的点都是内点。原点 $(0,0)\notin S$，故不考虑。
因此 $\mathring{S}=S$（$S$ 中的每个点都是内点，但注意 $S$ 不包含原点）。

**闭包**：$\overline{S}=S\cup S'=S\cup\{(x,y)\mid x^2+y^2=1\}\cup\{(0,0)\}=\{(x,y)\mid x^2+y^2\leq 1\}$（闭单位圆盘）。

**开/闭性判断**：
- $S$ 是开集吗？是的。$S$ 中每个点都是内点（如上所述）。$S$ 实际上就是去掉圆心后的开圆盘，可以看作开球 $B(0,1)$ 去掉一个点。但这不影响它是开集——因为 $0\notin S$，不需要检验 $0$。
- $S$ 是闭集吗？不是。因为 $S$ 不包含其所有聚点——例如 $(0,0)$ 和边界点都是 $S$ 的聚点但不属于 $S$。

故 $S$ 是开集，不是闭集。

</details>

---

**6.** 在 $\mathbb{R}$ 中构造一个集合 $S$，使得 $S'=\mathbb{N}$（即 $S$ 的聚点恰好是全体自然数）。

<details><summary>参考答案</summary>

**思路**：要让每个 $n\in\mathbb{N}$ 成为 $S$ 的聚点，需要在 $n$ 附近放置趋于 $n$ 的序列。同时要确保非自然数不是聚点。

**构造**：对每个 $n\in\mathbb{N}$，取序列 $x_k^{(n)}=n+\frac1k$（$k\in\mathbb{N}^+$）。令
$$S=\bigcup_{n\in\mathbb{N}}\left\{n+\frac1k\;\Big|\; k\in\mathbb{N}^+\right\}$$

即 $S=\{\dots,\; 1+\frac12,1+\frac13,\dots,\; 2+\frac12,2+\frac13,\dots,\; \dots\}$。

**验证**：
- 对任意固定的 $n_0\in\mathbb{N}$，$x_k=n_0+\frac1k\to n_0$，且 $x_k\in S\setminus\{n_0\}$，故 $n_0$ 是聚点。
- 对 $x\notin\mathbb{N}$，考虑两种情况：
  - 若 $x$ 在 $(n,n+1)$ 之间且不是任何 $n+1/k$：取 $\delta=\min\{|x-(n+\frac1k)|,\ |x-(n+1+\frac1k)|\}$ 的一个下界，但这样的构造比较繁琐。
  
  更简单的方法：注意到 $S$ 中的每个点都是孤立点（因为 $n+1/k$ 到相邻点 $n+1/(k+1)$ 有正距离），且 $S$ 的聚点只有 $\mathbb{N}$。由于 $S$ 没有其他子列收敛到非自然数，故 $S'=\mathbb{N}$。

**另一种构造**：$S=\bigcup_{n\in\mathbb{N}}\{n+\frac1{2^k}\mid k\in\mathbb{N}^+\}$，思路相同。

</details>
