# 03. 幂级数的收敛半径——Abel定理与Cauchy-Hadamard公式

> 所属章节：第十章 函数项级数  |  文件序号：03  |  难度：进阶
> 常见混淆点：Abel定理（幂级数版本）与Ch9-05的Abel判别法是两个不同的定理，仅共享Abel之名；Cauchy-Hadamard公式中的 $\varlimsup$ 即使普通极限不存在仍然存在，处理系数有缺项的幂级数时必须使用上极限而非普通极限

## 1. 学习目标与先修前置

### 学习目标
- 掌握幂级数的形式定义 $\sum_{n=0}^{\infty} a_n (x - x_0)^n$，理解其作为函数项级数的特殊地位
- 掌握Abel定理（幂级数版本）的完整陈述：收敛方向与发散方向
- 理解Abel定理的证明：利用收敛必要条件给出等比控制，再通过比较判别法
- 严格区分Abel定理（幂级数版本）与Abel判别法（数项级数版本，Ch9-05）
- 掌握收敛半径的形式化定义 $R = \sup\{|x| : \sum a_n x^n \text{ 收敛}\}$
- 掌握Cauchy-Hadamard公式 $R = 1 / \varlimsup \sqrt[n]{|a_n|}$ 及其推导
- 能处理系数有缺项的幂级数（某些 $a_n = 0$），正确使用上极限计算半径
- 理解Cauchy-Hadamard公式与比值判别法的关系：比值法更简便但有条件，上极限法更普适

### 先修知识
- 文件 02（第十章）：函数项级数的定义（定义10.1-10.3）、收敛域（定义10.3）——幂级数是函数项级数的特例
- 文件 01-07（第九章）：数项级数的所有判别法，特别是：
  - 根值判别法的上极限形式（定理9.9）——Cauchy-Hadamard公式的直接理论基础
  - 比值判别法（定理9.7）——在极限存在时可用于求收敛半径
  - 比较判别法（定理9.4）——Abel定理证明的核心工具
  - 级数收敛的必要条件（定理9.1）：$\sum a_n$ 收敛 $\Rightarrow$ $\lim a_n = 0$
  - Abel判别法（文件05，定理9.13）——需与此处的Abel定理严格区分
- 文件 08-11（第二章）：数列极限的定义、上极限的定义（定义9.7）
- 文件 11（第二章）：$\lim q^n$ 的结论（定理11.4）——等比级数的敛散性

---

## 2. 背景与应用场景

在文件02中，我们学习了函数项级数的一般理论。现在聚焦于一类特别重要的函数项级数——**幂级数**（power series）：

$$\sum_{n=0}^{\infty} a_n (x - x_0)^n = a_0 + a_1 (x - x_0) + a_2 (x - x_0)^2 + \cdots + a_n (x - x_0)^n + \cdots$$

其中 $\{a_n\}$ 是常系数（与 $x$ 无关），$x_0$ 称为**展开中心**。幂级数的每一项都是幂函数 $(x - x_0)^n$ 乘以系数 $a_n$。
为简化讨论，通常通过平移变换 $x - x_0 \mapsto x$ 将中心化为 $x_0 = 0$，简写为：

$$\sum_{n=0}^{\infty} a_n x^n = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n + \cdots$$

幂级数的重要性体现在以下方面：
- **Taylor展开**：任意光滑函数都可以在一点附近展开为幂级数（如 $e^x = \sum x^n/n!$，$\sin x = \sum (-1)^n x^{2n+1}/(2n+1)!$）
- **微分方程求解**：许多微分方程的解可以表示为幂级数（如Bessel方程、Legendre方程）
- **复分析**：幂级数是全纯函数的基本表示法
- **逼近论**：幂级数的部分和给出了函数的多项式逼近

首次遇到幂级数时，最核心的问题是：**对于给定的系数 $\{a_n\}$，幂级数在哪些 $x$ 处收敛？答案具有什么样的结构？**

从第九章的经验出发，固定 $x$ 后，$\sum a_n x^n$ 就是一个数项级数，可以用比值判别法或根值判别法分析。但需要注意：比值/根值的结果依赖于 $x$。具体地：

- 比值法：$\displaystyle\lim_{n\to\infty} \left|\frac{a_{n+1} x^{n+1}}{a_n x^n}\right| = |x| \cdot \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right|$
- 根值法：$\displaystyle\varlimsup_{n\to\infty} \sqrt[n]{|a_n x^n|} = |x| \cdot \varlimsup_{n\to\infty} \sqrt[n]{|a_n|}$

两种方法都揭示了一个重要事实：**收敛条件归结为 $|x|$ 与某个仅由系数 $\{a_n\}$ 决定的常数的比较**。这说明幂级数的收敛域具有非常特殊的结构——它必然是以原点为中心的对称区间。本文件将系统建立这一结论的严格理论基础。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 示例 |
|------|------|------|
| $\sum_{n=0}^{\infty} a_n x^n$ | 幂级数，$a_n$ 为系数 | $\sum_{n=0}^{\infty} \frac{x^n}{n!}$ |
| $\sum_{n=0}^{\infty} a_n (x - x_0)^n$ | 以 $x_0$ 为中心的幂级数 | $\sum_{n=0}^{\infty} \frac{(x-1)^n}{2^n}$ |
| $R$ | 收敛半径 | $R = 1$ 表示在 $(-1, 1)$ 内收敛 |
| $\varlimsup$ | 上极限（limsup） | $\varlimsup \sqrt[n]{|a_n|}$ |
| $a_n$ | 幂级数第 $n$ 项的系数 | $a_n = 1/n!$ |
| $x_0$ | 展开中心 | 通常简化为 $x_0 = 0$ |

### 3.1 幂级数的定义

**定义 10.6（幂级数）**：设 $\{a_n\}_{n=0}^{\infty}$ 是一个实数列。称形如

$$\sum_{n=0}^{\infty} a_n (x - x_0)^n$$

的函数项级数为**幂级数**（power series），其中 $x_0 \in \mathbb{R}$ 称为**展开中心**，$a_n$ 称为**系数**。特别地，当 $x_0 = 0$ 时简写为 $\sum_{n=0}^{\infty} a_n x^n$。

幂级数的部分和为 $S_N(x) = \sum_{n=0}^{N} a_n (x - x_0)^n$，它是一个 $N$ 次多项式。

**注意**：幂级数的定义域 $D = \mathbb{R}$（每个 $a_n (x - x_0)^n$ 作为单项式定义在整个实数轴上），但收敛域 $C$（定义10.3）是 $\mathbb{R}$ 的子集。

---

## 4. Abel定理（幂级数版本）

Abel定理是幂级数理论中第一个核心定理。它揭示了幂级数收敛域的本质结构：**收敛点都在一个以原点为中心的区间内**。这与第九章的Abel判别法是两个不同的定理，详见4.4节的说明。

### 4.1 定理陈述

**定理 10.3（Abel定理——幂级数版本）**：设 $\sum_{n=0}^{\infty} a_n x^n$ 为幂级数。

1. **收敛方向**：若存在 $x_0 \neq 0$ 使 $\sum_{n=0}^{\infty} a_n x_0^n$ 收敛，则对一切满足 $|x| < |x_0|$ 的 $x$，幂级数 $\sum_{n=0}^{\infty} a_n x^n$ **绝对收敛**。
2. **发散方向**：若存在 $x_0$ 使 $\sum_{n=0}^{\infty} a_n x_0^n$ 发散，则对一切满足 $|x| > |x_0|$ 的 $x$，幂级数 $\sum_{n=0}^{\infty} a_n x^n$ **发散**。

**直观理解**：收敛的"好消息"会从收敛点向内传播——$x_0$ 处收敛意味着更靠近原点的所有点都绝对收敛。反之，发散的"坏消息"会向外传播——$x_0$ 处发散意味着更远离原点的所有点都发散。

**推论（收敛域的对称性）**：幂级数 $\sum a_n x^n$ 的收敛域 $C$ 是以原点为中心的对称区间。即，若 $x_0 \in C$，则对一切 $|x| < |x_0|$ 有 $x \in C$；若 $x_0 \notin C$，则对一切 $|x| > |x_0|$ 有 $x \notin C$。收敛域具体是 $(-R, R)$、$[-R, R]$、$(-R, R]$ 或 $[-R, R)$ 中的一种（$R$ 为收敛半径）。

---

### 4.2 收敛方向的证明

设 $\sum_{n=0}^{\infty} a_n x_0^n$ 收敛（$x_0 \neq 0$）。对任意满足 $|x| < |x_0|$ 的 $x$，我们要证明 $\sum_{n=0}^{\infty} a_n x^n$ **绝对收敛**，即 $\sum_{n=0}^{\infty} |a_n x^n| < \infty$。

**第 1 步：有界性**。

由于 $\sum_{n=0}^{\infty} a_n x_0^n$ 收敛，由级数收敛的必要条件（定理9.1），其通项趋于零：

$$\lim_{n\to\infty} a_n x_0^n = 0$$

收敛数列必有界（数列极限的基本性质）：存在常数 $M > 0$，使得

$$|a_n x_0^n| \le M, \quad \forall n = 0, 1, 2, \dots$$

**第 2 步：构造等比控制**。

对任意满足 $|x| < |x_0|$ 的 $x$，写出 $|a_n x^n|$ 的估计：

$$|a_n x^n| = |a_n x_0^n| \cdot \left|\frac{x}{x_0}\right|^n \le M \cdot \left|\frac{x}{x_0}\right|^n$$

其中第一步利用了恒等式 $a_n x^n = a_n x_0^n \cdot (x/x_0)^n$，第二步代入了有界性条件。

记 $q = \left|\dfrac{x}{x_0}\right|$。由于 $|x| < |x_0|$，有 $0 \le q < 1$。

**第 3 步：比较判别法**。

由第 2 步的放缩：

$$\sum_{n=0}^{\infty} |a_n x^n| \le \sum_{n=0}^{\infty} M \cdot q^n = M \cdot \sum_{n=0}^{\infty} q^n$$

右端是公比为 $q$（$0 \le q < 1$）的等比级数，由定理9.2，该级数收敛于 $M/(1-q)$。

由比较判别法（定理9.4(1)）：$|a_n x^n| \le M q^n$ 且 $\sum M q^n$ 收敛 $\Rightarrow$ $\sum |a_n x^n|$ 收敛。

因此 $\sum a_n x^n$ 在 $|x| < |x_0|$ 上绝对收敛。证毕。

**证明的关键步骤**：收敛 $\Rightarrow$ 通项有界 $\Rightarrow$ 等比放缩 $\Rightarrow$ 比较判别法。三步缺一不可。

---

### 4.3 发散方向的证明

设 $\sum_{n=0}^{\infty} a_n x_0^n$ 发散。对任意满足 $|x| > |x_0|$ 的 $x$，我们要证明 $\sum a_n x^n$ 发散。

采用**反证法**。

假设 $\sum a_n x^n$ 在某个满足 $|x| > |x_0|$ 的点 $x$ 处收敛。由于 $|x_0| < |x|$，将Abel定理的收敛方向应用于"反向"：以 $x$ 为收敛点，则对一切满足 $|t| < |x|$ 的 $t$，$\sum a_n t^n$ 绝对收敛。

特别地，取 $t = x_0$。由于 $|x_0| < |x|$，由收敛方向知 $\sum a_n x_0^n$ **绝对收敛**（从而收敛），与假设矛盾！

因此假设不成立，$\sum a_n x^n$ 在一切 $|x| > |x_0|$ 处发散。证毕。

**注意**：发散方向的证明完全依赖于收敛方向——如果没有收敛方向的结论，无法进行反证。

---

### 4.4 与Abel判别法的区别

**必须明确区分**：本文件的**Abel定理**（定理10.3）与第九章文件05的**Abel判别法**（定理9.13）是两个完全不同的定理，仅共享数学家Niels Henrik Abel的名字。

| 性质 | Abel定理（本文件） | Abel判别法（Ch9-05） |
|------|-------------------|---------------------|
| 对象 | 幂级数 $\sum a_n x^n$ | 数项级数 $\sum a_n b_n$ |
| 条件 | 在某点 $x_0$ 收敛/发散 | $\{a_n\}$ 单调有界 + $\sum b_n$ 收敛 |
| 结论 | 更小/更大范围内的敛散性 | $\sum a_n b_n$ 收敛 |
| 证明方法 | 比较判别法 + 有界性 | 分部求和法（Abel变换） |

后续论及"Abel定理"时，如无特别说明，均指本文件的定理10.3。

---

## 5. 收敛半径的形式化定义

Abel定理揭示了幂级数收敛域的本质结构：所有收敛点构成一个以原点为中心的对称区间。这个区间的"宽度"由**收敛半径** $R$ 刻画。

### 5.1 定义

**定义 10.7（收敛半径）**：设 $\sum_{n=0}^{\infty} a_n x^n$ 为幂级数。定义

$$R = \sup\left\{ |x| \;\Big|\; \sum_{n=0}^{\infty} a_n x^n \text{ 收敛} \right\}$$

称 $R$ 为幂级数的**收敛半径**（radius of convergence）。约定：

- 若集合 $\{|x| : \sum a_n x^n \text{ 收敛}\}$ 无上界（即对一切实数 $x$ 级数收敛），则 $R = +\infty$；
- 若集合中仅 $x = 0$ 为收敛点（即 $\sup = 0$），则 $R = 0$。

**几何意义**：

根据 $R$ 的取值，收敛域的结构分为三种情形：

1. **$R = 0$**：幂级数仅在 $x = 0$ 处收敛（此时收敛域就是 $\{0\}$）。
2. **$R = +\infty$**：幂级数在全体实数上收敛（收敛域为 $\mathbb{R}$）。
3. **$0 < R < +\infty$**：幂级数在开区间 $(-R, R)$ 内**绝对收敛**，在 $(-R, R)$ 外**发散**。在端点 $x = \pm R$ 处，级数可能收敛也可能发散，需单独讨论。

**重要说明**：收敛半径 $R$ 的定义仅依赖于幂级数的系数，不依赖于端点的敛散性——端点是否需要包含在收敛域中，需要单独判断。

### 5.2 Abel定理保证收敛半径的well-defined性

定义10.7中的上确界一定存在（可能为 $+\infty$），这是数学上平凡的——任何实数集合都有上确界。但**Abel定理保证了收敛半径概念是有意义的**，即：收敛域确实具有"从原点向外辐射"的区间结构。

具体来说，Abel定理保证了以下性质：

- **对一切 $|x| < R$**：由 $R$ 的定义（上确界），存在收敛点 $x_1$ 满足 $|x| < |x_1| \le R$（或 $x_1$ 满足 $|x| < |x_1|$ 且 $|x_1|$ 任意接近 $R$）。由Abel定理的收敛方向，$\sum a_n x^n$ 绝对收敛。
- **对一切 $|x| > R$**：由 $R$ 的定义，不存在收敛点其绝对值大于 $R$。因此任何 $|x| > R$ 都是发散点。由Abel定理的发散方向，这进一步确认了发散性质。

**如果没有Abel定理**，收敛域可能是"散落"的——即可能有一些收敛点零散分布在数轴上，而中间夹着发散点。Abel定理排除了这种可能性，使收敛半径成为一个有意义的几何量。

### 5.3 三种情形的举例

**(1) $R = 0$ 的情形**：考虑 $\sum_{n=0}^{\infty} n! \, x^n$。对 $x \neq 0$，通项 $|n! x^n| \to \infty$（$n \to \infty$），不满足收敛的必要条件。仅 $x = 0$ 时级数收敛（和为 $0! \cdot 0^0 = 1$，约定 $0^0 = 1$）。收敛半径 $R = 0$。

**(2) $R = +\infty$ 的情形**：考虑 $\sum_{n=0}^{\infty} \frac{x^n}{n!}$。对任意固定的 $x$，由比值判别法：
$$\lim_{n\to\infty} \left|\frac{a_{n+1} x^{n+1}}{a_n x^n}\right| = \lim_{n\to\infty} \frac{|x|}{n+1} = 0 < 1$$
故级数对一切 $x$ 收敛。收敛半径 $R = +\infty$。

**(3) $0 < R < +\infty$ 的情形**：考虑 $\sum_{n=0}^{\infty} x^n$。由等比级数知识：
- $|x| < 1$ 时收敛于 $1/(1-x)$
- $|x| \ge 1$ 时发散
故收敛半径 $R = 1$，收敛域为 $(-1, 1)$（端点均发散）。

---

## 6. Cauchy-Hadamard公式

定义10.7给出了收敛半径的概念定义，但如何从系数 $\{a_n\}$ 直接计算出 $R$？Cauchy-Hadamard公式给出了答案。

### 6.1 定理陈述

**定理 10.4（Cauchy-Hadamard公式）**：设 $\sum_{n=0}^{\infty} a_n x^n$ 为幂级数，记

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|}$$

则收敛半径 $R$ 由下式给出：

$$R = \begin{cases}
\displaystyle \frac{1}{\rho}, & 0 < \rho < +\infty \\[10pt]
+\infty, & \rho = 0 \\[10pt]
0, & \rho = +\infty
\end{cases}$$

其中 $\varlimsup$ 表示上极限（定义9.7），并约定 $1/0 = +\infty$ 和 $1/(+\infty) = 0$。

---

### 6.2 证明（从根值判别法出发）

定理的证明基于第九章的根值判别法（定理9.9——上极限版本）。对每个固定的 $x$，$\sum_{n=0}^{\infty} a_n x^n$ 是数项级数，可应用根值判别法。

**情形 1：$0 < \rho < +\infty$**。

对固定的 $x$，计算通项 $a_n x^n$ 的 $n$ 次根的上极限：

$$\varlimsup_{n\to\infty} \sqrt[n]{|a_n x^n|} = \varlimsup_{n\to\infty} \left( \sqrt[n]{|a_n|} \cdot |x| \right) = |x| \cdot \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = |x| \cdot \rho$$

其中第二个等式利用了上极限的性质：若 $c \ge 0$ 为常数，则 $\varlimsup (c \cdot b_n) = c \cdot \varlimsup b_n$（当 $c = 0$ 时平凡；$c > 0$ 时由 $\sup$ 的定义可直接验证）。

由上极限版本的根值判别法（定理9.9）：

- 若 $|x| \cdot \rho < 1$，即 $|x| < 1/\rho$，则 $\sum a_n x^n$ **收敛**（实际上绝对收敛）。
- 若 $|x| \cdot \rho > 1$，即 $|x| > 1/\rho$，则 $\sum a_n x^n$ **发散**。
- 若 $|x| \cdot \rho = 1$，即 $|x| = 1/\rho$，无法判定（需要单独讨论端点）。

因此收敛半径 $R = 1/\rho = 1 / \varlimsup \sqrt[n]{|a_n|}$。

**情形 2：$\rho = 0$**。

此时 $\varlimsup \sqrt[n]{|a_n|} = 0$。对任意 $x \neq 0$：

$$\varlimsup_{n\to\infty} \sqrt[n]{|a_n x^n|} = |x| \cdot 0 = 0 < 1$$

由根值判别法（定理9.9(1)），$\sum a_n x^n$ 对一切 $x$ 收敛（绝对收敛）。故 $R = +\infty$。

**情形 3：$\rho = +\infty$**。

此时 $\varlimsup \sqrt[n]{|a_n|} = +\infty$。对任意 $x \neq 0$：

$$\varlimsup_{n\to\infty} \sqrt[n]{|a_n x^n|} = |x| \cdot (+\infty) = +\infty > 1$$

由根值判别法（定理9.9(2)），$\sum a_n x^n$ 在 $x \neq 0$ 时发散。仅 $x = 0$ 处级数收敛（所有项为零）。故 $R = 0$。

证毕。

**证明的实质**：Cauchy-Hadamard公式本质上就是**根值判别法在幂级数上的直接应用**。固定 $x$ 后，数项级数的根值判别法自动给出收敛条件 $|x| < 1/\varlimsup \sqrt[n]{|a_n|}$。

---

### 6.3 系数有缺项时的处理

当系数序列 $\{a_n\}$ 中有许多项为零时（即幂级数"缺项"），$\sqrt[n]{|a_n|}$ 的普通极限可能不存在，但**上极限**依然存在且有效。这正是Cauchy-Hadamard公式相对于比值判别法的核心优势所在。

**例**：考虑幂级数 $\sum_{n=0}^{\infty} a_n x^n$，其中

$$a_n = \begin{cases}
3^n, & n = 2k \text{（偶数）} \\
0, & n = 2k+1 \text{（奇数）}
\end{cases}$$

即级数为 $\sum_{k=0}^{\infty} 3^{2k} x^{2k} = \sum_{k=0}^{\infty} (3^2 x^2)^k = \sum_{k=0}^{\infty} (9x^2)^k$。

直接写出幂级数形式：
- $a_0 = 3^0 = 1$，$a_1 = 0$，$a_2 = 3^2 = 9$，$a_3 = 0$，$a_4 = 3^4 = 81$，$\dots$

计算 $\sqrt[n]{|a_n|}$：
- 当 $n$ 为奇数时（$a_n = 0$）：$\sqrt[n]{|a_n|} = 0$
- 当 $n$ 为偶数时（$n = 2k$）：$\sqrt[n]{|a_n|} = \sqrt[2k]{3^{2k}} = 3$

因此 $\sqrt[n]{|a_n|}$ 的子列极限为 $0$（奇数项）和 $3$（偶数项）。上极限取子列极限的最大值：

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = 3$$

由Cauchy-Hadamard公式，收敛半径 $R = 1/3$。

**验证**：原级数 $\sum_{k=0}^{\infty} (9x^2)^k$ 是公比为 $9x^2$ 的等比级数，当 $|9x^2| < 1$ 即 $|x| < 1/3$ 时收敛。这与 $R = 1/3$ 一致。

**为什么普通极限不存在**：$\sqrt[n]{|a_n|}$ 的奇数项为 $0$，偶数项为 $3$，两个子列收敛到不同的值，因此普通极限不存在。但上极限总是存在的（任何数列都有上极限，可为 $+\infty$）。

---

### 6.4 Cauchy-Hadamard公式的简化记法

为方便使用，Cauchy-Hadamard公式通常合并写为：

$$R = \frac{1}{\displaystyle\varlimsup_{n\to\infty} \sqrt[n]{|a_n|}}$$

并约定 $\dfrac{1}{0} = +\infty$，$\dfrac{1}{+\infty} = 0$。这种写法简洁地覆盖了三种情形。

---

## 7. 比值判别法与Cauchy-Hadamard公式的关系

当根值法中的普通极限 $\lim \sqrt[n]{|a_n|}$ 存在时，比值法可以作为Cauchy-Hadamard公式的一个替代（且计算更简便）。

### 7.1 比值法的可用条件

对于幂级数 $\sum a_n x^n$，若极限

$$\ell = \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right|$$

存在（允许 $\ell = +\infty$），则收敛半径也可由比值法求得：

$$R = \frac{1}{\ell}$$

并约定 $\ell = 0$ 时 $R = +\infty$，$\ell = +\infty$ 时 $R = 0$。

**推导**：固定 $x \neq 0$，对 $\sum a_n x^n$ 应用比值判别法（定理9.7）：

$$\lim_{n\to\infty} \left|\frac{a_{n+1} x^{n+1}}{a_n x^n}\right| = |x| \cdot \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = |x| \cdot \ell$$

收敛条件为 $|x| \cdot \ell < 1$，即 $|x| < 1/\ell$。因此 $R = 1/\ell$。

**注意**：比值法要求 $\lim |a_{n+1}/a_n|$ **存在**。若该极限不存在（例如系数振荡时），比值法失效，必须使用Cauchy-Hadamard公式（上极限法始终有效）。

### 7.2 两个命题的关系

关于两种方法的关系，有以下重要结论：

**命题**：若 $\lim_{n\to\infty} \left|\dfrac{a_{n+1}}{a_n}\right|$ 存在（为 $\ell$），则 $\varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \ell$（即Cauchy-Hadamard公式与比值法给出相同的 $R$）。

这一结论来自数学分析中的Cauchy命题：若 $\lim a_{n+1}/a_n = \ell$，则 $\lim \sqrt[n]{a_n} = \ell$（对正数列成立）。因此两种方法本质一致，但比值法仅适用于极限存在的情形，而上极限法则普适。

**何时使用比值法，何时使用上极限法？**

| 情形 | 推荐方法 | 示例 |
|------|---------|------|
| 系数含阶乘 $n!$ | 比值法（消去阶乘） | $\sum \frac{x^n}{n!}$，$\ell = 0$，$R = +\infty$ |
| 系数含指数 $a^n$ | 比值法或根值法 | $\sum \frac{x^n}{3^n}$，$\ell = 1/3$，$R = 3$ |
| 系数有缺项（大量 $a_n = 0$） | Cauchy-Hadamard公式（上极限） | $\sum 3^{2k} x^{2k}$，$\rho = 3$，$R = 1/3$ |
| 系数振荡（如 $[2+(-1)^n]^n$） | Cauchy-Hadamard公式（上极限） | $\sum [2+(-1)^n]^n x^n$，需上极限 |
| 系数为有理函数 | 比值法 | $\sum \frac{x^n}{n^2+1}$，$\ell = 1$，$R = 1$ |

**核心原则**：当极限 $\lim |a_{n+1}/a_n|$ 存在时，比值法计算更简便；当该极限不存在或不确定时，使用Cauchy-Hadamard公式（上极限法）。

---

## 8. 例题

### 例题 1：比值法与Cauchy-Hadamard公式的等价使用

求幂级数 $\displaystyle\sum_{n=1}^{\infty} \frac{3^n}{\sqrt{n}} x^n$ 的收敛半径，并讨论端点 $x = \pm\dfrac{1}{3}$ 处的敛散性，写出完整的收敛域。

**解**：

**第 1 步：判断使用哪种方法**。

系数 $a_n = \dfrac{3^n}{\sqrt{n}}$。极限 $\displaystyle\lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right|$ 存在，可用比值法或Cauchy-Hadamard公式。比值法计算更简便。

**第 2 a 步（比值法）**：

$$\lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{3^{n+1}/\sqrt{n+1}}{3^n/\sqrt{n}} = \lim_{n\to\infty} 3 \cdot \sqrt{\frac{n}{n+1}} = 3 \cdot 1 = 3$$

因此 $\ell = 3$，收敛半径 $R = 1/\ell = 1/3$。

**第 2 b 步（Cauchy-Hadamard公式验证）**：

$$\sqrt[n]{|a_n|} = \sqrt[n]{\frac{3^n}{\sqrt{n}}} = \frac{3}{n^{1/(2n)}}$$

由于 $\displaystyle\lim_{n\to\infty} n^{1/(2n)} = \lim_{n\to\infty} e^{\frac{\ln n}{2n}} = e^0 = 1$，故 $\rho = \varlimsup \sqrt[n]{|a_n|} = 3$，$R = 1/3$。两种方法结果一致。

**第 3 步：讨论端点**。

- 在 $x = \dfrac{1}{3}$ 处：级数化为 $\displaystyle\sum_{n=1}^{\infty} \frac{3^n}{\sqrt{n}} \cdot \left(\frac{1}{3}\right)^n = \sum_{n=1}^{\infty} \frac{1}{\sqrt{n}} = \sum_{n=1}^{\infty} \frac{1}{n^{1/2}}$

  这是 $p = 1/2 \le 1$ 的 $p$-级数，由定理9.5，该级数发散。

- 在 $x = -\dfrac{1}{3}$ 处：级数化为 $\displaystyle\sum_{n=1}^{\infty} \frac{3^n}{\sqrt{n}} \cdot \left(-\frac{1}{3}\right)^n = \sum_{n=1}^{\infty} \frac{(-1)^n}{\sqrt{n}}$

  这是交错级数，通项 $\dfrac{1}{\sqrt{n}}$ 单调递减趋于 $0$。由Leibniz判别法（定理9.11），该级数收敛（条件收敛）。

**第 4 步：写出收敛域**。

收敛域为 $[-\dfrac{1}{3}, \dfrac{1}{3})$。

---

### 例题 2：系数有缺项的幂级数

求幂级数 $\displaystyle\sum_{n=1}^{\infty} (-1)^n \frac{x^{2n}}{n \cdot 2^n}$ 的收敛半径。

**解**：

**第 1 步：化为标准形式**。

原级数只有 $x$ 的偶次幂。写成 $\sum a_n x^n$ 的形式：

$$a_n = \begin{cases}
(-1)^{n/2} \cdot \dfrac{1}{(n/2) \cdot 2^{n/2}}, & n = 2k \text{（偶数）} \\[6pt]
0, & n = 2k+1 \text{（奇数）}
\end{cases}$$

即 $a_{2k} = (-1)^k \cdot \dfrac{1}{k \cdot 2^k}$，$a_{2k+1} = 0$。

**第 2 步：用Cauchy-Hadamard公式**。

计算 $\sqrt[n]{|a_n|}$：

- 当 $n$ 为奇数时：$\sqrt[n]{|a_n|} = \sqrt[n]{0} = 0$。
- 当 $n = 2k$ 为偶数时：
  $$\sqrt[n]{|a_n|} = \sqrt[2k]{\frac{1}{k \cdot 2^k}} = \frac{1}{\sqrt[2k]{k} \cdot \sqrt{2}}$$

**第 3 步：求上极限**。

考察 $\sqrt[2k]{k}$ 的极限：

$$\lim_{k\to\infty} \sqrt[2k]{k} = \lim_{k\to\infty} e^{\frac{\ln k}{2k}} = e^0 = 1$$

因此当 $n = 2k$ 沿偶数趋于无穷时：

$$\lim_{k\to\infty} \sqrt[n]{|a_n|} = \frac{1}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}$$

奇数项子列的极限为 $0$。上极限取子列极限的最大值：

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \frac{1}{\sqrt{2}}$$

**第 4 步：计算收敛半径**。

$$R = \frac{1}{\rho} = \sqrt{2}$$

**验证**：将原级数改写为 $\displaystyle\sum_{k=1}^{\infty} (-1)^k \frac{(x^2)^k}{k \cdot 2^k}$，令 $t = x^2$，则化为 $\displaystyle\sum_{k=1}^{\infty} (-1)^k \frac{t^k}{k \cdot 2^k}$。对 $t$ 的级数，由比值法得收敛半径 $R_t = 2$。因此 $x^2 < 2$ 即 $|x| < \sqrt{2}$ 时收敛。与原结果一致。

---

### 例题 3：$R = 0$ 和 $R = +\infty$ 的情形

用Cauchy-Hadamard公式求下列幂级数的收敛半径：

(1) $\displaystyle\sum_{n=1}^{\infty} n^n x^n$

(2) $\displaystyle\sum_{n=1}^{\infty} \frac{1}{n^n} x^n$

**解**：

**(1)** $a_n = n^n$。

$$\sqrt[n]{|a_n|} = \sqrt[n]{n^n} = n$$

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \lim_{n\to\infty} n = +\infty$$

由Cauchy-Hadamard公式，$\rho = +\infty \Rightarrow R = 0$。

**含义**：该幂级数仅在 $x = 0$ 处收敛。对任意 $x \neq 0$，当 $n$ 充分大时 $|n^n x^n| = (n|x|)^n \to \infty$，通项不趋于零，级数发散。

---

**(2)** $a_n = \dfrac{1}{n^n}$。

$$\sqrt[n]{|a_n|} = \sqrt[n]{\frac{1}{n^n}} = \frac{1}{n}$$

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \lim_{n\to\infty} \frac{1}{n} = 0$$

由Cauchy-Hadamard公式，$\rho = 0 \Rightarrow R = +\infty$。

**含义**：该幂级数对一切实数 $x$ 收敛（实际上对一切 $x$ 绝对收敛）。这是因为分母 $n^n$ 衰减极快，比任何指数衰减 $r^n$ 都快，即使是很大的 $x$ 也无法阻止通量趋于零。

---

## 9. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 幂级数 $\sum a_n x^n$ 在 $x = x_0$ 收敛，则对 $|x| \ge |x_0|$ 也收敛 | 收敛方向只能推出 $|x| < |x_0|$ 时收敛。对 $|x| > |x_0|$ 的情况不能下结论（实际上可能发散，也可能发散方向告诉我们发散）——反例：$\sum (-1)^n x^n/n$ 在 $x=1$ 收敛，在 $x=2$ 发散 |
| Abel定理（幂级数版本）和Abel判别法（数项级数版本）是同一回事 | 两者完全不同。幂级数版本的Abel定理研究 $\sum a_n x^n$ 的收敛域结构；数项级数版本的Abel判别法研究 $\sum a_n b_n$ 的敛散性。仅共享Abel之名 |
| Cauchy-Hadamard公式中可以用 $\lim \sqrt[n]{|a_n|}$ 代替 $\varlimsup \sqrt[n]{|a_n|}$ | 仅当普通极限存在时可以。对于系数有缺项或振荡的情况，普通极限不存在，必须使用上极限。用普通极限会得到错误结果 |
| 比值法求得的 $R$ 与Cauchy-Hadamard公式求得的 $R$ 一定相等 | 当 $\lim |a_{n+1}/a_n|$ 存在时相等。但当该极限不存在时，比值法失效（不能使用），而Cauchy-Hadamard公式仍有效 |
| 收敛半径 $R$ 就是"级数收敛的最大区间半径" | 从定义上看 $R = \sup\{|x| : \text{收敛}\}$，这个上确界就是收敛区间的半径。但端点是否包含在收敛域中，需要单独判断，$R$ 本身不决定端点行为 |
| $\varlimsup \sqrt[n]{|a_n|} = \rho$ 中，$\rho$ 的取值仅与 $|a_n|$ 的"大小"有关 | 上极限捕捉的是系数序列的**子列极限的最大值**，而不是最大项。例如缺项级数中 $a_n$ 可能是 $0$（小）和 $3^n$（大）交替出现，上极限 $3$ 来自偶数子列 |

### 检查点

- [ ] 能否完整叙述Abel定理（定理10.3）的收敛方向和发散方向？
- [ ] 能否完成Abel定理收敛方向的完整证明（有界性 $\rightarrow$ 等比控制 $\rightarrow$ 比较判别法）？
- [ ] 能否说明Abel定理发散方向的证明为什么依赖于收敛方向（反证法）？
- [ ] 能否明确区分本文件的Abel定理与Ch9-05的Abel判别法？
- [ ] 能否写出收敛半径的形式化定义 $R = \sup\{|x| : \sum a_n x^n \text{ 收敛}\}$？
- [ ] 能否解释为什么Abel定理保证了收敛半径概念的well-defined性？
- [ ] 能否写出Cauchy-Hadamard公式（定理10.4）并完成从根值判别法的推导？
- [ ] 能否正确处理 $\rho = 0$ 和 $\rho = +\infty$ 两种特殊情形（$R$ 分别取 $+\infty$ 和 $0$）？
- [ ] 能否计算系数有缺项的幂级数的收敛半径（正确使用上极限）？
- [ ] 能否说明比值法求 $R$ 的条件和局限性，以及它与Cauchy-Hadamard公式的关系？
- [ ] 对于给定的幂级数，能否建立正确的计算流程：判断极限是否存在 $\rightarrow$ 选择比值法或Cauchy-Hadamard公式 $\rightarrow$ 计算 $R$ $\rightarrow$ 讨论端点？

---

## 练习题

### 基础巩固

**1.** 求下列幂级数的收敛半径，并讨论端点处的敛散性。

(1) $\displaystyle\sum_{n=1}^{\infty} \frac{x^n}{n^2}$

(2) $\displaystyle\sum_{n=1}^{\infty} \frac{2^n}{n} x^n$

(3) $\displaystyle\sum_{n=1}^{\infty} \frac{n!}{n^n} x^n$

<details><summary>参考答案</summary>

**(1)** $\displaystyle\sum_{n=1}^{\infty} \frac{x^n}{n^2}$

系数 $a_n = 1/n^2$。由比值法：
$$\ell = \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{1/(n+1)^2}{1/n^2} = \lim_{n\to\infty} \frac{n^2}{(n+1)^2} = 1$$
收敛半径 $R = 1/\ell = 1$。

端点讨论：
- $x = 1$：级数化为 $\sum 1/n^2$，$p = 2 > 1$ 的 $p$-级数，收敛（绝对收敛）。
- $x = -1$：级数化为 $\sum (-1)^n / n^2$，绝对收敛（因为 $|(-1)^n / n^2| = 1/n^2$ 收敛）。

因此收敛域为 $[-1, 1]$。

---

**(2)** $\displaystyle\sum_{n=1}^{\infty} \frac{2^n}{n} x^n$

系数 $a_n = 2^n / n$。由比值法：
$$\ell = \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{2^{n+1}/(n+1)}{2^n/n} = \lim_{n\to\infty} 2 \cdot \frac{n}{n+1} = 2$$
收敛半径 $R = 1/2$。

端点讨论：
- $x = 1/2$：级数化为 $\sum 2^n / n \cdot (1/2)^n = \sum 1/n$（调和级数），发散。
- $x = -1/2$：级数化为 $\sum 2^n / n \cdot (-1/2)^n = \sum (-1)^n / n$（交错调和级数），由Leibniz判别法，收敛（条件收敛）。

因此收敛域为 $[-1/2, 1/2)$。

---

**(3)** $\displaystyle\sum_{n=1}^{\infty} \frac{n!}{n^n} x^n$

系数 $a_n = n! / n^n$。由比值法：
$$\ell = \lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{(n+1)!/(n+1)^{n+1}}{n!/n^n} = \lim_{n\to\infty} \frac{(n+1)}{(n+1)^{n+1}} \cdot n^n = \lim_{n\to\infty} \frac{n^n}{(n+1)^n} = \lim_{n\to\infty} \frac{1}{(1+1/n)^n} = \frac{1}{e}$$
收敛半径 $R = e$。

端点讨论（$x = \pm e$）需要使用更精细的工具（如Raabe判别法），超出本课程范围。但可以判断：在 $|x| = e$ 处，通项 $\frac{n!}{n^n} \cdot e^n$ 的行为与Stirling公式有关，此处略去。

</details>

---

**2.** 用Cauchy-Hadamard公式求下列幂级数的收敛半径。

(1) $\displaystyle\sum_{n=1}^{\infty} \left(1 + \frac{1}{n}\right)^{n^2} x^n$

(2) $\displaystyle\sum_{n=1}^{\infty} \frac{[3 + (-1)^n]^n}{2^n} x^n$

<details><summary>参考答案</summary>

**(1)** $\displaystyle\sum_{n=1}^{\infty} \left(1 + \frac{1}{n}\right)^{n^2} x^n$

系数 $a_n = \left(1 + \frac{1}{n}\right)^{n^2}$。用Cauchy-Hadamard公式：
$$\sqrt[n]{|a_n|} = \sqrt[n]{\left(1 + \frac{1}{n}\right)^{n^2}} = \left(1 + \frac{1}{n}\right)^n$$

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \lim_{n\to\infty} \left(1 + \frac{1}{n}\right)^n = e$$

收敛半径 $R = 1/e$。

(教材中常用极限 $\lim (1+1/n)^n = e$，定理12.1)

---

**(2)** $\displaystyle\sum_{n=1}^{\infty} \frac{[3 + (-1)^n]^n}{2^n} x^n$

系数 $a_n = \dfrac{[3 + (-1)^n]^n}{2^n}$。计算 $n$ 次根：
$$\sqrt[n]{|a_n|} = \sqrt[n]{\frac{[3 + (-1)^n]^n}{2^n}} = \frac{3 + (-1)^n}{2}$$

$\sqrt[n]{|a_n|}$ 的取值依赖于 $n$ 的奇偶：
- 偶数 $n$：$3 + (-1)^n = 4$，$\sqrt[n]{|a_n|} = 2$
- 奇数 $n$：$3 + (-1)^n = 2$，$\sqrt[n]{|a_n|} = 1$

子列极限分别为 $2$ 和 $1$，上极限取最大值：
$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = 2$$

收敛半径 $R = 1/2$。

注意：普通极限 $\lim \sqrt[n]{|a_n|}$ 不存在（振荡于 $1$ 和 $2$ 之间），必须使用上极限版本的Cauchy-Hadamard公式。

</details>

---

### 迁移应用

**3.** 考虑幂级数 $\displaystyle\sum_{n=0}^{\infty} a_n x^n$，其中系数为

$$a_n = \begin{cases}
2^n, & n = 3k \text{（3的倍数）} \\
0, & \text{其他}
\end{cases}$$

(1) 写出该幂级数的具体形式（写出前若干非零项）。

(2) 用Cauchy-Hadamard公式求其收敛半径 $R$，写出 $\varlimsup \sqrt[n]{|a_n|}$ 的完整计算过程。

(3) 解释为什么比值判别法在此失效。

<details><summary>参考答案</summary>

**(1)** 幂级数的具体形式：

$$a_0 = 2^0 = 1, \quad a_1 = a_2 = 0, \quad a_3 = 2^3 = 8, \quad a_4 = a_5 = 0, \quad a_6 = 2^6 = 64, \dots$$

因此幂级数为：

$$1 + 8x^3 + 64x^6 + \cdots = \sum_{k=0}^{\infty} 2^{3k} x^{3k} = \sum_{k=0}^{\infty} (8x^3)^k$$

---

**(2)** 用Cauchy-Hadamard公式求收敛半径。

计算 $\sqrt[n]{|a_n|}$：
- 当 $n$ 为3的倍数时（$n = 3k$）：$\sqrt[n]{|a_n|} = \sqrt[3k]{2^{3k}} = 2$
- 当 $n$ 不是3的倍数时：$\sqrt[n]{|a_n|} = \sqrt[n]{0} = 0$

因此 $\sqrt[n]{|a_n|}$ 的子列极限为 $0$（非3倍数项）和 $2$（3倍数项）。上极限取最大值：

$$\rho = \varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = 2$$

收敛半径：

$$R = \frac{1}{\rho} = \frac{1}{2}$$

**验证**：原级数 $\sum_{k=0}^{\infty} (8x^3)^k$ 是公比为 $8x^3$ 的等比级数，收敛条件为 $|8x^3| < 1$，即 $|x| < 1/2$。与 $R = 1/2$ 一致。

---

**(3)** 比值判别法为什么失效？

考察比值 $\left|\dfrac{a_{n+1}}{a_n}\right|$：由于大量 $a_n = 0$，比值中分母为零的情况出现时（$a_n = 0$ 而 $a_{n+1} \neq 0$），比值无定义。即使考虑非零项的比值：

$$\frac{a_{3k+3}}{a_{3k}} = \frac{2^{3k+3}}{2^{3k}} = 8$$

但 $\dfrac{a_{3k+1}}{a_{3k}} = \dfrac{0}{2^{3k}} = 0$，$\dfrac{a_{3k+2}}{a_{3k+1}}$ 无定义。$\left|\dfrac{a_{n+1}}{a_n}\right|$ 不存在极限（甚至无法对整个数列定义）。因此比值判别法失效。

Cauchy-Hadamard公式（上极限法）不受此影响——它只关心非零系数的 $n$ 次根，忽略 $a_n = 0$ 的项（其贡献为 $0$，不影响上极限）。

</details>

---

**4.** （综合）设幂级数 $\sum_{n=0}^{\infty} a_n x^n$ 的收敛半径为 $R$（$0 < R < +\infty$）。证明：

(1) 若 $R > 0$，则级数在 $(-R, R)$ 内绝对收敛。

(2) 若 $R$ 为收敛半径，则 $\varlimsup_{n\to\infty} \sqrt[n]{|a_n|} = \dfrac{1}{R}$。

<details><summary>参考答案</summary>

**(1)** 级数在 $(-R, R)$ 内绝对收敛。

由收敛半径的定义 $R = \sup\{|x| : \sum a_n x^n \text{ 收敛}\}$，对任意 $x$ 满足 $|x| < R$，存在 $x_1$ 使得 $|x| < |x_1| \le R$ 且 $\sum a_n x_1^n$ 收敛（根据上确界的性质：对任意 $\varepsilon > 0$，存在收敛点 $x_1$ 使得 $|x_1| > R - \varepsilon$）。

由Abel定理的收敛方向（定理10.3(1)）：$\sum a_n x_1^n$ 收敛 $\Rightarrow$ 对一切 $|t| < |x_1|$，$\sum a_n t^n$ 绝对收敛。取 $t = x$ 即得 $\sum a_n x^n$ 绝对收敛。

---

**(2)** $\varlimsup \sqrt[n]{|a_n|} = 1/R$。

由Cauchy-Hadamard公式（定理10.4）：$R = 1 / \varlimsup \sqrt[n]{|a_n|}$。两边取倒数即得 $\varlimsup \sqrt[n]{|a_n|} = 1/R$。

注意：当 $R$ 已知时，这给出了从收敛半径反推系数的上极限的公式。

</details>
