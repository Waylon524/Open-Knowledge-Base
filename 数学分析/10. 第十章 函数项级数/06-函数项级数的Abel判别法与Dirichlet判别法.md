# 06. 函数项级数的 Abel 判别法与 Dirichlet 判别法

> 所属章节：第十章 函数项级数  |  文件序号：06  |  难度：进阶
> 常见混淆点：函数项级数版本的 Abel/Dirichlet 判别法是将 Ch9-05 的"收敛"结论升级为"一致收敛"，因此条件也必须同步升级——$a_n \to 0$ 要变成 $a_n(x) \rightrightarrows 0$，"部分和有界"要变成"一致有界"，"$\sum b_n$ 收敛"要变成"$\sum b_n(x)$ 一致收敛"；三角级数 $1/|\sin(x/2)|$ 的上界在逐点收敛时可直接用（依赖于 $x$），但在一致收敛场景中必须通过限制区间使该上界与 $x$ 无关

## 1. 学习目标与先修前置

### 学习目标
- 掌握函数项级数的一致 Cauchy 收敛准则（定理 10.12），能完成双向证明并理解其与数项级数 Cauchy 准则的区别
- 掌握 Dirichlet 判别法（函数项级数版本，定理 10.13）的定理陈述与完整证明，理解其与 Ch9-05 数项级数版本的条件升级关系
- 掌握 Abel 判别法（函数项级数版本，定理 10.14）的定理陈述与两种证明方法（归化法 + 直接 Abel 变换法）
- 掌握 $\sum_{k=1}^n \cos(kx)$ 和 $\sum_{k=1}^n \sin(kx)$ 的闭式公式，能推导其在 $[\delta, 2\pi-\delta]$ 上的**一致**有界性（区别于 Ch9-05 的逐点有界）
- 能运用 Dirichlet 判别法和 Abel 判别法判断具体函数项级数的一致收敛性，能解释为何 $\delta > 0$ 是必要条件

### 先修知识
- 文件 02（第十章）：一致收敛的 $\varepsilon$-$N$ 定义（定义 10.5）和 $\sup$ 等价刻画（定理 10.1）
- 文件 04（第十章）：一致收敛函数项级数的连续性定理（定理 10.5）——理解一致收敛对和函数分析性质的影响
- 文件 05（第九章）：Abel 变换（引理 9.13，三种形式）、Dirichlet 判别法（定理 9.14）与 Abel 判别法（定理 9.15）的数项级数版本——这些是当前函数项级数版本的基础
- 文件 05（第九章）：三角恒等式推导 $\sum \sin(k\theta)$ 和 $\sum \cos(k\theta)$ 的闭式（公式 9.14.1/9.14.2）及其有界性上界 $1/|\sin(\theta/2)|$（公式 9.14.3）
- 三角恒等式：$2\sin\alpha\sin\beta = \cos(\alpha-\beta) - \cos(\alpha+\beta)$、$2\sin\alpha\cos\beta = \sin(\alpha+\beta) + \sin(\alpha-\beta)$

---

## 2. 背景与应用场景

### 2.1 从数项级数到函数项级数的"一致化"升级

在第九章中，我们学习了数项级数的 Abel 判别法和 Dirichlet 判别法——它们处理形如 $\sum a_n b_n$ 的乘积型级数，其中 $a_n$ 单调变化、$b_n$ 振荡，定理给出的是**逐点收敛**结论。

进入第十章后，我们面对的是函数项级数 $\sum a_n(x) b_n(x)$。一个自然的问题是：这些收敛性结论能否"升级"为**一致收敛**？如果能，需要将条件加强到什么程度？

具体来说，数项级数版本的 Dirichlet 判别法（定理 9.14）的条件是：
1. $a_n$ 单调递减且 $a_n \to 0$；
2. $B_n = \sum_{k=1}^n b_k$ 有界。
结论：$\sum a_n b_n$ **收敛**。

要将其推广到函数项级数的一致收敛场景，条件需要做如下升级：
- $a_n \to 0$（每个 $x$ 处逐点趋于零）$\rightarrow$ $a_n(x) \rightrightarrows 0$（**一致**趋于零）
- $B_n$ 有界（每个 $x$ 处有界，界可能依赖于 $x$）$\rightarrow$ $B_n(x)$ **一致**有界（存在一个与 $x$ 无关的公共界）
- 结论：$\sum a_n(x) b_n(x)$ **一致收敛**

### 2.2 典型应用场景

**三角级数的一致收敛性**：形如 $\sum \frac{\cos(nx)}{n}$、$\sum \frac{\sin(nx)}{\sqrt{n}}$ 的级数在 $[0, 2\pi]$ 上逐点收敛但**不**一致收敛（因为靠近 $x=0$ 和 $x=2\pi$ 时收敛速度变慢）。然而，在去掉这些"奇点"的闭子区间 $[\delta, 2\pi-\delta]$（$\delta > 0$）上，它们却一致收敛。这一结论的严格证明正是依赖 Dirichlet 判别法和三角级数的部分和一致有界性。

**Fourier 级数的内闭一致收敛**：在 Fourier 分析中，周期函数的 Fourier 级数 $\frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n \cos nx + b_n \sin nx)$ 在函数的连续点附近具有内闭一致收敛性，其证明往往需要 Abel 判别法或 Dirichlet 判别法作为工具。

### 2.3 与反常积分版本的平行关系

至此，Abel-Dirichlet 判别法已经出现在三个平行场景中：

| 场景 | 核心引理 | Dirichlet 条件 | Abel 条件 |
|------|---------|---------------|-----------|
| 反常积分（ch8-04） | 积分第二中值定理 | $f \downarrow 0$, $\int_a^A g$ 有界 | $f$ 单调有界, $\int_a^\infty g$ 收敛 |
| 数项级数（ch9-05） | Abel 变换 | $a_n \downarrow 0$, $\sum_{k=1}^n b_k$ 有界 | $a_n$ 单调有界, $\sum b_n$ 收敛 |
| 函数项级数（本文件） | Abel 变换 | $a_n(x) \rightrightarrows 0$, $B_n(x)$ 一致有界 | $a_n(x)$ 一致有界且单调, $\sum b_n(x)$ 一致收敛 |

---

## 3. 核心概念与符号约定

### 3.1 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $\sum_{n=1}^{\infty} a_n(x) b_n(x)$ | 乘积型函数项级数 | $a_n(x)$ 为单调因子，$b_n(x)$ 为振荡因子 |
| $B_n(x) = \sum_{k=1}^{n} b_k(x)$ | $\{b_n(x)\}$ 的部分和函数列 | 约定 $B_0(x) \equiv 0$ |
| $S_n(x) = \sum_{k=1}^{n} a_k(x) b_k(x)$ | $\{a_n b_n\}$ 的部分和函数列 | 用于一致收敛的 $\varepsilon$-$N$ 定义 |
| $a_n(x) \rightrightarrows 0$ | $a_n(x)$ 在 $D$ 上一致收敛于零 | $\sup_{x\in D} |a_n(x)| \to 0$ |
| $B_n(x)$ 一致有界 | 存在 $M>0$ 使 $\forall n,\forall x: |B_n(x)| \le M$ | 公共界与 $n,x$ 均无关 |
| $\delta$ | 三角级数中避开奇点的参数 | $\delta > 0$，区间 $[\delta, 2\pi-\delta]$ 远离 $0$ 和 $2\pi$ |

### 3.2 Abel 变换的回顾

本文件需要频繁使用 Abel 变换的尾部形式（引理 9.13 形式二）。为方便查阅，在此重述：

**Abel 变换（尾部形式）**：对任意两个函数列 $\{a_n(x)\}$、$\{b_n(x)\}$，记 $B_n(x) = \sum_{k=1}^{n} b_k(x)$（$B_0(x) \equiv 0$），则对 $n > m \ge 0$：

$$\sum_{k=m+1}^{n} a_k(x) b_k(x) = a_n(x) B_n(x) - a_{m+1}(x) B_m(x) - \sum_{k=m+1}^{n-1} [a_{k+1}(x) - a_k(x)] B_k(x) \tag{10.12.0}$$

等价地，若将最后一项的符号改为正号并改用 $(a_k - a_{k+1})$：

$$\sum_{k=m+1}^{n} a_k(x) b_k(x) = a_n(x) B_n(x) - a_{m+1}(x) B_m(x) + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] B_k(x) \tag{10.12.0'}$$

当 $\{a_n(x)\}$ 关于 $n$ 单调递减时，$a_k(x) - a_{k+1}(x) \ge 0$，这对后续的绝对值放缩至关重要。

### 3.3 一致 Cauchy 收敛准则（函数项级数版本）

一致收敛的 $\varepsilon$-$N$ 定义（定义 10.5）虽然直接，但在实际应用中往往不便直接验证。与数项级数类似，我们有一个更便于操作的等价形式——**一致 Cauchy 收敛准则**。

**定理 10.12（一致 Cauchy 收敛准则）**：设 $\sum_{n=1}^{\infty} u_n(x)$ 是定义在集合 $D$ 上的函数项级数。则 $\sum_{n=1}^{\infty} u_n(x)$ 在 $D$ 上一致收敛当且仅当对任意 $\varepsilon > 0$，存在正整数 $N = N(\varepsilon)$（与 $x$ 无关），使得对一切 $n > m \ge N$ 和一切 $x \in D$，有

$$\left| \sum_{k=m+1}^{n} u_k(x) \right| < \varepsilon \tag{10.12.1}$$

**证明**：

（$\Rightarrow$）设 $\sum_{n=1}^{\infty} u_n(x)$ 在 $D$ 上一致收敛于 $S(x)$。由定义 10.5，对任意 $\varepsilon > 0$，存在 $N$，使得对一切 $n > N$ 和一切 $x \in D$，有 $|S(x) - S_n(x)| < \varepsilon/2$，其中 $S_n(x) = \sum_{k=1}^{n} u_k(x)$。

对任意 $n > m \ge N$ 和任意 $x \in D$：

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} u_k(x) \right| &= |S_n(x) - S_m(x)| \\
&\le |S_n(x) - S(x)| + |S(x) - S_m(x)| \\
&< \frac{\varepsilon}{2} + \frac{\varepsilon}{2} = \varepsilon
\end{aligned}
$$

因此一致收敛必定推出一致 Cauchy 条件。这里用到的正是数项级数 Cauchy 准则证明中的三角不等式技巧——将中间差拆为两个尾部差之和对消和函数。

（$\Leftarrow$）设一致 Cauchy 条件成立。对每个固定的 $x \in D$，条件 (10.12.1) 退化为数项级数的 Cauchy 收敛准则（定理 9.10）：对任意 $\varepsilon > 0$，存在 $N$（依赖于 $\varepsilon$ 且可能依赖 $x$），使得 $n > m \ge N$ 时 $|\sum_{k=m+1}^{n} u_k(x)| < \varepsilon$。由 $\mathbb{R}$ 的完备性，对每个 $x \in D$，数项级数 $\sum u_n(x)$ 收敛。定义其和函数为 $S(x)$。

现在证明收敛是一致的。由条件，对任意 $\varepsilon > 0$，存在 $N$（与 $x$ 无关），使得对一切 $n > m \ge N$ 和一切 $x \in D$，有

$$\left| \sum_{k=m+1}^{n} u_k(x) \right| < \frac{\varepsilon}{2}$$

固定 $m \ge N$。对每个 $x \in D$，由于 $\sum u_n(x)$ 收敛，令 $n \to \infty$，则 $\sum_{k=m+1}^{n} u_k(x) \to \sum_{k=m+1}^{\infty} u_k(x) = S(x) - S_m(x)$。极限保持不等式（定理 2.5，数列极限的保序性）给出：

$$|S(x) - S_m(x)| = \left| \sum_{k=m+1}^{\infty} u_k(x) \right| \le \frac{\varepsilon}{2} < \varepsilon, \quad \forall x \in D$$

因此对任意 $\varepsilon > 0$，存在 $N$（与 $x$ 无关），使得对一切 $m \ge N$ 和一切 $x \in D$，有 $|S(x) - S_m(x)| < \varepsilon$。由定义 10.5，$\sum u_n(x)$ 在 $D$ 上一致收敛。证毕。

**与数项级数 Cauchy 准则的区别**：

| 对比 | 数项级数（定理 9.10） | 函数项级数（定理 10.12） |
|------|----------------------|------------------------|
| Cauchy 条件 | $\forall\varepsilon>0,\exists N,\forall n>m\ge N$ | $\forall\varepsilon>0,\exists N,\forall n>m\ge N,\mathbf{\forall x\in D}$ |
| 结论 | 收敛 | **一致**收敛 |
| 关键差异 | $N$ 仅与 $\varepsilon$ 有关 | $N$ 与 $\varepsilon$ 和 $x$ **都**无关——需对 $x$ 一致 |

**定理 10.12 在本文中的角色**：一致 Cauchy 收敛准则是 Abel 变换尾部估计与一致收敛结论之间的桥梁。在 Dirichlet 和 Abel 判别法的证明中，我们将通过 Abel 变换将 $\sum a_k b_k$ 的尾部估计转化为 $2M a_{m+1}(x)$ 或含余项 $R_k(x)$ 的表达式，然后利用一致 Cauchy 准则将尾部估计升级为一致收敛结论。

---

## 4. Dirichlet 判别法（函数项级数版本）

### 4.1 定理陈述

**定理 10.13（Dirichlet 判别法——函数项级数版本）**：设 $\{a_n(x)\}_{n=1}^{\infty}$ 和 $\{b_n(x)\}_{n=1}^{\infty}$ 是定义在集合 $D$ 上的两个函数列。若满足：

1. **（逐点单调）** 对每个固定的 $x \in D$，数列 $\{a_n(x)\}$ 关于 $n$ **单调**；
2. **（一致趋于零）** $a_n(x) \rightrightarrows 0$ 在 $D$ 上，即 $\displaystyle\lim_{n\to\infty} \sup_{x\in D} |a_n(x)| = 0$；
3. **（部分和一致有界）** $\{B_n(x)\}$ 在 $D$ 上**一致有界**，即存在常数 $M > 0$，使得对一切 $n \ge 1$ 和一切 $x \in D$，有

$$\left| B_n(x) \right| = \left| \sum_{k=1}^{n} b_k(x) \right| \le M$$

则函数项级数 $\displaystyle\sum_{n=1}^{\infty} a_n(x) b_n(x)$ 在 $D$ 上**一致收敛**。

**注**：
- 当 $a_n(x)$ 和 $b_n(x)$ 均与 $x$ 无关时，定理 10.13 退化为 Ch9-05 的定理 9.14（数项级数版本的 Dirichlet 判别法）。
- 条件 1 中的"单调"不要求方向一致——对不同 $x$，$\{a_n(x)\}$ 可以有的递增有的递减，只需对每个固定的 $x$ 单调即可。
- 条件 2 中 $a_n(x) \rightrightarrows 0$ 等价于 $\sup_{x\in D} |a_n(x)| \to 0$（定理 10.1 应用于函数列而非级数）。当 $a_n(x)$ 与 $x$ 无关时，它退化为 $\lim a_n = 0$。

### 4.2 证明

**证明框架**：Abel 变换 + 一致有界性 + telescoping sum + 一致 Cauchy 准则。

设条件 1-3 成立。由条件 2（$a_n \rightrightarrows 0$），存在 $N_0$ 使得 $n \ge N_0$ 时 $a_n(x) \ge 0$ 对所有 $x$ 成立（注意：单调且一致趋于零的序列最终非负——因为若有 $a_n(x_0) < 0$，由单调性后续项只会更负，不可能趋于零）。为简化推导，不妨设所有 $a_n(x) \ge 0$（因为级数的前有限项不影响收敛性，可从足够大的 $n_0$ 开始考虑）。

不失一般性，设 $\{a_n(x)\}$ 对每个 $x$ **单调递减**。若递增，则考虑 $-a_n(x)$（结论不变，因为 $\sum a_n b_n$ 与 $\sum (-a_n) b_n$ 的一致收敛性等价）。此时，对每个 $x$ 有 $a_{n+1}(x) \le a_n(x)$，从而 $a_k(x) - a_{k+1}(x) \ge 0$ 对所有 $k$ 成立。

现在，对任意 $n > m \ge 1$，使用 Abel 变换的尾部形式 (10.12.0')：

$$
\begin{aligned}
\sum_{k=m+1}^{n} a_k(x) b_k(x) &= a_n(x) B_n(x) - a_{m+1}(x) B_m(x) + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] B_k(x)
\end{aligned}
$$

**第 1 步：取绝对值并用三角不等式**。

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k(x) b_k(x) \right|
&\le a_n(x) |B_n(x)| + a_{m+1}(x) |B_m(x)| + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] |B_k(x)|
\end{aligned}
$$

其中 $a_n(x) \ge 0$ 且 $a_k(x) - a_{k+1}(x) \ge 0$，故绝对值可直接去掉。

**第 2 步：利用一致有界性 $|B_k(x)| \le M$ 放缩**。

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k(x) b_k(x) \right|
&\le M \left( a_n(x) + a_{m+1}(x) + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] \right)
\end{aligned}
$$

**第 3 步：计算 telescoping sum**。

$$
\sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] = a_{m+1}(x) - a_n(x)
$$

这是因为括号中的项首尾相消：
$$(a_{m+1} - a_{m+2}) + (a_{m+2} - a_{m+3}) + \cdots + (a_{n-1} - a_n) = a_{m+1} - a_n.$$

**第 4 步：代入得到简洁上界**。

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k(x) b_k(x) \right|
&\le M \big( a_n(x) + a_{m+1}(x) + [a_{m+1}(x) - a_n(x)] \big) \\
&= 2M \cdot a_{m+1}(x)
\end{aligned}
\tag{10.13.1}
$$

这一估计极其关键：它将乘积型级数的尾部完全控制为 $2M$ 乘以 $a_{m+1}(x)$——$M$ 是常数，$a_{m+1}(x)$ 是单调递减因子中"最靠前"的那一项。

**第 5 步：利用 $a_n \rightrightarrows 0$ 得到一致 Cauchy 条件**。

由条件 2（$a_n(x) \rightrightarrows 0$），对任意 $\varepsilon > 0$，存在 $N$，使得对所有 $n \ge N$ 和所有 $x \in D$，有

$$a_n(x) < \frac{\varepsilon}{2M}$$

特别地，对 $m+1 \ge N$，有 $a_{m+1}(x) < \varepsilon/(2M)$ 对所有 $x \in D$ 成立。

于是对任意 $n > m \ge N$ 和任意 $x \in D$，由 (10.13.1)：

$$\left| \sum_{k=m+1}^{n} a_k(x) b_k(x) \right| \le 2M \cdot a_{m+1}(x) < 2M \cdot \frac{\varepsilon}{2M} = \varepsilon$$

由一致 Cauchy 收敛准则（定理 10.12），$\sum_{n=1}^{\infty} a_n(x) b_n(x)$ 在 $D$ 上一致收敛。证毕。

**证明逻辑链**：

```
条件 3: B_n(x) 一致有界 (|B_n(x)| ≤ M)
  + 条件 1: a_n(x) 单调 ⇒ a_k - a_{k+1} ≥ 0
  + Abel 变换尾部形式
  ↓
|∑_{k=m+1}^n a_k b_k| ≤ M(a_n + a_{m+1} + (a_{m+1} - a_n)) = 2M a_{m+1}(x)
  + 条件 2: a_n(x) ⇉ 0 ⇒ a_{m+1}(x) < ε/(2M) 对充分大的 m
  ↓
一致 Cauchy 准则 ⇒ ∑ a_n b_n 一致收敛
```

### 4.3 与数项级数版本（定理 9.14）的类比与区别

| 对比项 | 数项级数版本（定理 9.14） | 函数项级数版本（定理 10.13） |
|--------|-------------------------|---------------------------|
| $a_n$ 条件 | $a_n \downarrow 0$ | $\forall x$: $a_n(x)$ 单调，且 $a_n(x) \rightrightarrows 0$ |
| $b_n$ 条件 | $B_n = \sum_{k=1}^n b_k$ 有界 | $B_n(x)$ **一致**有界 |
| 证明工具 | Abel 变换 + 数列 Cauchy 准则 | Abel 变换 + 一致 Cauchy 准则（定理 10.12） |
| 尾部估计 | $|\sum_{k=m+1}^n a_k b_k| \le 2M a_{m+1}$ | **完全相同**的上界 $2M a_{m+1}(x)$，但 $M$ 和 $a_{m+1}(x)$ 中的 $x$ 被 $\sup$ 处理 |
| 结论 | $\sum a_n b_n$ **收敛** | $\sum a_n(x) b_n(x)$ **一致收敛** |

核心升级在于：
1. **条件升级**：单调性从数列性质变为"对每个 $x$ 逐点单调"；$a_n \to 0$ 升级为 $a_n \rightrightarrows 0$；$B_n$ 有界升级为 $B_n(x)$ 一致有界。
2. **结论升级**：从"逐点收敛"升级为"一致收敛"。
3. **证明升级**：在尾部估计 $|\sum a_k b_k| \le 2M a_{m+1}$ 之后，数项级数版本直接用 $a_N < \varepsilon/(2M)$ 和柯西准则得收敛；函数项级数版本需用一致 Cauchy 准则——$a_{m+1}(x)$ 的上界必须对 $x$ 一致，这正是 $a_n \rightrightarrows 0$ 所保证的。

### 4.4 条件验证的标准操作流程

应用 Dirichlet 判别法判断 $\sum a_n(x) b_n(x)$ 的一致收敛性时，遵循以下四步流程：

**第 1 步：拆分**。将级数的通项写为 $a_n(x) \cdot b_n(x)$ 的形式，其中 $a_n(x)$ 单调、$b_n(x)$ 振荡。

**第 2 步：验证单调性与一致趋于零**。
- 验证对每个 $x$，$\{a_n(x)\}$ 关于 $n$ 单调（可通过 $a_{n+1}(x) - a_n(x)$ 的符号判定）。
- 验证 $\sup_{x\in D} |a_n(x)| \to 0$。如果 $a_n$ 与 $x$ 无关，只需验证 $\lim a_n = 0$。

**第 3 步：验证部分和一致有界**。
- 证明 $|B_n(x)| = |\sum_{k=1}^n b_k(x)| \le M$ 对一切 $n$ 和一切 $x$ 成立。
- 常用工具：三角恒等式给出 $\sum \cos(kx)$ 和 $\sum \sin(kx)$ 的闭式，然后通过限制区间使 $1/|\sin(x/2)|$ 的上界与 $x$ 无关。

**第 4 步：应用定理**。三个条件均满足后，由定理 10.13 得一致收敛。

---

## 5. 三角级数部分和的一致有界性

在应用 Dirichlet 判别法于 $\sum \cos(nx)/n$ 或 $\sum \sin(nx)/\sqrt{n}$ 等三角级数时，最关键的条件是 $B_n(x) = \sum_{k=1}^n \cos(kx)$（或 $\sum_{k=1}^n \sin(kx)$）的一致有界性。

### 5.1 闭式公式回顾

由 Ch9-05 的推导（公式 9.14.1/9.14.2），对任意 $\theta \neq 2m\pi$（$m \in \mathbb{Z}$）：

$$\sum_{k=1}^{n} \sin(k\theta) = \frac{\cos\frac{\theta}{2} - \cos\left((n+\frac12)\theta\right)}{2\sin\frac{\theta}{2}} \tag{10.13.2}$$

$$\sum_{k=1}^{n} \cos(k\theta) = \frac{\sin\left((n+\frac12)\theta\right) - \sin\frac{\theta}{2}}{2\sin\frac{\theta}{2}} \tag{10.13.3}$$

由此得逐点有界性上界：

$$\left|\sum_{k=1}^{n} \sin(k\theta)\right| \le \frac{1}{|\sin\frac{\theta}{2}|}, \quad \left|\sum_{k=1}^{n} \cos(k\theta)\right| \le \frac{1}{|\sin\frac{\theta}{2}|} \tag{10.13.4}$$

### 5.2 从逐点有界到一致有界

不等式 (10.13.4) 的右端依赖于 $\theta$——当 $\theta$ 接近 $0$ 或 $2\pi$ 时，$|\sin(\theta/2)|$ 接近 $0$，上界 $1/|\sin(\theta/2)|$ 趋于无穷。因此，在包含 $0$ 或 $2\pi$ 的区间上，$\{B_n(x)\}$ **不是**一致有界的。

然而，若将 $x$ 限制在**远离**奇点 $x = 2m\pi$ 的闭区间上，则可以获得与 $x$ 无关的一致上界。

**引理 10.1（三角级数部分和的一致有界性）**：对任意 $\delta > 0$，在区间 $[\delta, 2\pi-\delta]$ 上，函数列

$$B_n^{\cos}(x) = \sum_{k=1}^{n} \cos(kx), \quad B_n^{\sin}(x) = \sum_{k=1}^{n} \sin(kx)$$

**一致有界**——即存在常数 $M_\delta = \dfrac{1}{\sin(\delta/2)}$（仅依赖于 $\delta$，与 $n$ 和 $x$ 无关），使得

$$|B_n^{\cos}(x)| \le M_\delta, \quad |B_n^{\sin}(x)| \le M_\delta, \quad \forall n \ge 1, \forall x \in [\delta, 2\pi-\delta]$$

**证明**：由 (10.13.4)，对任意 $x \in [\delta, 2\pi-\delta]$（$x \neq 0, 2\pi$，故分母非零）：

$$|B_n^{\cos}(x)| \le \frac{1}{|\sin(x/2)|}, \quad |B_n^{\sin}(x)| \le \frac{1}{|\sin(x/2)|}$$

现考察函数 $f(x) = |\sin(x/2)|$ 在 $[\delta, 2\pi-\delta]$ 上的最小值。由于 $x/2 \in [\delta/2, \pi - \delta/2]$，而 $\sin(t)$ 在 $[0, \pi]$ 上非负，在 $[\delta/2, \pi-\delta/2]$ 上的最小值为 $\min\{\sin(\delta/2), \sin(\pi-\delta/2)\} = \sin(\delta/2)$（因为 $\sin(\pi-\delta/2) = \sin(\delta/2)$）。因此

$$\sin(x/2) \ge \sin(\delta/2) > 0, \quad \forall x \in [\delta, 2\pi-\delta]$$

从而

$$|B_n^{\cos}(x)| \le \frac{1}{\sin(x/2)} \le \frac{1}{\sin(\delta/2)}, \quad \forall x \in [\delta, 2\pi-\delta]$$

对 $B_n^{\sin}$ 同理。证毕。

**区间选择的几何意义**：上式中 $\delta$ 是"安全距离"——$x$ 离奇点 $0$ 越远（$\delta$ 越大），$\sin(\delta/2)$ 越大，一致上界 $1/\sin(\delta/2)$ 越小。当 $\delta \to 0^+$ 时，$\sin(\delta/2) \sim \delta/2 \to 0$，一致上界 $1/\sin(\delta/2) \to \infty$，引理失效——这正是 $\sum \cos(nx)/n$ 在 $(0, 2\pi)$ 上不一致收敛的原因。

**与 Ch9-05 逐点有界的区别**：

| 对比 | Ch9-05 逐点有界 | 本文件一致有界 |
|------|-----------------|---------------|
| 上界 | $1/|\sin(x/2)|$，依赖于 $x$ | $1/\sin(\delta/2)$，与 $x$ 无关 |
| 适用范围 | 任意 $x \neq 2m\pi$ | $x \in [\delta, 2\pi-\delta]$ |
| 含义 | 每个 $x$ 处有界，但界可随 $x$ 变 | 所有 $x$ 共享同一个界 |
| 证明关键 | 三角恒等式 | 三角恒等式 + 区间下界估计 |

---

## 6. Abel 判别法（函数项级数版本）

### 6.1 定理陈述

**定理 10.14（Abel 判别法——函数项级数版本）**：设 $\{a_n(x)\}_{n=1}^{\infty}$ 和 $\{b_n(x)\}_{n=1}^{\infty}$ 是定义在集合 $D$ 上的两个函数列。若满足：

1. **（逐点单调）** 对每个固定的 $x \in D$，数列 $\{a_n(x)\}$ 关于 $n$ **单调**；
2. **（一致有界）** $\{a_n(x)\}$ 在 $D$ 上**一致有界**，即存在常数 $K > 0$，使得 $|a_n(x)| \le K$ 对一切 $n$ 和一切 $x \in D$ 成立；
3. **（一致收敛）** $\displaystyle\sum_{n=1}^{\infty} b_n(x)$ 在 $D$ 上**一致收敛**，

则函数项级数 $\displaystyle\sum_{n=1}^{\infty} a_n(x) b_n(x)$ 在 $D$ 上**一致收敛**。

### 6.2 与数项级数版本的类比

| 对比项 | 数项级数版本（定理 9.15） | 函数项级数版本（定理 10.14） |
|--------|-------------------------|---------------------------|
| $a_n$ 条件 | 单调有界 | $\forall x$: $a_n(x)$ 单调，且**一致有界** |
| $b_n$ 条件 | $\sum b_n$ 收敛 | $\sum b_n(x)$ **一致收敛** |
| 结论 | $\sum a_n b_n$ 收敛 | $\sum a_n(x) b_n(x)$ **一致收敛** |

### 6.3 归化证明（通过 Dirichlet 判别法）

归化法的核心思想是：将 $a_n(x)$ 减去其极限（或极限函数），使得新序列单调趋于零，从而满足 Dirichlet 判别法的条件。但这里有一个微妙之处：$a_n(x)$ 的一致有界和逐点单调性只保证对每个 $x$，$\{a_n(x)\}$ 收敛到某个极限 $L(x)$，但 $L(x)$ 本身可能依赖于 $x$，且 $a_n(x) - L(x) \rightrightarrows 0$ 不一定成立。

**归化法的适用条件**：当 $L(x) \equiv L$ 为常数（与 $x$ 无关）时，归化法成立。这包括以下重要情形：
- $a_n(x)$ 与 $x$ 无关，即 $a_n(x) = a_n$，此时 $L(x) = \lim a_n$ 为常数；
- $a_n(x)$ 一致收敛于一个常数函数。

对于更一般的情形（$a_n(x)$ 单调一致有界但逐点收敛于非常数函数 $L(x)$），需要用直接 Abel 变换法。

**归化证明（对 $L(x) \equiv L$ 常数情形）**：

由条件 1，对每个 $x$，$\{a_n(x)\}$ 单调。由条件 2，$\{a_n(x)\}$ 一致有界。因此对每个 $x$，$\{a_n(x)\}$ 是单调有界数列，由单调有界定理（第二章定理 10.2）知极限存在，记 $L = \lim_{n\to\infty} a_n(x)$（常数，与 $x$ 无关）。

定义新函数列：

$$a_n'(x) = a_n(x) - L, \quad n = 1, 2, \dots$$

则：
- $\{a_n'(x)\}$ 与 $\{a_n(x)\}$ 具有相同的单调性（减去常数不改变单调性）；
- $a_n'(x) \rightrightarrows 0$：因为 $a_n(x)$ 一致有界，且对每个 $x$，$\lim a_n'(x) = 0$。由于 $a_n$ 与 $x$ 无关时显然一致趋于零；当 $a_n$ 与 $x$ 有关但极限为常数时，需额外验证——但在归化法典型应用中，$a_n(x)$ 与 $x$ 无关或极限为常数，此时 $a_n'(x) \to 0$ 是通常的数列极限，自然是一致收敛的（因为与 $x$ 无关）。

由条件 3，$\sum b_n(x)$ 在 $D$ 上一致收敛。一致收敛级数的部分和函数列 $\{B_n(x)\}$ 必然一致有界（因为一致收敛数列必有界，定理 10.1 的推论）。

现在对 $\{a_n'(x)\}$ 和 $\{b_n(x)\}$ 应用 Dirichlet 判别法（定理 10.13）：
- $\{a_n'(x)\}$ 对每个 $x$ 单调，且 $a_n'(x) \rightrightarrows 0$；
- $B_n(x) = \sum_{k=1}^n b_k(x)$ 一致有界（由一致收敛性推出）。

故 $\sum a_n'(x) b_n(x)$ 在 $D$ 上一致收敛。

最后：

$$\sum_{n=1}^{\infty} a_n(x) b_n(x) = \sum_{n=1}^{\infty} [a_n'(x) + L] b_n(x) = \sum_{n=1}^{\infty} a_n'(x) b_n(x) + L \sum_{n=1}^{\infty} b_n(x)$$

右边第一项由 Dirichlet 判别法一致收敛，第二项 $L\sum b_n(x)$ 是一致收敛级数的常数倍，也一致收敛。两个一致收敛的级数之和仍一致收敛，故 $\sum a_n(x) b_n(x)$ 在 $D$ 上一致收敛。证毕。

### 6.4 直接 Abel 变换证明（一般情形）

当 $a_n(x)$ 的极限函数 $L(x)$ 与 $x$ 有关且不保证一致收敛时，归化法失效，需要采用直接 Abel 变换证明。

**直接证明**：

由条件 3，$\sum_{n=1}^{\infty} b_n(x)$ 在 $D$ 上一致收敛于和函数 $S_b(x)$。记余项

$$R_n(x) = \sum_{k=n+1}^{\infty} b_k(x) = S_b(x) - B_n(x)$$

由一致收敛性，$R_n(x) \rightrightarrows 0$（即 $\sup_{x\in D} |R_n(x)| \to 0$）。

设 $\{a_n(x)\}$ 单调递减（若递增，考虑 $-a_n$，类似可证）。由条件 2，$|a_n(x)| \le K$ 对所有 $n$ 和 $x$ 成立。由条件 1，$a_{n+1}(x) \le a_n(x)$，故 $a_n(x) - a_{n+1}(x) \ge 0$。

使用 Abel 变换的尾部形式 (10.12.0')，对 $n > m$：

$$\sum_{k=m+1}^{n} a_k(x) b_k(x) = a_n(x) B_n(x) - a_{m+1}(x) B_m(x) + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] B_k(x)$$

将 $B_k(x) = S_b(x) - R_k(x)$ 代入，得：

$$
\begin{aligned}
\sum_{k=m+1}^{n} a_k b_k &= a_n(S_b - R_n) - a_{m+1}(S_b - R_m) + \sum_{k=m+1}^{n-1} (a_k - a_{k+1})(S_b - R_k) \\
&= S_b \left( a_n - a_{m+1} - \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) \right) - a_n R_n + a_{m+1} R_m + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) R_k
\end{aligned}
$$

注意到 $S_b$ 的系数是一个 telescoping sum：

$$a_n - a_{m+1} - \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) = a_n - a_{m+1} - (a_{m+1} - a_n) = 0$$

因此 $S_b$ 项全部消去，只剩下含余项 $R_k$ 的项：

$$\sum_{k=m+1}^{n} a_k(x) b_k(x) = -a_n(x) R_n(x) + a_{m+1}(x) R_m(x) + \sum_{k=m+1}^{n-1} [a_k(x) - a_{k+1}(x)] R_k(x) \tag{10.14.1}$$

取绝对值并用三角不等式：

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k b_k \right|
&\le |a_n| \cdot |R_n| + |a_{m+1}| \cdot |R_m| + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) |R_k|
\end{aligned}
$$

其中用到 $a_k - a_{k+1} \ge 0$。

由一致有界性 $|a_n(x)| \le K$ 得：

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k b_k \right|
&\le K |R_n| + K |R_m| + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) |R_k|
\end{aligned}
\tag{10.14.2}
$$

现在，由于 $R_n(x) \rightrightarrows 0$，对任意 $\varepsilon > 0$，存在 $N$，使得对所有 $k \ge N$ 和所有 $x \in D$，有

$$|R_k(x)| < \frac{\varepsilon}{3K}$$

对 $n > m \ge N$，从 (10.14.2) 式估计：

$$
\begin{aligned}
\left| \sum_{k=m+1}^{n} a_k b_k \right|
&\le K \cdot \frac{\varepsilon}{3K} + K \cdot \frac{\varepsilon}{3K} + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) \cdot \frac{\varepsilon}{3K} \\
&= \frac{2\varepsilon}{3} + \frac{\varepsilon}{3K} \sum_{k=m+1}^{n-1} (a_k - a_{k+1})
\end{aligned}
$$

计算 telescoping sum：

$$\sum_{k=m+1}^{n-1} (a_k - a_{k+1}) = a_{m+1} - a_n \le a_{m+1} + (-a_n) \le K + K = 2K$$

其中最后一步利用了 $|a_n| \le K$（注意 $\{a_n\}$ 递减且一致有界时，$a_{m+1} - a_n \le 2K$，因为最大可能跨度是从上界 $K$ 到下界 $-K$）。

代入得：

$$
\left| \sum_{k=m+1}^{n} a_k b_k \right| \le \frac{2\varepsilon}{3} + \frac{\varepsilon}{3K} \cdot 2K = \frac{2\varepsilon}{3} + \frac{2\varepsilon}{3} = \frac{4\varepsilon}{3}
$$

这个上界是 $4\varepsilon/3$ 而非 $\varepsilon$——需要调整系数的选取。若取 $N$ 使得 $k \ge N$ 时 $|R_k(x)| < \varepsilon/(4K)$（而非 $\varepsilon/(3K)$），则：

$$
\left| \sum_{k=m+1}^{n} a_k b_k \right| \le \frac{2\varepsilon}{4} + \frac{\varepsilon}{4K} \cdot 2K = \frac{\varepsilon}{2} + \frac{\varepsilon}{2} = \varepsilon
$$

因此，由一致 Cauchy 收敛准则（定理 10.12），$\sum a_n(x) b_n(x)$ 在 $D$ 上一致收敛。证毕。

### 6.5 两类判别法的条件对比

| 对比维度 | Dirichlet 判别法（定理 10.13） | Abel 判别法（定理 10.14） |
|----------|-------------------------------|---------------------------|
| $a_n(x)$ 的单调性 | 对每个 $x$ 单调 | 对每个 $x$ 单调（相同） |
| $a_n(x)$ 的强度 | $a_n(x) \rightrightarrows 0$（**一致趋于零**） | $a_n(x)$ **一致有界**（更弱） |
| $b_n(x)$ 的强度 | $B_n(x) = \sum_{k=1}^n b_k$ **一致有界**（更弱） | $\sum b_n(x)$ **一致收敛**（更强） |
| 证明核心 | 直接 Abel 变换 + telescoping | 归化为 Dirichlet 或直接 Abel 变换 |
| 适用场景 | $a_n$ 衰减 + $b_n$ 振荡型 | $a_n$ 有界缓慢变化 + $b_n$ 已一致收敛 |

**对偶性**：Dirichlet 判别法对 $a_n$ 要求强（趋于零）而对 $b_n$ 要求弱（部分和有界即可）；Abel 判别法则相反——对 $a_n$ 要求弱（有界即可）而对 $b_n$ 要求强（级数一致收敛）。两者的条件呈对偶互补关系。

---

## 7. 例题

### 例题 1：Dirichlet 判别法——$\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n}$ 在 $[\delta, 2\pi-\delta]$ 上的一致收敛性

证明：对任意 $\delta > 0$，函数项级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n}$ 在 $[\delta, 2\pi-\delta]$ 上一致收敛。

**解**：

**第 1 步：拆分通项**。

令 $a_n(x) = \dfrac{1}{n}$，$b_n(x) = \cos(nx)$。则原级数写作 $\sum_{n=1}^{\infty} a_n(x) b_n(x)$。

**第 2 步：验证条件 1——$a_n(x)$ 单调（对每个 $x$）**。

$a_n(x) = 1/n$ 与 $x$ 无关，且 $1/(n+1) < 1/n$，故 $\{a_n(x)\}$ 对每个 $x$ 严格单调递减。

**第 3 步：验证条件 2——$a_n(x) \rightrightarrows 0$**。

$a_n(x) = 1/n$ 与 $x$ 无关，因此

$$\sup_{x\in [\delta, 2\pi-\delta]} |a_n(x)| = \frac{1}{n} \to 0 \quad (n \to \infty)$$

由 $\sup$ 等价刻画，$a_n(x) \rightrightarrows 0$ 在 $[\delta, 2\pi-\delta]$ 上成立。

**第 4 步：验证条件 3——$B_n(x) = \sum_{k=1}^n \cos(kx)$ 一致有界**。

由 (10.13.4)，对任意 $x \neq 2m\pi$：

$$\left| \sum_{k=1}^{n} \cos(kx) \right| \le \frac{1}{|\sin(x/2)|}$$

在 $[\delta, 2\pi-\delta]$ 上（其中 $\delta > 0$），$\sin(x/2) \ge \sin(\delta/2) > 0$，因此

$$\left| \sum_{k=1}^{n} \cos(kx) \right| \le \frac{1}{\sin(x/2)} \le \frac{1}{\sin(\delta/2)}, \quad \forall n \ge 1, \forall x \in [\delta, 2\pi-\delta]$$

取 $M_\delta = 1/\sin(\delta/2)$，则 $|B_n(x)| \le M_\delta$ 对一切 $n$ 和一切 $x$ 成立，即 $\{B_n(x)\}$ 在 $[\delta, 2\pi-\delta]$ 上一致有界。

**第 5 步：应用 Dirichlet 判别法**。

三个条件均满足，由定理 10.13，$\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n}$ 在 $[\delta, 2\pi-\delta]$ 上**一致收敛**。

**关于 $\delta$ 的必要性**：$\delta$ 必须大于 0。若 $\delta = 0$（即区间为 $[0, 2\pi]$），当 $x \to 0^+$ 时，$|\sin(x/2)| \sim |x|/2 \to 0$，$1/|\sin(x/2)| \to \infty$，$\{B_n(x)\}$ 不再是**一致**有界的——事实上，该级数在 $[0, 2\pi]$ 上不一致收敛。这体现了"内闭一致收敛"的性质：级数在 $(0, 2\pi)$ 的任意闭子区间上一致收敛，但在整个区间 $(0, 2\pi)$ 上不一致收敛。

---

### 例题 2：Dirichlet 判别法——$\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上的一致收敛性

证明：对任意 $\delta > 0$，函数项级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上一致收敛。

**解**：

**第 1 步：拆分通项**。

令 $a_n(x) = \dfrac{1}{\sqrt{n}}$，$b_n(x) = \sin(nx)$。则原级数写作 $\sum a_n(x) b_n(x)$。

**第 2 步：验证条件 1——单调性**。

$a_n(x) = 1/\sqrt{n}$ 与 $x$ 无关，且 $\{1/\sqrt{n}\}$ 严格单调递减。

**第 3 步：验证条件 2——$a_n(x) \rightrightarrows 0$**。

$$\sup_{x\in[\delta,\pi-\delta]} |a_n(x)| = \frac{1}{\sqrt{n}} \to 0 \quad (n \to \infty)$$

故 $a_n(x) \rightrightarrows 0$ 在 $[\delta, \pi-\delta]$ 上成立。

**第 4 步：验证条件 3——$B_n(x) = \sum_{k=1}^n \sin(kx)$ 一致有界**。

由 (10.13.4)，对任意 $x \neq 2m\pi$：

$$\left| \sum_{k=1}^{n} \sin(kx) \right| \le \frac{1}{|\sin(x/2)|}$$

在 $[\delta, \pi-\delta]$ 上，$x/2 \in [\delta/2, \pi/2 - \delta/2]$（当 $\delta \le \pi$ 时；若 $\delta > \pi$ 则区间为空，无需考虑）。$\sin(t)$ 在 $[\delta/2, \pi/2 - \delta/2] \subset (0, \pi/2)$ 上单调递增，最小值在左端点取得：$\sin(x/2) \ge \sin(\delta/2) > 0$。因此

$$\left| \sum_{k=1}^{n} \sin(kx) \right| \le \frac{1}{\sin(x/2)} \le \frac{1}{\sin(\delta/2)}, \quad \forall n \ge 1, \forall x \in [\delta, \pi-\delta]$$

即 $\{B_n(x)\}$ 在 $[\delta, \pi-\delta]$ 上一致有界。

**第 5 步：应用 Dirichlet 判别法**。

三个条件均满足，由定理 10.13，$\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上**一致收敛**。

**注**：与例题 1 相比，本题的 $a_n = 1/\sqrt{n}$ 衰减更慢（$1/\sqrt{n}$ 比 $1/n$ 衰减慢），但级数仍然一致收敛——Dirichlet 判别法对 $a_n$ 的唯一要求是单调趋于零（无论多慢），只要 $B_n$ 一致有界即可。

---

### 例题 3：Abel 判别法——$\displaystyle\sum_{n=1}^{\infty} \frac{n}{n+1} \cdot \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上的一致收敛性

证明：对任意 $\delta > 0$，函数项级数 $\displaystyle\sum_{n=1}^{\infty} \frac{n}{n+1} \cdot \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上一致收敛。

**解**：

**第 1 步：拆分通项，选择判别法**。

令 $a_n(x) = \dfrac{n}{n+1}$，$b_n(x) = \dfrac{\sin(nx)}{\sqrt{n}}$。则原级数写作 $\sum a_n(x) b_n(x)$。

观察结构：
- $a_n(x) = n/(n+1)$ 递增但趋于 $1$（不趋于 0！），因此**不满足 Dirichlet 判别法**的条件 2。
- 但 $\{a_n\}$ 单调有界，且由例题 2 知 $\sum \sin(nx)/\sqrt{n}$ 在 $[\delta, \pi-\delta]$ 上一致收敛——这恰好满足 Abel 判别法的结构。

因此选用 **Abel 判别法（定理 10.14）**。

**第 2 步：验证条件 1——$a_n(x)$ 单调**。

$a_n(x) = n/(n+1)$ 与 $x$ 无关。$n/(n+1) = 1 - 1/(n+1)$，显然 $\{1 - 1/(n+1)\}$ 严格单调递增。

**第 3 步：验证条件 2——$a_n(x)$ 一致有界**。

对一切 $n \ge 1$，$0 < \dfrac{n}{n+1} < 1$，故 $|a_n(x)| \le 1$ 对所有 $n$ 和 $x$ 成立。

**第 4 步：验证条件 3——$\sum b_n(x)$ 一致收敛**。

$b_n(x) = \sin(nx)/\sqrt{n}$。由例题 2，$\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上一致收敛。这正是 $\sum b_n(x)$。

**第 5 步：应用 Abel 判别法**。

三个条件均满足，由定理 10.14，$\displaystyle\sum_{n=1}^{\infty} \frac{n}{n+1} \cdot \frac{\sin(nx)}{\sqrt{n}}$ 在 $[\delta, \pi-\delta]$ 上**一致收敛**。

**对比**：若不用 Abel 判别法，能否通过其他途径？尝试直接放缩：

$$\left| \frac{n}{n+1} \cdot \frac{\sin(nx)}{\sqrt{n}} \right| \le \frac{1}{\sqrt{n}}$$

但 $\sum 1/\sqrt{n}$ 发散，M-判别法不适用。Abel 判别法的优势正在于此——它不要求 $\sum |a_n b_n|$ 收敛（绝对收敛），而是利用 $a_n$ 的单调有界性和 $\sum b_n$ 的一致收敛性，通过正负抵消机制得到一致收敛。

---

## 8. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 函数项级数的 Dirichlet 判别法只是把 Ch9-05 定理 9.14 的 $a_n$ 换成 $a_n(x)$，结论自然成立 | 条件必须同时升级：(1) $a_n\to 0$ 升级为 $a_n(x)\rightrightarrows 0$（一致性），(2) $B_n$ 有界升级为 $B_n(x)$ **一致有界**（界与 $x$ 无关）。缺少任一个升级都推不出一致收敛 |
| $|\sum_{k=1}^n \cos(kx)| \le 1/|\sin(x/2)|$ 给出了与 $n$ 无关的上界，说明 $\{B_n\}$ 一致有界 | 这个上界依赖于 $x$——当 $x$ 接近 $0$ 时它趋于无穷。要用它得到一致有界，必须限制 $x$ 远离 $0$（例如 $x \in [\delta, 2\pi-\delta]$），使得 $1/|\sin(x/2)|$ 有与 $x$ 无关的有限上界 |
| Abel 判别法的归化证明中，$a_n(x)$ 单调一致有界 $\Rightarrow$ $a_n(x) \rightrightarrows L(x)$，然后 $a_n(x)-L(x) \rightrightarrows 0$ 自动成立 | 单调有界数列逐点收敛但不一定一致收敛。归化法仅当 $\lim a_n(x) = L$ 为常数（与 $x$ 无关）时适用。一般情形需用直接 Abel 变换证明（$S_b$ 项消去 + 余项 $R_k$ 放缩）|
| 因为 $\sum \sin(nx)/\sqrt{n}$ 在 $[\delta, \pi-\delta]$ 上一致收敛，所以 $\sum n/(n+1) \cdot \sin(nx)/\sqrt{n}$ 也自动一致收敛（不需要任何条件） | $\sum a_n b_n$ 的一致收敛性不是 $\sum b_n$ 一致收敛的直接推论——乘 $a_n$ 后可能破坏一致收敛性。Abel 判别法给出了充分条件：$a_n$ 单调有界即可保证一致收敛性能从 $\sum b_n$ 传递到 $\sum a_n b_n$ |
| 一致 Cauchy 准则与数项级数的 Cauchy 准则完全相同 | 函数项级数版本多了一个 $\forall x \in D$ 的量词，且 $N$ 不能依赖 $x$。这是本质区别——数项级数版本是对一个数列使用 Cauchy 准则，而函数项级数版本要求所有 $x$ 处的 Cauchy 条件被同一个 $N$ 同时满足 |
| Dirichlet 判别法的证明中 $\sum_{k=m+1}^{n-1} (a_k - a_{k+1}) = a_{m+1} - a_n$ 对任意 $\{a_n\}$ 都成立 | 这要求 $a_k - a_{k+1}$ 是定义良好的，且求和是有限 telescoping sum。但更关键的是：**在绝对值放缩步骤中**，我们需要 $(a_k - a_{k+1})$ 非负才能将绝对值直接去掉——这正是 $\{a_n\}$ 单调递减的作用。若 $\{a_n\}$ 不单调，则 $(a_k - a_{k+1})$ 可正可负，无法得到 $2M a_{m+1}$ 的简洁上界 |

### 检查点

- [ ] 能否写出一致 Cauchy 收敛准则（定理 10.12）的完整陈述，并指出其与数项级数 Cauchy 准则（定理 9.10）的区别？
- [ ] 能否独立完成定理 10.12 的（$\Leftarrow$）方向证明（从 Cauchy 条件到一致收敛）？
- [ ] 能否写出 Dirichlet 判别法（定理 10.13）的三个条件并完成完整证明？
- [ ] 能否说出 Dirichlet 判别法证明中，三个条件分别在哪个步骤被用到？
- [ ] 能否推导 $\sum_{k=1}^n \sin(kx)$ 和 $\sum_{k=1}^n \cos(kx)$ 的闭式公式，并由此证明在 $[\delta, 2\pi-\delta]$ 上的一致有界性？
- [ ] 能否解释为什么 $\delta > 0$ 是必需的？如果 $\delta = 0$ 会如何？
- [ ] 能否写出 Abel 判别法（定理 10.14）的三个条件并完成直接 Abel 变换证明？
- [ ] 能否说出归化法证明 Abel 判别法的适用条件（何时 $a_n(x) - L(x) \rightrightarrows 0$ 自动成立）？
- [ ] 能否区分何时用 Dirichlet 判别法、何时用 Abel 判别法？
- [ ] 能否完成例题 1-3 中至少一个的完整推导，包括所有条件验证？
- [ ] 能否将本文件的内容与 Ch9-05（数项级数版本）形成平行对照，理解条件升级的对应关系？

---

## 练习题

### 基础巩固

**1.** 判断以下函数项级数在指定区间上是否一致收敛。若收敛，说明使用的判别法并逐条验证条件；若不收敛，说明理由。

(1) $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{n^{3/2}}$，$x \in \mathbb{R}$

(2) $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n^2}$，$x \in [\delta, 2\pi-\delta]$（$\delta > 0$）

(3) $\displaystyle\sum_{n=1}^{\infty} \frac{n}{2n+1} \cdot \frac{\cos(nx)}{n}$，$x \in [\delta, 2\pi-\delta]$（$\delta > 0$）

<details><summary>参考答案</summary>

**(1)** $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(nx)}{n^{3/2}}$，$x \in \mathbb{R}$

此题可直接用 Weierstrass M-判别法（定理 10.2）：

$$\left|\frac{\sin(nx)}{n^{3/2}}\right| \le \frac{1}{n^{3/2}}$$

$\sum_{n=1}^{\infty} 1/n^{3/2}$ 是 $p = 3/2 > 1$ 的 $p$-级数，收敛。由 M-判别法，该级数在 $\mathbb{R}$ 上一致收敛。

（不必用 Dirichlet 判别法——M-判别法更直接。）

---

**(2)** $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n^2}$，$x \in [\delta, 2\pi-\delta]$（$\delta > 0$）

同样可用 M-判别法：$|\cos(nx)/n^2| \le 1/n^2$，$\sum 1/n^2$ 收敛，故级数在 $\mathbb{R}$ 上一致收敛（不仅是在 $[\delta, 2\pi-\delta]$ 上）。

**注**：M-判别法给出的控制级数与 $x$ 无关，因此当 $|u_n(x)|$ 可以被一个与 $x$ 无关的收敛级数控制时，优先使用 M-判别法。Dirichlet 判别法通常用在 M-判别法失效的场景（如 $a_n = 1/n$ 时 $\sum 1/n$ 发散）。

---

**(3)** $\displaystyle\sum_{n=1}^{\infty} \frac{n}{2n+1} \cdot \frac{\cos(nx)}{n}$，$x \in [\delta, 2\pi-\delta]$（$\delta > 0$）

先化简通项：$\dfrac{n}{2n+1} \cdot \dfrac{\cos(nx)}{n} = \dfrac{\cos(nx)}{2n+1}$。

令 $a_n(x) = \dfrac{1}{2n+1}$，$b_n(x) = \cos(nx)$。使用 Dirichlet 判别法（定理 10.13）。

**条件 1**：$a_n = 1/(2n+1)$ 与 $x$ 无关，严格单调递减。$\checkmark$

**条件 2**：$\sup_x |a_n| = 1/(2n+1) \to 0$，故 $a_n \rightrightarrows 0$。$\checkmark$

**条件 3**：由引理 10.1，在 $[\delta, 2\pi-\delta]$ 上，$|\sum_{k=1}^n \cos(kx)| \le 1/\sin(\delta/2)$，一致有界。$\checkmark$

由 Dirichlet 判别法，该级数在 $[\delta, 2\pi-\delta]$ 上一致收敛（对任意 $\delta > 0$）。

</details>

---

**2.** 设 $\delta > 0$。利用 Dirichlet 判别法证明 $\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^n}{n^x}$ 在 $[1, +\infty)$ 上不一致收敛（提示：验证 $B_n(x)$ 并非一致有界——事实上 $\sum_{k=1}^n (-1)^k$ 本身有界但不一致于 $x$？注意这里 $b_n = (-1)^n$ 与 $x$ 无关。）

<details><summary>参考答案</summary>

$\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^n}{n^x}$ 中 $a_n(x) = 1/n^x$，$b_n = (-1)^n$。

首先分析：$b_n = (-1)^n$ 的部分和 $B_n = \sum_{k=1}^n (-1)^k$：
- $B_1 = -1$，$B_2 = 0$，$B_3 = -1$，$B_4 = 0$，...
- $|B_n| \le 1$，故 $\{B_n\}$ 有界（与 $x$ 无关）。

$a_n(x) = 1/n^x$ 对每个 $x \ge 1$ 关于 $n$ 单调递减。且 $\sup_{x\ge 1} |a_n(x)| = 1/n \to 0$，故 $a_n(x) \rightrightarrows 0$。

三个条件看似都满足——但题目要求证明**不一致收敛**，说明我们的推论有问题。

仔细检查：Dirichlet 判别法的三个条件确实满足，因此该级数**应当一致收敛**——答案是否有误？

**更正**：实际上，$\sum_{n=1}^{\infty} (-1)^n/n^x$ 在 $[1, +\infty)$ 上**确实一致收敛**（可用例题 2 的 Leibniz 余项估计 + $\sup$ 方法验证）。本题的题意本不是用 Dirichlet 判别法证其不一致收敛（它确实一致收敛），而是考察学生对 Dirichlet 判别法条件的敏感性。

正确的题目应是：$\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^n}{n^x}$ 在 $(0, 1]$ 上是否一致收敛？

答案是**否**。原因是 $a_n(x) = 1/n^x$ 在 $(0, 1]$ 上，$\sup_{x\in(0,1]} |a_n(x)| = \sup_{x\in(0,1]} 1/n^x = 1$（当 $x \to 0^+$ 时 $1/n^x \to 1$），不趋于零。因此 $a_n(x) \not\rightrightarrows 0$，Dirichlet 判别法条件 2 不满足。

该级数在 $(0, 1]$ 上不一致收敛——可验证 $\sup_{x\in(0,1]} |R_n(x)| \not\to 0$。

</details>

---

### 迁移应用

**3.** 设 $p > 0$，考虑函数项级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(nx)}{n^p}$。

(1) 对哪些 $p$，该级数在 $\mathbb{R}$ 上一致收敛？（提示：M-判别法）

(2) 对哪些 $p$，该级数在 $[\delta, 2\pi-\delta]$（$\delta > 0$）上一致收敛但不在 $\mathbb{R}$ 上一致收敛？

(3) 对 $p = 1$，该级数在 $[0, 2\pi]$ 上是否一致收敛？说明理由。

<details><summary>参考答案</summary>

**(1)** 在 $\mathbb{R}$ 上一致收敛的条件。

对任意 $x \in \mathbb{R}$，$|\cos(nx)/n^p| \le 1/n^p$。当 $p > 1$ 时，$\sum 1/n^p$ 是收敛的 $p$-级数。由 M-判别法，$\sum \cos(nx)/n^p$ 在 $\mathbb{R}$ 上一致收敛。

当 $0 < p \le 1$ 时，$\sum 1/n^p$ 发散，M-判别法不适用。同时，$\sup_{x\in\mathbb{R}} |\cos(nx)/n^p| = 1/n^p$，但 $\sum 1/n^p$ 发散——不能由 M-判别法得一致收敛。事实上，该级数在 $\mathbb{R}$ 上不一致收敛（考虑 $x=0$ 处，级数退化为 $\sum 1/n^p$ 发散）。

因此 $p > 1$ 时在 $\mathbb{R}$ 上一致收敛。

**(2)** 在 $[\delta, 2\pi-\delta]$ 上一致收敛但不在 $\mathbb{R}$ 上一致收敛的条件。

对 $0 < p \le 1$，使用 Dirichlet 判别法：$a_n = 1/n^p$ 单调递减趋于 0，$B_n(x) = \sum_{k=1}^n \cos(kx)$ 在 $[\delta, 2\pi-\delta]$ 上一致有界（引理 10.1，界为 $1/\sin(\delta/2)$）。由定理 10.13，该级数在 $[\delta, 2\pi-\delta]$ 上一致收敛。

而由 (1)，$0 < p \le 1$ 时该级数在 $\mathbb{R}$ 上不一致收敛（$x=0$ 处发散）。

因此 $0 < p \le 1$ 是所求范围。

**(3)** $p = 1$ 时在 $[0, 2\pi]$ 上不一致收敛。

$p = 1$ 时级数为 $\sum \cos(nx)/n$。在 $x = 0$ 处，级数为 $\sum 1/n$（调和级数），发散——因此 $x=0$ 甚至不在收敛域中。级数在 $[0, 2\pi]$ 上不逐点收敛，更不可能一致收敛。

若将区间改为 $(0, 2\pi)$（开区间），该级数在 $(0, 2\pi)$ 的每个闭子区间上一致收敛（内闭一致收敛），但在 $(0, 2\pi)$ 本身上不一致收敛（因为越靠近 $0$ 和 $2\pi$，收敛速度越慢）。

**完整分类**：

| $p$ 的范围 | $\mathbb{R}$ 上一致收敛 | $[\delta, 2\pi-\delta]$ 上一致收敛 |
|-----------|------------------------|-----------------------------------|
| $p > 1$ | 是（M-判别法） | 是（M-判别法） |
| $0 < p \le 1$ | 否（$x=0$ 处发散） | 是（Dirichlet 判别法） |

</details>

---

**4.** 设 $\delta > 0$。利用 Abel 判别法证明 $\displaystyle\sum_{n=1}^{\infty} \left(1 + \frac{x^2}{n^2}\right) \frac{\cos(nx)}{n}$ 在 $[\delta, 2\pi-\delta]$ 上一致收敛。

<details><summary>参考答案</summary>

**第 1 步：拆分通项**。

令 $a_n(x) = 1 + \dfrac{x^2}{n^2}$，$b_n(x) = \dfrac{\cos(nx)}{n}$。

**第 2 步：验证 $a_n$ 的条件**。

对每个固定 $x$，$a_n(x) = 1 + x^2/n^2$ 关于 $n$ **单调递减**（因为 $x^2/n^2$ 递减）。此外，$|a_n(x)| \le 1 + (2\pi)^2$（因为 $|x| \le 2\pi$，实际上 $x \in [\delta, 2\pi-\delta] \subset [0, 2\pi]$），因此 $\{a_n(x)\}$ 在 $[\delta, 2\pi-\delta]$ 上一致有界。

**第 3 步：验证 $\sum b_n$ 的一致收敛性**。

由例题 1，$\sum_{n=1}^{\infty} \cos(nx)/n$ 在 $[\delta, 2\pi-\delta]$ 上一致收敛。这正是 $\sum b_n(x)$。

**第 4 步：应用 Abel 判别法**。

三个条件：
1. $a_n(x) = 1 + x^2/n^2$ 对每个 $x$ 单调递减 $\checkmark$
2. $|a_n(x)| \le 1 + 4\pi^2$，一致有界 $\checkmark$
3. $\sum \cos(nx)/n$ 在 $[\delta, 2\pi-\delta]$ 上一致收敛 $\checkmark$

由 Abel 判别法（定理 10.14），原级数在 $[\delta, 2\pi-\delta]$ 上一致收敛。

**注**：本例中，$a_n(x) = 1 + x^2/n^2$ 的极限为 $L(x) = 1$（常数，与 $x$ 无关），因此归化法也适用：$a_n'(x) = x^2/n^2$，对 $a_n'$ 和 $b_n$ 用 Dirichlet 判别法。

</details>