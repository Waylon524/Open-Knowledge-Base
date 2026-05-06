# 05. Abel 判别法与 Dirichlet 判别法（级数版本）

> 所属章节：第九章 数项级数  |  文件序号：05  |  难度：进阶
> 常见混淆点：Dirichlet 判别法要求 $a_n$ 单调趋于零 + $\sum b_k$ 的部分和有界；Abel 判别法要求 $a_n$ 单调有界 + $\sum b_n$ 收敛——两者条件恰好"互换"，注意不要用混；对于 $\sum \sin(n\theta)/n$ 型级数，$\theta$ 是 $\pi$ 的有理数倍时可通过直接计算前若干项观察周期性证明部分和有界，但更一般地需用三角恒等式 $2\sin(\theta/2)\sum_{k=1}^n\sin(k\theta) = \cos(\theta/2) - \cos((n+1/2)\theta)$ 来建立有界性

## 1. 学习目标与先修前置

### 学习目标
- 掌握 Abel 变换（分部求和公式 / summation by parts）的推导与两种等价形式，理解其作为离散版本"分部积分"的本质
- 掌握 Dirichlet 判别法（级数版本）的定理陈述，能用 Abel 变换 + Cauchy 收敛准则完成证明
- 掌握 Abel 判别法（级数版本）的定理陈述，能通过归化为 Dirichlet 判别法完成证明
- 能用三角恒等式证明 $\sum_{k=1}^n \sin(k\theta)$ 和 $\sum_{k=1}^n \cos(k\theta)$ 的有界性，结合 Dirichlet 判别法判断 $\sum \sin(n\theta)/n$ 和 $\sum \cos(n\theta)/n$ 型级数的收敛性
- 能区分 Dirichlet 判别法与 Abel 判别法的适用条件，并对给定级数选择正确的判别法
- 能与反常积分版本的 Abel-Dirichlet 判别法（ch8-04）进行平行类比，理解级数版本与积分版本的条件差异

### 先修知识
- 文件 01（第九章）：级数的基本概念（定义 9.1）、部分和 $S_n$、级数收敛 $\iff$ $\lim S_n$ 存在且有限
- 文件 04（第九章）：Cauchy 收敛准则（定理 9.10）、Leibniz 判别法（定理 9.11）、绝对收敛与条件收敛的定义（定义 9.9、9.10）、绝对收敛 $\Rightarrow$ 收敛定理（定理 9.12）
- 文件 10（第二章）：单调有界定理（定理 10.2）——用于证明 $\sum a_n' b_n$ 中 $a_n'$ 的单调性以及存在极限 $L$
- 文件 02（第九章）：比较判别法（定理 9.4）、极限比较判别法（定理 9.6）、$p$-级数（定理 9.5）
- 文件 04（第八章）：反常积分的 Abel-Dirichlet 判别法（积分第二中值定理、Dirichlet 判别法、Abel 判别法）——用于平行类比
- 三角恒等式：$2\sin\alpha\sin\beta = \cos(\alpha-\beta) - \cos(\alpha+\beta)$、$2\sin\alpha\cos\beta = \sin(\alpha+\beta) + \sin(\alpha-\beta)$

---

## 2. 背景与应用场景

### 2.1 Leibniz 判别法的局限

在文件 04 中，我们学习了 Leibniz 判别法——它完美解决了**交错级数**（$\sum (-1)^{n+1} a_n$）的收敛性问题。然而，许多重要的级数并不是交错级数：

$$\sum_{n=1}^{\infty} \frac{\sin n}{n}, \quad \sum_{n=1}^{\infty} \frac{\cos n}{n}, \quad \sum_{n=1}^{\infty} \frac{(-1)^n \sin n}{\sqrt{n}}$$

这些级数的通项正负交替**不严格**（$\sin n$ 的符号没有固定的周期性模式），Leibniz 判别法中的 $(-1)^{n+1}$ 结构不再存在。同时，这些级数也不是正项级数，比较判别法、比值判别法、根值判别法均不能直接使用。

然而，这些级数有一个共同的特点——它们都可以写成**乘积形式**：

$$\sum_{n=1}^{\infty} a_n b_n, \quad \text{其中 } a_n \text{ 单调递减趋于零}, \quad b_n \text{ 振荡（有界部分和）}$$

### 2.2 核心思想：分部求和与离散分部积分

处理这种"单调因子 $\times$ 振荡因子"乘积级数的核心工具是 **Abel 变换**——它是分部积分公式在离散情形下的类比：

| 连续版本（分部积分） | 离散版本（Abel 变换） |
|---------------------|---------------------|
| $\int_a^b u\,dv = uv\big|_a^b - \int_a^b v\,du$ | $\sum_{k=m+1}^n a_k b_k = a_n B_n - a_{m+1} B_m - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) B_k$ |

在连续场合，分部积分将一个积分转化为另一个"可能更容易处理"的积分；在离散场合，Abel 变换将 $\sum a_k b_k$ 转化为 $\sum (a_{k+1} - a_k) B_k$——如果 $a_k$ 单调（差分为定号）且 $B_k$ 有界，这个变换后的级数更容易估计。

### 2.3 与反常积分版本的平行关系

在 ch8-04 中，我们学习了反常积分的 Abel-Dirichlet 判别法，其逻辑链条为：

```
积分第二中值定理 → Dirichlet 判别法（∫f g 收敛）→ Abel 判别法（∫f g 收敛）
```

级数版本完全平行，只是将**积分第二中值定理**替换为**Abel 变换（分部求和公式）**：

```
Abel 变换 → Dirichlet 判别法（∑ a_n b_n 收敛）→ Abel 判别法（∑ a_n b_n 收敛）
```

两类判别法的条件结构完全一致：Dirichlet 判别法要求 $f$（或 $a_n$）单调趋于零、$\int g$（或 $\sum b_k$）的部分和有界；Abel 判别法则将"趋于零"放松为"有界"，将"有界"加强为"收敛"。

---

## 3. 核心概念与符号约定

### 3.1 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $\{a_n\}$ | 单调因子（通常单调递减） | 用于控制振荡幅度 |
| $\{b_n\}$ | 振荡因子（正负交替或更一般地变号） | 提供正负抵消 |
| $B_n = \sum_{k=1}^{n} b_k$ | $\{b_n\}$ 的部分和数列 | 约定 $B_0 = 0$ |
| $A_k = \sum_{i=1}^{k} a_i$ | $\{a_n\}$ 的部分和数列 | 用于第二种形式的 Abel 变换 |
| $a_k - a_{k+1}$ | $\{a_n\}$ 的差分（**递减时**非负） | 出现在 Abel 变换的尾部 |
| $b_k - b_{k+1}$ | $\{b_n\}$ 的差分 | 出现在 Abel 变换的另一种形式 |
| $L = \lim a_n$ | $\{a_n\}$ 的极限（Abel 判别法中有限） | $a_n' = a_n - L$ 趋于零 |
| $\theta$ | 三角级数中的频率参数 | 出现在 $\sum \sin(n\theta)/n$ 中 |
| $M$ | 有界性常数 | $|B_n| \le M$ |

### 3.2 关键概念

#### 3.2.1 乘积型级数

在实际问题中，级数常常以两个数列的乘积形式出现：

$$\sum_{n=1}^{\infty} a_n b_n$$

其中 $\{a_n\}$ 和 $\{b_n\}$ 分别扮演不同的角色：
- $\{a_n\}$ 通常是单调的（控制量级）
- $\{b_n\}$ 通常是振荡的（控制符号交替）

如果 $\{a_n\}$ 本身单调递减趋于零且 $\{b_n\}$ 的部分和有界，那么即使 $\sum |a_n b_n|$ 发散，正负抵消也可能使 $\sum a_n b_n$ 收敛——这正是 Dirichlet 判别法的核心思想。

#### 3.2.2 与积分版本的术语对应

| 级数版本 | 反常积分版本（ch8-04） |
|---------|----------------------|
| 数列 $\{a_n\}$ | 函数 $f(x)$ |
| 数列 $\{b_n\}$ | 函数 $g(x)$ |
| 部分和 $B_n = \sum_{k=1}^n b_k$ | 变上限积分 $G(A) = \int_a^A g(x)\,dx$ |
| 部分和有界 | 变上限积分一致有界 |
| Abel 变换 | 积分第二中值定理 |
| 单调递减 $\to 0$ | 单调趋于零 |

---

## 4. 原理与方法

### 4.1 Abel 变换（分部求和公式 / Summation by Parts）

**引理 9.13（Abel 变换）**：设 $\{a_n\}$ 和 $\{b_n\}$ 为两个数列，记 $B_n = \sum_{k=1}^{n} b_k$（约定 $B_0 = 0$）。则对任意 $n \ge 1$，有以下两个等价形式：

**形式一（常用形式）**：
$$\sum_{k=1}^{n} a_k b_k = a_n B_n + \sum_{k=1}^{n-1} (a_k - a_{k+1}) B_k \tag{9.13.1}$$

**形式二（尾部形式）**：对 $n > m \ge 0$，
$$\sum_{k=m+1}^{n} a_k b_k = a_n B_n - a_{m+1} B_m - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) B_k \tag{9.13.2}$$

**形式三（对偶形式）**：记 $A_n = \sum_{i=1}^{n} a_i$，则
$$\sum_{k=1}^{n} a_k b_k = A_n b_{n+1} + \sum_{k=1}^{n} A_k (b_k - b_{k+1}) \tag{9.13.3}$$

**证明（形式一）**：
从 $b_k = B_k - B_{k-1}$ 出发：

$$
\begin{aligned}
\sum_{k=1}^{n} a_k b_k &= \sum_{k=1}^{n} a_k (B_k - B_{k-1}) \\
&= \sum_{k=1}^{n} a_k B_k - \sum_{k=1}^{n} a_k B_{k-1} \\
&= \left( a_n B_n + \sum_{k=1}^{n-1} a_k B_k \right) - \left( \sum_{k=2}^{n} a_k B_{k-1} + a_1 B_0 \right) \\
&= a_n B_n + \sum_{k=1}^{n-1} a_k B_k - \sum_{k=1}^{n-1} a_{k+1} B_k \quad (\text{将第二和中的指标 } k \mapsto k+1) \\
&= a_n B_n + \sum_{k=1}^{n-1} (a_k - a_{k+1}) B_k
\end{aligned}
$$

这里最后一步用到了 $B_0 = 0$ 的约定。证毕。

**（形式三）证明**：
令 $A_k = \sum_{i=1}^k a_i$，约定 $A_0 = 0$。则 $a_k = A_k - A_{k-1}$。

$$
\begin{aligned}
\sum_{k=1}^{n} a_k b_k &= \sum_{k=1}^{n} (A_k - A_{k-1}) b_k \\
&= \sum_{k=1}^{n} A_k b_k - \sum_{k=1}^{n} A_{k-1} b_k \\
&= \sum_{k=1}^{n} A_k b_k - \sum_{k=0}^{n-1} A_k b_{k+1} \quad (\text{第二和指标 } k-1 \mapsto k) \\
&= A_n b_n + \sum_{k=1}^{n-1} A_k b_k - \sum_{k=1}^{n-1} A_k b_{k+1} - A_0 b_1 \\
&= A_n b_n + \sum_{k=1}^{n-1} A_k (b_k - b_{k+1}) + A_n (b_{n+1} - b_n) - A_n (b_{n+1} - b_n) \\
\end{aligned}
$$

另一种更直接的推导：

$$
\begin{aligned}
\sum_{k=1}^{n} A_k (b_k - b_{k+1}) &= \sum_{k=1}^{n} A_k b_k - \sum_{k=1}^{n} A_k b_{k+1} \\
&= \sum_{k=1}^{n} A_k b_k - \sum_{k=2}^{n+1} A_{k-1} b_k \\
&= A_1 b_1 + \sum_{k=2}^{n} (A_k - A_{k-1}) b_k - A_n b_{n+1} \\
&= a_1 b_1 + \sum_{k=2}^{n} a_k b_k - A_n b_{n+1} \\
&= \sum_{k=1}^{n} a_k b_k - A_n b_{n+1}
\end{aligned}
$$

移项即得 $\sum_{k=1}^{n} a_k b_k = A_n b_{n+1} + \sum_{k=1}^{n} A_k (b_k - b_{k+1})$。证毕。

**形式二（尾部形式）的推导**：
对 $n > m$，应用形式一在 $[1, n]$ 和 $[1, m]$ 上并相减：

$$
\begin{aligned}
\sum_{k=m+1}^{n} a_k b_k &= \left[ a_n B_n + \sum_{k=1}^{n-1} (a_k - a_{k+1}) B_k \right] - \left[ a_m B_m + \sum_{k=1}^{m-1} (a_k - a_{k+1}) B_k \right] \\
&= a_n B_n - a_m B_m + \sum_{k=m}^{n-1} (a_k - a_{k+1}) B_k
\end{aligned}
$$

当需要从 $a_{m+1}$ 开始时，可调整为：
$$\sum_{k=m+1}^{n} a_k b_k = a_n B_n - a_{m+1} B_m - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) B_k$$

这两种尾部形式在使用时根据便利性选择。

**分部积分类比**：形式三 $\sum a b = A_n b_{n+1} + \sum A (b_k - b_{k+1})$ 与分部积分公式 $\int u\,dv = uv - \int v\,du$ 对比如下：
- $a_k$（差分得到 $A_k$） $\leftrightarrow$ $dv$（积分得到 $v$）
- $b_k$（差分 $b_k - b_{k+1}$） $\leftrightarrow$ $u$（微分得到 $du$）
- 离散求和 $\sum$ $\leftrightarrow$ 连续积分 $\int$
- 边界项 $A_n b_{n+1}$ $\leftrightarrow$ 边界项 $uv\big|_a^b$

---

### 4.2 Dirichlet 判别法（级数版本）

#### 4.2.1 定理陈述

**定理 9.14（Dirichlet 判别法 / Dirichlet Test）**：设 $\{a_n\}$ 和 $\{b_n\}$ 为两个数列。若满足：

1. $\{a_n\}$ **单调递减**且 $\displaystyle\lim_{n\to\infty} a_n = 0$；
2. 数列 $\{B_n\}$ **有界**，其中 $B_n = \sum_{k=1}^{n} b_k$（即 $\{b_n\}$ 的部分和数列有界），

则级数 $\displaystyle\sum_{n=1}^{\infty} a_n b_n$ **收敛**。

**注**：
- 若 $\{a_n\}$ 单调**递增**趋于 0（即从负方向趋于 0），可令 $a_n' = -a_n$，则 $\{a_n'\}$ 单调递减趋于 0，且 $\sum a_n b_n = -\sum a_n' b_n$——敛散性相同。因此"单调递减"的条件不是本质限制。
- 与反常积分版本（ch8-04 定理 8.9）对比：积分版要求 $f(x)$ 在 $[a, +\infty)$ 上单调且 $\lim_{x\to\infty} f(x) = 0$、$G(A) = \int_a^A g$ 一致有界；级数版完全平行。

#### 4.2.2 证明

设 $M$ 为 $\{B_n\}$ 的一个上界，即 $|B_n| \le M$ 对所有 $n \ge 1$ 成立（由条件 2 可知这样的 $M$ 存在）。由条件 1，$\{a_n\}$ 单调递减趋于 0，因此 $a_n \ge 0$ 对所有 $n$ 成立（因为递减且极限为 0 意味着所有项非负）。

我们要证明 $\sum a_n b_n$ 收敛。使用 Cauchy 收敛准则（定理 9.10）：对任意 $\varepsilon > 0$，寻找 $N$ 使得对所有 $n > m \ge N$ 有 $\left|\sum_{k=m+1}^{n} a_k b_k\right| < \varepsilon$。

**第 1 步：应用 Abel 变换于尾部。**

利用 Abel 变换的尾部形式（引理 9.13 形式二），对 $n > m$：

$$\sum_{k=m+1}^{n} a_k b_k = a_n B_n - a_{m+1} B_m + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) B_k$$

注意这里同乘以 $B_k$ 的是 $(a_k - a_{k+1})$ 而非 $(a_{k+1} - a_k)$，且符号为正——这是由形式一直接导出的结果。

**第 2 步：对每一项取绝对值并用三角不等式放缩。**

$$
\begin{aligned}
\left|\sum_{k=m+1}^{n} a_k b_k\right|
&\le |a_n B_n| + |a_{m+1} B_m| + \sum_{k=m+1}^{n-1} |(a_k - a_{k+1}) B_k| \\
&\le a_n |B_n| + a_{m+1} |B_m| + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) |B_k|
\end{aligned}
$$

在第三项中，由于 $\{a_n\}$ 递减，$a_k - a_{k+1} \ge 0$，因此绝对值可去掉。而 $a_n \ge 0$。

**第 3 步：利用有界性 $|B_k| \le M$ 进行整体放缩。**

$$
\begin{aligned}
\left|\sum_{k=m+1}^{n} a_k b_k\right|
&\le M \left( a_n + a_{m+1} + \sum_{k=m+1}^{n-1} (a_k - a_{k+1}) \right)
\end{aligned}
$$

**第 4 步：计算求和项——这是一个 telescoping sum。**

$$
\sum_{k=m+1}^{n-1} (a_k - a_{k+1}) = (a_{m+1} - a_{m+2}) + (a_{m+2} - a_{m+3}) + \cdots + (a_{n-1} - a_n) = a_{m+1} - a_n
$$

代入：

$$
\begin{aligned}
\left|\sum_{k=m+1}^{n} a_k b_k\right|
&\le M \big( a_n + a_{m+1} + (a_{m+1} - a_n) \big) \\
&= M (2a_{m+1}) = 2M a_{m+1}
\end{aligned}
$$

**第 5 步：用 $a_n \to 0$ 控制 $a_{m+1}$。**

由 $\lim a_n = 0$，对给定的 $\varepsilon > 0$，存在 $N$ 使得当 $n \ge N$ 时 $a_n < \dfrac{\varepsilon}{2M}$。由于 $\{a_n\}$ 递减，对任意 $m \ge N$，有 $a_{m+1} \le a_N < \dfrac{\varepsilon}{2M}$。

因此当 $n > m \ge N$ 时：

$$\left|\sum_{k=m+1}^{n} a_k b_k\right| \le 2M a_{m+1} < 2M \cdot \frac{\varepsilon}{2M} = \varepsilon$$

由 Cauchy 收敛准则，级数 $\sum_{n=1}^{\infty} a_n b_n$ 收敛。证毕。

**证明逻辑总结**：

```
条件 2: B_n 有界 (|B_n| ≤ M)
  + 条件 1: a_n 递减 → 0
  + Abel 变换
  ↓
|∑ a_k b_k| ≤ M (a_n + a_{m+1} + (a_{m+1} - a_n)) = 2M a_{m+1}
  + a_n → 0  ⇒  a_{m+1} < ε/(2M) 对充分大的 m
  ↓
Cauchy 准则 ⇒ ∑ a_n b_n 收敛
```

#### 4.2.3 关键工具：$\sum \sin(k\theta)$ 和 $\sum \cos(k\theta)$ 的有界性

应用 Dirichlet 判别法的核心是验证条件 2——$\{B_n\}$ 有界。对于 $\sum \sin(n\theta)/n$ 和 $\sum \cos(n\theta)/n$ 型级数，$b_n = \sin(n\theta)$ 或 $b_n = \cos(n\theta)$。我们需证明 $\sum_{k=1}^n \sin(k\theta)$ 和 $\sum_{k=1}^n \cos(k\theta)$ 有界。

**公式推导**：利用积化和差恒等式。

对 $\sin$ 求和：
$$
\begin{aligned}
2\sin\frac{\theta}{2} \sum_{k=1}^{n} \sin(k\theta) &= \sum_{k=1}^{n} 2\sin\frac{\theta}{2}\sin(k\theta) \\
&= \sum_{k=1}^{n} \left[ \cos\left((k-\tfrac12)\theta\right) - \cos\left((k+\tfrac12)\theta\right) \right] \\
&= \cos\frac{\theta}{2} - \cos\left((n+\tfrac12)\theta\right)
\end{aligned}
$$

因此，当 $\theta \neq 2m\pi$（$m \in \mathbb{Z}$）时：

$$\sum_{k=1}^{n} \sin(k\theta) = \frac{\cos\frac{\theta}{2} - \cos\left((n+\frac12)\theta\right)}{2\sin\frac{\theta}{2}} \tag{9.14.1}$$

同理，对 $\cos$ 求和：
$$
\begin{aligned}
2\sin\frac{\theta}{2} \sum_{k=1}^{n} \cos(k\theta) &= \sum_{k=1}^{n} 2\sin\frac{\theta}{2}\cos(k\theta) \\
&= \sum_{k=1}^{n} \left[ \sin\left((k+\tfrac12)\theta\right) - \sin\left((k-\tfrac12)\theta\right) \right] \\
&= \sin\left((n+\tfrac12)\theta\right) - \sin\frac{\theta}{2}
\end{aligned}
$$

因此，当 $\theta \neq 2m\pi$ 时：

$$\sum_{k=1}^{n} \cos(k\theta) = \frac{\sin\left((n+\frac12)\theta\right) - \sin\frac{\theta}{2}}{2\sin\frac{\theta}{2}} \tag{9.14.2}$$

**有界性结论**：从 (9.14.1) 和 (9.14.2) 可得：

$$\left|\sum_{k=1}^{n} \sin(k\theta)\right| \le \frac{1}{|\sin\frac{\theta}{2}|}, \quad \left|\sum_{k=1}^{n} \cos(k\theta)\right| \le \frac{1}{|\sin\frac{\theta}{2}|} \tag{9.14.3}$$

右端与 $n$ 无关，因此 $\{B_n\}$ 有界。

当 $\theta = 2m\pi$ 时，$\sin(k\theta) = \sin(2mk\pi) = 0$，$\cos(k\theta) = \cos(2mk\pi) = 1$。此时：
- $\sum \sin(k\theta) = 0$（平凡有界）
- $\sum \cos(k\theta) = n$（**无界**！——此时 $b_n = 1$，$\sum b_n$ 本身发散，不能使用 Dirichlet 判别法）

因此 Dirichlet 判别法适用于 $\sum \sin(n\theta)/n$ 对所有 $\theta$ 成立，但对 $\sum \cos(n\theta)/n$ 仅当 $\theta \neq 2m\pi$ 时成立。

**特例验证**：当 $\theta = \pi/4$ 时，$\sin(\pi/8) = \sqrt{2-\sqrt{2}}/2$，上界为 $1/\sin(\pi/8) \approx 2.613$。当 $\theta = \pi/3$ 时，$\sin(\pi/6) = 1/2$，上界为 $1/(1/2) = 2$。这与前序文件中通过直接计算前若干项观察周期性得到的有界性结论一致。

---

### 4.3 Abel 判别法（级数版本）

#### 4.3.1 定理陈述

**定理 9.15（Abel 判别法 / Abel Test）**：设 $\{a_n\}$ 和 $\{b_n\}$ 为两个数列。若满足：

1. $\{a_n\}$ **单调且有界**；
2. 级数 $\displaystyle\sum_{n=1}^{\infty} b_n$ **收敛**，

则级数 $\displaystyle\sum_{n=1}^{\infty} a_n b_n$ **收敛**。

**注**：
- 与 Dirichlet 判别法对比：Dirichlet 要求 $a_n \to 0$（比"有界"强）且 $B_n$ 有界（比"$\sum b_n$ 收敛"弱）；Abel 恰好互换——$a_n$ 只须有界（不必趋于 0），但 $\sum b_n$ 须收敛（比"有界"强）。两者是对偶关系。
- 与反常积分版本（ch8-04 定理 8.10）对比：条件结构完全一致。

#### 4.3.2 证明（归化法——通过 Dirichlet 判别法）

这是最简洁的证法，将 Abel 判别法归化为已经证明的 Dirichlet 判别法。

**证明**：

由条件 1，$\{a_n\}$ 单调且有界。由单调有界定理（第二章定理 10.2），$\{a_n\}$ 存在极限，记

$$L = \lim_{n\to\infty} a_n$$

（$L$ 为有限实数。）定义新数列 $\{a_n'\}$：

$$a_n' = a_n - L, \quad n = 1, 2, 3, \dots$$

则：
- $\{a_n'\}$ 与 $\{a_n\}$ 具有相同的单调性（减去常数不改变单调性）
- $\displaystyle\lim_{n\to\infty} a_n' = \lim_{n\to\infty} (a_n - L) = L - L = 0$

由条件 2，$\sum b_n$ 收敛。收敛级数的部分和数列 $\{B_n\}$ 必然有界（收敛数列必有界）。因此 $\{B_n\}$ 有界。

现在对数列 $\{a_n'\}$ 和 $\{b_n\}$ 应用 Dirichlet 判别法（定理 9.14）：
- $\{a_n'\}$ 单调且 $\lim a_n' = 0$ ✓
- $\{B_n\}$ 有界 ✓

故级数 $\sum_{n=1}^{\infty} a_n' b_n$ 收敛。

最后：

$$\sum_{n=1}^{\infty} a_n b_n = \sum_{n=1}^{\infty} (a_n' + L) b_n = \sum_{n=1}^{\infty} a_n' b_n + L \sum_{n=1}^{\infty} b_n$$

其中 $\sum a_n' b_n$ 收敛（由 Dirichlet 判别法），$L\sum b_n$ 也收敛（常数倍收敛级数）。两个收敛级数之和收敛，故 $\sum a_n b_n$ 收敛。证毕。

#### 4.3.3 直接证明（通过 Abel 变换 + Cauchy 收敛准则）

除归化法外，也可直接使用 Abel 变换证明。这里给出一种与 Dirichlet 判别法对称的论证。

**直接证明**：

由条件 2，$\sum b_n$ 收敛。记 $S = \sum_{n=1}^{\infty} b_n$，部分和 $B_n = \sum_{k=1}^{n} b_k$，余项 $R_n = S - B_n = \sum_{k=n+1}^{\infty} b_k$。由收敛性知 $\lim_{n\to\infty} R_n = 0$。

由条件 1，$\{a_n\}$ 单调有界，记其上界为 $M$（即 $|a_n| \le M$ 对所有 $n$ 成立），且极限 $L = \lim a_n$ 存在。

选取 Abel 变换的展开形式。注意到 $b_k = (B_k - B_{k-1})$，对 $n > m$ 有：

$$\sum_{k=m+1}^{n} a_k b_k = a_n B_n - a_{m+1} B_m - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) B_k$$

将 $B_k = S - R_k$ 代入：

$$
\begin{aligned}
\sum_{k=m+1}^{n} a_k b_k &= a_n (S - R_n) - a_{m+1} (S - R_m) - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k)(S - R_k) \\
&= S \left( a_n - a_{m+1} - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) \right) - a_n R_n + a_{m+1} R_m + \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) R_k
\end{aligned}
$$

注意到 $\sum_{k=m+1}^{n-1} (a_{k+1} - a_k) = a_n - a_{m+1}$（telescoping sum），因此方括号中的项为零：

$$a_n - a_{m+1} - \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) = a_n - a_{m+1} - (a_n - a_{m+1}) = 0$$

于是：

$$\sum_{k=m+1}^{n} a_k b_k = - a_n R_n + a_{m+1} R_m + \sum_{k=m+1}^{n-1} (a_{k+1} - a_k) R_k$$

取绝对值，利用 $|a_k| \le M$ 以及 $a_{k+1} - a_k$ 在 $\{a_n\}$ 单调时定号的特性：

$$
\begin{aligned}
\left|\sum_{k=m+1}^{n} a_k b_k\right|
&\le M|R_n| + M|R_m| + \sum_{k=m+1}^{n-1} |a_{k+1} - a_k| \cdot |R_k|
\end{aligned}
$$

由于 $R_k \to 0$，对任意 $\varepsilon > 0$，存在 $N_1$ 使得 $k \ge N_1$ 时 $|R_k| < \varepsilon/(3M)$。取 $N \ge N_1$，则对 $n > m \ge N$：

$$
\begin{aligned}
\left|\sum_{k=m+1}^{n} a_k b_k\right|
&\le M \cdot \frac{\varepsilon}{3M} + M \cdot \frac{\varepsilon}{3M} + \frac{\varepsilon}{3M} \sum_{k=m+1}^{n-1} |a_{k+1} - a_k| \\
&= \frac{2\varepsilon}{3} + \frac{\varepsilon}{3M} \sum_{k=m+1}^{n-1} |a_{k+1} - a_k|
\end{aligned}
$$

最后一项中，若 $\{a_n\}$ 单调递减，则 $a_{k+1} - a_k \le 0$，$|a_{k+1} - a_k| = a_k - a_{k+1}$，求和得 $a_{m+1} - a_n$。由有界性，$a_{m+1} - a_n \le 2M$。因此：

$$\sum_{k=m+1}^{n-1} |a_{k+1} - a_k| \le 2M$$

代入得：

$$\left|\sum_{k=m+1}^{n} a_k b_k\right| \le \frac{2\varepsilon}{3} + \frac{\varepsilon}{3M} \cdot 2M = \frac{2\varepsilon}{3} + \frac{2\varepsilon}{3} = \frac{4\varepsilon}{3}$$

这个上界未能压缩到 $\varepsilon$ 以内——调整系数的选取可以解决（例如选取 $N$ 使得 $|R_k| < \varepsilon/(4M)$ 即可得到 $\left|\sum a_k b_k\right| < \varepsilon$）。详细的系数调整留给读者作为练习。该直接证明展示了与归化法不同的视角——通过余项 $R_k$ 的衰减而非 $B_k$ 的有界性来估计尾部。

**实际使用时，我们推荐归化法**，因其简洁且无需调整系数。

#### 4.3.4 两类判别法的条件对比

| 对比维度 | Dirichlet 判别法（定理 9.14） | Abel 判别法（定理 9.15） |
|----------|-----------------------------|-------------------------|
| 对 $a_n$ 的要求 | 单调，**趋于零** | 单调，**有界**（不必趋于 0） |
| 对 $b_n$ 的要求 | $\sum_{k=1}^{n} b_k$ **有界** | $\sum_{n=1}^{\infty} b_n$ **收敛** |
| 谁更强 | $a_n \to 0$ 比有界强；但 $B_n$ 有界比 $\sum b_n$ 收敛弱 | $a_n$ 有界比趋于零弱；但 $\sum b_n$ 收敛比有界强 |
| 适用场景 | $a_n$ 衰减 + $b_n$ 振荡 | $a_n$ 有界变化 + $b_n$ 收敛 |
| 证明核心工具 | Abel 变换 + Cauchy 准则 | 归化为 Dirichlet（或直接 Abel 变换） |

### 4.4 级数版本与积分版本的条件对比

| 对比维度 | 级数版本 | 反常积分版本（ch8-04） |
|---------|---------|---------------------|
| 核心引理 | Abel 变换（引理 9.13） | 积分第二中值定理 |
| Dirichlet 判据 | $a_n \downarrow 0$，$\sum b_k$ 部分和有界 | $f(x) \downarrow 0$，$\int_a^A g$ 一致有界 |
| Abel 判据 | $a_n$ 单调有界，$\sum b_n$ 收敛 | $f$ 单调有界，$\int_a^\infty g$ 收敛 |
| Cauchy 准则 | 级数 Cauchy 准则（定理 9.10） | 反常积分 Cauchy 准则（定理 8.6） |
| 应用示例 | $\sum \frac{\sin n}{n}$ | $\int_1^\infty \frac{\sin x}{x}\,dx$ |

---

## 5. 例题

### 例题 1：Dirichlet 判别法（$\sum \sin(n\theta)/n$ 型）

判断级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(n\pi/4)}{n}$ 的敛散性。

**解**：

**第 1 步：观察结构，选取判别法**。

级数通项为 $\frac{\sin(n\pi/4)}{n}$，可写为 $a_n b_n$ 的形式，其中 $a_n = \frac{1}{n}$，$b_n = \sin(n\pi/4)$。

- 它不是交错级数（$\sin(n\pi/4)$ 不严格按 $(-1)^{n+1}$ 的规律交替）
- 它不是正项级数（$\sin(n\pi/4)$ 可取正、负、零）
- 但 $a_n = 1/n$ 单调递减趋于零，且 $b_n = \sin(n\pi/4)$ 的**部分和**有界——符合 Dirichlet 判别法的结构

因此选用 **Dirichlet 判别法（定理 9.14）**。

**第 2 步：验证条件 1——$\{a_n\}$ 单调递减趋于零**。

$a_n = \frac{1}{n}$（$n \ge 1$）。显然：
- $a_{n+1} = \frac{1}{n+1} < \frac{1}{n} = a_n$，故 $\{a_n\}$ 严格单调递减
- $\displaystyle\lim_{n\to\infty} \frac{1}{n} = 0$
- 条件满足 $\checkmark$

**第 3 步：验证条件 2——$\{B_n\} = \{\sum_{k=1}^n \sin(k\pi/4)\}$ 有界**。

**方法一（直接计算 + 周期性）**：

$\sin(k\pi/4)$ 的周期为 $2\pi/(\pi/4) = 8$。计算前 8 项：

| $k$ | $k\pi/4$ | $\sin(k\pi/4)$ | $B_k = \sum_{i=1}^k \sin(i\pi/4)$ |
|-----|----------|----------------|-----------------------------------|
| 1 | $\pi/4$ | $\sqrt{2}/2$ | $\sqrt{2}/2$ |
| 2 | $\pi/2$ | $1$ | $\sqrt{2}/2 + 1$ |
| 3 | $3\pi/4$ | $\sqrt{2}/2$ | $\sqrt{2} + 1$ |
| 4 | $\pi$ | $0$ | $\sqrt{2} + 1$ |
| 5 | $5\pi/4$ | $-\sqrt{2}/2$ | $\sqrt{2}/2 + 1$ |
| 6 | $3\pi/2$ | $-1$ | $\sqrt{2}/2$ |
| 7 | $7\pi/4$ | $-\sqrt{2}/2$ | $0$ |
| 8 | $2\pi$ | $0$ | $0$ |

$B_8 = 0$，且由周期性，$B_{8+r} = B_8 + \sum_{i=1}^r \sin(i\pi/4) = B_r$。故 $\{B_n\}$ 以 8 为周期，只取有限多个值：

$$\{0,\ \frac{\sqrt{2}}{2},\ \frac{\sqrt{2}}{2}+1,\ \sqrt{2}+1\}$$

因此 $\{B_n\}$ 有界，$|B_n| \le \sqrt{2} + 1 \approx 2.414$。$\checkmark$

**方法二（三角恒等式法）**：

用公式 (9.14.1)：$\theta = \pi/4$，$\sin(\theta/2) = \sin(\pi/8) = \sqrt{2-\sqrt{2}}/2$。

$$\left|\sum_{k=1}^{n} \sin(k\pi/4)\right| \le \frac{1}{\sin(\pi/8)} = \frac{2}{\sqrt{2-\sqrt{2}}} \approx 2.613$$

此上界虽然比方法一稍宽松，但适用于所有 $\theta$（无需周期性），是更一般的方法。

**第 4 步：应用 Dirichlet 判别法，得出结论**。

两个条件均满足，由 Dirichlet 判别法（定理 9.14），级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(n\pi/4)}{n}$ **收敛**。

**注意**：该级数**不绝对收敛**（$\sum |\sin(n\pi/4)|/n$ 发散，因为 $|\sin(n\pi/4)|$ 对约一半的 $n$ 取正值，与 $\sum 1/n$ 同量级），因此属于**条件收敛**。

---

### 例题 2：Dirichlet 判别法（$\sum \cos(n\theta)/n$ 型）

判断级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(n\pi/3)}{n}$ 的敛散性。

**解**：

**第 1 步：观察结构，选取判别法**。

令 $a_n = \frac{1}{n}$，$b_n = \cos(n\pi/3)$。$a_n \downarrow 0$，需要验证 $\{\sum \cos(k\pi/3)\}$ 有界。

**第 2 步：验证条件 1——$\{a_n\}$ 单调递减趋于零**。

同例题 1，$\{1/n\}$ 严格单调递减趋于 0。$\checkmark$

**第 3 步：验证条件 2——$\{B_n\} = \{\sum_{k=1}^n \cos(k\pi/3)\}$ 有界**。

**方法一（直接计算 + 周期性）**：

$\cos(k\pi/3)$ 的周期为 $2\pi/(\pi/3) = 6$。计算前 6 项：

| $k$ | $k\pi/3$ | $\cos(k\pi/3)$ | $B_k = \sum_{i=1}^k \cos(i\pi/3)$ |
|-----|----------|----------------|-----------------------------------|
| 1 | $\pi/3$ | $1/2$ | $1/2$ |
| 2 | $2\pi/3$ | $-1/2$ | $0$ |
| 3 | $\pi$ | $-1$ | $-1$ |
| 4 | $4\pi/3$ | $-1/2$ | $-3/2$ |
| 5 | $5\pi/3$ | $1/2$ | $-1$ |
| 6 | $2\pi$ | $1$ | $0$ |

$B_6 = 0$。由周期性，$\{B_n\}$ 以 6 为周期，取值均在 $[-3/2, 1/2]$ 内，故有界：$|B_n| \le \frac{3}{2}$。$\checkmark$

**方法二（三角恒等式法）**：

用公式 (9.14.2)：$\theta = \pi/3$，$\sin(\theta/2) = \sin(\pi/6) = 1/2$。

$$\left|\sum_{k=1}^{n} \cos(k\pi/3)\right| \le \frac{1}{\sin(\pi/6)} = \frac{1}{1/2} = 2$$

**第 4 步：应用 Dirichlet 判别法，得出结论**。

两个条件均满足，由 Dirichlet 判别法，级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(n\pi/3)}{n}$ **收敛**（条件收敛）。

---

### 例题 3：Abel 判别法

设 $a_n = 1 + \dfrac{1}{n}$，$b_n = \dfrac{(-1)^n}{n^2}$。判断级数 $\displaystyle\sum_{n=1}^{\infty} a_n b_n = \sum_{n=1}^{\infty} \left(1 + \frac{1}{n}\right)\frac{(-1)^n}{n^2}$ 的敛散性。

**解**：

**第 1 步：观察结构，选取判别法**。

级数通项为 $(1 + 1/n)\cdot (-1)^n/n^2$。

先看能否用 Dirichlet 判别法：$a_n = 1 + 1/n$ 单调递减但**不趋于零**（它趋于 1），不符合 Dirichlet 的条件 1。

再看 $\{b_n\} = \{(-1)^n/n^2\}$ 的级数：$\sum (-1)^n/n^2$ 是 $p = 2 > 1$ 的 $p$-级数取绝对值后收敛，故绝对收敛。

这恰好符合 **Abel 判别法（定理 9.15）** 的结构——$a_n$ 单调有界、$\sum b_n$ 收敛。

**第 2 步：验证条件 1——$\{a_n\}$ 单调且有界**。

$a_n = 1 + 1/n$（$n \ge 1$）。

- **单调性**：$\frac{1}{n+1} < \frac{1}{n}$，故 $a_{n+1} = 1 + \frac{1}{n+1} < 1 + \frac{1}{n} = a_n$。$\{a_n\}$ 严格单调递减。$\checkmark$

- **有界性**：对 $n \ge 1$，$1 < 1 + \frac{1}{n} \le 2$（$n = 1$ 时最大为 2）。因此 $\{a_n\}$ 有上下界，$|a_n| \le 2$。$\checkmark$

**第 3 步：验证条件 2——$\sum b_n$ 收敛**。

$$\sum_{n=1}^{\infty} b_n = \sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}$$

考察绝对值级数 $\sum |b_n| = \sum 1/n^2$。这是 $p$-级数，$p = 2 > 1$，由定理 9.5 知 $\sum 1/n^2$ 收敛。故 $\sum (-1)^n/n^2$ **绝对收敛**（定义 9.9），由定理 9.12 知它亦收敛。$\checkmark$

**第 4 步：应用 Abel 判别法，得出结论**。

两个条件均满足，由 Abel 判别法（定理 9.15），级数 $\displaystyle\sum_{n=1}^{\infty} \left(1 + \frac{1}{n}\right)\frac{(-1)^n}{n^2}$ **收敛**。

**验证**（利用绝对收敛 + 放缩法作为补充检查）：由于 $|(1+1/n)\cdot(-1)^n/n^2| \le 2/n^2$ 且 $\sum 2/n^2$ 收敛，由比较判别法知原级数**绝对收敛**——这一论证无需 Abel 判别法，但仅当 $a_n$ 本身乘上 $|b_n|$ 仍可被比较级数控制时才可行。当 $a_n$ 不再有界、或 $b_n$ 不绝对收敛时，绝对收敛放缩法失效，必须依赖 Abel 判别法。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| Dirichlet 判别法和 Abel 判别法都可以处理 $\sum (-1)^{n+1} a_n$ 型交错级数 | 对于交错级数，Leibniz 判别法已经足够且条件更简单（只需 $a_n \downarrow 0$）。Abel/Dirichlet 判别法用于更一般的乘积型级数 $\sum a_n b_n$，$b_n$ 不限于 $(-1)^{n+1}$ |
| 用 Dirichlet 判别法时，$a_n$ 必须递减到 0 且 $b_n$ 自身有界（即 $|b_n| \le M$） | 条件 2 要求的是 $\{B_n\} = \{\sum_{k=1}^n b_k\}$ 有界而非 $\{b_n\}$ 有界。这比"$b_n$ 有界"更宽松——$b_n$ 本身可以无界，只要其部分和数列有界即可 |
| Abel 判别法和 Dirichlet 判别法是完全独立的两个定理 | Abel 判别法可以通过构造 $a_n' = a_n - L$ 归化为 Dirichlet 判别法。二者共享 Abel 变换作为核心工具，条件形成对偶关系 |
| 只需 $|B_n| \le M$ 和 $a_n \to 0$ 就自动得到 $|\sum a_k b_k| \le 2M a_N$ | 推导中关键步用到了 $a_k - a_{k+1} \ge 0$（即 $\{a_n\}$ 单调递减）。若 $\{a_n\}$ 不单调，则 $(a_k - a_{k+1})$ 可正可负，telescoping sum 不再有 $a_{m+1} - a_n$ 的简洁形式 |
| 对任意 $\theta$，$\sum \cos(n\theta)/n$ 都收敛 | 当 $\theta = 2m\pi$ 时，$\cos(n\theta) = 1$，级数退化为 $\sum 1/n$（调和级数），发散！因为此时 $\sum \cos(k\theta) = n$ 无界，Dirichlet 判别法的条件 2 不满足 |
| Abel 变换中 $\sum_{k=1}^{n-1} (a_k - a_{k+1}) B_k$ 是一个与 $\sum a_k b_k$ 类似的新级数，没有简化问题 | 关键区别在于：$(a_k - a_{k+1})$ 是定号的（由单调性保证），且 telescoping sum $\sum (a_k - a_{k+1}) = a_1 - a_n$ 可以被控制，而 $B_k$ 有界——因此整个和式可以用 $\sup|B_k| \cdot (a_1 - a_n)$ 简单上界，实现了"化乘积和为上界乘积" |

### 检查点

- [ ] 能否写出 Abel 变换的两种等价形式（求和公式 9.13.1 和 9.13.3）并完成推导？
- [ ] 能否完整叙述 Dirichlet 判别法（定理 9.14）的两个条件并用 Abel 变换 + Cauchy 收敛准则给出证明？
- [ ] 能否完整叙述 Abel 判别法（定理 9.15）的两个条件并用归化法（通过 Dirichlet 判别法）给出证明？
- [ ] 能否说明为什么在 Dirichlet 判别法的证明中，$\{a_n\}$ 的单调性是不可或缺的？
- [ ] 能否推导并证明 $\sum_{k=1}^n \sin(k\theta)$ 和 $\sum_{k=1}^n \cos(k\theta)$ 的有界性公式（9.14.1/9.14.2）？
- [ ] 能否指出 $\theta$ 取何值时 $\sum \cos(n\theta)/n$ 不收敛，并解释原因？
- [ ] 能否判断级数 $\sum \sin(n\theta)/n$ 对任意 $\theta$ 是否都收敛？
- [ ] 能否区分 Dirichlet 判别法和 Abel 判别法的适用条件差异？
- [ ] 能用 Dirichlet 判别法判断 $\sum \sin(n\pi/4)/n$ 的敛散性吗？
- [ ] 能用 Abel 判别法判断 $\sum (1+1/n)(-1)^n/n^2$ 的敛散性吗？
- [ ] 能否说出级数版本的 Abel-Dirichlet 判别法与反常积分版本（ch8-04）的平行关系？
- [ ] 能否举出一个必须用 Dirichlet 判别法（而不能用其他任何前序判别法）的例子？

---

## 练习题

### 基础巩固

**1.** 完成 Abel 变换的推导填空：

(1) 已知 $B_n = \sum_{k=1}^n b_k$（$B_0 = 0$）。证明：
$$\sum_{k=1}^{n} a_k b_k = a_n B_n + \sum_{k=1}^{n-1} (a_k - a_{k+1}) B_k$$

(2) 已知 $A_n = \sum_{k=1}^n a_k$（$A_0 = 0$）。证明：
$$\sum_{k=1}^{n} a_k b_k = A_n b_{n+1} + \sum_{k=1}^{n} A_k (b_k - b_{k+1})$$

<details><summary>参考答案</summary>

**(1)** 见引理 9.13 中形式一的推导。核心步骤：
- 将 $b_k$ 写为 $B_k - B_{k-1}$
- 展开乘积和，重新组合指标
- 合并同类项得 $(a_k - a_{k+1})B_k$

**(2)** 见引理 9.13 中形式三的推导。核心步骤：
- 将 $a_k$ 写为 $A_k - A_{k-1}$
- 展开乘积和，重新组合指标
- 合并同类项得 $A_k(b_k - b_{k+1})$

两种形式互为对偶，区别在于："求和 $a$ 差分 $b$"（形式三） vs "求和 $b$ 差分 $a$"（形式一）。

</details>

---

**2.** 判断下列级数的敛散性。说明使用的判别法并逐条验证条件。

(1) $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(n\pi/6)}{n}$

(2) $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(n\pi/2)}{n}$

(3) $\displaystyle\sum_{n=1}^{\infty} \left(2 + \frac{(-1)^n}{n}\right) \cdot \frac{(-1)^n}{n^{3/2}}$

<details><summary>参考答案</summary>

**(1)** $\displaystyle\sum_{n=1}^{\infty} \frac{\sin(n\pi/6)}{n}$

令 $a_n = 1/n$，$b_n = \sin(n\pi/6)$。使用 Dirichlet 判别法。

**条件 1**：$\{a_n\} = \{1/n\}$ 递减趋于 0。$\checkmark$

**条件 2**：$\sum_{k=1}^n \sin(k\pi/6)$ 有界。

用三角恒等式法：$\theta = \pi/6$，$\sin(\theta/2) = \sin(\pi/12) = \sqrt{2-\sqrt{3}}/2$。

$$\left|\sum_{k=1}^n \sin(k\pi/6)\right| \le \frac{1}{\sin(\pi/12)} = \frac{2}{\sqrt{2-\sqrt{3}}}$$

右端有限（与 $n$ 无关），故有界。$\checkmark$

由 Dirichlet 判别法，级数收敛（条件收敛）。

---

**(2)** $\displaystyle\sum_{n=1}^{\infty} \frac{\cos(n\pi/2)}{n}$

令 $a_n = 1/n$，$b_n = \cos(n\pi/2)$。尝试使用 Dirichlet 判别法。

计算 $\cos(n\pi/2)$ 的前几项：$\cos(\pi/2) = 0$，$\cos(\pi) = -1$，$\cos(3\pi/2) = 0$，$\cos(2\pi) = 1$，$\cos(5\pi/2) = 0$，$\cos(3\pi) = -1$，$\dots$

部分和 $B_n = \sum_{k=1}^n \cos(k\pi/2)$：
$B_1 = 0$，$B_2 = -1$，$B_3 = -1$，$B_4 = 0$，$B_5 = 0$，$B_6 = -1$，$\dots$

$\{B_n\}$ 以 4 为周期，仅取值 $0$ 和 $-1$，故有界（$|B_n| \le 1$）。$\checkmark$

条件 1 显然满足。由 Dirichlet 判别法，级数收敛。

---

**(3)** $\displaystyle\sum_{n=1}^{\infty} \left(2 + \frac{(-1)^n}{n}\right) \cdot \frac{(-1)^n}{n^{3/2}}$

令 $a_n = 2 + (-1)^n/n$，$b_n = (-1)^n/n^{3/2}$。尝试用 Abel 判别法。

**条件 1**：$\{a_n\}$ 是否单调？
$a_1 = 2 - 1 = 1$，$a_2 = 2 + 1/2 = 2.5$，$a_3 = 2 - 1/3 \approx 1.667$，$a_4 = 2 + 1/4 = 2.25$，$\dots$
$\{a_n\}$ 在 $1$ 到 $2.5$ 之间振荡，**不单调**。Abel 判别法不适用。

改用绝对收敛放缩法：
$$|a_n b_n| = \left|\left(2 + \frac{(-1)^n}{n}\right)\frac{(-1)^n}{n^{3/2}}\right| \le \left(2 + \frac{1}{n}\right)\frac{1}{n^{3/2}} \le \frac{3}{n^{3/2}} \quad (n \ge 1)$$

$\sum 3/n^{3/2}$ 是 $p = 3/2 > 1$ 的 $p$-级数，收敛。由比较判别法，$\sum |a_n b_n|$ 收敛，故原级数**绝对收敛**（从而收敛）。

**注**：本题说明并非所有乘积型级数都必须依赖 Abel 判别法——当 $b_n$ 本身可被 $p$-级数控制时，绝对收敛放缩法更简单直接。

</details>

---

**3.** 利用三角恒等式证明：对任意 $\theta \neq 2m\pi$（$m \in \mathbb{Z}$），有

(1) $\displaystyle\sum_{k=1}^n \sin(k\theta) = \frac{\cos\frac{\theta}{2} - \cos\left((n+\frac12)\theta\right)}{2\sin\frac{\theta}{2}}$

(2) $\displaystyle\sum_{k=1}^n \cos(k\theta) = \frac{\sin\left((n+\frac12)\theta\right) - \sin\frac{\theta}{2}}{2\sin\frac{\theta}{2}}$

并由此推导上界 $\left|\sum_{k=1}^n \sin(k\theta)\right| \le \frac{1}{|\sin(\theta/2)|}$，$\left|\sum_{k=1}^n \cos(k\theta)\right| \le \frac{1}{|\sin(\theta/2)|}$。

<details><summary>参考答案</summary>

**(1)** 利用积化和差：$2\sin\frac{\theta}{2}\sin(k\theta) = \cos\left((k-\frac12)\theta\right) - \cos\left((k+\frac12)\theta\right)$。

对 $k = 1$ 到 $n$ 求和，右端构成 telescoping sum：
$$\sum_{k=1}^n \left[\cos\left((k-\tfrac12)\theta\right) - \cos\left((k+\tfrac12)\theta\right)\right] = \cos\frac{\theta}{2} - \cos\left((n+\tfrac12)\theta\right)$$

两边除以 $2\sin\frac{\theta}{2}$ 即得公式。

**(2)** 利用积化和差：$2\sin\frac{\theta}{2}\cos(k\theta) = \sin\left((k+\frac12)\theta\right) - \sin\left((k-\frac12)\theta\right)$。

求和得：
$$\sum_{k=1}^n \left[\sin\left((k+\tfrac12)\theta\right) - \sin\left((k-\tfrac12)\theta\right)\right] = \sin\left((n+\tfrac12)\theta\right) - \sin\frac{\theta}{2}$$

两边除以 $2\sin\frac{\theta}{2}$ 即得公式。

**上界推导**：
$$\left|\sum_{k=1}^n \sin(k\theta)\right| = \frac{|\cos\frac{\theta}{2} - \cos((n+\frac12)\theta)|}{2|\sin\frac{\theta}{2}|} \le \frac{|\cos\frac{\theta}{2}| + |\cos((n+\frac12)\theta)|}{2|\sin\frac{\theta}{2}|} \le \frac{1+1}{2|\sin\frac{\theta}{2}|} = \frac{1}{|\sin\frac{\theta}{2}|}$$

同理 $\left|\sum_{k=1}^n \cos(k\theta)\right| \le \frac{1}{|\sin(\theta/2)|}$。

</details>

---

### 迁移应用

**4.** 判断级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin n}{n}$ 的敛散性。（提示：$\theta = 1$ 弧度，不是 $\pi$ 的有理数倍，因此不能通过有限周期直接计算——需要使用三角恒等式法。）

<details><summary>参考答案</summary>

令 $a_n = 1/n$，$b_n = \sin n$。使用 Dirichlet 判别法。

**条件 1**：$\{1/n\}$ 递减趋于 0。$\checkmark$

**条件 2**：需证明 $\sum_{k=1}^n \sin k$ 有界。

这里 $\theta = 1$（弧度），$\sin(1/2) \neq 0$。用三角恒等式法：

$$\left|\sum_{k=1}^n \sin k\right| = \left|\frac{\cos\frac12 - \cos\left(n+\frac12\right)}{2\sin\frac12}\right| \le \frac{1 + 1}{2\sin\frac12} = \frac{1}{\sin\frac12}$$

$\sin\frac12 \approx 0.479$，故上界约为 $2.086$。$\{B_n\}$ 有界。$\checkmark$

由 Dirichlet 判别法，$\displaystyle\sum_{n=1}^{\infty} \frac{\sin n}{n}$ **收敛**（条件收敛）。

**注**：本题中 $\theta = 1$ 不是 $\pi$ 的有理数倍，$\sin n$ 不呈现周期性，无法通过计算前若干项观察周期来证明有界性。正是三角恒等式法的通用性使之成为解决此类问题的标准工具。

</details>

---

**5.** 设 $p > 0$，讨论级数 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin n}{n^p}$ 的敛散性，并进行收敛类型分类。

<details><summary>参考答案</summary>

令 $a_n = 1/n^p$，$b_n = \sin n$。

**$p > 1$ 情形**：

$$\sum_{n=1}^{\infty} \left|\frac{\sin n}{n^p}\right| \le \sum_{n=1}^{\infty} \frac{1}{n^p}$$

右端是 $p > 1$ 的 $p$-级数，收敛。由比较判别法，原级数**绝对收敛**。

**$0 < p \le 1$ 情形**：

$\{a_n\} = \{1/n^p\}$ 递减趋于 0。$\sum_{k=1}^n \sin k$ 有界（如上题，上界 $1/\sin\frac12$）。由 Dirichlet 判别法，级数**收敛**。

考察绝对值级数 $\sum |\sin n|/n^p$。对 $n \ge 1$，$|\sin n| \ge \sin^2 n = \frac{1-\cos 2n}{2}$。

$$\sum \frac{|\sin n|}{n^p} \ge \frac12 \sum \frac{1}{n^p} - \frac12 \sum \frac{\cos 2n}{n^p}$$

第一项是 $p \le 1$ 的 $p$-级数，发散。第二项由 Dirichlet 判别法收敛（$1/n^p$ 递减趋于 0，$\sum \cos 2n$ 部分和有界）。发散级数减去收敛级数仍发散，故绝对值级数发散。

因此当 $0 < p \le 1$ 时，原级数**条件收敛**。

**$p \le 0$ 情形**：

通项 $\sin n / n^p$ 不趋于零（$\sin n$ 振荡不趋于零，$n^{-p} \ge 1$），由级数收敛的必要条件（定理 9.1 逆否命题），级数**发散**。

**完整分类**：

| $p$ 的范围 | 敛散性 | 收敛类型 |
|-----------|-------|---------|
| $p > 1$ | 收敛 | **绝对收敛** |
| $0 < p \le 1$ | 收敛 | **条件收敛** |
| $p \le 0$ | 发散 | — |

该分类与反常积分 $\int_1^\infty \frac{\sin x}{x^p}\,dx$ 的完整分类（ch8-04 例题）完全一致，体现了级数与积分在 Abel-Dirichlet 理论框架下的高度统一性。

</details>

---

**6.** 判断级数 $\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^n}{n}\left(1 + \frac{1}{\sqrt{n}}\right)$ 的敛散性并说明理由。

<details><summary>参考答案</summary>

方法一：注意到 $\frac{(-1)^n}{n}\left(1 + \frac{1}{\sqrt{n}}\right) = \frac{(-1)^n}{n} + \frac{(-1)^n}{n^{3/2}}$。两个级数分别收敛：
- $\sum (-1)^n/n$：交错调和级数，由 Leibniz 判别法收敛（条件收敛）
- $\sum (-1)^n/n^{3/2}$：绝对值级数 $\sum 1/n^{3/2}$ 收敛（$p$-级数，$p = 3/2 > 1$），故绝对收敛

两收敛级数之和收敛，因此原级数收敛。

方法二（Abel 判别法）：令 $a_n = 1 + 1/\sqrt{n}$，$b_n = (-1)^n/n$。
- $\{a_n\}$ 单调递减且有界（$1 < a_n \le 2$）$\checkmark$
- $\sum b_n = \sum (-1)^n/n$ 收敛（Leibniz）$\checkmark$
由 Abel 判别法，原级数收敛。

</details>
