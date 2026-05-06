# 07. 级数的 Cauchy 乘积

> 所属章节：第九章 数项级数  |  文件序号：07  |  难度：进阶
> 常见混淆点：Cauchy 乘积定义为 $c_n = \sum_{k=0}^n a_k b_{n-k}$ 而非 $c_n = a_n b_n$（初学者常误以为乘积的系数就是对应项直接相乘）；Mertens 定理要求至少一个级数绝对收敛——当两个级数都是条件收敛时，Cauchy 乘积可能发散（反例：$a_n = b_n = (-1)^n/\sqrt{n+1}$ 的 Cauchy 乘积通项不趋于零）

## 1. 学习目标与先修前置

### 学习目标
- 掌握两个级数 Cauchy 乘积的正式定义：$c_n = \sum_{k=0}^n a_k b_{n-k}$，理解其与多项式乘法（系数卷积）的类比
- 理解 Cauchy 乘积的几何意义——三角形部分和与矩形部分和的关系：$C_n = \sum_{k=0}^n a_k B_{n-k}$
- 掌握 Mertens 定理（定理 9.19）的完整陈述与证明，理解"至少一个级数绝对收敛"这一条件为何不可或缺
- 能验证给定级数对是否满足 Mertens 定理的条件，并确认 Cauchy 乘积收敛
- 能给出一个反例说明两个条件收敛级数的 Cauchy 乘积可能发散
- 熟悉常见的应用场景：指数函数级数乘积（验证 $e^x \cdot e^{-x} = 1$ 在级数乘法下的意义）

### 先修知识
- 文件 01（第九章）：级数的基本概念（定义 9.1、9.2）、部分和 $S_n$、级数收敛 $\iff$ $\lim S_n$ 存在且有限
- 文件 04（第九章）：绝对收敛的定义（定义 9.9）、条件收敛与绝对收敛的关系（定理 9.12）
- 文件 05（第九章）：Abel 变换（引理 9.13，形式一与形式三）——Mertens 证明的核心工具
- 文件 01（第八章）：无穷限反常积分的概念——用于类比理解卷积求和
- 二项式定理：$(x+y)^n = \sum_{k=0}^n \binom{n}{k} x^k y^{n-k}$——用于验证指数函数乘积的 Cauchy 系数
- 读者应熟悉双重求和的换序技巧以及数列极限的基本运算

---

## 2. 背景与应用场景

### 2.1 问题的提出

前六文构建了数项级数的完整判别体系：从比较判别法到比值/根值判别法、从交错级数到 Abel-Dirichlet 判别法、从积分判别法到绝对收敛理论。现在面对一个自然的问题：

**给定两个级数 $\sum a_n$ 和 $\sum b_n$（分别收敛于 $A$ 和 $B$），我们能否构造它们的"乘积"——即一个以某种方式组合 $a_n$ 和 $b_n$ 的级数——使之收敛于 $AB$？**

这个问题之所以重要，是因为在数学分析、复变函数、生成函数等领域，我们经常需要对级数进行乘法运算。例如：
- 幂级数的乘法：$(\sum a_n x^n)(\sum b_n x^n) = \sum c_n x^n$
- 指数函数的性质：$e^x \cdot e^y = e^{x+y}$ 的级数证明
- 生成函数中序列的卷积

### 2.2 从多项式乘法出发

两个多项式的乘法是最直观的出发点。设：

$$P(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_m x^m$$
$$Q(x) = b_0 + b_1 x + b_2 x^2 + \cdots + b_n x^n$$

乘积 $P(x)Q(x)$ 中 $x^k$ 的系数是：

$$c_k = \sum_{i+j=k} a_i b_j = a_0 b_k + a_1 b_{k-1} + \cdots + a_k b_0$$

这个系数 $c_k$ 由 $a_i$ 和 $b_j$ 的**所有满足 $i+j = k$ 的组合**相加得到——这正是离散卷积。

对于无穷级数，完全相同的思路给出了 **Cauchy 乘积** 的定义。

### 2.3 为什么需要绝对收敛条件？

Cauchy（柯西）本人曾认为：若 $\sum a_n$ 和 $\sum b_n$ 都收敛，则其 Cauchy 乘积也收敛且和为 $AB$。但这被后来的反例证伪。事实上，当两个级数都是**条件收敛**时，Cauchy 乘积可能发散。

一个经典反例是：

$$a_n = b_n = \frac{(-1)^n}{\sqrt{n+1}} \quad (n \ge 0)$$

两个级数均满足 Leibniz 判别法的条件（交替递减趋于零），故收敛（条件收敛）。但它们的 Cauchy 乘积通项：

$$c_n = \sum_{k=0}^n \frac{(-1)^k}{\sqrt{k+1}} \cdot \frac{(-1)^{n-k}}{\sqrt{n-k+1}} = (-1)^n \sum_{k=0}^n \frac{1}{\sqrt{(k+1)(n-k+1)}}$$

对每个 $k$，由 AM-GM 不等式：
$$\sqrt{(k+1)(n-k+1)} \le \frac{(k+1)+(n-k+1)}{2} = \frac{n+2}{2}$$

因此：
$$\frac{1}{\sqrt{(k+1)(n-k+1)}} \ge \frac{2}{n+2}$$

于是：
$$|c_n| \ge (n+1) \cdot \frac{2}{n+2} \to 2 \quad (n \to \infty)$$

通项 $c_n$ **不趋于零**，由必要条件（定理 9.1 逆否命题），Cauchy 乘积 $\sum c_n$ **发散**。

这一反例表明，要使 Cauchy 乘积收敛，**必须对级数的收敛性施加额外的条件**。Mertens 定理给出了一个漂亮的结果：只要其中一个级数**绝对收敛**，Cauchy 乘积就一定收敛于 $AB$。

---

## 3. 核心概念与符号约定

### 3.1 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $\sum_{n=0}^\infty a_n$，$\sum_{n=0}^\infty b_n$ | 两个级数（约定从 $n=0$ 开始） | 本章此前从 $n=1$ 开始，Cauchy 乘积从 $n=0$ 开始更自然 |
| $\{A_n\}$ | $A_n = \sum_{k=0}^n a_k$（$\{a_n\}$ 的部分和数列） | $A_n \to A$ |
| $\{B_n\}$ | $B_n = \sum_{k=0}^n b_k$（$\{b_n\}$ 的部分和数列） | $B_n \to B$ |
| $\{C_n\}$ | $C_n = \sum_{k=0}^n c_k$（Cauchy 乘积的部分和数列） | 希望 $C_n \to AB$ |
| $\alpha_n$ | $\alpha_n = A - A_n = \sum_{k=n+1}^\infty a_k$（$\{A_n\}$ 的余项） | $\alpha_n \to 0$ |
| $\beta_n$ | $\beta_n = B - B_n = \sum_{k=n+1}^\infty b_k$（$\{B_n\}$ 的余项） | $\beta_n \to 0$ |
| $c_n = \sum_{k=0}^n a_k b_{n-k}$ | Cauchy 乘积的通项（卷积形式） | 第 $n$ 项 |

### 3.2 关于指标约定的说明

本章此前各文件使用 $n=1$ 作为级数的起始指标。但在 Cauchy 乘积的定义和 Mertens 定理的证明中，从 $n=0$ 开始更为自然（与幂级数乘法、卷积的定义一致）。本章中约定：

- 所有级数从 $n=0$ 开始：$\sum_{n=0}^\infty a_n$，$\sum_{n=0}^\infty b_n$
- 部分和 $A_n = \sum_{k=0}^n a_k$，$B_n = \sum_{k=0}^n b_k$
- 若级数原本从 $n=1$ 开始，可设置 $a_0 = 0$ 做平凡延拓，不影响敛散性

---

## 4. 原理与方法

### 4.1 Cauchy 乘积的正式定义

**定义 9.18（Cauchy 乘积 / Cauchy Product）**：设 $\sum_{n=0}^\infty a_n$ 和 $\sum_{n=0}^\infty b_n$ 是两个级数。定义它们的 Cauchy 乘积为级数 $\sum_{n=0}^\infty c_n$，其中

$$\boxed{c_n = \sum_{k=0}^n a_k b_{n-k} = a_0 b_n + a_1 b_{n-1} + a_2 b_{n-2} + \cdots + a_n b_0} \tag{9.18.1}$$

即 $c_n$ 是 **$\{a_n\}$ 和 $\{b_n\}$ 的离散卷积**。

**与多项式乘法的类比**：形式幂级数 $(\sum a_n x^n)(\sum b_n x^n)$ 中 $x^n$ 的系数正是 $c_n$。如果将级数视为"形式幂级数在 $x=1$ 处的取值"，则 Cauchy 乘积自然对应于幂级数的乘法。

### 4.2 Cauchy 乘积的部分和表达式

为研究 Cauchy 乘积的收敛性，需找出其部分和 $C_n = \sum_{k=0}^n c_k$ 与 $\{a_n\}$、$\{b_n\}$ 的部分和之间的关系。

**命题 9.18.1（Cauchy 乘积的部分和）**：设 $C_n = \sum_{k=0}^n c_k$ 为 Cauchy 乘积的部分和，$A_n = \sum_{k=0}^n a_k$，$B_n = \sum_{k=0}^n b_k$。则

$$\boxed{C_n = \sum_{k=0}^n a_k B_{n-k} = \sum_{k=0}^n A_{n-k} b_k} \tag{9.18.2}$$

**证明**：

$$C_n = \sum_{k=0}^n c_k = \sum_{k=0}^n \sum_{i=0}^k a_i b_{k-i}$$

交换求和次序（先固定 $a_i$，再对 $k$ 从 $i$ 到 $n$ 求和）：

$$
\begin{aligned}
C_n &= \sum_{i=0}^n a_i \sum_{k=i}^n b_{k-i} \\
    &= \sum_{i=0}^n a_i \sum_{j=0}^{n-i} b_j \quad (\text{令 } j = k-i) \\
    &= \sum_{i=0}^n a_i B_{n-i}
\end{aligned}
$$

这给出了第一个表达式。对称地（先固定 $b_j$）：

$$
C_n = \sum_{j=0}^n b_j \sum_{k=j}^n a_{k-j} = \sum_{j=0}^n b_j \sum_{i=0}^{n-j} a_i = \sum_{j=0}^n A_{n-j} b_j
$$

证毕。

**几何直观**：$C_n = \sum_{i=0}^n a_i B_{n-i}$ 意味着 Cauchy 乘积的 $n$ 阶部分和由 $\{a_i\}$ 与 $\{B_j\}$ 的卷积给出。如果将 $A_n B_n = \sum_{i=0}^n \sum_{j=0}^n a_i b_j$ 视为 $i$-$j$ 平面上的 $(n+1)\times(n+1)$ 矩形区域，则 $C_n$ 对应于三角形区域 $\{i+j \le n\}$ 上的求和。两者的差（矩形减三角形）恰好是"交叉项"。

### 4.3 基本例子

**例 1（等比级数的 Cauchy 积）**：设 $a_n = b_n = \left(\frac12\right)^n$（$n \ge 0$），求 Cauchy 乘积 $\sum c_n$ 的通项并验证其收敛性。

**解**：

**第 1 步：计算 $c_n$**。

$$c_n = \sum_{k=0}^n \left(\frac12\right)^k \left(\frac12\right)^{n-k} = \sum_{k=0}^n \left(\frac12\right)^n = (n+1)\left(\frac12\right)^n$$

**第 2 步：求 $\sum a_n$ 和 $\sum b_n$ 的和**。

$\sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \left(\frac12\right)^n = \frac{1}{1-1/2} = 2$

同理 $\sum_{n=0}^\infty b_n = 2$，因此 $A = B = 2$。

**第 3 步：检验 Cauchy 乘积的收敛性**。

对 $\sum c_n = \sum (n+1)(1/2)^n$ 用比值判别法（定理 9.7）：

$$\lim_{n\to\infty} \frac{c_{n+1}}{c_n} = \lim_{n\to\infty} \frac{n+2}{n+1} \cdot \frac12 = \frac12 < 1$$

故级数 $\sum c_n$ 收敛。计算其和：

$$\sum_{n=0}^\infty (n+1)\left(\frac12\right)^n = \frac{1}{(1-1/2)^2} = 4$$

恰好等于 $AB = 2 \times 2 = 4$。在此例中，$\sum a_n$ 和 $\sum b_n$ 均为绝对收敛（等比级数 $|r| < 1$），Cauchy 乘积确实收敛于 $AB$。

**例 2（条件收敛的反例——回顾）**：

$a_n = b_n = \frac{(-1)^n}{\sqrt{n+1}}$（$n \ge 0$）。

- $\sum a_n$ 和 $\sum b_n$ 均为条件收敛（Leibniz 判别法，$a_n$ 递减趋于零）
- $c_n = (-1)^n \sum_{k=0}^n \frac{1}{\sqrt{(k+1)(n-k+1)}}$
- $|c_n| \ge (n+1) \cdot \frac{2}{n+2} \to 2 \neq 0$
- 因此 $\sum c_n$ 发散（通项不趋于零）

此例表明，当两个级数均仅为条件收敛时，Cauchy 乘积可能发散。这引出 Mertens 定理——只需一个级数绝对收敛，即可保证 Cauchy 乘积的收敛性。

---

### 4.4 Mertens 定理

#### 4.4.1 定理陈述

**定理 9.19（Mertens 定理 / Mertens' Theorem）**：设级数 $\sum_{n=0}^\infty a_n$ **绝对收敛**于 $A$，级数 $\sum_{n=0}^\infty b_n$ **收敛**于 $B$。则它们的 Cauchy 乘积 $\sum_{n=0}^\infty c_n$（其中 $c_n = \sum_{k=0}^n a_k b_{n-k}$）**收敛**，且

$$\boxed{\sum_{n=0}^\infty c_n = A \cdot B}$$

**注**：
- 定理条件是不对称的：**只要求其中一个级数绝对收敛**，另一个可以是条件收敛。这恰好覆盖了例 1（两个绝对收敛）和例 2 的反面情形（若将其中一个改为绝对收敛，则乘积必然收敛）。
- Mertens 定理是分析学中"绝对收敛级数可以像有限和那样运算"这一基本原则的体现。

#### 4.4.2 证明

**记号**：
- $A_n = \sum_{k=0}^n a_k$，$A = \lim A_n = \sum_{k=0}^\infty a_k$
- $B_n = \sum_{k=0}^n b_k$，$B = \lim B_n = \sum_{k=0}^\infty b_k$
- $C_n = \sum_{k=0}^n c_k$，目标：证明 $\lim C_n = AB$
- $\alpha_n = A - A_n = \sum_{k=n+1}^\infty a_k$（$\{A_n\}$ 的余项），$\alpha_n \to 0$
- $\beta_n = B - B_n = \sum_{k=n+1}^\infty b_k$（$\{B_n\}$ 的余项），$\beta_n \to 0$

**第 1 步：用 $\{a_k\}$ 和 $\{B_k\}$ 表示 $C_n$。**

由命题 9.18.1：
$$C_n = \sum_{k=0}^n a_k B_{n-k} \tag{4.4.1}$$

**第 2 步：将 $C_n - AB$ 分解为两部分。**

将 $B_{n-k} = B - \beta_{n-k}$ 代入 (4.4.1)：

$$
\begin{aligned}
C_n &= \sum_{k=0}^n a_k (B - \beta_{n-k}) \\
    &= B \sum_{k=0}^n a_k - \sum_{k=0}^n a_k \beta_{n-k} \\
    &= B A_n - \sum_{k=0}^n a_k \beta_{n-k}
\end{aligned}
$$

因此：
$$\boxed{C_n - AB = B(A_n - A) - \sum_{k=0}^n a_k \beta_{n-k} = -B \cdot \alpha_n - \sum_{k=0}^n a_k \beta_{n-k}} \tag{4.4.2}$$

**第 3 步：分析 $(4.4.2)$ 中各项的极限。**

第一项：由于 $\sum a_n$ 收敛（绝对收敛蕴含收敛），$A_n \to A$，故 $\alpha_n \to 0$，因此 $-B \cdot \alpha_n \to 0$。

第二项：令 $T_n = \sum_{k=0}^n a_k \beta_{n-k}$。我们需要证明 $T_n \to 0$。

**第 4 步：用 Abel 变换分析 $T_n$。**

由于 $\beta_k \to 0$（$B_k \to B$），对于任意给定的 $\varepsilon > 0$，存在 $M \in \mathbb{N}$，使得当 $k \ge M$ 时 $|\beta_k| < \varepsilon$。

将 $T_n$ 的求和指标反转：令 $j = n - k$，则当 $k$ 从 $0$ 到 $n$ 变化时，$j$ 从 $n$ 到 $0$ 变化：

$$T_n = \sum_{k=0}^n a_k \beta_{n-k} = \sum_{j=0}^n a_{n-j} \beta_j \tag{4.4.3}$$

现在对 $n > 2M$ 的情况进行分析。将 $(4.4.3)$ 中的求和拆分为两部分：

$$T_n = \underbrace{\sum_{j=0}^{M-1} a_{n-j} \beta_j}_{T_n^{(1)}} + \underbrace{\sum_{j=M}^n a_{n-j} \beta_j}_{T_n^{(2)}} \tag{4.4.4}$$

**对 $T_n^{(2)}$ 的估计**（$\beta_j$ 小的部分）：当 $j \ge M$ 时 $|\beta_j| < \varepsilon$，因此

$$|T_n^{(2)}| \le \sum_{j=M}^n |a_{n-j}| \cdot |\beta_j| \le \varepsilon \sum_{j=M}^n |a_{n-j}| = \varepsilon \sum_{k=0}^{n-M} |a_k| \le \varepsilon \sum_{k=0}^\infty |a_k| \tag{4.4.5}$$

这里 $k = n-j$，所以 $\sum_{j=M}^n |a_{n-j}| = \sum_{k=0}^{n-M} |a_k| \le \sum_{k=0}^\infty |a_k|$。记 $S_a = \sum_{k=0}^\infty |a_k|$（$\sum a_n$ 绝对收敛，故 $S_a$ 有限）。于是：

$$|T_n^{(2)}| \le \varepsilon \, S_a \tag{4.4.5'}$$

**对 $T_n^{(1)}$ 的估计**（指标 $j$ 固定为 $0$ 到 $M-1$）：这部分只有 $M$ 项。由于 $\sum a_n$ 收敛，通项 $a_n \to 0$。因此当 $n$ 充分大时，对 $j = 0, 1, \dots, M-1$ 有 $|a_{n-j}|$ 可以任意小。

具体地，$\{\beta_j\}$ 为收敛数列，必有界：设 $|\beta_j| \le R$ 对所有 $j$ 成立（例如取 $R = \sup |\beta_j|$ 或 $R = |\beta_0| + 1$ 等）。由于 $a_k \to 0$，存在 $N_1$ 使得对所有 $k \ge N_1$ 有 $|a_k| < \varepsilon / (MR)$。

取 $n \ge N_1 + M - 1$，则对 $j = 0, 1, \dots, M-1$，有 $n-j \ge N_1$，从而 $|a_{n-j}| < \varepsilon / (MR)$。因此：

$$|T_n^{(1)}| \le \sum_{j=0}^{M-1} |a_{n-j}| \cdot |\beta_j| \le \sum_{j=0}^{M-1} \frac{\varepsilon}{MR} \cdot R = M \cdot \frac{\varepsilon}{M} = \varepsilon \tag{4.4.6}$$

**第 5 步：综合估计。**

取 $N = \max\{N_1 + M - 1,\ 2M\}$。则对任意 $n \ge N$：

$$|T_n| \le |T_n^{(1)}| + |T_n^{(2)}| \le \varepsilon + \varepsilon \, S_a = (1 + S_a)\varepsilon \tag{4.4.7}$$

因此 $T_n \to 0$（因为 $\varepsilon$ 任意小）。

**第 6 步：回到 $(4.4.2)$ 完成证明。**

由 $(4.4.2)$：

$$|C_n - AB| \le |B| \cdot |\alpha_n| + |T_n|$$

由于 $\alpha_n \to 0$，存在 $N_2$ 使得 $n \ge N_2$ 时 $|\alpha_n| < \varepsilon/(|B|+1)$（避免 $B=0$ 时分母为零）。因此对 $n \ge \max\{N, N_2\}$：

$$|C_n - AB| \le |B| \cdot \frac{\varepsilon}{|B|+1} + (1+S_a)\varepsilon \le \varepsilon + (1+S_a)\varepsilon = (2 + S_a)\varepsilon$$

由于 $\varepsilon > 0$ 是任意的，$\lim_{n\to\infty} C_n = AB$。即 Cauchy 乘积收敛于 $AB$。证毕。

---

#### 4.4.3 证明逻辑总结

```
已知: ∑|a_n| 收敛 (∑ a_n = A 绝对收敛), ∑ b_n 收敛于 B

第一步: C_n = ∑_{k=0}^n a_k B_{n-k}           (部分和的卷积表达式)

第二步: C_n - AB = -B·α_n - ∑_{k=0}^n a_k β_{n-k}    (分解误差)

第三步: 分析 ∑_{k=0}^n a_k β_{n-k}
  ├── β_j → 0 ⇒ 存在 M 使 j ≥ M 时 |β_j| < ε
  ├── 将指标反转: ∑_{j=0}^n a_{n-j} β_j
  ├── 拆为 j<M 和 j≥M 两部分:
  │   ├── j ≥ M: |β_j| < ε, 用绝对有界控制: ≤ ε·∑|a_k|
  │   └── j < M: a_k → 0, 用通项衰减控制: ≤ ε
  └── 两部分均 < C·ε ⇒ T_n → 0

第四步: α_n → 0 ⇒ B·α_n → 0

结论: C_n → AB
```

#### 4.4.4 关于证明中使用的 Abel 变换

上述证明在第 4 步中将 $T_n = \sum_{k=0}^n a_k \beta_{n-k}$ 拆分为"$\beta_j$ 小的部分"和"$a_{n-j}$ 小的部分"，这种**分拆尾部求和指标**的技巧与引理 9.13（Abel 变换）的证明思路完全一致。Abel 变换的核心思想是：将 $\sum a_k b_k$ 转化为 $\sum (a_k - a_{k+1}) B_k$，从而利用 $a_k$ 的单调性和 $B_k$ 的有界性进行估计。在本证明中，虽然没有显式使用 Abel 变换的公式，但**分拆求和范围、分别利用不同条件控制不同部分**的方法论与 Abel 变换殊途同归。

Mertens 定理还有另一种直接使用 Abel 变换公式的证明路线，感兴趣的同学可以尝试将其作为拓展练习。

---

### 4.5 Mertens 定理的应用

#### 4.5.1 应用的一般步骤

使用 Mertens 定理判断 Cauchy 乘积收敛性的步骤：

1. **识别级数对**：确认 $\sum a_n$ 和 $\sum b_n$ 分别收敛于 $A$ 和 $B$
2. **验证绝对收敛条件**：检查哪个级数绝对收敛（至少一个必须绝对收敛）
3. **应用 Mertens 定理**：若满足条件，则 Cauchy 乘积 $\sum c_n$ 收敛于 $AB$
4. **计算乘积和**（如需）：利用已知的和 $A$、$B$ 相乘得到 $AB$

---

#### 4.5.2 应用一：幂级数与等比级数的乘积

**例 3**：设 $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{2^n}$（公比 $r=1/2$，绝对收敛于 $2$），$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{n+1}$（条件收敛于 $\ln 2$）。判断 Cauchy 乘积 $\sum c_n$ 的敛散性，若收敛求其和。

**解**：

**第 1 步：识别收敛类型**。

$\sum a_n = \sum 1/2^n$ 是等比级数，$|1/2| < 1$，**绝对收敛**于 $A = 2$。

$\sum b_n = \sum (-1)^n/(n+1)$ 是交错调和级数，由 Leibniz 判别法收敛（条件收敛）于 $B = \ln 2$。

**第 2 步：应用 Mertens 定理**。

$\sum a_n$ 绝对收敛，$\sum b_n$ 收敛，满足 Mertens 定理的条件。故 Cauchy 乘积 $\sum c_n$ 收敛于 $AB = 2 \ln 2$。

**第 3 步（选做）：写出 $c_n$ 的表达式**。

$$c_n = \sum_{k=0}^n \frac{1}{2^k} \cdot \frac{(-1)^{n-k}}{n-k+1}$$

此表达式没有简单的闭式形式，但由 Mertens 定理，我们知道 $\sum_{n=0}^\infty c_n = 2\ln 2$。

**对比**：若 $\sum a_n$ 和 $\sum b_n$ 均为条件收敛（如例 2 中的 $(-1)^n/\sqrt{n+1}$），Cauchy 乘积可能发散。Mertens 定理的绝对收敛条件恰好排除了这种危险情形。

---

#### 4.5.3 应用二：指数函数级数的乘积——一个优美的验证

**例 4**：设 $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{n!}$（绝对收敛于 $e$），$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{n!}$（绝对收敛于 $e^{-1}$）。求 Cauchy 乘积的通项 $c_n$ 并验证 $\sum c_n = 1$。

**解**：

**第 1 步：验证条件**。

$\sum 1/n!$ 和 $\sum (-1)^n/n!$ 均绝对收敛（比值判别法：$\lim (n!/(n+1)!) = 0 < 1$），因此 Mertens 定理的条件满足。

**第 2 步：计算 $c_n$**。

$$c_n = \sum_{k=0}^n \frac{1}{k!} \cdot \frac{(-1)^{n-k}}{(n-k)!}$$

提出 $1/n!$：

$$c_n = \frac{1}{n!} \sum_{k=0}^n \frac{n!}{k!(n-k)!} \cdot 1^k \cdot (-1)^{n-k}$$

应用二项式定理：$\sum_{k=0}^n \binom{n}{k} x^k y^{n-k} = (x+y)^n$。这里 $x = 1$，$y = -1$：

$$\sum_{k=0}^n \binom{n}{k} \cdot 1^k \cdot (-1)^{n-k} = (1-1)^n = 0^n$$

注意 $0^n$ 表示：
- 当 $n = 0$ 时，$0^0$ 理解为 $1$（二项式定理求和 $k=0$ 到 $0$，结果为 $1$）
- 当 $n \ge 1$ 时，$0^n = 0$

因此：
$$c_n = \begin{cases}
1, & n = 0 \\
0, & n \ge 1
\end{cases}$$

**第 3 步：验证 $\sum c_n = 1 = e \cdot e^{-1}$**。

$$\sum_{n=0}^\infty c_n = 1 + 0 + 0 + \cdots = 1$$

而 $AB = e \cdot e^{-1} = 1$，两者相等。Mertens 定理的结论成立。

**意义**：此例展示了 Mertens 定理在形式幂级数运算中的核心作用——它保证了指数函数的级数定义与指数运算律 $e^x \cdot e^{-x} = 1$ 的一致性。更一般地，绝对收敛幂级数的乘法对应于函数的逐点乘法，这正是生成函数理论和复分析的基础。

---

#### 4.5.4 应用三：$p$-级数与交错调和级数的乘积

**例 5**：设 $\displaystyle \sum_{n=1}^\infty a_n = \sum_{n=1}^\infty \frac{1}{n^2}$（$p=2 > 1$，绝对收敛于 $\pi^2/6$），$\displaystyle \sum_{n=1}^\infty b_n = \sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}$（条件收敛于 $\ln 2$）。将 $\{a_n\}$ 和 $\{b_n\}$ 延拓到 $n=0$ 令 $a_0 = b_0 = 0$，判断 Cauchy 乘积的敛散性。

**解**：

**第 1 步：延拓到 $n=0$ 起始**。

定义 $a_0 = 0$，对 $n \ge 1$ 有 $a_n = 1/n^2$；$b_0 = 0$，对 $n \ge 1$ 有 $b_n = (-1)^{n+1}/n$。增加首项 $0$ 不影响级数的和（收敛性不变）。

**第 2 步：验证条件**。

$\sum_{n=0}^\infty a_n = \sum_{n=1}^\infty 1/n^2$ 是 $p=2 > 1$ 的 $p$-级数，**绝对收敛**于 $A = \pi^2/6$（定理 9.5）。

$\sum_{n=0}^\infty b_n = \sum_{n=1}^\infty (-1)^{n+1}/n$ 是交错调和级数，由 Leibniz 判别法**条件收敛**于 $B = \ln 2$。

由 Mertens 定理（$\sum a_n$ 绝对收敛），Cauchy 乘积收敛于 $AB = (\pi^2/6) \ln 2$。

**第 3 步：写出 $c_n$ 的表达式**。

由于 $a_0 = b_0 = 0$，对 $n \ge 1$：

$$c_n = \sum_{k=0}^n a_k b_{n-k} = \sum_{k=1}^{n-1} \frac{1}{k^2} \cdot \frac{(-1)^{n-k+1}}{n-k}$$

注意当 $k=0$ 或 $k=n$ 时 $a_k b_{n-k} = 0$（因为 $a_0 = b_0 = 0$），因此求和实际上从 $k=1$ 到 $n-1$。

Mertens 定理保证了 $\sum_{n=0}^\infty c_n = (\pi^2/6) \ln 2$，尽管 $c_n$ 本身没有简单的闭式形式。

---

### 4.6 Mertens 定理的条件分析

Mertens 定理要求"至少一个级数绝对收敛"，这个条件是否可以进一步削弱？答案是否定的——例 2（两个条件收敛级数的 Cauchy 乘积发散）表明该条件是**最优的**。

更精确地说，对于 Cauchy 乘积的收敛性，存在以下谱系：

| $\sum a_n$ | $\sum b_n$ | Cauchy 乘积 $\sum c_n$ | 依据 |
|-----------|-----------|----------------------|------|
| 绝对收敛 | 绝对收敛 | 收敛（于 $AB$） | Mertens 定理 |
| 绝对收敛 | 条件收敛 | 收敛（于 $AB$） | Mertens 定理 |
| 绝对收敛 | 发散 | 不确定 | 无一般结论 |
| 条件收敛 | 条件收敛 | **可能发散** | 例 2 反例 |
| 条件收敛 | 条件收敛 | 也可能收敛 | 需更精细条件（如 Abel-Dirichlet 型） |

---

## 5. 例题

### 例 1：验证 Cauchy 乘积定义

给定 $a_n = \left(\frac13\right)^n$，$b_n = \left(-\frac12\right)^n$（$n \ge 0$）。求 Cauchy 乘积的通项 $c_n$，并判断 $\sum c_n$ 的收敛性。

**解**：

**第 1 步：写出 $c_n$ 的表达式**。

$$c_n = \sum_{k=0}^n \left(\frac13\right)^k \left(-\frac12\right)^{n-k} = \left(-\frac12\right)^n \sum_{k=0}^n \left(\frac13\right)^k \left(-\frac12\right)^{-k} = \left(-\frac12\right)^n \sum_{k=0}^n \left(-\frac{2}{3}\right)^k$$

**第 2 步：计算等比数例求和**。

$$\sum_{k=0}^n \left(-\frac{2}{3}\right)^k = \frac{1 - (-2/3)^{n+1}}{1 - (-2/3)} = \frac{1 - (-2/3)^{n+1}}{5/3} = \frac{3}{5}\left(1 - \left(-\frac{2}{3}\right)^{n+1}\right)$$

因此：
$$c_n = \left(-\frac12\right)^n \cdot \frac{3}{5}\left(1 - \left(-\frac{2}{3}\right)^{n+1}\right) = \frac{3}{5}\left[\left(-\frac12\right)^n - \left(-\frac12\right)^n \cdot \left(-\frac{2}{3}\right)^{n+1}\right]$$

化简第二项：
$$\left(-\frac12\right)^n \cdot \left(-\frac{2}{3}\right)^{n+1} = \left(-\frac12\right)^n \cdot \left(-\frac{2}{3}\right)^n \cdot \left(-\frac{2}{3}\right) = \left(\frac13\right)^n \cdot \left(-\frac{2}{3}\right)$$

因此：
$$c_n = \frac{3}{5}\left[\left(-\frac12\right)^n + \frac{2}{3}\left(\frac13\right)^n\right]$$

**第 3 步：判断 $\sum c_n$ 的收敛性**。

$\sum c_n$ 是两个等比级数的线性组合：
- $\sum (-\frac12)^n$：公比 $r = -\frac12$，$|r| = \frac12 < 1$，收敛
- $\sum (\frac13)^n$：公比 $r = \frac13$，$|r| = \frac13 < 1$，收敛

因此 $\sum c_n$ 收敛。计算其和：
$$\sum_{n=0}^\infty c_n = \frac{3}{5}\left(\frac{1}{1 - (-1/2)} + \frac{2}{3} \cdot \frac{1}{1 - 1/3}\right) = \frac{3}{5}\left(\frac{1}{3/2} + \frac{2}{3} \cdot \frac{1}{2/3}\right)$$

$$= \frac{3}{5}\left(\frac{2}{3} + \frac{2}{3} \cdot \frac{3}{2}\right) = \frac{3}{5}\left(\frac{2}{3} + 1\right) = \frac{3}{5} \cdot \frac{5}{3} = 1$$

**验证**：$\sum a_n = \frac{1}{1-1/3} = \frac{3}{2}$，$\sum b_n = \frac{1}{1-(-1/2)} = \frac{2}{3}$，$AB = \frac{3}{2} \cdot \frac{2}{3} = 1$。Mertens 定理成立。

---

### 例 2：识别 Mertens 定理的适用性

判断以下各对级数的 Cauchy 乘积是否收敛（可用 Mertens 定理的尽量用 Mertens 定理）：

(1) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{(n+1)^3}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{\sqrt[3]{n+1}}$

(2) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{(-1)^n}{\sqrt{n+1}}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{n+1}$

(3) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{(n+1)^2}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{1}{(n+1)^3}$

**解**：

**(1)** $\sum a_n = \sum 1/(n+1)^3$ 是 $p=3 > 1$ 的 $p$-级数，**绝对收敛**。

$\sum b_n = \sum (-1)^n/\sqrt[3]{n+1}$：检查 Leibniz 判别法的条件：
- 通项 $1/\sqrt[3]{n+1}$ 递减趋于 $0$，交错项 $(-1)^n$ → Leibniz 判别法适用
- 因此 $\sum b_n$ **条件收敛**

由 Mertens 定理（$\sum a_n$ 绝对收敛），Cauchy 乘积 **收敛**于 $AB$。

**(2)** $\sum a_n = \sum (-1)^n/\sqrt{n+1}$ 条件收敛（Leibniz 判别法，已见例 2）。

$\sum b_n = \sum (-1)^n/(n+1)$ 也是条件收敛（交错调和级数）。

两个级数均仅为条件收敛，Mertens 定理**不适用**。实际上，由例 2 的推理（$a_n$ 的量级与 $(-1)^n/\sqrt{n}$ 类似），可以猜想 Cauchy 乘积可能发散。事实上，更精细的分析（超出本文件范围）表明该 Cauchy 乘积确实发散。

**(3)** $\sum a_n = \sum 1/(n+1)^2$ 是 $p=2 > 1$ 的 $p$-级数，**绝对收敛**。

$\sum b_n = \sum 1/(n+1)^3$ 是 $p=3 > 1$ 的 $p$-级数，**绝对收敛**。

两个级数均绝对收敛，Mertens 定理适用。Cauchy 乘积 **收敛**于 $AB = (\pi^2/6) \cdot \zeta(3)$（其中 $\zeta(3)$ 是 Apéry 常数，约 $1.2020569$）。

---

### 例 3：从定义出发验证 Mertens 定理

取 $a_n = \frac{1}{2^n}$，$b_n = 1$（$n \ge 0$）。直接计算 Cauchy 乘积 $\sum c_n$ 并验证 $\sum c_n = (\sum a_n)(\sum b_n)$。

**解**：

**第 1 步：直接计算**。

$\sum_{n=0}^\infty a_n = \sum 1/2^n = 2$（绝对收敛），$\sum_{n=0}^\infty b_n = \sum 1$ **发散**（通项不趋于零）。

因此 Mertens 定理的条件不满足。但为理解 Cauchy 乘积的构造，仍可计算 $c_n$：

$$c_n = \sum_{k=0}^n \frac{1}{2^k} \cdot 1 = \sum_{k=0}^n \frac{1}{2^k} = 2\left(1 - \frac{1}{2^{n+1}}\right)$$

部分和 $C_n = \sum_{k=0}^n c_k = \sum_{k=0}^n 2\left(1 - \frac{1}{2^{k+1}}\right) = 2(n+1) - 2\sum_{k=0}^n \frac{1}{2^{k+1}} = 2(n+1) - 2\left(1 - \frac{1}{2^{n+1}}\right)$。

当 $n \to \infty$ 时 $C_n \to \infty$，Cauchy 乘积发散。这符合预期——即使形式上计算 $\sum b_n$ 发散，Cauchy 乘积也无意义。

**第 2 步：修正——选择一个收敛的 $\sum b_n$**。

取 $b_n = \left(\frac13\right)^n$，则 $\sum a_n = 2$，$\sum b_n = \frac{1}{1-1/3} = \frac{3}{2}$，均绝对收敛。

计算 $c_n$ 用与例题 1 类似的方法（留作练习），可得 $\sum c_n = 3 = 2 \times \frac{3}{2}$。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| Cauchy 乘积的通项是 $c_n = a_n b_n$（对应项直接相乘） | $c_n = \sum_{k=0}^n a_k b_{n-k}$，是离散卷积而非逐项乘积。直接相乘会丢失 $i+j=n$ 的所有交叉项 |
| 若 $\sum a_n$ 和 $\sum b_n$ 都收敛，则 Cauchy 乘积必收敛 | 这是 Cauchy 最初犯的错误。两个条件收敛级数的 Cauchy 乘积可能发散（例：$a_n = b_n = (-1)^n/\sqrt{n+1}$） |
| 只需 $\sum a_n$ 绝对收敛，无论 $\sum b_n$ 是否收敛，Cauchy 乘积必收敛 | Mertens 定理要求 $\sum b_n$ 也收敛（至少条件收敛）。如果 $\sum b_n$ 发散，Cauchy 乘积不一定收敛 |
| Mertens 定理的证明复杂且无法理解 | 证明的核心思想很直观：将误差 $\sum a_k \beta_{n-k}$ 拆成两部分——"$\beta$ 小"的部分用绝对收敛控制，"$a$ 小"的部分用通项趋于零控制 |
| Cauchy 乘积的部分和 $C_n$ 与 $A_n B_n$ 相同 | $A_n B_n = \sum_{i=0}^n \sum_{j=0}^n a_i b_j$（矩形区域），而 $C_n = \sum_{i=0}^n \sum_{j=0}^{n-i} a_i b_j$（三角形区域）。两者之差恰好是"交叉项" $i+j > n$ 的部分 |
| $c_n = \sum_{k=0}^n a_k b_{n-k}$ 的求和范围是 $k=0$ 到 $n$，每次都要算 $n+1$ 项，太复杂 | 确实，$c_n$ 的计算量随 $n$ 增大而增大。但 Mertens 定理的好处在于：**我们不需要直接计算 $c_n$ 的具体值**——定理保证了 Cauchy 乘积收敛于 $AB$ |

### 检查点

- [ ] 能否写出 $\sum a_n$ 和 $\sum b_n$ 的 Cauchy 乘积的通项 $c_n$ 的表达式？
- [ ] 能否解释为什么 $c_n$ 的定义与多项式乘法的系数相一致？
- [ ] 能否通过交换求和次序证明 $C_n = \sum_{k=0}^n a_k B_{n-k}$？
- [ ] 能否完整叙述 Mertens 定理（定理 9.19）的条件和结论？
- [ ] 能否给出一个反例说明两个条件收敛级数的 Cauchy 乘积可能发散？
- [ ] 在 Mertens 定理的证明中，能否解释 $T_n = \sum a_k \beta_{n-k}$ 中的分拆技巧？
- [ ] 能否应用 Mertens 定理判断 $\sum 1/n^2 \times \sum (-1)^n/n$ 的 Cauchy 乘积收敛性？
- [ ] 能否用二项式定理验证 $e^x$ 和 $e^{-x}$ 的 Cauchy 乘积为常数 $1$？
- [ ] 能否指出 Mertens 定理中绝对收敛条件的作用——它为什么不可或缺？

---

## 练习题

### 基础巩固

**1.** 设 $a_n = \left(\frac14\right)^n$，$b_n = \left(\frac12\right)^n$（$n \ge 0$）。求 Cauchy 乘积的通项 $c_n$，并验证 $\sum c_n = (\sum a_n)(\sum b_n)$。

<details><summary>参考答案</summary>

**第 1 步：计算 $c_n$**。

$$c_n = \sum_{k=0}^n \left(\frac14\right)^k \left(\frac12\right)^{n-k} = \left(\frac12\right)^n \sum_{k=0}^n \left(\frac14\right)^k \left(\frac12\right)^{-k} = \left(\frac12\right)^n \sum_{k=0}^n \left(\frac12\right)^k$$

$$= \left(\frac12\right)^n \cdot \frac{1 - (1/2)^{n+1}}{1 - 1/2} = \left(\frac12\right)^n \cdot 2\left(1 - \frac{1}{2^{n+1}}\right) = 2\left(\frac12\right)^n - 2\left(\frac12\right)^{2n+1}$$

$$= 2\left(\frac12\right)^n - \left(\frac14\right)^n$$

**第 2 步：验证 $\sum c_n$ 收敛**。

$\sum c_n = 2\sum (1/2)^n - \sum (1/4)^n$。两个等比级数均收敛（$|r| < 1$），故 $\sum c_n$ 收敛。

**第 3 步：计算和**。

$$\sum_{n=0}^\infty c_n = 2 \cdot \frac{1}{1-1/2} - \frac{1}{1-1/4} = 2 \cdot 2 - \frac{1}{3/4} = 4 - \frac{4}{3} = \frac{8}{3}$$

**验证**：$\sum a_n = \frac{1}{1-1/4} = \frac{4}{3}$，$\sum b_n = \frac{1}{1-1/2} = 2$，$AB = \frac{4}{3} \cdot 2 = \frac{8}{3}$。结果一致。

</details>

---

**2.** 判断以下各对级数的 Cauchy 乘积是否收敛（可用 Mertens 定理的尽量用）：

(1) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{(n+1)^{3/2}}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{n+1}$

(2) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{1}{3^n}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{3^n}$

(3) $\displaystyle \sum_{n=0}^\infty a_n = \sum_{n=0}^\infty \frac{(-1)^n}{\ln(n+2)}$，$\displaystyle \sum_{n=0}^\infty b_n = \sum_{n=0}^\infty \frac{(-1)^n}{n+1}$

<details><summary>参考答案</summary>

**(1)** $\sum a_n = \sum 1/(n+1)^{3/2}$ 是 $p = 3/2 > 1$ 的 $p$-级数，**绝对收敛**。

$\sum b_n = \sum (-1)^n/(n+1)$ 是交错调和级数，由 Leibniz 判别法知**条件收敛**。

由 Mertens 定理（$\sum a_n$ 绝对收敛），Cauchy 乘积**收敛**于 $AB$。其中 $A = \zeta(3/2) = \sum_{n=1}^\infty 1/n^{3/2}$（约 $2.612$），$B = \ln 2$。

**(2)** $\sum a_n = \sum (1/3)^n$ 是等比级数，$|r| = 1/3 < 1$，**绝对收敛**于 $A = \frac{1}{1-1/3} = \frac{3}{2}$。

$\sum b_n = \sum (-1/3)^n$ 也是等比级数，$|r| = 1/3 < 1$，**绝对收敛**于 $B = \frac{1}{1-(-1/3)} = \frac{3}{4}$。

由 Mertens 定理，Cauchy 乘积**收敛**于 $AB = \frac{3}{2} \cdot \frac{3}{4} = \frac{9}{8}$。

（实际上，不借助 Mertens 定理也可计算：$c_n = \sum_{k=0}^n (1/3)^k(-1/3)^{n-k} = (n+1)(-1)^n(1/3)^n$，再用等比级数求和可得 $\sum c_n = 9/8$。）

**(3)** $\sum a_n = \sum (-1)^n/\ln(n+2)$：通项趋于零（$\ln(n+2) \to \infty$），且递减（$\ln$ 函数递增，倒数递减），由 Leibniz 判别法知**条件收敛**。

$\sum b_n = \sum (-1)^n/(n+1)$：**条件收敛**（交错调和级数）。

两个级数均仅为条件收敛，Mertens 定理**不适用**。事实上，通过对通项量级的分析（$1/\ln(n+2)$ 的衰减速度慢于 $1/\sqrt{n+1}$），可以证明 Cauchy 乘积发散。此例说明：即使两个条件收敛级数的通项都趋于零，Cauchy 乘积仍然可能发散。

</details>

---

**3.** 设 $a_n = 1$ 对一切 $n \ge 0$，$b_n = \frac{(-1)^n}{2^n}$ 对一切 $n \ge 0$。

(1) $\sum a_n$ 是否收敛？
(2) 能否应用 Mertens 定理？
(3) 直接计算 Cauchy 乘积并判断 $\sum c_n$ 的敛散性。

<details><summary>参考答案</summary>

**(1)** $\sum_{n=0}^\infty a_n = \sum 1$ 的部分和 $S_n = n+1 \to \infty$，**发散**。不满足 Mertens 定理的条件（两个级数都必须收敛）。

**(2)** 不能。Mertens 定理要求 $\sum a_n$ 和 $\sum b_n$ 都收敛（至少条件收敛），但 $\sum a_n$ 发散。

**(3)** 直接计算 $c_n$：

$$c_n = \sum_{k=0}^n 1 \cdot \frac{(-1)^{n-k}}{2^{n-k}} = \sum_{k=0}^n \left(-\frac12\right)^{n-k} = \sum_{j=0}^n \left(-\frac12\right)^j = \frac{1 - (-1/2)^{n+1}}{1 - (-1/2)} = \frac{2}{3}\left(1 - \left(-\frac12\right)^{n+1}\right)$$

部分和：
$$C_n = \sum_{k=0}^n c_k = \frac{2}{3}\sum_{k=0}^n \left(1 - \left(-\frac12\right)^{k+1}\right) = \frac{2}{3}\left(n+1 - \sum_{k=0}^n \left(-\frac12\right)^{k+1}\right)$$

$$\sum_{k=0}^n \left(-\frac12\right)^{k+1} = \left(-\frac12\right)\frac{1-(-1/2)^{n+1}}{1-(-1/2)} = -\frac12 \cdot \frac{2}{3}\left(1-(-1/2)^{n+1}\right) = -\frac13\left(1-(-1/2)^{n+1}\right)$$

因此：
$$C_n = \frac{2}{3}\left(n+1 + \frac13\left(1-(-1/2)^{n+1}\right)\right) \to \infty \quad (n \to \infty)$$

$\sum c_n$ **发散**（部分和趋于无穷）。这说明当其中一个级数发散时，Cauchy 乘积也可能发散。

</details>

---

### 迁移应用

**4.** 证明：若 $\sum a_n$ 和 $\sum b_n$ 均绝对收敛，则 $\sum c_n$（Cauchy 乘积）也绝对收敛，且 $\sum |c_n| \le (\sum |a_n|)(\sum |b_n|)$。

<details><summary>参考答案</summary>

**证明**：

**第 1 步：写出 $c_n$ 的表达式**。
$$c_n = \sum_{k=0}^n a_k b_{n-k}$$

**第 2 步：对 $|c_n|$ 用三角不等式**。
$$|c_n| \le \sum_{k=0}^n |a_k| \cdot |b_{n-k}|$$

**第 3 步：考虑 $\sum_{n=0}^\infty |c_n|$ 的部分和**。

对部分和 $S_N = \sum_{n=0}^N |c_n|$：
$$S_N \le \sum_{n=0}^N \sum_{k=0}^n |a_k| \cdot |b_{n-k}|$$

交换求和次序（先对 $k$ 求和）：
$$S_N \le \sum_{k=0}^N |a_k| \sum_{n=k}^N |b_{n-k}| = \sum_{k=0}^N |a_k| \sum_{j=0}^{N-k} |b_j| \le \sum_{k=0}^N |a_k| \sum_{j=0}^\infty |b_j| = \left(\sum_{j=0}^\infty |b_j|\right) \sum_{k=0}^N |a_k|$$

**第 4 步：取极限**。

$\sum_{k=0}^N |a_k| \le \sum_{k=0}^\infty |a_k|$，所以：
$$S_N \le \left(\sum_{j=0}^\infty |b_j|\right)\left(\sum_{k=0}^\infty |a_k|\right)$$

$S_N$ 单调递增且有上界，故 $\lim_{N\to\infty} S_N$ 存在且有限，即 $\sum |c_n|$ 收敛。且：
$$\sum_{n=0}^\infty |c_n| \le \left(\sum_{n=0}^\infty |a_n|\right)\left(\sum_{n=0}^\infty |b_n|\right)$$

证毕。

**注**：这个不等式表明，当两个级数都绝对收敛时，Cauchy 乘积不仅收敛，而且绝对收敛——这比 Mertens 定理的结论更强（Mertens 只要求一个绝对收敛，结论为普通收敛）。

</details>

---

**5.** （Cauchy 乘积的结合律）设 $\sum a_n$、$\sum b_n$、$\sum c_n$ 均绝对收敛。证明：$(\sum a_n \star \sum b_n) \star \sum c_n = \sum a_n \star (\sum b_n \star \sum c_n)$，其中 $\star$ 表示 Cauchy 乘积运算。

<details><summary>参考答案</summary>

**分析**：Cauchy 乘积的结合律是指，先计算 $\sum a_n$ 和 $\sum b_n$ 的 Cauchy 乘积得到 $\sum d_n$（其中 $d_n = \sum_{k=0}^n a_k b_{n-k}$），再与 $\sum c_n$ 做 Cauchy 乘积得到 $\sum e_n$，与先计算 $\sum b_n$ 和 $\sum c_n$ 的 Cauchy 乘积再与 $\sum a_n$ 计算的结果相同。

**证明**：

**第 1 步：定义两个 Cauchy 乘积**。

设 $d_n = \sum_{i=0}^n a_i b_{n-i}$（$\sum a$ 与 $\sum b$ 的 Cauchy 积通项）。

则 $(\sum a \star \sum b) \star \sum c$ 的通项为：
$$e_n = \sum_{j=0}^n d_j c_{n-j} = \sum_{j=0}^n \left(\sum_{i=0}^j a_i b_{j-i}\right) c_{n-j}$$

交换求和次序（先对 $i$ 求和，再对 $j \ge i$ 求和）：
$$e_n = \sum_{i=0}^n a_i \sum_{j=i}^n b_{j-i} c_{n-j} = \sum_{i=0}^n a_i \sum_{k=0}^{n-i} b_k c_{n-i-k} \quad (\text{令 } k = j-i)$$

**第 2 步：定义另一种结合方式**。

设 $f_n = \sum_{k=0}^n b_k c_{n-k}$（$\sum b$ 与 $\sum c$ 的 Cauchy 积通项）。

则 $\sum a \star (\sum b \star \sum c)$ 的通项为：
$$g_n = \sum_{i=0}^n a_i f_{n-i} = \sum_{i=0}^n a_i \left(\sum_{k=0}^{n-i} b_k c_{n-i-k}\right)$$

**第 3 步：比较 $e_n$ 与 $g_n$**。

$$e_n = \sum_{i=0}^n a_i \sum_{k=0}^{n-i} b_k c_{n-i-k} = g_n$$

两者完全相等！因此 Cauchy 乘积的结合律对**通项成立**（而不需要收敛性条件）。

**第 4 步：收敛性**。

由练习题 4 的结论，两个绝对收敛级数的 Cauchy 乘积绝对收敛。因此：
- $\sum d_n = \sum a \star \sum b$ 绝对收敛
- $(\sum a \star \sum b) \star \sum c$ 绝对收敛（从而普通收敛）
- 同理 $\sum a \star (\sum b \star \sum c)$ 绝对收敛

因此结合律在收敛和的意义上也成立。

**注**：结合律的证明中，通项 $e_n = g_n$ 的等式是**纯代数**的（不依赖任何收敛性条件），但等式两端的级数 $\sum e_n$ 和 $\sum g_n$ 的收敛性则依赖于 $\sum a$、$\sum b$、$\sum c$ 的绝对收敛性（由练习题 4 保证）。

</details>

---

**6.** （思考题）设 $a_n = \frac{(-1)^n}{n+1}$，$b_n = \frac{(-1)^n}{(n+1)^2}$。$\sum a_n$ 条件收敛，$\sum b_n$ 绝对收敛。由 Mertens 定理，Cauchy 乘积收敛。试计算 $c_0, c_1, c_2$ 并观察其数值。

<details><summary>参考答案</summary>

**第 1 步：写出通项**。

$$c_n = \sum_{k=0}^n \frac{(-1)^k}{k+1} \cdot \frac{(-1)^{n-k}}{(n-k+1)^2} = (-1)^n \sum_{k=0}^n \frac{1}{(k+1)(n-k+1)^2}$$

**第 2 步：计算 $c_0, c_1, c_2$**。

$n=0$：
$$c_0 = (-1)^0 \cdot \frac{1}{1 \cdot 1^2} = 1$$

$n=1$：
$$c_1 = (-1)^1 \left(\frac{1}{1 \cdot 2^2} + \frac{1}{2 \cdot 1^2}\right) = -\left(\frac14 + \frac12\right) = -\frac34$$

$n=2$：
$$c_2 = (-1)^2 \left(\frac{1}{1 \cdot 3^2} + \frac{1}{2 \cdot 2^2} + \frac{1}{3 \cdot 1^2}\right) = \frac19 + \frac18 + \frac13$$

通分 $72$：
$$\frac{8}{72} + \frac{9}{72} + \frac{24}{72} = \frac{41}{72}$$

因此 $c_0 = 1$，$c_1 = -\frac34 = -0.75$，$c_2 = \frac{41}{72} \approx 0.5694$。

**第 3 步：观察**。

$c_n$ 的绝对值似乎在减小（$1, 0.75, 0.5694$），但 $c_n$ 的正负交替（由 $(-1)^n$ 决定）。Mertens 定理保证 $\sum_{n=0}^\infty c_n$ 收敛于 $AB = (\ln 2) \cdot (\pi^2/6) \approx 0.6931 \times 1.6449 \approx 1.140$。

</details>

