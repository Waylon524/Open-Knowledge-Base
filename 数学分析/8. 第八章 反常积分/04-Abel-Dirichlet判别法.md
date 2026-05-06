# 04. Abel-Dirichlet 判别法

> 所属章节：第八章 反常积分  |  文件序号：04  |  难度：进阶
> 常见混淆点：Dirichlet 判别法要求 f **单调趋于零** + ∫g 的部分积分**一致有界**；Abel 判别法要求 f **单调有界** + ∫g **收敛**。两者条件恰好"互换"——不要用混；条件收敛证明中 |sin x| 放缩的**方向是向下而非向上**（要找下界而非上界），利用 |sin x| ≥ sin²x 将绝对值积分拆为发散 + 收敛两部分

## 1. 学习目标与先修前置

### 学习目标
- 掌握积分第二中值定理（Bonnet 形式），理解其作为 Dirichlet/Abel 判别法证明核心工具的角色
- 掌握 Dirichlet 判别法的定理陈述（无穷限 + 瑕积分两种情形），能独立完成其 Cauchy 准则证明
- 掌握 Abel 判别法的定理陈述（无穷限 + 瑕积分两种情形），能通过归化为 Dirichlet 判别法完成证明
- 能熟练运用 Dirichlet 判别法验证 ∫ sin x/x dx 等经典条件收敛积分的收敛性
- 能完成含参积分 I(p) = ∫_1^∞ sin x/x^p dx 的完整分类讨论（p>1 绝对收敛；0<p≤1 条件收敛；p≤0 发散）
- 掌握条件收敛证明的 sin²x 拆分技术：|sin x|/x^p ≥ sin²x/x^p = (1−cos 2x)/(2x^p)，拆为发散项 + 收敛项
- 能运用 Abel 判别法验证 f 单调有界 + ∫g 收敛 → ∫fg 收敛的完整流程

### 先修知识
- 文件 03（第八章）：Cauchy 收敛准则（定理 8.6/8.7）、绝对收敛与条件收敛的定义（定义 8.3/8.4）、绝对收敛 ⇒ 收敛定理（定理 8.8）
- 文件 02（第八章）：比较判别法（定理 8.3/8.4）——用于判定 p>1 情形的绝对收敛
- 文件 01（第八章）：p-积分的敛散性结论（定理 8.1/8.2）
- 三角不等式与三角恒等式：|sin x| ≥ sin²x，sin²x = (1 − cos 2x)/2
- 极限的有界性：若 lim_{x→∞} F(x) 存在有限，则 F(x) 有界（常用于证明 ∫g 收敛 ⇒ ∫g 的部分积分有界）

---

## 2. 背景与应用场景

### 2.1 比较判别法的"天花板"

ch8-02 的比较判别法和 ch8-03 的绝对收敛判定法解决了一大类问题——只要能够找到合适的非负比较函数，就能判断敛散性。然而，这类方法有一个根本局限：

> 比较判别法要求被积函数**非负**，所以只能处理 **∫|f| 的收敛性**。由定理 8.8，如果 ∫|f| 收敛（绝对收敛），则 ∫f 收敛。但**如果 ∫|f| 发散**，比较判别法对 ∫f 本身的收敛性无法给出任何信息。

考虑下面的积分：

$$\int_1^{+\infty} \frac{\sin x}{x}\,dx$$

这个积分中的 sin x 正负交替，被积函数既不是非负的，也不绝对收敛（|sin x|/x ~ 1/x 发散）。但它在某种意义下"应该"收敛——因为正负部分互相抵消，使得净面积有限。这个直觉是否正确？如果是，如何证明？

### 2.2 核心思想：利用乘积结构

许多涉及变号函数的反常积分具有乘积形式：

$$\int_a^{+\infty} f(x)g(x)\,dx$$

其中：
- $f(x)$ 是单调的（控制"衰减"或"增长"）
- $g(x)$ 是振荡的（提供正负交替）

Dirichlet 判别法和 Abel 判别法专门针对这种"单调因子 × 振荡因子"的结构。它们的核心思想是：**即使 |f(x)g(x)| 的积分发散（非绝对收敛），正负交替的振荡也能使净积分收敛**——条件是振荡的幅度（由 f 控制）最终趋于零。

### 2.3 知识图谱

```
ch8-01: 反常积分的定义（极限语言）
ch8-02: 比较判别法（非负函数）
ch8-03: Cauchy 准则 + 绝对/条件收敛的定义
       ↓
ch8-04: 积分第二中值定理（关键引理）
       → Dirichlet 判别法（f→0, ∫g 有界 ⇒ ∫fg 收敛）
       → Abel 判别法（f 有界, ∫g 收敛 ⇒ ∫fg 收敛）
       → 条件收敛的证明技术（sin²x 拆分）
       → 含参积分 I(p) 的完整讨论
```

---

## 3. 核心概念与符号约定

### 3.1 两类判别法的条件对比

| 对比维度 | Dirichlet 判别法 | Abel 判别法 |
|----------|-----------------|-------------|
| 对 $f$ 的要求 | 单调，$\lim_{x\to+\infty} f(x)=0$ | 单调，**有界**（不必趋于 0） |
| 对 $\int g$ 的要求 | $\int_a^A g(x)\,dx$ **一致有界** | $\int_a^{+\infty} g(x)\,dx$ **收敛** |
| 适用范围 | $f$ 衰减 + $g$ 振荡 | $f$ 有界变化 + $g$ 收敛 |
| 证明工具 | 积分第二中值定理 + Cauchy 准则 | 归化为 Dirichlet 判别法 |

### 3.2 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $f(x)$ | "衰减因子"或"有界因子" | 单调函数 |
| $g(x)$ | "振荡因子"或"被乘的积分" | 与 f 相乘构成被积函数 |
| $G(u) = \int_a^u g(x)\,dx$ | $g$ 的变上限积分 | 在 Dirichlet 条件中一致有界 |
| $\xi$ | 积分第二中值定理中的中间点 | $\xi \in [A_1, A_2]$ |
| $M$ | $|G(u)| \le M$ 中的上界 | Dirichlet 条件中的有界常数 |
| $K$ | $|f(x)| \le K$ 中的上界 | Abel 条件中的有界常数 |
| $L = \lim_{x\to\infty} f(x)$ | 单调有界函数 $f$ 的极限 | 用于将 Abel 归化为 Dirichlet |
| $\sin^2 x = (1-\cos 2x)/2$ | 三角恒等式 | 条件收敛证明的核心技巧 |

---

## 4. 原理与方法

### 4.1 积分第二中值定理（第二积分中值定理）

证明 Dirichlet 判别法和 Abel 判别法需要借助一个重要的引理——积分第二中值定理（Second Mean Value Theorem for Integrals, Bonnet 形式）。

**引理 8.1（积分第二中值定理 / Bonnet 形式）**：设 $f$ 在 $[a, b]$ 上单调，$g$ 在 $[a, b]$ 上 Riemann 可积。则存在 $\xi \in [a, b]$，使得：

$$\int_a^b f(x)g(x)\,dx = f(a)\int_a^{\xi} g(x)\,dx + f(b)\int_{\xi}^b g(x)\,dx$$

**特别地**：
- 若 $f$ 在 $[a, b]$ 上非负单调递减，则 $\exists \xi \in [a, b]$：
  $$\int_a^b f(x)g(x)\,dx = f(a)\int_a^{\xi} g(x)\,dx$$
- 若 $f$ 在 $[a, b]$ 上非负单调递增，则 $\exists \xi \in [a, b]$：
  $$\int_a^b f(x)g(x)\,dx = f(b)\int_{\xi}^b g(x)\,dx$$

**证明思路**：利用 Riemann 积分定义和 Abel 变换（分部求和公式）可证。此处从略，作为已知引理直接使用。

**直观理解**：该定理说，当 $f$ 单调时，$\int_a^b fg$ 可以被简化为 $f$ 在端点处的值乘以 $g$ 在某个子区间上的积分。这个"简化"是将 Dirichlet 判别法的证明从 $\int fg$ 的 Cauchy 条件转化为 $f$ 的**端点值**控制的关键。

---

### 4.2 Dirichlet 判别法

#### 4.2.1 无穷限情形的定理陈述与证明

**定理 8.9（Dirichlet 判别法——无穷限反常积分）**：设 $f(x), g(x)$ 在 $[a, +\infty)$ 上有定义，且对任意 $b > a$，$f$ 和 $g$ 在 $[a, b]$ 上 Riemann 可积。若：

1. **$f(x)$ 单调趋于零**：$f$ 在 $[a, +\infty)$ 上单调，且 $\displaystyle\lim_{x\to +\infty} f(x) = 0$；
2. **$\int_a^A g(x)\,dx$ 一致有界**：存在 $M > 0$，使得对一切 $A \ge a$：
   $$\left|\int_a^A g(x)\,dx\right| \le M$$

则反常积分 $\displaystyle\int_a^{+\infty} f(x)g(x)\,dx$ **收敛**。

**证明**：

我们使用 Cauchy 收敛准则（定理 8.6）来完成证明。对任意 $\varepsilon > 0$，需要找到 $X \ge a$，使得对所有 $A_2 > A_1 \ge X$ 有：

$$\left|\int_{A_1}^{A_2} f(x)g(x)\,dx\right| < \varepsilon$$

---

**第 1 步：利用 $f\to 0$ 控制端点值**。

由于 $\lim_{x\to +\infty} f(x) = 0$，对给定的 $\varepsilon > 0$，存在 $X \ge a$，使得对所有 $x \ge X$：

$$|f(x)| < \frac{\varepsilon}{4M}$$

---

**第 2 步：在 $[A_1, A_2]$ 上应用积分第二中值定理**。

取任意 $A_2 > A_1 \ge X$。由于 $f$ 在 $[A_1, A_2]$ 上单调，由引理 8.1（Bonnet 形式），存在 $\xi \in [A_1, A_2]$，使得：

$$\int_{A_1}^{A_2} f(x)g(x)\,dx = f(A_1)\int_{A_1}^{\xi} g(x)\,dx + f(A_2)\int_{\xi}^{A_2} g(x)\,dx$$

---

**第 3 步：用 $G(u) = \int_a^u g$ 的有界性控制中间积分**。

记 $G(u) = \displaystyle\int_a^u g(x)\,dx$。由条件 2，$|G(u)| \le M$ 对一切 $u \ge a$ 成立。

于是：

$$\left|\int_{A_1}^{\xi} g(x)\,dx\right| = |G(\xi) - G(A_1)| \le |G(\xi)| + |G(A_1)| \le 2M$$

$$\left|\int_{\xi}^{A_2} g(x)\,dx\right| = |G(A_2) - G(\xi)| \le |G(A_2)| + |G(\xi)| \le 2M$$

---

**第 4 步：合并放缩**。

$$\begin{aligned}
\left|\int_{A_1}^{A_2} f(x)g(x)\,dx\right|
&\le |f(A_1)|\cdot\left|\int_{A_1}^{\xi} g\right| + |f(A_2)|\cdot\left|\int_{\xi}^{A_2} g\right| \\
&\le |f(A_1)|\cdot 2M + |f(A_2)|\cdot 2M \\
&= 2M\big(|f(A_1)| + |f(A_2)|\big)
\end{aligned}$$

由于 $A_1, A_2 \ge X$，由第 1 步知 $|f(A_1)|, |f(A_2)| < \varepsilon/(4M)$。代入：

$$\left|\int_{A_1}^{A_2} f(x)g(x)\,dx\right| < 2M\left(\frac{\varepsilon}{4M} + \frac{\varepsilon}{4M}\right) = \varepsilon$$

由 Cauchy 收敛准则（定理 8.6），$\displaystyle\int_a^{+\infty} f(x)g(x)\,dx$ 收敛。$\square$

---

#### 4.2.2 瑕积分情形的 Dirichlet 判别法

**定理 8.9'（Dirichlet 判别法——瑕积分情形）**：设 $x = a$ 是瑕点，$f(x), g(x)$ 在 $(a, b]$ 上有定义，且对任意 $\eta > 0$，$f$ 和 $g$ 在 $[a+\eta, b]$ 上 Riemann 可积。若：

1. **$f(x)$ 在 $(a, b]$ 上单调**，且 $\displaystyle\lim_{x\to a^+} f(x) = 0$；
2. **$\displaystyle\int_{a+\eta}^b g(x)\,dx$ 一致有界**：存在 $M > 0$，使得对一切 $\eta \in (0, b-a)$：
   $$\left|\int_{a+\eta}^b g(x)\,dx\right| \le M$$

则瑕积分 $\displaystyle\int_a^b f(x)g(x)\,dx$ **收敛**。

**证明思路**：与定理 8.9 的证明完全相同，只需将 Cauchy 准则替换为定理 8.7（瑕积分版本的 Cauchy 准则），并在区间 $[a+\eta_1, a+\eta_2]$ 上应用积分第二中值定理即可。$\square$

---

### 4.3 Abel 判别法

#### 4.3.1 无穷限情形的定理陈述与证明

**定理 8.10（Abel 判别法——无穷限反常积分）**：设 $f(x), g(x)$ 在 $[a, +\infty)$ 上有定义，且对任意 $b > a$，$f$ 和 $g$ 在 $[a, b]$ 上 Riemann 可积。若：

1. **$f(x)$ 单调有界**：$f$ 在 $[a, +\infty)$ 上单调，且存在 $K > 0$ 使得 $|f(x)| \le K$；
2. **$\displaystyle\int_a^{+\infty} g(x)\,dx$ 收敛**。

则反常积分 $\displaystyle\int_a^{+\infty} f(x)g(x)\,dx$ **收敛**。

**证明**（归化为 Dirichlet 判别法）：

**第 1 步：提取 $f$ 的极限**。

由于 $f$ 单调有界，由单调有界定理，极限 $L = \displaystyle\lim_{x\to +\infty} f(x)$ 存在且有限（$|L| \le K$）。

定义 $\tilde{f}(x) = f(x) - L$，则 $\tilde{f}$ 与 $f$ 具有相同的单调性，且：

$$\lim_{x\to +\infty} \tilde{f}(x) = 0$$

---

**第 2 步：验证 $\tilde{f}$ 满足 Dirichlet 条件**。

$\tilde{f}$ 在 $[a, +\infty)$ 上单调且趋于零。接下来验证 $\int_a^A g(x)\,dx$ 的一致有界性：

由于 $\int_a^{+\infty} g(x)\,dx$ 收敛，记 $I = \displaystyle\lim_{A\to +\infty} \int_a^A g(x)\,dx$。则 $\int_a^A g(x)\,dx \to I$ 当 $A\to +\infty$。收敛数列必有界——因此存在 $M > 0$，使得对所有 $A \ge a$：

$$\left|\int_a^A g(x)\,dx\right| \le M$$

（例如可取 $M = \sup_{A\ge a} |\int_a^A g| + 1$，该上确界有限因为收敛数列有界。）

---

**第 3 步：应用 Dirichlet 判别法**。

$\tilde{f}$ 单调趋于零，$\int_a^A g$ 一致有界，由定理 8.9（Dirichlet 判别法）：

$$\int_a^{+\infty} \tilde{f}(x)g(x)\,dx \quad\text{收敛}$$

---

**第 4 步：合并**。

$$\int_a^{+\infty} f(x)g(x)\,dx = \int_a^{+\infty} (\tilde{f}(x) + L)g(x)\,dx = \int_a^{+\infty} \tilde{f}(x)g(x)\,dx + L\int_a^{+\infty} g(x)\,dx$$

右端第一项收敛（由第 3 步），第二项 $L\int_a^{+\infty} g(x)\,dx$ 也收敛（因为 $\int g$ 收敛且 $L$ 是常数）。两项之和收敛，故 $\int_a^{+\infty} f(x)g(x)\,dx$ 收敛。$\square$

---

#### 4.3.2 瑕积分情形的 Abel 判别法

**定理 8.10'（Abel 判别法——瑕积分情形）**：设 $x = a$ 是瑕点，$f(x), g(x)$ 在 $(a, b]$ 上有定义，且对任意 $\eta > 0$，$f$ 和 $g$ 在 $[a+\eta, b]$ 上 Riemann 可积。若：

1. **$f(x)$ 在 $(a, b]$ 上单调有界**；
2. **瑕积分 $\displaystyle\int_a^b g(x)\,dx$ 收敛**。

则瑕积分 $\displaystyle\int_a^b f(x)g(x)\,dx$ **收敛**。

**证明思路**：与定理 8.10 的归化思路相同，将 $f$ 减去其在 $a^+$ 处的极限（由单调有界性保证极限存在），将 $\tilde{f}$ 部分用 Dirichlet 判别法处理。$\square$

---

### 4.4 条件收敛的证明技术：sin²x 拆分法

绝对值积分 $\int |f|$ 发散的证明通常比收敛性证明更困难。对于形如 $\int_1^{+\infty} \frac{\sin x}{x^p}\,dx$（$0 < p \le 1$）的条件收敛积分，有一个经典技巧。

#### 4.4.1 核心思想

要证明 $\int_1^{+\infty} |f(x)|\,dx$ 发散，如果直接放缩 $|f|$ 找不到合适的发散比较函数，可以通过三角恒等式将 $|f|$ **放大为**某发散积分的下界。

关键技术：利用三角恒等式 $\sin^2 x = (1 - \cos 2x)/2$，以及 $|\sin x| \ge \sin^2 x$（因为 $0 \le \sin^2 x \le |\sin x| \le 1$）。

#### 4.4.2 标准流程

以 $I(p) = \int_1^{+\infty} \frac{\sin x}{x^p}\,dx$（$0 < p \le 1$）为例，证明 $\int |\sin x/x^p|$ 发散：

**第 1 步：建立下界**。

$$|\sin x| \ge \sin^2 x = \frac{1 - \cos 2x}{2}, \quad \forall x \in \mathbb{R}$$

因此：

$$\int_1^{+\infty} \frac{|\sin x|}{x^p}\,dx \ge \int_1^{+\infty} \frac{1 - \cos 2x}{2x^p}\,dx$$

**第 2 步：拆分为两项**。

$$\int_1^{+\infty} \frac{1 - \cos 2x}{2x^p}\,dx = \underbrace{\int_1^{+\infty} \frac{1}{2x^p}\,dx}_{:=A} - \underbrace{\int_1^{+\infty} \frac{\cos 2x}{2x^p}\,dx}_{:=B}$$

**第 3 步：分析 A 的敛散性**。

$A = \frac12\int_1^{+\infty} \frac{1}{x^p}\,dx$。由定理 8.1（无穷限 p-积分），当 $p \le 1$ 时该积分 **发散**（到 $+\infty$）。

**第 4 步：分析 B 的敛散性**。

$B = \frac12\int_1^{+\infty} \frac{\cos 2x}{x^p}\,dx$。这里 $\cos 2x$ 是振荡函数。用 Dirichlet 判别法：

- $f(x) = \frac{1}{x^p}$：当 $p > 0$ 时在 $[1, +\infty)$ 上单调递减，且 $\lim_{x\to +\infty} 1/x^p = 0$
- $\int_1^A \cos 2x\,dx = \left[\frac{\sin 2x}{2}\right]_1^A = \frac{\sin 2A - \sin 2}{2}$，故 $|\int_1^A \cos 2x\,dx| \le 1$ 一致有界

由定理 8.9（Dirichlet 判别法），$B$ **收敛**。

**第 5 步：由"发散 − 收敛 = 发散"得出结论**。

由 $A$ 发散、$B$ 收敛，得 $A - B$ 发散（因为若 $A - B$ 收敛，则 $A = (A - B) + B$ 也收敛，矛盾）。因此：

$$\int_1^{+\infty} \frac{|\sin x|}{x^p}\,dx \ge A - B = +\infty$$

即 $\int |\sin x/x^p|$ 发散。结合由 Dirichlet 判别法得到的收敛性，$\int_1^{+\infty} \sin x/x^p\,dx$ **条件收敛**。

---

### 4.5 两类判别法的对比总结

| 判别法 | $f$ 的条件 | $\int g$ 的条件 | 证明策略 |
|--------|-----------|----------------|---------|
| Dirichlet | 单调 $\to 0$ | $\int_a^A g$ 一致有界 | 第二中值定理 + Cauchy 准则 |
| Abel | 单调有界 | $\int_a^{+\infty} g$ 收敛 | $f = \tilde{f} + L$，归化为 Dirichlet |

**实际应用中的选择策略**：

- 如果振荡因子 $g(x)$ 本身的原函数有**一致有界但不一定收敛**（如 $\sin x$、$\cos x$、$\sin(2x)$ 等），使用 **Dirichlet 判别法**；$f$ 必须**趋于零**。
- 如果 $\int g$ 已知收敛（如已由 Dirichlet 判别法证明过的 $\int \sin x/x\,dx$），但需要在被积函数上再乘以一个**单调有界**的因子 $f$ 时，使用 **Abel 判别法**；$f$ 不必趋于零，只需**有界**。
- 如果 $f$ 既趋于零又有界（如 $f(x) = 1/x^p$，$p > 0$），两种判别法理论上都适用，但通常看 $g$ 的条件哪个更容易验证。

---

## 5. 例题

### 例 1：Dirichlet 判别法的完整应用

考虑反常积分 $\displaystyle\int_1^{+\infty} \frac{\sin x}{x}\,dx$。

**(1)** 叙述 Dirichlet 判别法。
**(2)** 验证 $\int_1^A \sin x\,dx$ 的一致有界性。
**(3)** 利用 Dirichlet 判别法证明该积分收敛。
**(4)** 判断该积分的收敛类型。

---

**解 (1)**：

由定理 8.9（Dirichlet 判别法——无穷限情形）：

设 $f(x), g(x)$ 在 $[a, +\infty)$ 上定义，若：
- $f(x)$ 在 $[a, +\infty)$ 上单调，且 $\lim_{x\to +\infty} f(x) = 0$；
- 存在 $M > 0$ 使 $\left|\int_a^A g(x)\,dx\right| \le M$ 对所有 $A \ge a$ 成立，

则 $\displaystyle\int_a^{+\infty} f(x)g(x)\,dx$ 收敛。

---

**解 (2)**：

计算定积分：

$$\int_1^A \sin x\,dx = \big[-\cos x\big]_1^A = (-\cos A) - (-\cos 1) = \cos 1 - \cos A$$

用绝对值不等式放缩：

$$\left|\int_1^A \sin x\,dx\right| = |\cos 1 - \cos A| \le |\cos 1| + |\cos A| \le 1 + 1 = 2$$

因此 $\left|\int_1^A \sin x\,dx\right| \le 2$ 对一切 $A \ge 1$ 成立，一致有界（$M = 2$）。

---

**解 (3)**：

将原积分写为 $\displaystyle\int_1^{+\infty} f(x)g(x)\,dx$ 的形式，取：

$$f(x) = \frac{1}{x}, \quad g(x) = \sin x$$

对 $f$ 的验证：
- $f(x) = 1/x$ 在 $[1, +\infty)$ 上**严格单调递减**（$f'(x) = -1/x^2 < 0$）。
- $\displaystyle\lim_{x\to +\infty} \frac{1}{x} = 0$。✓

对 $g$ 的验证：
- 由 (2)，$\left|\int_1^A \sin x\,dx\right| \le 2$ 对所有 $A \ge 1$ 成立。✓

由定理 8.9（Dirichlet 判别法），$\displaystyle\int_1^{+\infty} \frac{\sin x}{x}\,dx$ **收敛**。

（注：该积分称为 **Dirichlet 积分**，其精确值为 $\pi/2$，但此处只需证明收敛性。）

---

**解 (4)**：

由 (3) 知积分收敛。接下来判断是否绝对收敛。

考虑绝对值的积分 $\displaystyle\int_1^{+\infty} \frac{|\sin x|}{x}\,dx$。

利用 $|\sin x| \ge \sin^2 x = (1 - \cos 2x)/2$（对所有 $x$ 成立）：

$$\int_1^{+\infty} \frac{|\sin x|}{x}\,dx \ge \int_1^{+\infty} \frac{1 - \cos 2x}{2x}\,dx = \underbrace{\int_1^{+\infty} \frac{1}{2x}\,dx}_{A} - \underbrace{\int_1^{+\infty} \frac{\cos 2x}{2x}\,dx}_{B}$$

- $A = \frac12\int_1^{+\infty} \frac{1}{x}\,dx$：$p = 1$ 的 p-积分，由定理 8.1，**发散**。
- $B$：由 Dirichlet 判别法，$1/(2x)$ 单调减趋于 0，$\int_1^A \cos 2x\,dx = \frac{\sin 2A - \sin 2}{2}$ 有界 $\le 1$，故 $B$ **收敛**。

由"发散 − 收敛 = 发散"，$\int_1^{+\infty} \frac{|\sin x|}{x}\,dx$ 发散，即原积分不绝对收敛。

由定义 8.4（条件收敛）：积分收敛但 $\int |f|$ 发散 ⇒ $\displaystyle\int_1^{+\infty} \frac{\sin x}{x}\,dx$ **条件收敛**。

---

### 例 2：含参积分 $I(p) = \displaystyle\int_1^{+\infty} \frac{\sin x}{x^p}\,dx$ 的完整分类讨论

对不同的实数 $p$，判断 $I(p)$ 的敛散性及收敛类型。

---

**解：**

#### 情形 1：$p > 1$ —— 绝对收敛

**第 1 步：取绝对值**。

$$\left|\frac{\sin x}{x^p}\right| = \frac{|\sin x|}{x^p} \le \frac{1}{x^p}, \quad \forall x \ge 1$$

**第 2 步：比较判别法**。

$\int_1^{+\infty} \frac{1}{x^p}\,dx$ 是无穷限 p-积分。由定理 8.1，当 $p > 1$ 时收敛。

由定理 8.3（比较判别法——收敛方向），$0 \le |\sin x|/x^p \le 1/x^p$ 且 $\int 1/x^p$ 收敛，故 $\int |\sin x/x^p|$ 收敛。

**第 3 步：得出结论**。

由定义 8.3，$I(p)$ 绝对收敛（由定理 8.8，因此也收敛）。

---

#### 情形 2：$0 < p \le 1$ —— 条件收敛

**步骤 A：证明 $I(p)$ 收敛**。

将 $I(p)$ 写为 $\int_1^{+\infty} f(x)g(x)\,dx$ 的形式，取：

$$f(x) = \frac{1}{x^p}, \quad g(x) = \sin x$$

验证 Dirichlet 条件：

- $f(x) = 1/x^p$ 在 $[1, +\infty)$ 上单调递减（$p > 0$ 时 $f'(x) = -p/x^{p+1} < 0$）。
- $\lim_{x\to +\infty} 1/x^p = 0$。
- 由例 1 (2)，$\left|\int_1^A \sin x\,dx\right| \le 2$ 对所有 $A \ge 1$ 成立。

由定理 8.9（Dirichlet 判别法），$I(p)$ **收敛**。

---

**步骤 B：证明 $I(p)$ 非绝对收敛**。

考虑 $\displaystyle\int_1^{+\infty} \frac{|\sin x|}{x^p}\,dx$，证明其发散。

**第 1 步：建立下界**。

$$|\sin x| \ge \sin^2 x = \frac{1 - \cos 2x}{2}$$

因此：

$$\int_1^{+\infty} \frac{|\sin x|}{x^p}\,dx \ge \int_1^{+\infty} \frac{1 - \cos 2x}{2x^p}\,dx = \frac12\int_1^{+\infty} \frac{1}{x^p}\,dx - \frac12\int_1^{+\infty} \frac{\cos 2x}{x^p}\,dx$$

**第 2 步：分析两项**。

第一项 $A = \frac12\int_1^{+\infty} \frac{1}{x^p}\,dx$：这是 p-积分乘以 $1/2$。当 $p \le 1$ 时，由定理 8.1，$A$ **发散**（到 $+\infty$）。

第二项 $B = \frac12\int_1^{+\infty} \frac{\cos 2x}{x^p}\,dx$：应用 Dirichlet 判别法。
- $f_B(x) = 1/(2x^p)$ 在 $[1, +\infty)$ 上单调递减且 $\to 0$（$p > 0$）。
- $\int_1^A \cos 2x\,dx = \frac{\sin 2A - \sin 2}{2}$，故 $\left|\int_1^A \cos 2x\,dx\right| \le 1$，一致有界。

由定理 8.9，$B$ **收敛**。

**第 3 步：综合**。

$\int |\sin x/x^p| \ge A - B$。$A$ 发散到 $+\infty$，$B$ 收敛到有限值，故 $A - B$ 发散到 $+\infty$。因此 $\int |\sin x/x^p|$ 发散，$I(p)$ 非绝对收敛。

---

**步骤 C：判断收敛类型**。

由步骤 A（收敛）和步骤 B（非绝对收敛），由定义 8.4，$I(p)$ 在 $0 < p \le 1$ 时 **条件收敛**。

---

#### 情形 3：$p = 0$ —— 振荡发散

$$I(0) = \int_1^{+\infty} \sin x\,dx = \lim_{b\to +\infty} \big[-\cos x\big]_1^b = \lim_{b\to +\infty} (\cos 1 - \cos b)$$

$\cos b$ 在 $b\to +\infty$ 时**振荡**（不趋于任何极限），故该极限不存在。$I(0)$ **发散**。

---

#### 情形 4：$p < 0$ —— 振幅增长导致发散

记 $q = -p > 0$，则：

$$I(p) = \int_1^{+\infty} x^q \sin x\,dx$$

当 $x\to +\infty$ 时，振幅 $x^q \to +\infty$，被积函数的振荡幅度无限增大。用 Cauchy 准则验证：

取 $\varepsilon = 1$。对任意大的 $X$，总能找到 $A_1 = 2n\pi$、$A_2 = (2n+1)\pi$（其中 $n$ 充分大使 $A_1 \ge X$），此时 $\sin x$ 在 $[A_1, A_2]$ 上非负，且：

$$\begin{aligned}
\int_{A_1}^{A_2} x^q \sin x\,dx
&\ge \int_{2n\pi}^{(2n+1)\pi} (2n\pi)^q \sin x\,dx \quad (\text{在 } [2n\pi, (2n+1)\pi] \text{ 上 } \sin x \ge 0,\ x \ge 2n\pi)\\
&= (2n\pi)^q \int_{2n\pi}^{(2n+1)\pi} \sin x\,dx \\
&= (2n\pi)^q \cdot \big[-\cos x\big]_{2n\pi}^{(2n+1)\pi} \\
&= (2n\pi)^q \cdot (-\cos(2n+1)\pi + \cos 2n\pi) \\
&= (2n\pi)^q \cdot (-(-1) + 1) = 2(2n\pi)^q \to +\infty \ (\text{当 } n\to\infty)
\end{aligned}$$

因此对任意 $X$，都存在 $A_2 > A_1 \ge X$ 使 $\left|\int_{A_1}^{A_2} x^q \sin x\,dx\right| > 1$，不满足 Cauchy 准则（定理 8.6）。$I(p)$ 在 $p < 0$ 时 **发散**。

---

**分类汇总表**：

| $p$ 的范围 | 敛散性 | 收敛类型 | 依据 |
|-----------|--------|---------|------|
| $p > 1$ | 收敛 | **绝对收敛** | $|\sin x/x^p| \le 1/x^p$，$p>1$ 时 p-积分收敛 |
| $0 < p \le 1$ | 收敛 | **条件收敛** | Dirichlet 判别法 + sin²x 拆分 |
| $p = 0$ | 发散 | — | $\int \sin x\,dx$ 振荡 |
| $p < 0$ | 发散 | — | 振幅 $x^{-p}$ 无限增长 |

---

### 例 3：Abel 判别法的应用

已知 $\displaystyle\int_1^{+\infty} \frac{\sin x}{x}\,dx$ 收敛（例 1 结论）。利用 Abel 判别法证明：

$$\int_1^{+\infty} \frac{\sin x}{x+1}\,dx$$

收敛。

**解**：

将原积分写为 $\int_1^{+\infty} f(x)g(x)\,dx$ 的形式，取：

$$f(x) = \frac{x}{x+1}, \quad g(x) = \frac{\sin x}{x}$$

注意 $\frac{\sin x}{x} \cdot \frac{x}{x+1} = \frac{\sin x}{x+1}$，确实匹配。

---

**第 1 步：验证 $f(x) = x/(x+1)$ 单调有界**。

改写 $f(x)$：

$$f(x) = \frac{x}{x+1} = 1 - \frac{1}{x+1}$$

求导验证单调性：

$$f'(x) = \frac{(x+1)\cdot 1 - x\cdot 1}{(x+1)^2} = \frac{1}{(x+1)^2} > 0, \quad \forall x \ge 1$$

因此 $f(x)$ 在 $[1, +\infty)$ 上严格单调递增。

有界性：对 $x \ge 1$，

$$0 < \frac{1}{x+1} \le \frac12 \quad\Rightarrow\quad \frac12 \le 1 - \frac{1}{x+1} < 1$$

因此 $|f(x)| < 1$，$f$ 有界（可取 $K = 1$）。✓

---

**第 2 步：验证 $\displaystyle\int_1^{+\infty} g(x)\,dx = \int_1^{+\infty} \frac{\sin x}{x}\,dx$ 收敛**。

由例 1 (3) 的结论，该积分收敛。✓

---

**第 3 步：应用 Abel 判别法**。

$f(x)$ 在 $[1, +\infty)$ 上单调有界，$\int_1^{+\infty} g(x)\,dx$ 收敛。由定理 8.10（Abel 判别法）：

$$\int_1^{+\infty} f(x)g(x)\,dx = \int_1^{+\infty} \frac{\sin x}{x+1}\,dx$$

**收敛**。

---

**评注**：本例展示了 Abel 判别法的典型使用场景——**在已知某个积分 $\int g$ 收敛的基础上，乘以一个单调有界因子 $f$ 后，收敛性保持不变**。这与 Dirichlet 判别法形成对比：Dirichlet 需要 $f$ 趋于零和 $\int g$ 的有界性（而非收敛性），条件更"宽松"地适用于 $\int g$ 本身可能发散（如 $\int \sin x\,dx$ 发散但有界）的情形。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| Dirichlet 判别法和 Abel 判别法条件相同，只是名字不同 | 两者条件**恰好互换**：Dirichlet 要求 $f$ **单调趋于零** + $\int g$ 的变上限积分**一致有界**；Abel 要求 $f$ **单调有界** + $\int g$ **收敛**。前者对 $\int g$ 条件更弱（有界而非收敛），但对 $f$ 条件更强（趋于 0 而非有界） |
| 想证 $\int fg$ 收敛，可以直接对 $|fg|$ 用比较判别法 | 比较判别法要求被积函数**非负**。$fg$ 振荡变号时不满足前提。正确做法是对 $f$ 和 $g$ 分别提取单调性和有界性，用 Dirichlet/Abel 判别法处理 |
| $|\sin x| \ge \sin^2 x$ 的放缩方向搞反 | $|\sin x| \ge \sin^2 x$ 是正确的，因为 $0 \le \sin^2 x \le |\sin x| \le 1$。这个"往下放缩"的目的是**找 $\int |f|$ 的下界**，证其发散（而非收敛） |
| 证明条件收敛时，只需证收敛即可，$\int |f|$ 发散的部分不重要 | 条件收敛要求**两者都验证**：(1) $\int f$ 收敛；(2) $\int |f|$ 发散。缺一不可。仅验证收敛只能得到"收敛"，无法区分绝对收敛与条件收敛 |
| Abel 判别法的证明必须用积分第二中值定理重新证 | Abel 判别法可以通过归化为 Dirichlet 判别法来证明，无需从零开始。方法是构造 $\tilde{f} = f - L$（$L$ 为 $f$ 的极限），$\tilde{f}$ 满足 Dirichlet 条件 |
| $\int \sin x/x^p\,dx$ 在 $p<0$ 时发散，因为被积函数趋于无穷 | 准确地说是**振幅无限增长**导致发散。$\sin x$ 本身在 $[-1,1]$ 间振荡，但乘以 $x^{|p|}$ 后振幅趋于无穷，正负部分的绝对面积都无限增大，无法抵消 |
| Dirichlet 判别法中 $f$ 单调且趋于 0 就够，不需要检查单调性方向 | 单调性（递增或递减均可）是为了应用积分第二中值定理。第二中值定理要求 $f$ 单调（不一定非负）。$f$ 递减且趋于 0 是最常见的应用场景，但递增也成立 |

### 检查点

- [ ] 能否写出积分第二中值定理（Bonnet 形式）的完整陈述？
- [ ] 能否写出 Dirichlet 判别法（定理 8.9）的完整陈述并给出其基于 Cauchy 准则的证明？
- [ ] 能否写出 Abel 判别法（定理 8.10）的完整陈述并给出其通过归化为 Dirichlet 判别法的证明？
- [ ] 能否独立完成 $\int_1^{+\infty} \sin x/x\,dx$ 的：(a) 用 Dirichlet 判别法证收敛；(b) 用 sin²x 技巧证非绝对收敛？
- [ ] 能否说出 Dirichlet 和 Abel 判别法在条件上的核心差异？
- [ ] 对于 $I(p) = \int_1^{+\infty} \sin x/x^p\,dx$，能否说出每种 $p$ 对应的敛散性和收敛类型？
- [ ] 能否解释为什么 $\int_1^{+\infty} \sin x/(x+1)\,dx$ 可以用 Abel 判别法而 $\int_1^{+\infty} \sin x/x\,dx$ 不能用 Abel 判别法（只能用 Dirichlet）？
- [ ] 能否正确写出一致有界的定义（即 $\exists M: |\int_a^A g| \le M$ 对所有 $A$ 成立）？
- [ ] 能否完成 sin²x 拆分的关键推导：$|\sin x| \ge \sin^2 x = (1-\cos 2x)/2$ → 拆为两项 → 分别判定敛散性 → 综合结论？
- [ ] 能否解释 $\int \sin x/x\,dx$ 在 $p=0$ 时为何发散（振荡而非趋于无穷）与 $p<0$ 时为何发散（振幅增长）的差异？

---

## 练习题

### 基础巩固

**1.**（定理理解）判断下列说法是否正确，并说明理由：

(a) 若 $f$ 在 $[1, +\infty)$ 上单调且 $\lim_{x\to +\infty} f(x) = 0$，则 $\int_1^{+\infty} f(x)\sin x\,dx$ 必收敛。

(b) 若 $\int_1^{+\infty} g(x)\,dx$ 收敛且 $f(x)$ 有界，则 $\int_1^{+\infty} f(x)g(x)\,dx$ 必收敛。

(c) 若 $f$ 单调有界，$\left|\int_1^A g(x)\,dx\right| \le M$ 对所有 $A\ge 1$ 成立，则 $\int_1^{+\infty} f(x)g(x)\,dx$ 收敛。

<details><summary>参考答案</summary>

**(a) 不一定**。Dirichlet 判别法还需要验证 $\int_1^A \sin x\,dx$ 的一致有界性。$\sin x$ 满足这个条件（有界 $M=2$），因此**实际上** $\int_1^{+\infty} f(x)\sin x\,dx$ 确实收敛——但"必收敛"的结论依赖于 $f$ 单调趋于零 **且** $\int \sin$ 的有界性，不能仅凭 $f$ 的条件下结论。正确说法："若 $f$ 单调趋于零，则由 Dirichlet 判别法，$\int f \sin$ 收敛"——隐含地引用了 $\sin$ 的有界性。

**(b) 不一定**。Abel 判别法要求 $f$ **单调**有界，仅仅是"有界"（不单调）不够。反例：令 $f(x)$ 在 $[n, n+1]$ 上取 $1$ 和 $-1$ 交替（不单调但有界 $|f| \le 1$），$g(x)$ 取适当函数可使 $\int fg$ 发散。Abel 判别法中的单调性是使用积分第二中值定理的必要条件。

**(c) 不一定**。此时 $f$ 的条件更接近 Dirichlet 而非 Abel，但 Dirichlet 要求 $f\to 0$ 而非"有界"。若 $f$ 单调有界但不趋于零（如 $f(x) = 2 + 1/x$），$\int fg$ 不一定收敛。只有 $f$ 单调趋于零时才满足 Dirichlet 条件。

</details>

---

**2.**（Dirichlet 判别法验证）用 Dirichlet 判别法证明 $\displaystyle\int_1^{+\infty} \frac{\cos x}{x^{3/2}}\,dx$ 收敛，并判断其是否绝对收敛。

<details><summary>参考答案</summary>

**证明收敛（Dirichlet 判别法）**：

取 $f(x) = 1/x^{3/2}$，$g(x) = \cos x$。

- $f(x) = 1/x^{3/2}$ 在 $[1, +\infty)$ 上单调递减，且 $\lim_{x\to +\infty} 1/x^{3/2} = 0$。✓
- $\int_1^A \cos x\,dx = \sin A - \sin 1$，故 $\left|\int_1^A \cos x\,dx\right| \le |\sin A| + |\sin 1| \le 2$，一致有界。✓

由定理 8.9（Dirichlet 判别法），$\int_1^{+\infty} \cos x/x^{3/2}\,dx$ 收敛。

---

**判断绝对收敛**：

$$\left|\frac{\cos x}{x^{3/2}}\right| \le \frac{1}{x^{3/2}}$$

$\int_1^{+\infty} 1/x^{3/2}\,dx$ 是 p-积分，$p = 3/2 > 1$，由定理 8.1 收敛。

由定理 8.3（比较判别法收敛方向），$\int |\cos x/x^{3/2}|\,dx$ 收敛。

由定义 8.3，原积分**绝对收敛**。

**结论**：$\int_1^{+\infty} \cos x/x^{3/2}\,dx$ 绝对收敛（因此当然也收敛）。

</details>

---

### 迁移应用

**3.**（Abel 判别法应用）已知 $\displaystyle\int_1^{+\infty} \frac{\cos x}{x^{3/2}}\,dx$ 收敛（练习 2 结论）。证明 $\displaystyle\int_1^{+\infty} \frac{x+2}{x+1}\cdot\frac{\cos x}{x^{3/2}}\,dx$ 收敛。

<details><summary>参考答案</summary>

将原积分写为 $\int_1^{+\infty} f(x)g(x)\,dx$ 的形式，取：

$$f(x) = \frac{x+2}{x+1}, \quad g(x) = \frac{\cos x}{x^{3/2}}$$

**验证 $f$ 单调有界**：

$$f(x) = \frac{x+2}{x+1} = 1 + \frac{1}{x+1}$$

$$f'(x) = -\frac{1}{(x+1)^2} < 0, \quad \forall x \ge 1$$

故 $f$ 在 $[1, +\infty)$ 上**严格单调递减**。

有界性：对 $x \ge 1$，$0 < \frac{1}{x+1} \le \frac12$，故 $1 < f(x) \le \frac32$，$|f(x)| \le \frac32$。✓

**验证 $\int g$ 收敛**：

由练习 2 结论，$\int_1^{+\infty} \frac{\cos x}{x^{3/2}}\,dx$ 收敛。✓

**应用 Abel 判别法**：

由定理 8.10（Abel 判别法），$\int_1^{+\infty} f(x)g(x)\,dx$ 收敛。

即 $\displaystyle\int_1^{+\infty} \frac{x+2}{x+1}\cdot\frac{\cos x}{x^{3/2}}\,dx = \int_1^{+\infty} \frac{(x+2)\cos x}{(x+1)x^{3/2}}\,dx$ 收敛。

</details>

---

**4.**（综合——条件收敛证明）证明 $\displaystyle\int_1^{+\infty} \frac{\cos x}{x^{p}}\,dx$ 在 $0 < p \le 1$ 时条件收敛。

**提示**：模仿例 2 的步骤：先用 Dirichlet 判别法证明收敛，再用 $|\cos x| \ge \cos^2 x = (1 + \cos 2x)/2$ 证明非绝对收敛。

<details><summary>参考答案</summary>

**第 1 步：证明收敛（Dirichlet 判别法）**。

取 $f(x) = 1/x^p$，$g(x) = \cos x$。

- $f(x) = 1/x^p$ 在 $[1, +\infty)$ 上单调递减，且 $\lim_{x\to +\infty} 1/x^p = 0$（$p > 0$）。✓
- $\int_1^A \cos x\,dx = \sin A - \sin 1$，故 $\left|\int_1^A \cos x\,dx\right| \le 2$。✓

由定理 8.9（Dirichlet 判别法），$\int_1^{+\infty} \cos x/x^p\,dx$ **收敛**。

---

**第 2 步：证明 $\int |\cos x/x^p|$ 发散**。

利用 $|\cos x| \ge \cos^2 x = \dfrac{1 + \cos 2x}{2}$：

$$\int_1^{+\infty} \frac{|\cos x|}{x^p}\,dx \ge \int_1^{+\infty} \frac{1 + \cos 2x}{2x^p}\,dx = \underbrace{\int_1^{+\infty} \frac{1}{2x^p}\,dx}_{A} + \underbrace{\int_1^{+\infty} \frac{\cos 2x}{2x^p}\,dx}_{B}$$

---

**分析 $A$**：

$A = \frac12\int_1^{+\infty} 1/x^p\,dx$，$p \le 1$ 时发散。✓

**分析 $B$**：

应用 Dirichlet 判别法：
- $f_B(x) = 1/(2x^p)$ 单调减趋于 0（$p > 0$）。
- $\int_1^A \cos 2x\,dx = \frac{\sin 2A - \sin 2}{2}$，$\left|\int_1^A \cos 2x\,dx\right| \le 1$。

由 Dirichlet 判别法，$B$ **收敛**。

---

**第 3 步：综合**。

$\int |\cos x/x^p| \ge A + B$，其中 $A$ 发散（到 $+\infty$），$B$ 收敛。发散 + 收敛 = 发散，故 $\int |\cos x/x^p|$ 发散。

**结论**：$\int_1^{+\infty} \cos x/x^p\,dx$ 在 $0 < p \le 1$ 时**条件收敛**。

</details>
