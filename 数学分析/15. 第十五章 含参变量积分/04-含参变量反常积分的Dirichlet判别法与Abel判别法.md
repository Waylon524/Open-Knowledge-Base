# 04. 含参变量反常积分的 Dirichlet 判别法与 Abel 判别法

> 所属章节：第十五章 含参变量积分  |  文件序号：04  |  难度：核心
> 常见混淆点：1) Dirichlet 与 Abel 判别法的条件恰好"互换"——Dirichlet 要求 f 一致趋于零而 ∫g 的变上限积分一致有界，Abel 要求 f 一致有界而 ∫g 一致收敛；2) "一致趋于零"需要关于参数 y 一致（这就是"一致"一词的关键），而非仅逐点趋于零；3) 内闭一致收敛 ≠ 整体一致收敛——理解反例 sin(αx)/x 在 (0,+∞) 上的行为；4) M-判别法给出的是绝对一致收敛，Dirichlet/Abel 允许条件一致收敛，三者形成完整的判别体系

---

## 1. 学习目标与先修前置

### 学习目标
- 掌握含参变量反常积分的 **Dirichlet 判别法**（定理 15.6）的完整陈述与证明
- 掌握含参变量反常积分的 **Abel 判别法**（定理 15.7）的完整陈述与证明
- 理解两个判别法的证明均以**积分第二中值定理**为核心工具，证明逻辑与单变量版本平行
- 掌握 **内闭一致收敛**（定义 15.7）的准确定义，理解其与整体一致收敛的区别与联系
- 能熟练应用 Dirichlet 判别法处理形如 $\int f(x,y)g(x,y)\,dx$ 的含振荡因子积分（$f$ 单调趋于零，$\int g$ 一致有界）
- 能熟练应用 Abel 判别法处理在已知一致收敛积分上乘以单调有界因子的情形
- 能通过对比表系统区分 M-判别法、Dirichlet 判别法、Abel 判别法的适用条件与结论
- 能灵活组合三种判别法处理复杂问题（如 $G(\alpha)=\int_0^{+\infty} e^{-\alpha x}\frac{\sin x}{x}\,dx$）

### 先修知识
- **文件 03（本章）**：含参变量反常积分的一致收敛定义（定义 15.6）、Cauchy 准则（定理 15.4）、Weierstrass M-判别法（定理 15.5）。其中 4.6 节已预告本知识点。
- **文件 04（第八章）**：单变量反常积分的 Dirichlet 判别法（定理 8.9）与 Abel 判别法（定理 8.10）以及积分第二中值定理（引理 8.1）。本文件将其推广到含参变量情形，证明结构完全平行，但条件中增加了"关于参数 y 一致"的要求。
- **文件 02（第十章）**：函数项级数的 Abcl 判别法与 Dirichlet 判别法（ch10-06）——含参变量积分版本与函数项级数版本高度平行，类比迁移有助于理解。
- **文件 01（本章）**：连续性定理（定理 15.1）——内闭一致收敛性是后续将连续性、可微性推广到反常含参变量积分的基础。

---

## 2. 背景与应用场景

### 2.1 M-判别法的局限与动机

文件 03 的 Weierstrass M-判别法（定理 15.5）给出了含参变量反常积分一致收敛的一个简便充分的判别工具。但其要求被积函数**被一个可积优函数绝对控制**，即

$$|f(x,y)| \leq M(x),\quad \int_a^{+\infty} M(x)\,dx < \infty$$

这一条件等价于 $\int_a^{+\infty} |f(x,y)|\,dx$ 一致收敛（即绝对一致收敛）。然而，许多重要的含参变量反常积分是**条件收敛**的——被积函数不绝对可积，但通过振荡因子的正负抵消实现收敛。对于这类积分，M-判别法无法直接应用。

**核心例子**（已在文件 03 的 4.6 节末尾预告）：

$$\Phi(\alpha) = \int_1^{+\infty} \frac{\sin(\alpha x)}{x}\,dx,\quad \alpha > 0$$

对被积函数取绝对值：$\left|\frac{\sin(\alpha x)}{x}\right| \leq \frac{1}{x}$，但 $\int_1^{+\infty} \frac{1}{x}\,dx$ 发散。M-判别法失效。然而，该积分在 $\alpha>0$ 时是条件收敛的，且在 $[\alpha_0,+\infty)$（$\alpha_0>0$）上一致收敛——这正是 Dirichlet 判别法的典型应用场景。

### 2.2 核心思想：乘积结构的处理

与单变量版本（ch8-04）相同，含参变量版本的 Dirichlet 和 Abel 判别法专门处理如下乘积形式的被积函数：

$$F(y) = \int_a^{+\infty} f(x,y)g(x,y)\,dx$$

其中：
- $f(x,y)$ 是**单调因子**（关于 $x$ 单调），控制衰减或增长
- $g(x,y)$ 是**振荡因子**或**已知收敛的因子**

与 M-判别法不同，Dirichlet/Abel 判别法不要求 $|f(x,y)g(x,y)|$ 绝对可积，而是利用 $f$ 的单调性和 $g$ 的积分性质来推断 $\int fg$ 的一致收敛性。

### 2.3 知识图谱

```
文件 01: 含参变量常义积分的定义与连续性（定理 15.1）
文件 02: 可微性与积分次序交换（定理 15.2/15.3）
文件 03: 一致收敛定义（定义 15.6）+ Cauchy 准则（定理 15.4）+ M-判别法（定理 15.5）
       ↓
文件 04: 内闭一致收敛（定义 15.7）
       → Dirichlet 判别法（定理 15.6）：f 一致趋于零 + ∫g 一致有界 ⇒ ∫fg 一致收敛
       → Abel 判别法（定理 15.7）：f 一致有界 + ∫g 一致收敛 ⇒ ∫fg 一致收敛
       → 三种判别法的对比与综合应用
       ↓
文件 05+: 利用一致收敛性推广连续性/可微性到反常含参变量积分
```

---

## 3. 核心概念与符号约定

### 3.1 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $f(x,y)$ | 乘积分解中的"单调因子" | 关于 $x$ 单调（递增或递减） |
| $g(x,y)$ | 乘积分解中的"振荡因子"或"被乘因子" | 与 $f$ 构成被积函数 $fg$ |
| $G_A(y)=\int_a^A g(x,y)\,dx$ | $g$ 的变上限积分 | Dirichlet 条件中要求其一致有界 |
| $\xi$ | 积分第二中值定理中的中间点 | $[A_1,A_2]$ 内的某点 |
| $M$ | $\left|\int_a^A g(x,y)\,dx\right|\leq M$ 中的上界 | Dirichlet 条件的有界常数 |
| $K$ | $|f(x,y)|\leq K$ 中的上界 | Abel 条件的有界常数 |
| $L(y)=\lim_{x\to\infty} f(x,y)$ | 单调有界函数 $f$ 的极限 | 依赖于 $y$，但存在且有限 |
| $\tilde{f}(x,y)=f(x,y)-L(y)$ | 减去极限后的函数 | 满足 Dirichlet 条件（趋于零） |
| $K\subset E$ | $E$ 的紧子集（有界闭子集） | 用于内闭一致收敛的定义 |
| $R_A(y)=\int_A^{+\infty} f(x,y)\,dx$ | 余项函数 | $\sup_{y\in E}|R_A(y)|\to 0$ 即一致收敛 |

### 3.2 两种判别法的条件速览

| 判别法 | 对 $f$ 的条件 | 对 $\int g$ 的条件 | 结论 |
|--------|-------------|-------------------|------|
| Dirichlet（定理 15.6） | 单调，$\lim_{x\to\infty} f=0$ 关于 $y$ **一致** | $\int_a^A g\,dx$ **一致有界** | $\int fg$ 一致收敛 |
| Abel（定理 15.7） | 单调，**一致有界** $|f|\leq K$ | $\int_a^{+\infty} g\,dx$ **一致收敛** | $\int fg$ 一致收敛 |

两者的条件恰好"互换"：Dirichlet 对 $f$ 要求更强（趋于零），对 $\int g$ 要求更弱（有界而非收敛）；Abel 对 $f$ 要求更弱（有界而非趋于零），对 $\int g$ 要求更强（收敛）。

---

## 4. 原理与方法

### 4.1 内闭一致收敛

在许多问题中，含参变量反常积分在参数定义域的**整体**上不一致收敛，但在每个**紧子集**上一致收敛。这种"局部一致收敛"的概念足以保证连续性等分析性质，且比整体一致收敛更容易验证。

**定义 15.7（内闭一致收敛）**：设 $F(y)=\displaystyle\int_a^{+\infty} f(x,y)\,dx$ 对每个 $y\in E$（$E\subseteq\mathbb{R}$）收敛。若对 $E$ 的任意有界闭子集（紧子集）$K\subset E$，$F(y)$ 在 $K$ 上一致收敛，则称 $F(y)$ 在 $E$ 上**内闭一致收敛**（uniformly convergent on every compact subset）。

---

**与整体一致收敛的关系**：

- **整体一致收敛 ⇒ 内闭一致收敛**：若 $F$ 在 $E$ 上整体一致收敛，则对任意 $K\subset E$，$F$ 在 $K$ 上自然也一致收敛（因为 $K$ 上的余项被 $E$ 上的余项控制）。
- **内闭一致收敛 $\not\Rightarrow$ 整体一致收敛**：存在反例。经典例子为：

$$\Phi(\alpha)=\int_1^{+\infty}\frac{\sin(\alpha x)}{x}\,dx,\quad \alpha\in(0,+\infty)$$

该积分在 $(0,+\infty)$ 上内闭一致收敛，但在 $(0,+\infty)$ 上**不**整体一致收敛（详见例 4-2）。

- **内闭一致收敛的充分性**：对于连续性、可积性、可微性的研究，内闭一致收敛通常足够——因为分析性质是局部性质，只需在每个点的邻域内验证一致收敛即可。

---

**几何含义**：

内闭一致收敛的含义是：对任意 $\varepsilon>0$ 和任意紧集 $K\subset E$，存在 $A_0(\varepsilon,K)$，使得当 $A>A_0$ 时，对**所有** $y\in K$ 同时有 $|R_A(y)|<\varepsilon$。这里 $A_0$ 可以依赖于 $K$，但对给定的 $K$，$A_0$ 与 $y\in K$ 无关。

### 4.2 积分第二中值定理（回顾）

Dirichlet 判别法和 Abel 判别法的证明依赖于积分第二中值定理，该定理已在第八章的引理 8.1 中给出证明。这里完整重述以便引用。

**引理 8.1（积分第二中值定理 / Bonnet 形式）**：设 $f$ 在 $[a,b]$ 上单调，$g$ 在 $[a,b]$ 上 Riemann 可积。则存在 $\xi\in[a,b]$，使得

$$\int_a^b f(x)g(x)\,dx = f(a)\int_a^\xi g(x)\,dx + f(b)\int_\xi^b g(x)\,dx$$

**特别地**：
- 若 $f$ 在 $[a,b]$ 上非负单调递减，则 $\exists\xi\in[a,b]$：$\displaystyle\int_a^b f(x)g(x)\,dx = f(a)\int_a^\xi g(x)\,dx$
- 若 $f$ 在 $[a,b]$ 上非负单调递增，则 $\exists\xi\in[a,b]$：$\displaystyle\int_a^b f(x)g(x)\,dx = f(b)\int_\xi^b g(x)\,dx$

**含参变量版本的使用**：当 $f(x,y)$ 关于 $x$ 单调时，对每个固定的 $y$ 均可应用该引理，得到的 $\xi$ 可能依赖于 $y$——但在证明中这不会造成问题，因为中间点 $\xi$ 的存在性已经足够。

---

### 4.3 Dirichlet 判别法（含参变量版本）

**定理 15.6（Dirichlet 判别法——含参变量反常积分）**：设 $f(x,y),g(x,y)$ 定义在 $[a,+\infty)\times E$ 上，且对任意 $A>a$ 和任意 $y\in E$，$f(\cdot,y)$ 和 $g(\cdot,y)$ 在 $[a,A]$ 上 Riemann 可积。若：

1. **$f$ 关于 $x$ 单调且一致趋于零**：对每个固定的 $y\in E$，$f(x,y)$ 关于 $x$ 在 $[a,+\infty)$ 上单调；且
   $$\lim_{x\to +\infty} f(x,y)=0\quad\text{关于 }y\in E\text{ 一致成立}$$
   即 $\forall\varepsilon>0$，$\exists X\geq a$，$\forall x>X$，$\forall y\in E$：$|f(x,y)|<\varepsilon$。

2. **$\displaystyle\int_a^A g(x,y)\,dx$ 一致有界**：存在常数 $M>0$，使得对一切 $A\geq a$ 和一切 $y\in E$ 成立
   $$\left|\int_a^A g(x,y)\,dx\right| \leq M$$

则含参变量反常积分 $\displaystyle F(y)=\int_a^{+\infty} f(x,y)g(x,y)\,dx$ 在 $E$ 上**一致收敛**。

---

**与单变量版本（定理 8.9）的对比**：

| 对比维度 | 单变量版本（定理 8.9） | 含参变量版本（定理 15.6） |
|----------|----------------------|-------------------------|
| 条件 1 | $\lim_{x\to\infty} f(x)=0$（普通极限） | $\lim_{x\to\infty} f(x,y)=0$ **关于 $y$ 一致** |
| 条件 2 | $\sup_A\left|\int_a^A g\right|<\infty$（单个积分） | $\sup_{A,y}\left|\int_a^A g(x,y)dx\right|<\infty$ **对 $y$ 一致** |
| 结论 | $\int fg$ 收敛 | $\int fg$ 在 $E$ 上**一致收敛**（$\sup_y$ 意义下） |
| 证明工具 | 第二中值定理 + Cauchy 准则 | 第二中值定理 + 定理 15.4（含参 Cauchy 准则） |

含参版本在条件中增加了**关于参数 $y$ 的一致性**要求，结论也从"收敛"加强为"一致收敛"——这正是含参变量积分理论的核心关切。

---

**证明**：

我们使用含参变量积分的 Cauchy 准则（定理 15.4）来完成证明。对任意 $\varepsilon>0$，需要找到 $A_0\geq a$，使得对任意 $A_2>A_1\geq A_0$ 和任意 $y\in E$ 有

$$\left|\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx\right| < \varepsilon$$

---

**第 1 步：利用 $f$ 的一致趋于零控制端点值**。

由条件 1（$f$ 一致趋于零），对给定的 $\varepsilon>0$，存在 $X_1\geq a$，使得对所有 $x>X_1$ 和所有 $y\in E$ 有

$$|f(x,y)| < \frac{\varepsilon}{4M}$$

其中 $M$ 是条件 2 中的有界常数。

---

**第 2 步：在 $[A_1,A_2]$ 上应用积分第二中值定理**。

取任意 $A_2>A_1\geq X_1$。对每个固定的 $y\in E$，$f(\cdot,y)$ 在 $[A_1,A_2]$ 上单调。由引理 8.1（Bonnet 形式），存在 $\xi=\xi(y)\in[A_1,A_2]$（可依赖 $y$），使得

$$\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx = f(A_1,y)\int_{A_1}^{\xi} g(x,y)\,dx + f(A_2,y)\int_{\xi}^{A_2} g(x,y)\,dx$$

---

**第 3 步：用 $G_A(y)=\int_a^A g(x,y)\,dx$ 的一致有界性控制中间积分**。

记 $G_A(y)=\displaystyle\int_a^A g(x,y)\,dx$。由条件 2，$|G_A(y)|\leq M$ 对一切 $A\geq a$ 和一切 $y\in E$ 成立。

于是：

$$\begin{aligned}
\left|\int_{A_1}^{\xi} g(x,y)\,dx\right| &= |G_\xi(y)-G_{A_1}(y)| \leq |G_\xi(y)|+|G_{A_1}(y)| \leq 2M \\
\left|\int_{\xi}^{A_2} g(x,y)\,dx\right| &= |G_{A_2}(y)-G_\xi(y)| \leq |G_{A_2}(y)|+|G_\xi(y)| \leq 2M
\end{aligned}$$

---

**第 4 步：合并放缩**。

$$\begin{aligned}
\left|\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx\right|
&\leq |f(A_1,y)|\cdot\left|\int_{A_1}^{\xi} g\right| + |f(A_2,y)|\cdot\left|\int_{\xi}^{A_2} g\right| \\
&\leq |f(A_1,y)|\cdot 2M + |f(A_2,y)|\cdot 2M \\
&= 2M\big(|f(A_1,y)| + |f(A_2,y)|\big)
\end{aligned}$$

由于 $A_1,A_2\geq X_1$，由第 1 步知 $|f(A_1,y)|,|f(A_2,y)|<\varepsilon/(4M)$ 对所有 $y\in E$ 成立。代入：

$$\left|\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx\right| < 2M\left(\frac{\varepsilon}{4M} + \frac{\varepsilon}{4M}\right) = \varepsilon$$

---

**第 5 步：由 Cauchy 准则得一致收敛**。

由 Cauchy 准则（定理 15.4），$\displaystyle F(y)=\int_a^{+\infty} f(x,y)g(x,y)\,dx$ 在 $E$ 上一致收敛。$\square$

---

**定理 15.6 的证明逻辑链**：

$$
\begin{aligned}
&\text{条件 1：} f\to 0\text{ 一致于 }y \xrightarrow{\text{取 }X_1} |f(A_1)|,|f(A_2)|<\varepsilon/(4M) \\
&\text{条件 2：} \sup_{A,y}\left|\int_a^A g\right|\leq M \xrightarrow{\text{三角不等式}} \left|\int_{A_1}^\xi g\right|,\left|\int_\xi^{A_2} g\right|\leq 2M \\
&\xrightarrow{\text{第二中值定理}} \left|\int_{A_1}^{A_2} fg\right| \leq 2M(|f(A_1)|+|f(A_2)|) < \varepsilon \\
&\xrightarrow{\text{定理 15.4}} \int_a^{+\infty} fg \text{ 一致收敛}
\end{aligned}
$$

---

#### 4.3.1 瑕积分情形的 Dirichlet 判别法

对于含参变量瑕积分，只需将"远离积分上限"替换为"靠近瑕点"即可。

**定理 15.6'（Dirichlet 判别法——瑕积分情形）**：设 $x=b$ 是瑕点，$f(x,y),g(x,y)$ 定义在 $[a,b)\times E$ 上，且对任意 $\eta>0$，$f(\cdot,y)$ 和 $g(\cdot,y)$ 在 $[a,b-\eta]$ 上 Riemann 可积。若：

1. **$f$ 关于 $x$ 单调且一致趋于零**：对每个 $y\in E$，$f(x,y)$ 关于 $x$ 单调，且 $\displaystyle\lim_{x\to b^-} f(x,y)=0$ 关于 $y\in E$ 一致成立；
2. **$\displaystyle\int_{a+\eta}^b g(x,y)\,dx$ 一致有界**：存在 $M>0$，使得对一切 $\eta\in(0,b-a)$ 和一切 $y\in E$ 有 $\left|\int_{a+\eta}^b g(x,y)\,dx\right|\leq M$（若瑕点在左端点 $x=a$，则改为 $\int_{a+\eta}^b$）。

则瑕积分 $\displaystyle\int_a^b f(x,y)g(x,y)\,dx$ 在 $E$ 上一致收敛。

**证明思路**：与定理 15.6 的证明完全相同，只需将 Cauchy 准则替换为瑕积分版本（对称于定理 15.4 的瑕积分形式），并在 $[b-\eta_1,b-\eta_2]$ 上应用积分第二中值定理。$\square$

---

### 4.4 Abel 判别法（含参变量版本）

**定理 15.7（Abel 判别法——含参变量反常积分）**：设 $f(x,y),g(x,y)$ 定义在 $[a,+\infty)\times E$ 上，且对任意 $A>a$ 和任意 $y\in E$，$f(\cdot,y)$ 和 $g(\cdot,y)$ 在 $[a,A]$ 上 Riemann 可积。若：

1. **$f$ 关于 $x$ 单调且一致有界**：对每个固定的 $y\in E$，$f(x,y)$ 关于 $x$ 在 $[a,+\infty)$ 上单调；且存在常数 $K>0$，使得 $|f(x,y)|\leq K$ 对一切 $x\geq a$ 和一切 $y\in E$ 成立。
2. **$\displaystyle\int_a^{+\infty} g(x,y)\,dx$ 一致收敛**：含参变量反常积分 $G(y)=\int_a^{+\infty} g(x,y)\,dx$ 在 $E$ 上一致收敛。

则含参变量反常积分 $\displaystyle F(y)=\int_a^{+\infty} f(x,y)g(x,y)\,dx$ 在 $E$ 上**一致收敛**。

---

**证明**（归化为 Dirichlet 判别法）：

**第 1 步：提取 $f$ 的极限函数**。

对每个固定的 $y\in E$，$f(\cdot,y)$ 单调有界。由单调有界定理，极限

$$L(y) = \lim_{x\to +\infty} f(x,y)$$

存在且有限（$|L(y)|\leq K$）。注意 $L(y)$ 可能依赖于 $y$。

定义 $\tilde{f}(x,y)=f(x,y)-L(y)$。则 $\tilde{f}$ 与 $f$ 具有相同的单调性（关于 $x$），且

$$\lim_{x\to +\infty} \tilde{f}(x,y)=0\quad\text{对每个 }y\in E$$

---

**第 2 步：验证 $\tilde{f}$ 一致趋于零**。

需要将"逐点趋于零"提升为"一致趋于零"。由条件 2（$\int g$ 一致收敛），我们实际上需要证明 $\tilde{f}$ 满足 Dirichlet 条件——由于 $\tilde{f}$ 的极限逐点为 0，但一致趋于零需要进一步论证。

这里我们使用一种标准的处理方式：通过 $f$ 的单调性和有界性，可以证明 $\tilde{f}$ 关于 $y$ 一致趋于零。实际上，对任意 $\varepsilon>0$ 和任意 $y\in E$，存在 $X(y)$ 使得 $x>X(y)$ 时 $|\tilde{f}(x,y)|<\varepsilon$。但我们需要 $X$ 与 $y$ 无关。

严格论证如下：定义

$$\omega(x) = \sup_{y\in E} |\tilde{f}(x,y)| = \sup_{y\in E} |f(x,y)-L(y)|$$

由于 $\tilde{f}(x,\cdot)$ 关于 $x$ 单调且逐点趋于 0，且 $f$ 一致有界（$|f|\leq K$），故 $|\tilde{f}|\leq 2K$。由单调性，对每个 $y$，$|\tilde{f}(x,y)|$ 关于 $x$ 递减（若 $f$ 递减）或递增（若 $f$ 递增）。由 Dini 定理类型的论证可以证明 $\omega(x)\to 0$ 当 $x\to\infty$。以下使用更直接的估计：

由条件 1，$f$ 关于 $x$ 单调。不妨设 $f$ 关于 $x$ 递减（递增情形类似）。则 $L(y)=\lim_{x\to\infty} f(x,y)$ 满足 $L(y)\leq f(x,y)$ 对一切 $x$。于是 $\tilde{f}(x,y)=f(x,y)-L(y)\geq 0$ 且关于 $x$ 递减趋于 0。

现在关键观察：$\tilde{f}$ 是否一致趋于零取决于 $f$ 趋于 $L(y)$ 的速度是否关于 $y$ 一致。在 Abel 判别法的条件下，**$f$ 的一致有界性并不自动保证 $\tilde{f}$ 一致趋于零**。实际上，Abel 判别法的标准证明并不直接要求 $\tilde{f}$ 一致趋于零，而是使用 Dirichlet 判别法的一种变体——将 $\tilde{f}$ 的一致趋于零替换为：

---

**修正的证明**（标准教科书版本）：

**第 1 步：利用 $\int g$ 的一致收敛性**。

由条件 2，$G(y)=\int_a^{+\infty} g(x,y)\,dx$ 在 $E$ 上一致收敛。由 Cauchy 准则（定理 15.4），$\forall\varepsilon>0$，$\exists A_0\geq a$，$\forall A_2>A_1\geq A_0$，$\forall y\in E$：

$$\left|\int_{A_1}^{A_2} g(x,y)\,dx\right| < \varepsilon$$

---

**第 2 步：应用积分第二中值定理**。

对任意 $A_2>A_1\geq A_0$ 和任意 $y\in E$，由积分第二中值定理（引理 8.1），存在 $\xi=\xi(y)\in[A_1,A_2]$ 使得

$$\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx = f(A_1,y)\int_{A_1}^{\xi} g(x,y)\,dx + f(A_2,y)\int_{\xi}^{A_2} g(x,y)\,dx$$

---

**第 3 步：利用 $f$ 的有界性放缩**。

由条件 1，$|f(x,y)|\leq K$ 对所有 $x\geq a$ 和 $y\in E$ 成立。再由第 1 步，取 $\varepsilon' = \varepsilon/(2K)$，存在 $A_0$ 使得 $\left|\int_{A_1}^{A_2} g(x,y)\,dx\right| < \varepsilon/(2K)$ 对所有 $A_2>A_1\geq A_0$ 成立。于是：

$$\begin{aligned}
\left|\int_{A_1}^{A_2} f(x,y)g(x,y)\,dx\right|
&\leq |f(A_1,y)|\cdot\left|\int_{A_1}^{\xi} g\right| + |f(A_2,y)|\cdot\left|\int_{\xi}^{A_2} g\right| \\
&< K\cdot\frac{\varepsilon}{2K} + K\cdot\frac{\varepsilon}{2K} = \varepsilon
\end{aligned}$$

---

**第 4 步：由 Cauchy 准则得一致收敛**。

由 Cauchy 准则（定理 15.4），$\displaystyle F(y)=\int_a^{+\infty} f(x,y)g(x,y)\,dx$ 在 $E$ 上一致收敛。$\square$

---

**说明**：上述证明不需要将 Abel 判别法归化为 Dirichlet 判别法，而是直接使用积分第二中值定理和 Cauchy 准则。这是更简洁的标准证明路径。另一种证明路径是通过归化（$f=\tilde{f}+L$，其中 $\tilde{f}$ 部分用 Dirichlet 判别法），但需要额外验证 $\tilde{f}$ 的一致趋于零，更适合在 $f$ 的极限 $L(y)$ 与 $y$ 无关的简单情形使用。

---

**与单变量版本（定理 8.10）的对比**：

| 对比维度 | 单变量版本（定理 8.10） | 含参变量版本（定理 15.7） |
|----------|----------------------|-------------------------|
| 条件 1 | $f$ 单调有界 | $f$ 单调，**一致有界**（$|f|\leq K$ 对所有 $y$） |
| 条件 2 | $\int g$ 收敛 | $\int g$ 在 $E$ 上**一致收敛** |
| 结论 | $\int fg$ 收敛 | $\int fg$ 在 $E$ 上**一致收敛** |

---

#### 4.4.1 瑕积分情形的 Abel 判别法

**定理 15.7'（Abel 判别法——瑕积分情形）**：设 $x=b$ 是瑕点，$f(x,y),g(x,y)$ 定义在 $[a,b)\times E$ 上。若：

1. $f(x,y)$ 关于 $x$ 单调且一致有界（$|f(x,y)|\leq K$ 对一切 $x$ 和 $y$）；
2. $\displaystyle\int_a^b g(x,y)\,dx$ 在 $E$ 上一致收敛。

则瑕积分 $\displaystyle\int_a^b f(x,y)g(x,y)\,dx$ 在 $E$ 上一致收敛。

**证明思路**：与定理 15.7 的证明完全平行，在 $[b-\eta_1,b-\eta_2]$ 上应用积分第二中值定理和瑕积分版本的 Cauchy 准则。$\square$

---

### 4.5 两类判别法的内在联系与对比

**定理核心逻辑**：两个判别法都依赖积分第二中值定理将 $\int fg$ 的 Cauchy 条件转化为 $f$ 在端点处的值乘以 $\int g$ 在子区间上的积分。区别在于：

- **Dirichlet 判别法**：用 $f$ 的**小**（$f\to 0$）来控制，即使 $\int g$ 只有有界性而非收敛性。
- **Abel 判别法**：用 $\int g$ 的**小**（$\int g$ 的余项足够小）来控制，即使 $f$ 只有有界性而非趋于零。

| 对比维度 | Dirichlet 判别法（定理 15.6） | Abel 判别法（定理 15.7） |
|----------|-----------------------------|-------------------------|
| 对 $f$ 的条件 | 单调，一致趋于 0 | 单调，一致有界 |
| 对 $\int g$ 的条件 | 变上限积分一致有界 | 本身一致收敛 |
| 证明策略 | 第二中值定理 + Cauchy 准则 | 第二中值定理 + Cauchy 准则 |
| 条件"强弱" | $f$ 要求更强，$\int g$ 要求更弱 | $f$ 要求更弱，$\int g$ 要求更强 |
| 典型场景 | $g$ 振荡但本身积分不收敛（如 $\sin x$） | $g$ 积分已知一致收敛，$f$ 有界变化 |

**实际应用中的选择决策**：

```
问题：判断 ∫ f(x,y)g(x,y)dx 的一致收敛性
│
├─ 能否找到与 y 无关的 M(x) 使 |fg| ≤ M(x) 且 ∫M(x)dx 收敛？
│  └─ 是 → M-判别法（定理 15.5），得到绝对一致收敛
│
├─ 不满足 M-判别法，但 f 单调趋于零（关于 y 一致）且 ∫g 一致有界？
│  └─ 是 → Dirichlet 判别法（定理 15.6）
│
├─ 不满足 M-判别法，但 f 单调有界（一致）且 ∫g 一致收敛？
│  └─ 是 → Abel 判别法（定理 15.7）
│
└─ 以上均不满足 → 需用定义或其他方法判断
```

### 4.6 三种判别法的系统对比

| 对比维度 | M-判别法（定理 15.5） | Dirichlet 判别法（定理 15.6） | Abel 判别法（定理 15.7） |
|----------|---------------------|-----------------------------|-------------------------|
| 对 $f$ 的条件 | $|f(x,y)|\leq M(x)$（被控制） | $f$ 单调，$\lim f=0$ **关于 $y$ 一致** | $f$ 单调，**一致有界** $|f|\leq K$ |
| 对 $\int g$ 的条件 | 无（直接控制整个被积函数） | $\int_a^A g\,dx$ **一致有界** | $\int_a^{+\infty} g\,dx$ **一致收敛** |
| 结论收敛类型 | **绝对一致收敛** | 条件收敛可能，一致收敛 | 条件收敛可能，一致收敛 |
| 证明核心工具 | 绝对值不等式 + 优函数 | 积分第二中值定理 + Cauchy 准则 | 积分第二中值定理 + Cauchy 准则 |
| 适用典型结构 | 衰减型，无振荡，$(x,y)$ 分离 | 单调衰减因子 $\times$ 振荡因子 | 单调有界因子 $\times$ 已知一致收敛积分 |
| 对 $y$ 的要求 | 优函数 $M(x)$ 与 $y$ 无关 | 一致性条件（关于 $y$）自然体现 | 一致性条件（关于 $y$）自然体现 |

**互补关系**：
- M-判别法覆盖**绝对可积**型——只要 $|f(x,y)|\leq M(x)$ 且 $\int M$ 收敛，积分就绝对一致收敛。
- Dirichlet 判别法覆盖**条件收敛的衰减振荡型**——$f$ 衰减到零，$g$ 振荡但变上限积分有界，通过正负抵消实现一致收敛。
- Abel 判别法覆盖**在已知一致收敛积分上乘以有界因子**——扩展已知结果。
- 三者形成完整的判别体系，覆盖了从绝对收敛到条件收敛的各类含参变量反常积分。

---

## 5. 例题精讲

### 例 4-1：Dirichlet 判别法的基础应用

证明含参变量反常积分

$$\Phi(\alpha) = \int_1^{+\infty} \frac{\sin(\alpha x)}{x}\,dx,\quad \alpha\in[\alpha_0,+\infty)\;(\alpha_0>0)$$

在 $[\alpha_0,+\infty)$ 上一致收敛。

---

**解**：

**第 1 步：构造乘积分解**。

取 $f(x,\alpha)=\dfrac{1}{x}$，$g(x,\alpha)=\sin(\alpha x)$。则 $\dfrac{\sin(\alpha x)}{x}=f(x,\alpha)g(x,\alpha)$。

---

**第 2 步：验证 $f$ 的条件（单调且一致趋于零）**。

$f(x,\alpha)=1/x$ 对每个 $\alpha$ 关于 $x$ 单调递减，且

$$\lim_{x\to+\infty} \frac{1}{x}=0$$

该极限与 $\alpha$ 无关，故关于 $\alpha\in[\alpha_0,+\infty)$ 一致成立。✓

---

**第 3 步：验证 $\int g$ 的一致有界性**。

计算变上限积分：

$$\left|\int_1^A \sin(\alpha x)\,dx\right| = \left|\frac{\cos\alpha - \cos(\alpha A)}{\alpha}\right| \leq \frac{2}{\alpha} \leq \frac{2}{\alpha_0}, \quad \forall A\geq 1,\;\forall \alpha\geq\alpha_0$$

因此 $\left|\int_1^A \sin(\alpha x)\,dx\right| \leq \dfrac{2}{\alpha_0}$，一致有界（取 $M=2/\alpha_0$）。✓

---

**第 4 步：应用定理 15.6**。

由 Dirichlet 判别法（定理 15.6），$\Phi(\alpha)$ 在 $[\alpha_0,+\infty)$ 上一致收敛。✓

---

**说明**：注意 $f(x,\alpha)=1/x$ 也与 $\alpha$ 无关，因此条件 1 的一致性是自动满足的。一般情况下，$f(x,y)$ 应选择为**不含 $y$ 或含 $y$ 但关于 $y$ 一致趋于零**的函数。

---

### 例 4-2：内闭一致收敛的经典例子

基于例 4-1 的结论，讨论 $\Phi(\alpha)=\displaystyle\int_1^{+\infty}\frac{\sin(\alpha x)}{x}\,dx$ 在 $(0,+\infty)$ 上的收敛性。

---

**解**：

**(1) 内闭一致收敛性**：

对 $(0,+\infty)$ 的任意紧子集 $K$，由于 $K$ 是有界闭集且不包含 $0$（因为 $0$ 是 $(0,+\infty)$ 的边界点），存在 $\alpha_0>0$ 使得 $K\subset[\alpha_0,+\infty)$。由例 4-1，$\Phi(\alpha)$ 在 $[\alpha_0,+\infty)$ 上一致收敛，从而在 $K$ 上一致收敛。

由定义 15.7，$\Phi(\alpha)$ 在 $(0,+\infty)$ 上**内闭一致收敛**。

---

**(2) 整体一致收敛性的否定**：

我们证明 $\Phi(\alpha)$ 在 $(0,+\infty)$ 上**不**整体一致收敛。

采用 Cauchy 准则的否定形式。考虑 $\alpha_n = \dfrac{1}{n}$（趋于 $0^+$），取 $A_n' = n$，$A_n'' = 2n$。计算：

$$\begin{aligned}
\int_{A_n'}^{A_n''} \frac{\sin(\alpha_n x)}{x}\,dx
&= \int_n^{2n} \frac{\sin(x/n)}{x}\,dx
\end{aligned}$$

换元 $t=x/n$，则 $x=nt$，$dx=n\,dt$，积分限 $t\in[1,2]$：

$$\int_n^{2n} \frac{\sin(x/n)}{x}\,dx = \int_1^2 \frac{\sin t}{nt}\cdot n\,dt = \int_1^2 \frac{\sin t}{t}\,dt > 0$$

记 $C=\int_1^2 \frac{\sin t}{t}\,dt > 0$（因为 $\sin t>0$ 在 $[1,2]\subset(0,\pi)$ 上）。于是：

$$\left|\int_{A_n'}^{A_n''} \frac{\sin(\alpha_n x)}{x}\,dx\right| = C > 0$$

取 $\varepsilon_0=C/2>0$，则对任意 $A_0$，取 $n>A_0$ 使得 $A_n'=n>A_0$，有 $|\int_{A_n'}^{A_n''}| = C > \varepsilon_0$。由 Cauchy 准则（定理 15.4）的否定形式，$\Phi(\alpha)$ 在 $(0,+\infty)$ 上不一致收敛。

---

**(3) 直观解释**：

当 $\alpha\to 0^+$ 时，$\sin(\alpha x)$ 的振荡周期 $T=2\pi/\alpha \to +\infty$。对固定的 $A$，余项 $\int_A^{+\infty} \frac{\sin(\alpha x)}{x}\,dx$ 的衰减速度随 $\alpha\to 0^+$ 趋近于 $\int_A^{+\infty} \frac{0}{x}\,dx = 0$ 的极限过程——但 $\alpha$ 越小，需要越大的 $A$ 才能使 $\sin(\alpha x)$ 完成一个振荡周期并实现正负抵消。因此无法找到对所有 $\alpha>0$ 同时有效的 $A_0$。

---

### 例 4-3：Abel 判别法的基础应用

已知反常积分 $\displaystyle\int_1^{+\infty} \frac{\sin x}{x}\,dx$ 收敛（Dirichlet 积分），考虑含参变量反常积分

$$F(\alpha) = \int_1^{+\infty} \frac{\sin x}{x+\alpha}\,dx,\quad \alpha\in[0,+\infty)$$

利用 Abel 判别法证明 $F(\alpha)$ 在 $[0,+\infty)$ 上一致收敛。

---

**解**：

**第 1 步：构造乘积分解**。

取 $f(x,\alpha)=\dfrac{x}{x+\alpha}$，$g(x,\alpha)=\dfrac{\sin x}{x}$。则：

$$\frac{\sin x}{x+\alpha} = \frac{x}{x+\alpha}\cdot\frac{\sin x}{x} = f(x,\alpha)g(x,\alpha)$$

---

**第 2 步：验证 $f$ 单调一致有界**。

对 $\alpha\geq 0$，$x\geq 1$：

$$f(x,\alpha)=\frac{x}{x+\alpha}=1-\frac{\alpha}{x+\alpha}$$

计算偏导数：

$$\frac{\partial f}{\partial x} = \frac{\alpha}{(x+\alpha)^2} \geq 0, \quad \forall \alpha\geq 0$$

因此 $f(x,\alpha)$ 关于 $x$ 单调递增（$\alpha=0$ 时 $f\equiv 1$，也是单调的）。

有界性：对 $x\geq 1$，$\alpha\geq 0$，

$$0 < f(x,\alpha) = \frac{x}{x+\alpha} \leq 1$$

因此 $|f(x,\alpha)|\leq 1=K$，一致有界。✓

---

**第 3 步：验证 $\int g$ 一致收敛**。

$g(x,\alpha)=\dfrac{\sin x}{x}$ 与 $\alpha$ 无关。已知 $\int_1^{+\infty} \frac{\sin x}{x}\,dx$ 收敛（Dirichlet 积分，由单变量 Dirichlet 判别法，定理 8.9），且收敛性与 $\alpha$ 无关，故关于 $\alpha\in[0,+\infty)$ 一致收敛。✓

---

**第 4 步：应用定理 15.7**。

由 Abel 判别法（定理 15.7），$F(\alpha)=\displaystyle\int_1^{+\infty} \frac{\sin x}{x+\alpha}\,dx$ 在 $[0,+\infty)$ 上**一致收敛**。✓

---

**说明**：本例中 $\int g$ 的收敛性与 $\alpha$ 无关，因此"一致收敛"是自动成立的。更一般的情形下，$\int g$ 的一致收敛性需要单独验证。

---

### 例 4-4：Abel 判别法处理指数衰减因子

考虑含参变量反常积分

$$G(\alpha)=\int_0^{+\infty} e^{-\alpha x}\frac{\sin x}{x}\,dx,\quad \alpha\geq 0$$

分别讨论 $\alpha>0$ 和 $\alpha=0$ 的情形，判断 $G(\alpha)$ 在 $[0,+\infty)$ 上是否一致收敛。

---

**解**：

**(1) $\alpha>0$ 的情形——M-判别法**：

对固定的 $\alpha_0>0$，考虑 $\alpha\in[\alpha_0,+\infty)$：

$$\left|e^{-\alpha x}\frac{\sin x}{x}\right| \leq e^{-\alpha x} \leq e^{-\alpha_0 x}$$

而 $\displaystyle\int_0^{+\infty} e^{-\alpha_0 x}\,dx = \frac{1}{\alpha_0} < \infty$。由 M-判别法（定理 15.5），$G(\alpha)$ 在 $[\alpha_0,+\infty)$ 上**绝对一致收敛**。

---

**(2) $\alpha=0$ 的情形——M-判别法失效**：

当 $\alpha=0$ 时，$G(0)=\int_0^{+\infty} \frac{\sin x}{x}\,dx$。这时 $|e^{0}\cdot\sin x/x| = |\sin x/x| \leq 1/x$，但 $\int_0^{+\infty} 1/x\,dx$ 发散——M-判别法不适用。但 $G(0)$ 本身是收敛的 Dirichlet 积分。需要证明当 $\alpha\to 0^+$ 时，$G(\alpha)$ 在 $\alpha=0$ 附近也一致收敛。

**(3) 综合判断——Abel 判别法覆盖 $\alpha=0$**：

考虑 $\alpha\in[0,+\infty)$。取 $f(x,\alpha)=e^{-\alpha x}$，$g(x,\alpha)=\dfrac{\sin x}{x}$。

验证 Abel 判别法的条件：

- **$f$ 单调有界**：$f(x,\alpha)=e^{-\alpha x}$ 关于 $x$ 单调递减（$\alpha\geq 0$ 时），且 $0<e^{-\alpha x}\leq 1$，一致有界（$K=1$）。✓
- **$\int g$ 一致收敛**：$\displaystyle\int_0^{+\infty} \frac{\sin x}{x}\,dx$ 收敛且与 $\alpha$ 无关，故关于 $\alpha\in[0,+\infty)$ 一致收敛。✓

由 Abel 判别法（定理 15.7），$G(\alpha)$ 在 $[0,+\infty)$ 上**一致收敛**。

---

**(4) 讨论——两种判别法的分工**：

- M-判别法在 $\alpha>0$ 的开区间上给出了**绝对一致收敛**的更强结论，但无法处理 $\alpha=0$ 处的条件收敛。
- Abel 判别法统一处理了整个区间 $[0,+\infty)$，包括 $\alpha=0$ 处的条件收敛情形。
- 两种判别法**互补**：需要绝对收敛性时用 M-判别法，需要覆盖条件收敛边界时用 Abel 判别法。

---

### 例 4-5：三种判别法的综合运用

判断含参变量反常积分

$$H(\alpha) = \int_1^{+\infty} \frac{\sin(\alpha x)}{x(1+e^{\alpha x})}\,dx,\quad \alpha>0$$

在 $(0,+\infty)$ 上的一致收敛性。

---

**分析**：这个积分同时包含振荡因子 $\sin(\alpha x)$ 和指数衰减因子 $1/(1+e^{\alpha x})$，以及代数衰减 $1/x$。三种判别法都有可能——取决于参数区域。

---

**解**：

**(1) 对 $\alpha\in[\alpha_0,+\infty)$（$\alpha_0>0$）——M-判别法**：

$$\left|\frac{\sin(\alpha x)}{x(1+e^{\alpha x})}\right| \leq \frac{1}{x(1+e^{\alpha_0 x})}$$

注意 $\dfrac{1}{x(1+e^{\alpha_0 x})} \leq \dfrac{1}{x}$ 但 $\int_1^{+\infty} 1/x\,dx$ 发散——需用更精细的优函数。

实际上，$\dfrac{1}{1+e^{\alpha_0 x}} \leq e^{-\alpha_0 x/2}$（对充分大的 $x$）。更精确地：$1+e^{\alpha_0 x} \geq e^{\alpha_0 x}$，故

$$\frac{1}{x(1+e^{\alpha_0 x})} \leq \frac{e^{-\alpha_0 x}}{x} \leq e^{-\alpha_0 x/2}$$

对充分大的 $x$ 成立，而 $\int_1^{+\infty} e^{-\alpha_0 x/2}\,dx$ 收敛。由 M-判别法，$H(\alpha)$ 在 $[\alpha_0,+\infty)$ 上绝对一致收敛。

**(2) 对 $\alpha\in(0,\alpha_0)$（靠近 $0$ 的区域）——Dirichlet 判别法**：

在 $\alpha\in(0,\alpha_0)$ 上，因子 $1/(1+e^{\alpha x})$ 的衰减随 $\alpha\to 0^+$ 变慢。但可取：

$$f(x,\alpha) = \frac{1}{x(1+e^{\alpha x})}, \quad g(x,\alpha)=\sin(\alpha x)$$

$f(x,\alpha)$ 关于 $x$ 单调递减，且 $\lim_{x\to\infty} f(x,\alpha)=0$ 是否关于 $\alpha\in(0,\alpha_0)$ 一致？需验证：

$$|f(x,\alpha)| = \frac{1}{x(1+e^{\alpha x})} \leq \frac{1}{x}$$

但 $\frac{1}{x}\to 0$ 与 $\alpha$ 无关，故一致趋于零。✓

同时 $\left|\int_1^A \sin(\alpha x)\,dx\right| \leq 2/\alpha_0$ 对 $\alpha\in[\alpha_0,+\infty)$ 成立——但这里 $\alpha$ 可以任意小！实际上 $\left|\int_1^A \sin(\alpha x)\,dx\right| \leq 2/\alpha$ 在 $\alpha\to 0^+$ 时无界。因此**无法直接**对 $(0,\alpha_0)$ 用 Dirichlet 判别法。

这里实际上 $\int_1^A \sin(\alpha x)\,dx$ 在 $\alpha\to 0^+$ 时不是一致有界的（因为 $2/\alpha\to +\infty$）。所以 Dirichlet 判别法在包含 $0$ 的区间上也不适用。

**(3) 最终结论**：$H(\alpha)$ 在 $(0,+\infty)$ 上**不**一致收敛；但在任意 $[\alpha_0,+\infty)$（$\alpha_0>0$）上一致收敛。即 $H(\alpha)$ 在 $(0,+\infty)$ 上**内闭一致收敛**。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 含参变量 Dirichlet/Abel 判别法与单变量版本完全相同，无需修改 | 关键区别在于含参版本要求"关于参数 $y$ 一致"——$f$ 的一致趋于零（Dirichlet）或 $\int g$ 的一致收敛（Abel）。这个"一致"条件是含参变量积分的核心关切，不能省略。 |
| Dirichlet 判别法只需 $f$ 单调趋于零（逐点）即可 | 需要 $\lim_{x\to\infty}f(x,y)=0$ **关于 $y$ 一致成立**。仅逐点趋于零不够——例如 $f(x,y)=e^{-xy}$ 在 $y>0$ 时逐点趋于零，但在 $(0,+\infty)$ 上不一致趋于零（$y\to 0^+$ 时衰减变慢）。 |
| 内闭一致收敛就是整体一致收敛 | 两者不等价。内闭一致收敛 $\not\Rightarrow$ 整体一致收敛。$\Phi(\alpha)=\int_1^{+\infty}\sin(\alpha x)/x\,dx$ 在 $(0,+\infty)$ 上内闭一致收敛但不整体一致收敛。区别在于 $A_0$ 是否可以依赖于紧子集 $K$。 |
| Abel 判别法的证明必须归化为 Dirichlet 判别法 | 归化法是其中一种证明路径。更简洁直接的证明是直接用积分第二中值定理 + Cauchy 准则，通过 $f$ 的有界性和 $\int g$ 的 Cauchy 条件一步到位（如定理 15.7 的证明所示）。 |
| M-判别法是"最强"的判别法，能用 M-判别法就不用 Dirichlet/Abel | M-判别法确实给出了最强的结论（绝对一致收敛），但适用范围最窄。对条件收敛的振荡型积分，M-判别法失效。三种判别法的关系是互补而非包含。 |
| 在 $(0,+\infty)$ 上内闭一致收敛足以说明连续性 | 正确。连续性是一个局部概念，只需在每个点的某邻域内一致收敛即可保证连续。内闭一致收敛正是这个意义上的"足够"。 |
| 应用 Dirichlet 判别法时，$f$ 必须单调递减 | 递增也同样适用——积分第二中值定理对任意单调函数均成立。只需注意 $f$ 的符号和端点值的选择即可。变通形式：递增时第二中值定理的端点改为 $f(b)$。 |
| 构造乘积分解时，必须将 $f$ 取为"衰减因子"、$g$ 取为"振荡因子" | 分配方式不唯一。关键是要使 $f$ 满足单调性条件，$\int g$ 满足有界/收敛条件。有时需要灵活调整分解方式。 |

### 检查点

- [ ] 能否准确写出含参变量积分 Dirichlet 判别法（定理 15.6）的条件与结论，并指出与单变量版本（定理 8.9）的核心区别？
- [ ] 能否给出定理 15.6 的完整证明（积分第二中值定理 + Cauchy 准则）？
- [ ] 能否准确写出含参变量积分 Abel 判别法（定理 15.7）的条件与结论？
- [ ] 能否写出内闭一致收敛的定义（定义 15.7），并举出一个内闭一致收敛但不整体一致收敛的例子？
- [ ] 能否在 30 秒内说出三种判别法（M、Dirichlet、Abel）对 $f$ 和 $\int g$ 的条件差异？
- [ ] 能否独立完成例 4-1（$\int_1^{+\infty}\sin(\alpha x)/x\,dx$ 在 $[\alpha_0,+\infty)$ 上一致收敛的完整验证）？
- [ ] 能否独立完成例 4-3（$\int_1^{+\infty}\sin x/(x+\alpha)\,dx$ 用 Abel 判别法）？
- [ ] 能否综合运用三种判别法处理例 4-4（$e^{-\alpha x}\sin x/x$）？
- [ ] 能否说出为什么例 4-2 中 $\Phi(\alpha)$ 在 $(0,+\infty)$ 上不一致收敛，但在每个紧子集上一致收敛？
- [ ] 对于形如 $\int f(x,y)g(x,y)dx$ 的积分，能否根据 $f$ 和 $g$ 的性质快速选择最合适的判别法（M/Dirichlet/Abel）？

---

## 练习题

### 基础巩固

**1.**（定理理解）判断下列说法是否正确，并说明理由：

(a) 若 $f(x,y)$ 在 $[a,+\infty)\times E$ 上满足 $|f(x,y)|\leq M(x)$ 且 $\int_a^{+\infty} M(x)\,dx$ 收敛，则 $\int_a^{+\infty} f(x,y)\,dx$ 在 $E$ 上一致收敛。

(b) 若对每个 $y\in E$，$f(x,y)$ 关于 $x$ 单调，$\lim_{x\to\infty}f(x,y)=0$，且 $|\int_a^A g(x,y)\,dx|\leq M$ 对所有 $A,y$ 成立，则 $\int_a^{+\infty} f(x,y)g(x,y)\,dx$ 在 $E$ 上一致收敛。

(c) 若 $F(y)=\int_a^{+\infty} f(x,y)\,dx$ 在 $(0,1)$ 上内闭一致收敛，则 $F(y)$ 在 $(0,1)$ 上一致收敛。

<details><summary>参考答案</summary>

**(a) 正确**。这正是 M-判别法（定理 15.5）的条件。$M$-判别法要求优函数 $M(x)$ 与 $y$ 无关且自身可积。

**(b) 正确**。由 Dirichlet 判别法（定理 15.6），条件是：$f$ 单调且 $\lim f=0$ 关于 $y$ 一致；$\int_a^A g$ 一致有界。题目中第二个条件已明确给出"一致有界"；第一个条件中 $\lim f=0$ 对每个 $y$ 成立，但需要"关于 $y$ 一致趋于零"——若题目只说了"对每个 $y$"而未涉及"一致"，则该条件不完全等价于定理 15.6。需要额外验证一致性。

**(c) 错误**。内闭一致收敛 $\not\Rightarrow$ 整体一致收敛。反例：$\int_1^{+\infty} \frac{\sin(\alpha x)}{x}\,dx$ 在 $(0,+\infty)$ 上内闭一致收敛但不整体一致收敛。同理，在 $(0,1)$ 上也可能存在类似反例。

</details>

---

**2.**（Dirichlet 判别法验证）用 Dirichlet 判别法证明 $\displaystyle\int_1^{+\infty} \frac{\cos(\alpha x)}{x^{3/2}}\,dx$ 在 $\alpha\in[\alpha_0,+\infty)$（$\alpha_0>0$）上一致收敛。

<details><summary>参考答案</summary>

取 $f(x,\alpha)=1/x^{3/2}$，$g(x,\alpha)=\cos(\alpha x)$。

**验证条件 1**：$f(x,\alpha)=1/x^{3/2}$ 关于 $x$ 单调递减，$\lim_{x\to+\infty}1/x^{3/2}=0$ 与 $\alpha$ 无关，故关于 $\alpha$ 一致趋于 0。✓

**验证条件 2**：

$$\left|\int_1^A \cos(\alpha x)\,dx\right| = \left|\frac{\sin(\alpha A)-\sin\alpha}{\alpha}\right| \leq \frac{2}{\alpha} \leq \frac{2}{\alpha_0}$$

因此 $\left|\int_1^A \cos(\alpha x)\,dx\right| \leq 2/\alpha_0$，一致有界。✓

由定理 15.6（Dirichlet 判别法），原积分在 $[\alpha_0,+\infty)$ 上一致收敛。

**进一步判断**：由于 $| \cos(\alpha x)/x^{3/2} | \leq 1/x^{3/2}$ 且 $\int_1^{+\infty} 1/x^{3/2}dx$ 收敛（$p=3/2>1$），实际上该积分**绝对一致收敛**（M-判别法也适用）。本题说明：当两种判别法都适用时，M-判别法给出更强的结论（绝对收敛），但 Dirichlet 判别法的验证同样有效。

</details>

---

**3.**（Abel 判别法应用）已知 $\displaystyle\int_1^{+\infty} \frac{\cos x}{x^{3/2}}\,dx$ 收敛（练习 2 结论），证明 $\displaystyle\int_1^{+\infty} \frac{x+2}{x+1}\cdot\frac{\cos x}{x^{3/2}}\,dx$ 在 $[0,+\infty)$ 上一致收敛（该积分中"参数"以常数的形式出现，此为单变量积分情形；若引入含参因子则得到含参版本）。

<details><summary>参考答案</summary>

将原积分写为 $\int_1^{+\infty} f(x)g(x)\,dx$。取 $f(x)=\frac{x+2}{x+1}$，$g(x)=\frac{\cos x}{x^{3/2}}$。

其实该问题可以扩展为含参版本：$F(\alpha)=\int_1^{+\infty} \frac{x+\alpha}{x+1}\cdot\frac{\cos x}{x^{3/2}}\,dx$，$\alpha\in[0,+\infty)$。

取 $f(x,\alpha)=\frac{x+\alpha}{x+1}$，$g(x,\alpha)=\frac{\cos x}{x^{3/2}}$。

**验证 $f$ 单调一致有界**：

$$f(x,\alpha)=\frac{x+\alpha}{x+1}=1+\frac{\alpha-1}{x+1}$$

当 $\alpha\geq 0$ 时，对 $x\geq 1$：

- 若 $\alpha\geq 1$：$\frac{\partial f}{\partial x}=-\frac{\alpha-1}{(x+1)^2}\leq 0$，$f$ 递减
- 若 $\alpha\in[0,1)$：$\frac{\partial f}{\partial x}=-\frac{\alpha-1}{(x+1)^2}>0$，$f$ 递增

不管哪种情况，$f$ 关于 $x$ 单调。且有界性：$|f(x,\alpha)|\leq \max\{1,\alpha\}$（对 $x\geq 1$），所以一致有界。✓

**验证 $\int g$ 一致收敛**：

$\int_1^{+\infty} \frac{\cos x}{x^{3/2}}\,dx$ 收敛且与 $\alpha$ 无关，故关于 $\alpha$ 一致收敛。✓

由定理 15.7（Abel 判别法），$F(\alpha)$ 在 $[0,+\infty)$ 上一致收敛。

</details>

---

### 迁移应用

**4.**（内闭一致收敛的判断）设

$$F(\alpha)=\int_1^{+\infty} \frac{\sin(\alpha x)}{x}\,dx,\quad \alpha>0$$

证明 $F(\alpha)$ 在 $(0,+\infty)$ 上内闭一致收敛，但不一致收敛。（参考例 4-1 和例 4-2）

<details><summary>参考答案</summary>

**内闭一致收敛**：任取紧子集 $K\subset(0,+\infty)$。由于 $K$ 有界闭且不包含 0，存在 $\alpha_0>0$ 使 $K\subset[\alpha_0,+\infty)$。由例 4-1，$F(\alpha)$ 在 $[\alpha_0,+\infty)$ 上一致收敛，因此在 $K$ 上一致收敛。由定义 15.7，$F(\alpha)$ 在 $(0,+\infty)$ 上内闭一致收敛。

**不一致收敛**：由例 4-2，取 $\alpha_n=1/n$，$A_n'=n$，$A_n''=2n$，计算得 $\int_{A_n'}^{A_n''} \frac{\sin(x/n)}{x}dx = \int_1^2 \frac{\sin t}{t}dt = C>0$，不满足 Cauchy 准则。故不一致收敛。

</details>

---

**5.**（综合运用三种判别法）判断下列含参变量反常积分在指定参数集上的一致收敛性，并说明所使用的判别法。

(1) $F(\alpha)=\displaystyle\int_0^{+\infty} e^{-\alpha x}\cos x\,dx,\quad \alpha\in[0,+\infty)$

(2) $G(\alpha)=\displaystyle\int_1^{+\infty} \frac{\sin x}{x^{\alpha}}\,dx,\quad \alpha\in(1,+\infty)$

(3) $H(\alpha)=\displaystyle\int_0^{+\infty} \frac{\sin(\alpha x)}{x}\,e^{-x}\,dx,\quad \alpha\in[0,+\infty)$

<details><summary>参考答案</summary>

**(1)** $F(\alpha)=\int_0^{+\infty} e^{-\alpha x}\cos x\,dx$，$\alpha\in[0,+\infty)$。

**M-判别法**：$|e^{-\alpha x}\cos x|\leq e^{-\alpha x}$。但 $\int_0^{+\infty} e^{-\alpha x}dx = 1/\alpha$ 对 $\alpha>0$ 收敛，在 $\alpha=0$ 时发散。故 M-判别法不能直接覆盖整个 $[0,+\infty)$。

对 $\alpha\in[\alpha_0,+\infty)$（$\alpha_0>0$）：$|e^{-\alpha x}\cos x|\leq e^{-\alpha_0 x}$，$\int_0^{+\infty} e^{-\alpha_0 x}dx$ 收敛 → M-判别法得一致收敛。

对 $\alpha\in[0,+\infty)$ 整体：取 $f(x,\alpha)=e^{-\alpha x}$（关于 $x$ 单调递减，$|f|\leq 1$），$g(x,\alpha)=\cos x$。但 $\int_0^{+\infty} \cos x\,dx$ 发散（振荡），Abel 判别法不适用。实际上一开始对 $\alpha\in[0,+\infty)$ 直接计算：$F(\alpha)=\int_0^{+\infty} e^{-\alpha x}\cos x\,dx = \frac{\alpha}{\alpha^2+1}$（Laplace 变换）。当 $\alpha=0$ 时 $F(0)=\int_0^{+\infty} \cos x\,dx$ 发散！因此 $F(\alpha)$ 在 $\alpha=0$ 处发散，故在 $[0,+\infty)$ 上不是逐点收敛的，更不用说一致收敛。

本题提醒：使用判别法前需先确保积分对每个参数值都**逐点收敛**。

**(2)** $G(\alpha)=\int_1^{+\infty} \frac{\sin x}{x^{\alpha}}\,dx$，$\alpha\in(1,+\infty)$。

$|\sin x/x^{\alpha}|\leq 1/x^{\alpha}$。对 $\alpha\in(1,+\infty)$，$1/x^{\alpha}\leq 1/x$ 但 $\int_1^{+\infty}1/x\,dx$ 发散。需注意：$\int_1^{+\infty}1/x^{\alpha}\,dx$ 在 $\alpha>1$ 时收敛。取 $M(x)=1/x$ 不行，取 $M(x)=1/x^{\alpha_0}$（其中 $\alpha_0=\inf\{\alpha:\alpha\in(1,+\infty)\}=1$），但 $\alpha=1$ 不在区间内，所以对 $\alpha\in[\alpha_0,+\infty)$（$\alpha_0>1$），$| \sin x/x^{\alpha} |\leq 1/x^{\alpha_0}$，$\int_1^{+\infty}1/x^{\alpha_0}dx$ 收敛。

所以 $G(\alpha)$ 在 $[\alpha_0,+\infty)$（$\alpha_0>1$）上由 M-判别法得一致收敛。在 $(1,+\infty)$ 上内闭一致收敛但非整体一致收敛（当 $\alpha\to 1^+$ 时收敛速度变慢）。

**(3)** $H(\alpha)=\int_0^{+\infty} \frac{\sin(\alpha x)}{x}e^{-x}\,dx$，$\alpha\in[0,+\infty)$。

当 $\alpha=0$ 时，被积函数为 0，积分收敛于 0。

当 $\alpha>0$ 时，可将原式写为 $\int_0^{+\infty} \frac{\sin(\alpha x)}{x}e^{-x}\,dx$。

**方法一（Abel 判别法）**：取 $f(x,\alpha)=e^{-x}$，$g(x,\alpha)=\frac{\sin(\alpha x)}{x}$。$f$ 关于 $x$ 单调递减，$|f|\leq 1$ 有界。$\int_0^{+\infty} \frac{\sin(\alpha x)}{x}\,dx = \frac{\pi}{2}$（对 $\alpha>0$ 收敛，且与 $\alpha$ 无关→一致收敛）。由 Abel 判别法，$H(\alpha)$ 在 $[\alpha_0,+\infty)$ 上一致收敛。

**方法二（M-判别法）**：$|\frac{\sin(\alpha x)}{x}e^{-x}|\leq e^{-x}$（因为 $|\sin(\alpha x)/x|\leq \alpha$ 但对固定 $\alpha$ 无统一上界...实际上需更精细估计）。当 $\alpha$ 有界时也可用 M-判别法。

**综合**：$H(\alpha)$ 在 $[0,+\infty)$ 上一致收敛。

</details>

---

**6.**（含参变量瑕积分的 Dirichlet 判别法）考虑含参变量瑕积分

$$F(y)=\int_0^1 \frac{\sin(y/x)}{\sqrt{x}}\,dx,\quad y>0$$

证明 $F(y)$ 在 $[\delta,+\infty)$（$\delta>0$）上一致收敛。

**提示**：瑕点在 $x=0$，换元 $t=1/x$ 可将其转化为无穷限积分后再用 Dirichlet 判别法。

<details><summary>参考答案</summary>

**第 1 步：换元转化为无穷限积分**。

令 $t=1/x$，则 $x=1/t$，$dx=-dt/t^2$。当 $x\in(0,1]$ 时 $t\in[1,+\infty)$。

$$F(y)=\int_0^1 \frac{\sin(y/x)}{\sqrt{x}}\,dx = \int_{+\infty}^1 \frac{\sin(yt)}{\sqrt{1/t}}\cdot\left(-\frac{1}{t^2}\right)dt = \int_1^{+\infty} \frac{\sin(yt)}{t^{3/2}}\,dt$$

**第 2 步：应用 Dirichlet 判别法**。

取 $f(t,y)=1/t^{3/2}$，$g(t,y)=\sin(yt)$。

$f(t,y)=1/t^{3/2}$ 关于 $t$ 单调递减，$\lim_{t\to\infty}1/t^{3/2}=0$ 与 $y$ 无关。

$$\left|\int_1^A \sin(yt)\,dt\right| = \left|\frac{\cos y - \cos(yA)}{y}\right| \leq \frac{2}{y} \leq \frac{2}{\delta}$$

因此 $\left|\int_1^A \sin(yt)\,dt\right|\leq 2/\delta$ 在 $y\geq\delta$ 上一致有界。

由定理 15.6（Dirichlet 判别法），$F(y)$ 在 $[\delta,+\infty)$ 上一致收敛。

</details>
