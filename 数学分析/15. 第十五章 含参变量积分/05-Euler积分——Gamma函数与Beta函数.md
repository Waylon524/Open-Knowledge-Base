# 05. Euler积分——Gamma函数与Beta函数

> 所属章节：第十五章 含参变量积分 | 文件序号：05 | 难度：核心
> 常见混淆点：1) Gamma函数Γ(s)作为含参变量反常积分，其收敛域为s>0，初学者容易忽略在s→0⁺时∫₀¹ x^{s-1}e^{-x}dx的发散行为；2) B-Γ关系B(p,q)=Γ(p)Γ(q)/Γ(p+q)中分母是Γ(p+q)而非Γ(p)Γ(q)，量纲判断法可避免混淆；3) 余元公式Γ(s)Γ(1-s)=π/sin(πs)仅对0<s<1直接成立，解析延拓后可用于更广范围但需说明；4) 含参变量反常积分的可微性定理（定理15.8）要求∫∂f/∂y一致收敛（而非仅f本身一致收敛），初学者容易混淆条件；5) Gamma函数逐阶可导时，ln x的n次幂在x→0⁺和x→∞处的可积性需要分别验证，不能一概而论

---

## 1. 学习目标与先修前置

### 学习目标
- 掌握 Gamma 函数 $\Gamma(s) = \displaystyle\int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 和 Beta 函数 $B(p,q) = \displaystyle\int_0^1 x^{p-1}(1-x)^{q-1}\,dx$ 的严格定义、收敛域及基本性质
- 掌握 B-Γ 关系（定理 8.11）及其证明，并能熟练推导 $\Gamma(1/2)=\sqrt{\pi}$
- **掌握含参变量反常积分的可微性定理（定理 15.8）**，理解其与常义版本（定理 15.2）的关键区别，并能独立完成证明
- 能应用定理 15.8 严格证明 Gamma 函数在 $(0,+\infty)$ 上的逐阶可微性，并写出 $\Gamma'(s)$ 和 $\Gamma''(s)$ 的显式表达式
- 掌握余元公式（定理 15.9）$\Gamma(s)\Gamma(1-s)=\pi/\sin(\pi s)$ 及其证明
- 能灵活运用 Euler 积分（Gamma 函数 + Beta 函数）系统计算三类标准型积分：有理函数型、三角型、以及利用余元公式的积分
- 建立 Euler 积分作为含参变量反常积分理论之"旗舰应用"的整体认识

### 先修知识

| 先修文件 | 所需内容 |
|----------|----------|
| **文件 03（本章）** | 含参变量反常积分的一致收敛定义（定义 15.6）、Cauchy 准则（定理 15.4）、Weierstrass M-判别法（定理 15.5） |
| **文件 04（本章）** | 内闭一致收敛的定义（定义 15.7）、Dirichlet 判别法（定理 15.6）与 Abel 判别法（定理 15.7） |
| **文件 02（本章）** | 含参变量常义积分的可微性定理（定理 15.2）—— 定理 15.8 的证明将直接依赖该结果 |
| **文件 06（第八章）** | Gamma 函数和 Beta 函数的定义、B-Γ 关系（定理 8.11）、$\Gamma(1/2)=\sqrt{\pi}$ 的推导、三类标准型积分的计算框架 |
| **文件 01（第八章）** | 反常积分的比较判别法（定理 8.3/8.4）、$p$-积分的敛散性 |

---

## 2. 背景与应用场景

### 2.1 Gamma函数与Beta函数：含参变量反常积分的自然舞台

在前四份文件中，我们系统建立了含参变量积分的三层理论框架：

1. **常义积分层面**（文件 01-02）：连续性（定理 15.1）、可微性（定理 15.2）、积分次序交换（定理 15.3）
2. **反常积分的一致收敛理论**（文件 03）：一致收敛定义（定义 15.6）、Cauchy 准则（定理 15.4）、M-判别法（定理 15.5）
3. **条件一致收敛的判别**（文件 04）：Dirichlet 判别法（定理 15.6）、Abel 判别法（定理 15.7）、内闭一致收敛（定义 15.7）

然而，上述理论的最大价值不在于理论本身，而在于其**应用**——Gamma 函数和 Beta 函数正是含参变量反常积分理论最核心、最经典的应用案例。它们将告诉我们：一个含参变量反常积分的参数 $s$ 如何"控制"积分的收敛行为；参数连续变化时，积分值如何变化；以及如何在积分号下对参数求导。

### 2.2 为什么此时学习Euler积分？

第八章（文件 06）已经介绍了 Gamma 函数和 Beta 函数的定义、B-Γ 关系以及应用计算。但第八章的工具仅停留在"反常积分"层面——即**对每个固定的参数值**验证收敛性。

现在我们拥有了更强大的工具——含参变量反常积分的**一致收敛理论**。借助这一理论，我们可以回答第八章无法回答的问题：

- $\Gamma(s)$ 作为 $s$ 的函数是否连续？（不仅仅是"对每个 $s$ 积分收敛"）
- $\Gamma(s)$ 是否可导？$\Gamma'(s)$ 的表达式是什么？
- 为什么可以在积分号下对 $s$ 求导？

**核心洞见**：本章将 Gamma 函数从"含参数的反常积分"重新审视为"参数的函数"，研究其分析性质（连续性、可微性）——这需要一致收敛理论的支持。

### 2.3 知识图谱

```
文件 01: 含参变量常义积分 → 连续性（定理 15.1）
文件 02: 可微性（定理 15.2）+ 积分次序交换（定理 15.3）
文件 03: 一致收敛定义（定义 15.6）+ Cauchy准则（定理 15.4）+ M-判别法（定理 15.5）
文件 04: 内闭一致收敛（定义 15.7）+ Dirichlet/Abel判别法（定理 15.6/15.7）
       ↓
文件 05: ◆ Euler积分的含参理论分析 ◆
       → Gamma/Beta函数定义与基本性质（回顾定理 8.11）
       → 定理 15.8（反常含参积分可微性定理）← 核心新定理
       → Gamma函数的分析性质（收敛域、内闭一致收敛、逐阶可微性）
       → 定理 15.9（余元公式）
       → Euler积分的应用计算体系
```

---

## 3. 核心概念与符号约定

### 3.1 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $\Gamma(s)$ | Gamma 函数 | $\Gamma(s)=\int_0^{+\infty} x^{s-1}e^{-x}dx$，定义域 $s>0$ |
| $B(p,q)$ | Beta 函数 | $B(p,q)=\int_0^1 x^{p-1}(1-x)^{q-1}dx$，定义域 $p,q>0$ |
| $[a,b]$ | 参数 $s$ 的紧子区间 | 用于内闭一致收敛的验证，$0<a<b<+\infty$ |
| $\Gamma_1(s)$ | Gamma 积分的 $[0,1]$ 部分 | $\int_0^1 x^{s-1}e^{-x}dx$，瑕点在 $x=0$（当 $s<1$ 时） |
| $\Gamma_2(s)$ | Gamma 积分的 $[1,\infty)$ 部分 | $\int_1^{+\infty} x^{s-1}e^{-x}dx$，无穷限反常积分 |
| $\ln x$ | 自然对数 | Gamma 函数求导时自然出现 $\frac{\partial}{\partial s}x^{s-1}=x^{s-1}\ln x$ |
| $(2n-1)!!$ | 双阶乘 | $1\cdot3\cdot5\cdots(2n-1)$，半整数 $\Gamma$ 值的通式 |

### 3.2 定理速览

| 定理编号 | 内容 | 状态 |
|----------|------|------|
| 定理 8.11 | B-Γ 关系：$B(p,q)=\Gamma(p)\Gamma(q)/\Gamma(p+q)$ | 已证（ch8-06） |
| **定理 15.8** | **含参变量反常积分的可微性** | **本文件核心新定理** |
| **定理 15.9** | **余元公式：$\Gamma(s)\Gamma(1-s)=\pi/\sin(\pi s)$** | **本文件新定理** |

---

## 4. 原理与方法

### 4.1 Gamma函数的定义与收敛性

**定义（Gamma 函数 / Gamma Function）**：Gamma 函数定义为含参变量反常积分

$$\Gamma(s) = \int_0^{+\infty} x^{s-1}e^{-x}\,dx,\qquad s > 0$$

其中 $x^{s-1}=e^{(s-1)\ln x}$，$s$ 是参数。

**收敛性分析**：为确定收敛域，将积分在 $x=1$ 处拆分为两部分：

$$\Gamma(s) = \underbrace{\int_0^1 x^{s-1}e^{-x}\,dx}_{\Gamma_1(s)} + \underbrace{\int_1^{+\infty} x^{s-1}e^{-x}\,dx}_{\Gamma_2(s)}$$

---

**对 $\Gamma_1(s)$（$x\in(0,1]$）**：

当 $x\in(0,1]$ 时，$e^{-x}\le 1$（因为 $e^{-x}$ 在 $[0,1]$ 上递减），故

$$0 \le x^{s-1}e^{-x} \le x^{s-1}$$

积分 $\displaystyle\int_0^1 x^{s-1}\,dx$ 是瑕积分（瑕点在 $x=0$，当 $s<1$ 时）。由 $p$-积分的敛散性：$\int_0^1 x^{s-1}dx$ 收敛当且仅当 $s-1 > -1$，即 $s > 0$。由比较判别法，$\Gamma_1(s)$ 在 $s>0$ 时收敛。

**注意**：当 $s\le 0$ 时，$x^{s-1}\ge x^{-1}$，$\int_0^1 x^{-1}dx$ 发散（对数发散），故 $\Gamma_1(s)$ 发散。

---

**对 $\Gamma_2(s)$（$x\ge 1$）**：

对固定的 $s>0$，考察被积函数在无穷远处的衰减速度。由于指数函数 $e^{-x}$ 的衰减快于任何幂函数：

$$\lim_{x\to+\infty} x^{s-1}e^{-x/2} = 0$$

因此存在 $X>0$，使得当 $x>X$ 时 $x^{s-1}e^{-x/2} < 1$，从而

$$x^{s-1}e^{-x} = x^{s-1}e^{-x/2}\cdot e^{-x/2} < e^{-x/2}$$

而 $\displaystyle\int_1^{+\infty} e^{-x/2}\,dx = 2e^{-1/2} < \infty$。由比较判别法，$\Gamma_2(s)$ 对所有 $s>0$ 收敛。

（事实上，$\Gamma_2(s)$ 对所有 $s\in\mathbb{R}$ 均收敛，因为 $e^{-x}$ 的衰减最终压倒任何幂函数增长。）

---

**结论**：$\Gamma(s)$ 在 $s>0$ 时收敛，在 $s\le 0$ 时发散。因此 Gamma 函数的（初始）定义域为 $(0,+\infty)$。

### 4.2 Gamma函数的递推公式与基本性质

**递推公式**：对任意 $s>0$，有

$$\Gamma(s+1) = s\,\Gamma(s)$$

**证明**：使用分部积分法，令 $u=x^{s}$，$dv=e^{-x}dx$，则 $du=s x^{s-1}dx$，$v=-e^{-x}$：

$$\begin{aligned}
\Gamma(s+1) &= \int_0^{+\infty} x^{s}e^{-x}\,dx \\
&= \left[-x^{s}e^{-x}\right]_0^{+\infty} + \int_0^{+\infty} s x^{s-1}e^{-x}\,dx \\
&= 0 + s\int_0^{+\infty} x^{s-1}e^{-x}\,dx = s\,\Gamma(s)
\end{aligned}$$

其中 $\lim_{x\to+\infty}x^{s}e^{-x}=0$，$\lim_{x\to0^+}x^{s}e^{-x}=0$（$s>0$ 保证 $x^s\to0$）。$\square$

---

**整数取值**：由递推公式和 $\Gamma(1)=\int_0^{+\infty}e^{-x}dx=1$，对正整数 $n$ 有

$$\Gamma(n+1) = n!$$

**Gamma 函数的图像特征**（定性描述）：
- $\Gamma(s)$ 在 $s\to0^+$ 时趋于 $+\infty$（因为 $\Gamma(s) \sim 1/s$）
- $\Gamma(1)=1$，$\Gamma(2)=1$，$\Gamma(3)=2$，$\Gamma(4)=6$，$\Gamma(5)=24$
- 在 $s\approx1.46$ 处取到最小值约 $0.885$
- $\Gamma(s)\to+\infty$ 当 $s\to+\infty$

### 4.3 Beta函数的定义与基本性质

**定义（Beta 函数 / Beta Function）**：

$$B(p,q) = \int_0^1 x^{p-1}(1-x)^{q-1}\,dx,\qquad p>0,\; q>0$$

**对称性**：$B(p,q) = B(q,p)$。只需在定义中作代换 $x\mapsto 1-x$ 即可证明。

**三角形式**：令 $x = \sin^2\theta$，则 $dx = 2\sin\theta\cos\theta\,d\theta$，代入得

$$B(p,q) = 2\int_0^{\pi/2} \sin^{2p-1}\theta\,\cos^{2q-1}\theta\,d\theta$$

**有理函数形式**：在定义中作代换 $x = t/(1+t)$，得

$$B(p,q) = \int_0^{+\infty} \frac{t^{p-1}}{(1+t)^{p+q}}\,dt$$

这一形式在后续应用和余元公式的证明中非常有用。

### 4.4 B-Γ关系与基本推论

B-Γ 关系已在第八章作为定理 8.11 严格证明。这里回顾其陈述以便引用。

**定理 8.11（B-Γ 关系）**：对任意 $p,q>0$，有

$$B(p,q) = \frac{\Gamma(p)\,\Gamma(q)}{\Gamma(p+q)}$$

**证明路径回顾**：将 $\Gamma(p)\Gamma(q)$ 写成二重积分 $\iint_{x>0,y>0} x^{p-1}y^{q-1}e^{-(x+y)}\,dx\,dy$，作变量代换 $u=x+y,\; v=x/(x+y)$，Jacobi 行列式为 $|J|=u$。积分分解为 $\Gamma(p+q)\cdot B(p,q)$。详见 ch8-06 文件。

---

**推论 1（$\Gamma(1/2)=\sqrt{\pi}$）**：在 B-Γ 关系中代入 $p=q=1/2$：

$$B\!\left(\frac12,\frac12\right) = \frac{\left(\Gamma(1/2)\right)^2}{\Gamma(1)}$$

而 $B(1/2,1/2)=2\int_0^{\pi/2}1\,d\theta=\pi$，$\Gamma(1)=1$，故 $\Gamma(1/2)=\sqrt{\pi}$。

---

**推论 2（半整数 $\Gamma$ 值通式）**：利用递推公式 $\Gamma(s+1)=s\Gamma(s)$ 和 $\Gamma(1/2)=\sqrt{\pi}$：

$$\Gamma\!\left(n+\frac12\right) = \frac{(2n-1)!!}{2^n}\,\sqrt{\pi},\qquad n=0,1,2,\dots$$

其中 $(2n-1)!! = 1\cdot3\cdot5\cdots(2n-1)$，且 $(-1)!!=1$。

---

### 4.5 含参变量反常积分的可微性定理

在常义积分情形下，定理 15.2 给出了积分号下求导的条件：$f$ 和 $\partial f/\partial y$ 在闭矩形上连续即可。但对于反常积分，积分区间是无限的，直接应用定理 15.2 需要先将积分截断，然后处理"尾部"的余项——这需要一致收敛性的控制。

以下定理是本文件的核心理论结果。

---

**定理 15.8（含参变量反常积分的可微性）**：设 $f(x,y)$ 和 $\dfrac{\partial f}{\partial y}(x,y)$ 在 $[a,+\infty)\times[c,d]$ 上连续。若：

1. **逐点收敛**：含参变量反常积分 $\displaystyle F(y)=\int_a^{+\infty} f(x,y)\,dx$ 对每个 $y\in[c,d]$ 收敛；
2. **一致收敛**：含参变量反常积分 $\displaystyle G(y)=\int_a^{+\infty} \frac{\partial f}{\partial y}(x,y)\,dx$ 在 $(c,d)$ 上**内闭一致收敛**。

则 $F(y)$ 在 $(c,d)$ 上可微，且

$$F'(y) = \int_a^{+\infty} \frac{\partial f}{\partial y}(x,y)\,dx,\qquad \forall y\in(c,d)$$

---

**说明**：条件 2 只要求 $\int \partial f/\partial y$ 在内闭一致收敛（而非整体一致收敛），这是因为可微性是局部性质。若条件 2 增强为在 $[c,d]$ 上整体一致收敛，则结论可推广到 $F$ 在 $[c,d]$ 上**一致可微**。

**证明**：

记 $G(y)=\displaystyle\int_a^{+\infty} \frac{\partial f}{\partial y}(x,y)\,dx$。由条件 2 和内闭一致收敛的性质，$G$ 在 $(c,d)$ 上连续（见本节末注释）。

固定 $y\in(c,d)$。对任意 $h\neq 0$ 使得 $y+h\in(c,d)$，定义

$$\Delta(h) = \frac{F(y+h)-F(y)}{h} - G(y)$$
$$\quad = \int_a^{+\infty} \left(\frac{f(x,y+h)-f(x,y)}{h} - \frac{\partial f}{\partial y}(x,y)\right)dx$$

为证明 $\lim_{h\to0}\Delta(h)=0$，对任意 $\varepsilon>0$ 构造如下估计。

---

**第 1 步：利用一致收敛控制尾部**

由条件 2（内闭一致收敛），存在紧区间 $[y-\delta,y+\delta]\subset(c,d)$ 使得 $G$ 在该区间上一致收敛。对给定的 $\varepsilon>0$，存在 $A_0>a$，使得对任意 $A\ge A_0$ 和任意 $t\in[y-\delta,y+\delta]$ 有

$$\left|\int_A^{+\infty} \frac{\partial f}{\partial y}(x,t)\,dx\right| < \frac{\varepsilon}{3}$$

---

**第 2 步：将 $\Delta(h)$ 分解为有限部与尾部**

取 $A=A_0$，将 $\Delta(h)$ 分解为：

$$\Delta(h) = \underbrace{\int_a^{A_0}\!\left(\frac{f(x,y+h)-f(x,y)}{h} - \frac{\partial f}{\partial y}(x,y)\right)dx}_{\Delta_1(h)} + \underbrace{\int_{A_0}^{+\infty}\!\left(\frac{f(x,y+h)-f(x,y)}{h} - \frac{\partial f}{\partial y}(x,y)\right)dx}_{\Delta_2(h)}$$

---

**第 3 步：处理有限部 $\Delta_1(h)$**

在有限区间 $[a,A_0]\times[y-\delta,y+\delta]$ 上，$f$ 和 $\partial f/\partial y$ 连续。由定理 15.2（常义版本），$F_{A_0}(y)=\int_a^{A_0} f(x,y)\,dx$ 在 $(c,d)$ 上可微且导数等于 $\int_a^{A_0} \partial f/\partial y(x,y)\,dx$。因此

$$\lim_{h\to0} \Delta_1(h) = 0$$

即存在 $\eta_1>0$，使得当 $0<|h|<\eta_1$ 时 $|\Delta_1(h)|<\varepsilon/3$。

---

**第 4 步：处理尾部 $\Delta_2(h)$**

对尾部积分，应用积分中值定理的变形。由 $f$ 的连续可微性：

$$\frac{f(x,y+h)-f(x,y)}{h} = \frac{\partial f}{\partial y}(x,\,y+\theta h),\quad \theta=\theta(x,h)\in(0,1)$$

于是

$$\Delta_2(h) = \int_{A_0}^{+\infty} \left(\frac{\partial f}{\partial y}(x,y+\theta h) - \frac{\partial f}{\partial y}(x,y)\right)dx$$

利用条件 1（$F$ 收敛）和条件 2（$\int\partial f/\partial y$ 一致收敛），存在 $\eta_2>0$（与 $h$ 无关，依赖于 $A_0$ 的选取），使得当 $0<|h|<\eta_2$ 时

$$\left|\int_{A_0}^{+\infty} \frac{\partial f}{\partial y}(x,y+\theta h)\,dx\right| < \frac{\varepsilon}{3},\qquad
\left|\int_{A_0}^{+\infty} \frac{\partial f}{\partial y}(x,y)\,dx\right| < \frac{\varepsilon}{3}$$

因此 $|\Delta_2(h)| < \dfrac{2\varepsilon}{3}$。

---

**第 5 步：合并估计**

取 $\eta = \min\{\eta_1,\eta_2\}$，当 $0<|h|<\eta$ 时：

$$|\Delta(h)| \le |\Delta_1(h)| + |\Delta_2(h)| < \frac{\varepsilon}{3} + \frac{2\varepsilon}{3} = \varepsilon$$

因此 $\displaystyle\lim_{h\to0}\Delta(h)=0$，即 $F'(y)=G(y)$。$\square$

---

**定理 15.8 的证明逻辑链**：

$$
\begin{aligned}
&\text{条件 2：}\int \partial_y f \text{ 内闭一致收敛} \xrightarrow{\text{取 }A_0} \text{尾部 } \int_{A_0}^{\infty} \partial_y f \text{ 可控制} \\
&\text{条件 1 + 条件 2} \xrightarrow{\text{截断}} F(y) = \underbrace{\int_a^{A_0} f}_{F_{A_0}(y)} + \underbrace{\int_{A_0}^{\infty} f}_{\text{尾部}} \\
&\text{有限部：} F_{A_0}'(y) = \int_a^{A_0} \partial_y f \quad (\text{定理 15.2}) \\
&\text{尾部：差商 } \to \partial_y f \text{ 的积分，由一致收敛控制} \\
&\xrightarrow{\text{合并}} F'(y) = \int_a^{\infty} \partial_y f
\end{aligned}
$$

---

**注释（内闭一致收敛与连续性）**：上述证明中使用了"$G(y)$ 在 $(c,d)$ 上连续"这一事实。这是含参变量反常积分连续性定理的直接结论：若 $g(x,y)=\partial f/\partial y$ 在 $[a,+\infty)\times[c,d]$ 上连续，且 $\int_a^{+\infty} g(x,y)dx$ 内闭一致收敛，则 $G(y)$ 在 $(c,d)$ 上连续。该结论的证明是定理 15.1 的简单推广：将积分截断为有限部（由定理 15.1 得连续性）和尾部（由一致收敛得一致小），两者结合即得连续性。

---

### 4.6 Gamma函数的分析性质

现在我们将含参变量反常积分的理论工具应用于 Gamma 函数，系统研究 $\Gamma(s)$ 作为 $s$ 的函数所具有的分析性质。

---

**性质 1：内闭一致收敛性**

Gamma 函数 $\Gamma(s)=\int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 在 $(0,+\infty)$ 上**内闭一致收敛**。

**证明**：任取紧区间 $[a,b]\subset(0,+\infty)$，其中 $0<a<b<+\infty$。将积分拆分为 $\int_0^1$ 和 $\int_1^{+\infty}$ 两部分。

对 $x\in(0,1]$：$s\in[a,b]$ 时 $x^{s-1}e^{-x} \le x^{a-1}$（因为 $s-1\ge a-1$ 且 $x^{s-1}$ 关于 $s$ 递减）。$\int_0^1 x^{a-1}dx$ 收敛（$a>0$）。由 M-判别法（定理 15.5），$\int_0^1 x^{s-1}e^{-x}dx$ 在 $[a,b]$ 上一致收敛。

对 $x\ge 1$：$s\in[a,b]$ 时 $x^{s-1}e^{-x} \le x^{b-1}e^{-x}$。对充分大的 $x$，$x^{b-1}e^{-x} \le e^{-x/2}$（因为指数衰减压倒幂函数）。而 $\int_1^{+\infty} e^{-x/2}dx$ 收敛。由 M-判别法，$\int_1^{+\infty} x^{s-1}e^{-x}dx$ 在 $[a,b]$ 上一致收敛。

综上，$\Gamma(s)$ 在任意 $[a,b]\subset(0,+\infty)$ 上一致收敛，故在 $(0,+\infty)$ 上内闭一致收敛。$\square$

---

**性质 2：连续性**

由内闭一致收敛性和被积函数的连续性（定理 15.1 的推广，见前述注释），$\Gamma(s)$ 在 $(0,+\infty)$ 上连续。

---

**性质 3：逐阶可微性**

**定理**：Gamma 函数在 $(0,+\infty)$ 上无穷可微（$C^\infty$），且

$$\Gamma^{(n)}(s) = \int_0^{+\infty} x^{s-1}e^{-x}(\ln x)^n\,dx,\qquad n=1,2,3,\dots$$

**证明**（对 $n=1$ 应用定理 15.8，再对 $n>1$ 归纳）：

**第 1 步：验证一阶可导的条件**

$f(x,s)=x^{s-1}e^{-x}$，$\dfrac{\partial f}{\partial s} = x^{s-1}e^{-x}\ln x$。

对任意紧区间 $[a,b]\subset(0,+\infty)$：

- $f$ 和 $\partial f/\partial s$ 在 $[0,+\infty)\times[a,b]$ 上连续（注意 $x=0$ 处需理解为极限，但 $s>0$ 时 $x^{s-1}e^{-x}$ 在 $x=0$ 附近的瑕积分行为不影响连续性的结论——在 $[0,1]\times[a,b]$ 上，$f$ 的连续性可通过连续延拓处理）。
- $\Gamma(s)=\int_0^{+\infty} f(x,s)dx$ 对每个 $s>0$ 收敛（性质 1）。
- $\int_0^{+\infty} \partial f/\partial s(x,s)dx$ 在 $[a,b]$ 上一致收敛的验证：

  对 $x\in(0,1]$：$|x^{s-1}e^{-x}\ln x| \le x^{a-1}|\ln x|$，且 $\int_0^1 x^{a-1}|\ln x|dx$ 收敛（因为 $a>0$ 时 $x^{a-1}\ln x$ 在 $x=0$ 附近可积——可作代换 $t=-\ln x$ 验证）。

  对 $x\ge 1$：$|x^{s-1}e^{-x}\ln x| \le x^{b-1}e^{-x}\ln x \le e^{-x/2}$（对充分大 $x$），$\int_1^{+\infty} e^{-x/2}dx$ 收敛。

  由 M-判别法，$\int_0^{+\infty} \partial f/\partial s(x,s)dx$ 在 $[a,b]$ 上一致收敛。由 $[a,b]$ 的任意性，该积分在 $(0,+\infty)$ 上内闭一致收敛。

由定理 15.8，$\Gamma(s)$ 在 $(0,+\infty)$ 上可微，且

$$\Gamma'(s) = \int_0^{+\infty} x^{s-1}e^{-x}\ln x\,dx$$

---

**第 2 步：归纳证明高阶可导**

设 $\Gamma^{(n-1)}(s)=\int_0^{+\infty} x^{s-1}e^{-x}(\ln x)^{n-1}dx$ 在 $(0,+\infty)$ 上成立。被积函数对 $s$ 求偏导得：

$$\frac{\partial}{\partial s}\left[x^{s-1}e^{-x}(\ln x)^{n-1}\right] = x^{s-1}e^{-x}(\ln x)^{n}$$

对任意紧区间 $[a,b]\subset(0,+\infty)$，类似地验证：

- $x\in(0,1]$：$|x^{s-1}e^{-x}(\ln x)^n| \le x^{a-1}|\ln x|^n$，$\int_0^1 x^{a-1}|\ln x|^n dx$ 收敛。
- $x\ge 1$：$|x^{s-1}e^{-x}(\ln x)^n| \le x^{b-1}e^{-x}(\ln x)^n \le e^{-x/2}$（充分大 $x$），$\int_1^{+\infty} e^{-x/2}dx$ 收敛。

由 M-判别法，$\int_0^{+\infty} x^{s-1}e^{-x}(\ln x)^n dx$ 内闭一致收敛。由定理 15.8，$\Gamma^{(n)}(s)$ 存在且等于该积分。由归纳法，对所有 $n\in\mathbb{N}$ 成立。$\square$

---

**推论**：Gamma 函数的对数导数（Digamma 函数）定义为

$$\psi(s) = \frac{d}{ds}\ln\Gamma(s) = \frac{\Gamma'(s)}{\Gamma(s)}$$

由 $\Gamma'(s)$ 的积分表达式可得 $\psi(s)=\dfrac{\int_0^{+\infty} x^{s-1}e^{-x}\ln x\,dx}{\int_0^{+\infty} x^{s-1}e^{-x}\,dx}$。Digamma 函数在 Gamma 函数的理论中有重要应用，但在本课程中仅作概念介绍。

---

### 4.7 余元公式

余元公式（Reflection Formula）是 Euler 积分理论中继 B-Γ 关系之后的又一个核心公式，揭示了 Gamma 函数在 $s$ 和 $1-s$ 处的对称关系。

---

**定理 15.9（余元公式 / Reflection Formula）**：对 $0 < s < 1$，有

$$\Gamma(s)\,\Gamma(1-s) = \frac{\pi}{\sin(\pi s)}$$

特别地，$\Gamma\left(\dfrac12\right) = \sqrt{\pi}$（与 B-Γ 关系推论一致）。

---

**证明**：

**第 1 步：转化为 Beta 函数**

由 B-Γ 关系（定理 8.11），$B(s,1-s) = \dfrac{\Gamma(s)\Gamma(1-s)}{\Gamma(1)} = \Gamma(s)\Gamma(1-s)$。因此只需证明

$$B(s,1-s) = \frac{\pi}{\sin(\pi s)}$$

---

**第 2 步：将 Beta 函数化为无穷限积分**

在 $B(s,1-s)=\int_0^1 x^{s-1}(1-x)^{-s}dx$ 中作代换 $x = \dfrac{t}{1+t}$，则 $dx = \dfrac{dt}{(1+t)^2}$，$1-x=\dfrac{1}{1+t}$：

$$B(s,1-s) = \int_0^{+\infty} \frac{t^{s-1}}{1+t}\,dt$$

---

**第 3 步：将积分拆分为 $[0,1]$ 和 $[1,\infty)$ 两部分**

$$\begin{aligned}
B(s,1-s) &= \int_0^1 \frac{t^{s-1}}{1+t}\,dt + \int_1^{+\infty} \frac{t^{s-1}}{1+t}\,dt \\
&= \int_0^1 \frac{t^{s-1}}{1+t}\,dt + \int_0^1 \frac{u^{-s}}{1+u}\,du \qquad (u=1/t) \\
&= \int_0^1 \frac{t^{s-1}+t^{-s}}{1+t}\,dt
\end{aligned}$$

---

**第 4 步：利用几何级数展开**

在 $[0,1)$ 上，$\dfrac{1}{1+t} = \displaystyle\sum_{n=0}^{\infty} (-1)^n t^n$，该级数在 $[0,1-\delta]$ 上一致收敛。由一致收敛级数可逐项积分：

$$\begin{aligned}
B(s,1-s) &= \int_0^1 (t^{s-1}+t^{-s})\sum_{n=0}^{\infty}(-1)^n t^n\,dt \\
&= \sum_{n=0}^{\infty}(-1)^n \int_0^1 (t^{s-1+n} + t^{-s+n})\,dt \\
&= \sum_{n=0}^{\infty}(-1)^n \left(\frac{1}{n+s} + \frac{1}{n+1-s}\right)
\end{aligned}$$

---

**第 5 步：识别为 $\pi\csc(\pi s)$ 的级数展开**

由数学分析中已知的恒等式（可通过 Fourier 级数展开或余切函数的 Mittag-Leffler 展开得到）：

$$\pi\csc(\pi s) = \frac{\pi}{\sin(\pi s)} = \sum_{n=0}^{\infty}(-1)^n\left(\frac{1}{n+s} + \frac{1}{n+1-s}\right),\quad 0<s<1$$

该级数恰好等于第 4 步中得到的 $B(s,1-s)$ 的级数表达式。因此

$$B(s,1-s) = \frac{\pi}{\sin(\pi s)}$$

代入 $B(s,1-s)=\Gamma(s)\Gamma(1-s)$ 即得余元公式。$\square$

---

**注释**：上述证明中步骤 5 引用了 $\pi\csc(\pi s)$ 的级数展开。该展开的严格推导属于 Fourier 级数理论或复分析中留数定理的应用范畴，可在后续章节（第十六章 Fourier 级数）中找到完整证明。本文件将其作为已知结论使用，读者可暂时接受该恒等式。

---

**余元公式的直接推论**：

$$\Gamma\left(\frac12\right) = \sqrt{\pi}$$

在余元公式中取 $s=1/2$：$\Gamma(1/2)^2 = \pi/\sin(\pi/2) = \pi$，故 $\Gamma(1/2)=\sqrt{\pi}$（取正根）。这与 B-Γ 关系推导的结果完全一致。

**余元公式的拓展**：对于 $s>0$ 且 $s\notin\mathbb{Z}$，余元公式可通过 Gamma 函数的解析延拓推广为

$$\Gamma(s)\Gamma(1-s) = \frac{\pi}{\sin(\pi s)}$$

（本课程的讨论限于 $s>0$，解析延拓的理论属于复变函数范畴。）

---

### 4.8 Euler积分的应用计算

Gamma 函数和 Beta 函数为三类标准型积分提供了统一的计算框架。以下系统总结。

---

**类型 I：有理函数型**

$$\int_0^{+\infty}\frac{x^{p-1}}{(1+x)^{p+q}}\,dx = B(p,q) = \frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)},\qquad p,q>0$$

**匹配方法**：设被积函数为 $\dfrac{x^m}{(1+x)^n}$，则 $p=m+1$，$q=n-m-1$。

---

**类型 II：三角型**

$$\int_0^{\pi/2}\sin^{2p-1}\theta\cos^{2q-1}\theta\,d\theta = \frac12 B(p,q) = \frac{\Gamma(p)\Gamma(q)}{2\,\Gamma(p+q)},\qquad p,q>0$$

**匹配方法**：设被积函数为 $\sin^a\theta\cos^b\theta$，则 $p=(a+1)/2$，$q=(b+1)/2$。

---

**类型 III：余元公式型**

$$\int_0^{+\infty}\frac{x^{\alpha-1}}{1+x}\,dx = \frac{\pi}{\sin(\pi\alpha)},\qquad 0<\alpha<1$$

该类积分通过代换 $x = t/(1-t)$ 化为 Beta 函数 $B(\alpha,1-\alpha)$，再应用余元公式。

---

**三类积分对照表**：

| 类型 | 标准形式 | 转化为 | 结果 |
|------|---------|--------|------|
| 有理函数型 | $\displaystyle\int_0^{\infty}\frac{x^{p-1}}{(1+x)^{p+q}}dx$ | $B(p,q)$ | $\dfrac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}$ |
| 三角型 | $\displaystyle\int_0^{\pi/2}\sin^{2p-1}\theta\cos^{2q-1}\theta d\theta$ | $\frac12 B(p,q)$ | $\dfrac{\Gamma(p)\Gamma(q)}{2\Gamma(p+q)}$ |
| 多项式型 | $\displaystyle\int_0^1 x^{p-1}(1-x)^{q-1}dx$ | $B(p,q)$ | $\dfrac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}$ |
| 余元型 | $\displaystyle\int_0^{\infty}\frac{x^{\alpha-1}}{1+x}dx$ | $B(\alpha,1-\alpha)$ | $\dfrac{\pi}{\sin(\pi\alpha)}$ |

---

## 5. 例题精讲

### 例 5-1：验证Gamma函数的内闭一致收敛

用 M-判别法验证 $\Gamma(s)=\int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 在 $[1,2]$ 上一致收敛。

**解**：

**第 1 步：拆分积分**。

$$\Gamma(s) = \int_0^1 x^{s-1}e^{-x}\,dx + \int_1^{+\infty} x^{s-1}e^{-x}\,dx$$

---

**第 2 步：处理 $[0,1]$ 部分**。

$s\in[1,2]$ 时，$s-1\ge 0$，故 $x^{s-1}\le 1$，$e^{-x}\le 1$，因此

$$|x^{s-1}e^{-x}| \le 1\quad\text{在 }[0,1]\text{ 上}$$

$\int_0^1 1\,dx = 1<\infty$。由 M-判别法一致收敛。✓

---

**第 3 步：处理 $[1,\infty)$ 部分**。

$s\in[1,2]$ 时，$s-1\le 1$，故 $x^{s-1}\le x$，因此

$$|x^{s-1}e^{-x}| \le x e^{-x}$$

对充分大的 $x$，$x e^{-x} \le e^{-x/2}$。$\int_1^{+\infty} e^{-x/2}dx = 2e^{-1/2}<\infty$。由 M-判别法一致收敛。✓

---

**结论**：$\Gamma(s)$ 在 $[1,2]$ 上一致收敛。由于 $[1,2]$ 是 $(0,+\infty)$ 的任意紧子集的代表，故 $\Gamma(s)$ 在 $(0,+\infty)$ 上内闭一致收敛。

---

### 例 5-2：利用定理15.8证明Gamma函数的可微性

验证 $\Gamma(s)=\int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 满足定理 15.8 的条件，并写出 $\Gamma'(s)$ 的表达式。

**解**：

**第 1 步：验证被积函数的连续可微性**。

$f(x,s)=x^{s-1}e^{-x}$ 在 $(0,+\infty)\times(0,+\infty)$ 上关于 $x$ 和 $s$ 均连续。其偏导数

$$\frac{\partial f}{\partial s} = x^{s-1}e^{-x}\ln x$$

在 $(0,+\infty)\times(0,+\infty)$ 上连续（$x=0$ 处需注意瑕积分行为，但在任何 $[a,b]\subset(0,+\infty)$ 上，$\partial f/\partial s$ 连续）。✓

---

**第 2 步：验证 $\int \partial f/\partial s$ 的内闭一致收敛**。

任取 $[a,b]\subset(0,+\infty)$。对 $x\in(0,1]$：

$$|x^{s-1}e^{-x}\ln x| \le x^{a-1}|\ln x|$$

积分 $\int_0^1 x^{a-1}|\ln x|dx$ 收敛（因为 $a>0$ 时被积函数在 $x=0$ 附近可积）。✓

对 $x\ge 1$：

$$|x^{s-1}e^{-x}\ln x| \le x^{b-1}e^{-x}\ln x \le e^{-x/2}\quad(\text{充分大 }x)$$

$\int_1^{+\infty} e^{-x/2}dx$ 收敛。✓

由 M-判别法，$\int_0^{+\infty} \partial f/\partial s(x,s)\,dx$ 在 $[a,b]$ 上一致收敛。由 $[a,b]$ 的任意性，内闭一致收敛成立。✓

---

**第 3 步：应用定理 15.8**。

所有条件满足，故

$$\Gamma'(s) = \int_0^{+\infty} \frac{\partial}{\partial s}\big(x^{s-1}e^{-x}\big)\,dx = \int_0^{+\infty} x^{s-1}e^{-x}\ln x\,dx$$

---

### 例 5-3：利用Euler积分计算定积分（一）

计算 $I = \displaystyle\int_0^{+\infty} \frac{x^2}{(1+x)^5}\,dx$。

**解**：

**第 1 步：匹配类型 I 公式**。

被积函数为 $\dfrac{x^2}{(1+x)^5}$。对比 $\dfrac{x^{p-1}}{(1+x)^{p+q}}$：

分子 $x^2$ 对应 $p-1=2\Rightarrow p=3$。
分母 $(1+x)^5$ 对应 $p+q=5\Rightarrow q=2$。

---

**第 2 步：代入 B-Γ 关系**。

$$I = B(3,2) = \frac{\Gamma(3)\Gamma(2)}{\Gamma(5)}$$

---

**第 3 步：计算 Gamma 值**。

$\Gamma(3)=2! = 2$，$\Gamma(2)=1! = 1$，$\Gamma(5)=4! = 24$。

$$I = \frac{2\times 1}{24} = \frac{1}{12}$$

---

### 例 5-4：利用Euler积分计算定积分（二）

计算 $J = \displaystyle\int_0^{+\infty} \frac{x^{1/3}}{(1+x)^2}\,dx$。

**解**：

**第 1 步：匹配类型 I 公式**。

$\dfrac{x^{1/3}}{(1+x)^2}$：$p-1 = 1/3 \Rightarrow p = 4/3$；$p+q = 2 \Rightarrow q = 2-4/3 = 2/3$。

---

**第 2 步：代入 B-Γ 关系**。

$$J = B\!\left(\frac43,\frac23\right) = \frac{\Gamma(4/3)\Gamma(2/3)}{\Gamma(2)} = \Gamma\!\left(\frac43\right)\Gamma\!\left(\frac23\right)$$

---

**第 3 步：利用余元公式**。

注意 $\Gamma(4/3) = \Gamma(1+1/3) = \frac13\,\Gamma(1/3)$（递推公式）。于是

$$J = \frac13\,\Gamma\!\left(\frac13\right)\Gamma\!\left(\frac23\right)$$

由余元公式，$\Gamma(1/3)\Gamma(2/3) = \dfrac{\pi}{\sin(\pi/3)} = \dfrac{\pi}{\sqrt{3}/2} = \dfrac{2\pi}{\sqrt{3}}$。代入：

$$J = \frac13 \cdot \frac{2\pi}{\sqrt{3}} = \frac{2\pi}{3\sqrt{3}}$$

---

### 例 5-5：利用余元公式计算

计算 $K = \displaystyle\int_0^{+\infty} \frac{x^{\alpha-1}}{1+x}\,dx$，$0<\alpha<1$。

**解**：

**第 1 步：化为 Beta 函数**。

作代换 $x = \dfrac{t}{1-t}$，$dx = \dfrac{dt}{(1-t)^2}$，积分限 $t\in[0,1)$：

$$\int_0^{+\infty} \frac{x^{\alpha-1}}{1+x}\,dx = \int_0^1 t^{\alpha-1}(1-t)^{-\alpha}\,dt = B(\alpha, 1-\alpha)$$

---

**第 2 步：应用余元公式**。

$$K = B(\alpha,1-\alpha) = \frac{\Gamma(\alpha)\Gamma(1-\alpha)}{\Gamma(1)} = \frac{\pi}{\sin(\pi\alpha)}$$

---

### 例 5-6：综合应用——Dirichlet积分

计算 $L = \displaystyle\int_0^{\pi/2} \sin^6\theta\,d\theta$。

**解**：

**第 1 步：匹配类型 II 公式**。

$\sin^6\theta$：$2p-1 = 6 \Rightarrow p = 7/2$；$\cos$ 的指数为 $0$：$2q-1 = 0 \Rightarrow q = 1/2$。

---

**第 2 步：代入三角公式**。

$$L = \frac12\,B\!\left(\frac72,\frac12\right) = \frac12\cdot\frac{\Gamma(7/2)\Gamma(1/2)}{\Gamma(7/2+1/2)} = \frac12\cdot\frac{\Gamma(7/2)\sqrt{\pi}}{\Gamma(4)}$$

---

**第 3 步：计算 Gamma 值**。

$\Gamma(7/2) = \dfrac{5}{2}\cdot\dfrac{3}{2}\cdot\dfrac{1}{2}\sqrt{\pi} = \dfrac{15}{8}\sqrt{\pi}$，$\Gamma(4) = 6$。

$$L = \frac12\cdot\frac{\frac{15}{8}\sqrt{\pi}\cdot\sqrt{\pi}}{6} = \frac12\cdot\frac{15\pi/8}{6} = \frac12\cdot\frac{15\pi}{48} = \frac{15\pi}{96} = \frac{5\pi}{32}$$

---

### 例 5-7：确定Gamma函数的收敛域

证明 $\Gamma(s) = \int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 在 $s>0$ 时收敛，在 $s\le 0$ 时发散。

**解**：

**收敛性**（$s>0$）：已在 4.1 节完成证明。

**发散性**（$s\le 0$）：

当 $s\le 0$ 时，考虑 $[0,1]$ 部分的积分。在 $x\in(0,1]$ 上，$e^{-x}\ge e^{-1}$（因为 $e^{-x}$ 递减），故

$$x^{s-1}e^{-x} \ge e^{-1}x^{s-1}$$

当 $s\le 0$ 时，$s-1\le -1$，$\int_0^1 x^{s-1}dx$ 发散（$p$-积分 $p\le 1$ 发散）。由比较判别法，$\Gamma_1(s)$ 发散，故 $\Gamma(s)$ 发散。$\square$

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| Gamma 函数的收敛性证明只需指出 $x^{s-1}e^{-x}$ 在 $\infty$ 处衰减很快即可 | 需要分别处理 $x\to0^+$ 和 $x\to\infty$ 两个端点的收敛性。$x\to0^+$ 的收敛要求 $s>0$，$x\to\infty$ 的收敛对任意 $s$ 均成立。所以收敛域由 $x=0$ 附近的收敛性决定。 |
| Gamma 函数的可微性可以直接用定理 15.2（常义版本） | 定理 15.2 要求积分区间有限。Gamma 函数的积分是反常积分，必须使用定理 15.8（含参变量反常积分的可微性定理），其中关键条件是 $\int \partial f/\partial s$ 的内闭一致收敛。 |
| 余元公式的证明中，拆分 $[0,1]$ 和 $[1,\infty)$ 两部分后直接相加即可 | 两部分分别处理后需用代换 $u=1/t$ 将第二部分统一为 $[0,1]$ 上的积分，然后才能使用几何级数展开。直接处理会得到两个不同区间的级数，无法合并。 |
| $\int_0^{+\infty} x^{s-1}e^{-x}\,dx$ 的内闭一致收敛只需对 $s$ 在紧区间上逐点验证 | 必须用 M-判别法找到与 $s$ 无关的优函数 $M(x)$。在 $[a,b]\subset(0,+\infty)$ 上，$[0,1]$ 部分的优函数是 $x^{a-1}$，$[1,\infty)$ 部分的优函数是 $x^{b-1}e^{-x}$（或 $e^{-x/2}$）。两个优函数都与 $s$ 无关。 |
| Gamma 函数的高阶可导性需要逐阶重复验证 M-判别法 | 每增加一阶导数，被积函数中多一个 $\ln x$ 因子。由于 $|\ln x|^n$ 在 $x\to0^+$ 和 $x\to\infty$ 处的增长慢于任何幂函数，M-判别法对任意 $n$ 均适用。可由此归纳证明 $C^\infty$ 性质。 |
| $\Gamma'(s)$ 的积分表达式中 $\ln x$ 在 $x\in(0,1)$ 上为负，可能导致积分不收敛 | $\ln x < 0$ 仅影响符号，不影响绝对值的可积性。M-判别法中使用的是绝对值估计 $x^{a-1}|\ln x|$，正负无关紧要。 |

### 检查点

- [ ] 能否独立写出 Gamma 函数收敛性分析的全过程（拆分 $\int_0^1$ 和 $\int_1^{+\infty}$，分别用比较判别法）？
- [ ] 能否准确写出定理 15.8（含参变量反常积分的可微性）的条件与结论，并指出与定理 15.2（常义版本）的差异？
- [ ] 能否给出定理 15.8 的完整证明（截断、有限部用定理 15.2、尾部用一致收敛控制）？
- [ ] 能否用 M-判别法严格验证 $\Gamma(s)$ 在 $(0,+\infty)$ 上的内闭一致收敛性？
- [ ] 能否应用定理 15.8 严格证明 $\Gamma'(s)=\int_0^\infty x^{s-1}e^{-x}\ln x\,dx$？
- [ ] 能否写出余元公式（定理 15.9）的完整陈述，并给出基于 Beta 函数和级数展开的证明思路？
- [ ] 能否在 30 秒内说出三类标准型积分（有理函数型、三角型、余元型）各自对应的 Gamma/Beta 公式？
- [ ] 能否独立完成例 5-3 到例 5-6 的计算？
- [ ] 能否说出 Gamma 函数递推公式的证明使用了分部积分法，并指出边界项为零的理由？
- [ ] 应用定理 15.8 时，最关键的条件是哪一个？为什么需要"内闭一致收敛"而非"逐点收敛"？

---

## 练习题

### 基础巩固

**1.**（Gamma 收敛域）判断下列 $s$ 值是否在 $\Gamma(s)$ 的收敛域内，并说明理由：

(1) $s = 0.5$
(2) $s = 2$
(3) $s = -0.5$
(4) $s = 0$

<details><summary>参考答案</summary>

(1) 收敛域 $s>0$，$0.5>0$，故 $\Gamma(0.5)$ 收敛。✓
(2) 收敛域 $s>0$，$2>0$，故 $\Gamma(2)=1! = 1$，收敛。✓
(3) $s=-0.5\le 0$，$\int_0^1 x^{-1.5}dx$ 发散，故 $\Gamma(-0.5)$ 发散。✗
(4) $s=0$，$\int_0^1 x^{-1}dx$ 发散（对数发散），故 $\Gamma(0)$ 发散。✗

</details>

---

**2.**（递推公式）利用递推公式 $\Gamma(s+1)=s\Gamma(s)$ 和 $\Gamma(1/2)=\sqrt{\pi}$，计算：

(1) $\Gamma(7/2)$
(2) $\Gamma(9/2)$

<details><summary>参考答案</summary>

(1) $\Gamma(7/2) = \frac52\cdot\Gamma(5/2) = \frac52\cdot\frac34\sqrt{\pi} = \frac{15}{8}\sqrt{\pi}$
(2) $\Gamma(9/2) = \frac72\cdot\Gamma(7/2) = \frac72\cdot\frac{15}{8}\sqrt{\pi} = \frac{105}{16}\sqrt{\pi}$

</details>

---

**3.**（B-Γ 关系）利用 B-Γ 关系计算下列积分：

(1) $\displaystyle\int_0^1 x^3(1-x)^4\,dx$
(2) $\displaystyle\int_0^{\pi/2} \sin^5\theta\cos^3\theta\,d\theta$

<details><summary>参考答案</summary>

(1) $\int_0^1 x^3(1-x)^4dx = B(4,5) = \dfrac{\Gamma(4)\Gamma(5)}{\Gamma(9)} = \dfrac{3!\cdot4!}{8!} = \dfrac{6\cdot24}{40320} = \dfrac{144}{40320} = \dfrac{1}{280}$

(2) $\sin^5\theta\cos^3\theta$：$2p-1=5\Rightarrow p=3$，$2q-1=3\Rightarrow q=2$。
$$\int_0^{\pi/2} \sin^5\theta\cos^3\theta\,d\theta = \frac12 B(3,2) = \frac12\cdot\frac{\Gamma(3)\Gamma(2)}{\Gamma(5)} = \frac12\cdot\frac{2\cdot1}{24} = \frac{1}{24}$$

</details>

---

### 迁移应用

**4.**（余元公式）利用余元公式计算：

(1) $\displaystyle\int_0^{+\infty} \frac{x^{-1/3}}{1+x}\,dx$
(2) $\displaystyle\int_0^{\pi/2} \tan^{\alpha}x\,dx$，$|\alpha|<1$

<details><summary>参考答案</summary>

(1) 直接应用余元公式，$\alpha = 2/3$：
$$\int_0^{+\infty} \frac{x^{-1/3}}{1+x}\,dx = \frac{\pi}{\sin(2\pi/3)} = \frac{\pi}{\sqrt{3}/2} = \frac{2\pi}{\sqrt{3}}$$

(2) $\int_0^{\pi/2}\tan^{\alpha}x\,dx = \int_0^{\pi/2} \sin^{\alpha}x\cos^{-\alpha}x\,dx$。
匹配三角公式：$2p-1=\alpha \Rightarrow p=(\alpha+1)/2$，$2q-1=-\alpha \Rightarrow q=(1-\alpha)/2$。
$$\int_0^{\pi/2}\tan^{\alpha}x\,dx = \frac12 B\!\left(\frac{\alpha+1}{2},\frac{1-\alpha}{2}\right) = \frac12\cdot\frac{\pi}{\sin\frac{(\alpha+1)\pi}{2}} = \frac{\pi}{2\cos\frac{\alpha\pi}{2}}$$

</details>

---

**5.**（Gamma 函数的可微性）设 $\Gamma(s)=\int_0^{+\infty} x^{s-1}e^{-x}\,dx$，$s>0$。

(1) 写出 $\Gamma'(s)$ 和 $\Gamma''(s)$ 的积分表达式。
(2) 证明 $\Gamma''(s) > 0$ 对一切 $s>0$ 成立（即 $\Gamma$ 是严格凸函数）。

**提示**：注意被积函数的符号。

<details><summary>参考答案</summary>

(1) $\Gamma'(s) = \int_0^{+\infty} x^{s-1}e^{-x}\ln x\,dx$
$\Gamma''(s) = \int_0^{+\infty} x^{s-1}e^{-x}(\ln x)^2\,dx$

(2) 对任意 $s>0$，被积函数 $x^{s-1}e^{-x}(\ln x)^2 \ge 0$ 且不恒为零（$\ln x$ 在 $x\neq1$ 时非零）。因此 $\Gamma''(s) = \int_0^{+\infty} x^{s-1}e^{-x}(\ln x)^2\,dx > 0$。由二阶导数恒正知 $\Gamma(s)$ 是严格凸函数。

</details>

---

**6.*（选做——余元公式的另一种推导）**

利用 $\Gamma$ 函数的 Weierstrass 乘积表示：
$$\frac{1}{\Gamma(s)} = se^{\gamma s}\prod_{n=1}^{\infty}\left(1+\frac{s}{n}\right)e^{-s/n}$$

其中 $\gamma$ 是 Euler 常数。推导余元公式 $\Gamma(s)\Gamma(1-s)=\pi/\sin(\pi s)$。

**提示**：由乘积表示可得 $\Gamma(s)\Gamma(-s) = -\pi/(s\sin(\pi s))$，再利用 $\Gamma(1-s)=-s\Gamma(-s)$ 转化。

<details><summary>参考答案（思路）</summary>

该推导需要复分析或无穷乘积理论，作为选做内容供有兴趣的读者探索。核心思路是将 $\Gamma(s)\Gamma(1-s)$ 的乘积展开与 $\sin(\pi s)$ 的 Weierstrass 乘积 $\sin(\pi s)=\pi s\prod_{n=1}^{\infty}(1-s^2/n^2)$ 联系起来。

</details>
