# 05. 函数的幂级数展开——Taylor级数与Maclaurin级数

> 所属章节：第十章 函数项级数  |  文件序号：05  |  难度：进阶
> 常见混淆点：Taylor多项式（有限和+余项）与Taylor级数（无穷级数）是两个不同的概念——前者是多项式逼近，后者是幂级数展开，两者通过 $R_n(x)\to 0$ 的极限过程连接；"Taylor级数收敛"与"Taylor级数收敛于 $f(x)$"是两个不同的问题——前者是幂级数作为级数的收敛性，后者是级数的和是否恰好等于 $f(x)$

## 1. 学习目标与先修前置

### 学习目标
- 理解Taylor多项式（Ch5-03）与Taylor级数（本文件）的本质区别与联系
- 掌握Taylor级数和Maclaurin级数的形式化定义：$f(x) = \displaystyle\sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n$
- 掌握Taylor级数收敛于原函数的充要条件：$\displaystyle\lim_{n\to\infty} R_n(x) = 0$
- 能利用Lagrange余项证明 $e^x$、$\sin x$、$\cos x$ 的Maclaurin级数在 $\mathbb{R}$ 上收敛于原函数
- 掌握 $\displaystyle\lim_{n\to\infty} \frac{x^n}{n!} = 0$ 的严格证明
- 熟记五个基本初等函数的Maclaurin级数，能写出通项公式和收敛半径
- 掌握广义二项式定理：$(1+x)^\alpha = \displaystyle\sum_{n=0}^{\infty} \binom{\alpha}{n} x^n$ 及收敛半径 $R=1$
- 掌握 $\sqrt[n]{n!} \to \infty$ 的证明

### 先修知识
- 文件 03（第十章）：幂级数的定义（定义10.6）、Cauchy-Hadamard公式（定理10.4）——Taylor级数是一种特殊的幂级数
- 文件 04（第十章）：逐项积分定理（定理10.6）、逐项求导定理（定理10.7）——从已知级数通过积分/求导推出新级数的核心工具
- 文件 03（第五章）：Taylor公式——有限Taylor多项式的定义（定理5.5、定理5.6）、Lagrange余项公式、$e^x$、$\sin x$、$\cos x$、$\ln(1+x)$ 的有限Taylor展开

---

## 2. 背景与应用场景

在第五章中，我们学习了**Taylor多项式**：对在 $x=a$ 处 $n$ 阶可导的函数 $f(x)$，有

$$f(x) = P_n(x) + R_n(x)$$

其中

$$P_n(x) = \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k$$

是 $n$ 次Taylor多项式，$R_n(x)$ 是余项（Peano型或Lagrange型）。

Taylor多项式给出了函数在 $x=a$ 附近的多项式逼近。一个自然的问题是：**如果让 $n$ 趋于无穷，会得到什么？**

从形式上看，令 $n \to \infty$，Taylor多项式 $P_n(x)$ 变成一个无穷级数：

$$\sum_{k=0}^{\infty} \frac{f^{(k)}(a)}{k!} (x-a)^k = \frac{f(a)}{0!} + \frac{f'(a)}{1!} (x-a) + \frac{f''(a)}{2!} (x-a)^2 + \cdots$$

这是一个**幂级数**——其通项是系数 $a_k = f^{(k)}(a)/k!$ 乘以 $(x-a)^k$。这个幂级数称为 $f(x)$ 的**Taylor级数**。

但这里有两个根本问题需要解决：

1. **收敛性**：这个幂级数是否收敛？如果收敛，其收敛域是什么？
2. **展开的合理性**：即使幂级数收敛，其和是否恰好等于 $f(x)$？还是收敛到别的函数？

回答这两个问题，正是本文件的核心任务。

从应用角度看，函数的幂级数展开是分析学中最强大的工具之一：
- **数值计算**：用幂级数的部分和近似计算函数值（如 $e^x$、$\sin x$、$\ln(1+x)$）
- **微分方程**：将未知函数展开为幂级数求解（幂级数解法）
- **积分计算**：对被积函数展开后逐项积分
- **极限计算**：用幂级数展开处理 $\frac{0}{0}$ 型极限

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 示例 |
|------|------|------|
| $T_f(x;a)$ | 函数 $f$ 在 $x=a$ 处的Taylor级数 | $T_{\sin}(x;0) = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{(2n+1)!}$ |
| $P_n(x)$ | $n$ 次Taylor多项式（有限和） | $P_3(x) = x - \frac{x^3}{6}$（$\sin x$） |
| $R_n(x)$ | Taylor余项 | $R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} x^{n+1}$（Lagrange型） |
| $\binom{\alpha}{n}$ | 广义二项式系数（$\alpha \in \mathbb{R}$） | $\binom{1/2}{2} = \frac{(1/2)(-1/2)}{2!} = -\frac{1}{8}$ |
| $\mathcal{C}^\infty(a,b)$ | 在 $(a,b)$ 上无穷次可微的函数类 | $e^x \in \mathcal{C}^\infty(\mathbb{R})$ |

### 3.1 从Taylor多项式到Taylor级数——形式上的极限

回顾Taylor多项式（定理5.6，Lagrange余项形式）：

$$f(x) = \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k + R_n(x)$$

其中余项

$$R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} (x-a)^{n+1}, \quad \xi \text{ 介于 } x \text{ 与 } a \text{ 之间}$$

若我们固定 $x$，考察当 $n$ 增大时的情况：
- 多项式部分 $P_n(x)$ 的项数不断增加，越来越"长"
- 余项 $R_n(x)$ 是 $n+1$ 次项，随 $n$ 增大而改变

如果 $\lim_{n\to\infty} R_n(x) = 0$，则取极限得到：

$$f(x) = \lim_{n\to\infty} \left[ \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k + R_n(x) \right] = \sum_{k=0}^{\infty} \frac{f^{(k)}(a)}{k!} (x-a)^k$$

这个无穷级数就是 $f(x)$ 在 $x=a$ 处的**Taylor级数**。

### 3.2 Taylor级数与Maclaurin级数的形式化定义

**定义 10.9（Taylor级数与Maclaurin级数）**：设 $f(x)$ 在 $x=a$ 处无穷次可导（即对一切 $n \in \mathbb{N}$，$f^{(n)}(a)$ 存在）。

1. 称幂级数

   $$T_f(x;a) = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n$$

   为 $f(x)$ 在 $x=a$ 处的 **Taylor级数**（Taylor series）。系数 $a_n = f^{(n)}(a)/n!$ 称为 **Taylor系数**。

2. 当 $a = 0$ 时，称

   $$T_f(x;0) = \sum_{n=0}^{\infty} \frac{f^{(n)}(0)}{n!} x^n$$

   为 $f(x)$ 的 **Maclaurin级数**（Maclaurin series）。

**注**：Taylor级数的定义只需要 $f(x)$ 在 $x=a$ 处无穷次可导。但定义本身只给出了一个幂级数，并没有说这个级数是否收敛，更没有说它是否收敛于 $f(x)$。这两个问题需要单独判断。

### 3.3 Taylor级数收敛于原函数的充要条件

Taylor级数 $T_f(x;a)$ 是一个幂级数，它自然有一个收敛半径 $R$（由Cauchy-Hadamard公式或比值法确定）。在收敛区间 $(-R, R)$ 内，它收敛于某个和函数 $S(x)$。

但关键在于：**$S(x)$ 是否等于 $f(x)$？**

**定理 10.10（Taylor级数收敛于原函数的判定）**：设 $f(x)$ 在 $x=a$ 的某邻域内无穷次可导，$R_n(x)$ 为 Taylor 公式的 Lagrange 余项。则 Taylor 级数 $\sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n$ 在区间 $I$ 上收敛于 $f(x)$ 当且仅当对每个 $x \in I$，

$$\lim_{n\to\infty} R_n(x) = 0$$

**证明**：

对每个固定的 $x \in I$，将 Taylor 公式写作：

$$f(x) - \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k = R_n(x)$$

其中 $R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} (x-a)^{n+1}$，$\xi$ 介于 $x$ 与 $a$ 之间。

Taylor 级数的部分和 $S_n(x) = \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k$ 正是 $n$ 次 Taylor 多项式。因此：

$$f(x) - S_n(x) = R_n(x)$$

级数 $\sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n$ 收敛于 $f(x)$ 的定义是 $\lim_{n\to\infty} S_n(x) = f(x)$，而

$$\lim_{n\to\infty} S_n(x) = f(x) \iff \lim_{n\to\infty} [f(x) - S_n(x)] = 0 \iff \lim_{n\to\infty} R_n(x) = 0$$

证毕。

**注**：定理 10.10 给出了一个**充要条件**，但实际使用中，直接证明 $\lim R_n(x) = 0$ 通常依赖于对 $f^{(n+1)}(\xi)$ 的有界性估计。

### 3.4 关键引理：$\displaystyle\lim_{n\to\infty} \frac{x^n}{n!} = 0$

在后续用 Lagrange 余项证明 $e^x$、$\sin x$、$\cos x$ 的 Maclaurin 级数在 $\mathbb{R}$ 上收敛时，我们需要反复用到这一事实。

**引理 10.1**：对任意固定的实数 $x$，有

$$\lim_{n\to\infty} \frac{x^n}{n!} = 0$$

**证明**：

**第 1 步：将问题转化为控制项**。

若 $x = 0$，结论显然成立（数列恒为 $0$）。下设 $x \neq 0$。取正整数 $N$ 满足 $N > 2|x|$。例如，$N = \lfloor 2|x| \rfloor + 1$。

**第 2 步：对 $n > N$ 进行递推放缩**。

对任意 $n > N$，将 $\frac{|x|^n}{n!}$ 分解为前 $N$ 项和后 $n-N$ 项的乘积：

$$\frac{|x|^n}{n!} = \underbrace{\frac{|x|^N}{N!}}_{\text{常数 }C} \cdot \underbrace{\frac{|x|}{N+1} \cdot \frac{|x|}{N+2} \cdots \frac{|x|}{n}}_{\text{每个因子都小于 }1/2}$$

由于对任意 $k > N$ 有 $N > 2|x|$，故 $k > 2|x|$，从而 $\frac{|x|}{k} < \frac{1}{2}$。因此：

$$\frac{|x|^n}{n!} \le C \cdot \left(\frac{1}{2}\right)^{n-N}, \quad \forall n > N$$

其中 $C = \frac{|x|^N}{N!}$ 是与 $n$ 无关的常数。

**第 3 步：用夹逼定理**。

对 $n > N$ 有：

$$0 \le \left|\frac{x^n}{n!}\right| = \frac{|x|^n}{n!} \le C \cdot \left(\frac{1}{2}\right)^{n-N}$$

当 $n \to \infty$ 时，$\left(\frac{1}{2}\right)^{n-N} \to 0$（公比 $< 1$ 的等比数列趋于 $0$）。由夹逼定理：

$$\lim_{n\to\infty} \frac{x^n}{n!} = 0$$

证毕。

**注**：这个证明的精髓在于"取 $N > 2|x|$"——一旦 $n$ 超过 $2|x|$，后续每个因子 $\frac{|x|}{k}$ 都小于 $1/2$，从而使尾部被等比数列控制。

---

## 4. 初等函数的Maclaurin级数展开

本节系统推导五个基本初等函数的Maclaurin级数。对每个函数，我们给出其 Maclaurin 级数的具体形式、收敛半径，并利用定理 10.10 证明级数确实收敛于原函数。

### 4.1 指数函数 $e^x$

**展开式**：

$$e^x = \sum_{n=0}^{\infty} \frac{x^n}{n!} = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots + \frac{x^n}{n!} + \cdots, \quad x \in \mathbb{R} \ (R = +\infty)$$

**推导**：

**第 1 步：计算 Taylor 系数**。

对 $f(x) = e^x$，有 $f^{(n)}(x) = e^x$ 对一切 $n \ge 0$ 成立。因此：

$$f^{(n)}(0) = e^0 = 1, \quad \forall n \ge 0$$

Maclaurin 系数为：

$$a_n = \frac{f^{(n)}(0)}{n!} = \frac{1}{n!}$$

**第 2 步：写出 Maclaurin 级数**。

$$\sum_{n=0}^{\infty} \frac{f^{(n)}(0)}{n!} x^n = \sum_{n=0}^{\infty} \frac{x^n}{n!}$$

**第 3 步：求收敛半径**。

用比值法：

$$\lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{1/(n+1)!}{1/n!} = \lim_{n\to\infty} \frac{1}{n+1} = 0$$

因此收敛半径 $R = \infty$。

**第 4 步：证明级数收敛于 $e^x$**。

对任意固定的 $x \in \mathbb{R}$，Lagrange 余项为：

$$R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} x^{n+1} = \frac{e^{\xi}}{(n+1)!} x^{n+1}, \quad \xi \text{ 介于 } 0 \text{ 与 } x \text{ 之间}$$

由于 $\xi$ 介于 $0$ 与 $x$ 之间，$|\xi| \le |x|$，故 $e^{\xi} \le e^{|x|}$。因此：

$$|R_n(x)| \le \frac{e^{|x|} |x|^{n+1}}{(n+1)!}$$

由引理 10.1，$\lim_{n\to\infty} \frac{|x|^{n+1}}{(n+1)!} = 0$。因此 $\lim_{n\to\infty} |R_n(x)| = 0$，即 $\lim_{n\to\infty} R_n(x) = 0$。

由定理 10.10，Taylor 级数在 $x$ 处收敛于 $e^x$。由于 $x$ 是任意的，该结论对一切 $x \in \mathbb{R}$ 成立。

---

### 4.2 正弦函数 $\sin x$

**展开式**：

$$\sin x = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{(2n+1)!} = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \frac{x^7}{7!} + \cdots, \quad x \in \mathbb{R} \ (R = +\infty)$$

**推导**：

**第 1 步：计算 Taylor 系数**。

$f(x) = \sin x$ 的导数呈周期性变化：

$$f(x) = \sin x, \quad f'(x) = \cos x, \quad f''(x) = -\sin x, \quad f'''(x) = -\cos x, \quad f^{(4)}(x) = \sin x, \cdots$$

周期为 4。在 $x = 0$ 处取值：

$$f(0) = 0, \quad f'(0) = 1, \quad f''(0) = 0, \quad f'''(0) = -1, \quad f^{(4)}(0) = 0, \quad f^{(5)}(0) = 1, \cdots$$

一般地：

$$f^{(n)}(0) = \begin{cases}
0, & n = 2k \text{（偶数阶）} \\
(-1)^k, & n = 2k+1 \text{（奇数阶）}
\end{cases}$$

**第 2 步：写出 Maclaurin 级数**。

偶次项系数为零，仅奇次项非零。令 $n = 2k+1$：

$$a_{2k+1} = \frac{f^{(2k+1)}(0)}{(2k+1)!} = \frac{(-1)^k}{(2k+1)!}$$

因此：

$$\sum_{n=0}^{\infty} \frac{f^{(n)}(0)}{n!} x^n = \sum_{k=0}^{\infty} (-1)^k \frac{x^{2k+1}}{(2k+1)!}$$

**第 3 步：求收敛半径**。

考虑级数 $\sum_{k=0}^{\infty} \frac{x^{2k+1}}{(2k+1)!}$。令 $a_k = \frac{1}{(2k+1)!}$，则：

$$\lim_{k\to\infty} \frac{a_{k+1}}{a_k} = \lim_{k\to\infty} \frac{1/(2k+3)!}{1/(2k+1)!} = \lim_{k\to\infty} \frac{1}{(2k+2)(2k+3)} = 0$$

所以 $\sum \frac{x^{2k+1}}{(2k+1)!}$ 的收敛半径为 $+\infty$，乘以有界项 $(-1)^k$ 不改变收敛性，故 $\sin x$ 的 Maclaurin 级数收敛半径 $R = +\infty$。

**第 4 步：证明级数收敛于 $\sin x$**。

对任意固定的 $x \in \mathbb{R}$，Lagrange 余项为：

$$R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} x^{n+1}, \quad \xi \text{ 介于 } 0 \text{ 与 } x \text{ 之间}$$

由于 $\sin x$ 和 $\cos x$ 的任意阶导数都是 $\pm\sin x$ 或 $\pm\cos x$，其绝对值 $\le 1$。因此 $|f^{(n+1)}(\xi)| \le 1$，从而：

$$|R_n(x)| \le \frac{|x|^{n+1}}{(n+1)!}$$

由引理 10.1，$\lim_{n\to\infty} \frac{|x|^{n+1}}{(n+1)!} = 0$，故 $\lim_{n\to\infty} R_n(x) = 0$。

由定理 10.10，Taylor 级数在任意 $x \in \mathbb{R}$ 处收敛于 $\sin x$。

---

### 4.3 余弦函数 $\cos x$

**展开式**：

$$\cos x = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n}}{(2n)! = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \frac{x^6}{6!} + \cdots, \quad x \in \mathbb{R} \ (R = +\infty)$$

**推导**：

**第 1 步：计算 Taylor 系数**。

$f(x) = \cos x$ 的导数同样呈周期为 4 的循环：

$$f(0) = 1, \quad f'(0) = 0, \quad f''(0) = -1, \quad f'''(0) = 0, \quad f^{(4)}(0) = 1, \cdots$$

一般地：

$$f^{(n)}(0) = \begin{cases}
(-1)^{k}, & n = 2k \text{（偶数阶）} \\
0, & n = 2k+1 \text{（奇数阶）}
\end{cases}$$

**第 2 步：写出 Maclaurin 级数**。

奇次项系数为零，仅偶次项非零。令 $n = 2k$：

$$a_{2k} = \frac{f^{(2k)}(0)}{(2k)!} = \frac{(-1)^k}{(2k)!}$$

因此：

$$\sum_{n=0}^{\infty} \frac{f^{(n)}(0)}{n!} x^n = \sum_{k=0}^{\infty} (-1)^k \frac{x^{2k}}{(2k)!}$$

**第 3 步：收敛半径和收敛性**。

与 $\sin x$ 完全类似的分析可得 $R = +\infty$。Lagrange 余项同样满足 $|R_n(x)| \le \frac{|x|^{n+1}}{(n+1)!}$，故 $\lim_{n\to\infty} R_n(x) = 0$。

由定理 10.10，级数在 $\mathbb{R}$ 上收敛于 $\cos x$。

**验证**：对 $\sin x$ 的展开式逐项求导：

$$\frac{d}{dx} \sin x = \frac{d}{dx} \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{(2n+1)!} = \sum_{n=0}^{\infty} (-1)^n \frac{(2n+1)x^{2n}}{(2n+1)!} = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n}}{(2n)! = \cos x$$

这与 $\frac{d}{dx} \sin x = \cos x$ 一致，验证了展开式的正确性。

---

### 4.4 对数函数 $\ln(1+x)$

**展开式**：

$$\ln(1+x) = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n} = x - \frac{x^2}{2} + \frac{x^3}{3} - \frac{x^4}{4} + \cdots, \quad -1 < x \le 1 \ (R = 1)$$

**推导（方法一：直接求导计算Taylor系数）**：

**第 1 步：计算各阶导数**。

$f(x) = \ln(1+x)$，定义域 $x > -1$。

$$f'(x) = \frac{1}{1+x}, \quad f''(x) = -\frac{1}{(1+x)^2}, \quad f'''(x) = \frac{2!}{(1+x)^3}, \quad f^{(4)}(x) = -\frac{3!}{(1+x)^4}, \cdots$$

一般地，对 $n \ge 1$：

$$f^{(n)}(x) = (-1)^{n-1} \frac{(n-1)!}{(1+x)^n}$$

**第 2 步：在 $x = 0$ 处取值**。

$$f(0) = \ln 1 = 0, \quad f^{(n)}(0) = (-1)^{n-1} (n-1)! \ (n \ge 1)$$

**第 3 步：写出 Maclaurin 级数**。

$$a_0 = \frac{f(0)}{0!} = 0, \quad a_n = \frac{f^{(n)}(0)}{n!} = \frac{(-1)^{n-1} (n-1)!}{n!} = \frac{(-1)^{n-1}}{n} \ (n \ge 1)$$

因此：

$$\sum_{n=0}^{\infty} \frac{f^{(n)}(0)}{n!} x^n = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n}$$

**第 4 步：求收敛半径**。

用比值法（从 $n=1$ 开始）：

$$\lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{1/(n+1)}{1/n} = \lim_{n\to\infty} \frac{n}{n+1} = 1$$

因此收敛半径 $R = 1$。在端点处：
- $x = 1$：级数 $\sum_{n=1}^{\infty} (-1)^{n-1} \frac{1}{n}$ 是交错调和级数，由 Leibniz 判别法知收敛。
- $x = -1$：级数 $\sum_{n=1}^{\infty} (-1)^{n-1} \frac{(-1)^n}{n} = \sum_{n=1}^{\infty} -\frac{1}{n}$ 是负的调和级数，发散。

收敛域为 $(-1, 1]$。

**推导（方法二：从等比级数逐项积分）**：

已知 $\frac{1}{1+t} = \sum_{n=0}^{\infty} (-1)^n t^n$ 对 $|t| < 1$ 成立（等比级数）。将两端从 $0$ 到 $x$ 积分：

$$\int_0^x \frac{1}{1+t}\,dt = \sum_{n=0}^{\infty} (-1)^n \int_0^x t^n\,dt$$

左边 = $\ln(1+x)$，右边 = $\sum_{n=0}^{\infty} (-1)^n \frac{x^{n+1}}{n+1} = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n}$。

注意：逐项积分的合法性由定理 10.6（逐项积分）保证——等比级数在 $|t| \le r < 1$ 上内闭一致收敛。

---

### 4.5 反正切函数 $\arctan x$

**展开式**：

$$\arctan x = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{2n+1} = x - \frac{x^3}{3} + \frac{x^5}{5} - \frac{x^7}{7} + \cdots, \quad -1 \le x \le 1 \ (R = 1)$$

**推导**：

**第 1 步：从 $\frac{1}{1+x^2}$ 的展开出发**。

已知 $\frac{1}{1+u} = \sum_{n=0}^{\infty} (-1)^n u^n$ 对 $|u| < 1$ 成立。令 $u = x^2$，得：

$$\frac{1}{1+x^2} = \sum_{n=0}^{\infty} (-1)^n x^{2n}, \quad |x| < 1$$

**第 2 步：逐项积分**。

将两端从 $0$ 到 $x$ 积分：

$$\int_0^x \frac{1}{1+t^2}\,dt = \sum_{n=0}^{\infty} (-1)^n \int_0^x t^{2n}\,dt$$

左边 = $\arctan x$（因为 $\frac{d}{dx} \arctan x = \frac{1}{1+x^2}$ 且 $\arctan 0 = 0$）。

右边 = $\sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{2n+1}$。

**第 3 步：收敛域**。

用比值法。级数 $\sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{2n+1}$ 的通项为 $u_n(x) = (-1)^n \frac{x^{2n+1}}{2n+1}$。考察绝对值级数：

$$\lim_{n\to\infty} \left|\frac{u_{n+1}(x)}{u_n(x)}\right| = \lim_{n\to\infty} \frac{|x|^{2n+3}/(2n+3)}{|x|^{2n+1}/(2n+1)} = |x|^2 \lim_{n\to\infty} \frac{2n+1}{2n+3} = |x|^2$$

当 $|x| < 1$ 时，$|x|^2 < 1$，绝对收敛；当 $|x| > 1$ 时，$|x|^2 > 1$，发散。故 $R = 1$。

端点处：
- $x = 1$：级数 $\sum_{n=0}^{\infty} (-1)^n \frac{1}{2n+1} = 1 - \frac{1}{3} + \frac{1}{5} - \frac{1}{7} + \cdots$ 是交错级数，通项单调趋于 $0$，由 Leibniz 判别法知收敛。此时 $\arctan 1 = \frac{\pi}{4}$，所以 $\frac{\pi}{4} = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1}$——这就是 Leibniz 公式。
- $x = -1$：同理收敛。

收敛域为 $[-1, 1]$。

### 4.6 标准展开汇总表

| 函数 | Maclaurin 级数 | 收敛半径 | 收敛域 |
|------|---------------|---------|-------|
| $e^x$ | $\displaystyle\sum_{n=0}^{\infty} \frac{x^n}{n!}$ | $+\infty$ | $\mathbb{R}$ |
| $\sin x$ | $\displaystyle\sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{(2n+1)!}$ | $+\infty$ | $\mathbb{R}$ |
| $\cos x$ | $\displaystyle\sum_{n=0}^{\infty} (-1)^n \frac{x^{2n}}{(2n)!} | $+\infty$ | $\mathbb{R}$ |
| $\ln(1+x)$ | $\displaystyle\sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n}$ | $1$ | $(-1, 1]$ |
| $\arctan x$ | $\displaystyle\sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{2n+1}$ | $1$ | $[-1, 1]$ |
| $\frac{1}{1-x}$ | $\displaystyle\sum_{n=0}^{\infty} x^n$ | $1$ | $(-1, 1)$ |
| $\frac{1}{1+x}$ | $\displaystyle\sum_{n=0}^{\infty} (-1)^n x^n$ | $1$ | $(-1, 1)$ |

注：最后两行来自等比级数（文件 03 的基础知识），在此一并列出便于后续使用。

---

## 5. 广义二项式定理

### 5.1 定义与定理

**定理 10.11（广义二项式定理）**：设 $\alpha \in \mathbb{R}$ 为任意实数，定义**广义二项式系数**（generalized binomial coefficient）：

$$\binom{\alpha}{n} = \frac{\alpha(\alpha-1)(\alpha-2)\cdots(\alpha-n+1)}{n!}, \quad n \in \mathbb{N}$$

并约定 $\binom{\alpha}{0} = 1$。则函数 $(1+x)^\alpha$ 的 Maclaurin 级数为：

$$(1+x)^\alpha = \sum_{n=0}^{\infty} \binom{\alpha}{n} x^n, \quad |x| < 1$$

当 $\alpha \notin \mathbb{N}$（非负整数）时，收敛半径 $R = 1$。当 $\alpha \in \mathbb{N}$ 时，级数退化为有限和（多项式），此时 $R = +\infty$。

**注**：当 $\alpha = m$ 为正整数时，对 $n > m$，因式 $\alpha - n + 1$ 中会出现 $0$，故 $\binom{m}{n} = 0$ 对 $n > m$ 成立，级数退化为普通二项式定理的有限和：

$$(1+x)^m = \sum_{n=0}^{m} \binom{m}{n} x^n$$

### 5.2 定理的推导（Taylor系数法）

**第 1 步：计算 $(1+x)^\alpha$ 的各阶导数**。

设 $f(x) = (1+x)^\alpha$（$\alpha$ 为任意实数，$x > -1$ 以保证 $(1+x)^\alpha$ 有定义）。

$$f'(x) = \alpha(1+x)^{\alpha-1}$$
$$f''(x) = \alpha(\alpha-1)(1+x)^{\alpha-2}$$
$$f'''(x) = \alpha(\alpha-1)(\alpha-2)(1+x)^{\alpha-3}$$

一般地，对 $n \ge 1$：

$$f^{(n)}(x) = \alpha(\alpha-1)(\alpha-2)\cdots(\alpha-n+1)(1+x)^{\alpha-n}$$

**第 2 步：在 $x = 0$ 处取值**。

$$f(0) = 1^\alpha = 1$$
$$f^{(n)}(0) = \alpha(\alpha-1)(\alpha-2)\cdots(\alpha-n+1)$$

**第 3 步：写出 Taylor 系数**。

$$a_0 = \frac{f(0)}{0!} = 1 = \binom{\alpha}{0}$$

对 $n \ge 1$：

$$a_n = \frac{f^{(n)}(0)}{n!} = \frac{\alpha(\alpha-1)(\alpha-2)\cdots(\alpha-n+1)}{n!} = \binom{\alpha}{n}$$

因此 Maclaurin 级数为 $\sum_{n=0}^{\infty} \binom{\alpha}{n} x^n$。

**第 4 步：求收敛半径**。

对 $n \ge 1$，计算相邻系数的比值：

$$\left|\frac{a_{n+1}}{a_n}\right| = \left|\frac{\binom{\alpha}{n+1}}{\binom{\alpha}{n}}\right| = \left|\frac{\alpha(\alpha-1)\cdots(\alpha-n+1)(\alpha-n)}{(n+1)!} \cdot \frac{n!}{\alpha(\alpha-1)\cdots(\alpha-n+1)}\right| = \left|\frac{\alpha-n}{n+1}\right|$$

当 $n \to \infty$ 时：

$$\lim_{n\to\infty} \left|\frac{a_{n+1}}{a_n}\right| = \lim_{n\to\infty} \frac{|\alpha-n|}{n+1} = \lim_{n\to\infty} \frac{n-\alpha}{n+1} = 1$$

由比值法，收敛半径 $R = 1$。

**第 5 步：收敛性证明**（Lagrange余项法）。

此处不展开详细的 Lagrange 余项证明（涉及对 $(1+\xi)^{\alpha-n-1}$ 的估计），接受结论：对 $|x| < 1$，级数收敛于 $(1+x)^\alpha$。

### 5.3 $\alpha = 1/2$ 特例

最常用的特例之一是 $\alpha = \frac{1}{2}$，即 $\sqrt{1+x}$ 的展开。

**广义二项式系数**：

$$a_0 = \binom{1/2}{0} = 1$$

$$a_1 = \binom{1/2}{1} = \frac{1/2}{1!} = \frac{1}{2}$$

$$a_2 = \binom{1/2}{2} = \frac{(1/2)(-1/2)}{2!} = -\frac{1}{8}$$

$$a_3 = \binom{1/2}{3} = \frac{(1/2)(-1/2)(-3/2)}{3!} = -\frac{1}{16}$$

$$a_4 = \binom{1/2}{4} = \frac{(1/2)(-1/2)(-3/2)(-5/2)}{4!} = -\frac{15/16}{24} = -\frac{5}{128}$$

通项公式：对 $n \ge 1$，

$$a_n = \frac{(1/2)(-1/2)(-3/2)\cdots((3-2n)/2)}{n!} = (-1)^{n-1} \frac{3 \cdot 5 \cdots (2n-3)}{2^n n!}$$

因此：

$$\sqrt{1+x} = 1 + \frac{1}{2}x - \frac{1}{8}x^2 - \frac{1}{16}x^3 - \frac{5}{128}x^4 + \cdots, \quad |x| < 1$$

### 5.4 其他常见特例

$\alpha = -1$：
$$(1+x)^{-1} = \sum_{n=0}^{\infty} \binom{-1}{n} x^n = \sum_{n=0}^{\infty} (-1)^n x^n = 1 - x + x^2 - x^3 + \cdots$$

验证：$\binom{-1}{n} = \frac{(-1)(-2)\cdots(-n)}{n!} = (-1)^n$，与等比级数一致。

$\alpha = -\frac{1}{2}$：

$$\frac{1}{\sqrt{1+x}} = \sum_{n=0}^{\infty} \binom{-1/2}{n} x^n = 1 - \frac{1}{2}x + \frac{3}{8}x^2 - \frac{5}{16}x^3 + \cdots$$

### 5.5 广义二项式系数与普通二项式系数的关系

当 $\alpha = m$ 为正整数时：

$$\binom{m}{n} = \frac{m(m-1)\cdots(m-n+1)}{n!} = \frac{m!}{n!(m-n)!}$$

这就是中学阶段学过的普通二项式系数。广义二项式系数是它的自然推广——将分子中的正整数 $m$ 替换为任意实数 $\alpha$。

---

## 6. 技术引理：$\sqrt[n]{n!} \to \infty$

在处理含 $n!$ 的系数时（如 $\sum \frac{x^n}{n!}$），Cauchy-Hadamard 公式需要计算 $\varlimsup \sqrt[n]{|a_n|}$。对于 $a_n = 1/n!$，我们需要 $\sqrt[n]{n!}$ 的极限信息。

**引理 10.2**：
$$\lim_{n\to\infty} \sqrt[n]{n!} = +\infty$$

**证明**：采用下界放缩法。

**第 1 步：将 $n!$ 拆分为前后两半**。

对 $n \ge 2$，$n! = 1 \times 2 \times \cdots \times n$。将其分为前 $\lfloor n/2 \rfloor$ 个因子和后 $n - \lfloor n/2 \rfloor$ 个因子。

为简化，考虑 $n$ 为偶数的情形（$n = 2k$）：

$$(2k)! = 1 \times 2 \times \cdots \times k \times (k+1) \times \cdots \times (2k)$$

后 $k$ 个因子中，每个都 $\ge k$，因此：

$$(2k)! \ge (k+1)(k+2)\cdots(2k) \ge k \cdot k \cdots k = k^k$$

**第 2 步：推广到一般 $n$**。

对任意 $n \ge 2$，$n!$ 的后 $\lfloor n/2 \rfloor$ 个因子的最小值是 $\lceil n/2 \rceil$，因此：

$$n! \ge \left(\frac{n}{2}\right)^{\lfloor n/2 \rfloor} \ge \left(\frac{n}{2}\right)^{n/2 - 1}$$

对 $n \ge 4$，更简洁的估计：

$$n! \ge \left(\frac{n}{2}\right)^{n/2}$$

（后一半因子每个 $\ge \frac{n}{2}$，共 $\frac{n}{2}$ 个）。严格推导：对 $n \ge 2$，$\lfloor n/2 \rfloor \ge 1$，$n!$ 的后 $\lfloor n/2 \rfloor$ 个因子分别为 $\lfloor n/2 \rfloor + 1, \lfloor n/2 \rfloor + 2, \ldots, n$，每个都 $\ge \lceil n/2 \rceil \ge n/2$。

**第 3 步：开 $n$ 次方**。

$$\sqrt[n]{n!} \ge \sqrt[n]{\left(\frac{n}{2}\right)^{n/2}} = \left(\frac{n}{2}\right)^{1/2} = \sqrt{\frac{n}{2}}$$

**第 4 步：取极限**。

当 $n \to \infty$ 时，$\sqrt{\frac{n}{2}} \to \infty$，因此 $\sqrt[n]{n!} \to \infty$。

证毕。

**应用**：用 Cauchy-Hadamard 公式求 $\sum \frac{x^n}{n!}$ 的收敛半径。

$$\varlimsup_{n\to\infty} \sqrt[n]{\left|\frac{1}{n!}\right|} = \frac{1}{\lim_{n\to\infty} \sqrt[n]{n!}} = \frac{1}{\infty} = 0$$

因此 $R = 1/0 = +\infty$，与比值法结论一致。

---

## 7. 例题

### 例题 1：利用已知展开求新级数

求函数 $f(x) = \frac{x}{1-x^2}$ 的 Maclaurin 级数，并求收敛半径。

**解**：

**第 1 步：将函数分解为已知展开的组合**。

已知 $\frac{1}{1-u} = \sum_{n=0}^{\infty} u^n$ 对 $|u| < 1$ 成立。令 $u = x^2$，得：

$$\frac{1}{1-x^2} = \sum_{n=0}^{\infty} (x^2)^n = \sum_{n=0}^{\infty} x^{2n}, \quad |x| < 1$$

**第 2 步：乘以 $x$**。

$$f(x) = x \cdot \frac{1}{1-x^2} = x \cdot \sum_{n=0}^{\infty} x^{2n} = \sum_{n=0}^{\infty} x^{2n+1}$$

**第 3 步：确定收敛半径**。

$\frac{1}{1-x^2}$ 的收敛半径为 $R=1$（在 $|x|<1$ 内收敛）。乘以 $x$ 不改变收敛半径，故 $f(x)$ 的收敛半径也是 $R=1$。

**验证**：$f(x) = \frac{x}{1-x^2}$ 的奇点是 $x = \pm 1$（分母为零），因此确实 $R=1$。

**第 4 步：写出最终结果**。

$$\frac{x}{1-x^2} = \sum_{n=0}^{\infty} x^{2n+1} = x + x^3 + x^5 + x^7 + \cdots, \quad |x| < 1$$

---

### 例题 2：通过逐项积分求 $\ln(1+x)$ 的展开

已知 $\frac{1}{1+x} = \sum_{n=0}^{\infty} (-1)^n x^n$（$|x| < 1$），用逐项积分法求 $\ln(1+x)$ 的 Maclaurin 级数。

**解**：

**第 1 步：验证逐项积分的条件**。

对任意 $r \in (0, 1)$，在 $[0, r]$ 上，$\sum_{n=0}^{\infty} (-1)^n t^n$ 满足：
- 每项 $u_n(t) = (-1)^n t^n$ 连续；
- $|u_n(t)| \le r^n$，$\sum r^n$ 收敛（等比级数，公比 $r < 1$）。

由 Weierstrass M-判别法，级数在 $[0, r]$ 上一致收敛。

**第 2 步：应用逐项积分定理（定理 10.6）**。

对 $x \in (0, r)$：

$$\int_0^x \frac{1}{1+t}\,dt = \sum_{n=0}^{\infty} \int_0^x (-1)^n t^n\,dt$$

**第 3 步：计算两边**。

左边：$\int_0^x \frac{1}{1+t}\,dt = \left[\ln(1+t)\right]_0^x = \ln(1+x)$

右边：$\sum_{n=0}^{\infty} (-1)^n \int_0^x t^n\,dt = \sum_{n=0}^{\infty} (-1)^n \frac{x^{n+1}}{n+1} = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n}$

**第 4 步：结论**。

$$\ln(1+x) = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{x^n}{n}, \quad |x| < 1$$

（端点 $x = 1$ 处的收敛性需单独判断——由 Leibniz 判别法知收敛。）

---

### 例题 3：利用广义二项式定理展开

求 $f(x) = \sqrt{1-x}$ 的 Maclaurin 级数（写出前四项及通项）。

**解**：

**第 1 步：将函数形式与广义二项式定理匹配**。

$$f(x) = \sqrt{1-x} = (1 + (-x))^{1/2}$$

令 $u = -x$，$\alpha = 1/2$，则 $f(x) = (1+u)^{1/2}$。

**第 2 步：写出广义二项式展开**。

由定理 10.11：

$$(1+u)^{1/2} = \sum_{n=0}^{\infty} \binom{1/2}{n} u^n, \quad |u| < 1$$

**第 3 步：代入 $u = -x$ 并写出前四项**。

$$\sqrt{1-x} = \sum_{n=0}^{\infty} \binom{1/2}{n} (-x)^n = \sum_{n=0}^{\infty} \binom{1/2}{n} (-1)^n x^n$$

计算前四项的系数：

- $n = 0$：$\binom{1/2}{0} (-1)^0 = 1$
- $n = 1$：$\binom{1/2}{1} (-1)^1 = \frac{1}{2} \cdot (-1) = -\frac{1}{2}$
- $n = 2$：$\binom{1/2}{2} (-1)^2 = \left(-\frac{1}{8}\right) \cdot 1 = -\frac{1}{8}$
- $n = 3$：$\binom{1/2}{3} (-1)^3 = \frac{1}{16} \cdot (-1) = -\frac{1}{16}$

**第 4 步：写出通项**。

对 $n \ge 1$：

$$a_n = \binom{1/2}{n} (-1)^n = \frac{(1/2)(-1/2)(-3/2)\cdots((3-2n)/2)}{n!} \cdot (-1)^n$$
$$= \frac{(-1)^{n-1} \cdot 1 \cdot 3 \cdot 5 \cdots (2n-3)}{2^n n!} \cdot (-1)^n$$
$$= -\frac{1 \cdot 3 \cdot 5 \cdots (2n-3)}{2^n n!}$$

**结果**：

$$\sqrt{1-x} = 1 - \frac{1}{2}x - \frac{1}{8}x^2 - \frac{1}{16}x^3 - \frac{5}{128}x^4 - \cdots, \quad |x| < 1$$

---

## 8. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 只要函数无穷次可导，其Taylor级数就一定收敛于该函数 | 无穷次可导只能保证Taylor级数作为形式幂级数存在，但不能保证收敛性，更不能保证收敛于原函数。经典反例：$f(x) = e^{-1/x^2}$（$x \neq 0$），$f(0) = 0$，在 $x=0$ 处无穷次可导，但所有导数均为 $0$，Taylor级数为 $0 \neq f(x)$（$x \neq 0$） |
| Taylor级数收敛 $\iff$ Taylor级数收敛于原函数 | 两个概念不同。例如 $f(x) = e^{-1/x^2}$ 的Taylor级数是 $0$，它在全实数轴上收敛，但仅收敛于 $0$ 而非原函数。收敛于原函数需要额外验证 $\lim R_n(x) = 0$ |
| $\lim_{n\to\infty} \frac{x^n}{n!} = 0$ 的证明可以用 $n! \sim \sqrt{2\pi n}(n/e)^n$（Stirling公式） | Stirling公式只是给出了更精确的估计，但基础证明只需要初等放缩：取 $N > 2|x|$，则 $\frac{|x|^n}{n!} \le \frac{|x|^N}{N!} \cdot (\frac{1}{2})^{n-N}$ |
| $(1+x)^\alpha$ 的Maclaurin级数收敛半径总是 $R=1$ | 当 $\alpha$ 为非负整数时，级数退化为有限多项式，此时 $R = +\infty$。只有非整数的 $\alpha$ 才满足 $R=1$ |
| $\ln(1+x)$ 的麦克劳林级数在 $x = -1$ 处发散，所以 $\ln 0$ 无定义 | 发散的原因是右端级数 $\sum (-1)^{n-1} (-1)^n/n = -\sum 1/n$ 发散，而左端 $\ln 0$ 本身趋向 $-\infty$，两者都发散到无穷——逻辑上一致，但 $\ln(1+x)$ 在 $x=-1$ 处无定义 |
| 广义二项式系数 $\binom{\alpha}{n}$ 对 $\alpha \notin \mathbb{N}$ 保持交替符号 | 符号取决于 $\alpha$ 的取值。例如 $\alpha = 1/2$ 时，$a_1 > 0$（正），$a_2 < 0$（负），之后交替——因为分子中每个新增因子都增加一个负号。对 $\alpha > 0$，符号规律需要具体分析 |

### 检查点

- [ ] 能否准确说出Taylor级数与Taylor多项式的区别？能否写出两者的形式化定义？
- [ ] 能否写出定理10.10（Taylor级数收敛于原函数的充要条件）并给出证明？
- [ ] 能否独立完成 $\lim_{n\to\infty} x^n/n! = 0$ 的完整证明（取 $N > 2|x|$ 的放缩法）？
- [ ] 能否默写出 $e^x$、$\sin x$、$\cos x$、$\ln(1+x)$、$\arctan x$ 的Maclaurin级数及收敛半径？
- [ ] 能否用两种方法推导 $\ln(1+x)$ 的Maclaurin级数（直接求Taylor系数 + 逐项积分）？
- [ ] 能否写出广义二项式系数的定义，并计算 $\binom{1/2}{3}$ 和 $\binom{-1/3}{2}$？
- [ ] 能否推导 $(1+x)^\alpha$ 的Maclaurin级数并确定收敛半径？
- [ ] 能否证明 $\lim_{n\to\infty} \sqrt[n]{n!} = +\infty$？
- [ ] 给定一个函数（如 $\frac{x}{1+x^2}$），能否通过已知展开的变量替换或代数操作求出其Maclaurin级数？

---

## 练习题

### 基础巩固

**1.** 求下列函数的Maclaurin级数（利用已知展开进行变量替换、逐项积分或逐项求导）：

(1) $f(x) = e^{-x^2}$
(2) $f(x) = \frac{x}{1+x^2}$
(3) $f(x) = \ln(1-x)$

<details><summary>参考答案</summary>

**(1)** $f(x) = e^{-x^2}$

已知 $e^u = \sum_{n=0}^{\infty} \frac{u^n}{n!}$ 对 $u \in \mathbb{R}$ 成立。令 $u = -x^2$：

$$e^{-x^2} = \sum_{n=0}^{\infty} \frac{(-x^2)^n}{n!} = \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n}}{n!}, \quad x \in \mathbb{R}$$

收敛半径 $R = \infty$（因为 $e^u$ 的展开对一切 $u \in \mathbb{R}$ 成立）。

**(2)** $f(x) = \frac{x}{1+x^2}$

已知 $\frac{1}{1+u} = \sum_{n=0}^{\infty} (-1)^n u^n$ 对 $|u| < 1$ 成立。令 $u = x^2$：

$$\frac{1}{1+x^2} = \sum_{n=0}^{\infty} (-1)^n x^{2n}, \quad |x| < 1$$

乘以 $x$：

$$f(x) = x \cdot \frac{1}{1+x^2} = \sum_{n=0}^{\infty} (-1)^n x^{2n+1}, \quad |x| < 1$$

收敛半径 $R = 1$（奇点在 $x = \pm i$ 处）。

**(3)** $f(x) = \ln(1-x)$

方法一（直接求Taylor系数）：

$f(x) = \ln(1-x)$，$f'(x) = -\frac{1}{1-x}$，$f''(x) = -\frac{1}{(1-x)^2}$。

$f(0) = 0$，$f^{(n)}(0) = -(n-1)!$ 对 $n \ge 1$。

$$a_n = \frac{f^{(n)}(0)}{n!} = -\frac{1}{n}, \quad n \ge 1$$

因此 $\ln(1-x) = -\sum_{n=1}^{\infty} \frac{x^n}{n}$，$|x| < 1$。

方法二（从已知展开出发）：

已知 $\ln(1+u) = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{u^n}{n}$，令 $u = -x$：

$$\ln(1-x) = \sum_{n=1}^{\infty} (-1)^{n-1} \frac{(-x)^n}{n} = \sum_{n=1}^{\infty} (-1)^{n-1} (-1)^n \frac{x^n}{n} = -\sum_{n=1}^{\infty} \frac{x^n}{n}$$

两种方法结果一致。收敛半径 $R = 1$，收敛域 $[-1, 1)$。

</details>

---

**2.** 用广义二项式定理展开 $f(x) = \frac{1}{\sqrt{1-x}}$，写出前四项和通项公式。

<details><summary>参考答案</summary>

**第 1 步：匹配形式**。

$$\frac{1}{\sqrt{1-x}} = (1-x)^{-1/2} = (1+u)^{-1/2}, \quad u = -x, \ \alpha = -\frac{1}{2}$$

**第 2 步：计算广义二项式系数**。

$$a_0 = \binom{-1/2}{0} = 1$$

$$a_1 = \binom{-1/2}{1} (-1)^1 = \left(-\frac{1}{2}\right)(-1) = \frac{1}{2}$$

$$a_2 = \binom{-1/2}{2} (-1)^2 = \frac{(-1/2)(-3/2)}{2!} \cdot 1 = \frac{3}{8}$$

$$a_3 = \binom{-1/2}{3} (-1)^3 = \frac{(-1/2)(-3/2)(-5/2)}{3!} \cdot (-1) = \frac{5}{16}$$

**第 3 步：通项公式**。

对 $n \ge 1$：

$$a_n = \binom{-1/2}{n} (-1)^n = \frac{(-1/2)(-3/2)(-5/2)\cdots(-(2n-1)/2)}{n!} \cdot (-1)^n$$

注意分子有 $n$ 个因子，每个负号：

$$\binom{-1/2}{n} = \frac{(-1)^n \cdot 1 \cdot 3 \cdot 5 \cdots (2n-1)}{2^n n!}$$

再乘以 $(-1)^n$ 得：

$$a_n = \frac{1 \cdot 3 \cdot 5 \cdots (2n-1)}{2^n n!}$$

**结果**：

$$\frac{1}{\sqrt{1-x}} = 1 + \frac{1}{2}x + \frac{3}{8}x^2 + \frac{5}{16}x^3 + \cdots + \frac{1 \cdot 3 \cdot 5 \cdots (2n-1)}{2^n n!} x^n + \cdots, \quad |x| < 1$$

</details>

---

### 迁移应用

**3.** 利用已知展开求 $\displaystyle\frac{1}{1+x^2}$ 的 Maclaurin 级数，然后通过逐项积分推导 $\arctan x$ 的展开式。利用这个展开式计算 $\pi$ 的近似值（提示：$\arctan 1 = \pi/4$）。

<details><summary>参考答案</summary>

**第 1 步：求 $\frac{1}{1+x^2}$ 的展开**。

由 $\frac{1}{1+u} = \sum_{n=0}^{\infty} (-1)^n u^n$ 对 $|u| < 1$ 成立，令 $u = x^2$：

$$\frac{1}{1+x^2} = \sum_{n=0}^{\infty} (-1)^n x^{2n}, \quad |x| < 1$$

**第 2 步：逐项积分**。

对任意 $r \in (0, 1)$，在 $[0, r]$ 上级数一致收敛（M-判别法：$|(-1)^n t^{2n}| \le r^{2n}$）。由定理 10.6：

$$\int_0^x \frac{1}{1+t^2}\,dt = \sum_{n=0}^{\infty} (-1)^n \int_0^x t^{2n}\,dt$$

左边 $= \arctan x$，右边 $= \sum_{n=0}^{\infty} (-1)^n \frac{x^{2n+1}}{2n+1}$。

**第 3 步：计算 $\pi$ 的近似值**。

取 $x = 1$：

$$\frac{\pi}{4} = \arctan 1 = \sum_{n=0}^{\infty} (-1)^n \frac{1}{2n+1} = 1 - \frac{1}{3} + \frac{1}{5} - \frac{1}{7} + \cdots$$

于是 $\pi = 4\left(1 - \frac{1}{3} + \frac{1}{5} - \frac{1}{7} + \cdots\right)$。

**第 4 步：数值举例**。

取前 4 项：$\pi \approx 4\left(1 - \frac{1}{3} + \frac{1}{5} - \frac{1}{7}\right) = 4 \times \frac{76}{105} \approx 2.8952$（误差较大）。

取前 10 项：$\pi \approx 4\sum_{n=0}^{9} \frac{(-1)^n}{2n+1} \approx 3.0418$。

这个级数收敛速度较慢（因为它是交错调和级数型的，误差约 $4/(2N+1)$），需要很多项才能达到高精度。实际计算 $\pi$ 常用收敛更快的 Machin 公式等。

</details>

---

**4.** 设 $f(x) = \frac{x}{(1-x)^2}$，求 $f(x)$ 的 Maclaurin 级数。

<details><summary>参考答案</summary>

**方法一（代数分解 + 已知展开）**：

$\frac{1}{(1-x)^2}$ 可以从两个角度得到：

角度 1：对 $\frac{1}{1-x} = \sum_{n=0}^{\infty} x^n$ 逐项求导。由定理 10.7，对 $|x| < 1$：

$$\frac{d}{dx} \frac{1}{1-x} = \frac{1}{(1-x)^2} = \frac{d}{dx} \sum_{n=0}^{\infty} x^n = \sum_{n=1}^{\infty} n x^{n-1} = \sum_{n=0}^{\infty} (n+1) x^n$$

角度 2：用已知展开 $\frac{1}{(1-x)^2} = \sum_{n=0}^{\infty} (n+1) x^n$（可用 Cauchy 乘积验证：$\frac{1}{1-x} \cdot \frac{1}{1-x} = (\sum x^n)(\sum x^n) = \sum (n+1)x^n$）。

**第 2 步：乘以 $x$**。

$$f(x) = x \cdot \frac{1}{(1-x)^2} = x \cdot \sum_{n=0}^{\infty} (n+1) x^n = \sum_{n=0}^{\infty} (n+1) x^{n+1} = \sum_{n=1}^{\infty} n x^{n}$$

**结果**：

$$\frac{x}{(1-x)^2} = \sum_{n=1}^{\infty} n x^n = x + 2x^2 + 3x^3 + 4x^4 + \cdots, \quad |x| < 1$$

**验证**：$x = \frac{1}{2}$ 时，左边 $\frac{1/2}{(1-1/2)^2} = 2$，右边 $\sum_{n=1}^{\infty} n/2^n = 2$（这是已知的数项级数结论）。

</details>
