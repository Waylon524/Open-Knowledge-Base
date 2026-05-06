# 第十六章 Fourier级数

---

## §3 Fourier级数的分析性质——逐项积分与逐项求导

---

### 3.1 Fourier级数的逐项积分定理

在§1中我们通过Euler-Fourier公式计算了函数 $x^2$ 的Fourier系数，直接得到了其Fourier级数。一个自然的问题是：能否从 $x$ 的Fourier级数出发，通过逐项积分获得 $x^2$ 的展开？更一般地，对Fourier级数进行逐项积分是否总是合法的？

#### 定理 16.17（Fourier级数的逐项积分定理）

设 $f$ 是以 $2\pi$ 为周期的函数，在 $[-\pi,\pi]$ 上分段连续（至多有限个第一类间断点）。记其Fourier级数为
$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} \bigl(a_n \cos nx + b_n \sin nx\bigr).
$$

定义
$$
F(x) = \int_{0}^{x} \Bigl[f(t) - \frac{a_0}{2}\Bigr]\,dt,
$$
则对任意 $x \in [-\pi,\pi]$，成立
$$
\boxed{\;
F(x) = \sum_{n=1}^{\infty} \frac{a_n \sin nx + b_n(1 - \cos nx)}{n},
\;}
$$
且右端级数在 $[-\pi,\pi]$ 上 **绝对且一致收敛**。

换言之，$F$ 的Fourier级数可以通过对 $f$ 的Fourier级数（减去常数项后）逐项积分得到，且积分后的级数收敛性从逐点提升为一致收敛。

---

**证明** 分三步进行。

**第一步：$F$ 的周期性和连续性。** 由 $a_0$ 的定义，
$$
\int_{-\pi}^{\pi} \Bigl[f(t) - \frac{a_0}{2}\Bigr]\,dt = \int_{-\pi}^{\pi} f(t)\,dt - \pi a_0 = 0.
$$

因此
$$
\begin{aligned}
F(x+2\pi) - F(x) &= \int_{x}^{x+2\pi} \Bigl[f(t) - \frac{a_0}{2}\Bigr]\,dt \\
&= \int_{-\pi}^{\pi} \Bigl[f(t) - \frac{a_0}{2}\Bigr]\,dt = 0,
\end{aligned}
$$
故 $F$ 是以 $2\pi$ 为周期的 **连续函数**。特别地，$F(-\pi) = F(\pi)$。同时，在 $f$ 的连续点处 $F'(x) = f(x) - a_0/2$，在 $f$ 的间断点处 $F$ 仍连续（积分不依赖于单点值）。故 $F$ 在 $[-\pi,\pi]$ 上分段 $C^1$，满足Dirichlet-Jordan收敛定理的条件，其Fourier级数处处逐点收敛到 $F(x)$。

**第二步：计算 $F$ 的Fourier系数。** 记 $F$ 的Fourier系数为 $A_n, B_n$：
$$
A_n = \frac{1}{\pi}\int_{-\pi}^{\pi} F(x)\cos nx\,dx,\qquad
B_n = \frac{1}{\pi}\int_{-\pi}^{\pi} F(x)\sin nx\,dx.
$$

对 $n \ge 1$，利用分部积分（$F$ 连续且 $F(-\pi)=F(\pi)$，边界项为零）：
$$
\begin{aligned}
A_n &= \frac{1}{\pi}\int_{-\pi}^{\pi} F(x)\cos nx\,dx \\
&= \frac{1}{\pi}\Big( \big[F(x)\frac{\sin nx}{n}\big]_{-\pi}^{\pi} - \frac{1}{n}\int_{-\pi}^{\pi} F'(x)\sin nx\,dx \Big) \\
&= -\frac{1}{n\pi}\int_{-\pi}^{\pi} \big[f(x) - \frac{a_0}{2}\big]\sin nx\,dx \\
&= -\frac{1}{n\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx \qquad \big(\text{因为 } \tfrac{a_0}{2}\sin nx \text{ 的积分为零}\big) \\
&= -\frac{b_n}{n}.
\end{aligned}
$$

同理，
$$
\begin{aligned}
B_n &= \frac{1}{\pi}\int_{-\pi}^{\pi} F(x)\sin nx\,dx \\
&= \frac{1}{\pi}\Big( \big[-F(x)\frac{\cos nx}{n}\big]_{-\pi}^{\pi} + \frac{1}{n}\int_{-\pi}^{\pi} F'(x)\cos nx\,dx \Big) \\
&= \frac{1}{n\pi}\int_{-\pi}^{\pi} \big[f(x) - \frac{a_0}{2}\big]\cos nx\,dx \\
&= \frac{1}{n\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx - \frac{a_0}{2n\pi}\int_{-\pi}^{\pi}\cos nx\,dx \\
&= \frac{a_n}{n} \qquad \big(\text{因为 } \int_{-\pi}^{\pi}\cos nx\,dx = 0,\ n\ge 1\big).
\end{aligned}
$$

因此，对 $n \ge 1$ 有
$$
A_n = -\frac{b_n}{n},\qquad B_n = \frac{a_n}{n}.
$$

**第三步：确定常数项 $A_0$。** $F$ 的Fourier级数为
$$
F(x) \sim \frac{A_0}{2} + \sum_{n=1}^{\infty} (A_n\cos nx + B_n\sin nx)
    = \frac{A_0}{2} + \sum_{n=1}^{\infty} \Bigl(-\frac{b_n}{n}\cos nx + \frac{a_n}{n}\sin nx\Bigr).
$$

由 $F(0)=0$ 代入上式：
$$
0 = F(0) = \frac{A_0}{2} + \sum_{n=1}^{\infty} \Bigl(-\frac{b_n}{n}\Bigr) \quad\Longrightarrow\quad \frac{A_0}{2} = \sum_{n=1}^{\infty} \frac{b_n}{n}.
$$

因此
$$
\begin{aligned}
F(x) &= \frac{A_0}{2} + \sum_{n=1}^{\infty} \Bigl(-\frac{b_n}{n}\cos nx + \frac{a_n}{n}\sin nx\Bigr) \\
&= \sum_{n=1}^{\infty} \frac{b_n}{n} + \sum_{n=1}^{\infty} \Bigl(-\frac{b_n}{n}\cos nx + \frac{a_n}{n}\sin nx\Bigr) \\
&= \sum_{n=1}^{\infty} \frac{a_n\sin nx + b_n(1 - \cos nx)}{n}.
\end{aligned}
$$

这正是对 $f$ 的Fourier级数（减去常数项）从 $0$ 到 $x$ 逐项积分的结果：
$$
\int_{0}^{x} \Bigl[f(t) - \frac{a_0}{2}\Bigr]\,dt
\sim \sum_{n=1}^{\infty} \int_{0}^{x} \bigl(a_n\cos nt + b_n\sin nt\bigr)\,dt
= \sum_{n=1}^{\infty} \frac{a_n\sin nx + b_n(1 - \cos nx)}{n}.
$$

**一致收敛性**：对 $n \ge 1$，
$$
\Bigl|\frac{a_n\sin nx + b_n(1 - \cos nx)}{n}\Bigr| \le \frac{|a_n| + 2|b_n|}{n}.
$$

由Bessel不等式（定理16.15），$\sum (a_n^2 + b_n^2)$ 收敛，故 $\sum |a_n|/n$ 和 $\sum |b_n|/n$ 收敛（利用Cauchy-Schwarz不等式）：
$$
\sum_{n=1}^{\infty} \frac{|a_n|}{n} \le \Bigl(\sum_{n=1}^{\infty} a_n^2\Bigr)^{1/2} \Bigl(\sum_{n=1}^{\infty} \frac{1}{n^2}\Bigr)^{1/2} < \infty,
$$
对 $b_n$ 同理。因此 $\sum (|a_n| + 2|b_n|)/n$ 收敛，由Weierstrass M-判别法，右端级数在 $[-\pi,\pi]$ 上 **绝对且一致收敛**。$\square$

---

**注记 16.1（积分改善收敛性的本质）** 对比 $f$ 和 $F$ 的Fourier系数：
- $f$ 的系数：$a_n, b_n$ 满足 $a_n,b_n \to 0$（推论16.1），但衰减速度取决于 $f$ 的光滑程度。
- $F$ 的系数：$A_n = -b_n/n$, $B_n = a_n/n$，额外多出因子 $1/n$，因此衰减速度提升了一阶。——即使 $f$ 仅分段连续（系数 $O(1/n)$），$F$ 的系数也达到 $O(1/n^2)$，足以保证绝对一致收敛。

这揭示了Fourier分析中的一个重要原则：**积分改善收敛性，求导削弱收敛性**。

---

**推论 16.2（一般区间上的逐项积分）** 设 $f$ 满足定理16.17的条件，$a,b \in [-\pi,\pi]$，则
$$
\int_a^b f(x)\,dx = \frac{a_0(b-a)}{2} + \sum_{n=1}^{\infty} \int_a^b \bigl(a_n\cos nx + b_n\sin nx\bigr)\,dx,
$$
且右端级数一致收敛。

**证明** 由定理16.17，$F(x)$ 的级数在 $[-\pi,\pi]$ 上一致收敛，代入 $F(b)-F(a)$ 即得。$\square$

---

**例题 16.12（从 $x$ 的Fourier级数求 $x^2$ 的展开）** 利用逐项积分定理，由
$$
x \sim 2\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n}\sin nx,\quad x\in(-\pi,\pi)
$$
推导 $x^2$ 在 $[-\pi,\pi]$ 上的Fourier级数。

**解** 取 $f(x)=x$。在 $[-\pi,\pi]$ 上 $f$ 为奇函数，$a_0=0$，$a_n=0$，$b_n = \dfrac{2(-1)^{n+1}}{n}$。由定理16.17（此时 $a_0=0$，$f(t)-a_0/2 = f(t)$）：
$$
\int_0^x t\,dt = \frac{x^2}{2} = \sum_{n=1}^{\infty} \frac{b_n(1-\cos nx)}{n}
= \sum_{n=1}^{\infty} \frac{2(-1)^{n+1}}{n} \cdot \frac{1-\cos nx}{n}
= 2\sum_{n=1}^{\infty} \frac{(-1)^{n+1}(1-\cos nx)}{n^2}.
$$

因此
$$
x^2 = 4\sum_{n=1}^{\infty} \frac{(-1)^{n+1}(1-\cos nx)}{n^2}.
$$

分离常数项与余弦项：
$$
\begin{aligned}
x^2 &= 4\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n^2} - 4\sum_{n=1}^{\infty} \frac{(-1)^{n+1}\cos nx}{n^2} \\
&= 4\cdot\frac{\pi^2}{12} - 4\sum_{n=1}^{\infty} \frac{(-1)^{n+1}\cos nx}{n^2} \qquad \big(\sum_{n=1}^{\infty}\frac{(-1)^{n+1}}{n^2}=\frac{\pi^2}{12}\big) \\
&= \frac{\pi^2}{3} + 4\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}\cos nx.
\end{aligned}
$$

与§1例题16.2直接计算的结果完全一致。验证完成。$\square$

---

**注记 16.2** 本例完美展示了逐项积分定理的价值：通过积分一个已知函数的Fourier级数，可以间接获得新函数的展开，避免重新计算复杂的Fourier系数积分。这种 **间接展开法** 是Fourier分析中的重要技巧。

---

### 3.2 Fourier级数的逐项求导定理

与逐项积分相比，逐项求导的条件更为严格。直观上，求导会放大高频分量（导数使Fourier系数乘以 $n$），因此对函数的光滑性要求更高。

#### 定理 16.18（Fourier级数的逐项求导定理）

设 $f$ 是以 $2\pi$ 为周期的函数，满足：
1. $f$ 在 $\mathbb{R}$ 上 **连续**；
2. $f(-\pi) = f(\pi)$（由周期性自动满足，但此条件在分部积分中起关键作用）；
3. 导函数 $f'$ 在 $[-\pi,\pi]$ 上 **分段连续**（即至多有限个第一类间断点）。

记 $f$ 的Fourier级数为
$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} \bigl(a_n\cos nx + b_n\sin nx\bigr),
$$
则 $f'$ 的Fourier级数可由 $f$ 的Fourier级数 **逐项求导** 得到：
$$
\boxed{\;
f'(x) \sim \sum_{n=1}^{\infty} \bigl(-n a_n\sin nx + n b_n\cos nx\bigr),
\;}
$$
且右端级数在 $f'$ 的连续点处收敛到 $f'(x)$，在 $f'$ 的跳跃间断点处收敛到左右极限的算术平均值。

---

**证明** 设 $f'$ 的Fourier系数为 $a_n', b_n'$。利用分部积分：

对 $a_n'$（$n \ge 1$）：
$$
\begin{aligned}
a_n' &= \frac{1}{\pi}\int_{-\pi}^{\pi} f'(x)\cos nx\,dx \\
&= \frac{1}{\pi}\Big( \big[f(x)\cos nx\big]_{-\pi}^{\pi} + n\int_{-\pi}^{\pi} f(x)\sin nx\,dx \Big) \\
&= \frac{n}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx \qquad (\text{因为 }f(-\pi)=f(\pi)\text{ 且 }\cos n\pi = \cos(-n\pi)) \\
&= n b_n.
\end{aligned}
$$

对 $b_n'$（$n \ge 1$）：
$$
\begin{aligned}
b_n' &= \frac{1}{\pi}\int_{-\pi}^{\pi} f'(x)\sin nx\,dx \\
&= \frac{1}{\pi}\Big( \big[f(x)\sin nx\big]_{-\pi}^{\pi} - n\int_{-\pi}^{\pi} f(x)\cos nx\,dx \Big) \\
&= -\frac{n}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx \qquad (\text{因为 }\sin n\pi = \sin(-n\pi)=0) \\
&= -n a_n.
\end{aligned}
$$

因此 $f'$ 的Fourier级数为
$$
f'(x) \sim \frac{a_0'}{2} + \sum_{n=1}^{\infty} (a_n'\cos nx + b_n'\sin nx)
      = \frac{a_0'}{2} + \sum_{n=1}^{\infty} (n b_n\cos nx - n a_n\sin nx).
$$

而形式上对 $f$ 的Fourier级数逐项求导：
$$
\frac{d}{dx}\Bigl[ \frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n\cos nx + b_n\sin nx) \Bigr]
= \sum_{n=1}^{\infty} (-n a_n\sin nx + n b_n\cos nx),
$$
两者完全一致。

由条件 $f'$ 分段连续，其Fourier级数满足Dirichlet-Jordan收敛定理（定理16.13），故在 $f'$ 的连续点处收敛到 $f'(x)$，在跳跃间断点处收敛到左右极限的算术平均值。$\square$

---

**注记 16.3（逐项求导定理与定理16.11的关系）** 定理16.11（Fourier系数与导数系数的关系）已经给出了 $a_n' = n b_n$, $b_n' = -n a_n$，但该定理要求 $f$ 连续可微。定理16.18将条件放宽为 $f$ 连续且 $f'$ 分段连续，同时回答了"级数是否收敛"以及"收敛到什么值"的问题。

---

#### 条件 $f(-\pi) = f(\pi)$ 的必要性

回顾证明中 $a_n'$ 的计算：
$$
a_n' = \frac{1}{\pi}\Big( \big[f(x)\cos nx\big]_{-\pi}^{\pi} + n\int_{-\pi}^{\pi} f(x)\sin nx\,dx \Big).
$$

边界项的具体值为：
$$
\big[f(x)\cos nx\big]_{-\pi}^{\pi} = f(\pi)\cos n\pi - f(-\pi)\cos(-n\pi) = \big(f(\pi) - f(-\pi)\big)\cos n\pi.
$$

当 $f(-\pi) = f(\pi)$ 时，此项为零，得到 $a_n' = n b_n$。当 $f(-\pi) \neq f(\pi)$ 时，边界项非零，$a_n' \neq n b_n$，因此逐项求导不成立。

**实际后果**：以 $f(x) = x$ 在 $[-\pi,\pi]$ 上为例：
- $f(-\pi) = -\pi \neq \pi = f(\pi)$，故 $f$ 的周期延拓在端点 $x = \pm\pi$ 处有跳跃间断点。
- 若形式地对 $x \sim 2\sum (-1)^{n+1}\sin(nx)/n$ 逐项求导：
  $$
  \frac{d}{dx}\Bigl[ 2\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n}\sin nx \Bigr] = 2\sum_{n=1}^{\infty} (-1)^{n+1}\cos nx,
  $$
  右端级数通项 $\cos nx$ 不趋于零，级数 **处处发散**（对任何 $x$，通项不趋于零）。
- 但 $f'(x) = 1$ 的Fourier级数是合法的正弦级数（不含常数项时 $1$ 展开为 $1 \sim \frac{1}{2} + \frac{2}{\pi}\sum_{k=1}^{\infty}\frac{\sin(2k-1)x}{2k-1}$），与逐项求导的结果完全不同。

这个例子清楚地表明：**$f(-\pi) = f(\pi)$ 是逐项求导定理不可或缺的条件**。

---

**例题 16.13** 考察 $f(x) = x^2$ 在 $[-\pi,\pi]$ 上的逐项求导。

**解** $f(x) = x^2$ 在 $[-\pi,\pi]$ 上连续，$f(-\pi) = \pi^2 = f(\pi)$，且 $f'(x) = 2x$ 在 $[-\pi,\pi]$ 上连续（从而分段连续），故满足定理16.18的全部条件。

$f$ 的Fourier级数为（§1例题16.2）：
$$
x^2 = \frac{\pi^2}{3} + 4\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}\cos nx,\qquad x\in[-\pi,\pi].
$$

逐项求导：
$$
\frac{d}{dx}\Bigl[ \frac{\pi^2}{3} + 4\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}\cos nx \Bigr]
= 4\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}\cdot(-n\sin nx)
= 4\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n}\sin nx.
$$

另一方面，直接计算 $f'(x) = 2x$ 的Fourier系数（$2x$ 是奇函数，$a_n=0$）：
$$
b_n = \frac{1}{\pi}\int_{-\pi}^{\pi} 2x\sin nx\,dx
     = \frac{4}{\pi}\int_0^{\pi} x\sin nx\,dx
     = \frac{4}{\pi}\cdot\frac{-\pi(-1)^n}{n}
     = \frac{4(-1)^{n+1}}{n}.
$$

因此 $2x \sim 4\sum_{n=1}^{\infty} \dfrac{(-1)^{n+1}}{n}\sin nx$，与逐项求导结果完全一致。$\square$

---

**注记 16.4（逐项积分与逐项求导的条件对比）**

| 性质 | 逐项积分（定理16.17） | 逐项求导（定理16.18） |
|------|---------------------|---------------------|
| $f$ 的连续性 | 分段连续即可 | 必须连续（周期延拓后） |
| $f(-\pi)=f(\pi)$ | 自动满足（积分后） | 必须满足 |
| 导数要求 | 无 | $f'$ 分段连续 |
| 结果级数收敛性 | 一致收敛 | 逐点收敛 |
| 系数变化 | $\times 1/n$（改善） | $\times n$（恶化） |

这一对比清晰展示了Fourier分析中"积分改善收敛性、求导削弱收敛性"的普遍规律。

---

### 3.3 Fourier系数的唯一性定理

逐项积分定理和逐项求导定理回答了"如何对Fourier级数进行运算"的问题。而唯一性定理则回答了更基本的问题：**Fourier级数是否唯一决定函数？**

#### 定理 16.19（Fourier系数的唯一性定理）

设 $f, g$ 在 $[-\pi,\pi]$ 上连续，且它们的Fourier系数全部相等（即 $a_n(f) = a_n(g)$，$b_n(f) = b_n(g)$ 对一切 $n$ 成立），则
$$
\boxed{\; f \equiv g \quad \text{在 } [-\pi,\pi] \text{ 上恒等}. \;}
$$

换言之，连续函数由它的Fourier系数 **唯一确定**。

---

**证明** 令 $h(x) = f(x) - g(x)$。由于 $f,g$ 的Fourier系数对应相等，$h$ 的所有Fourier系数为零：
$$
a_n(h) = \frac{1}{\pi}\int_{-\pi}^{\pi} h(x)\cos nx\,dx = 0,\quad n = 0,1,2,\dots,
$$
$$
b_n(h) = \frac{1}{\pi}\int_{-\pi}^{\pi} h(x)\sin nx\,dx = 0,\quad n = 1,2,3,\dots.
$$

下面给出两种证明方法。

**方法一（利用Bessel不等式）**：由Bessel不等式（定理16.15）：
$$
\frac{a_0^2}{2} + \sum_{n=1}^{\infty} (a_n^2 + b_n^2) \le \frac{1}{\pi}\int_{-\pi}^{\pi} h^2(x)\,dx.
$$

由于 $h$ 的所有Fourier系数为零，左端为零，故
$$
0 \le \frac{1}{\pi}\int_{-\pi}^{\pi} h^2(x)\,dx \le 0 \quad\Longrightarrow\quad \int_{-\pi}^{\pi} h^2(x)\,dx = 0.
$$

$h$ 在 $[-\pi,\pi]$ 上连续且 $h^2(x) \ge 0$。若存在 $x_0$ 使得 $h(x_0) \neq 0$，由连续性存在邻域 $(x_0-\delta, x_0+\delta)$ 使得 $|h(x)| \ge |h(x_0)|/2 > 0$，从而 $\int_{-\pi}^{\pi} h^2(x)\,dx > 0$，矛盾。因此 $h \equiv 0$，即 $f \equiv g$。

**方法二（利用Fejér核与Cesàro平均）**：考虑 $h$ 的Fejér（算术平均）部分和：
$$
\sigma_N(x) = \frac{S_0(x) + S_1(x) + \cdots + S_{N-1}(x)}{N},
$$
其中 $S_n(x)$ 是 $h$ 的Fourier级数的 $n$ 阶部分和。Fejér定理指出：对连续函数 $h$，$\sigma_N(x)$ 在 $[-\pi,\pi]$ 上 **一致收敛** 到 $h(x)$。

但 $\sigma_N(x)$ 完全由 $h$ 的Fourier系数线性组合而成。由于 $h$ 的所有Fourier系数为零，必有 $\sigma_N(x) \equiv 0$ 对一切 $N$ 成立。因此一致极限 $h \equiv 0$。$\square$

---

**注记 16.5** 方法一（Bessel不等式法）简洁直接，但依赖于Parseval等式的完备性（即Bessel不等式中等号成立）；方法二（Fejér核法）不依赖 $L^2$ 理论，是更初等的证明。两种证明各有优势。

---

**推论 16.3** 设 $f$ 在 $[-\pi,\pi]$ 上连续。若 $f$ 的所有Fourier系数为零，则 $f \equiv 0$。

**证明** 在定理16.19中取 $g \equiv 0$ 即得。$\square$

---

**推论 16.4（正交系的完备性）** 设 $f$ 在 $[-\pi,\pi]$ 上连续。若 $f$ 与三角函数系中每个函数正交（即 $\int_{-\pi}^{\pi} f(x)\cos nx\,dx = \int_{-\pi}^{\pi} f(x)\sin nx\,dx = 0$ 对一切 $n$ 成立），则 $f \equiv 0$。

换言之，在连续函数空间中，三角函数系是 **完备正交系** ——除了零函数外，没有其他连续函数与所有三角函数正交。

---

**例题 16.14** 设 $f$ 在 $[-\pi,\pi]$ 上分段连续，且 $\displaystyle\int_{-\pi}^{\pi} f(t)\,dt = 0$。定义
$$
\Phi(x) = \int_{-\pi}^{x} f(t)\,dt.
$$
利用逐项积分定理和唯一性定理，证明 $\Phi$ 的Fourier级数可通过逐项积分得到，且一致收敛到 $\Phi(x)$。

**解** 由 $\int_{-\pi}^{\pi} f(t)\,dt = 0$ 知 $a_0 = 0$，故 $f$ 的Fourier级数不含常数项：
$$
f(x) \sim \sum_{n=1}^{\infty} (a_n\cos nx + b_n\sin nx).
$$

定义
$$
\Psi(x) = \int_{-\pi}^{x} f(t)\,dt.
$$
则 $\Psi$ 连续，$\Psi(-\pi)=0$，且 $\Psi(\pi)=\int_{-\pi}^{\pi} f(t)\,dt = 0$，故 $\Psi(-\pi) = \Psi(\pi)$。

对 $\Psi$ 应用逐项积分定理（取 $F(x) = \int_0^x [f(t)-a_0/2]\,dt$ 的形式稍作调整）。具体地，由于 $a_0 = 0$，定理16.17给出：
$$
\int_0^x f(t)\,dt = \sum_{n=1}^{\infty} \frac{a_n\sin nx + b_n(1 - \cos nx)}{n},
$$
且右端级数一致收敛。

而
$$
\Psi(x) = \int_{-\pi}^{x} f(t)\,dt = \int_0^{x} f(t)\,dt - \int_0^{-\pi} f(t)\,dt.
$$

记 $C = \int_0^{-\pi} f(t)\,dt$，则
$$
\Psi(x) = \sum_{n=1}^{\infty} \frac{a_n\sin nx + b_n(1 - \cos nx)}{n} - C.
$$

由 $\Psi(-\pi)=0$ 可确定常数 $C$，但更重要的是：由逐项积分定理，右端级数 **一致收敛** 到某个连续函数。由唯一性定理（定理16.19），这个一致收敛的级数必然等于 $\Psi(x)$ 本身。

因此 $\Phi$ 的Fourier级数可以通过对 $f$ 的Fourier级数逐项积分得到，且该级数在 $[-\pi,\pi]$ 上一致收敛。$\square$

---

### 3.4 综合例题

---

**例题 16.15** 设 $f(x)$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上分段光滑。利用逐项积分定理和逐项求导定理，判断以下命题的真伪：

**(1)** 若 $f$ 的Fourier级数在 $[-\pi,\pi]$ 上一致收敛，则 $f$ 连续。  
**(2)** 若 $f$ 连续且 $f(-\pi)=f(\pi)$，则 $f$ 的Fourier级数在 $[-\pi,\pi]$ 上一致收敛。  
**(3)** 若 $f$ 的Fourier级数逐项求导后一致收敛，则 $f'$ 连续。

**解**

**(1) 真。** 若Fourier级数一致收敛，则和函数（即 $f$ 在连续点处的值，修正后得到整体函数）连续，因为一致收敛的连续函数项级数的和函数连续。因此 $f$ 的周期延拓必须在 $\mathbb{R}$ 上连续。

**(2) 伪。** 反例：$f(x) = \begin{cases} \pi - x, & x \in (0,\pi], \\ x + \pi, & x \in [-\pi,0] \end{cases}$。$f$ 连续且 $f(-\pi)=f(\pi)$，但其Fourier级数（见§1例题16.7）为
$$
f(x) \sim \frac{\pi}{2} + \frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\cos(2k-1)x}{(2k-1)^2},
$$
该级数在 $[-\pi,\pi]$ 上绝对一致收敛（系数 $O(1/k^2)$），与(2)的结论并不矛盾——实际上该反例支持一致收敛，而非否定。更好的反例：取 $f$ 连续且 $f(-\pi)=f(\pi)$ 但不可微，其Fourier系数约为 $O(1/n)$，此时级数不一定一致收敛。实际上，若 $f$ 连续且分段光滑，但端点处导数有跳跃（如折线函数），系数衰减为 $O(1/n^2)$，级数一致收敛。因此(2)的"若则"方向需要进一步条件（如 $f$ 绝对连续或 $f'$ 平方可积）。命题本身过于宽泛。

**(3) 真。** 设 $\sum_{n=1}^{\infty} (-n a_n\sin nx + n b_n\cos nx)$ 一致收敛到某函数 $g(x)$。由于一致收敛的级数可以逐项积分，对 $g$ 的级数逐项积分再结合定理16.17，可得 $f$ 处处可微且 $f' = g$。而一致收敛级数的和函数连续，故 $f'$ 连续。$\square$

---

**例题 16.16** 利用逐项积分定理，求 $\displaystyle\sum_{n=1}^{\infty} \frac{\sin nx}{n}$ 的和函数。

**解** 考虑 $f(x) \sim \sum_{n=1}^{\infty} \frac{\sin nx}{n}$。这是奇函数形式的Fourier级数，对应某个奇函数 $f$。

由逐项积分定理的逆推：若 $f(x) = \sum_{n=1}^{\infty} b_n \sin nx$，则积分后得到：
$$
\int_0^x f(t)\,dt = \sum_{n=1}^{\infty} \frac{b_n(1 - \cos nx)}{n}.
$$

令 $b_n = 1/n$，则
$$
\int_0^x f(t)\,dt = \sum_{n=1}^{\infty} \frac{1 - \cos nx}{n^2}.
$$

右端级数在 $[-\pi,\pi]$ 上一致收敛。我们来识别左端的函数。计算右端的和函数：
$$
\sum_{n=1}^{\infty} \frac{1 - \cos nx}{n^2} = \sum_{n=1}^{\infty} \frac{1}{n^2} - \sum_{n=1}^{\infty} \frac{\cos nx}{n^2} = \frac{\pi^2}{6} - \sum_{n=1}^{\infty} \frac{\cos nx}{n^2}.
$$

而 $\sum_{n=1}^{\infty} \frac{\cos nx}{n^2}$ 是 $x^2$ 的Fourier级数的变形。由§1例题16.2：
$$
x^2 = \frac{\pi^2}{3} + 4\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}\cos nx,
$$
这里 $(-1)^n$ 表示不同区间的系数。对于 $[0,2\pi]$ 区间上的展开，有
$$
\frac{x^2}{2} - \pi x + \frac{\pi^2}{3} = \sum_{n=1}^{\infty} \frac{\cos nx}{n^2},\qquad x\in[0,2\pi].
$$

更直接地，在 $[-\pi,\pi]$ 上：
$$
\sum_{n=1}^{\infty} \frac{\cos nx}{n^2} = \frac{x^2}{4} - \frac{\pi x}{2} + \frac{\pi^2}{6},\qquad x\in[0,2\pi].
$$

代入得
$$
\int_0^x f(t)\,dt = \frac{\pi^2}{6} - \left(\frac{x^2}{4} - \frac{\pi x}{2} + \frac{\pi^2}{6}\right) = \frac{\pi x}{2} - \frac{x^2}{4}.
$$

因此
$$
f(x) = \frac{d}{dx}\left(\frac{\pi x}{2} - \frac{x^2}{4}\right) = \frac{\pi}{2} - \frac{x}{2},\qquad x\in(0,2\pi).
$$

即
$$
\boxed{\; \sum_{n=1}^{\infty} \frac{\sin nx}{n} = \frac{\pi - x}{2},\qquad x\in(0,2\pi). \;}
$$

在 $[-\pi,\pi]$ 上，该级数对应 $f(x) = -\frac{x}{2}$（调整常数项），这正是$f(x)=x$的Fourier级数减去因子。$\square$

**注**：本例展示了利用逐项积分定理 **求和函数** 的技巧——先积分消除分母中的 $n$，得到已知级数，再求导还原。

---

### 习题

1. 设 $f(x) = \begin{cases} 0, & -\pi \le x < 0, \\ 1, & 0 \le x \le \pi \end{cases}$，以 $2\pi$ 为周期延拓。
   (1) 求 $f$ 的Fourier级数；
   (2) 利用逐项积分定理，求 $\int_0^x f(t)\,dt$ 的Fourier级数展开；
   (3) 验证结果与直接计算Fourier系数的一致性。

2. 设 $f$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上分段光滑。证明：
   $$
   \int_{-\pi}^{\pi} [f(x)]^2\,dx = \frac{\pi a_0^2}{2} + \pi\sum_{n=1}^{\infty} (a_n^2 + b_n^2) - 2\pi\sum_{n=1}^{\infty} \frac{a_n b_n}{n}.
   $$
   （提示：考虑 $f$ 的Fourier级数与 $f'$ 的Fourier级数的Parseval等式关系。）

3. 利用逐项积分定理，从 $f(x) = x$ 的Fourier级数出发，求 $f(x) = x^3$ 在 $[-\pi,\pi]$ 上的Fourier级数。（提示：可能需要反复应用逐项积分定理。）

4. 判断以下函数是否满足逐项求导定理的条件：
   (a) $f(x) = |x|$ 在 $[-\pi,\pi]$ 上；
   (b) $f(x) = \sin x$ 以 $2\pi$ 为周期；
   (c) $f(x) = \begin{cases} 1, & -\pi < x < 0, \\ -1, & 0 < x < \pi \end{cases}$ 以 $2\pi$ 为周期延拓；
   (d) $f(x) = \pi^2 - x^2$ 在 $[-\pi,\pi]$ 上。

5. 设 $f$ 在 $[-\pi,\pi]$ 上连续，$f(-\pi)=f(\pi)$，且 $f'$ 分段连续。证明 $f$ 的Fourier级数在 $[-\pi,\pi]$ 上一致收敛。（提示：利用逐项求导定理和Weierstrass M-判别法，结合 $a_n',b_n'$ 的衰减性质。）

6. 利用唯一性定理证明：若 $f$ 在 $[-\pi,\pi]$ 上连续，且对一切 $n$ 有 $a_{2n}=b_{2n}=0$，则 $f(x+\pi) = -f(x)$。

7. 设 $f$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上光滑。利用逐项积分定理和逐项求导定理，证明：
   $$
   \sum_{n=1}^{\infty} n^2 (a_n^2 + b_n^2) = \frac{1}{\pi}\int_{-\pi}^{\pi} [f'(x)]^2\,dx.
   $$

8. **（选做）** 利用定理16.17和16.18，证明：
   $$
   \sum_{n=1}^{\infty} \frac{\cos nx}{n^2} = \frac{x^2}{4} - \frac{\pi x}{2} + \frac{\pi^2}{6},\qquad x\in[0,2\pi].
   $$
   并用此结果求 $\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^n}{n^2}$ 和 $\displaystyle\sum_{n=1}^{\infty} \frac{1}{n^2}$。

---

### 本章小结

本节建立了Fourier级数的三个核心分析性质及其相互关系：

| 主题 | 核心内容 | 对应定理 |
|------|----------|----------|
| 逐项积分定理 | 积分改善收敛性（系数 $\times 1/n$），结果一致收敛 | 定理16.17 |
| 逐项求导定理 | 求导削弱收敛性（系数 $\times n$），需 $f$ 连续且 $f(-\pi)=f(\pi)$ | 定理16.18 |
| 唯一性定理 | 连续函数由Fourier系数唯一确定 | 定理16.19 |
| 积分与求导的对比 | 积分无条件、求导有条件；积分改善收敛、求导恶化收敛 | 注记16.4 |

**核心思想**：Fourier级数的分析性质（逐项积分、逐项求导、唯一性）构成了Fourier分析理论的基础工具。它们不仅揭示了Fourier系数与函数光滑性之间的深刻联系，也为间接展开、求和函数等实际应用提供了严格的数学依据。积分"平滑"函数从而改善收敛性，而求导"尖锐"函数从而削弱收敛性——这一对偶关系贯穿整个Fourier分析。
