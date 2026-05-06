# 第十六章 Fourier级数

---

## §1 Fourier级数的基本概念与展开

---

### 1.1 三角函数系的正交性

**定义 16.1（三角函数系）** 函数集合

$$
\{1,\; \cos x,\; \cos 2x,\; \dots,\; \cos nx,\; \dots;\; \sin x,\; \sin 2x,\; \dots,\; \sin nx,\; \dots\}
$$

称为 **三角函数系**。这些函数都是以 $2\pi$ 为周期的周期函数。

**定义 16.2（函数的正交性）** 设 $f(x), g(x)$ 在区间 $[a, b]$ 上可积。若

$$
\int_a^b f(x)g(x)\,dx = 0,
$$

则称 $f$ 与 $g$ 在 $[a, b]$ 上 **正交**。

---

**定理 16.1（三角函数系的正交性）** 三角函数系在区间 $[-\pi, \pi]$ 上具有正交性，即对任意非负整数 $m, n$，以下关系成立：

**(1)** 余弦函数与余弦函数：

$$
\int_{-\pi}^{\pi} \cos mx \cos nx \, dx =
\begin{cases}
0, & m \neq n, \\[4pt]
\pi, & m = n \neq 0, \\[4pt]
2\pi, & m = n = 0.
\end{cases}
$$

**(2)** 正弦函数与正弦函数：

$$
\int_{-\pi}^{\pi} \sin mx \sin nx \, dx =
\begin{cases}
0, & m \neq n, \\[4pt]
\pi, & m = n \neq 0.
\end{cases}
$$

**(3)** 余弦函数与正弦函数：

$$
\int_{-\pi}^{\pi} \cos mx \sin nx \, dx = 0, \quad \forall\, m, n.
$$

**证明** 利用三角恒等式将乘积化为和差。

**(1)** 当 $m \neq n$ 时：

$$
\begin{aligned}
\int_{-\pi}^{\pi} \cos mx \cos nx \, dx
&= \frac{1}{2}\int_{-\pi}^{\pi} \bigl[\cos(m+n)x + \cos(m-n)x\bigr]\,dx \\
&= \frac{1}{2}\left[ \frac{\sin(m+n)x}{m+n} + \frac{\sin(m-n)x}{m-n} \right]_{-\pi}^{\pi} = 0.
\end{aligned}
$$

当 $m = n \neq 0$ 时：

$$
\begin{aligned}
\int_{-\pi}^{\pi} \cos^2 nx \, dx
&= \frac{1}{2}\int_{-\pi}^{\pi} (1 + \cos 2nx)\,dx \\
&= \frac{1}{2}\left[ x + \frac{\sin 2nx}{2n} \right]_{-\pi}^{\pi} = \pi.
\end{aligned}
$$

当 $m = n = 0$ 时：$\displaystyle\int_{-\pi}^{\pi} 1^2\,dx = 2\pi$。

**(2)** 当 $m \neq n$ 时：

$$
\begin{aligned}
\int_{-\pi}^{\pi} \sin mx \sin nx \, dx
&= -\frac{1}{2}\int_{-\pi}^{\pi} \bigl[\cos(m+n)x - \cos(m-n)x\bigr]\,dx \\
&= -\frac{1}{2}\left[ \frac{\sin(m+n)x}{m+n} - \frac{\sin(m-n)x}{m-n} \right]_{-\pi}^{\pi} = 0.
\end{aligned}
$$

当 $m = n \neq 0$ 时：

$$
\int_{-\pi}^{\pi} \sin^2 nx \, dx = \frac{1}{2}\int_{-\pi}^{\pi} (1 - \cos 2nx)\,dx = \pi.
$$

**(3)** 对任意 $m, n$：

$$
\begin{aligned}
\int_{-\pi}^{\pi} \cos mx \sin nx \, dx
&= \frac{1}{2}\int_{-\pi}^{\pi} \bigl[\sin(m+n)x + \sin(n-m)x\bigr]\,dx = 0.
\end{aligned}
$$

综上，三角函数系在 $[-\pi, \pi]$ 上两两正交。$\square$

---

### 1.2 周期函数的 Fourier 级数展开

**问题**：设 $f(x)$ 是以 $2\pi$ 为周期的函数（或定义在 $[-\pi, \pi]$ 上的函数），能否将其表示为三角函数系的线性组合？

**定义 16.3（Fourier 级数）** 设 $f(x)$ 在 $[-\pi, \pi]$ 上可积，称形式级数

$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} \bigl(a_n \cos nx + b_n \sin nx\bigr)
$$

为 $f$ 的 **Fourier 级数**，其中系数 $a_0, a_n, b_n$ 由以下 **Euler-Fourier 公式** 给出：

$$
\boxed{\;
a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \cos nx \, dx, \quad n = 0, 1, 2, \dots
\;}
$$

$$
\boxed{\;
b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \sin nx \, dx, \quad n = 1, 2, 3, \dots
\;}
$$

**注**：
1. $a_0$ 公式包含在 $a_n$ 中（取 $n=0$，$\cos 0 = 1$），即 $a_0 = \dfrac{1}{\pi} \displaystyle\int_{-\pi}^{\pi} f(x)\,dx$。
2. 记号 "$\sim$" 表示右端是 $f$ 的形式 Fourier 级数，其收敛性需要进一步讨论。
3. 常数项写作 $\dfrac{a_0}{2}$ 是为了使 $a_0$ 的公式与 $a_n$ 统一。

**公式推导（系数确定）** 假设 $f$ 可以表示为一致收敛的三角级数：

$$
f(x) = \frac{a_0}{2} + \sum_{k=1}^{\infty} \bigl(a_k \cos kx + b_k \sin kx\bigr).
$$

利用正交性求 $a_n$：两边同乘 $\cos nx$，在 $[-\pi, \pi]$ 上积分：

$$
\int_{-\pi}^{\pi} f(x) \cos nx \, dx
= \frac{a_0}{2} \int_{-\pi}^{\pi} \cos nx \, dx
+ \sum_{k=1}^{\infty} \Bigl(a_k \int_{-\pi}^{\pi} \cos kx \cos nx \, dx + b_k \int_{-\pi}^{\pi} \sin kx \cos nx \, dx\Bigr).
$$

由正交性：
- 当 $n = 0$ 时，右端仅 $\int_{-\pi}^{\pi} 1^2\,dx = 2\pi$ 非零，得 $a_0 = \dfrac{1}{\pi} \int_{-\pi}^{\pi} f(x)\,dx$。
- 当 $n \ge 1$ 时，右端仅 $k = n$ 的余弦项贡献 $\pi a_n$，得 $a_n = \dfrac{1}{\pi} \int_{-\pi}^{\pi} f(x)\cos nx \, dx$。

同理，两边同乘 $\sin nx$ 后积分，得 $b_n = \dfrac{1}{\pi} \int_{-\pi}^{\pi} f(x)\sin nx \, dx$。

---

**例题 16.1** 求函数 $f(x) = x$ 在 $[-\pi, \pi]$ 上的 Fourier 级数。

**解** $f(x) = x$ 是奇函数，$a_n = 0$（见 1.3 节奇偶性简化）。

计算正弦系数：

$$
\begin{aligned}
b_n &= \frac{1}{\pi} \int_{-\pi}^{\pi} x \sin nx \, dx
= \frac{2}{\pi} \int_{0}^{\pi} x \sin nx \, dx \quad (\text{奇函数} \times \text{奇函数} = \text{偶函数}) \\
&= \frac{2}{\pi} \left( \left[ -\frac{x \cos nx}{n} \right]_{0}^{\pi} + \frac{1}{n} \int_{0}^{\pi} \cos nx \, dx \right) \\
&= \frac{2}{\pi} \left( -\frac{\pi \cos n\pi}{n} + \frac{1}{n} \left[ \frac{\sin nx}{n} \right]_{0}^{\pi} \right) \\
&= \frac{2}{\pi} \cdot \frac{-\pi (-1)^n}{n} = \frac{2(-1)^{n+1}}{n}.
\end{aligned}
$$

因此

$$
x \sim 2 \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n} \sin nx, \quad x \in (-\pi, \pi).
$$

---

### 1.3 奇偶函数的 Fourier 级数简化

若 $f$ 具有奇偶性，Fourier 系数可大幅简化。

**定理 16.2（奇函数的 Fourier 级数）** 若 $f(x)$ 在 $[-\pi, \pi]$ 上是奇函数，则

$$
a_n = 0 \quad (n = 0, 1, 2, \dots), \qquad
b_n = \frac{2}{\pi} \int_{0}^{\pi} f(x) \sin nx \, dx.
$$

此时 Fourier 级数退化为 **正弦级数**：

$$
f(x) \sim \sum_{n=1}^{\infty} b_n \sin nx.
$$

**证明** 奇函数 $f(x)\cos nx$ 是奇函数（奇 $\times$ 偶 $=$ 奇），在对称区间上积分为零，故 $a_n = 0$。而 $f(x)\sin nx$ 是偶函数（奇 $\times$ 奇 $=$ 偶），故

$$
b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \sin nx \, dx = \frac{2}{\pi} \int_{0}^{\pi} f(x) \sin nx \, dx.
$$

$\square$

**定理 16.3（偶函数的 Fourier 级数）** 若 $f(x)$ 在 $[-\pi, \pi]$ 上是偶函数，则

$$
b_n = 0 \quad (n = 1, 2, \dots), \qquad
a_n = \frac{2}{\pi} \int_{0}^{\pi} f(x) \cos nx \, dx.
$$

此时 Fourier 级数退化为 **余弦级数**：

$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos nx.
$$

**证明** 偶函数 $f(x)\sin nx$ 是奇函数（偶 $\times$ 奇 $=$ 奇），积分为零，故 $b_n = 0$。而 $f(x)\cos nx$ 是偶函数（偶 $\times$ 偶 $=$ 偶），故

$$
a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \cos nx \, dx = \frac{2}{\pi} \int_{0}^{\pi} f(x) \cos nx \, dx.
$$

$\square$

---

**例题 16.2** 求 $f(x) = x^2$ 在 $[-\pi, \pi]$ 上的 Fourier 级数。

**解** $f(x) = x^2$ 是偶函数，故 $b_n = 0$。

计算 $a_0$：

$$
a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi} x^2 \, dx = \frac{2}{\pi} \int_{0}^{\pi} x^2 \, dx = \frac{2}{\pi} \cdot \frac{\pi^3}{3} = \frac{2\pi^2}{3}.
$$

计算 $a_n \; (n \ge 1)$：

$$
\begin{aligned}
a_n &= \frac{2}{\pi} \int_{0}^{\pi} x^2 \cos nx \, dx.
\end{aligned}
$$

使用分部积分：

$$
\begin{aligned}
\int_{0}^{\pi} x^2 \cos nx \, dx
&= \left[ \frac{x^2 \sin nx}{n} \right]_{0}^{\pi} - \frac{2}{n} \int_{0}^{\pi} x \sin nx \, dx \\
&= 0 - \frac{2}{n} \int_{0}^{\pi} x \sin nx \, dx.
\end{aligned}
$$

再分部积分：

$$
\begin{aligned}
\int_{0}^{\pi} x \sin nx \, dx
&= \left[ -\frac{x \cos nx}{n} \right]_{0}^{\pi} + \frac{1}{n} \int_{0}^{\pi} \cos nx \, dx \\
&= -\frac{\pi \cos n\pi}{n} + \frac{1}{n} \left[ \frac{\sin nx}{n} \right]_{0}^{\pi} \\
&= -\frac{\pi (-1)^n}{n}.
\end{aligned}
$$

因此

$$
\int_{0}^{\pi} x^2 \cos nx \, dx = -\frac{2}{n} \cdot \left( -\frac{\pi (-1)^n}{n} \right) = \frac{2\pi (-1)^n}{n^2},
$$

$$
a_n = \frac{2}{\pi} \cdot \frac{2\pi (-1)^n}{n^2} = \frac{4(-1)^n}{n^2}.
$$

故

$$
\boxed{\; x^2 = \frac{\pi^2}{3} + 4 \sum_{n=1}^{\infty} \frac{(-1)^n}{n^2} \cos nx, \quad x \in [-\pi, \pi]. \;}
$$

**注**：上式在 $x = \pi$ 处给出 $\pi^2 = \dfrac{\pi^2}{3} + 4 \displaystyle\sum_{n=1}^{\infty} \dfrac{(-1)^n}{n^2} (-1)^n$，整理可得

$$
\sum_{n=1}^{\infty} \frac{1}{n^2} = \frac{\pi^2}{6}.
$$

这是 Euler 著名结果的 Fourier 级数推导。

---

### 1.4 半幅展开（Half-Range Expansion）

实际问题中，函数往往只定义在半个区间 $[0, \pi]$ 上。我们可以通过 **奇延拓** 或 **偶延拓** 将其延拓到 $[-\pi, \pi]$ 上，然后展开成 Fourier 级数。

**定义 16.4（半幅展开）** 设 $f(x)$ 定义在 $[0, \pi]$ 上。

**(a) 余弦级数（偶延拓）**：将 $f$ 延拓为 $[-\pi, \pi]$ 上的偶函数 $\tilde{f}(x)$：

$$
\tilde{f}(x) = \begin{cases}
f(x), & x \in [0, \pi], \\
f(-x), & x \in [-\pi, 0).
\end{cases}
$$

则 $\tilde{f}$ 的 Fourier 级数退化为余弦级数，系数为

$$
\boxed{\;
a_n = \frac{2}{\pi} \int_{0}^{\pi} f(x) \cos nx \, dx, \quad n = 0, 1, 2, \dots
\;}
$$

$$
b_n = 0, \quad n = 1, 2, \dots
$$

在 $[0, \pi]$ 上，$\tilde{f}(x) = f(x)$，故

$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos nx, \quad x \in [0, \pi].
$$

**(b) 正弦级数（奇延拓）**：将 $f$ 延拓为 $[-\pi, \pi]$ 上的奇函数 $\tilde{f}(x)$：

$$
\tilde{f}(x) = \begin{cases}
f(x), & x \in [0, \pi], \\
-f(-x), & x \in (-\pi, 0), \\
0, & x = 0.
\end{cases}
$$

则 $\tilde{f}$ 的 Fourier 级数退化为正弦级数，系数为

$$
\boxed{\;
b_n = \frac{2}{\pi} \int_{0}^{\pi} f(x) \sin nx \, dx, \quad n = 1, 2, 3, \dots
\;}
$$

$$
a_n = 0, \quad n = 0, 1, 2, \dots
$$

在 $[0, \pi]$ 上，$\tilde{f}(x) = f(x)$，故

$$
f(x) \sim \sum_{n=1}^{\infty} b_n \sin nx, \quad x \in [0, \pi].
$$

---

**例题 16.3** 设 $f(x) = x$，$x \in [0, \pi]$，分别展开为余弦级数和正弦级数。

**解 (a) 余弦级数（偶延拓）**：

偶延拓后 $\tilde{f}(x) = |x|$ 在 $[-\pi, \pi]$ 上是偶函数，$b_n = 0$。

$$
a_0 = \frac{2}{\pi} \int_{0}^{\pi} x \, dx = \frac{2}{\pi} \cdot \frac{\pi^2}{2} = \pi.
$$

对 $n \ge 1$：

$$
\begin{aligned}
a_n &= \frac{2}{\pi} \int_{0}^{\pi} x \cos nx \, dx \\
&= \frac{2}{\pi} \left( \left[ \frac{x \sin nx}{n} \right]_{0}^{\pi} - \frac{1}{n} \int_{0}^{\pi} \sin nx \, dx \right) \\
&= \frac{2}{\pi} \left( 0 + \frac{1}{n} \left[ \frac{\cos nx}{n} \right]_{0}^{\pi} \right) \\
&= \frac{2}{\pi} \cdot \frac{\cos n\pi - 1}{n^2} = \frac{2\bigl[(-1)^n - 1\bigr]}{\pi n^2}.
\end{aligned}
$$

当 $n$ 为偶数时 $a_n = 0$；当 $n$ 为奇数时，记 $n = 2k-1$：

$$
a_{2k-1} = \frac{2}{\pi} \cdot \frac{-2}{(2k-1)^2} = -\frac{4}{\pi(2k-1)^2}.
$$

因此

$$
\boxed{\; x = \frac{\pi}{2} - \frac{4}{\pi} \sum_{k=1}^{\infty} \frac{\cos(2k-1)x}{(2k-1)^2}, \quad x \in [0, \pi]. \;}
$$

**(b) 正弦级数（奇延拓）**：

奇延拓后 $\tilde{f}(x) = x$ 在 $[-\pi, \pi]$ 上是奇函数，$a_n = 0$。由例题 16.1 结果：

$$
b_n = \frac{2}{\pi} \int_{0}^{\pi} x \sin nx \, dx = \frac{2(-1)^{n+1}}{n}.
$$

因此

$$
\boxed{\; x = 2 \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n} \sin nx, \quad x \in [0, \pi). \;}
$$

**注**：在端点 $x = \pi$ 处，正弦级数收敛到 $0$（奇延拓的间断点），而 $f(\pi) = \pi$，因此等式只在 $x \in [0, \pi)$ 成立。余弦级数在端点处无此跳跃现象。

---

### 1.5 一般区间上的 Fourier 展开

对于定义在长度为 $2\pi$ 的一般区间 $[a, a+2\pi]$ 上的函数，只需将积分区间平移。

**定理 16.4（一般区间的 Fourier 系数）** 设 $f(x)$ 在 $[a, a+2\pi]$ 上可积，则其 Fourier 级数为

$$
f(x) \sim \frac{a_0}{2} + \sum_{n=1}^{\infty} \bigl(a_n \cos nx + b_n \sin nx\bigr),
$$

其中

$$
\boxed{\;
a_n = \frac{1}{\pi} \int_{a}^{a+2\pi} f(x) \cos nx \, dx, \quad n = 0, 1, 2, \dots
\;}
$$

$$
\boxed{\;
b_n = \frac{1}{\pi} \int_{a}^{a+2\pi} f(x) \sin nx \, dx, \quad n = 1, 2, 3, \dots
\;}
$$

**证明** 正交性在任意长度为 $2\pi$ 的区间上成立（被积函数以 $2\pi$ 为周期）。系数公式的推导与标准区间完全相同。$\square$

---

**例题 16.4** 求 $f(x) = x$ 在 $[0, 2\pi]$ 上的 Fourier 级数。

**解** 取 $a = 0$，使用一般区间公式。

计算 $a_0$：

$$
a_0 = \frac{1}{\pi} \int_{0}^{2\pi} x \, dx = \frac{1}{\pi} \cdot \frac{(2\pi)^2}{2} = 2\pi.
$$

计算 $a_n \; (n \ge 1)$：

$$
\begin{aligned}
a_n &= \frac{1}{\pi} \int_{0}^{2\pi} x \cos nx \, dx \\
&= \frac{1}{\pi} \left( \left[ \frac{x \sin nx}{n} \right]_{0}^{2\pi} - \frac{1}{n} \int_{0}^{2\pi} \sin nx \, dx \right) \\
&= \frac{1}{\pi} \left( 0 + \frac{1}{n} \left[ \frac{\cos nx}{n} \right]_{0}^{2\pi} \right) \\
&= \frac{1}{\pi} \cdot \frac{1 - 1}{n^2} = 0.
\end{aligned}
$$

计算 $b_n$：

$$
\begin{aligned}
b_n &= \frac{1}{\pi} \int_{0}^{2\pi} x \sin nx \, dx \\
&= \frac{1}{\pi} \left( \left[ -\frac{x \cos nx}{n} \right]_{0}^{2\pi} + \frac{1}{n} \int_{0}^{2\pi} \cos nx \, dx \right) \\
&= \frac{1}{\pi} \left( -\frac{2\pi \cos(2n\pi)}{n} + \frac{1}{n} \left[ \frac{\sin nx}{n} \right]_{0}^{2\pi} \right) \\
&= \frac{1}{\pi} \left( -\frac{2\pi}{n} \right) = -\frac{2}{n}.
\end{aligned}
$$

因此

$$
\boxed{\; x = \pi - 2 \sum_{n=1}^{\infty} \frac{\sin nx}{n}, \quad x \in (0, 2\pi). \;}
$$

**对比**：与例 16.1 比较——同是 $f(x)=x$，在 $[-\pi,\pi]$ 上展开为正弦级数，在 $[0,2\pi]$ 上却出现了非零常数项 $\pi$。这是因为区间位置不同导致函数的平均分量不同。

---

### 1.6 周期函数的 Fourier 系数对称性

利用周期函数的对称性可以得到 Fourier 系数的特殊消零规律。

**定理 16.5（半周期平移的对称性）** 设 $f(x)$ 是以 $2\pi$ 为周期的可积函数，且满足

$$
f(x + \pi) = f(x), \quad \forall x \in [-\pi, \pi],
$$

则 $f$ 的 Fourier 级数中所有奇数阶系数为零，即

$$
a_{2n-1} = b_{2n-1} = 0, \quad n = 1, 2, 3, \dots
$$

换言之，Fourier 级数中只出现 $1, \cos 2x, \cos 4x, \dots, \sin 2x, \sin 4x, \dots$ 等偶次谐波。

**证明** 以 $a_{2n-1}$ 为例：

$$
\begin{aligned}
a_{2n-1} &= \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \cos(2n-1)x \, dx \\
&= \frac{1}{\pi} \int_{-\pi}^{0} f(x) \cos(2n-1)x \, dx + \frac{1}{\pi} \int_{0}^{\pi} f(x) \cos(2n-1)x \, dx.
\end{aligned}
$$

在第一项中令 $t = x + \pi$，则 $x = t - \pi$，$x \in [-\pi, 0]$ 对应 $t \in [0, \pi]$：

$$
\begin{aligned}
\int_{-\pi}^{0} f(x) \cos(2n-1)x \, dx
&= \int_{0}^{\pi} f(t - \pi) \cos\bigl[(2n-1)(t - \pi)\bigr] \, dt \\
&= \int_{0}^{\pi} f(t) \cos\bigl[(2n-1)t - (2n-1)\pi\bigr] \, dt \quad (\because f(t-\pi)=f(t)) \\
&= \int_{0}^{\pi} f(t) \bigl[ \cos(2n-1)t \cdot \underbrace{\cos(2n-1)\pi}_{= -1} + \sin(2n-1)t \cdot \underbrace{\sin(2n-1)\pi}_{= 0} \bigr] \, dt \\
&= -\int_{0}^{\pi} f(t) \cos(2n-1)t \, dt.
\end{aligned}
$$

因此

$$
a_{2n-1} = \frac{1}{\pi} \left( -\int_{0}^{\pi} f(x) \cos(2n-1)x \, dx + \int_{0}^{\pi} f(x) \cos(2n-1)x \, dx \right) = 0.
$$

对 $b_{2n-1}$ 完全同理：

$$
\begin{aligned}
b_{2n-1} &= \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \sin(2n-1)x \, dx \\
&= \frac{1}{\pi} \int_{-\pi}^{0} + \frac{1}{\pi} \int_{0}^{\pi}.
\end{aligned}
$$

令 $t = x + \pi$，利用 $\sin\bigl[(2n-1)(t - \pi)\bigr] = \sin\bigl[(2n-1)t - (2n-1)\pi\bigr] = \sin(2n-1)t \cdot \underbrace{\cos(2n-1)\pi}_{=-1} - \cos(2n-1)t \cdot \underbrace{\sin(2n-1)\pi}_{=0} = -\sin(2n-1)t$，

同样得到两项相消，$b_{2n-1} = 0$。$\square$

**几何意义**：条件 $f(x+\pi) = f(x)$ 表明函数以 $\pi$ 为周期（即半周期平移不变）。此时函数只包含偶次谐波分量，奇次谐波系数均为零。

---

### 1.7 Fourier 级数的收敛性（Dirichlet 条件）

**定理 16.6（Dirichlet 充分条件）** 设 $f(x)$ 是以 $2\pi$ 为周期的函数，在 $[-\pi, \pi]$ 上满足：
1. 至多只有有限个第一类间断点；
2. 至多只有有限个极值点（即分段单调）。

则 $f$ 的 Fourier 级数处处收敛，且

$$
\frac{a_0}{2} + \sum_{n=1}^{\infty} \bigl(a_n \cos nx + b_n \sin nx\bigr) = \frac{f(x^+) + f(x^-)}{2},
$$

其中 $f(x^+)$ 和 $f(x^-)$ 分别为 $f$ 在 $x$ 处的右极限和左极限。在连续点处，右端等于 $f(x)$；在跳跃间断点处，收敛到左右极限的算术平均值。

---

### 1.8 综合例题

**例题 16.5** 求 $f(x) = \begin{cases} 0, & -\pi < x < 0, \\ 1, & 0 < x < \pi \end{cases}$ 的 Fourier 级数。

**解** 计算系数：

$$
\begin{aligned}
a_0 &= \frac{1}{\pi} \int_{0}^{\pi} 1 \, dx = 1, \\[6pt]
a_n &= \frac{1}{\pi} \int_{0}^{\pi} \cos nx \, dx = \frac{1}{\pi} \left[ \frac{\sin nx}{n} \right]_{0}^{\pi} = 0, \quad n \ge 1, \\[6pt]
b_n &= \frac{1}{\pi} \int_{0}^{\pi} \sin nx \, dx = \frac{1}{\pi} \left[ -\frac{\cos nx}{n} \right]_{0}^{\pi} = \frac{1 - (-1)^n}{n\pi}.
\end{aligned}
$$

当 $n$ 为偶数时 $b_n = 0$；当 $n$ 为奇数，$n = 2k-1$ 时 $b_{2k-1} = \dfrac{2}{(2k-1)\pi}$。

因此

$$
f(x) \sim \frac{1}{2} + \frac{2}{\pi} \sum_{k=1}^{\infty} \frac{\sin(2k-1)x}{2k-1}.
$$

在 $x = 0$ 处，Fourier 级数收敛到 $\dfrac{f(0^+) + f(0^-)}{2} = \dfrac{1 + 0}{2} = \dfrac{1}{2}$，而级数本身在 $x=0$ 时每项为零，和为 $0$。这是因为 Dirichlet 条件在间断点处取算术平均值，但此处 $f$ 的 Fourier 级数实际收敛到 $\frac{1}{2}$，而级数的部分和在 $x=0$ 处确实为 $0$，这提醒我们 Fourier 级数的收敛行为在间断点处需要注意。更仔细的分析表明，上述级数在 $x=0$ 处收敛到 $\frac{1}{2}$ 是通过极限 $x\to 0^+$ 和 $x\to 0^-$ 实现的。

---

**例题 16.6** 利用 Fourier 级数求 $\displaystyle\sum_{n=1}^{\infty} \frac{1}{n^2}$。

**解** 在例题 16.2 中已得：

$$
x^2 = \frac{\pi^2}{3} + 4 \sum_{n=1}^{\infty} \frac{(-1)^n}{n^2} \cos nx, \quad x \in [-\pi, \pi].
$$

令 $x = \pi$：

$$
\pi^2 = \frac{\pi^2}{3} + 4 \sum_{n=1}^{\infty} \frac{(-1)^n}{n^2} \cos n\pi = \frac{\pi^2}{3} + 4 \sum_{n=1}^{\infty} \frac{(-1)^n}{n^2} (-1)^n = \frac{\pi^2}{3} + 4 \sum_{n=1}^{\infty} \frac{1}{n^2}.
$$

解得

$$
\boxed{\; \sum_{n=1}^{\infty} \frac{1}{n^2} = \frac{\pi^2}{6}. \;}
$$

这是数学史上著名的 Basel 问题（Euler, 1735）。

---

**例题 16.7** 求 $f(x) = \begin{cases} x + \pi, & -\pi \le x < 0, \\ \pi - x, & 0 \le x \le \pi \end{cases}$ 的 Fourier 级数。

**解** 注意到 $f$ 是偶函数（$f(-x) = f(x)$），故 $b_n = 0$。

$$
\begin{aligned}
a_0 &= \frac{2}{\pi} \int_{0}^{\pi} (\pi - x) \, dx = \frac{2}{\pi} \left[ \pi x - \frac{x^2}{2} \right]_{0}^{\pi} = \frac{2}{\pi} \cdot \frac{\pi^2}{2} = \pi. \\[6pt]
a_n &= \frac{2}{\pi} \int_{0}^{\pi} (\pi - x) \cos nx \, dx \\
&= \frac{2}{\pi} \left( \pi \int_{0}^{\pi} \cos nx \, dx - \int_{0}^{\pi} x \cos nx \, dx \right) \\
&= \frac{2}{\pi} \left( 0 - \frac{(-1)^n - 1}{n^2} \right) = \frac{2[1 - (-1)^n]}{\pi n^2}.
\end{aligned}
$$

当 $n$ 为偶数时 $a_n = 0$；$n = 2k-1$ 时 $a_{2k-1} = \dfrac{4}{\pi(2k-1)^2}$。

因此

$$
f(x) \sim \frac{\pi}{2} + \frac{4}{\pi} \sum_{k=1}^{\infty} \frac{\cos(2k-1)x}{(2k-1)^2}.
$$

**注**：此结果与例题 16.3(a) 中 $x$ 在 $[0, \pi]$ 上的余弦级数形式相同，这是因为两者偶延拓后在 $[0, \pi]$ 上取值不同但 $[-\pi, \pi]$ 上的 Fourier 级数结构相似。

---

### 习题

1. 验证函数系 $\{1, \cos x, \cos 2x, \cos 3x\}$ 在 $[-\pi, \pi]$ 上的正交性，写出所有非零内积值。

2. 求下列函数在 $[-\pi, \pi]$ 上的 Fourier 级数：
   (a) $f(x) = |x|$；
   (b) $f(x) = \begin{cases} -1, & -\pi < x < 0, \\ 1, & 0 < x < \pi. \end{cases}$

3. 将 $f(x) = \cos^2 x$ 在 $[-\pi, \pi]$ 上展开为 Fourier 级数，并与三角恒等式比较。

4. 将 $f(x) = \pi - x$ 在 $[0, \pi]$ 上分别展开为余弦级数和正弦级数。

5. 将 $f(x) = x(\pi - x)$ 在 $[0, \pi]$ 上展开为正弦级数，并用其结果求 $\displaystyle\sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{(2n-1)^3}$。

6. 求 $f(x) = e^x$ 在 $[-\pi, \pi]$ 上的 Fourier 级数。

7. 设 $f(x)$ 是以 $2\pi$ 为周期的函数，且满足 $f(x + \pi) = -f(x)$。证明 $f$ 的 Fourier 级数中所有偶数阶系数为零（即 $a_{2n} = b_{2n} = 0$）。

8. 利用 Fourier 级数求 $\displaystyle\sum_{n=1}^{\infty} \frac{1}{(2n-1)^2}$。（提示：用例 16.6 的结果）

---

### 本章小结

本节建立了 Fourier 级数的完整基础框架：

| 主题 | 核心内容 | 对应定义/定理 |
|------|----------|--------------|
| 正交性 | 三角函数系在 $[-\pi, \pi]$ 上两两正交 | 定义 16.1-16.2, 定理 16.1 |
| Fourier 系数 | Euler-Fourier 公式 $a_n, b_n$ | 定义 16.3 |
| 奇偶简化 | 奇→正弦级数，偶→余弦级数 | 定理 16.2-16.3 |
| 半幅展开 | 奇延拓/偶延拓后的半区间展开 | 定义 16.4 |
| 一般区间 | $[a, a+2\pi]$ 上的系数公式 | 定理 16.4 |
| 对称性 | 半周期平移导致奇次谐波归零 | 定理 16.5 |
| 收敛性 | Dirichlet 条件 | 定理 16.6 |

下一节将深入讨论 Fourier 级数的收敛性理论和 Gibbs 现象。
