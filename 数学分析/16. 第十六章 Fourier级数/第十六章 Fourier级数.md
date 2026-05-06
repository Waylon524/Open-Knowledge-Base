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


# 第十六章 Fourier级数

---

## §2 Fourier级数的逐点收敛与Dirichlet收敛定理

---

### 2.1 Dirichlet核与部分和的积分表示

在讨论Fourier级数的收敛性之前，我们需要将部分和 $S_n(x)$ 表示为便于分析极限行为的积分形式。核心工具是 **Dirichlet核**。

**定义 16.7（Dirichlet核）** 称
$$
D_n(t) = \frac12 + \sum_{k=1}^n \cos kt
$$
为 $n$ 阶 **Dirichlet核**。它是 $2\pi$ 为周期的偶函数。

---

**定理 16.7（Dirichlet核的闭式）** 对 $t \neq 2k\pi$（$k \in \mathbb{Z}$），有
$$
\boxed{\; D_n(t) = \frac{\sin\!\big(n+\frac12\big)t}{2\sin\frac{t}{2}} \;}
$$
在 $t = 2k\pi$ 处，$D_n(t)$ 按连续延拓取值 $D_n(2k\pi) = n + \frac12$。

**证明** 利用三角恒等式：
$$
2\sin\frac{t}{2} \cdot D_n(t) = 2\sin\frac{t}{2}\Big(\frac12 + \sum_{k=1}^n \cos kt\Big)
= \sin\frac{t}{2} + \sum_{k=1}^n 2\sin\frac{t}{2}\cos kt.
$$

由积化和差 $2\sin\frac{t}{2}\cos kt = \sin\!\big(k+\frac12\big)t - \sin\!\big(k-\frac12\big)t$，代入得：
$$
\begin{aligned}
2\sin\frac{t}{2}\cdot D_n(t) &= \sin\frac{t}{2} + \sum_{k=1}^n \Big[\sin\!\big(k+\frac12\big)t - \sin\!\big(k-\frac12\big)t\Big] \\
&= \sin\frac{t}{2} + \big[\sin\frac{3t}{2} - \sin\frac{t}{2}\big] + \big[\sin\frac{5t}{2} - \sin\frac{3t}{2}\big] + \cdots + \big[\sin\!\big(n+\frac12\big)t - \sin\!\big(n-\frac12\big)t\big] \\
&= \sin\!\big(n+\frac12\big)t.
\end{aligned}
$$

逐项相消后，仅剩 $\sin\big(n+\frac12\big)t$，故 $D_n(t) = \dfrac{\sin\big(n+\frac12\big)t}{2\sin\frac{t}{2}}$。

当 $t \to 0$ 时，利用 $\sin\frac{t}{2} \sim \frac{t}{2}$ 和 $\sin\big(n+\frac12\big)t \sim \big(n+\frac12\big)t$，得 $\displaystyle\lim_{t\to 0} D_n(t) = n + \frac12$。$\square$

**注**：Dirichlet核在 $t=0$ 处有最大值 $n+\frac12$，随着 $n$ 增大，主峰越来越高、越来越窄。这一性质对理解Fourier级数的收敛性和Gibbs现象至关重要。

---

**定理 16.8（部分和的卷积表示）** 设 $f$ 以 $2\pi$ 为周期且在 $[-\pi,\pi]$ 上可积，则 $f$ 的Fourier级数部分和 $S_n(x)$ 可表示为 $f$ 与 Dirichlet核的卷积：
$$
\boxed{\; S_n(x) = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x+t)D_n(t)\,dt = \frac{1}{\pi}\int_{-\pi}^{\pi} f(t)D_n(x-t)\,dt. \;}
$$

**证明** 将 Euler-Fourier 公式
$$
a_k = \frac{1}{\pi}\int_{-\pi}^{\pi} f(u)\cos ku\,du,\qquad
b_k = \frac{1}{\pi}\int_{-\pi}^{\pi} f(u)\sin ku\,du
$$
代入部分和 $S_n(x) = \dfrac{a_0}{2} + \displaystyle\sum_{k=1}^n (a_k\cos kx + b_k\sin kx)$：
$$
\begin{aligned}
S_n(x) &= \frac{1}{2\pi}\int_{-\pi}^{\pi} f(u)\,du + \sum_{k=1}^n \frac{1}{\pi}\int_{-\pi}^{\pi} f(u)\big[\cos kx\cos ku + \sin kx\sin ku\big]\,du \\
&= \frac{1}{\pi}\int_{-\pi}^{\pi} f(u)\Big[\frac12 + \sum_{k=1}^n \cos k(u-x)\Big]\,du \\
&= \frac{1}{\pi}\int_{-\pi}^{\pi} f(u)D_n(u-x)\,du.
\end{aligned}
$$

令 $t = u - x$，利用 $f$ 的周期性，积分区间平移后不变，得 $S_n(x) = \dfrac{1}{\pi}\displaystyle\int_{-\pi}^{\pi} f(x+t)D_n(t)\,dt$。$\square$

**引理 16.1** Dirichlet核的积分为
$$
\frac{1}{\pi}\int_0^{\pi} D_n(t)\,dt = \frac12.
$$

**证明** 由 $D_n(t) = \frac12 + \sum_{k=1}^n \cos kt$，逐项积分得 $\int_0^{\pi} \cos kt\,dt = 0$（$k \ge 1$），故 $\int_0^{\pi} D_n(t)\,dt = \frac{\pi}{2}$，除以 $\pi$ 即得 $\frac12$。$\square$

---

**定理 16.9（部分和与函数值的差积分）** 设 $f$ 以 $2\pi$ 为周期且在 $[-\pi,\pi]$ 上可积，$f(x^+)$ 和 $f(x^-)$ 分别表示 $f$ 在 $x$ 处的右极限和左极限，则
$$
S_n(x) - \frac{f(x^+)+f(x^-)}{2}
= \frac{1}{\pi}\int_0^{\pi} \frac{f(x+t)-f(x^+)}{2\sin\frac{t}{2}}\sin\!\Big(n+\frac12\Big)t\,dt
+ \frac{1}{\pi}\int_0^{\pi} \frac{f(x-t)-f(x^-)}{2\sin\frac{t}{2}}\sin\!\Big(n+\frac12\Big)t\,dt.
$$

**证明** 由定理16.8及 $D_n$ 为偶函数，将 $[-\pi,0]$ 段作代换 $t \to -t$：
$$
S_n(x) = \frac{1}{\pi}\int_0^{\pi} \big[f(x+t)+f(x-t)\big]D_n(t)\,dt.
$$

由引理16.1：
$$
\frac{f(x^+)+f(x^-)}{2} = \frac{1}{\pi}\int_0^{\pi} \big[f(x^+)+f(x^-)\big]D_n(t)\,dt.
$$

两式相减：
$$
S_n(x) - \frac{f(x^+)+f(x^-)}{2}
= \frac{1}{\pi}\int_0^{\pi} \Big[f(x+t)-f(x^+) + f(x-t)-f(x^-)\Big]D_n(t)\,dt.
$$

代入 $D_n(t) = \dfrac{\sin(n+\frac12)t}{2\sin\frac{t}{2}}$，即得所述表达式。$\square$

**注**：定理16.9将收敛性问题转化为形如 $\int_0^{\pi} \varphi(t)\sin\big(n+\frac12\big)t\,dt$ 的积分在 $n\to\infty$ 时的极限行为问题。这是利用 Riemann-Lebesgue 引理研究Fourier级数收敛性的关键步骤。

---

### 2.2 Riemann-Lebesgue引理与Fourier系数的衰减

**定理 16.10（Riemann-Lebesgue引理）** 设 $\varphi(u)$ 在 $[a,b]$ 上可积（或绝对可积），则
$$
\boxed{\;
\lim_{p\to\infty}\int_a^b \varphi(u)\sin pu\,du = 0,\qquad
\lim_{p\to\infty}\int_a^b \varphi(u)\cos pu\,du = 0.
\;}
$$

**证明思路** 分三步：
1. **阶梯函数情形**：若 $\varphi(u) \equiv c$ 在 $[\alpha,\beta]$ 上为常数，则
   $$
   \int_{\alpha}^{\beta} c\sin pu\,du = c\cdot\frac{\cos p\alpha - \cos p\beta}{p} \xrightarrow{p\to\infty} 0.
   $$
   （对 $\cos pu$ 同理。）

2. **可积函数逼近**：对任意可积函数 $\varphi$ 和任意 $\varepsilon>0$，存在阶梯函数 $g$ 使得 $\int_a^b |\varphi(u)-g(u)|\,du < \varepsilon$。

3. **三角不等式估计**：
   $$
   \Big|\int_a^b \varphi(u)\sin pu\,du\Big|
   \le \Big|\int_a^b g(u)\sin pu\,du\Big| + \int_a^b |\varphi(u)-g(u)|\,du.
   $$
   当 $p$ 足够大时第一项 $<\varepsilon$，第二项 $<\varepsilon$，故总积分 $<2\varepsilon$。$\square$

**推论 16.1（Fourier系数趋于零）** 若 $f$ 在 $[-\pi,\pi]$ 上可积，则其Fourier系数满足
$$
\lim_{n\to\infty} a_n = 0,\qquad \lim_{n\to\infty} b_n = 0.
$$

**证明** 在Riemann-Lebesgue引理中取 $\varphi(u)=f(u)/\pi$，$p=n$，$a=-\pi$，$b=\pi$，由 $a_n=\frac1\pi\int_{-\pi}^{\pi}f(u)\cos nu\,du$ 即得 $\lim a_n=0$；对 $b_n$ 同理。$\square$

---

**定理 16.11（Fourier系数与导数系数的关系）** 设 $f$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上连续且有连续导函数 $f'$。记 $f$ 的Fourier系数为 $a_n,b_n$，$f'$ 的Fourier系数为 $a_n',b_n'$，则
$$
\boxed{\;
a_n' = n b_n,\qquad b_n' = -n a_n,\qquad n=1,2,\dots
\;}
$$

**证明** 对 $a_n'$ 使用分部积分，利用周期性 $f(-\pi)=f(\pi)$ 消去边界项：
$$
\begin{aligned}
a_n' &= \frac{1}{\pi}\int_{-\pi}^{\pi} f'(x)\cos nx\,dx \\
&= \frac{1}{\pi}\Big(\big[f(x)\cos nx\big]_{-\pi}^{\pi} + n\int_{-\pi}^{\pi} f(x)\sin nx\,dx\Big) \\
&= \frac{n}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx = n b_n.
\end{aligned}
$$

对 $b_n'$ 同理：
$$
\begin{aligned}
b_n' &= \frac{1}{\pi}\int_{-\pi}^{\pi} f'(x)\sin nx\,dx \\
&= \frac{1}{\pi}\Big(\big[f(x)\sin nx\big]_{-\pi}^{\pi} - n\int_{-\pi}^{\pi} f(x)\cos nx\,dx\Big) \\
&= -\frac{n}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx = -n a_n.
\end{aligned}
$$
$\square$

---

**定理 16.12（Fourier系数的衰减速度）** 设 $f$ 以 $2\pi$ 为周期。
1. 若 $f$ 连续且 $f(-\pi)=f(\pi)$，则 $a_n,b_n = o(1/n)$；
2. 若 $f \in C^k$（$k$ 阶连续可导），则 $a_n,b_n = o(1/n^k)$。

**证明** （1）由 $f$ 连续且 $f(-\pi)=f(\pi)$，可定义周期函数 $f$ 处处连续。记 $F(x)=\int_0^x f(t)\,dt$ 虽然不一定有周期，但分部积分法可用标准方法处理。更直接的方式是利用定理16.11的逆推：

设 $f$ 连续且 $f(-\pi)=f(\pi)$，则 $f$ 可视为某个周期函数的原函数。考虑 $f$ 的Fourier系数 $a_n,b_n$。由定理16.11对 $f'$ 的系数关系（若 $f$ 可微），但此处 $f$ 不一定可微。不过，利用分部积分直接对 $a_n$ 做一次分部积分：
$$
a_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx
= \frac{1}{\pi}\Big(\big[f(x)\frac{\sin nx}{n}\big]_{-\pi}^{\pi} - \frac{1}{n}\int_{-\pi}^{\pi} f'(x)\sin nx\,dx\Big).
$$
若 $f$ 连续且 $f(-\pi)=f(\pi)$，边界项为零（$\sin(\pm n\pi)=0$）。被积函数中 $f'(x)$ 不一定存在，但若 $f$ 连续可微，则上式给出 $a_n = O(1/n)$。再由Riemann-Lebesgue引理（$f'$ 的系数趋于零）得 $n a_n \to 0$，即 $a_n = o(1/n)$。对 $b_n$ 同理。

（2）若 $f \in C^k$，反复利用分部积分 $k$ 次：
$$
a_n = \frac{(-1)^k}{\pi n^k}\int_{-\pi}^{\pi} f^{(k)}(x)\big(\text{$\cos nx$ 或 $\sin nx$}\big)\,dx + \text{（边界项）},
$$
其中边界项因 $f^{(j)}(-\pi)=f^{(j)}(\pi)$（$j=0,1,\dots,k-1$）而消去。由Riemann-Lebesgue引理，积分部分趋于零，故 $a_n,b_n = o(1/n^k)$。$\square$

**Fourier系数衰减与函数光滑性的对应关系**：

| 函数光滑程度 | Fourier系数衰减速度 | 级数收敛性 |
|-------------|---------------------|-----------|
| 可积（仅保证） | $a_n,b_n \to 0$ | 无具体衰减率 |
| 有界变差/分段连续 | $a_n,b_n = O(1/n)$ | 逐点收敛 |
| 连续且端点值相等 | $a_n,b_n = o(1/n)$ | 可能一致收敛 |
| $C^1$ 连续可微 | $a_n,b_n = o(1/n)$ | 绝对一致收敛 |
| $C^k$ 光滑 | $a_n,b_n = o(1/n^k)$ | 快速收敛 |
| 无穷次可微 | $a_n,b_n = o(1/n^m)$ 对任意 $m$ | 超快收敛 |
| 解析 | 指数衰减 | 极快收敛 |

这一关系是Fourier分析的核心结论之一，也是谱方法数值计算的理论基础——函数越光滑，Fourier级数收敛越快。

---

### 2.3 Fourier级数的收敛判别法

利用Dirichlet核的积分表示和Riemann-Lebesgue引理，可以建立Fourier级数收敛的充分条件。

**定理 16.13（Dirichlet-Jordan收敛定理）** 设 $f(x)$ 是以 $2\pi$ 为周期的函数，在 $[-\pi,\pi]$ 上满足：
1. 至多只有有限个第一类间断点；
2. **分段单调**，即函数可分成有限个单调区间（至多只有有限个极值点）。

则 $f$ 的Fourier级数处处收敛，且
$$
\boxed{\;
\frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n\cos nx + b_n\sin nx) = \frac{f(x^+) + f(x^-)}{2}.
\;}
$$

在连续点处，$f(x^+)=f(x^-)=f(x)$，故收敛到 $f(x)$；在跳跃间断点处，收敛到左右极限的算术平均值。

**注**：定理16.13与§1中的定理16.6（Dirichlet充分条件）表述完全一致。§1将之作为Fourier展开时的收敛性说明，§2正式赋予其完整的理论名称——Dirichlet-Jordan定理，以便与其他判别法区分。

**证明概要**：利用定理16.9将收敛判定转化为证明
$$
\lim_{n\to\infty} \frac{1}{\pi}\int_0^{\pi} \frac{f(x+t)-f(x^+)}{2\sin\frac{t}{2}}\sin\!\Big(n+\frac12\Big)t\,dt = 0
$$
以及右极限部分的对应项为零。在分段单调条件下，函数 $\dfrac{f(x+t)-f(x^+)}{\sin\frac{t}{2}}$ 在 $t=0$ 附近具有可控性，再结合Riemann-Lebesgue引理及积分第二中值定理可完成证明。$\square$

---

**定理 16.14（Dini-Lipschitz判别法）** 设 $f$ 以 $2\pi$ 为周期。若在点 $x$ 处存在常数 $\alpha \in (0,1]$ 和 $\delta > 0$，使得当 $|h| < \delta$ 时
$$
|f(x+h) - f(x)| \le L\,|h|^{\alpha},
$$
则 $f$ 的Fourier级数在 $x$ 处收敛到 $f(x)$。

换言之，若 $f$ 在 $x$ 附近满足 **Hölder条件**（$\alpha=1$ 时即为 **Lipschitz条件**），则Fourier级数在该点收敛。

**证明思路**：在定理16.9的表达式中，令 $\varphi(t) = \dfrac{f(x+t)-f(x^+)}{2\sin\frac{t}{2}}$。由Hölder条件，在 $t$ 充分小时 $|f(x+t)-f(x^+)| \le L t^{\alpha}$，且 $\sin\frac{t}{2} \sim \frac{t}{2}$，故 $|\varphi(t)| \le C t^{\alpha-1}$。由于 $\alpha > 0$，函数 $\varphi(t)$ 在 $t=0$ 附近可积（$\int_0^{\delta} t^{\alpha-1}\,dt < \infty$）。由Riemann-Lebesgue引理即得积分趋于零。$\square$

---

**两类判别法的本质区别**：

| 判别法 | 条件性质 | 覆盖范围 | 典型满足函数 | 典型不满足函数 |
|--------|----------|----------|-------------|---------------|
| Dirichlet-Jordan | 整体性（全区间分段单调） | 有界变差函数类 | $f_1(x)=1/\ln|x|$（单调连续） | $f_2(x)=x\cos(\pi/2x)$（无穷次振荡） |
| Dini-Lipschitz | 局部性（单点附近Hölder连续） | Hölder连续函数类 | $f_2(x)=x\cos(\pi/2x)$（Lipschitz） | $f_1(x)=1/\ln|x|$（非Hölder） |

两类判别法 **互不包含**，各自覆盖了不同的函数类。

---

**例题 16.8** 判断以下函数在 $x=0$ 处Fourier级数的收敛性：
$$
f_1(x)=\begin{cases}
\displaystyle\frac{1}{\ln\frac{|x|}{2\pi}}, & x\neq 0,\\[6pt]
0, & x=0,
\end{cases} \qquad
f_2(x)=\begin{cases}
x\cos\frac{\pi}{2x}, & x\neq 0,\\[4pt]
0, & x=0.
\end{cases}
$$

**解** 
- **$f_1$**：$f_1$ 在 $[-\pi,\pi]$ 上连续（$\lim_{x\to0}\frac{1}{\ln(|x|/2\pi)}=0$），且对 $x>0$ 有 $f_1'(x)=-\frac{1}{x(\ln(x/2\pi))^2}<0$，故 $f_1$ 单调递减，**满足Dirichlet-Jordan判别法**。但
  $$
  \frac{|f_1(0+u)-f_1(0^+)|}{u^{\alpha}} = \frac{1}{u^{\alpha}|\ln\frac{u}{2\pi}|} \xrightarrow{u\to0^+} \infty,
  $$
  故 $f_1$ **不满足Dini-Lipschitz判别法**。

- **$f_2$**：$|f_2(x)-f_2(0)| = |x\cos\frac{\pi}{2x}| \le |x|$，故 $f_2$ 在 $x=0$ 处满足Lipschitz条件（$\alpha=1$），**满足Dini-Lipschitz判别法**。但 $f_2$ 在 $x=0$ 的任何邻域内振荡无穷多次（导函数变号无穷次），不是分段单调函数，故 **不满足Dirichlet-Jordan判别法**。

两个函数均能以 $2\pi$ 为周期延拓，且 $f_1$ 和 $f_2$ 在 $x=0$ 处均有收敛的Fourier级数（$f_1$ 由Dirichlet-Jordan保证，$f_2$ 由Dini-Lipschitz保证）。这一例子充分说明了两类判别法的独立性。

---

### 2.4 Parseval等式与Bessel不等式

**定理 16.15（Bessel不等式）** 设 $f$ 在 $[-\pi,\pi]$ 上平方可积，其Fourier系数为 $a_n,b_n$，则
$$
\boxed{\;
\frac{a_0^2}{2} + \sum_{n=1}^{\infty} (a_n^2 + b_n^2) \le \frac{1}{\pi}\int_{-\pi}^{\pi} f^2(x)\,dx.
\;}
$$

**证明** 考虑三角多项式 $S_n(x) = \frac{a_0}{2} + \sum_{k=1}^n (a_k\cos kx + b_k\sin kx)$。计算 $\int_{-\pi}^{\pi} (f(x)-S_n(x))^2\,dx \ge 0$ 并展开，利用正交性化简即得。$\square$

Bessel不等式表明Fourier系数的平方和受到函数 $L^2$ 模的控制，从而 $\sum (a_n^2+b_n^2)$ 收敛。

---

**定理 16.16（Parseval等式）** 设 $f$ 在 $[-\pi,\pi]$ 上平方可积，则Bessel不等式中等号成立：
$$
\boxed{\;
\frac{1}{\pi}\int_{-\pi}^{\pi} f^2(x)\,dx = \frac{a_0^2}{2} + \sum_{n=1}^{\infty} (a_n^2 + b_n^2).
\;}
$$

**证明** 证明Parseval等式需要利用Fourier级数在 $L^2$ 意义下的收敛性（即均方收敛）。可以证明：
$$
\lim_{n\to\infty} \int_{-\pi}^{\pi} |f(x) - S_n(x)|^2\,dx = 0,
$$
展开后即得Parseval等式。$\square$

**特化形式**：
- **正弦级数**（$a_n=0$）：
  $$
  \frac{2}{\pi}\int_0^{\pi} f^2(x)\,dx = \sum_{n=1}^{\infty} b_n^2.
  $$
- **余弦级数**（$b_n=0$）：
  $$
  \frac{1}{\pi}\int_{-\pi}^{\pi} f^2(x)\,dx = \frac{a_0^2}{2} + \sum_{n=1}^{\infty} a_n^2.
  $$

---

**例题 16.9（用Parseval等式求 $\displaystyle\sum_{n=1}^{\infty}\frac{1}{n^2}$）** 利用Parseval等式重新推导Basel问题的结果。

**解** 取 $f(x)=x$，$x\in(-\pi,\pi)$。$f$ 是奇函数，其Fourier级数为（见§1例题16.1）：
$$
x \sim 2\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n}\sin nx,
$$
系数为 $a_n=0$，$b_n = \dfrac{2(-1)^{n+1}}{n}$。

代入Parseval等式（正弦级数形式）：
$$
\frac{2}{\pi}\int_0^{\pi} x^2\,dx = \sum_{n=1}^{\infty} b_n^2 = \sum_{n=1}^{\infty} \frac{4}{n^2}.
$$

计算左端：$\displaystyle\frac{2}{\pi}\int_0^{\pi} x^2\,dx = \frac{2}{\pi}\cdot\frac{\pi^3}{3} = \frac{2\pi^2}{3}$。

因此 $\displaystyle\frac{2\pi^2}{3} = 4\sum_{n=1}^{\infty}\frac{1}{n^2}$，解得
$$
\boxed{\; \sum_{n=1}^{\infty}\frac{1}{n^2} = \frac{\pi^2}{6}. \;}
$$

**注**：§1例题16.6通过代入 $x=\pi$ 的方式也得到了相同结果。Parseval等式提供了另一条途径，且更具一般性——它不仅适用于特定的 $x$ 值，而是利用了函数的整体 $L^2$ 模。

---

**例题 16.10** 利用Parseval等式求 $\displaystyle\sum_{n=1}^{\infty}\frac{1}{(2n-1)^2}$。

**解** 取 $f(x)=\begin{cases} -1, & -\pi < x < 0, \\ 1, & 0 < x < \pi \end{cases}$（方波函数）。其Fourier级数为（见§1例题16.5）：
$$
f(x) \sim \frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\sin(2k-1)x}{2k-1},
$$
系数为 $b_{2k-1} = \dfrac{4}{\pi(2k-1)}$，其余系数为零。

由Parseval等式：
$$
\frac{1}{\pi}\int_{-\pi}^{\pi} f^2(x)\,dx = \sum_{k=1}^{\infty} b_{2k-1}^2 = \frac{16}{\pi^2}\sum_{k=1}^{\infty}\frac{1}{(2k-1)^2}.
$$

左端：$\frac{1}{\pi}\int_{-\pi}^{\pi} 1\,dx = 2$。

因此 $\displaystyle 2 = \frac{16}{\pi^2}\sum_{k=1}^{\infty}\frac{1}{(2k-1)^2}$，解得
$$
\boxed{\; \sum_{k=1}^{\infty}\frac{1}{(2k-1)^2} = \frac{\pi^2}{8}. \;}
$$

---

### 2.5 Gibbs现象

当函数 $f$ 有跳跃间断点时，即使Fourier级数在每点都收敛（由Dirichlet-Jordan定理保证），在间断点附近也会出现一种特殊的振荡行为——Gibbs现象。

**定义 16.8（Gibbs现象）** 设 $f$ 在 $x_0$ 处有跳跃间断点，$S_n(x)$ 为其Fourier级数的部分和。Gibbs现象指以下特征：
1. 在间断点两侧，$S_n(x)$ 以振荡方式逼近 $f(x)$；
2. 紧邻间断点的位置，$S_n(x)$ 的值 **超出** $f$ 的极限值，形成"过冲"（overshoot）和"欠冲"（undershoot）；
3. 随着 $n$ 增大，过冲/欠冲的宽度趋于零，但 **高度趋近于一个非零常数**。

**数学本质**：Gibbs现象来源于Dirichlet核的旁瓣效应，是Fourier级数用连续函数（三角多项式）逼近间断函数时的固有特征。其本质是 **非一致收敛性** 的体现——在间断点附近，无论 $n$ 多大，总存在点使得 $|S_n(x)-f(x)|$ 不小于一个正常数。

---

**Gibbs现象的定量分析（方波函数）**

考虑方波函数：
$$
f(x)=\begin{cases}
-1, & -\pi < x < 0,\\
0, & x=0,\\
1, & 0 < x < \pi,
\end{cases}
$$
以 $2\pi$ 为周期延拓。其Fourier级数为：
$$
f(x) \sim \frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\sin(2k-1)x}{2k-1}.
$$

部分和 $S_n(x)$ 在 $x=0$ 右侧附近的行为由Dirichlet核决定。可以证明，过冲的极限值为：
$$
\lim_{n\to\infty} S_n\!\Big(\frac{\pi}{n}\Big) = \frac{2}{\pi}\int_0^{\pi} \frac{\sin t}{t}\,dt \approx 1.17898.
$$

跳跃度 $J = f(0^+)-f(0^-) = 1-(-1) = 2$，$f(0^+)=1$。过冲量相对于跳跃度的百分比：
$$
\frac{S_n(\pi/n) - f(0^+)}{J} \times 100\% \approx \frac{1.17898 - 1}{2} \times 100\% \approx 8.95\%.
$$

即过冲约占跳跃度的 **约9%**，这是一个与具体函数无关的普适常数（只要函数在该点有跳跃）。

---

**Gibbs现象与一致收敛性**

1. **是否随 $n\to\infty$ 消失？** 不会。过冲的位置 $x_n \sim \pi/n$ 趋于零（收缩到间断点附近），但过冲的振幅趋于正常数 $1.179$，**不会**趋于 $f(0^+)=1$。因此Gibbs现象不会消失。

2. **对逐点收敛的影响**：对任意固定的 $x\neq 0$，当 $n$ 足够大时，$|x| > \pi/n$，过冲区域已被压缩到比 $x$ 更靠近原点的位置，故 $S_n(x) \to f(x)$。因此 **逐点收敛性不受影响**。

3. **对一致收敛的影响**：由于存在过冲，$\displaystyle\sup_{x\in(-\pi,\pi)} |S_n(x)-f(x)|$ 不趋于零（始终存在约 $0.179$ 的偏差）。因此Fourier级数 **不可能在包含间断点的区间上一致收敛**。这从另一角度证实了定理16.13的结论——Dirichlet-Jordan条件只能保证逐点收敛，不能保证一致收敛。

---

**例题 16.11** 设
$$
f(x) = \begin{cases}
0, & -\pi \le x < 0, \\
x, & 0 \le x \le \pi.
\end{cases}
$$
（1）求其Fourier级数；（2）利用Dirichlet收敛定理写出在 $x=-\pi, 0, \pi/2, \pi$ 处的和函数值；（3）Fourier级数在 $[-\pi,\pi]$ 上是否一致收敛？

**解** （1）$f$ 在 $[-\pi,\pi]$ 上分段光滑，满足Dirichlet-Jordan条件。计算Fourier系数：
$$
a_0 = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\,dx = \frac{1}{\pi}\int_0^{\pi} x\,dx = \frac{\pi}{2}.
$$

对 $n\ge1$：
$$
\begin{aligned}
a_n &= \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx = \frac{1}{\pi}\int_0^{\pi} x\cos nx\,dx \\
&= \frac{1}{\pi}\Big(\Big[\frac{x\sin nx}{n}\Big]_0^{\pi} - \frac{1}{n}\int_0^{\pi} \sin nx\,dx\Big) \\
&= \frac{1}{\pi}\Big(0 + \frac{1}{n}\Big[\frac{\cos nx}{n}\Big]_0^{\pi}\Big) = \frac{\cos n\pi - 1}{\pi n^2} = \frac{(-1)^n-1}{\pi n^2}.
\end{aligned}
$$

$$
\begin{aligned}
b_n &= \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx = \frac{1}{\pi}\int_0^{\pi} x\sin nx\,dx \\
&= \frac{1}{\pi}\Big(\Big[-\frac{x\cos nx}{n}\Big]_0^{\pi} + \frac{1}{n}\int_0^{\pi} \cos nx\,dx\Big) \\
&= \frac{1}{\pi}\Big(-\frac{\pi(-1)^n}{n} + 0\Big) = \frac{(-1)^{n+1}}{n}.
\end{aligned}
$$

因此
$$
f(x) \sim \frac{\pi}{4} + \sum_{n=1}^{\infty}\Big(\frac{(-1)^n-1}{\pi n^2}\cos nx + \frac{(-1)^{n+1}}{n}\sin nx\Big).
$$

（2）由Dirichlet-Jordan收敛定理：
- $x=-\pi$（周期延拓后为跳跃间断点，左极限 $f(\pi^-)=\pi$，右极限 $f(-\pi^+)=0$）：收敛值 $\dfrac{\pi+0}{2} = \dfrac{\pi}{2}$。
- $x=0$（连续点，$f(0)=0$）：收敛值 $0$。
- $x=\pi/2$（连续点，$f(\pi/2)=\pi/2$）：收敛值 $\pi/2$。
- $x=\pi$（跳跃间断点，左极限 $f(\pi^-)=\pi$，右极限 $f(\pi^+)=f(-\pi^+)=0$）：收敛值 $\dfrac{\pi+0}{2} = \dfrac{\pi}{2}$。

（3）**不一致收敛**。$f$ 的周期延拓在 $x=\pm\pi$ 处有跳跃间断点。若Fourier级数在 $[-\pi,\pi]$ 上一致收敛，由于每项连续，和函数应连续，但 $f$ 在端点不连续，矛盾。故Fourier级数不可能在 $[-\pi,\pi]$ 上一致收敛。

实际上，对任意 $\delta>0$，级数在 $[-\pi+\delta,\pi-\delta]$ 上内闭一致收敛。端点附近的非一致收敛正是Gibbs现象的表现。

---

### 习题

1. 证明Dirichlet核 $D_n(t)$ 满足：
   （1）$\displaystyle\frac{1}{\pi}\int_{-\pi}^{\pi} D_n(t)\,dt = 1$；
   （2）$\displaystyle\int_{-\pi}^{\pi} |D_n(t)|\,dt \sim \frac{4}{\pi^2}\ln n$（$n\to\infty$，提示：利用 $|D_n(t)|$ 的估计）。

2. 利用Riemann-Lebesgue引理证明：
   $$
   \lim_{n\to\infty} \int_0^{\pi} \frac{\cos nx}{x}\,dx = 0.
   $$

3. 设 $f(x)$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上分段光滑。证明：若 $f$ 在 $x_0$ 处连续，则其Fourier级数在 $x_0$ 处收敛到 $f(x_0)$。

4. 判断以下函数在 $x=0$ 处满足哪种判别法的条件（Dirichlet-Jordan, Dini-Lipschitz, 或两者皆满足）：
   （1）$f(x) = |x|^{\frac12}$；
   （2）$f(x) = \begin{cases} x\sin\frac{1}{x}, & x\neq 0, \\ 0, & x=0; \end{cases}$
   （3）$f(x) = \sqrt[3]{x}$。

5. 设 $f(x) = \begin{cases} 1, & 0 < x < \pi, \\ 0, & \pi < x < 2\pi \end{cases}$，以 $2\pi$ 为周期延拓。
   （1）求其Fourier级数；
   （2）利用Dirichlet收敛定理写出在 $x=0,\pi$ 处的和函数值；
   （3）Fourier级数在何处有一致收敛性？何处没有？

6. 设 $f(x) = x^2$ 在 $(-\pi,\pi)$ 上，利用Parseval等式求 $\displaystyle\sum_{n=1}^{\infty}\frac{1}{n^4}$。（提示：$f$ 的Fourier系数见§1例题16.2）

7. 利用Fourier系数衰减理论说明：若 $f$ 以 $2\pi$ 为周期且 $f(x) = \sum_{n=1}^{\infty} \frac{\sin nx}{n}$，则 $f$ 是分段连续的，且其周期延拓在 $x=0$ 处有跳跃间断点。并判断该Fourier级数是否能一致收敛。

8. **（选做）** 设 $f$ 以 $2\pi$ 为周期，且 $f\in C^2$。证明其Fourier级数在 $\mathbb{R}$ 上绝对一致收敛。

---

### 本章小结

本节深入研究了Fourier级数的收敛性理论，是§1展开理论的延续与深化。

| 主题 | 核心内容 | 对应定义/定理 |
|------|----------|--------------|
| Dirichlet核 | 闭式 $D_n(t)=\frac{\sin(n+\frac12)t}{2\sin\frac{t}{2}}$ | 定义16.7, 定理16.7 |
| 卷积表示 | $S_n(x)=\frac1\pi\int_{-\pi}^\pi f(x+t)D_n(t)\,dt$ | 定理16.8 |
| Riemann-Lebesgue | $\int\varphi\sin pu\to0$ ($p\to\infty$) | 定理16.10 |
| 系数衰减 | $f\in C^k \implies a_n,b_n=o(1/n^k)$ | 定理16.11-16.12 |
| Dirichlet-Jordan | 分段单调 $\Rightarrow$ 逐点收敛 | 定理16.13 |
| Dini-Lipschitz | Hölder连续 $\Rightarrow$ 逐点收敛 | 定理16.14 |
| Parseval等式 | $L^2$ 模 = 系数平方和 | 定理16.15-16.16 |
| Gibbs现象 | 间断点处过冲约9%，非一致收敛 | 定义16.8 |

本节的核心思想是：**Fourier级数的收敛性分析远比幂级数复杂**，需要借助Dirichlet核、Riemann-Lebesgue引理和多种专门的收敛判别法。函数的光滑程度决定了Fourier系数衰减速度和级数收敛性质——这是Fourier分析中贯穿始终的核心线索。


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
