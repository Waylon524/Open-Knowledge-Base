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
