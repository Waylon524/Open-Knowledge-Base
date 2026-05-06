# 09. 调和函数与Laplace方程

> 所属章节：第十四章 曲线积分、曲面积分与场论  |  文件序号：09  |  难度：进阶
> 常见混淆点：误以为调和函数必须满足 $\nabla^2 u = 0$ 在整个空间处处成立而忽略奇点（如 $1/r$ 仅在 $\mathbb{R}^3\setminus\{0\}$ 上调和）；混淆方向导数 $\partial u/\partial \mathbf{n}$ 与偏导数 $\partial u/\partial x$——方向导数是沿任意方向的导数，偏导数仅是沿坐标轴方向的特例；混淆Green公式的环量形式与散度形式——两者数学上等价，但物理意义截然不同；误以为 $\oint_C (\partial u/\partial \mathbf{n})\,ds = 0$ 对任何函数都成立——仅当 $u$ 是调和函数时成立；遗忘逆命题证明中连续性条件的关键作用——需要 $\nabla^2 u$ 连续才能用反证法

---

## 1. 学习目标与先修前置

### 学习目标

- 掌握调和函数的定义及其与Laplace方程 $\nabla^2 u = 0$ 的关系
- 掌握方向导数的定义 $\partial u/\partial \mathbf{n} = \nabla u\cdot\mathbf{n}$ 及其与法向导数的关系
- 掌握Green公式的通量/散度形式（二维散度定理）及其与环量形式的等价性
- 理解并掌握调和函数的法向导数通量性质 $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_D \nabla^2 u\,dA$
- 掌握基本调和函数 $1/r$ 的梯度、Laplacian计算及其与Gauss公式的结合应用
- 理解有势无源场的势函数必为调和函数的结论

### 先修知识

- **文件08（第十四章）**：Nabla算子框架（梯度/散度/旋度），Laplace算子 $\nabla^2 = \nabla\cdot(\nabla)$ 的定义（第3.4节，定理14.21前），有势场定义与判定
- **文件03（第十四章）**：Green公式的环量形式（定理14.7），正向边界，标准四步工作流
- **文件06（第十四章）**：Gauss公式（定理14.14），散度的定义与物理意义，源/汇分类
- **第十二章**：偏导数计算、链式法则、乘积法则、混合偏导对称性（Clairaut定理）

---

## 2. Laplace算子复习与调和函数定义 (Type A-1)

### 2.1 Laplace算子的定义复习

从文件08（第3.4节）我们知道，Laplace算子（Laplacian）是梯度与散度的复合运算：

**定义 14.28（Laplace算子 / Laplacian）**：设 $u(x,y,z)$ 具有二阶连续偏导数。$u$ 的Laplacian定义为：
$$\boxed{\nabla^2 u = \nabla\cdot(\nabla u) = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}}$$

Laplace算子将**标量函数**映射为**标量函数**：输入 $u$（标量场），输出 $\nabla^2 u$（标量场）。

对于二维情形（仅依赖 $x,y$），Laplacian退化为：
$$\nabla^2 u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}$$

### 2.2 Laplace方程与调和函数

**定义 14.29（Laplace方程 / Laplace's Equation）**：方程
$$\boxed{\nabla^2 u = 0}$$
称为**Laplace方程**（或**调和方程**）。它是数学物理中最基本的偏微分方程之一，出现在静电场（电势满足 $\nabla^2 V = 0$）、稳态热传导（温度满足 $\nabla^2 T = 0$）、不可压缩无旋流体的速度势等众多物理场景中。

**定义 14.30（调和函数 / Harmonic Function）**：设 $\Omega \subset \mathbb{R}^n$ 是开区域，$u: \Omega \to \mathbb{R}$ 在 $\Omega$ 上具有二阶连续偏导数。若 $u$ 在 $\Omega$ 上处处满足Laplace方程：
$$\boxed{\nabla^2 u = 0 \quad (\forall\, \mathbf{x} \in \Omega)}$$
则称 $u$ 是 $\Omega$ 上的**调和函数**（harmonic function）。

**判断调和函数的标准流程**：计算 $\nabla^2 u$ 的表达式并检查它是否在定义域内恒为零。

### 2.3 典型调和函数示例

**例 1**（线性函数）：$u(x,y) = ax + by + c$（$a,b,c$ 为常数）
$$\frac{\partial^2 u}{\partial x^2} = 0,\quad \frac{\partial^2 u}{\partial y^2} = 0 \quad\Rightarrow\quad \nabla^2 u = 0$$
任何线性函数都是调和的。

---

**例 2**（二次调和多项式）：$u(x,y) = x^2 - y^2$
$$\frac{\partial^2 u}{\partial x^2} = 2,\quad \frac{\partial^2 u}{\partial y^2} = -2 \quad\Rightarrow\quad \nabla^2 u = 2 + (-2) = 0$$

**例 3**（三次调和多项式）：$u(x,y) = x^3 - 3xy^2$
$$\frac{\partial u}{\partial x} = 3x^2 - 3y^2,\quad \frac{\partial^2 u}{\partial x^2} = 6x$$
$$\frac{\partial u}{\partial y} = -6xy,\quad \frac{\partial^2 u}{\partial y^2} = -6x$$
$$\nabla^2 u = 6x + (-6x) = 0$$

**例 4**（指数-三角型）：$u(x,y) = e^x\sin y$
$$\frac{\partial u}{\partial x} = e^x\sin y,\quad \frac{\partial^2 u}{\partial x^2} = e^x\sin y$$
$$\frac{\partial u}{\partial y} = e^x\cos y,\quad \frac{\partial^2 u}{\partial y^2} = -e^x\sin y$$
$$\nabla^2 u = e^x\sin y + (-e^x\sin y) = 0$$

---

**例 5**（三维调和多项式）：$u(x,y,z) = x^2 + y^2 - 2z^2$
$$\frac{\partial^2 u}{\partial x^2} = 2,\quad \frac{\partial^2 u}{\partial y^2} = 2,\quad \frac{\partial^2 u}{\partial z^2} = -4$$
$$\nabla^2 u = 2 + 2 + (-4) = 0$$

**例 6**（对数调和函数）：$u(x,y) = \ln(x^2 + y^2)$，定义域 $\mathbb{R}^2 \setminus \{(0,0)\}$
$$\frac{\partial u}{\partial x} = \frac{2x}{x^2+y^2},\quad \frac{\partial^2 u}{\partial x^2} = \frac{2(y^2-x^2)}{(x^2+y^2)^2}$$
$$\frac{\partial u}{\partial y} = \frac{2y}{x^2+y^2},\quad \frac{\partial^2 u}{\partial y^2} = \frac{2(x^2-y^2)}{(x^2+y^2)^2}$$
$$\nabla^2 u = \frac{2(y^2-x^2)}{(x^2+y^2)^2} + \frac{2(x^2-y^2)}{(x^2+y^2)^2} = 0$$
注意：在原点的奇点不影响 $\mathbb{R}^2\setminus\{(0,0)\}$ 上的调和性。

### 2.4 调和函数的线性叠加性质

**定理 14.31（调和函数的线性叠加）**：设 $u_1, u_2$ 是区域 $\Omega$ 上的调和函数，$c_1, c_2$ 是任意实常数。则它们的线性组合 $u = c_1u_1 + c_2u_2$ 也是 $\Omega$ 上的调和函数。

**证明**：由Laplace算子的线性性质（微分算子的线性性）：
$$\nabla^2 (c_1u_1 + c_2u_2) = c_1\nabla^2 u_1 + c_2\nabla^2 u_2 = c_1\cdot 0 + c_2\cdot 0 = 0 \quad \blacksquare$$

**推论**：全体调和函数构成一个**向量空间**（在线性运算下封闭）。

### 2.5 二次调和多项式的系数条件

考虑一般二次多项式：
$$u(x,y) = ax^2 + bxy + cy^2 + dx + ey + f$$

计算Laplacian：
$$\frac{\partial^2 u}{\partial x^2} = 2a,\quad \frac{\partial^2 u}{\partial y^2} = 2c$$
$$\nabla^2 u = 2a + 2c = 2(a + c)$$

$u$ 为调和函数当且仅当 $a + c = 0$。条件与 $b, d, e, f$ 无关。因此，含混合项的二次调和多项式如 $u(x,y) = 2xy$（取 $a=c=0, b=2$）也是一个调和函数。

验证：$\partial^2(2xy)/\partial x^2 = 0$，$\partial^2(2xy)/\partial y^2 = 0$，所以 $\nabla^2(2xy) = 0$。

---

## 3. 方向导数的定义与公式 (Type A-2)

### 3.1 方向导数的概念

在分析沿曲面的法向通量时，我们需要计算函数沿任意方向的瞬时变化率——这就是**方向导数**。

**定义 14.32（方向导数 / Directional Derivative）**：设 $u(x,y,z)$ 在点 $M$ 处可微，$\mathbf{n} = (n_x, n_y, n_z)$ 是一个**单位向量**。$u$ 在 $M$ 点沿方向 $\mathbf{n}$ 的**方向导数**定义为：
$$\boxed{\frac{\partial u}{\partial \mathbf{n}} = \nabla u \cdot \mathbf{n} = \lim_{t\to 0} \frac{u(M + t\mathbf{n}) - u(M)}{t}}$$

方向导数度量了函数在给定方向上的瞬时变化率。

**法向导数**（normal derivative）是方向导数的特例——当 $\mathbf{n}$ 取为曲面（或曲线）的法向量时，方向导数 $\partial u/\partial \mathbf{n}$ 称为**法向导数**。这在后续讨论通量时极为重要。

### 3.2 展开公式

将方向导数用分量形式展开：
$$\boxed{\frac{\partial u}{\partial \mathbf{n}} = \frac{\partial u}{\partial x}\,n_x + \frac{\partial u}{\partial y}\,n_y + \frac{\partial u}{\partial z}\,n_z}$$

其中 $\mathbf{n} = (n_x, n_y, n_z) = (\cos\alpha,\; \cos\beta,\; \cos\gamma)$ 是单位向量，$\alpha,\beta,\gamma$ 是 $\mathbf{n}$ 与三个坐标轴正方向的夹角（方向余弦）。因此也可写为：
$$\frac{\partial u}{\partial \mathbf{n}} = \frac{\partial u}{\partial x}\cos\alpha + \frac{\partial u}{\partial y}\cos\beta + \frac{\partial u}{\partial z}\cos\gamma$$

### 3.3 计算示例

**例 7**：设 $u(x,y) = x^2 + y^2$，$\mathbf{n} = \left(\frac{\sqrt{2}}{2},\; \frac{\sqrt{2}}{2}\right)$，求 $\partial u/\partial \mathbf{n}$ 在点 $(1,0)$ 处的值。

**解**：
$$\nabla u = (2x,\; 2y),\quad \nabla u(1,0) = (2,\; 0)$$
$$\frac{\partial u}{\partial \mathbf{n}} = \nabla u\cdot\mathbf{n} = (2,\;0)\cdot\left(\frac{\sqrt{2}}{2},\;\frac{\sqrt{2}}{2}\right) = 2\cdot\frac{\sqrt{2}}{2} + 0\cdot\frac{\sqrt{2}}{2} = \sqrt{2}$$

**物理意义**：在点 $(1,0)$ 处沿 $45^\circ$ 方向，函数 $u$ 的变化率为 $\sqrt{2}$。

---

**例 8**：设 $u(x,y,z) = 1/r$，$r = \sqrt{x^2 + y^2 + z^2}$，求 $u$ 沿径向向外的方向导数（即 $\mathbf{n} = (x,y,z)/r$）。

**解**：
$$\nabla\left(\frac{1}{r}\right) = \left(-\frac{x}{r^3},\; -\frac{y}{r^3},\; -\frac{z}{r^3}\right)$$
$$\frac{\partial}{\partial \mathbf{n}}\left(\frac{1}{r}\right) = \nabla\left(\frac{1}{r}\right)\cdot\frac{(x,y,z)}{r} = -\frac{x^2 + y^2 + z^2}{r^4} = -\frac{r^2}{r^4} = -\frac{1}{r^2}$$

---

## 4. 二维散度定理——Green公式的通量形式 (Type A-3)

### 4.1 从环量到通量

文件03（定理14.7）给出了Green公式的**环量形式**（circulation form）：
$$\boxed{\oint_{\partial D} P\,dx + Q\,dy = \iint_D \left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right) dA}$$

该形式描述的是"边界上的环量 = 内部旋度（标量旋度）的累积"。但物理学中经常出现的另一个量是**通量**（flux）——向量场在法向上的积分。下面给出Green公式的等价形式——**通量形式**（flux form）。

### 4.2 通量形式的陈述

**定理 14.33（Green公式——通量/散度形式 / Flux Form / 2D Divergence Theorem）**：设 $D \subset \mathbb{R}^2$ 为有界闭区域，其边界 $\partial D$ 由有限条分段光滑的简单闭曲线组成，且取**正向**（行走时区域在左侧）。若向量场 $\mathbf{F}(x,y) = (P(x,y),\; Q(x,y))$ 在包含 $D$ 的开区域上具有连续偏导数，则：
$$\boxed{\oint_{\partial D} \mathbf{F}\cdot\mathbf{n}\,ds = \iint_D \nabla\cdot\mathbf{F}\,dA}$$

其中 $\mathbf{n}$ 是 $\partial D$ 上的**外侧单位法向量**，$ds$ 是弧长微元。

### 4.3 与环量形式的等价性证明

**证明思路**：将通量形式通过向量旋转转化为环量形式。

设曲线 $\partial D$ 的单位切向量为 $\boldsymbol{\tau} = (\tau_x, \tau_y)$，单位法向量为 $\mathbf{n} = (n_x, n_y)$，取正向（逆时针）时，法向量指向区域外部。几何关系为：
$$\begin{cases}
n_x = \tau_y \\
n_y = -\tau_x
\end{cases}$$

直观理解：将切向量逆时针旋转 $90^\circ$ 得到外法向量。

由弧长微元 $ds$ 与投影微元的关系：
$$\tau_x\,ds = dx,\quad \tau_y\,ds = dy$$

因此法向通量的被积表达式可作如下变换：
$$\begin{aligned}
\mathbf{F}\cdot\mathbf{n}\,ds &= (P n_x + Q n_y)\,ds \\
&= (P\tau_y - Q\tau_x)\,ds \quad (\text{代入 }n_x=\tau_y,\; n_y=-\tau_x) \\
&= P\,dy - Q\,dx
\end{aligned}$$

于是：
$$\oint_{\partial D} \mathbf{F}\cdot\mathbf{n}\,ds = \oint_{\partial D} P\,dy - Q\,dx = \oint_{\partial D} (-Q)\,dx + P\,dy$$

对最后一式应用Green公式的环量形式（取 $\widetilde{P} = -Q$，$\widetilde{Q} = P$）：
$$\begin{aligned}
\oint_{\partial D} (-Q)\,dx + P\,dy &= \iint_D \left( \frac{\partial P}{\partial x} - \frac{\partial(-Q)}{\partial y} \right) dA \\
&= \iint_D \left( \frac{\partial P}{\partial x} + \frac{\partial Q}{\partial y} \right) dA \\
&= \iint_D \nabla\cdot\mathbf{F}\,dA \quad \blacksquare
\end{aligned}$$

### 4.4 两种形式的对比

| 形式 | 被积式 | 边界 | 右端 | 物理含义 |
|:-----|:-------|:-----|:-----|:---------|
| 环量形式 | $P\,dx + Q\,dy = \mathbf{F}\cdot\boldsymbol{\tau}\,ds$ | 切向分量 | $\iint (\partial Q/\partial x - \partial P/\partial y)\,dA$ | 边界环量 = 内部旋度累积 |
| 散度形式 | $\mathbf{F}\cdot\mathbf{n}\,ds$ | 法向分量 | $\iint \nabla\cdot\mathbf{F}\,dA$ | 边界通量 = 内部散度累积 |

**核心统一视角**：两者都是"边界积分 = 内部导数累积"这一微积分基本定理思想在二维的体现——环量形式对应旋度，散度形式对应散度。从环量形式出发，将向量场旋转 $90^\circ$ 即可得到散度形式。

### 4.5 与三维Gauss公式的类比

注意到二维散度定理与三维Gauss公式（定理14.14，文件06）具有完全相同的结构：

| 公式 | 维度 | 边界积分 | 区域积分 |
|:-----|:-----|:---------|:---------|
| 二维散度定理 | 2D | $\displaystyle\oint_{\partial D} \mathbf{F}\cdot\mathbf{n}\,ds$ | $\displaystyle\iint_D \nabla\cdot\mathbf{F}\,dA$ |
| Gauss公式 | 3D | $\displaystyle\iint_{\partial\Omega} \mathbf{F}\cdot\mathbf{n}\,dS$ | $\displaystyle\iiint_\Omega \nabla\cdot\mathbf{F}\,dV$ |

二维散度定理正是三维Gauss公式在二维空间中的退化形式——将 $\mathbb{R}^3$ 中不依赖于 $z$ 的向量场投影到 $xy$ 平面即得。

---

## 5. 调和函数的法向导数通量性质 (Type A-4)

### 5.1 法向导数线积分与Laplacian的关系

**定理 14.34（调和函数的法向导数通量性质）**：设 $u(x,y)$ 在区域 $\Omega \subset \mathbb{R}^2$ 上具有二阶连续偏导数。$D \subset \Omega$ 是有界闭区域，其边界 $\partial D$ 是分段光滑简单闭曲线。则：
$$\boxed{\oint_{\partial D} \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_D \nabla^2 u\,dA}$$

其中 $\mathbf{n}$ 是 $\partial D$ 的外侧单位法向量。

**证明**：由方向导数定义（定义14.32），$\partial u/\partial \mathbf{n} = \nabla u\cdot\mathbf{n}$。以 $\mathbf{F} = \nabla u$ 代入二维散度定理（定理14.33）：
$$\oint_{\partial D} \frac{\partial u}{\partial \mathbf{n}}\,ds = \oint_{\partial D} \nabla u\cdot\mathbf{n}\,ds = \iint_D \nabla\cdot(\nabla u)\,dA = \iint_D \nabla^2 u\,dA \quad \blacksquare$$

这个简洁的证明直接利用了方向导数定义和二维散度定理，将法向导数的线积分与Laplacian的面积分联系起来。它也可以绕过二维散度定理，通过Green公式的环量形式配合几何关系得到相同结论（见下面的替代推导）。

**替代推导**（使用环量形式）：利用3.3节中的几何关系，将 $\partial u/\partial \mathbf{n}\,ds$ 展开：
$$\begin{aligned}
\frac{\partial u}{\partial \mathbf{n}}\,ds &= \left( \frac{\partial u}{\partial x}n_x + \frac{\partial u}{\partial y}n_y \right) ds \\
&= \left( \frac{\partial u}{\partial x}\tau_y - \frac{\partial u}{\partial y}\tau_x \right) ds \\
&= \frac{\partial u}{\partial x}\,dy - \frac{\partial u}{\partial y}\,dx
\end{aligned}$$

代入Green公式环量形式，取 $P = -\partial u/\partial y$，$Q = \partial u/\partial x$：
$$\begin{aligned}
\oint_{\partial D} \frac{\partial u}{\partial \mathbf{n}}\,ds &= \oint_{\partial D} \frac{\partial u}{\partial x}\,dy - \frac{\partial u}{\partial y}\,dx \\
&= \oint_{\partial D} \left(-\frac{\partial u}{\partial y}\right) dx + \left(\frac{\partial u}{\partial x}\right) dy \\
&= \iint_D \left[ \frac{\partial}{\partial x}\left(\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial y}\left(-\frac{\partial u}{\partial y}\right) \right] dA \\
&= \iint_D \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) dA = \iint_D \nabla^2 u\,dA
\end{aligned}$$

两种推导结果一致。通量形式更为直接，环量形式则展示了更深层的几何变换。

### 5.2 正向命题：调和函数的法向导数通量为零

**正向命题**：若 $u$ 是区域 $\Omega$ 上的调和函数，则对 $\Omega$ 内任意光滑封闭曲线 $C$，有：
$$\boxed{\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0}$$

**证明**：由定理14.34，对以 $C$ 为边界的区域 $D$：
$$\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_D \nabla^2 u\,dA$$
由于 $u$ 是调和函数，$\nabla^2 u = 0$ 在 $D$ 上处处成立，因此右端积分为零。$\blacksquare$

### 5.3 反向命题：法向导数通量恒为零蕴含调和性

**反向命题**：设 $u$ 在 $\Omega$ 上具有**二阶连续偏导数**，且对 $\Omega$ 内任意光滑封闭曲线 $C$ 均有 $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0$。则 $u$ 是 $\Omega$ 上的调和函数。

**证明**（反证法）：假设 $u$ 不是调和的，则存在点 $M_0 \in \Omega$ 使得 $\nabla^2 u(M_0) \neq 0$。不妨设 $\nabla^2 u(M_0) = \delta > 0$（若 $\delta < 0$，考虑 $-u$ 同理可证）。

由于 $\nabla^2 u$ 在 $\Omega$ 上连续（$u$ 具有二阶连续偏导数 $\Rightarrow$ $\nabla^2 u$ 连续），由连续函数的局部保号性，存在 $M_0$ 的一个邻域 $B(M_0, \varepsilon) \subset \Omega$，使得在所有 $M \in B(M_0, \varepsilon)$ 上都有：
$$\nabla^2 u(M) > \frac{\delta}{2} > 0$$

取以 $M_0$ 为圆心、半径为 $\varepsilon$ 的圆周 $C_\varepsilon$ 作为积分曲线，$D_\varepsilon$ 为所围圆盘。由定理14.34：
$$\oint_{C_\varepsilon} \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_{D_\varepsilon} \nabla^2 u\,dA > \iint_{D_\varepsilon} \frac{\delta}{2}\,dA = \frac{\delta}{2} \cdot \pi\varepsilon^2 > 0$$

这与"对任意光滑封闭曲线 $C$ 有 $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0$"的假设矛盾。因此 $\nabla^2 u \equiv 0$ 在 $\Omega$ 上处处成立，即 $u$ 是 $\Omega$ 上的调和函数。$\blacksquare$

**定理 14.34'（完整等价性）**：设 $u$ 在 $\Omega$ 上具有二阶连续偏导数，则：
$$u\text{ 是 }\Omega\text{ 上的调和函数} \quad\Longleftrightarrow\quad \text{对 }\Omega\text{ 内任意光滑封闭曲线 }C\text{ 有 } \oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0$$

正向由定理14.34直接推出，反向由上述反证法得到。$\blacksquare$

**关键条件说明**：反向证明依赖 $\nabla^2 u$ 的连续性。若 $\nabla^2 u$ 不连续，则不能使用局部保号性推出矛盾。这正是定理中"二阶连续偏导数"条件的必要性所在。

---

## 6. 基本调和函数 1/r (Type A-5)

### 6.1 三维径向调和函数 $1/r$

在三维空间中，最基本也是最重要的调和函数之一是 $1/r$（除原点外），其中 $r = \sqrt{x^2 + y^2 + z^2}$。它在静电学中代表点电荷的电势。

**引理 14.35（$1/r$ 的梯度）**：设 $r = \sqrt{x^2 + y^2 + z^2} > 0$，则：
$$\boxed{\nabla\left(\frac{1}{r}\right) = -\frac{(x,y,z)}{r^3}}$$

**证明**：$\displaystyle \frac{\partial r}{\partial x} = \frac{x}{r}$
$$\frac{\partial}{\partial x}\left(\frac{1}{r}\right) = \frac{\partial}{\partial x}(r^{-1}) = -r^{-2}\cdot\frac{\partial r}{\partial x} = -\frac{1}{r^2}\cdot\frac{x}{r} = -\frac{x}{r^3}$$

同理 $\displaystyle \frac{\partial}{\partial y}\left(\frac{1}{r}\right) = -\frac{y}{r^3}$，$\displaystyle \frac{\partial}{\partial z}\left(\frac{1}{r}\right) = -\frac{z}{r^3}$。合并即得梯度表达式。$\blacksquare$

### 6.2 验证 $1/r$ 是调和函数

**定理 14.36（$1/r$ 的调和性）**：函数 $u(x,y,z) = 1/r$ 在 $\mathbb{R}^3 \setminus \{(0,0,0)\}$ 上满足Laplace方程：
$$\boxed{\nabla^2\left(\frac{1}{r}\right) = 0 \quad (r > 0)}$$

**证明**：计算二阶偏导数。由 $\partial(1/r)/\partial x = -x/r^3$：
$$\begin{aligned}
\frac{\partial^2}{\partial x^2}\left(\frac{1}{r}\right) &= \frac{\partial}{\partial x}\left(-\frac{x}{r^3}\right) = -\frac{\partial}{\partial x}\big(x\cdot r^{-3}\big) \\
&= -\left( r^{-3} + x \cdot (-3)r^{-4}\cdot\frac{\partial r}{\partial x} \right) \quad (\text{乘积法则 + 链式法则}) \\
&= -\left( \frac{1}{r^3} - \frac{3x^2}{r^5} \right) = -\frac{1}{r^3} + \frac{3x^2}{r^5}
\end{aligned}$$

由对称性：
$$\frac{\partial^2}{\partial y^2}\left(\frac{1}{r}\right) = -\frac{1}{r^3} + \frac{3y^2}{r^5},\qquad
\frac{\partial^2}{\partial z^2}\left(\frac{1}{r}\right) = -\frac{1}{r^3} + \frac{3z^2}{r^5}$$

求和：
$$\begin{aligned}
\nabla^2\left(\frac{1}{r}\right) &= \left(-\frac{1}{r^3} + \frac{3x^2}{r^5}\right) + \left(-\frac{1}{r^3} + \frac{3y^2}{r^5}\right) + \left(-\frac{1}{r^3} + \frac{3z^2}{r^5}\right) \\
&= -\frac{3}{r^3} + \frac{3(x^2+y^2+z^2)}{r^5} \\
&= -\frac{3}{r^3} + \frac{3r^2}{r^5} = -\frac{3}{r^3} + \frac{3}{r^3} = 0 \quad \blacksquare
\end{aligned}$$

**重要说明**：$1/r$ 在原点 $(0,0,0)$ 处无定义（趋于无穷大），因此调和性仅适用于 $\mathbb{R}^3\setminus\{(0,0,0)\}$。原点称为**奇点**（singularity）。

### 6.3 结合Gauss公式——通量与半径无关

**例 9**：设 $\Sigma_R$ 是以原点为球心、$R$ 为半径的球面，取外侧法向量 $\mathbf{n}$。计算：
$$I_R = \iint_{\Sigma_R} \nabla\left(\frac{1}{r}\right)\cdot\mathbf{n}\,dS$$

**解**：

**Step 1：求被积函数在球面上的值**

球面上 $r = R$，外侧单位法向量 $\mathbf{n} = (x,y,z)/R$。由引理14.35：
$$\nabla(1/r) = -\frac{(x,y,z)}{r^3}$$

在球面上：
$$\nabla(1/r)\cdot\mathbf{n} = \left(-\frac{(x,y,z)}{R^3}\right)\cdot\left(\frac{(x,y,z)}{R}\right) = -\frac{x^2 + y^2 + z^2}{R^4} = -\frac{R^2}{R^4} = -\frac{1}{R^2}$$

**Step 2：计算通量**

球面面积为 $4\pi R^2$，因此：
$$\begin{aligned}
I_R &= \iint_{\Sigma_R} \left(-\frac{1}{R^2}\right) dS \\
&= -\frac{1}{R^2} \times (\text{球面面积}) = -\frac{1}{R^2} \times 4\pi R^2 = -4\pi
\end{aligned}$$

**结果**：$\boxed{I_R = -4\pi}$，**与 $R$ 无关**。

### 6.4 通量与半径无关的物理解释

为什么 $I_R$ 与 $R$ 无关？根本原因在于 $1/r$ 在 $\mathbb{R}^3\setminus\{0\}$ 上是调和的，而 $\Sigma_R$ 包围了原点处的奇点。

**分析**：
1. 在球体 $B_R = \{r \leq R\}$ 内部，$\nabla^2(1/r) = 0$ 在 $B_R\setminus\{0\}$ 上成立，但在原点处无定义。
2. 若取两个球面 $\Sigma_{R_1}$ 和 $\Sigma_{R_2}$（$R_1 > R_2 > 0$），考虑夹在两个球面之间的区域 $\Omega = \{R_2 \leq r \leq R_1\}$。在 $\Omega$ 上，$\nabla^2(1/r) = 0$ 处处成立，且 $\nabla(1/r)$ 在 $\Omega$ 上具有连续偏导数。
3. 对 $\Omega$ 应用Gauss公式（外侧为 $\Sigma_{R_1}$ 的外侧 + $\Sigma_{R_2}$ 的内侧）：
   $$\iint_{\Sigma_{R_1}} \nabla(1/r)\cdot\mathbf{n}\,dS + \iint_{\Sigma_{R_2}} \nabla(1/r)\cdot(-\mathbf{n})\,dS = \iiint_{\Omega} \nabla^2(1/r)\,dV = 0$$
   因此 $I_{R_1} = I_{R_2}$，即通量与半径无关。

4. 直接计算验证：$I_R = -4\pi$ 与 $R$ 无关，正是这一性质的体现。

### 6.5 物理类比——点电荷的静电场

在静电学中，位于原点的点电荷 $q$ 产生的电势为：
$$V(r) = \frac{1}{4\pi\varepsilon_0}\cdot\frac{q}{r}$$

电场强度 $\mathbf{E} = -\nabla V = -\dfrac{q}{4\pi\varepsilon_0}\nabla\left(\dfrac{1}{r}\right)$。穿过包围该点电荷的任意封闭曲面的电通量为：
$$\iint_{\Sigma} \mathbf{E}\cdot\mathbf{n}\,dS = -\frac{q}{4\pi\varepsilon_0} \iint_{\Sigma} \nabla\left(\frac{1}{r}\right)\cdot\mathbf{n}\,dS = \frac{q}{\varepsilon_0}$$

这正是Gauss定律——电通量只与曲面内部的电荷有关，与曲面的半径无关。

### 6.6 有势无源场的势函数是调和函数

**定理 14.37（有势无源场的势函数调和性）**：设 $\Omega \subset \mathbb{R}^3$ 是单连通区域，$\mathbf{F}$ 是 $\Omega$ 上的向量场，且同时满足：
1. **有势场**（保守场）：存在势函数 $U$ 使 $\mathbf{F} = \nabla U$（定义14.23）
2. **无源场**（无散场）：$\nabla\cdot\mathbf{F} = 0$ 在 $\Omega$ 上处处成立

则势函数 $U$ 是 $\Omega$ 上的调和函数。

**证明**：
$$\nabla^2 U = \nabla\cdot(\nabla U) = \nabla\cdot\mathbf{F} = 0 \quad \blacksquare$$

这是调和函数与有势场理论（文件08）之间的直接联系：有势无源场的势函数必为调和函数。

**例 10**：验证 $\mathbf{F} = (x,\; y,\; -2z)$ 在 $\mathbb{R}^3$ 上既是有势场又是无源场，求其势函数 $U$，并验证 $U$ 是调和函数。

**解**：

**Step 1：验证有势性**
$$\nabla\times\mathbf{F} = \left( \frac{\partial(-2z)}{\partial y} - \frac{\partial y}{\partial z},\; \frac{\partial x}{\partial z} - \frac{\partial(-2z)}{\partial x},\; \frac{\partial y}{\partial x} - \frac{\partial x}{\partial y} \right) = (0,0,0)$$
$\mathbb{R}^3$ 单连通，因此 $\mathbf{F}$ 是有势场。

**Step 2：验证无源性**
$$\nabla\cdot\mathbf{F} = \frac{\partial}{\partial x}(x) + \frac{\partial}{\partial y}(y) + \frac{\partial}{\partial z}(-2z) = 1 + 1 - 2 = 0$$

**Step 3：求势函数**（直接积分法，参见文件08第5.1节）
由 $\partial U/\partial x = x$：$U = \dfrac{x^2}{2} + \varphi(y,z)$
由 $\partial U/\partial y = y$：$\partial\varphi/\partial y = y$，$\varphi = \dfrac{y^2}{2} + \psi(z)$
由 $\partial U/\partial z = -2z$：$\psi'(z) = -2z$，$\psi(z) = -z^2$
$$U(x,y,z) = \frac{x^2}{2} + \frac{y^2}{2} - z^2$$

**Step 4：验证调和性**
$$\nabla^2 U = \frac{\partial^2}{\partial x^2}\left(\frac{x^2}{2}\right) + \frac{\partial^2}{\partial y^2}\left(\frac{y^2}{2}\right) + \frac{\partial^2}{\partial z^2}(-z^2) = 1 + 1 - 2 = 0$$
因此 $U$ 是 $\mathbb{R}^3$ 上的调和函数。$\checkmark$

---

## 7. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|:---------|:---------|
| 调和函数必须 $\nabla^2 u = 0$ 在全空间处处成立，包括奇点 | 调和性在定义域内讨论。$1/r$ 在 $r=0$ 处无定义，只需在 $\mathbb{R}^3\setminus\{0\}$ 上满足 $\nabla^2(1/r) = 0$ 即可称其为上的调和函数 |
| $\partial u/\partial \mathbf{n} = \nabla u\cdot\mathbf{n}$ 中的 $\mathbf{n}$ 可以是任意非零向量 | 方向导数定义要求 $\mathbf{n}$ 是**单位向量**。若 $\mathbf{n}$ 不是单位向量，$\nabla u\cdot\mathbf{n}$ 给出的是方向导数的 $\|\mathbf{n}\|$ 倍，而非方向导数本身 |
| Green公式的通量形式和环量形式是两个独立定理 | 两者等价，可以通过将被积向量场旋转 $90^\circ$ 相互转化。通量形式是环量形式的直接推论（取 $\widetilde{P} = -Q$，$\widetilde{Q} = P$） |
| $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0$ 对所有函数都成立 | 这仅对**调和函数**成立（或更一般地，当 $\nabla^2 u$ 在 $C$ 所围区域内积分为零时）。对一般函数，该积分等于 $\iint_D \nabla^2 u\,dA$ |
| $\nabla^2(1/r) = 0$ 的证明可以省略，因为"显然" | 该计算虽不复杂，但需要仔细运用链式法则和乘积法则。二阶偏导项需要三项相加，不可遗漏或计算错误 |
| 用反证法证明逆命题时，连续性条件可有可无 | 连续性至关重要。若 $\nabla^2 u$ 不连续，即使在某点 $\nabla^2 u(M_0) > 0$，也不一定存在一个邻域使 $\nabla^2 u$ 保持正号，反证法失效 |
| 有势场 $\mathbf{F} = \nabla U$ 且 $\nabla\cdot\mathbf{F} = 0$ 时，$U$ 一定是调和函数，反过来也对 | 反过来正确：若 $U$ 是调和函数，则 $\mathbf{F} = \nabla U$ 满足 $\nabla\cdot\mathbf{F} = 0$，但 $\mathbf{F}$ 不一定是有势场——$\nabla U$ 本身就是梯度的形式，自动满足有势性定义 |
| 二维散度定理与三维Gauss公式无关 | 二维散度定理是三维Gauss公式在二维空间中的自然退化。两者共享完全相同的结构：边界法向通量 = 内部散度的体积分 |

### 检查点

- [ ] 能否写出Laplace算子 $\nabla^2$ 的定义并准确计算几个函数的Laplacian？
- [ ] 能否给出调和函数的定义，并列举3-5个不同类型的调和函数示例？
- [ ] 能否独立验证 $u = e^x\sin y$、$u = x^3 - 3xy^2$、$u = x^2 + y^2 - 2z^2$ 是调和函数？
- [ ] 能否写出方向导数的定义 $\partial u/\partial \mathbf{n} = \nabla u\cdot\mathbf{n}$ 及其分量展开形式？
- [ ] 能否陈述Green公式的通量形式（二维散度定理）并完成其与环量形式的等价性证明？
- [ ] 能否推导 $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_D \nabla^2 u\,dA$（至少用一种方法）？
- [ ] 能否独立完成定理14.34'（正向与反向）的完整证明？
- [ ] 能否解释反向证明中连续性条件的作用？
- [ ] 能否计算 $\nabla(1/r)$ 和 $\nabla^2(1/r)$ 的完整表达式？
- [ ] 能否计算球面通量 $I_R = \iint_{\Sigma_R} \nabla(1/r)\cdot\mathbf{n}\,dS$ 并解释为何与 $R$ 无关？
- [ ] 能否说明有势无源场的势函数必为调和函数？

---

## 8. 练习题

### 基础巩固

**1.**（调和性判别）判断下列函数在指定区域上是否为调和函数，并写出 $\nabla^2 u$ 的完整计算过程：

**(1)** $u(x,y) = x^3 - 3xy^2$，在 $\mathbb{R}^2$ 上

**(2)** $u(x,y) = e^x \cos y$，在 $\mathbb{R}^2$ 上

**(3)** $u(x,y,z) = x^2 + y^2 - 2z^2$，在 $\mathbb{R}^3$ 上

**(4)** $u(x,y) = \ln(x^2 + y^2)$，在 $\mathbb{R}^2\setminus\{(0,0)\}$ 上

<details><summary>参考答案</summary>

**(1)** $\dfrac{\partial u}{\partial x} = 3x^2 - 3y^2$，$\dfrac{\partial^2 u}{\partial x^2} = 6x$
$\dfrac{\partial u}{\partial y} = -6xy$，$\dfrac{\partial^2 u}{\partial y^2} = -6x$
$\nabla^2 u = 6x + (-6x) = 0$，是调和函数。✓

**(2)** $\dfrac{\partial u}{\partial x} = e^x\cos y$，$\dfrac{\partial^2 u}{\partial x^2} = e^x\cos y$
$\dfrac{\partial u}{\partial y} = -e^x\sin y$，$\dfrac{\partial^2 u}{\partial y^2} = -e^x\cos y$
$\nabla^2 u = e^x\cos y + (-e^x\cos y) = 0$，是调和函数。✓

**(3)** $\dfrac{\partial^2 u}{\partial x^2} = 2$，$\dfrac{\partial^2 u}{\partial y^2} = 2$，$\dfrac{\partial^2 u}{\partial z^2} = -4$
$\nabla^2 u = 2 + 2 + (-4) = 0$，是调和函数。✓

**(4)** $\dfrac{\partial u}{\partial x} = \dfrac{2x}{x^2+y^2}$，$\dfrac{\partial^2 u}{\partial x^2} = \dfrac{2(y^2-x^2)}{(x^2+y^2)^2}$
$\dfrac{\partial u}{\partial y} = \dfrac{2y}{x^2+y^2}$，$\dfrac{\partial^2 u}{\partial y^2} = \dfrac{2(x^2-y^2)}{(x^2+y^2)^2}$
$\nabla^2 u = \dfrac{2(y^2-x^2)}{(x^2+y^2)^2} + \dfrac{2(x^2-y^2)}{(x^2+y^2)^2} = 0$，是调和函数。✓

</details>

---

**2.**（方向导数计算）计算下列函数在给定点沿指定方向的方向导数：

**(1)** $u(x,y) = x^2 + y^2$，点 $(1,2)$，方向 $\mathbf{n} = \left(\dfrac{3}{5},\; \dfrac{4}{5}\right)$

**(2)** $u(x,y) = xy$，点 $(2,1)$，方向 $\mathbf{n} = \left(\dfrac{\sqrt{2}}{2},\; \dfrac{\sqrt{2}}{2}\right)$

**(3)** $u(x,y,z) = xyz$，点 $(1,1,1)$，方向 $\mathbf{n} = \left(\dfrac{1}{\sqrt{3}},\; \dfrac{1}{\sqrt{3}},\; \dfrac{1}{\sqrt{3}}\right)$

<details><summary>参考答案</summary>

**(1)** $\nabla u = (2x,\; 2y)$，$\nabla u(1,2) = (2,\; 4)$
$\dfrac{\partial u}{\partial \mathbf{n}} = (2,4)\cdot\left(\dfrac{3}{5},\dfrac{4}{5}\right) = \dfrac{6}{5} + \dfrac{16}{5} = \dfrac{22}{5}$

**(2)** $\nabla u = (y,\; x)$，$\nabla u(2,1) = (1,\; 2)$
$\dfrac{\partial u}{\partial \mathbf{n}} = (1,2)\cdot\left(\dfrac{\sqrt{2}}{2},\dfrac{\sqrt{2}}{2}\right) = \dfrac{\sqrt{2}}{2} + \sqrt{2} = \dfrac{3\sqrt{2}}{2}$

**(3)** $\nabla u = (yz,\; xz,\; xy)$，$\nabla u(1,1,1) = (1,\;1,\;1)$
$\dfrac{\partial u}{\partial \mathbf{n}} = (1,1,1)\cdot\left(\dfrac{1}{\sqrt{3}},\dfrac{1}{\sqrt{3}},\dfrac{1}{\sqrt{3}}\right) = \dfrac{3}{\sqrt{3}} = \sqrt{3}$

</details>

---

**3.**（二维散度定理的应用）利用二维散度定理（Green公式的通量形式）将下列曲线积分转化为二重积分并计算：

$$\oint_C (x^2 + y^2)\mathbf{n}\,ds$$

其中 $C$ 是正向圆周 $x^2 + y^2 = 1$，$\mathbf{n}$ 是外侧单位法向量。（提示：$\mathbf{F} = (x^2+y^2, 0)$ 不是向量场——这里的被积表达式是向量值，需要取一个分量。实际上应先将 $\oint_C \mathbf{F}\cdot\mathbf{n}\,ds$ 转化为 $\iint_D \nabla\cdot\mathbf{F}\,dA$ 的形式求解。）

<details><summary>参考答案</summary>

实际上本题应进一步明确为标量积分的计算。通量形式适用于 $\oint_C \mathbf{F}\cdot\mathbf{n}\,ds$。例如取 $\mathbf{F} = (x^2, y^2)$：
$$\nabla\cdot\mathbf{F} = 2x + 2y$$
$$\oint_C (x^2, y^2)\cdot\mathbf{n}\,ds = \iint_D (2x + 2y)\,dA = \int_0^{2\pi}\int_0^1 (2r\cos\theta + 2r\sin\theta) r\,dr\,d\theta = 0$$

（由对称性，$\cos\theta$ 和 $\sin\theta$ 在一个完整周期上的积分为零。）

</details>

---

### 迁移应用

**4.**（调和多项式构造）求所有形如 $u(x,y) = ax^2 + bxy + cy^2 + dx + ey + f$ 且满足Laplace方程 $\nabla^2 u = 0$ 的调和多项式（确定系数 $a,b,c,d,e,f$ 应满足的条件）。在此基础上给出一个包含混合项（$xy$项）且非常数调和的二次多项式示例。

<details><summary>参考答案</summary>

计算Laplacian：
$$\frac{\partial^2 u}{\partial x^2} = 2a,\quad \frac{\partial^2 u}{\partial y^2} = 2c$$
$$\nabla^2 u = 2a + 2c = 2(a + c)$$

$\nabla^2 u = 0$ 当且仅当 $a + c = 0$。系数 $b, d, e, f$ 不受限制。

一个含混合项且非常数调和的示例：取 $a = 0, b = 1, c = 0, d = 0, e = 0, f = 0$：
$$u(x,y) = xy$$
验证：$\partial^2(xy)/\partial x^2 = 0$，$\partial^2(xy)/\partial y^2 = 0$，$\nabla^2(xy) = 0$。✓

另一个示例：$u(x,y) = 3xy + x - y$（也是调和的）。

</details>

---

**5.**（法向导数通量的正向应用）设 $u(x,y) = x^2 - y^2$，$C$ 是正向单位圆周 $x^2 + y^2 = 1$。计算：
$$I = \oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds$$
（直接使用定理14.34的结论，无需做积分计算）

<details><summary>参考答案</summary>

$u(x,y) = x^2 - y^2$ 是调和函数（$\nabla^2 u = 2 - 2 = 0$）。由定理14.34'的正向命题：
$$\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_D \nabla^2 u\,dA = \iint_D 0\,dA = 0$$

直接验证：$\partial u/\partial \mathbf{n} = \nabla u\cdot\mathbf{n} = (2x, -2y)\cdot(x,y) = 2x^2 - 2y^2$（单位圆上 $\mathbf{n} = (x,y)$）。
$$\oint_C (2x^2 - 2y^2)\,ds = \int_0^{2\pi} (2\cos^2 t - 2\sin^2 t)\,dt = \int_0^{2\pi} 2\cos 2t\,dt = 0 \quad \checkmark$$

</details>

---

**6.**（法向导数通量的反向命题应用）设 $u(x,y)$ 在 $\mathbb{R}^2$ 上具有二阶连续偏导数。已知对任意光滑封闭曲线 $C$ 有 $\displaystyle\oint_C \frac{\partial u}{\partial \mathbf{n}}\,ds = 0$，且 $\nabla^2 u$ 在 $\mathbb{R}^2$ 上连续。证明 $u$ 是调和函数。

<details><summary>参考答案</summary>

这正是定理14.34'的反向命题。用反证法：

假设 $u$ 不是调和的，则存在点 $M_0$ 使 $\nabla^2 u(M_0) \neq 0$。不妨设 $\nabla^2 u(M_0) = \delta > 0$（$\delta < 0$ 时同理）。

由 $\nabla^2 u$ 的连续性，存在 $\varepsilon > 0$ 使在圆盘 $D_\varepsilon: (x-x_0)^2 + (y-y_0)^2 \leq \varepsilon^2$ 上，$\nabla^2 u > \delta/2 > 0$。

取 $C_\varepsilon = \partial D_\varepsilon$ 为正向圆周，由定理14.34：
$$\oint_{C_\varepsilon} \frac{\partial u}{\partial \mathbf{n}}\,ds = \iint_{D_\varepsilon} \nabla^2 u\,dA > \iint_{D_\varepsilon} \frac{\delta}{2}\,dA = \frac{\delta}{2}\cdot\pi\varepsilon^2 > 0$$

与假设矛盾。因此 $\nabla^2 u \equiv 0$，$u$ 是调和函数。$\blacksquare$

</details>

---

**7.**（$1/r$ 的Gauss公式综合应用）设 $r = \sqrt{x^2 + y^2 + z^2}$，$\Sigma_R$ 是以原点为球心、$R$ 为半径的球面（取外侧）。

**(1)** 计算 $I_R = \displaystyle\iint_{\Sigma_R} \nabla\left(\frac{1}{r}\right)\cdot\mathbf{n}\,dS$

**(2)** 设 $R_1 > R_2 > 0$，利用Gauss公式解释为什么 $I_{R_1} = I_{R_2}$

**(3)** 计算 $\displaystyle\iint_{\Sigma_R} \frac{\partial}{\partial \mathbf{n}}\left(\frac{1}{r}\right)\,dS$

<details><summary>参考答案</summary>

**(1)** 由例9的计算：
$$\nabla(1/r)\cdot\mathbf{n} = -\frac{1}{R^2}$$
$$I_R = \iint_{\Sigma_R} \left(-\frac{1}{R^2}\right) dS = -\frac{1}{R^2} \cdot 4\pi R^2 = -4\pi$$

**(2)** 考虑夹在 $\Sigma_{R_1}$ 和 $\Sigma_{R_2}$ 之间的区域 $\Omega = \{R_2 \leq r \leq R_1\}$。在 $\Omega$ 上，$\nabla^2(1/r) = 0$ 处处成立，且 $\nabla(1/r)$ 连续可微。对 $\Omega$ 应用Gauss公式（外侧为 $\Sigma_{R_1}$ 外侧 + $\Sigma_{R_2}$ 内侧）：
$$\iint_{\Sigma_{R_1}} \nabla(1/r)\cdot\mathbf{n}\,dS + \iint_{\Sigma_{R_2}} \nabla(1/r)\cdot(-\mathbf{n})\,dS = \iiint_{\Omega} \nabla^2(1/r)\,dV = 0$$
因此 $I_{R_1} - I_{R_2} = 0$，即 $I_{R_1} = I_{R_2}$。

**(3)** 由方向导数定义，$\partial(1/r)/\partial \mathbf{n} = \nabla(1/r)\cdot\mathbf{n}$，因此：
$$\iint_{\Sigma_R} \frac{\partial}{\partial \mathbf{n}}\left(\frac{1}{r}\right) dS = I_R = -4\pi$$

</details>

---

**8.**（思考题——有势无源场与调和函数）设 $\mathbf{F} = (y^2 - x^2,\; 2xy,\; 0)$ 是定义在 $\mathbb{R}^3$ 上的向量场。

**(1)** 证明 $\mathbf{F}$ 是有势场，并求其势函数 $U$。

**(2)** 验证 $\mathbf{F}$ 是无源场（$\nabla\cdot\mathbf{F} = 0$）。

**(3)** 直接验证 $U$ 是调和函数。

<details><summary>参考答案</summary>

**(1)** 旋度检验：$P = y^2 - x^2$，$Q = 2xy$，$R = 0$。
$$\nabla\times\mathbf{F} = (0-0,\; 0-0,\; 2y - 2y) = (0,0,0)$$
$\mathbb{R}^3$ 单连通，故是有势场。

求势函数（直接积分法）：
$\partial U/\partial x = y^2 - x^2$ → $U = y^2 x - x^3/3 + \varphi(y,z)$
$\partial U/\partial y = 2xy + \partial\varphi/\partial y = 2xy$ → $\partial\varphi/\partial y = 0$ → $\varphi = \psi(z)$
$\partial U/\partial z = \psi'(z) = 0$ → $\psi(z) = C$
$$U(x,y,z) = xy^2 - \frac{x^3}{3} + C$$

**(2)** 验证无源性：
$$\nabla\cdot\mathbf{F} = \frac{\partial}{\partial x}(y^2-x^2) + \frac{\partial}{\partial y}(2xy) + \frac{\partial}{\partial z}(0) = -2x + 2x + 0 = 0$$

**(3)** 验证调和性：
$$\frac{\partial^2 U}{\partial x^2} = \frac{\partial}{\partial x}(y^2 - x^2) = -2x$$
$$\frac{\partial^2 U}{\partial y^2} = \frac{\partial}{\partial y}(2xy) = 2x$$
$$\frac{\partial^2 U}{\partial z^2} = 0$$
$$\nabla^2 U = -2x + 2x + 0 = 0$$
因此 $U$ 是调和函数。✓

</details>
