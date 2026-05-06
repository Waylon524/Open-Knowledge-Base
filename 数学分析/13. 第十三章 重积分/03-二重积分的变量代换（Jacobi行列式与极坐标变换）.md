# 03. 二重积分的变量代换（Jacobi行列式与极坐标变换）

> 所属章节：第十三章 重积分  |  文件序号：03  |  难度：进阶
> 常见混淆点：极坐标变换下面积元素是 $dxdy = r\,dr\,d\theta$ 而非 $dr\,d\theta$（漏掉因子 $r$ 是最常见错误）；变量代换公式中 Jacobi 行列式必须取绝对值 $|J|$（负的行列式对应定向反转，面积伸缩因子取正）

## 1. 学习目标与先修前置

### 学习目标
- 掌握 Jacobi 行列式的定义 $\displaystyle J = \frac{\partial(x,y)}{\partial(u,v)}$，能熟练计算 2×2 行列式展开并取绝对值
- 理解一般变量代换公式 $\displaystyle \iint_D f(x,y)\,dxdy = \iint_{D'} f(x(u,v), y(u,v))\,|J(u,v)|\,du\,dv$ 及其几何含义（面积元的伸缩因子）
- 掌握极坐标变换 $x = r\cos\theta,\ y = r\sin\theta$ 下 $dxdy = r\,dr\,d\theta$ 的两种推导方法（几何法与 Jacobi 法）
- 能将圆形、扇形区域上的二重积分通过极坐标变换化为累次积分并完成计算
- 能将平行四边形区域上的二重积分通过线性变量代换化为矩形区域上的积分并完成计算
- 能与一维换元公式 $\int_a^b f(x)\,dx = \int_\alpha^\beta f(\psi(t))\,\psi'(t)\,dt$ 进行类比，理解 $|J|$ 与 $|\psi'(t)|$ 的对应关系

### 先修知识
- 文件01（第十三章）：二重积分的定义、Riemann 和、面积元素 $dxdy$ 的概念
- 文件02（第十三章）：X-型/Y-型区域的不等式组描述、累次积分的计算
- 文件03（第七章）：一维定积分的换元公式（定理7.7 和 定理7.8）——用于类比理解变量代换
- 文件01（第十二章）：偏导数的定义与计算方法（$\frac{\partial f}{\partial x}$、$\frac{\partial f}{\partial y}$ 的计算）
- 读者应熟悉 2×2 行列式的展开：$\begin{vmatrix} a & b \\ c & d \end{vmatrix} = ad - bc$

---

## 2. 背景与应用场景

### 2.1 从一维换元到二维换元

在第七章中，定积分的换元公式（定理7.8）告诉我们：

$$\int_a^b f(x)\,dx = \int_\alpha^\beta f(\psi(t))\,\psi'(t)\,dt,\quad x = \psi(t)$$

其中 $\psi'(t)$ 的作用是**补偿变量替换导致的长度伸缩**。当 $t$ 变化 $dt$ 时，$x$ 的变化量为 $dx = \psi'(t)\,dt$，因此积分区间被"拉伸"或"压缩"了 $\psi'(t)$ 倍。

在二维情形中，当我们对二重积分 $\iint_D f(x,y)\,dxdy$ 做变量代换 $(x,y) = T(u,v)$ 时，面积元 $dxdy$ 会被拉伸或压缩。但这个"拉伸"不再是一个简单的导数因子，而是一个**行列式**——因为二维的伸缩涉及两个方向（$u$ 方向和 $v$ 方向）的联合效应。

### 2.2 极坐标变换的直观动机

许多区域在直角坐标系下描述非常复杂——例如圆盘 $x^2 + y^2 \leq R^2$ 在 X-型下需要分段：$x \in [-R, R]$，$y \in [-\sqrt{R^2-x^2}, \sqrt{R^2-x^2}]$，内层积分限含根号，计算繁琐。但在极坐标系下，同一个圆盘简化为 $0 \leq r \leq R$，$0 \leq \theta \leq 2\pi$，是一个完美的矩形区域。

同样，被积函数 $e^{-(x^2+y^2)}$ 在直角坐标系下几乎无法积分（它的原函数不是初等函数），但在极坐标系下变为 $e^{-r^2}$，配合 $r\,dr$ 可以轻松积分。

### 2.3 线性变量代换的几何动机

对于平行四边形区域，若通过线性变换将其"拉直"为矩形，二重积分就可以用 Fubini 定理直接计算。这本质上是把坐标轴方向调整到和平行四边形的边平行的方向，从而让积分区域变成"坐标轴对齐"的矩形。

---

## 3. 核心概念与符号约定

### 3.1 一维换元法回顾——从类比出发

**定理回顾（一维第二换元法，定理7.8）**：设 $x = \psi(t)$ 在 $[\alpha, \beta]$ 上可导且 $\psi'(t)$ 连续，$\psi(\alpha)=a$，$\psi(\beta)=b$，则：

$$\int_a^b f(x)\,dx = \int_\alpha^\beta f(\psi(t))\,\psi'(t)\,dt$$

**关键解读**：换元函数 $\psi$ 将 $t$ 轴上的区间 $[\alpha, \beta]$ 映射到 $x$ 轴上的区间 $[a, b]$。在 $t$ 轴上的一个微小长度 $dt$ 对应 $x$ 轴上的长度 $dx = \psi'(t)\,dt$。因此，$x$ 轴上的积分转换到 $t$ 轴上时，**长度元的变换系数**是 $|\psi'(t)|$（取绝对值是因为长度总是正的）。

二维的推广思路：在二重积分中，$uv$ 平面上的一个小矩形区域（面积 $du\,dv$）经过变换 $(x,y) = T(u,v)$ 映射到 $xy$ 平面上的一小片区域（面积 $dxdy$）。两个面积之间的比例系数就是 **Jacobi 行列式的绝对值** $|J(u,v)|$。

### 3.2 Jacobi 行列式的定义（Type A-3）

**定义 13.9（Jacobi 行列式 / Jacobian Determinant）**：设 $(x,y) = T(u,v)$ 是从 $uv$ 平面到 $xy$ 平面的可微变换，即：

$$x = x(u,v),\quad y = y(u,v)$$

则 $T$ 的 **Jacobi 行列式**（简称 Jacobi 或 Jacobian）定义为：

$$J(u,v) = \frac{\partial(x,y)}{\partial(u,v)} = \begin{vmatrix}
\displaystyle\frac{\partial x}{\partial u} & \displaystyle\frac{\partial x}{\partial v} \\[8pt]
\displaystyle\frac{\partial y}{\partial u} & \displaystyle\frac{\partial y}{\partial v}
\end{vmatrix}$$

展开这个 2×2 行列式：

$$J(u,v) = \frac{\partial x}{\partial u} \cdot \frac{\partial y}{\partial v} - \frac{\partial x}{\partial v} \cdot \frac{\partial y}{\partial u}$$

**计算 Jacobi 行列式的三步流程**：

1. **计算四个偏导数**：$\frac{\partial x}{\partial u}$、$\frac{\partial x}{\partial v}$、$\frac{\partial y}{\partial u}$、$\frac{\partial y}{\partial v}$
2. **排列为 2×2 矩阵**：第一行是 $x$ 对 $u,v$ 的偏导，第二行是 $y$ 对 $u,v$ 的偏导
3. **展开行列式**：$(\frac{\partial x}{\partial u})(\frac{\partial y}{\partial v}) - (\frac{\partial x}{\partial v})(\frac{\partial y}{\partial u})$
4. **取绝对值**得 $|J(u,v)|$

**记号说明**：$\displaystyle\frac{\partial(x,y)}{\partial(u,v)}$ 是 Jacobi 行列式的标准记号。若变换是 $(u,v) \to (x,y)$ 则记为 $\frac{\partial(x,y)}{\partial(u,v)}$；逆变换的 Jacobi 记为 $\frac{\partial(u,v)}{\partial(x,y)}$，且有关系 $\frac{\partial(x,y)}{\partial(u,v)} \cdot \frac{\partial(u,v)}{\partial(x,y)} = 1$。

**例**：对线性变换 $x = 2u + v,\ y = u - 3v$，计算 $J$。

解：
- $\frac{\partial x}{\partial u} = 2$，$\frac{\partial x}{\partial v} = 1$
- $\frac{\partial y}{\partial u} = 1$，$\frac{\partial y}{\partial v} = -3$

$$J = \begin{vmatrix} 2 & 1 \\ 1 & -3 \end{vmatrix} = 2 \times (-3) - 1 \times 1 = -6 - 1 = -7,\quad |J| = 7$$

### 3.3 极坐标变换的定义（Type A-1）

**极坐标变换**（polar coordinate transformation）是变量代换中最重要的一类特例：

$$\begin{cases}
x = r\cos\theta \\[2pt]
y = r\sin\theta
\end{cases}\quad (r \geq 0,\ 0 \leq \theta \leq 2\pi)$$

其中 $r$ 表示点到原点的距离，$\theta$ 表示从 $x$ 轴正方向逆时针旋转的角度。

**极坐标变换下区域描述的规则**：

| 直角坐标区域 | 极坐标描述 |
|:---|:---|
| 圆 $x^2 + y^2 \leq R^2$ | $0 \leq r \leq R$，$0 \leq \theta \leq 2\pi$ |
| 第一象限四分之一圆（$x \geq 0, y \geq 0$） | $0 \leq r \leq R$，$0 \leq \theta \leq \frac{\pi}{2}$ |
| 上半圆（$y \geq 0$） | $0 \leq r \leq R$，$0 \leq \theta \leq \pi$ |
| 扇形（从 $\theta=\alpha$ 到 $\theta=\beta$） | $0 \leq r \leq R$，$\alpha \leq \theta \leq \beta$ |
| 环形 $a^2 \leq x^2 + y^2 \leq b^2$ | $a \leq r \leq b$，$0 \leq \theta \leq 2\pi$ |

**极坐标变换的 Jacobi 行列式计算**：

$$\begin{aligned}
\frac{\partial x}{\partial r} &= \cos\theta, &\quad \frac{\partial x}{\partial \theta} &= -r\sin\theta \\[4pt]
\frac{\partial y}{\partial r} &= \sin\theta, &\quad \frac{\partial y}{\partial \theta} &= r\cos\theta
\end{aligned}$$

$$J = \frac{\partial(x,y)}{\partial(r,\theta)} = \begin{vmatrix}
\cos\theta & -r\sin\theta \\
\sin\theta & r\cos\theta
\end{vmatrix} = \cos\theta \cdot r\cos\theta - (-r\sin\theta) \cdot \sin\theta = r(\cos^2\theta + \sin^2\theta) = r$$

因此 **Jacobi 行列式的绝对值**为 $|J| = r$（因为 $r \geq 0$），面积元素满足：

$$\boxed{dxdy = r\,dr\,d\theta}$$

### 3.4 一般变量代换公式（Type A-2）

**定理 13.10（二重积分的变量代换公式 / Change of Variables for Double Integrals）**：设 $T: (u,v) \mapsto (x,y)$ 是从 $uv$ 平面上的有界闭区域 $D'$ 到 $xy$ 平面上的有界闭区域 $D$ 的一一对应、连续可微的变换，且 Jacobi 行列式 $J(u,v) = \frac{\partial(x,y)}{\partial(u,v)}$ 在 $D'$ 上恒不为零。则对在 $D$ 上连续的函数 $f(x,y)$，有：

$$\iint_D f(x,y)\,dxdy = \iint_{D'} f\big(x(u,v), y(u,v)\big) \cdot |J(u,v)| \,du\,dv$$

**与一维换元公式的对比**：

| 维度 | 换元公式 | 伸缩因子 |
|:---|:---|:---|
| 一维 | $\int_a^b f(x)dx = \int_\alpha^\beta f(\psi(t))\,\psi'(t)\,dt$ | $|\psi'(t)|$ |
| 二维 | $\iint_D f(x,y)dxdy = \iint_{D'} f(x(u,v), y(u,v))\,|J(u,v)|\,du\,dv$ | $|J(u,v)|$ |

**公式的几何直观**：在 $uv$ 平面上取一个小矩形，边长 $\Delta u$、$\Delta v$，面积为 $\Delta u \Delta v$。经过变换 $T$ 后，这个小矩形映射到 $xy$ 平面上近似为一个平行四边形，其面积约为 $|J(u,v)| \cdot \Delta u \Delta v$。因此 $|J(u,v)|$ 就是**面积伸缩因子**（area scaling factor）。

### 3.5 符号表

| 符号 | 含义 | 说明 |
|:---|:---|:---|
| $J(u,v)$ | Jacobi 行列式 $\partial(x,y)/\partial(u,v)$ | 行列式值（可正可负） |
| $|J(u,v)|$ | Jacobi 行列式的绝对值 | 面积伸缩因子，始终非负 |
| $\frac{\partial(x,y)}{\partial(u,v)}$ | Jacobi 行列式的标准记号 | 记号的分子分母分别对应像空间和原空间 |
| $r$ | 极径（点到原点的距离） | $r \geq 0$ |
| $\theta$ | 极角（从 $x$ 轴正方向逆时针） | $\theta \in [0, 2\pi]$ |
| $dxdy = r\,dr\,d\theta$ | 极坐标下的面积元素 | 通过 $|J| = r$ 得到 |
| $T$ | 变量代换映射 | $T(u,v) = (x(u,v), y(u,v))$ |

---

## 4. 原理与方法

### 4.1 Jacobi 行列式的几何含义（Type A-2 + Type A-3）

为什么面积伸缩因子是 Jacobi 行列式的绝对值？让我们从线性逼近的角度推导。

设在 $uv$ 平面上取一个微小矩形 $R_{ij}$，顶点为 $(u_i, v_j)$，边长 $\Delta u$、$\Delta v$。这个矩形在变换 $T$ 下的像 $T(R_{ij})$ 近似是一个平行四边形，其两条邻边由以下向量给出：

沿 $u$ 方向的边：
$$T(u_i+\Delta u, v_j) - T(u_i, v_j) \approx \left( \frac{\partial x}{\partial u}\Delta u,\ \frac{\partial y}{\partial u}\Delta u \right) = \left( \frac{\partial x}{\partial u},\ \frac{\partial y}{\partial u} \right) \Delta u$$

沿 $v$ 方向的边：
$$T(u_i, v_j+\Delta v) - T(u_i, v_j) \approx \left( \frac{\partial x}{\partial v}\Delta v,\ \frac{\partial y}{\partial v}\Delta v \right) = \left( \frac{\partial x}{\partial v},\ \frac{\partial y}{\partial v} \right) \Delta v$$

由这两个向量张成的平行四边形的面积为：

$$\text{Area} = \left| \det \begin{pmatrix}
\frac{\partial x}{\partial u}\Delta u & \frac{\partial x}{\partial v}\Delta v \\[4pt]
\frac{\partial y}{\partial u}\Delta u & \frac{\partial y}{\partial v}\Delta v
\end{pmatrix} \right| = \left| \det \begin{pmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\[4pt]
\frac{\partial y}{\partial u} & \frac{\partial y}{\partial v}
\end{pmatrix} \right| \cdot \Delta u \Delta v = |J(u,v)| \cdot \Delta u \Delta v$$

因此，$uv$ 平面上面积为 $\Delta u \Delta v$ 的小矩形，在 $xy$ 平面上的像面积约为 $|J(u,v)| \cdot \Delta u \Delta v$。当取极限 $\Delta u,\Delta v \to 0$ 时，这个近似成为精确的等式，从而 $d\sigma_{xy} = |J|\,d\sigma_{uv}$，即 $dxdy = |J|\,du\,dv$。

### 4.2 极坐标面积元：两种推导方法（Type A-1）

**方法一（几何法）**：

在极坐标系 $(r,\theta)$ 中，考虑由 $r$ 到 $r+dr$、$\theta$ 到 $\theta+d\theta$ 围成的一个"微扇形"区域。这个区域的形状近似为一个矩形：

- 径向边长：$dr$
- 弧向边长：$r\,d\theta$（半径为 $r$ 时，角度 $d\theta$ 对应的弧长为 $r\,d\theta$）

因此面积元为：

$$d\sigma = (dr) \cdot (r\,d\theta) = r\,dr\,d\theta$$

**方法二（Jacobi 法）**：

直接代入定理 13.10。对极坐标变换 $x = r\cos\theta$、$y = r\sin\theta$，已在 3.3 节中计算得到 $J = r$，$|J| = r$，所以 $dxdy = r\,dr\,d\theta$。

两种方法结果一致。几何法更直观，Jacobi 法更系统——对于更一般的曲线坐标（如椭圆坐标、抛物线坐标），几何法可能不再直观，而 Jacobi 法始终适用。

### 4.3 极坐标变换的完整流程（Type A-1）

**第 1 步：区域描述**——将 $D$ 用 $r,\theta$ 的不等式组描述为 $D'$

注意到 $D'$ 在 $(r,\theta)$ 平面上通常是一个矩形（或矩形的一部分），这是因为：
- $r$ 的范围由点到原点的距离决定：$r_{\min} \leq r \leq r_{\max}$
- $\theta$ 的范围由区域覆盖的角度决定：$\theta_{\min} \leq \theta \leq \theta_{\max}$

**第 2 步：被积函数转换**——将 $f(x,y)$ 写为 $r,\theta$ 的表达式：$f(r\cos\theta, r\sin\theta)$

**第 3 步：面积元素转换**——$dxdy = r\,dr\,d\theta$

**第 4 步：化为累次积分**——在 $D'$ 上写出累次积分：

$$\iint_D f(x,y)\,dxdy = \int_{\theta_{\min}}^{\theta_{\max}} \int_{r_{\min}(\theta)}^{r_{\max}(\theta)} f(r\cos\theta, r\sin\theta) \cdot r\,dr\,d\theta$$

**第 5 步：计算积分**——若 $r$ 的范围与 $\theta$ 无关（圆盘、环形等），则 $\theta$ 积分和 $r$ 积分可分离：

$$\iint_D f(x,y)\,dxdy = \int_{\theta_{\min}}^{\theta_{\max}} g(\theta)\,d\theta \cdot \int_{r_{\min}}^{r_{\max}} h(r)\,r\,dr$$

**关键技巧**：遇到 $re^{-r^2}$、$r \cdot \text{函数}$ 形式时，用换元法 $u = r^2$，$du = 2r\,dr$，即 $r\,dr = \frac{1}{2}du$。

### 4.4 一般变量代换的完整流程（Type A-2）

**第 1 步：选择代换**——根据区域的形状和被积函数的结构，选择合适的变量代换 $(x,y) = T(u,v)$

常见选择策略：
- 区域边界由直线 $ax+by = \text{常数}$ 构成 $\to$ 设 $u = ax+by$，$v = cx+dy$
- 区域是圆形/扇形 $\to$ 极坐标变换
- 被积函数含 $x^2+y^2$ $\to$ 极坐标变换

**第 2 步：映射区域**——将 $D$ 的边界曲线映射到 $uv$ 平面，得到新区域 $D'$ 的不等式组描述

当变换是线性时，直线映射为直线，求法为：
- 将 $x = x(u,v)$、$y = y(u,v)$ 代入原边界方程，化简得到 $u,v$ 的关系式

**第 3 步：计算 Jacobi 行列式**——按三步流程（3.2 节）计算 $J(u,v)$ 并取绝对值 $|J(u,v)|$

**第 4 步：转换被积函数**——将 $f(x,y)$ 写为 $f(x(u,v), y(u,v))$

**第 5 步：在新区域上积分**——在 $D'$ 上用 Fubini 定理计算 $\iint_{D'} f(x(u,v), y(u,v))\,|J(u,v)|\,du\,dv$

### 4.5 线性变换的特例简化

当变换 $T$ 是线性变换时：

$$\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} \begin{pmatrix} u \\ v \end{pmatrix} + \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$$

Jacobi 行列式就是系数矩阵的行列式（常数）：

$$J = \det \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$$

这是一个常数，可以提到积分号外面。**线性变换的特殊性**：区域内各点的面积伸缩因子处处相等。

---

## 5. 例题

### 例 1：极坐标变换——四分之一单位圆盘（Type A-1）

计算 $\displaystyle I = \iint_D y\,dxdy$，其中 $D$ 是由 $x^2 + y^2 \leq 1$、$x \geq 0$、$y \geq 0$ 所围成的第一象限四分之一单位圆盘。

**第 1 步：引入极坐标变换**

$$x = r\cos\theta,\quad y = r\sin\theta,\quad r \geq 0,\ 0 \leq \theta \leq 2\pi$$

**第 2 步：用 $r,\theta$ 描述区域 $D$**

$x^2 + y^2 \leq 1$ 变为 $r^2 \leq 1$，即 $0 \leq r \leq 1$。
$x \geq 0$ 要求 $\cos\theta \geq 0$，即 $\theta \in [-\frac{\pi}{2}, \frac{\pi}{2}]$。结合 $y \geq 0$ 要求 $\sin\theta \geq 0$，得 $\theta \in [0, \frac{\pi}{2}]$。

因此 $D$ 在 $(r,\theta)$ 平面上的像为：
$$D' = \big\{ (r,\theta) \mid 0 \leq r \leq 1,\ 0 \leq \theta \leq \tfrac{\pi}{2} \big\}$$

这是一个矩形区域。

**第 3 步：转换被积函数和面积元素**

被积函数：$y = r\sin\theta$
面积元素：$dxdy = r\,dr\,d\theta$

**第 4 步：写出极坐标下的累次积分**

$$\begin{aligned}
I &= \iint_D y\,dxdy = \iint_{D'} (r\sin\theta) \cdot r\,dr\,d\theta \\
  &= \int_{\theta=0}^{\pi/2} \int_{r=0}^{1} r^2 \sin\theta \,dr\,d\theta
\end{aligned}$$

**第 5 步：分离变量并计算**

被积函数 $r^2 \sin\theta$ 是 $r$ 的函数和 $\theta$ 的函数的乘积，且积分区域是矩形，因此积分可分离：

$$\begin{aligned}
I &= \left( \int_0^{\pi/2} \sin\theta\,d\theta \right) \cdot \left( \int_0^1 r^2\,dr \right) \\
  &= \left[ -\cos\theta \right]_0^{\pi/2} \cdot \left[ \frac{r^3}{3} \right]_0^1 \\
  &= \big( -\cos\frac{\pi}{2} - (-\cos 0) \big) \cdot \frac{1}{3} \\
  &= (0 + 1) \cdot \frac{1}{3} = \frac{1}{3}
\end{aligned}$$

因此 $\displaystyle \iint_D y\,dxdy = \frac{1}{3}$。

**验证——用直角坐标**（Fubini 定理，X-型区域）：

区域 $D$ 的 X-型描述：$D = \{(x,y) \mid 0 \leq x \leq 1,\ 0 \leq y \leq \sqrt{1-x^2}\}$

$$\begin{aligned}
I &= \int_0^1 dx \int_0^{\sqrt{1-x^2}} y\,dy \\
  &= \int_0^1 \left[ \frac{y^2}{2} \right]_0^{\sqrt{1-x^2}} dx \\
  &= \int_0^1 \frac{1-x^2}{2}\,dx = \frac{1}{2} \left[ x - \frac{x^3}{3} \right]_0^1 = \frac{1}{2} \cdot \frac{2}{3} = \frac{1}{3}
\end{aligned}$$

两种方法结果一致。$\checkmark$

---

### 例 2：线性变量代换——平行四边形区域（Type A-2 + Type A-3）

计算 $\displaystyle I = \iint_D (x + 2y)^2\,dxdy$，其中 $D$ 是由直线 $x+2y=2$、$x+2y=4$、$x-y=0$、$x-y=2$ 所围成的平行四边形区域。

**第 1 步：选择代换**

观察区域的边界线：
- $x + 2y = 2$ 和 $x + 2y = 4$ ——这是两条平行线
- $x - y = 0$ 和 $x - y = 2$ ——这是另外两条平行线

令 $u = x + 2y$，$v = x - y$。这样，$D$ 的边界在 $uv$ 平面上变为常数线，$D$ 映射为一个矩形。

**第 2 步：映射区域**

当 $(x,y)$ 在 $D$ 中时：
- $u = x+2y$ 的范围从 $2$ 到 $4$
- $v = x-y$ 的范围从 $0$ 到 $2$

因此 $D$ 在 $uv$ 平面上的像为：
$$D' = \big\{ (u,v) \mid 2 \leq u \leq 4,\ 0 \leq v \leq 2 \big\}$$

这是一个矩形。

**第 3 步：求逆变换并计算 Jacobi 行列式**

我们需要 $x$ 和 $y$ 关于 $u$、$v$ 的表达式。解方程组：

$$\begin{cases}
u = x + 2y \\
v = x - y
\end{cases}$$

从第二式得 $x = v + y$，代入第一式：
$$u = (v + y) + 2y = v + 3y \quad\Longrightarrow\quad y = \frac{u - v}{3}$$

代回 $x = v + y = v + \frac{u - v}{3} = \frac{3v + u - v}{3} = \frac{u + 2v}{3}$。

因此逆变换为：
$$x = \frac{u + 2v}{3},\quad y = \frac{u - v}{3}$$

计算 Jacobi 行列式（注意这里是从 $(u,v)$ 到 $(x,y)$，所以 $\frac{\partial(x,y)}{\partial(u,v)}$）：

$$\begin{aligned}
\frac{\partial x}{\partial u} &= \frac{1}{3}, &\quad \frac{\partial x}{\partial v} &= \frac{2}{3} \\[4pt]
\frac{\partial y}{\partial u} &= \frac{1}{3}, &\quad \frac{\partial y}{\partial v} &= -\frac{1}{3}
\end{aligned}$$

$$J(u,v) = \begin{vmatrix}
\frac{1}{3} & \frac{2}{3} \\[6pt]
\frac{1}{3} & -\frac{1}{3}
\end{vmatrix}
= \frac{1}{3} \cdot \left(-\frac{1}{3}\right) - \frac{2}{3} \cdot \frac{1}{3}
= -\frac{1}{9} - \frac{2}{9} = -\frac{1}{3}$$

$$|J(u,v)| = \frac{1}{3}$$

**第 4 步：转换被积函数**

$$(x + 2y)^2 = u^2$$

这里直接用 $u$ 表达，非常简洁——这正是选择 $u = x+2y$ 作为新变量的原因。

**第 5 步：在新区域上计算积分**

$$\begin{aligned}
I &= \iint_D (x+2y)^2\,dxdy = \iint_{D'} u^2 \cdot |J(u,v)| \,du\,dv \\
  &= \iint_{D'} u^2 \cdot \frac{1}{3}\,du\,dv \\
  &= \frac{1}{3} \int_{u=2}^{4} \int_{v=0}^{2} u^2 \,dv\,du
\end{aligned}$$

先对 $v$ 积分（被积函数 $u^2$ 与 $v$ 无关，$\int_0^2 dv = 2$）：

$$\begin{aligned}
I &= \frac{1}{3} \int_{2}^{4} u^2 \cdot 2\,du \\
  &= \frac{2}{3} \int_{2}^{4} u^2\,du = \frac{2}{3} \cdot \left[ \frac{u^3}{3} \right]_{2}^{4} \\
  &= \frac{2}{3} \cdot \left( \frac{64}{3} - \frac{8}{3} \right) = \frac{2}{3} \cdot \frac{56}{3} = \frac{112}{9}
\end{aligned}$$

因此 $\displaystyle \iint_D (x+2y)^2\,dxdy = \frac{112}{9}$。

**思考**：如果不使用变量代换，直接在直角坐标系下计算这个平行四边形区域上的二重积分，需要将区域分割为多个 X-型或 Y-型子区域，计算量远大于变量代换法。

---

### 例 3：极坐标变换——整圆盘与 Gauss 型积分（Type A-1）

计算 $\displaystyle I = \iint_D e^{-(x^2+y^2)}\,dxdy$，其中 $D$ 是由 $x^2 + y^2 \leq 4$ 所围成的圆盘（半径为 2）。

**第 1 步：引入极坐标变换**

$$x = r\cos\theta,\quad y = r\sin\theta$$

**第 2 步：用 $r,\theta$ 描述区域 $D$**

$x^2 + y^2 \leq 4$ 变为 $r^2 \leq 4$，即 $0 \leq r \leq 2$。$\theta$ 的取值范围覆盖整个圆周：$0 \leq \theta \leq 2\pi$。

$$D' = \big\{ (r,\theta) \mid 0 \leq r \leq 2,\ 0 \leq \theta \leq 2\pi \big\}$$

**第 3 步：转换被积函数和面积元素**

被积函数：$e^{-(x^2+y^2)} = e^{-r^2}$
面积元素：$dxdy = r\,dr\,d\theta$

**第 4 步：写出极坐标下的累次积分**

$$I = \int_{\theta=0}^{2\pi} \int_{r=0}^{2} e^{-r^2} \cdot r\,dr\,d\theta$$

**第 5 步：分离变量并计算**

被积函数 $e^{-r^2} \cdot r$ 与 $\theta$ 无关，积分区域为矩形，因此：

$$I = \left( \int_0^{2\pi} d\theta \right) \cdot \left( \int_0^2 r e^{-r^2}\,dr \right) = 2\pi \cdot \int_0^2 r e^{-r^2}\,dr$$

计算 $r$ 积分：作换元 $u = r^2$，则 $du = 2r\,dr$，$r\,dr = \frac{1}{2}du$。当 $r = 0$ 时 $u = 0$，$r = 2$ 时 $u = 4$。

$$\int_0^2 r e^{-r^2}\,dr = \int_0^4 e^{-u} \cdot \frac{1}{2}\,du = \frac{1}{2} \left[ -e^{-u} \right]_0^4 = \frac{1}{2}(-e^{-4} + 1) = \frac{1 - e^{-4}}{2}$$

因此：

$$I = 2\pi \cdot \frac{1 - e^{-4}}{2} = \pi (1 - e^{-4})$$

**重要注解**：函数 $e^{-x^2}$ 的原函数不是初等函数（无法用有限次基本初等函数组合表示），因此定积分 $\int_{-2}^{2} e^{-x^2}dx$ 在直角坐标系下无法直接通过 Newton-Leibniz 公式计算。但通过二重积分和极坐标变换，我们巧妙地求出了 $\iint_D e^{-(x^2+y^2)}dxdy$ 的值。这一技巧在概率论和数理统计中有着极其重要的应用——用于推导正态分布的概率密度函数的归一化常数 $\int_{-\infty}^{\infty} e^{-x^2}dx = \sqrt{\pi}$。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|:---|:---|
| 极坐标下面积元素就是 $dxdy = dr\,d\theta$（忘了 $r$） | 面积元素是 $dxdy = r\,dr\,d\theta$，因子 $r$ 来源于弧长 $r\,d\theta$ |
| 变量代换后 $dxdy = du\,dv$（忘了 Jacobi 行列式） | 必须乘以 $|J(u,v)|$：$dxdy = |J(u,v)|\,du\,dv$ |
| Jacobi 行列式 $J$ 直接就是伸缩因子 | 伸缩因子是 $|J|$（绝对值）。$J$ 可为负，代表定向反转，但面积不能为负 |
| 变量代换只须转换被积函数和面积元，不用管区域边界 | 区域 $D$ 的边界必须同步映射为 $D'$ 的边界。换元后在新区域 $D'$ 上积分 |
| 极坐标变换中 $r$ 的范围总是从 $0$ 开始 | 对于环形区域 $a \leq r \leq b$，$r$ 从 $a$ 开始。$r$ 的上下限取决于区域形状 |
| 线性变换的 Jacobi 行列式依赖于 $(u,v)$ 的位置 | 线性变换的 Jacobi 行列式是常数（系数矩阵的行列式），与位置无关 |
| 极坐标变换只适用于圆形区域 | 极坐标适用于任何在极坐标下描述更简便的区域，如心形线 $r = 1+\cos\theta$ 等 |

### 检查点

- [ ] 能否写出 2×2 Jacobi 行列式的定义公式（包括四个偏导数的排列规则）？
- [ ] 能否用几何法（扇形面积）和 Jacobi 法两种方式推导 $dxdy = r\,dr\,d\theta$？
- [ ] 能否完成极坐标变换的完整五步流程：区域描述、被积函数转换、面积元素转换、化累次积分、计算？
- [ ] 能否说出一般变量代换公式中 $|J(u,v)|$ 的几何意义？
- [ ] 对于线性变换 $u = ax+by$、$v = cx+dy$，能否说明为什么 $J(u,v)$ 是常数？
- [ ] 当极坐标下被积函数出现 $re^{-r^2}$ 时，能否用换元 $u = r^2$ 计算？
- [ ] 能否将定理 13.10（二维变量代换）与定理 7.8（一维换元）进行类比，说明 $|J|$ 与 $|\psi'|$ 的对应关系？
- [ ] 对于由直线围成的平行四边形区域，能否识别出合适的线性变量代换？

---

## 练习题

### 基础巩固

**1.**（Jacobi 行列式计算）对下列变换，计算 Jacobi 行列式 $J = \frac{\partial(x,y)}{\partial(u,v)}$ 及其绝对值。

**(a)** $x = 3u - v$，$y = u + 2v$

**(b)** $x = u^2 - v^2$，$y = 2uv$

**(c)** $x = e^u \cos v$，$y = e^u \sin v$

<details><summary>参考答案</summary>

**(a)** 
$$\frac{\partial x}{\partial u} = 3,\ \frac{\partial x}{\partial v} = -1,\ \frac{\partial y}{\partial u} = 1,\ \frac{\partial y}{\partial v} = 2$$
$$J = \begin{vmatrix} 3 & -1 \\ 1 & 2 \end{vmatrix} = 3 \times 2 - (-1) \times 1 = 6 + 1 = 7,\quad |J| = 7$$

**(b)**
$$\frac{\partial x}{\partial u} = 2u,\ \frac{\partial x}{\partial v} = -2v,\ \frac{\partial y}{\partial u} = 2v,\ \frac{\partial y}{\partial v} = 2u$$
$$J = \begin{vmatrix} 2u & -2v \\ 2v & 2u \end{vmatrix} = (2u)(2u) - (-2v)(2v) = 4u^2 + 4v^2 = 4(u^2+v^2)$$
$$|J| = 4(u^2+v^2)$$

**(c)**
$$\frac{\partial x}{\partial u} = e^u\cos v,\ \frac{\partial x}{\partial v} = -e^u\sin v,\ \frac{\partial y}{\partial u} = e^u\sin v,\ \frac{\partial y}{\partial v} = e^u\cos v$$
$$J = \begin{vmatrix} e^u\cos v & -e^u\sin v \\ e^u\sin v & e^u\cos v \end{vmatrix} = e^{2u}(\cos v \cdot \cos v - (-\sin v) \cdot \sin v) = e^{2u}(\cos^2v + \sin^2v) = e^{2u}$$
$$|J| = e^{2u}$$

</details>

---

**2.**（极坐标变换基础）利用极坐标变换计算下列二重积分：

**(a)** $\displaystyle I = \iint_D x^2\,dxdy$，其中 $D$ 是半径为 $R$ 的圆盘 $x^2 + y^2 \leq R^2$

**(b)** $\displaystyle I = \iint_D \sqrt{x^2 + y^2}\,dxdy$，其中 $D$ 是由 $x^2 + y^2 \leq 1$ 所围成的单位圆盘

**(c)** $\displaystyle I = \iint_D y\,dxdy$，其中 $D$ 是上半圆盘 $x^2 + y^2 \leq 1$，$y \geq 0$

<details><summary>参考答案</summary>

**(a)** 极坐标变换 $x = r\cos\theta$，$y = r\sin\theta$，$dxdy = r\,dr\,d\theta$。

区域 $D: 0 \leq r \leq R,\ 0 \leq \theta \leq 2\pi$。

被积函数 $x^2 = r^2\cos^2\theta$。

$$\begin{aligned}
I &= \int_0^{2\pi} \int_0^R r^2\cos^2\theta \cdot r\,dr\,d\theta \\
  &= \left( \int_0^{2\pi} \cos^2\theta\,d\theta \right) \left( \int_0^R r^3\,dr \right) \\
  &= \left( \int_0^{2\pi} \frac{1+\cos 2\theta}{2}\,d\theta \right) \cdot \left[ \frac{r^4}{4} \right]_0^R \\
  &= \frac{1}{2} \cdot 2\pi \cdot \frac{R^4}{4} = \frac{\pi R^4}{4}
\end{aligned}$$

**(b)** 极坐标变换，$x^2 + y^2 = r^2$，所以 $\sqrt{x^2 + y^2} = r$。

$$I = \int_0^{2\pi} \int_0^1 r \cdot r\,dr\,d\theta = 2\pi \int_0^1 r^2\,dr = 2\pi \cdot \frac{1}{3} = \frac{2\pi}{3}$$

**(c)** 极坐标变换，$y = r\sin\theta$。

上半圆盘：$0 \leq r \leq 1$，$0 \leq \theta \leq \pi$。

$$\begin{aligned}
I &= \int_0^{\pi} \int_0^1 (r\sin\theta) \cdot r\,dr\,d\theta \\
  &= \left( \int_0^{\pi} \sin\theta\,d\theta \right) \left( \int_0^1 r^2\,dr \right) \\
  &= \left[ -\cos\theta \right]_0^{\pi} \cdot \frac{1}{3} = (1 - (-1)) \cdot \frac{1}{3} = 2 \cdot \frac{1}{3} = \frac{2}{3}
\end{aligned}$$

</details>

---

### 迁移应用

**3.**（线性变量代换）计算 $\displaystyle I = \iint_D (x - y)^2 e^{x+y}\,dxdy$，其中 $D$ 是由 $x+y=1$、$x+y=3$、$x-y=-1$、$x-y=1$ 所围成的平行四边形区域。

<details><summary>参考答案</summary>

**第 1 步：选择代换。** 令 $u = x + y$，$v = x - y$。

**第 2 步：映射区域。**

$u$ 从 $1$ 到 $3$，$v$ 从 $-1$ 到 $1$。
$$D' = \big\{ (u,v) \mid 1 \leq u \leq 3,\ -1 \leq v \leq 1 \big\}$$

**第 3 步：求逆变换和 Jacobi 行列式。**

解方程组：
$$\begin{cases}
u = x + y \\
v = x - y
\end{cases}$$

两式相加：$u + v = 2x$，得 $x = \dfrac{u+v}{2}$。
两式相减：$u - v = 2y$，得 $y = \dfrac{u-v}{2}$。

计算偏导数：
$$\frac{\partial x}{\partial u} = \frac{1}{2},\ \frac{\partial x}{\partial v} = \frac{1}{2},\ \frac{\partial y}{\partial u} = \frac{1}{2},\ \frac{\partial y}{\partial v} = -\frac{1}{2}$$

$$J = \begin{vmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & -\frac{1}{2} \end{vmatrix}
= \frac{1}{2} \cdot \left(-\frac{1}{2}\right) - \frac{1}{2} \cdot \frac{1}{2}
= -\frac{1}{4} - \frac{1}{4} = -\frac{1}{2}$$

$$|J| = \frac{1}{2}$$

**第 4 步：转换被积函数。**

$(x-y)^2 e^{x+y} = v^2 e^u$

**第 5 步：在新区域上计算积分。**

$$\begin{aligned}
I &= \iint_{D'} v^2 e^u \cdot \frac{1}{2}\,du\,dv \\
  &= \frac{1}{2} \int_{u=1}^{3} \int_{v=-1}^{1} v^2 e^u \,dv\,du \\
  &= \frac{1}{2} \left( \int_{1}^{3} e^u\,du \right) \left( \int_{-1}^{1} v^2\,dv \right) \\
  &= \frac{1}{2} \left[ e^u \right]_{1}^{3} \cdot \left[ \frac{v^3}{3} \right]_{-1}^{1} \\
  &= \frac{1}{2} (e^3 - e) \cdot \left( \frac{1}{3} - \left(-\frac{1}{3}\right) \right) \\
  &= \frac{1}{2} (e^3 - e) \cdot \frac{2}{3} = \frac{e^3 - e}{3}
\end{aligned}$$

因此 $\displaystyle I = \frac{e^3 - e}{3}$。

</details>

---

**4.**（极坐标——环形区域）计算 $\displaystyle I = \iint_D \frac{1}{x^2 + y^2}\,dxdy$，其中 $D$ 是由 $x^2 + y^2 \geq 1$ 和 $x^2 + y^2 \leq e^2$ 所围成的环形区域（$e$ 是自然常数）。

<details><summary>参考答案</summary>

**第 1 步：极坐标变换。**

$x = r\cos\theta$，$y = r\sin\theta$，$dxdy = r\,dr\,d\theta$。

**第 2 步：区域描述。**

$x^2 + y^2 \geq 1$ 变为 $r \geq 1$；$x^2 + y^2 \leq e^2$ 变为 $r \leq e$。
$$\theta \in [0, 2\pi],\quad r \in [1, e]$$

**第 3 步：转换被积函数。**

$\dfrac{1}{x^2 + y^2} = \dfrac{1}{r^2}$

**第 4 步：写累次积分并计算。**

$$\begin{aligned}
I &= \int_0^{2\pi} \int_1^e \frac{1}{r^2} \cdot r\,dr\,d\theta \\
  &= \int_0^{2\pi} d\theta \cdot \int_1^e \frac{1}{r}\,dr \\
  &= 2\pi \cdot \left[ \ln r \right]_1^e \\
  &= 2\pi (\ln e - \ln 1) = 2\pi
\end{aligned}$$

因此 $\displaystyle I = 2\pi$。

**关键点**：这里 $r$ 从 $1$ 开始（不是从 $0$），因为区域是环形，内部被挖空。$r$ 的范围必须覆盖整个环形区域，即 $1 \leq r \leq e$。

</details>
