# 07. Stokes公式（斯托克斯公式）

> 所属章节：第十四章 曲线积分、曲面积分与场论  |  文件序号：07  |  难度：进阶
> 常见混淆点：混淆旋度 $\nabla\times\mathbf{F}$（向量）与散度 $\nabla\cdot\mathbf{F}$（标量）——旋度是向量场的"旋转强度"和方向，散度是"源/汇强度"，两者通过 Nabla 算子的不同运算（叉积 vs 点积）区分；忽略右手法则——曲面定向与边界曲线的方向必须通过右手法则协调，方向不一致时计算结果差一个负号

## 1. 学习目标与先修前置

### 学习目标
- 掌握三维旋度的完整定义 $\displaystyle \nabla\times\mathbf{F} = \left(\frac{\partial R}{\partial y}-\frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z}-\frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\right)$，理解行列式记忆法和 Nabla 叉积记法
- 理解旋度的物理含义——局部旋转的轴方向和强度，能够计算给定向量的旋度
- 掌握 Stokes 公式的完整陈述（有向投影形式和向量点积形式）、成立条件及两种形式的等价关系
- 理解右手法则：曲面定向与边界定向的对应关系，能在具体问题中应用
- 掌握 Stokes 公式退化到 Green 公式的完整推理——当曲面退化为 $xy$ 平面区域时，Stokes 公式给出 Green 公式
- 掌握利用 Stokes 公式将空间闭曲线上的第二型曲线积分转化为曲面积分的五步标准工作流

### 先修知识
- 文件03（第十四章）：Green 公式 $\oint_{\partial D} P\,dx+Q\,dy = \iint_D (\partial Q/\partial x-\partial P/\partial y)\,dA$ 及其二维旋度概念——Stokes 公式是 Green 公式从平面到空间曲面的推广
- 文件04（第十四章）：第二类曲面积分的通量定义（定义14.10）、有向投影微元 $dy\wedge dz, dz\wedge dx, dx\wedge dy$、侧反转变号（定理14.11）——Stokes 公式将曲线积分转化为曲面积分
- 文件05（第十四章）：第一类曲面积分的 $dS$ 公式 $dS = \sqrt{1+z_x^2+z_y^2}\,dx\,dy$——Stokes 公式的向量形式 $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ 通过第一类曲面积分计算
- 文件06（第十四章）：Nabla 算子 $\nabla = (\partial/\partial x,\;\partial/\partial y,\;\partial/\partial z)$ 的定义、散度 $\nabla\cdot\mathbf{F}$、Gauss 公式——三维向量微分的场论框架
- 文件02（第十四章）：第二型曲线积分的参数化计算方法、方向反转性质
- 文件12（第十二章）：偏导数计算能力——旋度计算的基础

---

## 2. 背景与应用场景

### 2.1 从二维 Green 公式到三维 Stokes 公式

在文件03（第十四章）中，我们学习了 Green 公式：
$$\oint_{\partial D} P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$$

这个公式建立了**平面闭曲线上的环量**与**内部旋度的面积分**之间的等式。Green 公式本质上是二维的——它要求积分曲线 $C$ 是平面上的闭曲线，所围区域 $D$ 是平面区域。

现在考虑三维空间中的情形：设 $\Sigma$ 是一张有向曲面（可能弯曲），$\partial\Sigma$ 是它的边界曲线（一条空间闭曲线）。能否将沿 $\partial\Sigma$ 的曲线积分转化为 $\Sigma$ 上的曲面积分？这正是 **Stokes 公式**要回答的问题：

$$\text{三维空间闭曲线上的环量} = \text{穿过曲面的旋度通量}$$

### 2.2 微积分基本定理的推广链

回顾整个曲线/曲面积分领域的核心思想：

| 维度 | 定理名称 | 边界积分 | 内部"导数"积累 |
|:-----|:---------|:---------|:---------------|
| 1D | 微积分基本定理 | $f(b)-f(a)$ | $\int_a^b f'(x)\,dx$ |
| 2D | Green 公式 | $\oint_{\partial D} P\,dx+Q\,dy$ | $\iint_D (\partial Q/\partial x-\partial P/\partial y)\,dA$ |
| 3D | **Stokes 公式** | $\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r}$ | $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ |
| 3D | Gauss 公式 | $\iint_{\partial\Omega} \mathbf{F}\cdot\mathbf{n}\,dS$ | $\iiint_\Omega \nabla\cdot\mathbf{F}\,dV$ |

每一行都体现了统一的模式：**边界上的积分等于内部某种"导数"的区域积分**。Stokes 公式和 Gauss 公式分别是 Green 公式在三维空间中的两个不同推广方向：

- **Stokes 公式**：将 Green 公式推广到**弯曲曲面**上，描述**环量**与**旋度通量**的关系
- **Gauss 公式**：将 Green 公式推广到**三维体积**上，描述**通量**与**散度积累**的关系

### 2.3 物理应用

- **电磁学（法拉第电磁感应定律）**：$\displaystyle \oint_{\partial\Sigma} \mathbf{E}\cdot d\mathbf{r} = -\frac{d}{dt}\iint_\Sigma \mathbf{B}\cdot\mathbf{n}\,dS$ ——电场沿闭回路的环量等于穿过该回路所围曲面的磁通量的变化率。这正是 Stokes 公式在电磁学中的直接应用。
- **流体力学**：涡线（vortex line）上的速度场的环量等于穿过涡管截面的涡量（vorticity，即速度场的旋度）的通量。
- **保守力场判断**：若 $\nabla\times\mathbf{F} = \mathbf{0}$ 在某单连通区域内处处成立，则 $\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = 0$ 对任何可收缩闭曲线成立，场是保守的。

---

## 3. 核心概念与符号约定

### 3.1 三维旋度（Curl）的完整定义（A-1）

在文件03（第十四章）中，我们接触过二维向量场 $(P,Q)$ 的**标量旋度**：$\partial Q/\partial x - \partial P/\partial y$，它度量了平面场的局部旋转强度。对于三维向量场 $\mathbf{F} = (P,Q,R)$，旋度不再是一个标量，而是一个**向量**——它既指示旋转的轴方向，也指示旋转的强度。

**定义 14.15（三维旋度 / Curl in Three Dimensions）**：设 $\mathbf{F}(x,y,z) = (P(x,y,z),\; Q(x,y,z),\; R(x,y,z))$ 是定义在区域 $\Omega \subset \mathbb{R}^3$ 上的向量场，$P,Q,R$ 均具有一阶连续偏导数。定义 $\mathbf{F}$ 的**旋度**为：
$$\boxed{\nabla\times\mathbf{F} = \left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z} - \frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right)}$$

**三种等价的表示形式**：

| 形式 | 表达式 | 适用场景 |
|:-----|:-------|:---------|
| 分量展开 | $\displaystyle \left( \frac{\partial R}{\partial y}-\frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z}-\frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y} \right)$ | 具体数值计算 |
| 行列式形式 | $\displaystyle \nabla\times\mathbf{F} = \begin{vmatrix} \mathbf{i} & \mathbf{j} & \mathbf{k} \\[2pt] \dfrac{\partial}{\partial x} & \dfrac{\partial}{\partial y} & \dfrac{\partial}{\partial z} \\[6pt] P & Q & R \end{vmatrix}$ | 记忆与推导 |
| Nabla 叉积 | $\displaystyle \nabla\times\mathbf{F} = \left(\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z}\right) \times (P,Q,R)$ | 简洁记号 |

**行列式形式的展开验证**：
$$\begin{aligned}
\nabla\times\mathbf{F} &= \begin{vmatrix}
\mathbf{i} & \mathbf{j} & \mathbf{k} \\
\partial_x & \partial_y & \partial_z \\
P & Q & R
\end{vmatrix} \\
&= \mathbf{i}\left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z} \right) - \mathbf{j}\left( \frac{\partial R}{\partial x} - \frac{\partial P}{\partial z} \right) + \mathbf{k}\left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) \\
&= \left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z} - \frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right)
\end{aligned}$$

其中 $\partial_x = \partial/\partial x$，$\partial_y = \partial/\partial y$，$\partial_z = \partial/\partial z$。

**与文件03中二维旋度的联系**：当 $\mathbf{F}$ 的 $z$ 分量 $R \equiv 0$ 且 $P,Q$ 与 $z$ 无关时（即 $\mathbf{F}$ 退化为二维场）：
$$\nabla\times\mathbf{F} = \left( 0,\; 0,\; \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right)$$

旋度的前两个分量为零，仅剩下 $z$ 分量——这正是文件03中定义的二维标量旋度。因此，二维旋度是三维旋度的 $z$ 分量，三维旋度是二维旋度的完整向量推广。

#### 旋度计算示例

**例 A**（考试题一）：$\mathbf{F} = (y, z, x)$，即 $P = y$，$Q = z$，$R = x$。

使用分量公式逐步计算：
- 第一分量：$\partial R/\partial y - \partial Q/\partial z = \partial(x)/\partial y - \partial(z)/\partial z = 0 - 1 = -1$
- 第二分量：$\partial P/\partial z - \partial R/\partial x = \partial(y)/\partial z - \partial(x)/\partial x = 0 - 1 = -1$
- 第三分量：$\partial Q/\partial x - \partial P/\partial y = \partial(z)/\partial x - \partial(y)/\partial y = 0 - 1 = -1$

因此 $\nabla\times\mathbf{F} = (-1, -1, -1)$。

**例 B**（考试题二）：$\mathbf{F} = (3z, 5x, -2y)$，即 $P = 3z$，$Q = 5x$，$R = -2y$。
- 第一分量：$\partial R/\partial y - \partial Q/\partial z = \partial(-2y)/\partial y - \partial(5x)/\partial z = -2 - 0 = -2$
- 第二分量：$\partial P/\partial z - \partial R/\partial x = \partial(3z)/\partial z - \partial(-2y)/\partial x = 3 - 0 = 3$
- 第三分量：$\partial Q/\partial x - \partial P/\partial y = \partial(5x)/\partial x - \partial(3z)/\partial y = 5 - 0 = 5$

因此 $\nabla\times\mathbf{F} = (-2, 3, 5)$。

**例 C**：$\mathbf{F} = (y^2, -xz, \sin z)$，即 $P = y^2$，$Q = -xz$，$R = \sin z$。
- 第一分量：$\partial(\sin z)/\partial y - \partial(-xz)/\partial z = 0 - (-x) = x$
- 第二分量：$\partial(y^2)/\partial z - \partial(\sin z)/\partial x = 0 - 0 = 0$
- 第三分量：$\partial(-xz)/\partial x - \partial(y^2)/\partial y = (-z) - (2y) = -z - 2y$

因此 $\nabla\times\mathbf{F} = (x,\; 0,\; -z-2y)$。

#### 旋度的物理含义

在文件06（第十四章）中，散度 $\nabla\cdot\mathbf{F}$ 描述了场在某点的"源/汇强度"（通量体密度）。旋度 $\nabla\times\mathbf{F}$ 则描述了场在某点的**旋转强度**（环量体密度），具体表现为：

- **方向**：$\nabla\times\mathbf{F}$ 的方向指示了局部旋转的**轴方向**——将一微小测试轮（风车）放入场中，旋转轴的方向就是旋度的方向
- **大小**：$\|\nabla\times\mathbf{F}\|$ 度量了旋转的**强度**——旋度越大，旋转越快
- **无旋场**：若 $\nabla\times\mathbf{F} = \mathbf{0}$ 处处成立，则场为无旋场（保守场），任何闭曲线上的环量为零

旋度与散度的对比：

| 比较项 | 散度 $\nabla\cdot\mathbf{F}$ | 旋度 $\nabla\times\mathbf{F}$ |
|:-------|:----------------------------|:------------------------------|
| 运算类型 | Nabla 点积 | Nabla 叉积 |
| 结果类型 | 标量 | 向量 |
| 物理含义 | 源/汇强度（通量体密度） | 旋转强度（环量体密度） |
| 正/负含义 | 正=源，负=汇，零=无散 | 方向由旋度向量指示 |
| 对应公式 | Gauss 公式（通量型） | Stokes 公式（环量型） |

#### 拉普拉斯算子的引入（为将来准备）

虽然本节不需要，但值得注意的是一个重要的复合运算：梯度的散度——**拉普拉斯算子**（Laplacian）：
$$\nabla^2 f = \nabla\cdot(\nabla f) = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} + \frac{\partial^2 f}{\partial z^2}$$

此外，一个重要的恒等式（将在后续场论部分证明）：**旋度的散度恒为零**。
$$\nabla\cdot(\nabla\times\mathbf{F}) = 0$$

即任何向量场的旋度都是无散场。这为电磁学中 $\nabla\cdot\mathbf{B} = 0$（磁场无散）提供了理论基础。

### 3.2 Stokes 公式的正式陈述（A-2）

有了三维旋度的定义，我们可以正式陈述 Stokes 公式。

**定理 14.16（Stokes 公式 / Stokes' Theorem）**：设 $\Sigma \subset \mathbb{R}^3$ 为分片光滑的**有向曲面**，$\partial\Sigma$ 为 $\Sigma$ 的边界曲线（分段光滑闭曲线），其方向与 $\Sigma$ 的定向满足**右手法则**（见 3.3 节）。若向量场 $\mathbf{F}(x,y,z) = (P(x,y,z),\; Q(x,y,z),\; R(x,y,z))$ 的三个分量 $P,Q,R$ 在包含 $\Sigma$ 的某个开区域上具有**连续偏导数**，则：

**有向投影形式**：
$$\boxed{\oint_{\partial\Sigma} P\,dx + Q\,dy + R\,dz = \iint_{\Sigma} \left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z} \right) dy\wedge dz + \left( \frac{\partial P}{\partial z} - \frac{\partial R}{\partial x} \right) dz\wedge dx + \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dx\wedge dy}$$

**向量点积形式**（最简洁）：
$$\boxed{\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS}$$

其中 $\mathbf{n}$ 是 $\Sigma$ 上的单位法向量（指向与 $\partial\Sigma$ 的方向满足右手法则），$dS$ 是无向面积微元。

**行列式展开形式**（便于记忆）：
$$\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} \begin{vmatrix}
\mathbf{n} & dS \\[2pt]
\nabla & \times \\
\mathbf{F} &
\end{vmatrix} \cdot \mathbf{n}\,dS = \iint_{\Sigma} \begin{vmatrix}
\mathbf{i} & \mathbf{j} & \mathbf{k} \\
\partial_x & \partial_y & \partial_z \\
P & Q & R
\end{vmatrix} \cdot \mathbf{n}\,dS$$

#### 条件分析

**条件 1：$\Sigma$ 是分片光滑的有向曲面**

- 分片光滑意味着 $\Sigma$ 可以由若干片光滑曲面拼接而成（如立方体表面有 6 个平面片，拼接处可以有棱）
- 定向（有向）意味着 $\Sigma$ 的每一点指定了法向量方向（上侧/下侧、外侧/内侧）

**条件 2：$\partial\Sigma$ 是分段光滑闭曲线，方向与 $\Sigma$ 满足右手法则**

- 边界曲线必须闭合——如果 $\Sigma$ 是封闭曲面（如球面），则 $\partial\Sigma$ 为空集，Stokes 公式右端为零（没有边界曲线积分）
- 方向必须协调——这是应用 Stokes 公式最容易出错的地方（详见 3.3 节）

**条件 3：$P,Q,R$ 在包含 $\Sigma$ 的开区域上有连续偏导数**

- 连续偏导数保证了旋度 $\nabla\times\mathbf{F}$ 有定义且可积
- 若 $\mathbf{F}$ 在 $\Sigma$ 上有奇点（如分母为零的点），需挖洞处理——与 Green 公式和 Gauss 公式的条件类似

### 3.3 右手法则：曲面定向与边界定向的对应关系（A-3）

Stokes 公式中，曲面 $\Sigma$ 的定向（法向量指向）与边界曲线 $\partial\Sigma$ 的方向必须满足某种对应关系，否则公式右端会差一个负号。这一对应关系由**右手法则**描述。

**右手法则（Right-Hand Rule）**：将右手的拇指指向曲面 $\Sigma$ 的单位法向量 $\mathbf{n}$ 方向，则其余四指弯曲的方向就是边界曲线 $\partial\Sigma$ 的**正方向**（也称为 $\partial\Sigma$ 的**正向**）。

**等价描述方法**：
- 沿 $\partial\Sigma$ 的正方向行走时，曲面 $\Sigma$ 始终位于行走方向的**左侧**
- 当 $\Sigma$ 是 $xy$ 平面上的区域取上侧时，$\partial\Sigma$ 的正向即**逆时针方向**（这与文件03中平面正向边界的定义一致——"区域在左侧"）

**示例说明**（考试题一）：设 $\Sigma$ 是平面 $x+y+z=1$ 在第一卦限的三角形，取上侧（法向量指向 $z$ 正方向一侧）。由右手法则：拇指指向上侧法向量方向，四指弯曲方向即边界曲线的正方向——从 $z$ 轴正向看去，边界为逆时针方向。题目中给出的"从 $z$ 轴正向看去取逆时针方向"正好与此一致。

**示例说明**（考试题二）：设 $\Sigma$ 是平面 $z=y+3$ 上被圆柱 $x^2+y^2=1$ 截得的椭圆盘。曲线 $L$ 从 $z$ 轴正向看去取逆时针方向。由右手法则：拇指指向曲面法向量方向，四指指向边界方向（逆时针）。因此需要检查法向量的指向是否与右手法则一致。

具体操作：当边界取逆时针时，由右手法则，拇指应指向 $z$ 轴正向（与逆时针四指相对应的拇指方向）。因此曲面应取上侧（$n_z > 0$）。平面 $z=y+3$ 的法向量为 $\nabla(z-y-3) = (0,-1,1)$，其 $n_z = 1 > 0$，符合要求。

**右手法则总结**：
- $\mathbf{n}$ 指向拇指方向 $\iff$ $\partial\Sigma$ 的正向为四指弯曲方向
- 若 $\Sigma$ 取上侧（$\mathbf{n}$ 指向 $z$ 正方向），则 $\partial\Sigma$ 为逆时针（从 $z$ 正向看）
- 若 $\Sigma$ 取下侧（$\mathbf{n}$ 指向 $z$ 负方向），则 $\partial\Sigma$ 为顺时针（从 $z$ 正向看）

### 3.4 Stokes 公式与 Green 公式、Gauss 公式的关系

| 比较项 | Green 公式（文件03） | Stokes 公式（本文件） | Gauss 公式（文件06） |
|:-------|:---------------------|:---------------------|:--------------------|
| 边界积分 | $\oint_{\partial D} P\,dx+Q\,dy$ | $\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r}$ | $\iint_{\partial\Omega} \mathbf{F}\cdot\mathbf{n}\,dS$ |
| 内部积分 | $\iint_D (\partial Q/\partial x-\partial P/\partial y)\,dA$ | $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ | $\iiint_\Omega \nabla\cdot\mathbf{F}\,dV$ |
| 内部"导数" | 二维标量旋度 | 三维向量旋度 | 散度 |
| 区域类型 | 平面区域 $D$ | 空间曲面 $\Sigma$ | 空间区域 $\Omega$ |
| 边界类型 | 平面闭曲线 | 空间闭曲线 | 封闭曲面 |
| 推广方向 | — | 平面 $\to$ 空间曲面 | 二维 $\to$ 三维体积 |

**关键理解**：
- Stokes 将 Green 从**平面闭曲线**推广到**空间任意有向曲面的边界**
- Gauss 则将 Green 从**二维区域**推广到**三维体积**
- 三者分别描述了"环量 = 旋度通量"（Stokes）和"通量 = 散度积累"（Gauss），是微积分基本定理在三维空间中不同侧面的体现

### 3.5 符号表

| 符号 | 含义 | 说明 |
|:-----|:-----|:------|
| $\nabla\times\mathbf{F}$ | 向量场 $\mathbf{F}$ 的旋度 | 三维向量，度量局部旋转的轴和强度 |
| $\nabla\cdot\mathbf{F}$ | 向量场 $\mathbf{F}$ 的散度（文件06） | 标量，度量源/汇强度 |
| $\partial\Sigma$ | 有向曲面 $\Sigma$ 的边界曲线 | 方向与 $\Sigma$ 满足右手法则 |
| $\displaystyle\oint_{\partial\Sigma}$ | 沿边界闭曲线的第二型曲线积分 | $\partial\Sigma$ 取正向 |
| $(\nabla\times\mathbf{F})\cdot\mathbf{n}$ | 旋度在法线方向上的投影 | Stokes 公式右端被积函数 |
| $\mathbf{r} = (x,y,z)$ | 位置向量 | $d\mathbf{r} = (dx,dy,dz)$ |
| $\cos\alpha,\cos\beta,\cos\gamma$ | 法向量的方向余弦 | $\mathbf{n} = (\cos\alpha,\cos\beta,\cos\gamma)$ |

---

## 4. 原理与方法

### 4.1 Stokes 公式退化到 Green 公式（A-4）

Stokes 公式是 Green 公式的自然推广。为了深刻理解这一关系，我们展示从 Stokes 公式到 Green 公式的完整退化过程。

**退化设定**：设 $\Sigma$ 是 $xy$ 平面上的一个有界区域 $D$，取**上侧**（法向量 $\mathbf{n} = (0,0,1)$，$dS = dx\,dy$）。$\partial\Sigma = C$ 是 $D$ 的正向边界（逆时针方向）。向量场 $\mathbf{F} = (P(x,y,z), Q(x,y,z), R(x,y,z))$。

**退化推导**：

**第 1 步**：写出 Stokes 公式的向量形式：
$$\oint_C \mathbf{F}\cdot d\mathbf{r} = \iint_D (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$$

**第 2 步**：代入 $\mathbf{n} = (0,0,1)$：
$$(\nabla\times\mathbf{F})\cdot\mathbf{n} = \left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z} - \frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) \cdot (0,0,1) = \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}$$

**第 3 步**：代入 $\mathbf{F}\cdot d\mathbf{r}$。在 $z=0$ 平面上，$dz = 0$，因此：
$$\mathbf{F}\cdot d\mathbf{r} = P\,dx + Q\,dy + R\,dz = P\,dx + Q\,dy$$

**第 4 步**：代入 $dS = dx\,dy$（因为 $\Sigma$ 是平面区域取上侧，$dS = dx\,dy$）。

**第 5 步**：合并结果：
$$\boxed{\oint_C P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dx\,dy}$$

这正是 **Green 公式**（定理 14.7）。因此，**Green 公式是 Stokes 公式在曲面退化为 $xy$ 平面区域时的特例**。

**物理理解**：当曲面 $\Sigma$ 是平面区域时，旋度在法线方向的分量恰好退化为二维标量旋度 $\partial Q/\partial x - \partial P/\partial y$，边界曲线积分退化为平面闭曲线积分。一切"三维"信息消失，Stokes 公式回归到二维的 Green 公式。

**重要启示**：这一定位关系告诉我们，任何适用于 Green 公式的题目都可以看作 Stokes 公式的特例。反过来，Stokes 公式是 Green 公式在空间弯曲曲面上的推广——当你将弯曲曲面"压平"到坐标平面上时，Stokes 公式就变成了 Green 公式。

### 4.2 Stokes 公式应用标准工作流（A-5）

应用 Stokes 公式将空间闭曲线上的第二型曲线积分转化为曲面积分，遵循以下五步工作流。

**第 1 步：计算旋度 $\nabla\times\mathbf{F}$**

由定义 14.15，计算向量场 $\mathbf{F} = (P,Q,R)$ 的旋度：
$$\nabla\times\mathbf{F} = \left( \frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z},\; \frac{\partial P}{\partial z} - \frac{\partial R}{\partial x},\; \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right)$$

使用行列式形式或分量公式均可。

**第 2 步：确定曲面 $\Sigma$ 及其法向量 $\mathbf{n}$（右手法则验证）**

- 选择以给定闭曲线 $C$ 为边界的曲面 $\Sigma$（自由选择——只要 $C = \partial\Sigma$ 且 $\Sigma$ 光滑即可）
- 确定 $\Sigma$ 的定向：使 $\Sigma$ 的法向量指向与 $C$ 的正向满足右手法则
- 计算 $\Sigma$ 的单位法向量 $\mathbf{n}$

**第 3 步：计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$**

将旋度向量投影到法线方向，得到标量被积函数。

**第 4 步：写出曲面积分**

将原曲线积分转化为第一类曲面积分：
$$\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$$

**第 5 步：计算曲面面积微元 $dS$ 并求值**

- 若 $\Sigma$ 是显式曲面 $z = f(x,y)$，则 $dS = \sqrt{1 + f_x^2 + f_y^2}\;dx\,dy$
- 若 $\Sigma$ 是平面，常用其面积（当 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$ 为常数时直接面积 $\times$ 常数）
- 转换为投影区域上的二重积分并计算

**工作流速查表**：

| 步骤 | 操作 | 示例（考试题一） |
|:-----|:-----|:----------------|
| Step 1 | 计算 $\nabla\times\mathbf{F}$ | $\nabla\times(y,z,x) = (-1,-1,-1)$ |
| Step 2 | 确定 $\Sigma$ 和 $\mathbf{n}$（右手法则） | $x+y+z=1$，上侧 $\mathbf{n}=(1,1,1)/\sqrt{3}$ |
| Step 3 | 计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$ | $(-1,-1,-1)\cdot(1,1,1)/\sqrt{3} = -\sqrt{3}$ |
| Step 4 | 写出 $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ | $\iint_\Sigma (-\sqrt{3})\,dS$ |
| Step 5 | 计算 $dS$ 并求值 | $dS=\sqrt{3}\,dx\,dy$，$I=-3\times\text{Area}(D)=-3/2$ |

---

## 5. 例题

### 例 1：平面三角形上的 Stokes 公式

设向量场 $\mathbf{F}(x,y,z) = (y,\; z,\; x)$。$L$ 是平面 $x+y+z=1$ 在第一卦限中被坐标平面截得的三角形边界，从 $z$ 轴正向看去取逆时针方向。计算曲线积分：

$$I = \oint_L y\,dx + z\,dy + x\,dz$$

**解**：

#### 方法一：Stokes 公式法

**Step 1：计算旋度**

由定义 14.15，对 $\mathbf{F} = (y,z,x)$（$P=y,\; Q=z,\; R=x$）：
$$\nabla\times\mathbf{F} = \left( \frac{\partial x}{\partial y} - \frac{\partial z}{\partial z},\; \frac{\partial y}{\partial z} - \frac{\partial x}{\partial x},\; \frac{\partial z}{\partial x} - \frac{\partial y}{\partial y} \right) = (0-1,\; 0-1,\; 0-1) = (-1,\; -1,\; -1)$$

**Step 2：确定曲面和法向量**

取 $\Sigma$ 为 $L$ 所围的平面三角形区域 $x+y+z=1$（$x,y,z \geq 0$），取上侧。

由 3.1 节"法向量"：$F(x,y,z)=x+y+z-1=0$，$\nabla F = (1,1,1)$。

单位法向量（上侧，$n_z = 1/\sqrt{3} > 0$）：
$$\mathbf{n} = \frac{(1,1,1)}{\sqrt{3}}$$

右手法则验证：拇指指向上侧法向量 $(1,1,1)$ 方向，四指弯曲方向为逆时针（从 $z$ 正向看），与题目给定的 $L$ 方向一致。$\checkmark$

**Step 3：计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$**

$$(\nabla\times\mathbf{F})\cdot\mathbf{n} = (-1,-1,-1)\cdot\frac{(1,1,1)}{\sqrt{3}} = -\frac{1+1+1}{\sqrt{3}} = -\frac{3}{\sqrt{3}} = -\sqrt{3}$$

**Step 4：写出曲面积分**

由 Stokes 公式：
$$I = \iint_{\Sigma} (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS = \iint_{\Sigma} (-\sqrt{3})\,dS$$

**Step 5：计算 $dS$ 并求值**

$\Sigma$ 是平面 $z = f(x,y) = 1 - x - y$，$(x,y) \in D$，其中 $D$ 是 $xy$ 平面上的三角形：
$$D = \{(x,y) \mid x \geq 0,\; y \geq 0,\; x+y \leq 1\}$$

$f_x = -1$，$f_y = -1$，因此：
$$dS = \sqrt{1 + f_x^2 + f_y^2}\;dx\,dy = \sqrt{1 + 1 + 1}\;dx\,dy = \sqrt{3}\;dx\,dy$$

$$\begin{aligned}
I &= \iint_{\Sigma} (-\sqrt{3})\,dS = \iint_D (-\sqrt{3})\cdot\sqrt{3}\;dx\,dy \\
  &= -3 \iint_D dx\,dy = -3 \times \text{Area}(D)
\end{aligned}$$

三角形 $D$ 的面积 $= \frac12 \times 1 \times 1 = \frac12$，因此：
$$I = -3 \times \frac12 = -\frac{3}{2}$$

**结果**：$\boxed{I = -\dfrac{3}{2}}$。

#### 方法二：直接参数化（验证）

三角形 $L$ 由三条边 $L_1$（A$\to$B）、$L_2$（B$\to$C）、$L_3$（C$\to$A）组成，其中 $A(1,0,0)$、$B(0,1,0)$、$C(0,0,1)$。

**$L_1$（A$\to$B）**：在 $z=0$ 平面上从 $(1,0,0)$ 到 $(0,1,0)$。
参数化：$x = 1-t$，$y = t$，$z = 0$，$t$ 从 $0$ 到 $1$。
$dx = -dt$，$dy = dt$，$dz = 0$。
$$\begin{aligned}
\int_{L_1} y\,dx + z\,dy + x\,dz &= \int_0^1 [t\cdot(-dt) + 0\cdot dt + (1-t)\cdot 0] \\
&= -\int_0^1 t\,dt = -\frac12
\end{aligned}$$

**$L_2$（B$\to$C）**：在 $x=0$ 平面上从 $(0,1,0)$ 到 $(0,0,1)$。
参数化：$x = 0$，$y = 1-t$，$z = t$，$t$ 从 $0$ 到 $1$。
$dx = 0$，$dy = -dt$，$dz = dt$。
$$\begin{aligned}
\int_{L_2} y\,dx + z\,dy + x\,dz &= \int_0^1 [ (1-t)\cdot 0 + t\cdot(-dt) + 0\cdot dt] \\
&= -\int_0^1 t\,dt = -\frac12
\end{aligned}$$

**$L_3$（C$\to$A）**：在 $y=0$ 平面上从 $(0,0,1)$ 到 $(1,0,0)$。
参数化：$x = t$，$y = 0$，$z = 1-t$，$t$ 从 $0$ 到 $1$。
$dx = dt$，$dy = 0$，$dz = -dt$。
$$\begin{aligned}
\int_{L_3} y\,dx + z\,dy + x\,dz &= \int_0^1 [ 0\cdot dt + (1-t)\cdot 0 + t\cdot(-dt)] \\
&= -\int_0^1 t\,dt = -\frac12
\end{aligned}$$

**求和**：$I = -\frac12 - \frac12 - \frac12 = -\frac32$。与 Stokes 公式结果一致。$\checkmark$

---

### 例 2：倾斜椭圆盘上的 Stokes 公式

设向量场 $\mathbf{F}(x,y,z) = (3z,\; 5x,\; -2y)$。$L$ 是圆柱面 $x^2+y^2=1$ 与平面 $z = y+3$ 的交线（椭圆），从 $z$ 轴正向看去取逆时针方向。计算曲线积分：

$$I = \oint_L 3z\,dx + 5x\,dy - 2y\,dz$$

**解**：

**Step 1：计算旋度**

$\mathbf{F} = (3z, 5x, -2y)$，即 $P = 3z$，$Q = 5x$，$R = -2y$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(-2y)}{\partial y} - \frac{\partial(5x)}{\partial z},\; \frac{\partial(3z)}{\partial z} - \frac{\partial(-2y)}{\partial x},\; \frac{\partial(5x)}{\partial x} - \frac{\partial(3z)}{\partial y} \right) \\
&= (-2 - 0,\; 3 - 0,\; 5 - 0) = (-2,\; 3,\; 5)
\end{aligned}$$

**Step 2：确定曲面和法向量**

取 $\Sigma$ 为 $L$ 所围的平面区域——平面 $z = y+3$ 上 $x^2+y^2 \leq 1$ 截得的椭圆盘。

由右手法则：$L$ 从 $z$ 正向看为逆时针，因此拇指应指向 $z$ 正向，曲面取上侧。检查平面 $z = y+3$ 的法向量：将平面写为 $F(x,y,z) = z - y - 3 = 0$，则 $\nabla F = (0, -1, 1)$。其 $n_z = 1 > 0$，指向上侧。符合右手法则。$\checkmark$

单位法向量：
$$\mathbf{n} = \frac{(0,-1,1)}{\sqrt{0^2+(-1)^2+1^2}} = \left(0,\; -\frac{1}{\sqrt{2}},\; \frac{1}{\sqrt{2}}\right)$$

**Step 3：计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$**

$$\begin{aligned}
(\nabla\times\mathbf{F})\cdot\mathbf{n} &= (-2,\; 3,\; 5)\cdot\left(0,\; -\frac{1}{\sqrt{2}},\; \frac{1}{\sqrt{2}}\right) \\
&= (-2)\cdot 0 + 3\cdot\left(-\frac{1}{\sqrt{2}}\right) + 5\cdot\frac{1}{\sqrt{2}} \\
&= -\frac{3}{\sqrt{2}} + \frac{5}{\sqrt{2}} = \frac{2}{\sqrt{2}} = \sqrt{2}
\end{aligned}$$

**Step 4：写出曲面积分**

由 Stokes 公式：
$$I = \iint_{\Sigma} (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS = \iint_{\Sigma} \sqrt{2}\;dS = \sqrt{2} \times \text{Area}(\Sigma)$$

**Step 5：计算面积并求值**

需要计算椭圆盘 $\Sigma$ 的面积。$\Sigma$ 是平面 $z=y+3$ 上的区域，其投影到 $xy$ 平面为单位圆盘 $x^2+y^2 \leq 1$。

$z = y+3$，$z_x = 0$，$z_y = 1$，因此：
$$dS = \sqrt{1 + 0^2 + 1^2}\;dx\,dy = \sqrt{2}\;dx\,dy$$

椭圆盘面积：
$$\text{Area}(\Sigma) = \iint_\Sigma 1\,dS = \iint_{D} \sqrt{2}\;dx\,dy = \sqrt{2} \times \text{Area}(D) = \sqrt{2} \times \pi \times 1^2 = \sqrt{2}\,\pi$$

因此：
$$I = \sqrt{2} \times \sqrt{2}\,\pi = 2\pi$$

**结果**：$\boxed{I = 2\pi}$。

**关键看点**：$(\nabla\times\mathbf{F})\cdot\mathbf{n} = \sqrt{2}$ 和 $dS = \sqrt{2}\,dx\,dy$ 中的 $\sqrt{2}$ 因子相乘得到 $2$，使得最终积分简化为 $2 \times \text{Area}(D) = 2\pi$。这展示了在工作流中正确应用 $dS$ 公式的重要性。

---

### 例 3：Stokes 公式退化到 Green 公式

设向量场 $\mathbf{F} = (y,\; x^2,\; z)$，$L$ 是圆周 $x^2+y^2=R^2$（$R>0$），$z=0$，从 $z$ 轴正向看去取逆时针方向。计算曲线积分：

$$I = \oint_L y\,dx + x^2\,dy + z\,dz$$

**解法一：Stokes 公式法（3D 视角）**

**Step 1：计算旋度**

$\mathbf{F} = (y, x^2, z)$，即 $P = y$，$Q = x^2$，$R = z$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial z}{\partial y} - \frac{\partial(x^2)}{\partial z},\; \frac{\partial y}{\partial z} - \frac{\partial z}{\partial x},\; \frac{\partial(x^2)}{\partial x} - \frac{\partial y}{\partial y} \right) \\
&= (0 - 0,\; 0 - 0,\; 2x - 1) = (0,\; 0,\; 2x - 1)
\end{aligned}$$

**Step 2：确定曲面和法向量**

取 $\Sigma$ 为 $L$ 所围的圆盘 $x^2+y^2 \leq R^2$，位于 $z=0$ 平面上，取上侧。法向量 $\mathbf{n} = (0,0,1)$。

右手法则：上侧时边界为逆时针，与题目给定方向一致。$\checkmark$

**Step 3：计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$**

$$(\nabla\times\mathbf{F})\cdot\mathbf{n} = (0,0,2x-1)\cdot(0,0,1) = 2x - 1$$

**Step 4：写出曲面积分**

由 Stokes 公式：
$$I = \iint_{\Sigma} (2x - 1)\,dS$$

**Step 5：计算曲面积分**

在 $z=0$ 平面上，$dS = dx\,dy$（因为 $z_x = z_y = 0$，$\sqrt{1+0+0}=1$）。$D$ 是圆盘 $x^2+y^2 \leq R^2$。

$$\begin{aligned}
I &= \iint_D (2x - 1)\,dx\,dy \\
  &= 2\iint_D x\,dx\,dy - \iint_D 1\,dx\,dy
\end{aligned}$$

第一项：由对称性，$x$ 在圆盘上的积分为零（$x$ 是奇函数，区域关于 $x$ 轴对称）。因此 $\iint_D x\,dx\,dy = 0$。

第二项：$\iint_D 1\,dx\,dy = \text{Area}(D) = \pi R^2$。

因此：
$$I = 2 \times 0 - \pi R^2 = -\pi R^2$$

**结果**：$\boxed{I = -\pi R^2}$。

**解法二：Green 公式法（退化方法，验证）**

在 $z=0$ 平面上，$dz = 0$，因此 $z\,dz = 0$。曲线积分退化为：
$$I = \oint_L y\,dx + x^2\,dy$$

这正是 $xy$ 平面上的第二型闭曲线积分。由 Green 公式（定理 14.7），$P = y$，$Q = x^2$：
$$\frac{\partial Q}{\partial x} = 2x,\quad \frac{\partial P}{\partial y} = 1$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 2x - 1$$

因此：
$$I = \iint_D (2x - 1)\,dA = 2\iint_D x\,dA - \pi R^2 = 0 - \pi R^2 = -\pi R^2$$

两种方法结果一致。$\checkmark$

**退化确认**：本题中，曲面 $\Sigma$ 正是 $xy$ 平面上的区域 $D$，$\mathbf{n} = (0,0,1)$，Stokes 公式直接退化为了 Green 公式。特别注意到 $\nabla\times\mathbf{F} = (0,0,2x-1)$ 的前两个分量为零——退化的标志性特征。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|:---------|:---------|
| 旋度 $\nabla\times\mathbf{F}$ 和散度 $\nabla\cdot\mathbf{F}$ 只是记号不同，本质是一回事 | 旋度是向量（Nabla 叉积结果），散度是标量（Nabla 点积结果）。旋度度量旋转，散度度量源汇。只有三维以上空间才有旋度向量，散度在各维度都存在 |
| Stokes 公式的右端 $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ 中的 $\mathbf{n}\,dS$ 直接写成 $dx\,dy$ 等而不考虑定向 | $\mathbf{n}\,dS$ 转化为有向投影时需注意：$n_z dS = dx \wedge dy$，符号由法向量指向决定（上侧正、下侧负） |
| 任意给定空间闭曲线，都能直接用 Stokes 公式 | 需要满足：(1) 存在以该曲线为边界的可定向曲面；(2) $P,Q,R$ 在曲面上有连续偏导数。若曲线是空间结（如 trefoil knot）或场有奇点，需谨慎处理 |
| 使用 Stokes 公式时，曲面定向和边界方向可以任意选择 | 两者必须通过右手法则协调。如果方向不一致，需要在 Stokes 公式右端加负号：$\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = -\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ |
| Stokes 公式和 Green 公式没有关系 | Green 公式是 Stokes 公式在曲面为 $xy$ 平面区域时的特例。当 $\Sigma$ 退化为平面区域 $D$ 且 $\mathbf{n} = (0,0,1)$ 时，$(\nabla\times\mathbf{F})\cdot\mathbf{n} = \partial Q/\partial x - \partial P/\partial y$，Stokes 公式直接退化到 Green 公式 |
| 只有向量点积形式 $\oint \mathbf{F}\cdot d\mathbf{r} = \iint (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ 一种写法 | 还有有向投影形式，将 $\mathbf{n}\,dS$ 展开为 $dy\wedge dz,\; dz\wedge dx,\; dx\wedge dy$ 三项，便于具体计算 |
| 用 Stokes 公式计算时，旋度一定是常数才能用 | 旋度可以是位置函数，通过曲面积分计算即可。如例 3 中 $(\nabla\times\mathbf{F})\cdot\mathbf{n} = 2x-1$ 不是常数，但通过对称性简化计算 |
| Stokes 公式只在曲面是平面时成立 | Stokes 公式适用于任何分片光滑的有向曲面，包括弯曲曲面（如球冠、抛物面等）。只是在弯曲曲面上计算 $dS$ 更复杂 |

### 检查点

- [ ] 能否写出三维旋度的完整定义式（三个分量）和行列式记忆形式？
- [ ] 能否用自己的语言完整叙述 Stokes 公式（包括有向投影形式和向量点积形式）以及所有条件？
- [ ] 能否解释右手法则——给定曲面法向量方向，如何确定边界曲线的正方向？
- [ ] 能否完成 Stokes 公式五步工作流的每一步：计算旋度 $\to$ 确定曲面和法向量 $\to$ 计算 $(\nabla\times\mathbf{F})\cdot\mathbf{n}$ $\to$ 写出曲面积分 $\to$ 计算 $dS$ 并求值？
- [ ] 能否推导 Stokes 公式退化到 Green 公式的完整过程？
- [ ] 能否区分 Stokes 公式（环量型）、Green 公式（平面环量型）和 Gauss 公式（通量型）三者的联系与区别？
- [ ] 给定一个具体的空间闭曲线和向量场，能否判断是直接用参数化方法还是用 Stokes 公式更简便？
- [ ] 如果用 Stoke公式且曲面定向与边界方向不一致（违反了右手法则），公式如何修正？
- [ ] 能否说明二维标量旋度 $\partial Q/\partial x - \partial P/\partial y$ 与三维旋度 $\nabla\times\mathbf{F}$ 的关系？

---

## 练习题

### 基础巩固

**1.**（旋度计算）计算下列向量场的旋度 $\nabla\times\mathbf{F}$：

**(1)** $\mathbf{F} = (x^2,\; y^2,\; z^2)$

**(2)** $\mathbf{F} = (yz,\; zx,\; xy)$

**(3)** $\mathbf{F} = (e^x,\; \sin y,\; \cos z)$

**(4)** $\mathbf{F} = (2x + y,\; x - 3y,\; z^2)$

<details><summary>参考答案</summary>

**(1)** $P = x^2$，$Q = y^2$，$R = z^2$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(z^2)}{\partial y} - \frac{\partial(y^2)}{\partial z},\; \frac{\partial(x^2)}{\partial z} - \frac{\partial(z^2)}{\partial x},\; \frac{\partial(y^2)}{\partial x} - \frac{\partial(x^2)}{\partial y} \right) \\
&= (0 - 0,\; 0 - 0,\; 0 - 0) = (0,0,0)
\end{aligned}$$

$\mathbf{F} = (x^2, y^2, z^2)$ 的旋度为零向量——这是一个无旋场（保守场）。

**(2)** $P = yz$，$Q = zx$，$R = xy$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(xy)}{\partial y} - \frac{\partial(zx)}{\partial z},\; \frac{\partial(yz)}{\partial z} - \frac{\partial(xy)}{\partial x},\; \frac{\partial(zx)}{\partial x} - \frac{\partial(yz)}{\partial y} \right) \\
&= (x - x,\; y - y,\; z - z) = (0,0,0)
\end{aligned}$$

也是一个无旋场。

**(3)** $P = e^x$，$Q = \sin y$，$R = \cos z$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(\cos z)}{\partial y} - \frac{\partial(\sin y)}{\partial z},\; \frac{\partial(e^x)}{\partial z} - \frac{\partial(\cos z)}{\partial x},\; \frac{\partial(\sin y)}{\partial x} - \frac{\partial(e^x)}{\partial y} \right) \\
&= (0 - 0,\; 0 - 0,\; 0 - 0) = (0,0,0)
\end{aligned}$$

三个分量各自只依赖一个变量，所有交叉偏导数均为零，旋度为零向量。

**(4)** $P = 2x + y$，$Q = x - 3y$，$R = z^2$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(z^2)}{\partial y} - \frac{\partial(x-3y)}{\partial z},\; \frac{\partial(2x+y)}{\partial z} - \frac{\partial(z^2)}{\partial x},\; \frac{\partial(x-3y)}{\partial x} - \frac{\partial(2x+y)}{\partial y} \right) \\
&= (0 - 0,\; 0 - 0,\; 1 - 1) = (0,0,0)
\end{aligned}$$

也是无旋场。

**总结**：以上四个场均为无旋场（旋度为零）。这听起来也许令人惊讶，但正是因为它们都是"良态"的保守场。在实际问题中，旋度非零的场往往具有某种"旋转"结构，如 $\mathbf{F} = (-y, x, 0)$ 的旋度为 $(0,0,2) \neq \mathbf{0}$。

</details>

---

**2.**（条件判断）判断下列曲线积分能否直接用 Stokes 公式转化为曲面积分。若能，写出 $\nabla\times\mathbf{F}$；若不能，说明原因。

**(1)** $\displaystyle \oint_L x\,dx + y\,dy + z\,dz$，$L$ 是空间任意分段光滑闭曲线

**(2)** $\displaystyle \oint_L \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy + z\,dz$，$L$ 是 $xy$ 平面上的圆周 $x^2+y^2=1$（从 $z$ 正向看逆时针）

**(3)** $\displaystyle \int_L y\,dx + z\,dy + x\,dz$，$L$ 是从 $A(0,0,0)$ 到 $B(1,1,1)$ 的直线段

<details><summary>参考答案</summary>

**(1)** 可以。$\mathbf{F} = (x,y,z)$，$P=x$，$Q=y$，$R=z$。

$$\nabla\times\mathbf{F} = \left( \frac{\partial z}{\partial y} - \frac{\partial y}{\partial z},\; \frac{\partial x}{\partial z} - \frac{\partial z}{\partial x},\; \frac{\partial y}{\partial x} - \frac{\partial x}{\partial y} \right) = (0,0,0)$$

因此由 Stokes 公式，$\oint_L x\,dx+y\,dy+z\,dz = \iint_\Sigma \mathbf{0}\,dS = 0$ 对任何闭曲线 $L$ 成立。$\checkmark$

注意：旋度为零说明 $\mathbf{F}$ 是保守场，其势函数为 $f(x,y,z) = \frac12(x^2+y^2+z^2)$。沿闭曲线积分恒为零。

**(2)** **可以，但需注意奇点问题**。$\mathbf{F} = \left( \dfrac{-y}{x^2+y^2},\; \dfrac{x}{x^2+y^2},\; z \right)$。

计算旋度：在定义域内（排除 $z$ 轴 $x=y=0$ 上的点），
$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}\left(\frac{x}{x^2+y^2}\right) = \frac{(x^2+y^2) - x(2x)}{(x^2+y^2)^2} = \frac{y^2 - x^2}{(x^2+y^2)^2}$$

$$\frac{\partial P}{\partial y} = \frac{\partial}{\partial y}\left(\frac{-y}{x^2+y^2}\right) = -\frac{(x^2+y^2) - y(2y)}{(x^2+y^2)^2} = -\frac{x^2 - y^2}{(x^2+y^2)^2} = \frac{y^2 - x^2}{(x^2+y^2)^2}$$

所以 $\partial Q/\partial x - \partial P/\partial y = 0$。还需计算其他分量——但关键是：$P$ 和 $Q$ 在 $L$ 所围圆的内部有奇点 $(0,0)$，而 $L$ 是 $xy$ 平面上的圆周 $x^2+y^2=1$。如果选取以 $L$ 为边界的曲面为圆盘，则奇点位于曲面上，不满足连续性条件。但 $z\,dz$ 项没有问题。**需要选择避开奇点的曲面**（如螺旋曲面等）。

更简单的做法是直接利用退化：在 $z=0$ 平面上，$dz=0$，积分退化为二维问题。对 $P = -y/(x^2+y^2)$，$Q = x/(x^2+y^2)$，由 Green 公式的复连通版本（文件03练习题6），在 $x^2+y^2=1$ 上的积分值为 $2\pi$。

**(3)** **不能**。Stokes 公式要求积分曲线是**封闭**的（$\partial\Sigma$ 必须是闭曲线）。这里 $L$ 是开放线段，没有围成曲面，无法使用 Stokes 公式。应使用第二型曲线积分的参数化方法直接计算。

</details>

---

**3.**（Stokes 公式正向应用）利用 Stokes 公式计算下列曲线积分：

**(1)** $\displaystyle \oint_L y\,dx + z\,dy + x\,dz$，$L$ 是四个点 $(1,0,0)$、$(0,1,0)$、$(0,0,1)$ 连成的三角形边界（方向合适）

**(2)** $\displaystyle \oint_L (y-z)\,dx + (z-x)\,dy + (x-y)\,dz$，$L$ 是椭圆 $x^2+y^2=1$，$z=2$（从 $z$ 正向看逆时针）

<details><summary>参考答案</summary>

**(1)** 本题与例1相同，$\mathbf{F} = (y,z,x)$。由例1，$\nabla\times\mathbf{F} = (-1,-1,-1)$。取 $\Sigma$ 为三角形平面 $x+y+z=1$（上侧），$\mathbf{n} = (1,1,1)/\sqrt{3}$。

$$(\nabla\times\mathbf{F})\cdot\mathbf{n} = (-1,-1,-1)\cdot(1,1,1)/\sqrt{3} = -\sqrt{3}$$

$dS = \sqrt{3}\,dx\,dy$，投影区域 $D$ 为 $xy$ 平面三角形 $x \geq 0, y \geq 0, x+y \leq 1$，面积 $= 1/2$。

$$\oint_L = \iint_{\Sigma} (-\sqrt{3})\,dS = -3 \times \frac12 = -\frac32$$

**(2)** $\mathbf{F} = (y-z,\; z-x,\; x-y)$，即 $P = y-z$，$Q = z-x$，$R = x-y$。

计算旋度：
$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(x-y)}{\partial y} - \frac{\partial(z-x)}{\partial z},\; \frac{\partial(y-z)}{\partial z} - \frac{\partial(x-y)}{\partial x},\; \frac{\partial(z-x)}{\partial x} - \frac{\partial(y-z)}{\partial y} \right) \\
&= (1 - 1,\; 1 - 1,\; 1 - 1) = (0,0,0)
\end{aligned}$$

旋度为零，因此由 Stokes 公式：
$$\oint_L (y-z)\,dx + (z-x)\,dy + (x-y)\,dz = \iint_{\Sigma} (0,0,0)\cdot\mathbf{n}\,dS = 0$$

无论 $L$ 是什么闭曲线，该积分恒为零。这是保守场的特征——事实上，$\mathbf{F}$ 的势函数可验证为 $f(x,y,z) = xy + yz + zx$ 或类似形式。

</details>

---

### 迁移应用

**4.**（右手法则验证）设 $\mathbf{F} = (y, -x, z^2)$，$\Sigma$ 是上半球面 $z = \sqrt{1 - x^2 - y^2}$（$x^2+y^2 \leq 1$），取下侧。$L$ 是 $\Sigma$ 的边界（即 $xy$ 平面上的圆周 $x^2+y^2=1$，$z=0$）。

**(1)** 对于 $\Sigma$ 取下侧，由右手法则确定 $L$ 的正向。

**(2)** 计算 $\nabla\times\mathbf{F}$。

**(3)** 利用 Stokes 公式计算 $\displaystyle \oint_L \mathbf{F}\cdot d\mathbf{r}$。

<details><summary>参考答案</summary>

**(1)** 右手法则：拇指指向法向量方向，四指弯曲方向为边界正方向。

$\Sigma$ 取下侧：法向量指向 $z$ 轴负方向。因此拇指指向 $z$ 负方向，四指弯曲方向为**顺时针**（从 $z$ 轴正向看去）。

所以 $L$ 的正向为顺时针（从 $z$ 正向看）。即：沿 $L$ 行走时，曲面 $\Sigma$ 在左侧。

**(2)** $\mathbf{F} = (y, -x, z^2)$，$P = y$，$Q = -x$，$R = z^2$。

$$\begin{aligned}
\nabla\times\mathbf{F} &= \left( \frac{\partial(z^2)}{\partial y} - \frac{\partial(-x)}{\partial z},\; \frac{\partial y}{\partial z} - \frac{\partial(z^2)}{\partial x},\; \frac{\partial(-x)}{\partial x} - \frac{\partial y}{\partial y} \right) \\
&= (0 - 0,\; 0 - 0,\; -1 - 1) = (0,0,-2)
\end{aligned}$$

**(3)** 由 Stokes 公式：
$$\oint_L \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$$

$\nabla\times\mathbf{F} = (0,0,-2)$。$\Sigma$ 是上半球面取下侧。

上半球面的法向量（上侧时为 $(x,y,z)/\sqrt{x^2+y^2+z^2} = (x,y,z)$，因为 $x^2+y^2+z^2=1$）。取下侧时，$\mathbf{n} = -(x,y,z) = (-x,-y,-z)$。

因此：
$$(\nabla\times\mathbf{F})\cdot\mathbf{n} = (0,0,-2)\cdot(-x,-y,-z) = (-2)(-z) = 2z$$

代入 Stokes 公式：
$$\oint_L \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} 2z\,dS$$

计算曲面积分（上半球面取下侧）。用球坐标或投影法：上半球面 $z = \sqrt{1-x^2-y^2}$，$dS = \frac{1}{\sqrt{1-x^2-y^2}}\,dx\,dy$。

$$\begin{aligned}
\iint_{\Sigma} 2z\,dS &= \iint_D 2\sqrt{1-x^2-y^2} \cdot \frac{1}{\sqrt{1-x^2-y^2}}\,dx\,dy \\
&= \iint_D 2\,dx\,dy = 2 \times \pi \times 1^2 = 2\pi
\end{aligned}$$

**结果**：$\boxed{\displaystyle \oint_L \mathbf{F}\cdot d\mathbf{r} = 2\pi}$。

**验证**：直接参数化圆周 $L$（顺时针方向）：
$x = \cos t$，$y = -\sin t$，$t$ 从 $0$ 到 $2\pi$（顺时针参数化）。
$dx = -\sin t\,dt$，$dy = -\cos t\,dt$，$dz = 0$。
$$\begin{aligned}
\oint_L y\,dx - x\,dy + z^2\,dz &= \int_0^{2\pi} [(-\sin t)(-\sin t) - (\cos t)(-\cos t) + 0]\,dt \\
&= \int_0^{2\pi} (\sin^2 t + \cos^2 t)\,dt = \int_0^{2\pi} 1\,dt = 2\pi \quad \checkmark
\end{aligned}$$

</details>

---

**5.**（思考题——利用旋度的恒等式）证明：若 $\mathbf{F} = \nabla f$ 是某个标量函数 $f$ 的梯度场，则 $\nabla\times\mathbf{F} = \mathbf{0}$。并由此说明保守力场的闭曲线积分恒为零与 Stokes 公式的一致性。

<details><summary>参考答案</summary>

**证明**：设 $\mathbf{F} = \nabla f = \left( \dfrac{\partial f}{\partial x},\; \dfrac{\partial f}{\partial y},\; \dfrac{\partial f}{\partial z} \right)$，即 $P = \partial f/\partial x$，$Q = \partial f/\partial y$，$R = \partial f/\partial z$。

计算 $\nabla\times\mathbf{F}$ 的第一分量：
$$\frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z} = \frac{\partial}{\partial y}\left(\frac{\partial f}{\partial z}\right) - \frac{\partial}{\partial z}\left(\frac{\partial f}{\partial y}\right) = \frac{\partial^2 f}{\partial y\partial z} - \frac{\partial^2 f}{\partial z\partial y}$$

由混合偏导数的对称性（Clairaut 定理），当 $f$ 具有连续的二阶偏导数时，$\partial^2 f/\partial y\partial z = \partial^2 f/\partial z\partial y$，因此：
$$\frac{\partial R}{\partial y} - \frac{\partial Q}{\partial z} = 0$$

同理可证第二分量和第三分量也为零。所以 $\nabla\times(\nabla f) = \mathbf{0}$。

**物理解释**：梯度场 $\nabla f$ 的旋度恒为零——梯度场总是无旋场（保守场）。这与我们之前的理解一致：在保守力场中，沿任何闭曲线做功为零。

**与 Stokes 公式的一致性**：若 $\nabla\times\mathbf{F} = \mathbf{0}$，则由 Stokes 公式，对任何以 $\Sigma$ 为边界的闭曲线 $\partial\Sigma$：
$$\oint_{\partial\Sigma} \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma} \mathbf{0}\cdot\mathbf{n}\,dS = 0$$

因此任何保守力场在闭曲线上的环量恒为零。$\checkmark$

**重要恒等式**：$\nabla\times(\nabla f) \equiv \mathbf{0}$ 是向量分析中的一个基本恒等式，称为"梯度的旋度恒为零"。它与另一个重要恒等式 $\nabla\cdot(\nabla\times\mathbf{F}) \equiv 0$（"旋度的散度恒为零"）一起构成了向量分析的基础。

</details>

---

**6.**（思考题——Stokes 公式与曲面的选择无关性）设 $\mathbf{F}$ 和曲线 $L$ 固定。$\Sigma_1$ 和 $\Sigma_2$ 是以 $L$ 为边界的两个不同的有向曲面（方向与 $L$ 满足右手法则）。由 Stokes 公式：
$$\iint_{\Sigma_1} (\nabla\times\mathbf{F})\cdot\mathbf{n}_1\,dS = \oint_L \mathbf{F}\cdot d\mathbf{r} = \iint_{\Sigma_2} (\nabla\times\mathbf{F})\cdot\mathbf{n}_2\,dS$$

这说明 $\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ 的值只依赖于边界曲线 $L$，而**不依赖于 $\Sigma$ 的具体形状**（只要 $L = \partial\Sigma$ 且方向一致）。利用这一性质，解释为什么在例 2 中我们可以取 $\Sigma$ 为椭圆盘（而不必取圆柱面本身）。

<details><summary>参考答案</summary>

**解答**：在例 2 中，$L$ 是圆柱面 $x^2+y^2=1$ 与平面 $z=y+3$ 的交线。$L$ 作为椭圆盘 $\Sigma$（平面 $z=y+3$ 上 $x^2+y^2 \leq 1$ 部分）的边界，也可以作为其他曲面的边界（比如圆柱面的侧面从 $z=0$ 到 $z=y+3$ 的部分与底面 $z=0$ 的组合）。

但是根据 Stokes 公式的边界决定性质，无论选择哪一个以 $L$ 为边界的曲面，$\iint_\Sigma (\nabla\times\mathbf{F})\cdot\mathbf{n}\,dS$ 的值都等于 $\oint_L \mathbf{F}\cdot d\mathbf{r}$，因此与曲面的选择无关。

这就是为什么我们可以自由选择**最简单**的以 $L$ 为边界的曲面——在例 2 中，椭圆盘是最简单的选择（平面，计算 $dS$ 最方便）。

**条件提醒**：这一性质成立的前提是 $\nabla\times\mathbf{F}$ 在两个曲面所夹的区域上处处有定义且连续。如果场在某个曲面附近有奇点，则不能任意选择曲面。

**实际意义**：在应用 Stokes 公式时，我们总是选择**最容易计算曲面积分**的有向曲面。通常的选择策略是：
- 优先选择平面（$dS$ 计算最简单）
- 其次选择与坐标平面平行的平面（投影法和 $dS$ 最简化）
- 避免选择弯曲复杂的曲面（除非必要）

</details>