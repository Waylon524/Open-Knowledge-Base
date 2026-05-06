# 03. Green公式（格林公式）

> 所属章节：第十四章 曲线积分、曲面积分与场论  |  文件序号：03  |  难度：基础
> 常见混淆点：将正向边界误认为"逆时针"而忽略"区域始终在左侧"的完整含义——对于带洞区域（复连通），外边界逆时针、内边界顺时针；混淆Green公式的旋度形式与散度形式——$\oint_C P\,dx+Q\,dy = \iint_D (\partial Q/\partial x-\partial P/\partial y)\,dA$ 描述的是环量（旋度），而非通量

## 1. 学习目标与先修前置

### 学习目标
- 掌握Green公式（Green's Theorem）的完整陈述：$\displaystyle \oint_{\partial D} P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$
- 理解正向边界（positive orientation）的准确定义——沿边界行走时区域 $D$ 始终在左侧
- 理解Green公式的条件要求（$C$ 分段光滑简单闭曲线，$P,Q$ 有连续偏导数），能够判断公式是否适用
- 掌握利用Green公式将封闭曲线上的第二型曲线积分转化为二重积分的计算方法
- 掌握从Green公式推导平面区域面积公式 $A = \dfrac{1}{2}\oint_C x\,dy - y\,dx$ 并用于计算椭圆等规则图形的面积
- 理解Green公式与第二型曲线积分方向反转性质的联系——顺时针边界如何修正公式
- 了解Green公式的物理含义：环量等于旋度在区域上的累积

### 先修知识
- 文件02（第十四章）：第二型曲线积分的定义（定义14.3）、参数化转化公式（定理14.4）、方向反转变号性质——Green公式建立了第二型闭曲线积分与二重积分之间的桥梁
- 文件01（第十四章）：第一型曲线积分的概念与参数化方法——对比理解两类曲线积分的区别
- 文件01-03（第十三章）：二重积分的计算（X-型区域/Y-型区域的积分限确定、极坐标变换）——Green公式的右端是二重积分
- 文件03（第十三章）：二重积分的极坐标变换——用于圆形区域上的Green公式计算
- 文件03（第十二章）：偏导数的计算——$\partial Q/\partial x$ 和 $\partial P/\partial y$ 的计算是应用Green公式的前提

---

## 2. 背景与应用场景

### 2.1 从第二型曲线积分到二重积分的转化

在文件02（第十四章）中，我们学习了第二型曲线积分 $\int_L P\,dx + Q\,dy$ 的计算方法——通过参数化转化为定积分。这种方法对任意曲线（开放或封闭）都适用，但计算量较大，且每次都要写出参数方程并代入。

在自然科学和工程中，我们经常遇到**封闭曲线**上的第二型曲线积分：

$$\oint_C P\,dx + Q\,dy$$

其中 $C$ 是一条封闭的有向曲线，$\oint$ 表示沿封闭曲线的积分（circle indicates closed curve）。例如：

- **环量（circulation）**：流体沿封闭回路的速度场积分 $\oint_C \mathbf{v} \cdot d\mathbf{r}$，描述流体在该回路中的旋转强度
- **电动势（EMF）**：电场沿闭合回路的积分 $\oint_C \mathbf{E} \cdot d\mathbf{r}$，在静电场中为零，在变化磁场中不为零（法拉第电磁感应定律）
- **功的净效果**：在保守力场中 $\oint_C \mathbf{F} \cdot d\mathbf{r} = 0$，但在非保守力场中通常不为零

**核心问题**：对封闭曲线上的第二型曲线积分，是否存在一种方法将其转化为该曲线所围区域上的二重积分，从而简化计算？

答案是肯定的——这就是**Green公式**。

### 2.2 Green公式的物理直觉

考虑一个向量场 $\mathbf{F} = (P, Q)$，沿一条很小的封闭曲线 $C$ 的环量 $\oint_C P\,dx + Q\,dy$ 描述了"场绕这个回路旋转"的强度。

在回路内部，场的旋转（旋度）由 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$ 度量。Green公式的核心思想是：

$$\text{边界上的总环量} = \text{内部所有旋度的总和}$$

这类似于一维微积分基本定理：

$$\int_a^b f'(x)\,dx = f(b) - f(a)$$

在一维中，区间 $[a,b]$ 上的导数累积等于边界点（端点）上的函数差值。在二维中，区域 $D$ 上的旋度（导数）累积等于边界 $\partial D$ 上的环量。因此，Green公式可以看作**微积分基本定理在二维的推广**。

### 2.3 从路径依赖性到Green公式的动机

文件02（第十四章）的例1中，我们观察到沿不同路径从 $(0,0)$ 到 $(1,1)$ 的线积分结果不同（$\frac13$ 和 $\frac12$）。但对**封闭**曲线（起点=终点），如果线积分的值不为零，说明场不是保守的，而非零的闭曲线积分值是衡量场"非保守程度"的重要指标。

Green公式提供了计算闭曲线积分的有力工具，同时也为判断保守力场（$\partial P/\partial y = \partial Q/\partial x$）提供了理论基础——当 $\partial Q/\partial x - \partial P/\partial y = 0$ 时，闭曲线积分恒为零，场是保守的。

### 2.4 Green公式的应用领域

| 应用领域 | 具体问题 | Green公式的作用 |
|:---------|:---------|:----------------|
| 力学 | 计算力沿闭合路径所做的净功 | 转换为区域内角动量的累积 |
| 流体力学 | 计算速度场的环量 | 转换为涡量的面积分 |
| 电磁学 | 计算电动势 | 法拉第定律的数学基础 |
| 几何学 | 计算不规则图形的面积 | $A = \frac12\oint_C x\,dy - y\,dx$ |
| 大地测量 | 计算多边形区域的面积 | 离散化的Green公式（鞋带公式） |

---

## 3. 核心概念与符号约定

### 3.1 简单闭曲线与正向边界

在正式陈述Green公式之前，我们需要精确定义两类关键概念。

**定义 14.5（简单闭曲线 / Simple Closed Curve）**：若曲线 $C$ 的参数方程 $\mathbf{r}(t),\; a \leq t \leq b$ 满足：
- $\mathbf{r}(a) = \mathbf{r}(b)$（闭合，起点等于终点）
- 对任意 $t_1, t_2 \in (a, b)$，若 $t_1 \neq t_2$ 则 $\mathbf{r}(t_1) \neq \mathbf{r}(t_2)$（不自交）
则称 $C$ 是**简单闭曲线**（simple closed curve）。

简单闭曲线将平面划分为三个互不相交的部分：内部（有界区域）、外部（无界区域）和曲线本身（Jordan曲线定理，此处不证明）。

**定义 14.6（正向边界 / Positive Orientation）**：设 $C$ 是区域 $D$ 的边界曲线 $C = \partial D$。若沿 $C$ 的指定方向行走时，区域 $D$ **始终位于行走方向的左侧**，则称该方向为 $C$ 的**正向**（positive orientation），记作 $C^+$ 或简写为 $C$（默认即正向）。

**直观理解**：
- 对于**无洞区域**（单连通区域），正向就是**逆时针方向**（counterclockwise）
- 对于**有洞区域**（复连通区域），外边界取逆时针（正向），内边界（洞的边界）取**顺时针**（也是正向——因为行走时区域在左侧）

**为什么是"左侧"而非"逆时针"**：
"左侧规则"是正确定义，而"逆时针"只在单连通无洞区域的特例中与左侧规则一致。对于复连通区域，如果简单地认为正向=逆时针，会在内边界上出错。

**本节的约定**：在本节中，除非特别说明，所有曲线均取正向。

### 3.2 Green公式的正式陈述

**定理 14.7（Green公式 / Green's Theorem）**：设 $D \subset \mathbb{R}^2$ 为有界闭区域，其边界 $\partial D$ 由有限条分段光滑的简单闭曲线组成，且取**正向**（即沿边界行走时区域 $D$ 始终在左侧）。若函数 $P(x,y)$ 和 $Q(x,y)$ 在包含 $D$ 的某个开区域上具有**连续的偏导数**，则：

$$\boxed{\oint_{\partial D} P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA}$$

其中：
- $\partial D$ — 区域 $D$ 的**正向边界**
- $dA$ — 面积微元（可写为 $dx\,dy$ 或 $d\sigma$）
- $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$ — **旋度**（curl）的二维版本，也称**标量旋度**

**公式的另一种常见写法**（向量形式）：

$$\oint_{\partial D} \mathbf{F} \cdot d\mathbf{r} = \iint_D (\nabla \times \mathbf{F}) \cdot \mathbf{k}\,dA$$

其中 $\mathbf{F} = (P, Q)$，$\nabla \times \mathbf{F} = \left(0, 0, \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right)$，$\mathbf{k} = (0,0,1)$。

### 3.3 条件分析：为什么需要这些条件？

Green公式的三个核心条件各有其必要性：

**条件 1：$D$ 是有界闭区域，$\partial D$ 是分段光滑简单闭曲线**

- **有界**：区域面积有限，二重积分有意义
- **闭区域**：包含边界，边界上的点属于 $D$
- **分段光滑**：边界可以包含角点（如三角形的顶点、矩形的角），但每条分段光滑
- **简单闭曲线**：不自交，确保 $D$ 有明确的"内部"

**反例**：若边界自交（如"8"字形曲线），自交点处区域 $D$ 的定义出现歧义——曲线两侧的区域不再是简单分立的"内部"和"外部"，Green公式不能直接应用。

**条件 2：$C = \partial D$ 取正向**

- 若曲线取反向（顺时针），公式变为：
  $$\oint_{\partial D^-} P\,dx + Q\,dy = -\iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$$
  这是因为第二型曲线积分方向反转变号（文件02第十四章，3.2节），而右端二重积分与边界方向无关。

**条件 3：$P$ 和 $Q$ 在包含 $D$ 的开区域上有连续偏导数**

- **连续偏导数存在**是二重积分 $\iint_D (\partial Q/\partial x - \partial P/\partial y)\,dA$ 有定义的保证
- 若 $P$ 或 $Q$ 在 $D$ 内某点处无定义或偏导数不连续（例如分母为零的点），则需要将该点挖掉，对复连通区域应用Green公式

### 3.4 旋度的直观含义

表达式 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$ 在物理学中称为向量场 $(P,Q)$ 的**标量旋度**（scalar curl），也称为**二维旋度**。

**物理直觉**：在场中放一个很小的"风车"（测试轮），其旋转趋势由旋度度量：
- 旋度 $> 0$：场有逆时针旋转的趋势
- 旋度 $< 0$：场有顺时针旋转的趋势
- 旋度 $= 0$：场无旋转（保守场/无旋场）

**与Green公式的关系**：旋度在区域上的累积等于边界上的环量——即"总旋转 = 边界上的净环量"。

### 3.5 平面区域面积公式的推导

利用Green公式可以推导出一个用线积分计算平面区域面积的重要公式。

设区域 $D$ 的面积为 $A(D)$。我们希望找到 $P, Q$ 使得 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y} \equiv 1$，这样由Green公式：

$$A(D) = \iint_D 1\,dA = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA = \oint_{\partial D} P\,dx + Q\,dy$$

有多种选择满足 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y} = 1$：

| 选择 | $P$ | $Q$ | $\partial Q/\partial x - \partial P/\partial y$ | 面积公式 |
|:-----|:----|:----|:----------------------------------------------|:---------|
| 1 | $-y$ | $0$ | $0 - (-1) = 1$ | $A = -\oint_C y\,dx$ |
| 2 | $0$ | $x$ | $1 - 0 = 1$ | $A = \oint_C x\,dy$ |
| 3 | $-\frac{y}{2}$ | $\frac{x}{2}$ | $\frac12 - (-\frac12) = 1$ | $A = \frac12\oint_C x\,dy - y\,dx$ |

**最常用的对称形式**：取 $P = -\dfrac{y}{2}$，$Q = \dfrac{x}{2}$，得到：

$$\boxed{A(D) = \frac{1}{2} \oint_{\partial D} x\,dy - y\,dx}$$

**验证**：此时 $\partial Q/\partial x = 1/2$，$\partial P/\partial y = -1/2$，差值为 $1/2 - (-1/2) = 1$，代入Green公式即得面积。

**物理含义**：$x\,dy - y\,dx = \left| \begin{array}{cc} x & y \\ dx & dy \end{array} \right|$ 是位置向量 $\mathbf{r} = (x,y)$ 与有向线微元 $d\mathbf{r} = (dx,dy)$ 所张平行四边形的有向面积（的2倍）。沿封闭曲线积分就得到整个区域的总面积。

### 3.6 符号表

| 符号 | 含义 | 说明 |
|:-----|:-----|:-----|
| $\oint_C$ | 沿封闭曲线 $C$ 的积分 | 表示起点=终点的线积分 |
| $\partial D$ | 区域 $D$ 的边界 | 取正向（区域在左侧） |
| $\partial D^-$ | 区域 $D$ 的反向边界 | 取反向（区域在右侧） |
| $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$ | 二维向量场的标量旋度 | 度量场的局部旋转强度 |
| $C^+$ | 曲线 $C$ 的正向 | 单连通区域下=逆时针方向 |
| $C^-$ | 曲线 $C$ 的反向 | 与正相反的向 |
| $A(D)$ | 平面区域 $D$ 的面积 | 面积公式推导的结果 |
| $\mathbf{k}$ | $z$ 轴正方向的单位向量 | $(0,0,1)$ |

---

## 4. 原理与方法

### 4.1 Green公式的证明思路（对X-型区域的简单情形）

虽然完整的证明需要处理一般区域，但我们可以通过X-型区域和Y-型区域的特殊情形来理解核心思想。这个证明也展示了如何将边界线积分转化为区域上的二重积分。

**情形 1：$D$ 为X-型区域**

设 $D$ 可表示为：
$$D = \{(x,y) \mid a \leq x \leq b,\; y_1(x) \leq y \leq y_2(x)\}$$

其中 $y_1(x), y_2(x)$ 在 $[a,b]$ 上连续，且 $y_1(x) \leq y_2(x)$。

考虑只含有 $P\,dx$ 项的Green公式（即先证 $\oint_{\partial D} P\,dx = -\iint_D \frac{\partial P}{\partial y}\,dA$）。

**第 1 步：计算右端二重积分**

$$\begin{aligned}
-\iint_D \frac{\partial P}{\partial y}\,dA &= -\int_a^b \int_{y_1(x)}^{y_2(x)} \frac{\partial P}{\partial y}(x,y)\,dy\,dx \\
&= -\int_a^b \left[ P(x,y) \right]_{y=y_1(x)}^{y=y_2(x)} dx \quad (\text{对 } y \text{ 积分，由微积分基本定理}) \\
&= -\int_a^b \big[ P(x, y_2(x)) - P(x, y_1(x)) \big]\,dx \\
&= \int_a^b P(x, y_1(x))\,dx - \int_a^b P(x, y_2(x))\,dx
\end{aligned}$$

**第 2 步：计算左端线积分**

边界 $\partial D$ 由四部分组成：
- 下边界 $C_1$：$y = y_1(x)$，$x$ 从 $a$ 到 $b$（方向：左→右）
- 右边界 $C_2$：$x = b$，$y$ 从 $y_1(b)$ 到 $y_2(b)$（方向：下→上）
- 上边界 $C_3$：$y = y_2(x)$，$x$ 从 $b$ 到 $a$（方向：右→左，与 $C_1$ 相反）
- 左边界 $C_4$：$x = a$，$y$ 从 $y_2(a)$ 到 $y_1(a)$（方向：上→下）

沿正向边界 $\partial D$ 的线积分 $\oint P\,dx$ 中，在 $C_2$ 和 $C_4$ 上 $dx = 0$（因为 $x$ 为常数），所以：
$$\oint_{\partial D} P\,dx = \int_{C_1} P\,dx + \int_{C_3} P\,dx$$

计算 $C_1$（$x$ 从 $a$ 到 $b$，$y = y_1(x)$，$dx = dx$）：
$$\int_{C_1} P\,dx = \int_a^b P(x, y_1(x))\,dx$$

计算 $C_3$（$x$ 从 $b$ 到 $a$，$y = y_2(x)$，$dx = dx$）：
$$\int_{C_3} P\,dx = \int_b^a P(x, y_2(x))\,dx = -\int_a^b P(x, y_2(x))\,dx$$

因此：
$$\oint_{\partial D} P\,dx = \int_a^b P(x, y_1(x))\,dx - \int_a^b P(x, y_2(x))\,dx$$

**第 3 步：比较左右两端**

左端 $\oint_{\partial D} P\,dx$ 的表达式与右端 $-\iint_D \frac{\partial P}{\partial y}\,dA$ 的表达式完全相同。因此：
$$\oint_{\partial D} P\,dx = -\iint_D \frac{\partial P}{\partial y}\,dA$$

**情形 2：$D$ 为Y-型区域**

类似地，将 $D$ 表示为 $D = \{(x,y) \mid c \leq y \leq d,\; x_1(y) \leq x \leq x_2(y)\}$，可以证明：
$$\oint_{\partial D} Q\,dy = \iint_D \frac{\partial Q}{\partial x}\,dA$$

**合并**：将两个等式相加，得到：
$$\oint_{\partial D} (P\,dx + Q\,dy) = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$$

**一般区域的推广**：对于不是同时满足X-型和Y-型的一般区域，可以将其分割为有限个X-型或Y-型子区域的并，在每个子区域上应用上述证明，然后利用边界线积分的可加性（相邻边界方向相反，相互抵消），最终得到整体区域上的Green公式。

### 4.2 Green公式应用的标准工作流

**第 1 步：确认条件**
- 积分曲线 $C$ 是封闭的简单曲线吗？（若未封闭，不能直接使用Green公式）
- $C$ 是分段光滑的吗？（有角点不影响，但需分段处理）
- $C$ 取正向了吗？（若为顺时针，需要添加负号）
- $P, Q$ 在 $D$ 内有连续的偏导数吗？（若有奇点，需挖洞处理）

**第 2 步：识别 $P$ 和 $Q$**
- 从被积表达式 $P\,dx + Q\,dy$ 中读出 $P(x,y)$ 和 $Q(x,y)$
- 注意：$dx$ 的系数是 $P$，$dy$ 的系数是 $Q$

**第 3 步：计算旋度**
- 计算 $\dfrac{\partial Q}{\partial x}$ 和 $\dfrac{\partial P}{\partial y}$
- 计算 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$

**第 4 步：转化为二重积分**
- 写出 $\displaystyle \oint_C P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$
- 确定区域 $D$（$C$ 所围成的区域）
- 选择适当的积分方法（X-型/Y-型积分限确定，或极坐标变换）

**第 5 步：计算二重积分**
- 用第十三章的方法计算二重积分

### 4.3 方向修正：顺时针曲线的处理

若曲线 $C$ 取**顺时针**方向（记为 $C^-$），则不能直接套用Green公式。处理方式有两种：

**方法一（直接加负号）**：
$$\oint_{C^-} P\,dx + Q\,dy = -\iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$$

理由：第二型曲线积分方向反转变号（文件02第十四章，3.2节），而正向曲线 $C$ 的Green公式为 $\oint_C P\,dx + Q\,dy = \iint_D (\partial Q/\partial x - \partial P/\partial y)\,dA$，所以 $C^-$ 上的积分等于正向积分的相反数。

**方法二（直接使用反向参数化）**：
将曲线参数化时取 $t$ 递增方向对应顺时针，然后直接代入定理14.4计算——这时Green公式中的右端会自然产生一个负号，因为$dx, dy$的符号发生了变化。

**实用建议**：始终先确认曲线的方向。如果曲线是逆时针（正向），直接使用Green公式；如果曲线是顺时针，在公式右端添加负号。

### 4.4 Green公式与保守力场的联系

在文件02（第十四章）的迁移应用练习题5中，我们发现力场 $\mathbf{F} = (y, x)$ 沿不同路径的积分值相同，而 $\mathbf{F} = (0, x^2)$ 则不同。

利用Green公式可以给出路径无关性的充要条件：

**定理（保守力场的旋度条件）**：设 $D$ 是单连通区域，$\mathbf{F} = (P,Q)$ 在 $D$ 上有连续偏导数。则 $\int_L \mathbf{F} \cdot d\mathbf{r}$ 在 $D$ 内与路径无关当且仅当：
$$\frac{\partial P}{\partial y} = \frac{\partial Q}{\partial x} \quad \text{在 } D \text{ 内处处成立}$$
或等价地，标量旋度为零：
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 0$$

**解释**：若旋度处处为零，则由Green公式，对 $D$ 内任何闭曲线 $C$，$\oint_C P\,dx + Q\,dy = 0$。而路径无关等价于任意闭曲线积分为零（文件02第十四章中已讨论）。在本节中，我们仅指出这一联系——保守力场的完整理论将在后续章节（如势函数、全微分条件）中系统展开。

---

## 5. 例题

### 例 1：三角形区域上的Green公式应用（对应考试题 Q1）

设曲线 $C$ 是以 $A(0,0)$、$B(1,0)$、$C(1,1)$ 为顶点的三角形的**正向**边界。计算曲线积分：

$$I = \oint_C (xy)\,dx + (x^2)\,dy$$

**解**：

**第 1 步：确认条件和方向**

- 曲线 $C$ 是三角形边界——分段光滑简单闭曲线。$\checkmark$
- 正向：三角形边界取逆时针方向（$A \to B \to C \to A$）。验证：沿此方向行走时三角形区域在左侧。$\checkmark$
- 被积函数 $P(x,y) = xy$，$Q(x,y) = x^2$ 在全平面有连续偏导数。$\checkmark$
- 条件满足，可以应用Green公式。

**第 2 步：计算旋度**

$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}(x^2) = 2x$$
$$\frac{\partial P}{\partial y} = \frac{\partial}{\partial y}(xy) = x$$

因此：
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 2x - x = x$$

**第 3 步：确定区域 $D$ 并建立二重积分**

$D$ 是顶点为 $A(0,0), B(1,0), C(1,1)$ 的三角形区域。由Green公式：

$$I = \iint_D x\,dx\,dy$$

**第 4 步：化为累次积分并计算**

$D$ 是X-型区域：
$$0 \leq x \leq 1,\quad 0 \leq y \leq x$$

$$\begin{aligned}
I &= \iint_D x\,dx\,dy \\
  &= \int_{0}^{1} \int_{0}^{x} x\,dy\,dx \\
  &= \int_{0}^{1} x \left( \int_{0}^{x} dy \right) dx \\
  &= \int_{0}^{1} x \cdot \big[ y \big]_{y=0}^{y=x} \,dx \\
  &= \int_{0}^{1} x \cdot x\,dx \\
  &= \int_{0}^{1} x^2\,dx \\
  &= \left[ \frac{x^3}{3} \right]_{0}^{1} \\
  &= \frac{1}{3}
\end{aligned}$$

因此 $\displaystyle I = \frac{1}{3}$。

**第 5 步：验证（方向反转检查）**

若 $C$ 取顺时针方向（$A \to C \to B \to A$），则结果为 $I_{C^-} = -\dfrac{1}{3}$。这是因为：
- 由第二型曲线积分方向反转变号性质（文件02第十四章）：$\oint_{C^-} = -\oint_C$
- 直接代入修正后的Green公式：$\oint_{C^-} P\,dx + Q\,dy = -\iint_D x\,dA = -\dfrac{1}{3}$

两种推理结果一致。$\checkmark$

**思考**：如果直接用参数化方法沿三条边分别计算，结果也是 $\frac13$ 吗？读者可以自行验证——这相当于同时验证了Green公式的正确性。

---

### 例 2：用线积分计算椭圆面积（对应考试题 Q2）

**(1)** 取 $P = -y$，$Q = x$，计算 $\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}$。

**解**：

$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}(x) = 1,\qquad \frac{\partial P}{\partial y} = \frac{\partial}{\partial y}(-y) = -1$$

$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 1 - (-1) = 2$$

---

**(2)** 利用Green公式，将 $\displaystyle \oint_C -y\,dx + x\,dy$ 转化为二重积分，并说明与 $C$ 所围区域 $D$ 的面积之间的关系。

**解**：

由Green公式：

$$\oint_C -y\,dx + x\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA = \iint_D 2\,dA = 2 \times A(D)$$

其中 $D$ 为 $C$ 所围成的区域。因此：

$$\oint_C -y\,dx + x\,dy = 2 \cdot A(D)$$

即该线积分等于 $C$ 所围区域面积的 2 倍。

---

**(3)** 由 (2) 推导出用线积分表示平面区域面积的公式。

**解**：

由 (2) 的结果 $2A(D) = \oint_C -y\,dx + x\,dy$，解得：

$$A(D) = \frac{1}{2} \oint_C (-y)\,dx + x\,dy = \frac{1}{2} \oint_C x\,dy - y\,dx$$

这就是用线积分计算平面区域面积的公式。注意这与我们在 3.5 节中选取 $P = -y/2,\; Q = x/2$ 得到的公式 $A = \frac12 \oint x\,dy - y\,dx$ 完全一致（因为 $(-y)\,dx + x\,dy = x\,dy - y\,dx$）。

---

**(4)** 参数化椭圆 $\dfrac{x^2}{4} + \dfrac{y^2}{9} = 1$（取正向），利用上题公式计算该椭圆所围成的面积。

**解**：

**第 1 步：参数化椭圆（正向——逆时针）**

椭圆的标准参数方程：
$$\begin{cases}
x = 2\cos t \\
y = 3\sin t
\end{cases},\quad 0 \leq t \leq 2\pi$$

验证正向：$t=0$ 时 $(x,y) = (2,0)$；$t$ 增加时 $y = 3\sin t > 0$（对小的 $t > 0$），点沿逆时针方向运动。$\checkmark$

**第 2 步：计算有向投影微元**

$$dx = -2\sin t\,dt,\qquad dy = 3\cos t\,dt$$

**第 3 步：代入面积公式**

$$\begin{aligned}
A &= \frac{1}{2} \oint_C x\,dy - y\,dx \\
  &= \frac{1}{2} \int_{0}^{2\pi} \big[ (2\cos t)(3\cos t\,dt) - (3\sin t)(-2\sin t\,dt) \big] \\
  &= \frac{1}{2} \int_{0}^{2\pi} (6\cos^2 t + 6\sin^2 t)\,dt \\
  &= \frac{1}{2} \int_{0}^{2\pi} 6(\cos^2 t + \sin^2 t)\,dt \\
  &= \frac{1}{2} \int_{0}^{2\pi} 6\,dt \\
  &= \frac{1}{2} \cdot 6 \cdot 2\pi \\
  &= 6\pi
\end{aligned}$$

**第 4 步：验证**

椭圆面积公式为 $A = \pi ab$，其中 $a$ 为半长轴，$b$ 为半短轴。本题中 $a = 2$（$x$ 方向半轴），$b = 3$（$y$ 方向半轴），所以：
$$A = \pi \cdot 2 \cdot 3 = 6\pi$$

与线积分计算结果一致。$\checkmark$

**关键要点**：面积公式 $A = \frac12 \oint_C x\,dy - y\,dx$ 适用于任何用正向简单闭曲线围成的区域，无论其形状是否规则。对于像椭圆这样的规则图形，参数化代入计算即可；对于多边形区域，可以离散化为鞋带公式。

---

### 例 3：圆周上的Green公式与参数化验证（对应考试题 Q3）

设 $P(x,y) = y$，$Q(x,y) = x^2$，曲线 $C$ 为正向圆周 $x^2 + y^2 = R^2$（$R > 0$）。

计算 $\displaystyle I = \oint_C y\,dx + x^2\,dy$：

**(1)** 利用Green公式计算；

**(2)** 直接参数化圆周逐项计算，验证结果与 (1) 一致；

**(3)** 若 $C$ 改为顺时针，结果如何变化？

---

#### (1) 利用Green公式计算

**第 1 步：确认条件**

- $C$ 是半径为 $R$ 的圆周——分段光滑简单闭曲线。$\checkmark$
- $C$ 取正向（逆时针）。$\checkmark$
- $P(x,y) = y$，$Q(x,y) = x^2$ 在全平面有连续偏导数。$\checkmark$
- 条件满足。

**第 2 步：计算旋度**

$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}(x^2) = 2x,\qquad \frac{\partial P}{\partial y} = \frac{\partial}{\partial y}(y) = 1$$

$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 2x - 1$$

**第 3 步：转化为二重积分**

由Green公式：
$$I = \iint_D (2x - 1)\,dA$$

其中 $D$ 是圆盘 $x^2 + y^2 \leq R^2$。

**第 4 步：用极坐标变换计算二重积分**

采用极坐标变换（参见第十三章文件03）：
$$\begin{cases}
x = r\cos\theta \\
y = r\sin\theta
\end{cases},\quad 0 \leq r \leq R,\quad 0 \leq \theta \leq 2\pi,\quad dA = r\,dr\,d\theta$$

$$\begin{aligned}
I &= \int_{0}^{2\pi} \int_{0}^{R} (2r\cos\theta - 1)\, r\,dr\,d\theta \\
  &= \int_{0}^{2\pi} \int_{0}^{R} (2r^2\cos\theta - r)\,dr\,d\theta \\
  &= \int_{0}^{2\pi} \left[ \frac{2r^3}{3}\cos\theta - \frac{r^2}{2} \right]_{r=0}^{r=R} d\theta \\
  &= \int_{0}^{2\pi} \left( \frac{2R^3}{3}\cos\theta - \frac{R^2}{2} \right) d\theta \\
  &= \left[ \frac{2R^3}{3}\sin\theta - \frac{R^2}{2}\theta \right]_{0}^{2\pi} \\
  &= \left( \frac{2R^3}{3}\cdot 0 - \frac{R^2}{2}\cdot 2\pi \right) - \left( \frac{2R^3}{3}\cdot 0 - \frac{R^2}{2}\cdot 0 \right) \\
  &= -\pi R^2
\end{aligned}$$

因此 $\displaystyle I = -\pi R^2$。

**结果分析**：$I$ 为负值，说明向量场 $\mathbf{F} = (y, x^2)$ 沿逆时针圆周的环量为负——场整体有顺时针旋转的趋势。

---

#### (2) 直接参数化验证

**第 1 步：参数化圆周（正向——逆时针）**

$$\begin{cases}
x = R\cos t \\
y = R\sin t
\end{cases},\quad 0 \leq t \leq 2\pi$$

验证正向：$t=0$ 时 $(x,y)=(R,0)$；$t$ 略大于 $0$ 时 $y = R\sin t > 0$，点向上运动——逆时针方向。$\checkmark$

**第 2 步：计算有向投影微元**

$$dx = -R\sin t\,dt,\qquad dy = R\cos t\,dt$$

**第 3 步：代入被积函数并展开**

$$\begin{aligned}
I &= \oint_C y\,dx + x^2\,dy \\
  &= \int_{0}^{2\pi} \big[ (R\sin t)(-R\sin t) + (R^2\cos^2 t)(R\cos t) \big]\,dt \\
  &= \int_{0}^{2\pi} \big( -R^2\sin^2 t + R^3\cos^3 t \big)\,dt \\
  &= -R^2 \int_{0}^{2\pi} \sin^2 t\,dt + R^3 \int_{0}^{2\pi} \cos^3 t\,dt
\end{aligned}$$

**第 4 步：计算两个积分**

**积分一：** $\displaystyle \int_0^{2\pi} \sin^2 t\,dt$

利用降幂公式 $\sin^2 t = \dfrac{1 - \cos 2t}{2}$：

$$\begin{aligned}
\int_0^{2\pi} \sin^2 t\,dt &= \int_0^{2\pi} \frac{1 - \cos 2t}{2}\,dt \\
&= \frac{1}{2} \int_0^{2\pi} 1\,dt - \frac{1}{2} \int_0^{2\pi} \cos 2t\,dt \\
&= \frac{1}{2} \cdot 2\pi - \frac{1}{2} \cdot 0 \quad (\text{余弦函数在一个完整周期上的积分为零}) \\
&= \pi
\end{aligned}$$

**积分二：** $\displaystyle \int_0^{2\pi} \cos^3 t\,dt$

利用 $\cos^3 t = \cos t\,(1 - \sin^2 t)$，令 $u = \sin t$，则 $du = \cos t\,dt$。当 $t$ 从 $0$ 到 $2\pi$ 时，$u$ 从 $0$ 变回 $0$（$\sin t$ 在一个完整周期内从 $0 \to 1 \to 0 \to -1 \to 0$），因此周期积分为零。更直接的方法：$\cos^3 t$ 是偶函数，在 $[0,2\pi]$ 上对称区间的积分值为零。

$$\int_0^{2\pi} \cos^3 t\,dt = 0$$

**第 5 步：合并结果**

$$\begin{aligned}
I &= -R^2 \cdot \pi + R^3 \cdot 0 \\
  &= -\pi R^2
\end{aligned}$$

**结论**：直接参数化计算的结果为 $-\pi R^2$，与 (1) 中利用Green公式所得结果完全一致。验证了Green公式的正确性。$\checkmark$

---

#### (3) 顺时针方向的影响

若将 $C$ 改为顺时针方向（记为 $C^-$），结果如何？

**分析**：由第二型曲线积分的方向反转变号性质（文件02第十四章，3.2节），当曲线方向反转时，积分值变号：

$$\oint_{C^-} y\,dx + x^2\,dy = -\oint_{C} y\,dx + x^2\,dy = -(-\pi R^2) = \pi R^2$$

**与Green公式的一致性**：Green公式要求曲线取**正向**（逆时针）。对于顺时针曲线，修正后的公式为：

$$\oint_{C^-} P\,dx + Q\,dy = -\iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA$$

代入 $\partial Q/\partial x - \partial P/\partial y = 2x - 1$：

$$\oint_{C^-} y\,dx + x^2\,dy = -\iint_D (2x - 1)\,dA = -(-\pi R^2) = \pi R^2$$

两种方法结果一致。$\checkmark$

**关键理解**：Green公式与方向反转性质不矛盾——前者提供了闭曲线积分的计算工具，后者描述了方向变化时积分的变换规律。两者结合使用可以处理任意方向的闭曲线积分。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|:---------|:---------|
| 正向边界就是逆时针方向，对任何区域都成立 | 正向边界的准确定义是"区域始终在左侧"。对无洞区域确实等价于逆时针，但对有洞的复连通区域，外边界逆时针、内边界顺时针才是正向 |
| Green公式中右端的二重积分被积函数是 $\partial P/\partial y - \partial Q/\partial x$ | 正确形式是 $\partial Q/\partial x - \partial P/\partial y$。顺序不能颠倒——线积分中 $P$ 与 $dx$ 配对，$Q$ 与 $dy$ 配对，偏导差是 $Q$ 对 $x$ 减 $P$ 对 $y$ |
| 任何封闭曲线都可用Green公式 | 需要满足条件：(1) 简单闭曲线（不自交）；(2) 分段光滑；(3) $P,Q$ 在区域内连续偏导。若 $D$ 内有奇点（如分母为零的点），需要挖洞处理 |
| Green公式只能用于正向曲线，反向曲线不能用 | 反向曲线可以用修正公式 $\oint_{C^-} = -\iint_D(\cdots)\,dA$，或者先按正向计算再取负 |
| 面积公式 $A = \frac12\oint x\,dy - y\,dx$ 对所有曲线都成立 | 公式要求 $C$ 取正向（逆时针）。若 $C$ 为顺时针，$A = -\frac12\oint x\,dy - y\,dx$（因为积分值变号） |
| Green公式的右端二重积分可以任意选择积分次序 | 二重积分的计算确实可以用任意次序，但需正确确定积分限。用X-型或Y-型描述区域时要注意上下限对应正确的边界函数 |
| 只有 $P$ 和 $Q$ 有定义就能用Green公式 | 不仅要有定义，还需要 $P$ 和 $Q$ 有**连续的偏导数**。如果偏导数不连续（如 $P = y/(x^2+y^2)$ 在原点无定义），不能直接应用 |

### 检查点

- [ ] 能否用自己的语言完整叙述Green公式，包括所有条件要求（$C$ 简单闭曲线、分段光滑、正向、$P,Q$ 有连续偏导数）？
- [ ] 能否准确写出 $\partial Q/\partial x - \partial P/\partial y$ 的表达式并解释其物理含义（二维旋度）？
- [ ] 给定一个具体的闭曲线和向量场 $(P,Q)$，能否按照五步工作流应用Green公式？
- [ ] 能否说明正向边界的完整定义（"区域在左侧"而非仅仅"逆时针"）？
- [ ] 能否从Green公式推导出面积公式 $A = \frac12\oint x\,dy - y\,dx$？
- [ ] 能否用面积公式计算椭圆、圆等规则图形的面积？
- [ ] 如果曲线取顺时针方向，Green公式如何修正？
- [ ] 能否通过直接参数化验证Green公式的计算结果（如例3）？
- [ ] 能否说明Green公式与文件02中方向反转性质的关系？
- [ ] 能否识别什么时候不能使用Green公式（边界自交、有奇点等）？

---

## 练习题

### 基础巩固

**1.**（Green公式的条件判断）判断下列曲线积分能否直接用Green公式（将闭曲线积分转化为二重积分）。若能，写出 $\partial Q/\partial x - \partial P/\partial y$；若不能，说明原因。

**(1)** $\displaystyle \oint_C xy\,dx + (x+y)\,dy$，$C$ 是以 $(0,0), (2,0), (0,2)$ 为顶点的三角形正向边界

**(2)** $\displaystyle \oint_C \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy$，$C$ 是圆周 $x^2+y^2 = 1$（正向）

**(3)** $\displaystyle \int_L y\,dx - x\,dy$，$L$ 是从 $(1,0)$ 到 $(0,1)$ 的直线段（不是封闭曲线）

<details><summary>参考答案</summary>

**(1)** 可以。$P = xy$，$Q = x + y$。
$$\frac{\partial Q}{\partial x} = 1,\quad \frac{\partial P}{\partial y} = x$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 1 - x$$
三角形顶点 $(0,0), (2,0), (0,2)$，正向即为逆时针。$\checkmark$

**(2)** **不能直接应用**。虽然 $C$ 是简单闭曲线，但 $P = \dfrac{-y}{x^2+y^2}$ 和 $Q = \dfrac{x}{x^2+y^2}$ 在原点 $(0,0)$ 无定义（分母为零），而原点在 $C$ 所围的圆盘 $x^2+y^2 \leq 1$ 内部。因此 $P, Q$ 在包含 $D$ 的开区域上不是连续可微的，不满足Green公式的条件。需要将原点挖掉（挖去一个小圆盘），然后对复连通区域使用Green公式（涉及两个边界的线积分之和）。

**(3)** **不能**。Green公式要求积分曲线是**封闭**的简单闭曲线。这里 $L$ 是开放的直线段，没有围成区域，因此无法使用Green公式。应使用文件02（第十四章）中第二型曲线积分的参数化方法直接计算。

</details>

---

**2.**（直接应用Green公式）利用Green公式计算下列闭曲线积分：

**(1)** $\displaystyle \oint_C (x^2 - y)\,dx + (y^2 + x)\,dy$，$C$ 是以 $(0,0), (1,0), (1,1), (0,1)$ 为顶点的正方形正向边界

**(2)** $\displaystyle \oint_C y\,dx + x\,dy$，$C$ 是圆周 $x^2 + y^2 = 1$（正向）

<details><summary>参考答案</summary>

**(1)** $P = x^2 - y$，$Q = y^2 + x$。

$$\frac{\partial Q}{\partial x} = 1,\quad \frac{\partial P}{\partial y} = -1$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 1 - (-1) = 2$$

$D$ 是正方形 $0 \leq x \leq 1,\; 0 \leq y \leq 1$。

由Green公式：
$$\begin{aligned}
\oint_C (x^2 - y)\,dx + (y^2 + x)\,dy &= \iint_D 2\,dA \\
&= 2 \times \text{Area}(D) = 2 \times (1 \times 1) = 2
\end{aligned}$$

因此积分值为 $2$。

**(2)** $P = y$，$Q = x$。

$$\frac{\partial Q}{\partial x} = 1,\quad \frac{\partial P}{\partial y} = 1$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 1 - 1 = 0$$

由Green公式：
$$\oint_C y\,dx + x\,dy = \iint_D 0\,dA = 0$$

因此对任何闭曲线 $C$（只要 $C$ 围成的区域内 $P, Q$ 连续可微），该积分恒为零。$\checkmark$

注意：由Green公式，$\partial Q/\partial x = \partial P/\partial y$ 等价于该场是保守场——线积分与路径无关，闭曲线积分为零。

</details>

---

**3.**（面积公式应用）利用面积公式 $A = \dfrac{1}{2}\oint_C x\,dy - y\,dx$ 计算下列图形的面积：

**(1)** 圆心在原点、半径为 $r$ 的圆

**(2)** 椭圆 $\dfrac{x^2}{a^2} + \dfrac{y^2}{b^2} = 1$

**(3)** 顶点为 $(0,0), (a,0), (a,b), (0,b)$ 的矩形（$a,b > 0$），取正向边界

<details><summary>参考答案</summary>

**(1)** 圆的参数化：$x = r\cos t$，$y = r\sin t$，$0 \leq t \leq 2\pi$（正向：逆时针）。

$dx = -r\sin t\,dt$，$dy = r\cos t\,dt$。

$$\begin{aligned}
A &= \frac{1}{2} \oint_C x\,dy - y\,dx \\
  &= \frac{1}{2} \int_0^{2\pi} \big[ (r\cos t)(r\cos t) - (r\sin t)(-r\sin t) \big]\,dt \\
  &= \frac{1}{2} \int_0^{2\pi} (r^2\cos^2 t + r^2\sin^2 t)\,dt \\
  &= \frac{1}{2} \int_0^{2\pi} r^2\,dt \\
  &= \frac{1}{2} \cdot r^2 \cdot 2\pi = \pi r^2 \quad \checkmark
\end{aligned}$$

**(2)** 椭圆的参数化：$x = a\cos t$，$y = b\sin t$，$0 \leq t \leq 2\pi$（正向：逆时针）。

$dx = -a\sin t\,dt$，$dy = b\cos t\,dt$。

$$\begin{aligned}
A &= \frac{1}{2} \oint_C x\,dy - y\,dx \\
  &= \frac{1}{2} \int_0^{2\pi} \big[ (a\cos t)(b\cos t) - (b\sin t)(-a\sin t) \big]\,dt \\
  &= \frac{1}{2} \int_0^{2\pi} (ab\cos^2 t + ab\sin^2 t)\,dt \\
  &= \frac{1}{2} \int_0^{2\pi} ab\,dt \\
  &= \frac{1}{2} \cdot ab \cdot 2\pi = \pi ab \quad \checkmark
\end{aligned}$$

**(3)** 矩形边界由四条直线段组成，正向（逆时针）为：$(0,0) \to (a,0) \to (a,b) \to (0,b) \to (0,0)$。

将矩形边界 $C = C_1 \cup C_2 \cup C_3 \cup C_4$ 分段处理：

- $C_1$：$y = 0$，$x$ 从 $0$ 到 $a$，$dy = 0$
- $C_2$：$x = a$，$y$ 从 $0$ 到 $b$，$dx = 0$
- $C_3$：$y = b$，$x$ 从 $a$ 到 $0$，$dy = 0$
- $C_4$：$x = 0$，$y$ 从 $b$ 到 $0$，$dx = 0$

用面积公式 $A = \frac12 \oint_C x\,dy - y\,dx$ 计算：

**$C_1$ 上**（$y=0$，$dy=0$，$x$ 从 $0$ 到 $a$）：
$$\int_{C_1} x\,dy - y\,dx = \int_{C_1} x\cdot 0 - 0\cdot dx = 0$$

**$C_2$ 上**（$x=a$，$dx=0$，$y$ 从 $0$ 到 $b$）：
$$\int_{C_2} x\,dy - y\,dx = \int_0^b a\,dy - y\cdot 0 = a\cdot b$$

**$C_3$ 上**（$y=b$，$dy=0$，$x$ 从 $a$ 到 $0$）：
$$\int_{C_3} x\,dy - y\,dx = \int_C^0 x\cdot 0 - b\cdot dx = -\int_0^a (-b)\,dx = \int_a^0 (-b)\,dx = b\cdot a$$

等等，让我重新计算 $C_3$：

在 $C_3$ 上：$y = b$，$dy = 0$，$x$ 从 $a$ 到 $0$（注意方向），$dx = dx$。

$$\int_{C_3} x\,dy - y\,dx = \int_{x=a}^{x=0} \big[ x\cdot 0 - b\,dx \big] = \int_a^0 (-b)\,dx = -b \cdot (0 - a) = ab$$

或者参数化：$x = t$，$y = b$，$t$ 从 $a$ 到 $0$，$dx = dt$，$dy = 0$：
$$\int_{C_3} = \int_a^0 \big[ t\cdot 0 - b\cdot 1 \big]\,dt = \int_a^0 (-b)\,dt = -b(0-a) = ab$$

**$C_4$ 上**（$x=0$，$dx=0$，$y$ 从 $b$ 到 $0$）：
$$\int_{C_4} x\,dy - y\,dx = \int_b^0 0\cdot dy - y\cdot 0 = 0$$

**求和**：
$$\oint_C x\,dy - y\,dx = 0 + ab + ab + 0 = 2ab$$

因此：
$$A = \frac{1}{2} \cdot 2ab = ab \quad \checkmark$$

即矩形面积为长 $a$ 宽 $b$ 的乘积。验证了面积公式对多边形区域也适用。

**拓展**：对于任意多边形，将面积公式离散化就得到"鞋带公式"（shoelace formula），这是大地测量中计算多边形面积的标准方法。本题中矩形是一个特例。

</details>

---

### 迁移应用

**4.**（Green公式与方向反转）设 $C$ 为正向圆周 $x^2 + y^2 = 4$，计算曲线积分：

$$I = \oint_C (x^2 + y)\,dx + (y^2 - x)\,dy$$

**(1)** 利用Green公式直接计算。

**(2)** 若将 $C$ 改为顺时针方向，结果如何？

<details><summary>参考答案</summary>

**(1)** $P = x^2 + y$，$Q = y^2 - x$。

$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}(y^2 - x) = -1$$
$$\frac{\partial P}{\partial y} = \frac{\partial}{\partial y}(x^2 + y) = 1$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = -1 - 1 = -2$$

由Green公式：
$$I = \iint_D (-2)\,dA = -2 \times \text{Area}(D)$$

$D$ 是半径为 $2$ 的圆盘，面积 $= \pi \cdot 2^2 = 4\pi$。

因此：
$$I = -2 \times 4\pi = -8\pi$$

**(2)** 若 $C$ 改为顺时针方向（$C^-$），由方向反转性质：
$$I_{C^-} = -I = -(-8\pi) = 8\pi$$

**验证**：用修正后的Green公式：
$$\oint_{C^-} P\,dx + Q\,dy = -\iint_D (\partial Q/\partial x - \partial P/\partial y)\,dA = -\iint_D (-2)\,dA = 2 \times 4\pi = 8\pi$$

结果一致。$\checkmark$

</details>

---

**5.**（Green公式的路径无关性应用）已知向量场 $\mathbf{F}(x,y) = (y, x)$ 满足 $\partial Q/\partial x = \partial P/\partial y$（即旋度为零）。

**(1)** 利用Green公式证明：对 $\mathbb{R}^2$ 中的任何简单闭曲线 $C$，有 $\displaystyle \oint_C y\,dx + x\,dy = 0$。

**(2)** 由此证明：$\displaystyle \int_{L} y\,dx + x\,dy$ 与路径无关（即从点 $A$ 到点 $B$ 沿任意分段光滑曲线的积分值都相等）。

<details><summary>参考答案</summary>

**(1)** 对于 $\mathbf{F} = (y, x)$：
$$P = y,\quad Q = x$$
$$\frac{\partial Q}{\partial x} = 1,\quad \frac{\partial P}{\partial y} = 1$$
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = 1 - 1 = 0$$

设 $C$ 为 $\mathbb{R}^2$ 中任意一条正向简单闭曲线，$D$ 为 $C$ 所围的区域。由Green公式：
$$\oint_C y\,dx + x\,dy = \iint_D 0\,dA = 0$$

由于 $\partial Q/\partial x - \partial P/\partial y \equiv 0$ 在全平面成立，这一结论对 $\mathbb{R}^2$ 中任何简单闭曲线 $C$ 都成立。

**(2)** 证明路径无关：设 $\Gamma_1$ 和 $\Gamma_2$ 是 $\mathbb{R}^2$ 中从 $A$ 到 $B$ 的两条不同路径。考虑封闭曲线 $C = \Gamma_1 \cup \Gamma_2^-$（即从 $A$ 沿 $\Gamma_1$ 到 $B$，再沿 $\Gamma_2$ 的反向从 $B$ 回到 $A$）。

由 (1)：
$$\oint_C y\,dx + x\,dy = 0$$

由曲线可加性：
$$\int_{\Gamma_1} y\,dx + x\,dy + \int_{\Gamma_2^-} y\,dx + x\,dy = 0$$

由方向反转性质：
$$\int_{\Gamma_1} y\,dx + x\,dy - \int_{\Gamma_2} y\,dx + x\,dy = 0$$

因此：
$$\int_{\Gamma_1} y\,dx + x\,dy = \int_{\Gamma_2} y\,dx + x\,dy$$

即积分值与路径无关。$\checkmark$

**联系**：这个结论与文件02（第十四章）的迁移应用练习题5一致——那里我们通过具体计算验证了沿两条不同路径的积分值相等（均等于1）。现在Green公式从理论上给出了路径无关性的条件：$\partial Q/\partial x = \partial P/\partial y$。

</details>

---

**6.**（思考题——复连通区域上的Green公式）设 $D$ 是圆环区域 $1 \leq x^2 + y^2 \leq 4$，其边界由外圆周 $C_1: x^2 + y^2 = 4$（正向）和内圆周 $C_2: x^2 + y^2 = 1$（正向）组成。计算：

$$\oint_{C_1} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy + \oint_{C_2} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy$$

其中 $C_1$ 取逆时针方向，$C_2$ 取顺时针方向（注意：对圆环区域，$C_2$ 的正向是顺时针）。

<details><summary>参考答案</summary>

这是复连通区域上的Green公式问题。在圆环区域 $D$ 上，$P = \dfrac{-y}{x^2+y^2}$，$Q = \dfrac{x}{x^2+y^2}$。计算旋度：

$$\frac{\partial Q}{\partial x} = \frac{\partial}{\partial x}\left( \frac{x}{x^2+y^2} \right) = \frac{(x^2+y^2) - x(2x)}{(x^2+y^2)^2} = \frac{y^2 - x^2}{(x^2+y^2)^2}$$

$$\frac{\partial P}{\partial y} = \frac{\partial}{\partial y}\left( \frac{-y}{x^2+y^2} \right) = -\frac{(x^2+y^2) - y(2y)}{(x^2+y^2)^2} = -\frac{x^2 - y^2}{(x^2+y^2)^2} = \frac{y^2 - x^2}{(x^2+y^2)^2}$$

因此：
$$\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = \frac{y^2 - x^2}{(x^2+y^2)^2} - \frac{y^2 - x^2}{(x^2+y^2)^2} = 0$$

但注意：$P, Q$ 在原点 $(0,0)$ 无定义，而原点不在 $D$ 内（因为 $D$ 是 $1 \leq x^2+y^2 \leq 4$ 的圆环，原点被挖掉了）。在 $D$ 内，$P, Q$ 有连续的偏导数。

由复连通区域上的Green公式：
$$\oint_{\partial D} P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dA = \iint_D 0\,dA = 0$$

而 $\partial D = C_1 \cup C_2$，其中 $C_1$ 取逆时针，$C_2$ 取顺时针（正向——区域 $D$ 在左侧）。

因此：
$$\oint_{C_1} P\,dx + Q\,dy + \oint_{C_2} P\,dx + Q\,dy = 0$$

即：
$$\oint_{C_1} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy = -\oint_{C_2} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy$$

这个结果很有意思：如果单独看 $C_1$ 上的积分，$P, Q$ 在 $C_1$ 围成的圆盘内有奇点（原点），不能直接应用Green公式。但是通过挖去内圆盘 $C_2$，我们可以在圆环区域 $D$ 上安全地使用Green公式，得到两个边界积分之和为零。

**实际计算**：通过直接参数化可以验证（读者可自行计算）：
$$\oint_{C_1} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy = 2\pi$$
$$\oint_{C_2} \frac{-y}{x^2+y^2}\,dx + \frac{x}{x^2+y^2}\,dy = -2\pi$$

验证了和为 $0$。注意 $C_2$ 上的积分值为 $-2\pi$ 是因为 $C_2$ 取顺时针方向，而 $2\pi$ 是逆时针圆周上的值。

</details>

