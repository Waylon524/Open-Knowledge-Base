# 08. 多元函数的Taylor公式

> 所属章节：第十二章 多元函数的微分学  |  文件序号：08  |  难度：综合
> 常见混淆点：1) 带Lagrange余项与带Peano余项的Taylor公式的适用条件不同——前者需要整体区间上的高阶可微性而后者只需一点处的；2) 乘法展开时需根据"总次数"筛选保留的项，不能简单保留所有项——$x^3y$ 的总次数为4，若展开到3阶则应舍弃

## 1. 学习目标与先修前置

### 学习目标
- 掌握辅助函数 $\varphi(t)=f(\mathbf{P}_0+t\mathbf{h})$ 构造法，理解多元Taylor公式向一元的化归推导
- 掌握带Lagrange余项的一阶Taylor公式（多元微分中值定理），并能用它证明 $\nabla f\equiv\mathbf{0}\Rightarrow f\equiv\text{常数}$
- 掌握带Lagrange余项的二阶Taylor公式（偏导数在中间点取值的形式）
- 掌握多重指标记号 $\alpha=(\alpha_1,\alpha_2)$、$|\alpha|$、$\alpha!$、$\partial^\alpha f$ 及n阶Taylor展开的通项公式
- 掌握算子形式 $(\Delta x\partial_x+\Delta y\partial_y)^kf$ 及其与多重指标的关系
- 掌握两个一元展开式相乘得到二元展开的乘法技术及总次数筛选规则
- 掌握二阶Taylor近似计算的标准化四步工作流程
- 掌握利用正定Hessian矩阵的 $\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}\geq m\|\mathbf{h}\|^2$ 性质和 $o(\|\mathbf{h}\|^2)$ 处理技巧证明极值充分条件

### 先修知识
- 文件 01（第十二章）：全微分 $df=f_xdx+f_ydy$ 与一阶线性近似 $f\approx f_0+f_x\Delta x+f_y\Delta y$
- 文件 03（第十二章）：全微分的定义（定义12.3），可微性的判别
- 文件 04（第十二章）：多元复合函数的链式法则（定理12.7），抽象复合函数的偏导记号
- 文件 06（第十二章）：Hessian矩阵的定义（定义12.8），带Peano余项的二阶Taylor公式（证明思路中已给出向量形式），Hessian判别法（定理12.12）
- 文件 03（第五章）：一元函数的泰勒公式（定理5.5带Peano余项，定理5.6带Lagrange余项），常见麦克劳林展开式
- 文件 01（第五章）：拉格朗日中值定理（定理5.3）：$f(b)-f(a)=f'(\xi)(b-a)$

---

## 2. 背景与应用场景

在文件06中，我们为了证明Hessian判别法（定理12.12），已经使用了如下二阶Taylor展开：
$$f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)=\nabla f\cdot\mathbf{h}+\frac12\mathbf{h}^\mathsf{T}\mathbf{H}f\,\mathbf{h}+o(\|\mathbf{h}\|^2)$$

但那个公式是作为"已知结论"直接引用的。本章的目标是：**系统建立多元函数的Taylor公式理论**，包括：

- **一般n阶公式**：从一元到多元的完整推广
- **两种余项形式**：Peano余项（局部估计）和Lagrange余项（区间估计）
- **辅助函数构造法**：将多元问题化归为一元问题
- **算子逼近技术**：$(\Delta x\partial_x+\Delta y\partial_y)^k$ 的高效计算方法

**实际应用**：
- **数值近似**：用二阶Taylor公式比一阶线性近似获得更高的计算精度（如 $8.96^{2.03}$ 的近似）
- **极值判定**：利用二阶Taylor展开证明极值的充分条件——这是Hessian判别法的理论基础
- **误差分析**：Lagrange余项给出近似计算中的误差上界
- **函数逼近**：将复杂函数用多项式近似，是数值分析和科学计算的基础

**知识链位置**：本章前7个文件依次建立了偏导数、方向导数、全微分、链式法则、隐函数、无条件极值和条件极值的理论体系。多元Taylor公式是这些工具的**顶峰**——它将一阶近似（全微分）推广到任意高阶，并为极值理论提供了严格的理论基础。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $\mathbf{P}_0=(x_0,y_0)$ | 展开中心点 | 向量表示法 |
| $\mathbf{h}=(\Delta x,\Delta y)$ | 自变量增量向量 | $\|\mathbf{h}\|=\sqrt{\Delta x^2+\Delta y^2}$ |
| $\varphi(t)=f(\mathbf{P}_0+t\mathbf{h})$ | 辅助函数 | 将多元问题化归为一元 |
| $\alpha=(\alpha_1,\alpha_2)$ | 多重指标 | $\alpha_1,\alpha_2$ 为非负整数 |
| $|\alpha|=\alpha_1+\alpha_2$ | 多重指标的阶 | 即总求导次数 |
| $\alpha!=\alpha_1!\alpha_2!$ | 多重指标的阶乘 | 各分量阶乘之积 |
| $\partial^\alpha f=\partial_x^{\alpha_1}\partial_y^{\alpha_2}f$ | 多重指标偏导 | 对 $x$ 求 $\alpha_1$ 次、对 $y$ 求 $\alpha_2$ 次 |
| $\mathbf{h}^\alpha=(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}$ | 增量的幂 | 各分量增量幂之积 |
| $(\Delta x\partial_x+\Delta y\partial_y)^kf$ | 算子形式的 $k$ 阶偏导组合 | 展开后即用二项式定理 |
| $\rho=\|\mathbf{h}\|$ | 增量向量的模长 | 用于Peano余项 $o(\rho^n)$ |
| $m=\lambda_{\min}(\mathbf{H})$ | Hessian矩阵最小特征值 | 用于正定矩阵的下界估计 |

### 3.1 一元Taylor公式回顾

为方便类比,先回顾一元函数的Taylor公式（来自第五章文件03）：

**带Peano余项（定理5.5）**：若 $f$ 在 $a$ 处 $n$ 阶可导，则当 $x\to a$ 时：
$$f(x)=\sum_{k=0}^n\frac{f^{(k)}(a)}{k!}(x-a)^k+o\!\big((x-a)^n\big)$$

**带Lagrange余项（定理5.6）**：若 $f$ 在 $[a,x]$ 上有连续 $n$ 阶导数且在 $(a,x)$ 内有 $n+1$ 阶导数，则存在 $\xi$ 介于 $a$ 与 $x$ 之间：
$$f(x)=\sum_{k=0}^n\frac{f^{(k)}(a)}{k!}(x-a)^k+\frac{f^{(n+1)}(\xi)}{(n+1)!}(x-a)^{n+1}$$

**核心理念**：泰勒多项式通过匹配点 $a$ 处的各阶导数值来逼近原函数，余项描述了近似的误差。

### 3.2 多重指标记号（A2核心）

为了简洁地表示多元函数的n阶Taylor公式，我们需要引入多重指标（multi-index）记号。

**定义（多重指标）**：设 $\alpha=(\alpha_1,\alpha_2)$，其中 $\alpha_1,\alpha_2$ 均为非负整数。定义：
- **阶**：$|\alpha|=\alpha_1+\alpha_2$
- **阶乘**：$\alpha!=\alpha_1!\alpha_2!$
- **偏导算子**：$\partial^\alpha f=\displaystyle\frac{\partial^{|\alpha|}f}{\partial x^{\alpha_1}\partial y^{\alpha_2}}=\partial_x^{\alpha_1}\partial_y^{\alpha_2}f$
- **幂**：$\mathbf{h}^\alpha=(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}$

**示例**：
| $\alpha$ | $|\alpha|$ | $\alpha!$ | $\partial^\alpha f$ | $\mathbf{h}^\alpha$ |
|----------|:---------:|:---------:|:-----------------:|:-----------------:|
| $(0,0)$ | $0$ | $0!0!=1$ | $f$ | $1$ |
| $(1,0)$ | $1$ | $1!0!=1$ | $f_x$ | $\Delta x$ |
| $(0,1)$ | $1$ | $0!1!=1$ | $f_y$ | $\Delta y$ |
| $(2,0)$ | $2$ | $2!0!=2$ | $f_{xx}$ | $\Delta x^2$ |
| $(1,1)$ | $2$ | $1!1!=1$ | $f_{xy}$ | $\Delta x\Delta y$ |
| $(0,2)$ | $2$ | $0!2!=2$ | $f_{yy}$ | $\Delta y^2$ |
| $(2,1)$ | $3$ | $2!1!=2$ | $f_{xxy}$ | $\Delta x^2\Delta y$ |

---

## 4. 原理与方法

### 4.1 核心方法：辅助函数构造法

将多元Taylor公式化归为一元Taylor公式的核心工具是**辅助函数**。

**构造**：设 $f(x,y)$ 在点 $\mathbf{P}_0=(x_0,y_0)$ 的某邻域内有足够高阶的连续偏导。令 $\mathbf{h}=(\Delta x,\Delta y)$，定义关于 $t$ 的一元函数：
$$\varphi(t)=f(x_0+t\Delta x,\;y_0+t\Delta y),\quad t\in[0,1]$$

则 $\varphi(0)=f(x_0,y_0)$，$\varphi(1)=f(x_0+\Delta x,y_0+\Delta y)$。我们的目标是将 $\varphi(1)$ 在 $t=0$ 处展开。

**链式法则求 $\varphi(t)$ 的各阶导数**：

一阶导数：
$$\varphi'(t)=\frac{\partial f}{\partial x}\cdot\Delta x+\frac{\partial f}{\partial y}\cdot\Delta y=f_x\Delta x+f_y\Delta y$$

二阶导数（对 $\varphi'(t)$ 再次使用链式法则）：
$$
\begin{aligned}
\varphi''(t)&=\frac{\partial}{\partial x}(f_x\Delta x+f_y\Delta y)\cdot\Delta x+\frac{\partial}{\partial y}(f_x\Delta x+f_y\Delta y)\cdot\Delta y \\
&=(f_{xx}\Delta x+f_{xy}\Delta y)\Delta x+(f_{yx}\Delta x+f_{yy}\Delta y)\Delta y \\
&=f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2
\end{aligned}
$$

这里利用了 $f_{xy}=f_{yx}$（$f$ 二阶连续可微时成立）。

三阶导数：
$$
\varphi'''(t)=f_{xxx}\Delta x^3+3f_{xxy}\Delta x^2\Delta y+3f_{xyy}\Delta x\Delta y^2+f_{yyy}\Delta y^3
$$

**算子记法**：上述模式可统一写作算子形式：
$$\varphi^{(k)}(t)=\big(\Delta x\partial_x+\Delta y\partial_y\big)^k f(x_0+t\Delta x,\;y_0+t\Delta y)$$

其中 $(\Delta x\partial_x+\Delta y\partial_y)^k$ 是一个微分算子，按二项式定理展开：
$$(\Delta x\partial_x+\Delta y\partial_y)^k=\sum_{j=0}^k\binom{k}{j}(\Delta x)^{k-j}(\Delta y)^j\partial_x^{k-j}\partial_y^j$$

**验证**：$k=1$ 时得 $\Delta x\partial_x+\Delta y\partial_y$，$k=2$ 时得 $(\Delta x)^2\partial_x^2+2\Delta x\Delta y\partial_x\partial_y+(\Delta y)^2\partial_y^2$，与上述手动计算一致。

---

### 4.2 带Lagrange余项的一阶Taylor公式（多元微分中值定理）（A1）

**定理 12.13（多元微分中值定理）**：设 $f(x,y)$ 在凸区域 $D\subset\mathbb{R}^2$ 上具有一阶连续偏导数（$f\in C^1$）。则对 $D$ 中任意两点 $\mathbf{P}_0=(x_0,y_0)$ 和 $\mathbf{P}_0+\mathbf{h}=(x_0+\Delta x,y_0+\Delta y)$，存在 $\theta\in(0,1)$，使得
$$\boxed{f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)=f_x(\mathbf{P}_0+\theta\mathbf{h})\Delta x+f_y(\mathbf{P}_0+\theta\mathbf{h})\Delta y}$$

其中 $\mathbf{P}_0+\theta\mathbf{h}$ 表示连接 $\mathbf{P}_0$ 和 $\mathbf{P}_0+\mathbf{h}$ 的线段上的某一点。

**证明**：

**第 1 步：构造辅助函数**。令 $\varphi(t)=f(x_0+t\Delta x,\;y_0+t\Delta y)$，$t\in[0,1]$。

**第 2 步：求 $\varphi$ 的一阶导数**。由链式法则：
$$\varphi'(t)=f_x(x_0+t\Delta x,\;y_0+t\Delta y)\Delta x+f_y(x_0+t\Delta x,\;y_0+t\Delta y)\Delta y$$

**第 3 步：对一元函数 $\varphi$ 应用拉格朗日中值定理**。$\varphi$ 在 $[0,1]$ 上连续且在 $(0,1)$ 内可导（因为 $f\in C^1$），由一元拉格朗日中值定理（定理5.3），存在 $\theta\in(0,1)$ 使得：
$$\varphi(1)-\varphi(0)=\varphi'(\theta)$$

**第 4 步：代回原变量**。$\varphi(1)-\varphi(0)=f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)$，而：
$$\varphi'(\theta)=f_x(\mathbf{P}_0+\theta\mathbf{h})\Delta x+f_y(\mathbf{P}_0+\theta\mathbf{h})\Delta y$$

因此原式成立。证毕。

**与一元拉格朗日中值定理的类比**：

| 对比项 | 一元 | 多元 |
|--------|:----:|:----:|
| 公式 | $f(b)-f(a)=f'(\xi)(b-a)$ | $f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)=\nabla f(\mathbf{P}_0+\theta\mathbf{h})\cdot\mathbf{h}$ |
| 中间点 | $\xi$ 在 $a$ 与 $b$ 之间 | $\mathbf{P}_0+\theta\mathbf{h}$ 在线段上 |
| 应用 | 估计函数增量 | 估计函数的全增量 |

**推论（梯度为零则函数为常数）**：若 $f$ 在凸区域 $D$ 上具有一阶连续偏导数，且 $\nabla f\equiv\mathbf{0}$ 在 $D$ 上处处成立，则 $f$ 在 $D$ 上恒为常数。

**证明**：

任取 $D$ 中两点 $\mathbf{A},\mathbf{B}$。由于 $D$ 是凸区域，连接 $\mathbf{A}$ 和 $\mathbf{B}$ 的线段 $\mathbf{A}+t(\mathbf{B}-\mathbf{A})$（$t\in[0,1]$）完全包含在 $D$ 内。

由定理12.13（微分中值定理），存在 $\theta\in(0,1)$ 使得：
$$f(\mathbf{B})-f(\mathbf{A})=\nabla f(\mathbf{A}+\theta(\mathbf{B}-\mathbf{A}))\cdot(\mathbf{B}-\mathbf{A})$$

而题设 $\nabla f\equiv\mathbf{0}$，故右端 $=\mathbf{0}\cdot(\mathbf{B}-\mathbf{A})=0$。因此 $f(\mathbf{B})=f(\mathbf{A})$。

由 $\mathbf{A},\mathbf{B}$ 的任意性，$f$ 在 $D$ 上为常数。证毕。

**注**：凸区域条件不可缺少。若 $D$ 不是凸的（例如环形区域），可能存在两点使得连接线段不完全在 $D$ 内，此时定理12.13不能直接应用。

---

### 4.3 带Lagrange余项的二阶Taylor公式（A1）

**定理 12.14（带Lagrange余项的二阶Taylor公式）**：设 $f(x,y)$ 在点 $\mathbf{P}_0=(x_0,y_0)$ 的某邻域内具有二阶连续偏导数（$f\in C^2$）。则对充分小的 $\mathbf{h}=(\Delta x,\Delta y)$，存在 $\theta\in(0,1)$，使得
$$\boxed{f(\mathbf{P}_0+\mathbf{h})=f(\mathbf{P}_0)+\nabla f(\mathbf{P}_0)\cdot\mathbf{h}+\frac12\big(f_{xx}(\xi)\Delta x^2+2f_{xy}(\xi)\Delta x\Delta y+f_{yy}(\xi)\Delta y^2\big)}$$

其中 $\xi=\mathbf{P}_0+\theta\mathbf{h}=(x_0+\theta\Delta x,\;y_0+\theta\Delta y)$ 是线段上的中间点。

**证明**（使用辅助函数法）：

**第 1 步**：令 $\varphi(t)=f(x_0+t\Delta x,\;y_0+t\Delta y)$，$t\in[0,1]$。则 $\varphi$ 是 $t$ 的一元函数且二阶可导。

**第 2 步**：计算 $\varphi'(t)$ 和 $\varphi''(t)$（由4.1节）：
$$\varphi'(t)=f_x\Delta x+f_y\Delta y$$
$$\varphi''(t)=f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2$$

其中各偏导均在 $(x_0+t\Delta x,\;y_0+t\Delta y)$ 处取值。

**第 3 步**：对 $\varphi$ 应用一元带Lagrange余项的一阶Taylor公式（即定理5.6中 $n=1$ 的情形）：
$$\varphi(1)=\varphi(0)+\varphi'(0)\cdot1+\frac12\varphi''(\theta)\cdot1^2,\quad\theta\in(0,1)$$

**第 4 步**：代回原变量。$\varphi(1)=f(\mathbf{P}_0+\mathbf{h})$，$\varphi(0)=f(\mathbf{P}_0)$，$\varphi'(0)=\nabla f(\mathbf{P}_0)\cdot\mathbf{h}$，而 $\varphi''(\theta)=f_{xx}(\xi)\Delta x^2+2f_{xy}(\xi)\Delta x\Delta y+f_{yy}(\xi)\Delta y^2$（其中 $\xi=\mathbf{P}_0+\theta\mathbf{h}$）。代入即得定理。证毕。

**Peano余项与Lagrange余项的对比**：

| 对比项 | Peano余项 | Lagrange余项 |
|--------|:---------:|:-----------:|
| 公式 | $o(\|\mathbf{h}\|^2)$ | $\frac12(f_{xx}(\xi)\Delta x^2+2f_{xy}(\xi)\Delta x\Delta y+f_{yy}(\xi)\Delta y^2)$ |
| 条件 | $f$ 在 $\mathbf{P}_0$ 处二阶可导 | $f$ 在邻域内二阶连续可导 |
| 信息量 | 只知趋于零的速度 | 给出余项具体表达式 |
| 适用场景 | 极限计算、局部分析 | 误差估计、区间上的不等式证明 |
| 文件来源 | 文件06已作为工具使用 | 本文件首次系统建立 |

---

### 4.4 一般n阶Taylor公式与多重指标记号（A2）

利用辅助函数 $\varphi(t)=f(\mathbf{P}_0+t\mathbf{h})$ 和一元Taylor公式，我们可以将多元Taylor公式推广到任意 $n$ 阶。

**定理 12.15（带Peano余项的n阶Taylor公式）**：设 $f(x,y)$ 在点 $\mathbf{P}_0=(x_0,y_0)$ 处 $n$ 阶可微。则当 $\mathbf{h}\to\mathbf{0}$ 时，有
$$\boxed{f(\mathbf{P}_0+\mathbf{h})=\sum_{k=0}^n\frac{1}{k!}\big(\Delta x\partial_x+\Delta y\partial_y\big)^k f(\mathbf{P}_0)+o(\|\mathbf{h}\|^n)}$$

**定理 12.16（带Lagrange余项的n阶Taylor公式）**：设 $f(x,y)$ 在 $\mathbf{P}_0$ 的某邻域内具有 $n+1$ 阶连续偏导数。则存在 $\theta\in(0,1)$，使得
$$\boxed{f(\mathbf{P}_0+\mathbf{h})=\sum_{k=0}^n\frac{1}{k!}\big(\Delta x\partial_x+\Delta y\partial_y\big)^k f(\mathbf{P}_0)+R_n}$$

其中Lagrange余项为：
$$R_n=\frac{1}{(n+1)!}\big(\Delta x\partial_x+\Delta y\partial_y\big)^{n+1}f(\mathbf{P}_0+\theta\mathbf{h})$$

**证明概要**：
对 $\varphi(t)=f(\mathbf{P}_0+t\mathbf{h})$ 应用一元Taylor公式。由4.1节知 $\varphi^{(k)}(0)=(\Delta x\partial_x+\Delta y\partial_y)^k f(\mathbf{P}_0)$。代入一元带Peano或Lagrange余项的Taylor公式即得。证毕。

#### 4.4.1 用多重指标表示的通项公式

算子形式虽然简洁，但在具体展开时仍需要计算二项式系数。使用多重指标记号可以给出更对称的形式。

将二项式定理应用于算子：
$$(\Delta x\partial_x+\Delta y\partial_y)^k=\sum_{\alpha_1+\alpha_2=k}\frac{k!}{\alpha_1!\alpha_2!}(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}\partial_x^{\alpha_1}\partial_y^{\alpha_2}$$

用多重指标 $\alpha=(\alpha_1,\alpha_2)$ 重写为：
$$(\Delta x\partial_x+\Delta y\partial_y)^k=\sum_{|\alpha|=k}\frac{k!}{\alpha!}\,\mathbf{h}^\alpha\,\partial^\alpha$$

因此Taylor公式的通项为：
$$\frac{1}{k!}(\Delta x\partial_x+\Delta y\partial_y)^k f(\mathbf{P}_0)=\frac{1}{k!}\sum_{|\alpha|=k}\frac{k!}{\alpha!}\,\partial^\alpha f(\mathbf{P}_0)\,\mathbf{h}^\alpha=\sum_{|\alpha|=k}\frac{1}{\alpha!}\,\partial^\alpha f(\mathbf{P}_0)\,\mathbf{h}^\alpha$$

**用多重指标表示的n阶Taylor公式**（最简洁的形式）：

**Peano余项**：
$$\boxed{f(\mathbf{P}_0+\mathbf{h})=\sum_{|\alpha|\leq n}\frac{1}{\alpha!}\,\partial^\alpha f(\mathbf{P}_0)\,\mathbf{h}^\alpha+o(\|\mathbf{h}\|^n)}$$

**Lagrange余项**：
$$f(\mathbf{P}_0+\mathbf{h})=\sum_{|\alpha|\leq n}\frac{1}{\alpha!}\,\partial^\alpha f(\mathbf{P}_0)\,\mathbf{h}^\alpha+\sum_{|\alpha|=n+1}\frac{1}{\alpha!}\,\partial^\alpha f(\mathbf{P}_0+\theta\mathbf{h})\,\mathbf{h}^\alpha$$

**示例**：$n=2$ 时展开（展开到 $|\alpha|\leq2$ 的所有项）：
$$
\begin{aligned}
f(\mathbf{P}_0+\mathbf{h})&=f(\mathbf{P}_0) \\
&\quad+\big(f_x\Delta x+f_y\Delta y\big) \\
&\quad+\frac12\big(f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2\big) \\
&\quad+o(\|\mathbf{h}\|^2)
\end{aligned}
$$

这正是我们在文件06中使用的公式，现在可以看到它是一般n阶公式的 $n=2$ 特例。

#### 4.4.2 多重指标求和表（$n=3$ 时所有项）

| $|\alpha|$ | 所有 $\alpha$ | 项数 | 通项 $\frac1{\alpha!}\partial^\alpha f\,\mathbf{h}^\alpha$ |
|:---------:|:-------------:|:----:|:-------------------------------------------------------:|
| 0 | $(0,0)$ | 1 | $f$ |
| 1 | $(1,0),(0,1)$ | 2 | $f_x\Delta x+f_y\Delta y$ |
| 2 | $(2,0),(1,1),(0,2)$ | 3 | $\frac12f_{xx}\Delta x^2+f_{xy}\Delta x\Delta y+\frac12f_{yy}\Delta y^2$ |
| 3 | $(3,0),(2,1),(1,2),(0,3)$ | 4 | $\frac16f_{xxx}\Delta x^3+\frac12f_{xxy}\Delta x^2\Delta y+\frac12f_{xyy}\Delta x\Delta y^2+\frac16f_{yyy}\Delta y^3$ |

**项数规律**：$|\alpha|=k$ 时的项数为 $k+1$（因为 $\alpha_1$ 可取 $0,1,\ldots,k$，$\alpha_2=k-\alpha_1$）。

---

### 4.5 乘法展开技术（A3）

对于形如 $f(x,y)=g(x)h(y)$ 的乘积函数（其中 $g$ 只依赖 $x$，$h$ 只依赖 $y$），可以利用一元展开式相乘得到二元展开式。

**核心原理**：若
$$g(x)=\sum_{i=0}^m a_i x^i+o(x^m),\quad h(y)=\sum_{j=0}^n b_j y^j+o(y^n)$$

则乘积为：
$$g(x)h(y)=\sum_{i=0}^m\sum_{j=0}^n a_i b_j x^i y^j+\text{高阶项}$$

关键在于：**每个单项 $x^iy^j$ 的总次数**（total degree）定义为 $i+j$。

**总次数筛选规则**：若要将 $f(x,y)$ 展开到 $n$ 阶（即所有总次数 $\leq n$ 的项），则只需保留：
- 所有满足 $i+j\leq n$ 的项 $a_i b_j x^i y^j$

对于 $i+j>n$ 的项，它们属于 $o(\rho^n)$（其中 $\rho=\sqrt{x^2+y^2}$），应合并到余项中。

**Peano余项的规范写法**：
- 展开到 $n$ 阶（即所有总次数不超过 $n$ 的项全部保留时），余项写作 $o(\rho^n)$
- 其中 $\rho=\sqrt{x^2+y^2}$ 是向量的模长
- 注意 $o(\rho^n)$ 与 $o((x^2+y^2)^{n/2})$ 是等价的

**无穷小量乘积的阶数判断**：
- 若 $A=o(x^p)$ 且 $B=o(y^q)$，则 $A\cdot B$ 作为 $(x,y)\to(0,0)$ 时的无穷小，其阶数至少为 $p+q$（但具体阶数还需看 $x$ 和 $y$ 的相对大小）
- 实际操作中：直接在展开式中丢弃所有总次数 $>n$ 的项，将它们的和统一纳入 $o(\rho^n)$

**重要安全准则**：
- 若 $g(x)$ 展开到 $x^m$ 且 $h(y)$ 展开到 $y^n$，则乘积展开中能保留的项的最大总次数为 $\min(m,n)$（考虑 $x$ 和 $y$ 对称的情形）或更一般地为 $m+n$（但实际操作时按总次数筛选）
- 当两个展开式中都有 $o$ 余项时，它们的乘积属于更高阶的无穷小，直接归入最终余项

---

### 4.6 二阶Taylor近似计算的标准化流程（A4）

利用二阶Taylor公式作近似计算时，可按以下四步法系统操作：

| 步骤 | 操作 | 说明 |
|------|------|------|
| 第1步：选函数与基点 | 构造 $f(x,y)$，选择展开点 $(x_0,y_0)$ 使各阶偏导易算且 $\Delta x,\Delta y$ 充分小 | 与一阶线性近似的第一步相同 |
| 第2步：计算偏导数值 | 计算 $f(x_0,y_0)$、$f_x$、$f_y$、$f_{xx}$、$f_{xy}$、$f_{yy}$ 在展开点的值 | 注意指数型函数要化成 $e^{b\ln a}$ 处理 |
| 第3步：写出二阶展开式 | $$f(x_0+\Delta x,y_0+\Delta y)\approx f(x_0,y_0)+f_x\Delta x+f_y\Delta y+\frac12(f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2)$$ | 使用Peano余项形式的二阶公式，忽略余项 |
| 第4步：代入计算 | 将 $\Delta x=x-x_0$，$\Delta y=y-y_0$ 代入展开式 | 得到数值近似结果 |

**与一阶线性近似的对比**：
- 一阶近似：$f\approx f_0+f_x\Delta x+f_y\Delta y$，误差 $O(\|\mathbf{h}\|^2)$
- 二阶近似：增加 $\frac12(f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2)$ 项，误差 $O(\|\mathbf{h}\|^3)$
- 因此对于足够小的 $\|\mathbf{h}\|$，二阶近似比一阶近似精度更高

**指数型函数 $a^b$ 的偏导计算方法**：对于 $f(x,y)=u(x,y)^{v(x,y)}$，取自然对数转化为指数形式：
$$f(x,y)=e^{v(x,y)\ln u(x,y)}$$
然后使用链式法则和乘积法则求偏导。特别地，当 $u$ 和 $v$ 均为简单形式时，可直接使用导数公式：
$$f_x=v\cdot u^{v-1}\cdot u_x+u^v\ln u\cdot v_x,\qquad f_y=v\cdot u^{v-1}\cdot u_y+u^v\ln u\cdot v_y$$

---

### 4.7 利用Taylor公式证明极值充分条件（A5）

在文件06中，我们使用配方法证明了Hessian判别法（定理12.12）。现在利用二阶Taylor公式和正定矩阵的性质，给出另一条证明路径。

**定理 12.17（极值的二阶充分条件——正定矩阵方法）**：设 $f(x,y)$ 在点 $\mathbf{P}_0=(x_0,y_0)$ 的某邻域内具有二阶连续偏导数，且：
1. $\nabla f(\mathbf{P}_0)=\mathbf{0}$（驻点条件）
2. Hessian矩阵 $\mathbf{H}= \begin{pmatrix} f_{xx}(\mathbf{P}_0) & f_{xy}(\mathbf{P}_0) \\ f_{xy}(\mathbf{P}_0) & f_{yy}(\mathbf{P}_0) \end{pmatrix}$ 正定

则 $f$ 在 $\mathbf{P}_0$ 处取得**严格局部极小值**。

**证明**：

**第 1 步：写出驻点处的二阶Taylor展开**。由定理12.15（带Peano余项的二阶Taylor公式），由于 $\nabla f(\mathbf{P}_0)=\mathbf{0}$，一阶项消失：
$$f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)=\frac12\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}+o(\|\mathbf{h}\|^2)$$

其中 $\mathbf{H}=\mathbf{H}f(\mathbf{P}_0)$ 是 $\mathbf{P}_0$ 处的Hessian矩阵。

**第 2 步：利用正定矩阵的下界性质**。$\mathbf{H}$ 正定 $\Rightarrow$ $\mathbf{H}$ 的所有特征值 $\lambda_1,\lambda_2>0$。记最小特征值为 $m=\min\{\lambda_1,\lambda_2\}>0$。则对任意非零向量 $\mathbf{h}=(h_1,h_2)$ 有：
$$\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}\geq m\|\mathbf{h}\|^2$$

**为什么这个不等式成立？** 对于实对称矩阵 $\mathbf{H}$，存在正交矩阵 $\mathbf{Q}$ 使得 $\mathbf{Q}^\mathsf{T}\mathbf{H}\mathbf{Q}=\text{diag}(\lambda_1,\lambda_2)$。令 $\mathbf{y}=\mathbf{Q}^\mathsf{T}\mathbf{h}$，则 $\|\mathbf{y}\|=\|\mathbf{h}\|$，且：
$$\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}=(\mathbf{Q}\mathbf{y})^\mathsf{T}\mathbf{H}(\mathbf{Q}\mathbf{y})=\mathbf{y}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{H}\mathbf{Q}\mathbf{y}=\lambda_1 y_1^2+\lambda_2 y_2^2\geq m(y_1^2+y_2^2)=m\|\mathbf{h}\|^2$$

**第 3 步：处理 $o(\|\mathbf{h}\|^2)$ 项**。将余项写作：
$$o(\|\mathbf{h}\|^2)=\varepsilon(\mathbf{h})\|\mathbf{h}\|^2$$
其中 $\varepsilon(\mathbf{h})\to 0$ 当 $\mathbf{h}\to\mathbf{0}$（由 $o$ 记号的定义：$\varepsilon(\mathbf{h})=\frac{o(\|\mathbf{h}\|^2)}{\|\mathbf{h}\|^2}$）。

**第 4 步：合并并进行 $\varepsilon$-$\delta$ 论证**。
$$
\begin{aligned}
f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0) &=\frac12\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}+\varepsilon(\mathbf{h})\|\mathbf{h}\|^2 \\
&\geq\frac12 m\|\mathbf{h}\|^2+\varepsilon(\mathbf{h})\|\mathbf{h}\|^2 \\
&=\|\mathbf{h}\|^2\left(\frac{m}{2}+\varepsilon(\mathbf{h})\right)
\end{aligned}
$$

由于 $\varepsilon(\mathbf{h})\to0$，存在 $\delta>0$ 使得当 $0<\|\mathbf{h}\|<\delta$ 时，$|\varepsilon(\mathbf{h})|<\dfrac{m}{4}$。于是：
$$\frac{m}{2}+\varepsilon(\mathbf{h})>\frac{m}{2}-\frac{m}{4}=\frac{m}{4}>0$$

因此当 $0<\|\mathbf{h}\|<\delta$ 时：
$$f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)\geq\|\mathbf{h}\|^2\cdot\frac{m}{4}>0$$
即 $f(\mathbf{P}_0+\mathbf{h})>f(\mathbf{P}_0)$。

**第 5 步：结论**。存在 $\delta>0$，使得对所有满足 $0<\|\mathbf{h}\|<\delta$ 的 $\mathbf{h}$ 有 $f(\mathbf{P}_0+\mathbf{h})>f(\mathbf{P}_0)$。由严格局部极小值的定义（定义12.6），$f$ 在 $\mathbf{P}_0$ 处取得严格局部极小值。证毕。

**注**：
- 若 $\mathbf{H}$ 负定（即 $-\mathbf{H}$ 正定），类似可证 $f$ 在 $\mathbf{P}_0$ 处取得严格局部极大值
- 若 $\mathbf{H}$ 不定（特征值异号），则 $\mathbf{P}_0$ 是鞍点——沿特征向量方向，一个方向上升、另一个方向下降
- 该证明与文件06中的配方法证明本质等价：配方法对应 $2\times2$ 矩阵正定的判别条件 $H>0,A>0$，而这里用特征值给出了更一般的论证

**两种证明路径的对比**：

| 对比项 | 配方法（文件06） | 正定矩阵法（本文件） |
|--------|:--------------:|:-----------------:|
| 核心工具 | 二次型配方 $Q=A(\Delta x+\frac{B}{A}\Delta y)^2+\frac{H}{A}\Delta y^2$ | 特征值下界 $\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}\geq m\|\mathbf{h}\|^2$ |
| 适用范围 | 仅 $2\times2$ 矩阵 | 可推广到 $\mathbb{R}^n$ |
| 证明风格 | 初等代数 | 线性代数 |
| 条件判断 | $H=AC-B^2>0$ 且 $A>0$ | 特征值均 $>0$（等价于顺序主子式 $>0$） |

---

## 5. 例题

### 例题 1：多元微分中值定理与常数命题（A1）

设 $f(x,y)=\ln(1+x^2+y^2)$。

(1) 应用多元微分中值定理（定理12.13），写出 $f(1,2)-f(0,0)$ 的表达式。
(2) 利用多元微分中值定理证明：若可微函数 $g(x,y)$ 在凸区域 $D$ 上满足 $\dfrac{\partial g}{\partial x}\equiv0$ 且 $\dfrac{\partial g}{\partial y}\equiv0$，则 $g$ 在 $D$ 上为常数。

**解 (1)**：

令 $\mathbf{P}_0=(0,0)$，$\mathbf{h}=(1,2)$。由定理12.13，存在 $\theta\in(0,1)$ 使得：
$$f(1,2)-f(0,0)=f_x(\theta,2\theta)\cdot1+f_y(\theta,2\theta)\cdot2$$

计算：
$$f_x=\frac{2x}{1+x^2+y^2},\qquad f_y=\frac{2y}{1+x^2+y^2}$$

代入：
$$f_x(\theta,2\theta)=\frac{2\theta}{1+\theta^2+4\theta^2}=\frac{2\theta}{1+5\theta^2},\quad f_y(\theta,2\theta)=\frac{4\theta}{1+5\theta^2}$$

因此：
$$f(1,2)-f(0,0)=\frac{2\theta}{1+5\theta^2}\cdot1+\frac{4\theta}{1+5\theta^2}\cdot2=\frac{10\theta}{1+5\theta^2}$$

其中 $\theta\in(0,1)$ 是某个特定值。

**验证**：$f(1,2)-f(0,0)=\ln(1+1+4)-\ln1=\ln6\approx1.7918$。而 $\frac{10\theta}{1+5\theta^2}$ 在 $\theta\in(0,1)$ 上的取值范围为 $(0,\frac{10}{6})\approx(0,1.6667)$ ——不包含 $\ln6\approx1.7918$。这不矛盾吗？

注意：多元微分中值定理要求 $f$ 在 $\mathbf{P}_0$ 和 $\mathbf{P}_0+\mathbf{h}$ 的**连线**上一阶连续可导。$f(x,y)=\ln(1+x^2+y^2)$ 确实在 $\mathbb{R}^2$ 上 $C^1$。但 $(0,0)$ 到 $(1,2)$ 的连线上点均为 $(t,2t)$，代入 $f$ 得 $\varphi(t)=\ln(1+5t^2)$。

用一元中值定理：$\varphi(1)-\varphi(0)=\varphi'(\theta)$，$\varphi'(\theta)=\frac{10\theta}{1+5\theta^2}$。

方程 $\frac{10\theta}{1+5\theta^2}=\ln6$ 的解：$10\theta=\ln6(1+5\theta^2)$，即 $5\ln6\cdot\theta^2-10\theta+\ln6=0$。解得 $\theta\approx0.224$ 或 $\theta\approx0.892$，均在 $(0,1)$ 内，均可作为中间点。因此上述计算成立。

**解 (2)**：

若 $\dfrac{\partial g}{\partial x}\equiv0$ 且 $\dfrac{\partial g}{\partial y}\equiv0$，则 $\nabla g\equiv\mathbf{0}$。

任取 $D$ 中两点 $\mathbf{A},\mathbf{B}$，由凸性，线段 $\mathbf{A}+t(\mathbf{B}-\mathbf{A})$（$t\in[0,1]$）包含在 $D$ 内。由定理12.13，存在 $\theta\in(0,1)$ 使：
$$g(\mathbf{B})-g(\mathbf{A})=\nabla g(\mathbf{A}+\theta(\mathbf{B}-\mathbf{A}))\cdot(\mathbf{B}-\mathbf{A})=\mathbf{0}\cdot(\mathbf{B}-\mathbf{A})=0$$

因此 $g(\mathbf{B})=g(\mathbf{A})$，$g$ 在 $D$ 上为常数。证毕。

---

### 例题 2：$e^{x+y}$ 的n阶Taylor展开通项公式（A2）

求 $f(x,y)=e^{x+y}$ 在 $(0,0)$ 处的n阶Taylor展开式（写出通项公式），并写出Lagrange余项。

**解**：

**第 1 步：计算各阶偏导**。对任意多重指标 $\alpha=(\alpha_1,\alpha_2)$，$|\alpha|=k$：
$$\partial^\alpha f=\frac{\partial^k}{\partial x^{\alpha_1}\partial y^{\alpha_2}}e^{x+y}=e^{x+y}$$

（因为 $e^{x+y}$ 对 $x$ 和 $y$ 的任何偏导都等于自身。）

因此在 $(0,0)$ 处：$\partial^\alpha f(0,0)=e^{0}=1$ 对所有 $\alpha$ 成立。

**第 2 步：写出多重指标形式的Taylor公式**。由定理12.16：
$$f(\mathbf{h})=\sum_{|\alpha|\leq n}\frac{1}{\alpha!}\,\partial^\alpha f(\mathbf{0})\,\mathbf{h}^\alpha+R_n$$

代入 $\partial^\alpha f(\mathbf{0})=1$，$\mathbf{h}=(\Delta x,\Delta y)$：
$$f(\Delta x,\Delta y)=\sum_{|\alpha|\leq n}\frac{(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}}{\alpha_1!\alpha_2!}+R_n$$

**第 3 步：用二项式定理化简**。注意求和 $\sum_{|\alpha|\leq n}\frac{(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}}{\alpha_1!\alpha_2!}$ 恰好是 $e^{\Delta x+\Delta y}$ 的Taylor展开。我们来验证：

对固定 $k=|\alpha|$，$\sum_{\alpha_1+\alpha_2=k}\frac{(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}}{\alpha_1!\alpha_2!}=\frac{(\Delta x+\Delta y)^k}{k!}$（由二项式定理）。

因此：
$$\sum_{|\alpha|\leq n}\frac{(\Delta x)^{\alpha_1}(\Delta y)^{\alpha_2}}{\alpha_1!\alpha_2!}=\sum_{k=0}^n\frac{(\Delta x+\Delta y)^k}{k!}$$

所以：
$$e^{\Delta x+\Delta y}=\sum_{k=0}^n\frac{(\Delta x+\Delta y)^k}{k!}+R_n$$

**第 4 步：写出Lagrange余项**。由于 $f$ 的所有 $n+1$ 阶偏导均为 $e^{x+y}$，由定理12.16：
$$R_n=\frac{1}{(n+1)!}(\Delta x\partial_x+\Delta y\partial_y)^{n+1}e^{\theta(\Delta x+\Delta y)},\quad\theta\in(0,1)$$

计算算子作用：
$$(\Delta x\partial_x+\Delta y\partial_y)^{n+1}e^{\theta(\Delta x+\Delta y)}=(\Delta x+\Delta y)^{n+1}e^{\theta(\Delta x+\Delta y)}$$

因此：
$$\boxed{e^{\Delta x+\Delta y}=\sum_{k=0}^n\frac{(\Delta x+\Delta y)^k}{k!}+\frac{(\Delta x+\Delta y)^{n+1}}{(n+1)!}e^{\theta(\Delta x+\Delta y)},\quad\theta\in(0,1)}$$

**换回 $x,y$ 记号**：令 $\Delta x=x$，$\Delta y=y$，以 $(0,0)$ 为展开点：
$$e^{x+y}=\sum_{k=0}^n\frac{(x+y)^k}{k!}+\frac{(x+y)^{n+1}}{(n+1)!}e^{\theta(x+y)}$$

**验证**：当 $n=2$ 时，$e^{x+y}=1+(x+y)+\frac{(x+y)^2}{2}+R_2$，展开得：
$$e^{x+y}=1+x+y+\frac{x^2}{2}+xy+\frac{y^2}{2}+R_2$$
这与直接用二阶多元Taylor公式展开 $e^{x+y}$ 的结果一致。

---

### 例题 3：乘法展开技术——$\sin x\ln(1+y)$ 的三阶展开（A3）

求 $f(x,y)=\sin x\ln(1+y)$ 在 $(0,0)$ 处的三阶Taylor展开式（展开到所有总次数不超过3的项），写出Peano余项。

**解**：

**第 1 步：写出给定的一元Maclaurin展开式**。
$$\sin x=x-\frac{x^3}{6}+o(x^3)$$
$$\ln(1+y)=y-\frac{y^2}{2}+\frac{y^3}{3}+o(y^3)$$

**第 2 步：逐项相乘**。将两个展开式做乘法（暂不考虑余项），展开到所有 $x^py^q$ 且 $p+q\leq3$ 的项：

$$
\begin{aligned}
\sin x\ln(1+y)&=\left(x-\frac{x^3}{6}+o(x^3)\right)\left(y-\frac{y^2}{2}+\frac{y^3}{3}+o(y^3)\right) \\
&=x\cdot y + x\cdot\left(-\frac{y^2}{2}\right) + x\cdot\frac{y^3}{3} \\
&\quad+\left(-\frac{x^3}{6}\right)\cdot y + \left(-\frac{x^3}{6}\right)\cdot\left(-\frac{y^2}{2}\right) + \left(-\frac{x^3}{6}\right)\cdot\frac{y^3}{3} \\
&\quad+ \text{余项项的乘积（更高阶的无穷小）}
\end{aligned}
$$

**第 3 步：按总次数筛选**。

逐项检查每个乘积的总次数 $p+q$：

| 项 | 来源 | 总次数 $p+q$ | 是否保留（$p+q\leq3$） |
|:---|:----:|:-----------:|:---------------------:|
| $xy$ | $x\cdot y$ | $2$ | 保留 |
| $x\cdot(-\frac{y^2}{2})=-\frac12xy^2$ | $x\cdot(-\frac{y^2}{2})$ | $3$ | 保留 |
| $x\cdot\frac{y^3}{3}=\frac13xy^3$ | $x\cdot\frac{y^3}{3}$ | $4$ | 舍弃$\to o(\rho^3)$ |
| $(-\frac{x^3}{6})\cdot y=-\frac16x^3y$ | $(-\frac{x^3}{6})\cdot y$ | $4$ | 舍弃$\to o(\rho^3)$ |
| $(-\frac{x^3}{6})\cdot(-\frac{y^2}{2})=\frac1{12}x^3y^2$ | $(-\frac{x^3}{6})\cdot(-\frac{y^2}{2})$ | $5$ | 舍弃$\to o(\rho^3)$ |
| $(-\frac{x^3}{6})\cdot\frac{y^3}{3}=-\frac1{18}x^3y^3$ | $(-\frac{x^3}{6})\cdot\frac{y^3}{3}$ | $6$ | 舍弃$\to o(\rho^3)$ |

同理，余项与任何项的乘积均属于更高阶的无穷小，统一归入 $o(\rho^3)$。

**第 4 步：写出结果**。
$$\boxed{\sin x\ln(1+y)=xy-\frac12xy^2+o(\rho^3)}$$

其中 $\rho=\sqrt{x^2+y^2}$，$o(\rho^3)$ 表示比 $\rho^3$ 更高阶的无穷小量。

**验证**：展开式中没有 $x^2y$ 项（因为 $\sin x$ 展开中没有 $x^2$ 项），也没有 $y^3$ 项（因为 $x^0$ 项为零），这与直接分析一致。保留的三阶项仅为 $-\frac12xy^2$。

**安全检查**：总次数不超过3的所有可能项为：
$$1,\;x,\;y,\;x^2,\;xy,\;y^2,\;x^3,\;x^2y,\;xy^2,\;y^3$$

在本题中：
- 常数项和 $x$ 项：$\sin x$ 无常数项，$\ln(1+y)$ 的常数项为0 → 无
- $y$ 项：$\sin x$ 展开第一项为 $x$，与 $\ln(1+y)$ 的常数项0相乘得0 → 无
- $x^2$ 项：$\sin x$ 无 $x^2$ 项 → 无
- $xy$ 项：保留 ✓
- $y^2$ 项：$\sin x$ 展开第一项为 $x$，与 $\ln(1+y)$ 的 $-\frac{y^2}{2}$ 相乘得 $-\frac12xy^2$（非纯 $y^2$）→ 无纯 $y^2$ 项
- $x^3$ 项：$\sin x$ 的 $-\frac{x^3}{6}$ 与 $\ln(1+y)$ 的常数项0相乘 → 无
- $x^2y$ 项：$\sin x$ 无 $x^2$ 项 → 无
- $xy^2$ 项：$-\frac12xy^2$ ✓
- $y^3$ 项：$x\cdot\frac{y^3}{3}$ 得 $\frac13xy^3$（总次数4，舍弃）

因此最终结果正确。

---

### 例题 4：二阶Taylor近似计算 $8.96^{2.03}$（A4）

利用二阶Taylor公式近似计算 $8.96^{2.03}$。

**解**：

**第 1 步：选函数与基点**。令 $f(x,y)=(9+x)^{2+y}$，选择展开点 $(x_0,y_0)=(0,0)$，则：
$$f(0,0)=9^2=81$$
此时 $\Delta x=8.96-9=-0.04$，$\Delta y=2.03-2=0.03$，增量足够小。

**第 2 步：计算各阶偏导在 $(0,0)$ 处的值**。

先将 $f(x,y)=(9+x)^{2+y}$ 改写为指数形式以便求导：
$$f(x,y)=e^{(2+y)\ln(9+x)}$$

**一阶偏导**：
$$f_x=e^{(2+y)\ln(9+x)}\cdot\frac{2+y}{9+x}=(9+x)^{2+y}\cdot\frac{2+y}{9+x}$$
$$f_x(0,0)=81\cdot\frac{2}{9}=18$$

$$f_y=e^{(2+y)\ln(9+x)}\cdot\ln(9+x)=(9+x)^{2+y}\ln(9+x)$$
$$f_y(0,0)=81\cdot\ln9=81\times2.1972=177.9732$$

**二阶偏导**：

$f_{xx}$：对 $f_x=(9+x)^{2+y}\cdot\dfrac{2+y}{9+x}$ 求 $x$ 的偏导。注意 $f_x$ 可写为 $(2+y)(9+x)^{1+y}$，于是：
$$f_{xx}=(2+y)(1+y)(9+x)^{y}$$
$$f_{xx}(0,0)=2\cdot1\cdot9^0=2$$

$f_{xy}$：对 $f_x=(2+y)(9+x)^{1+y}$ 求 $y$ 的偏导。
$$f_{xy}=1\cdot(9+x)^{1+y}+(2+y)(9+x)^{1+y}\ln(9+x)$$
（使用乘积法则：$\frac{\partial}{\partial y}(2+y)=1$，$\frac{\partial}{\partial y}(9+x)^{1+y}=(9+x)^{1+y}\ln(9+x)$）
$$f_{xy}=(9+x)^{1+y}\big(1+(2+y)\ln(9+x)\big)$$
$$f_{xy}(0,0)=9^1\cdot\big(1+2\ln9\big)=9\times(1+2\times2.1972)=9\times5.3944=48.5496$$

$f_{yy}$：对 $f_y=(9+x)^{2+y}\ln(9+x)$ 求 $y$ 的偏导。
$$f_{yy}=(9+x)^{2+y}\big(\ln(9+x)\big)^2$$
$$f_{yy}(0,0)=81\times(\ln9)^2=81\times4.8283=391.0923$$

**第 3 步：写出二阶Taylor展开式**（使用Peano余项形式，略去余项）：
$$
\begin{aligned}
f(\Delta x,\Delta y)&\approx f(0,0)+f_x\Delta x+f_y\Delta y \\
&\quad+\frac12\big(f_{xx}\Delta x^2+2f_{xy}\Delta x\Delta y+f_{yy}\Delta y^2\big)
\end{aligned}
$$

代入数值：
$$
\begin{aligned}
8.96^{2.03}&\approx81+18\times(-0.04)+177.9732\times0.03 \\
&\quad+\frac12\big[2\times(-0.04)^2+2\times48.5496\times(-0.04)\times0.03+391.0923\times0.03^2\big]
\end{aligned}
$$

**第 4 步：逐项计算**。

线性部分：
$$81-0.72+5.339196=85.619196$$

计算二阶项的各分量：
$$f_{xx}\Delta x^2=2\times0.0016=0.0032$$
$$2f_{xy}\Delta x\Delta y=2\times48.5496\times(-0.0012)=-0.116519$$
$$f_{yy}\Delta y^2=391.0923\times0.0009=0.351983$$

二阶项之和：
$$\frac12(0.0032-0.116519+0.351983)=\frac12\times0.238664=0.119332$$

最终近似值：
$$85.619196+0.119332=85.738528\approx\boxed{85.74}$$

**对比一阶线性近似**（仅使用前两项）：
$$81-0.72+5.339196=85.6192$$
二阶近似的修正项 $+0.1193$ 使结果更精确。

**真实值验证**：$8.96^{2.03}=e^{2.03\ln8.96}$。$\ln8.96\approx2.192926$，$2.03\times2.192926\approx4.45164$，$e^{4.45164}\approx85.78$。二阶近似的绝对误差约为 $0.04$，相对误差约 $0.05\%$，精度已相当高。

---

### 例题 5：利用正定Hessian证明极值充分条件（A5）

设 $f(x,y)=x^2+4xy+5y^2-2x-6y$。

(1) 求 $f$ 的驻点。
(2) 计算Hessian矩阵并判断其正定性。
(3) 利用二阶Taylor公式和正定矩阵方法证明该驻点是严格局部极小值。

**解**：

**(1) 求驻点**：
$$f_x=2x+4y-2=0,\qquad f_y=4x+10y-6=0$$

解方程组：由第一式得 $2x=2-4y$，$x=1-2y$。代入第二式：
$$4(1-2y)+10y-6=4-8y+10y-6=2y-2=0\Rightarrow y=1$$
$x=1-2= -1$。唯一驻点：$\mathbf{P}_0=(-1,1)$。

**(2) 计算Hessian矩阵**：
$$f_{xx}=2,\quad f_{xy}=4,\quad f_{yy}=10$$
$$\mathbf{H}=\begin{pmatrix}2&4\\4&10\end{pmatrix}$$

判断正定性：顺序主子式 $D_1=2>0$，$D_2=2\times10-4^2=20-16=4>0$。因此 $\mathbf{H}$ 正定。最小特征值可由特征方程 $\det(\mathbf{H}-\lambda\mathbf{I})=(2-\lambda)(10-\lambda)-16=0$ 解得：$\lambda^2-12\lambda+4=0$，$\lambda=6\pm\sqrt{32}$，最小特征值 $m=6-4\sqrt{2}\approx0.343$。

**(3) 用Taylor公式证明严格局部极小值**：

由于 $f\in C^\infty$（多项式函数），可在 $\mathbf{P}_0$ 处展开。在 $\mathbf{P}_0$ 处 $\nabla f=0$，由定理12.15（二阶Peano余项）：
$$f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)=\frac12\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}+o(\|\mathbf{h}\|^2)$$

由正定性：
$$\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}\geq m\|\mathbf{h}\|^2,\quad m=6-4\sqrt{2}>0$$

将余项写作 $o(\|\mathbf{h}\|^2)=\varepsilon(\mathbf{h})\|\mathbf{h}\|^2$，其中 $\varepsilon(\mathbf{h})\to0$（当 $\mathbf{h}\to\mathbf{0}$）。于是：
$$
\begin{aligned}
f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0) &\geq \frac12 m\|\mathbf{h}\|^2+\varepsilon(\mathbf{h})\|\mathbf{h}\|^2 \\
&= \|\mathbf{h}\|^2\left(\frac{m}{2}+\varepsilon(\mathbf{h})\right)
\end{aligned}
$$

取 $\delta>0$ 使得 $0<\|\mathbf{h}\|<\delta$ 时 $|\varepsilon(\mathbf{h})|<\dfrac{m}{4}$。则对 $0<\|\mathbf{h}\|<\delta$：
$$f(\mathbf{P}_0+\mathbf{h})-f(\mathbf{P}_0)\geq \|\mathbf{h}\|^2\cdot\frac{m}{4}>0$$

因此 $f$ 在 $\mathbf{P}_0$ 处取得严格局部极小值。

**计算极小值**：
$$f(-1,1)=(-1)^2+4(-1)(1)+5(1)^2-2(-1)-6(1)=1-4+5+2-6=-2$$

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 多元Taylor公式的Lagrange余项与一元完全一样，只需把 $f^{(n+1)}(\xi)$ 换成 $\partial^\alpha f(\xi)$ | Lagrange余项在多元中是所有 $|\alpha|=n+1$ 项的加权和：$\sum_{|\alpha|=n+1}\frac{1}{\alpha!}\partial^\alpha f(\xi)\mathbf{h}^\alpha$，不是单个偏导项 |
| 辅助函数 $\varphi(t)=f(x_0+t\Delta x,y_0+t\Delta y)$ 的各阶导数就是对应偏导的简单替换 | $\varphi^{(k)}(t)$ 需要通过算子 $(\Delta x\partial_x+\Delta y\partial_y)^kf$ 展开，涉及多个偏导项的组合（二项式系数），不是单一偏导 |
| 两个一元展开式相乘时，保留所有项即可 | 必须按"总次数"（各单项中 $x$ 和 $y$ 的指数之和）筛选。总次数 $>n$ 的项应归入 $o(\rho^n)$ |
| $a^b$ 型函数的偏导可直接套用 $(x^a)'=ax^{a-1}$ | $a^b$ 中底数和指数都是变量，需转化为 $e^{b\ln a}$ 后用链式法则求导；或先固定一变量当常数再对另一变量求导 |
| 正定矩阵证明极值时可以直接用 $h^\mathsf{T}Hh > 0$ 得出结论 | 必须处理 $o(\|h\|^2)$ 项：写成 $\varepsilon(h)\|h\|^2$ 后通过 $\varepsilon(h)\to0$ 和特征值下界结合，用 $\varepsilon$-$\delta$ 论证才能得到严格不等式 |
| 总次数 $p+q$ 等于单项 $x^py^q$ 的"阶数"，与一元展开中 $x^k$ 的阶数 $k$ 含义相同 | 在多元中，$x^py^q$ 作为 $x,y\to0$ 时的无穷小，其阶数依赖于路径。$x^2$ 沿 $x=y$ 方向是 $O(\rho^2)$，但 $xy$ 沿同一方向也是 $O(\rho^2)$。$p+q$ 是"总次数"而非严格阶数，但在各方向同阶的假设下，$p+q$ 越大项越小 |
| 梯度为零的常数证明不需要凸区域条件 | 凸区域条件保证了连接任意两点的线段完全在区域内，从而中值定理可应用。若 $D$ 非凸，可能存在两点间的线段跑出 $D$，此时定理12.13不适用 |

### 检查点

- [ ] 能否构造辅助函数 $\varphi(t)=f(\mathbf{P}_0+t\mathbf{h})$ 并用链式法则求其前二阶导数？
- [ ] 能否写出带Lagrange余项的一阶Taylor公式（多元微分中值定理）并给出完整证明？
- [ ] 能否用微分中值定理证明 $\nabla f\equiv\mathbf{0}\Rightarrow f\equiv\text{常数}$？
- [ ] 能否用多重指标记号写出n阶Taylor公式的通项 $\sum_{|\alpha|\leq n}\frac{1}{\alpha!}\partial^\alpha f(\mathbf{P}_0)\mathbf{h}^\alpha$？
- [ ] 能否解释算子 $(\Delta x\partial_x+\Delta y\partial_y)^k$ 的展开与多重指标求和的关系？
- [ ] 给定两个一元展开式，能否通过乘法得到二元展开式并按总次数筛选保留项？
- [ ] 能否按四步法利用二阶Taylor公式进行数值近似计算？
- [ ] 能否用 $\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}\geq m\|\mathbf{h}\|^2$ 和 $o(\|\mathbf{h}\|^2)=\varepsilon(\mathbf{h})\|\mathbf{h}\|^2$ 证明正定Hessian处的严格局部极小性？

---

## 练习题

### 基础巩固

**1.** 设 $f(x,y)=x^2y+2xy^2$，$\mathbf{P}_0=(1,1)$，$\mathbf{h}=(\Delta x,\Delta y)$。

(1) 构造辅助函数 $\varphi(t)=f(1+t\Delta x,1+t\Delta y)$，求 $\varphi'(0)$ 和 $\varphi''(0)$。
(2) 写出 $f$ 在 $(1,1)$ 处的带Peano余项的二阶Taylor公式。

<details><summary>参考答案</summary>

**(1)** 构造 $\varphi(t)=f(1+t\Delta x,1+t\Delta y)$。

求一阶偏导：$f_x=2xy+2y^2$，$f_y=x^2+4xy$。

$$f_x(1,1)=2+2=4,\quad f_y(1,1)=1+4=5$$
$$\varphi'(0)=f_x(1,1)\Delta x+f_y(1,1)\Delta y=4\Delta x+5\Delta y$$

二阶偏导：$f_{xx}=2y$，$f_{xy}=2x+4y$，$f_{yy}=4x$。

$$f_{xx}(1,1)=2,\quad f_{xy}(1,1)=2+4=6,\quad f_{yy}(1,1)=4$$
$$\varphi''(0)=2\Delta x^2+2\cdot6\Delta x\Delta y+4\Delta y^2=2\Delta x^2+12\Delta x\Delta y+4\Delta y^2$$

**(2)** 二阶Taylor公式（Peano余项）：
$$f(1+\Delta x,1+\Delta y)=f(1,1)+(4\Delta x+5\Delta y)+\frac12(2\Delta x^2+12\Delta x\Delta y+4\Delta y^2)+o(\|\mathbf{h}\|^2)$$

其中 $f(1,1)=1+2=3$。整理得：
$$f(1+\Delta x,1+\Delta y)=3+4\Delta x+5\Delta y+\Delta x^2+6\Delta x\Delta y+2\Delta y^2+o(\Delta x^2+\Delta y^2)$$

</details>

---

**2.** 利用多元微分中值定理（一阶Lagrange余项）证明：若 $f(x,y)$ 在凸区域 $D$ 上满足 $|f_x(x,y)|\leq M$ 和 $|f_y(x,y)|\leq M$ 对一切 $(x,y)\in D$ 成立，则对 $D$ 中任意两点 $A,B$，有 $|f(B)-f(A)|\leq 2M\|B-A\|$。

<details><summary>参考答案</summary>

设 $A=(x_1,y_1)$，$B=(x_2,y_2)$。由凸性，连接 $A$ 和 $B$ 的线段包含在 $D$ 内。令 $\mathbf{h}=B-A=(\Delta x,\Delta y)$。

由定理12.13（多元微分中值定理），存在 $\theta\in(0,1)$ 使得：
$$f(B)-f(A)=f_x(A+\theta\mathbf{h})\Delta x+f_y(A+\theta\mathbf{h})\Delta y$$

取绝对值，利用三角不等式和题设 $|f_x|,|f_y|\leq M$：
$$
\begin{aligned}
|f(B)-f(A)| &\leq |f_x(A+\theta\mathbf{h})|\cdot|\Delta x|+|f_y(A+\theta\mathbf{h})|\cdot|\Delta y| \\
&\leq M(|\Delta x|+|\Delta y|)
\end{aligned}
$$

由 $|\Delta x|+|\Delta y|\leq\sqrt{2}\sqrt{\Delta x^2+\Delta y^2}=\sqrt{2}\|\mathbf{h}\|$（柯西-施瓦茨不等式），得：
$$|f(B)-f(A)|\leq M\sqrt{2}\|\mathbf{h}\|<\sqrt{2}M\|\mathbf{h}\|$$

实际上因为题目要求证明 $|f(B)-f(A)|\leq 2M\|B-A\|$，而 $\sqrt{2}<2$，所以成立。更紧的界是 $|f(B)-f(A)|\leq \sqrt{2}M\|B-A\|$。

</details>

---

### 迁移应用

**3.** 利用乘法展开技术，求 $f(x,y)=e^x\cos y$ 在 $(0,0)$ 处的三阶Taylor展开式（展开到所有总次数不超过3的项），写出Peano余项。

（提示：$e^x=1+x+\dfrac{x^2}{2}+\dfrac{x^3}{6}+o(x^3)$，$\cos y=1-\dfrac{y^2}{2}+\dfrac{y^4}{24}+o(y^4)$。）

<details><summary>参考答案</summary>

**第1步**：写出给定的一元展开式。
$$e^x=1+x+\frac{x^2}{2}+\frac{x^3}{6}+o(x^3)$$
$$\cos y=1-\frac{y^2}{2}+\frac{y^4}{24}+o(y^4)$$

**第2步**：逐项相乘并筛选总次数 $p+q\leq3$ 的项。

$$
\begin{aligned}
e^x\cos y &=\left(1+x+\frac{x^2}{2}+\frac{x^3}{6}+o(x^3)\right)\left(1-\frac{y^2}{2}+o(y^4)\right) \\
&=1\cdot1 + 1\cdot\left(-\frac{y^2}{2}\right) \\
&\quad+ x\cdot1 + x\cdot\left(-\frac{y^2}{2}\right) \\
&\quad+ \frac{x^2}{2}\cdot1 \\
&\quad+ \frac{x^3}{6}\cdot1 \\
&\quad+ \text{更高阶项（总次数$>3$，舍弃）}
\end{aligned}
$$

按总次数逐项检查：
| 项 | 内容 | 总次数 | 保留？ |
|:---|:----:|:-----:|:------:|
| $1\cdot1$ | $1$ | $0$ | 保留 |
| $1\cdot(-\frac{y^2}{2})$ | $-\frac12y^2$ | $2$ | 保留 |
| $1\cdot\frac{y^4}{24}$ | $\frac1{24}y^4$ | $4$ | 舍弃 |
| $x\cdot1$ | $x$ | $1$ | 保留 |
| $x\cdot(-\frac{y^2}{2})$ | $-\frac12xy^2$ | $3$ | **保留** |
| $x\cdot\frac{y^4}{24}$ | $\frac1{24}xy^4$ | $5$ | 舍弃 |
| $\frac{x^2}{2}\cdot1$ | $\frac12x^2$ | $2$ | 保留 |
| $\frac{x^2}{2}\cdot(-\frac{y^2}{2})$ | $-\frac14x^2y^2$ | $4$ | 舍弃 |
| $\frac{x^3}{6}\cdot1$ | $\frac16x^3$ | $3$ | 保留 |
| $\frac{x^3}{6}\cdot(-\frac{y^2}{2})$ | $-\frac1{12}x^3y^2$ | $5$ | 舍弃 |

**第3步**：写出结果。
$$\boxed{e^x\cos y=1+x+\frac12x^2+\frac16x^3-\frac12y^2-\frac12xy^2+o(\rho^3)}$$

其中 $\rho=\sqrt{x^2+y^2}$。

**验证**：总次数不超过3的项共有6项。注意 $\cos y$ 的展开式中 $y^4$ 项虽是一元4阶，但总次数为4，应舍弃。同时 $x^3$ 项总次数为3，保留。

</details>
