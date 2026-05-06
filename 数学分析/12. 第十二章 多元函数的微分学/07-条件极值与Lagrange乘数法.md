# 07. 条件极值与Lagrange乘数法

> 所属章节：第十二章 多元函数的微分学  |  文件序号：07  |  难度：综合
> 常见混淆点：1) 条件极值点处 $\nabla f$ 不一定为零向量——梯度为零是**无条件极值**的必要条件；条件极值点处 $\nabla f$ 与约束函数的梯度 $\nabla g$ **平行**；2) Lagrange 乘子 $\lambda$ 不是目标函数的变量——它是辅助变量，引入乘子是为了将约束极值化为等价的无条件驻点问题

## 1. 学习目标与先修前置

### 学习目标
- 理解条件极值问题的精确定义，区分条件极值与无条件极值
- 掌握代入消元法——当约束可显式解出某变量时，将条件极值化为一元无条件极值
- 理解 Lagrange 乘数法的核心思想：在约束极值点处 $\nabla f \parallel \nabla g$
- 能从一个约束的 Lagrange 乘数法写出方程组并求解（二元函数和三元函数情形）
- 能从多个约束的 Lagrange 乘数法写出方程组并求解（引入多个乘子）
- 掌握 Lagrange 方程组的标准求解策略：消去乘子法、分类讨论法、利用约束几何意义
- 理解条件极值的存在性论证：连续函数在有界闭集（紧集）上必有最值
- 掌握 Lagrange 乘数法的标准化五步求解流程

### 先修知识
- 文件 02（第十二章）：梯度 $\nabla f$ 的定义（定义 12.5）、梯度与等值线垂直的性质（定理 12.4）
- 文件 03（第十二章）：可微性与全微分（定义 12.3），用于方向导数的推导和 Taylor 展开
- 文件 04（第十二章）：多元复合函数的链式法则（定理 12.7），用于隐函数代入后的求导
- 文件 05（第十二章）：隐函数存在定理（定理 12.9），用于理解约束 $g(x,y)=0$ 在局部确定隐函数
- 文件 06（第十二章）：无条件极值的定义（定义 12.6）、极值的必要条件 $\nabla f=0$（定理 12.11）
- 一元函数 Fermat 引理（第五章）：可导函数在极值点处导数为零

---

## 2. 背景与应用场景

**核心问题**：在文件 06 中，我们学习了如何寻找多元函数的无条件极值——自变量 $x,y$ 在整个 $\mathbb{R}^2$ 上自由变化。但现实世界中的优化问题几乎总是带有**约束条件**：

- **成本最小化**：给定体积 $V$，设计一个圆柱形罐子使其表面积最小。目标函数 $S(r,h)=2\pi r^2+2\pi rh$，约束条件 $\pi r^2h=V$——半径 $r$ 和高度 $h$ 不能自由选择，必须满足体积约束。
- **投资组合优化**：在预期收益不低于某水平的前提下，最小化投资组合的风险。资产权重之和必须为 1（约束条件），不能自由分配。
- **最大熵原理**：在给定均值的条件下，寻找使信息熵最大的概率分布。概率之和为 1 且均值固定，这是带有两个约束的优化问题。
- **物理中的约束运动**：质点在曲面 $g(x,y,z)=0$ 上运动，求势能 $V(x,y,z)$ 的极小值——质点位置被限制在曲面上，不能自由移动。

**条件极值与无条件极值的对比**：

| 对比项 | 无条件极值 | 条件极值 |
|--------|-----------|----------|
| 自变量变化范围 | $\mathbb{R}^n$ 中的开集 | 约束集合 $\{x: g_i(x)=0\}$ 的子集 |
| 极值处梯度 | $\nabla f = \mathbf{0}$ | $\nabla f$ 与 $\nabla g$ 平行（一个约束时） |
| 典型方法 | 计算 $\nabla f=0$ 得驻点 | Lagrange 乘数法或代入消元法 |

**与一元函数的类比**：一元函数的条件极值就是"在满足某个条件 $g(x)=0$ 的前提下求 $f(x)$ 的极值"。当只有一个自变量时，条件 $g(x)=0$ 通常直接确定了 $x$ 的值（如果 $g$ 是单射），因此没有自由度。多元函数则不同——约束 $g(x,y)=0$ 通常是一条曲线，自变量的自由度从 2 降为 1，仍然有优化空间。

**本章的知识链位置**：文件 06 研究了**无条件**极值的 Hessian 判别法。本文件研究的是**条件极值**——自变量必须满足一个或多个约束方程。条件极值是多元微积分中最具实用价值的部分之一，支撑着工程优化、经济学和物理学中的大量应用。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 说明 |
|------|------|------|
| $f(\mathbf{x})$ | 目标函数（需优化的函数） | 可以是二元 $f(x,y)$ 或三元 $f(x,y,z)$ |
| $g(\mathbf{x})=0$ | 约束条件（方程形式） | $g(x,y)=0$ 为平面曲线，$g(x,y,z)=0$ 为空间曲面 |
| $\nabla f$ | 目标函数的梯度 | 指向 $f$ 增长最快的方向 |
| $\nabla g$ | 约束函数的梯度 | 垂直于约束曲面 |
| $\lambda$ | Lagrange 乘子（一个约束时） | 满足 $\nabla f = \lambda\nabla g$ 的标量 |
| $\lambda,\mu$ | Lagrange 乘子（两个约束时） | 满足 $\nabla f = \lambda\nabla g_1 + \mu\nabla g_2$ |
| $L$ | Lagrange 函数 | $L = f - \sum \lambda_i g_i$，其驻点对应条件极值的候选点 |
| $S$ | 可行域（约束集合） | $S = \{\mathbf{x}: g_i(\mathbf{x})=0, i=1,\ldots,m\}$ |

### 3.1 条件极值的定义（A-1）

**定义 12.13（条件极值）**：设 $f(\mathbf{x})$ 是定义在开集 $D \subseteq \mathbb{R}^n$ 上的函数，$g_1(\mathbf{x}),\ldots,g_m(\mathbf{x})$ 是 $D$ 上的函数（$m < n$）。考虑约束集合
$$S = \{\mathbf{x} \in D : g_1(\mathbf{x}) = 0,\; g_2(\mathbf{x}) = 0,\; \ldots,\; g_m(\mathbf{x}) = 0\}$$

若存在点 $\mathbf{x}_0 \in S$ 及其邻域 $U$，使得对所有 $\mathbf{x} \in S \cap U$ 有：
- $f(\mathbf{x}) \leq f(\mathbf{x}_0)$，则称 $\mathbf{x}_0$ 为 $f$ 在约束 $S$ 下的**条件极大值点**（constrained maximum point），$f(\mathbf{x}_0)$ 为条件极大值；
- $f(\mathbf{x}) \geq f(\mathbf{x}_0)$，则称 $\mathbf{x}_0$ 为 $f$ 在约束 $S$ 下的**条件极小值点**（constrained minimum point），$f(\mathbf{x}_0)$ 为条件极小值。

统称为**条件极值**（constrained extremum）或**条件极值点**。

**与无条件极值的本质区别**：
- 无条件极值：在 $\mathbf{x}_0$ 的**整个邻域**内比较 $f$ 的值——自变量可以在任意方向上自由移动。
- 条件极值：仅在**约束集合 $S$ 上**比较 $f$ 的值——自变量被限制在 $S$ 上，只能沿约束曲面的切方向移动。
- 因此，条件极值点处 $\nabla f$ **不一定为零**——因为自变量被限制住了，$f$ 可能在 $S$ 上达到最大/最小，但沿垂直于 $S$ 的方向仍有上升空间。

**几何直观**：将 $z=f(x,y)$ 看作一张曲面，约束 $g(x,y)=0$ 是 $xy$ 平面上的曲线。条件极值相当于在这条曲线上方 "架起" 的曲面弧线中找最高点和最低点——不是在整个曲面上找，而是在曲面上方的某条"路径"上找。

### 3.2 一个约束时的梯度平行关系（A-3）

条件极值理论的核心几何直观来自以下观察。

设 $f(x,y)$ 在约束 $g(x,y)=0$ 下在 $P_0=(x_0,y_0)$ 处取得条件极值。考虑约束集合 $S = \{(x,y): g(x,y)=0\}$，它是 $xy$ 平面中的一条曲线。在 $P_0$ 处，曲线 $S$ 有切向量 $\mathbf{v}$，且 $\nabla g(P_0)$ 垂直于 $S$（由定理 12.4，梯度与等值线正交），因此 $\nabla g(P_0) \cdot \mathbf{v} = 0$。

现取 $S$ 上的任意光滑参数曲线 $\gamma(t) = (x(t), y(t))$，满足 $\gamma(0) = P_0$，$\gamma'(0) = \mathbf{v}$。由于 $\gamma(t)$ 始终位于 $S$ 上，故 $g(\gamma(t)) \equiv 0$，链式法则给出：
$$0 = \frac{d}{dt}g(\gamma(t))\big|_{t=0} = \nabla g(P_0) \cdot \gamma'(0) = \nabla g(P_0) \cdot \mathbf{v}$$

现在考察 $f$ 在该曲线上的限制 $h(t) = f(\gamma(t))$。由于 $P_0$ 是 $f$ 在 $S$ 上的条件极值点，$h(t)$ 作为一元函数在 $t=0$ 处取得极值。由 Fermat 引理，$h'(0) = 0$：
$$0 = h'(0) = \frac{d}{dt}f(\gamma(t))\big|_{t=0} = \nabla f(P_0) \cdot \gamma'(0) = \nabla f(P_0) \cdot \mathbf{v}$$

因此 $\nabla f(P_0)$ 也与 $S$ 的任意切向量 $\mathbf{v}$ 垂直。由于在二维平面上，与 $S$ 的所有切向量垂直的向量必与 $\nabla g(P_0)$ 平行（两者都平行于 $S$ 的法线方向），我们得到核心结论：

**定理 12.14（条件极值的梯度必要条件）**：设 $f$ 和 $g$ 在 $P_0$ 处可微，$g(P_0)=0$，且 $P_0$ 是 $f$ 在约束 $g=0$ 下的条件极值点。若 $\nabla g(P_0) \neq \mathbf{0}$，则存在实数 $\lambda$ 使得
$$\boxed{\nabla f(P_0) = \lambda \nabla g(P_0)}$$

**说明**：
- $\nabla g(P_0) \neq \mathbf{0}$ 确保约束曲线在 $P_0$ 处有确定的法线方向。若 $\nabla g(P_0)=\mathbf{0}$，则该点是 $g$ 的奇点，约束曲线的行为可能异常。
- $\lambda$ 称为 **Lagrange 乘子**（Lagrange multiplier），其值在求解过程中确定。
- 定理 12.14 是**必要条件**——满足 $\nabla f = \lambda \nabla g$ 的点是候选点，不一定是真正的极值点。

**多约束情形的推广**：当有两个约束 $g_1=0, g_2=0$ 时（以三元函数为例），$S$ 是两条曲面的交线（一条空间曲线）。$S$ 的切方向 $\mathbf{v}$ 同时垂直于 $\nabla g_1$ 和 $\nabla g_2$。条件极值点处 $\nabla f \cdot \mathbf{v} = 0$ 对任意切向量 $\mathbf{v}$ 成立，这意味着 $\nabla f$ 位于 $\nabla g_1$ 和 $\nabla g_2$ 张成的平面中。因此存在 $\lambda,\mu$ 使得：
$$\boxed{\nabla f = \lambda \nabla g_1 + \mu \nabla g_2}$$

**几何解释图示**：
```
一个约束 (g(x,y)=0)：              两个约束 (g₁=0, g₂=0)：
                                     
   ∇g          ∇f                   ∇g₁  ∇g₂
    ↑           ↑                    ↖   ↗
     \         /                      \ /
      \ 平行  /                        × ← ∇f 在张成空间中
       \     /                        / \
        \   /                        /   \
         \ /                      ∇f 在 ∇g₁,∇g₂ 张成的平面中
          P₀
          |
         g=0  (曲线)
```

---

## 4. 原理与方法

### 4.1 代入消元法（A-2）

代入消元法是最直观的条件极值求解方法——当约束 $g(x,y)=0$ 能显式解出某个变量时（如 $y = \varphi(x)$），将其代入目标函数，将条件极值转化为一元函数的无条件极值。

**方法步骤**：

| 步骤 | 操作 | 说明 |
|------|------|------|
| 第 1 步：解约束 | 从 $g(x,y)=0$ 中解出 $y = \varphi(x)$（或 $x = \psi(y)$） | 需要 $g$ 在局部可显式解出某变量 |
| 第 2 步：代入目标 | 构造一元函数 $\Phi(x) = f(x, \varphi(x))$ | 条件极值问题化为一元无条件极值 |
| 第 3 步：求驻点 | 解 $\Phi'(x) = 0$，得候选点的 $x$ 坐标 | 应用一元 Fermat 引理 |
| 第 4 步：判定类型 | 用 $\Phi''(x)$ 的符号判定极大/极小 | $\Phi''(x)<0$ 为极大，$\Phi''(x)>0$ 为极小 |

**条件极值 $f(x,y)$ 在约束 $g(x,y)=0$ 下的代入消元法流程图**：
```
g(x,y)=0  →  y = φ(x)  →  Φ(x)=f(x,φ(x))  →  Φ'(x)=0  →  候选点
                                        ↓
                                   Φ''(x)>0 → 极小
                                   Φ''(x)<0 → 极大
```

**代入消元法的适用条件与局限性**：
- **适用**：约束 $g(x,y)=0$ 可显式解出 $y=\varphi(x)$ 或 $x=\psi(y)$。
- **局限性**：
  - 很多约束无法全局解出某个变量（如 $x^2+y^2=1$ 给出 $y=\pm\sqrt{1-x^2}$，需分段处理）。
  - 对于三元以上或两个以上约束，代入变得极其繁琐甚至不可行。
  - Lagrange 乘数法能统一处理这些复杂情况。

**例（代入消元法演示）**：求 $f(x,y)=x^2-y^2$ 在约束 $x+y=2$ 下的条件极值。

**解**：由 $x+y=2$ 得 $y=2-x$。代入目标函数：
$$\Phi(x) = f(x,2-x) = x^2 - (2-x)^2 = x^2 - (4-4x+x^2) = 4x-4$$

求导：$\Phi'(x) = 4$，令 $\Phi'(x)=0$ 得 $4=0$，无解。

这说明该条件极值问题**没有驻点**——在直线 $y=2-x$ 上，$f(x,y)=4x-4$ 是 $x$ 的线性函数（无界）。实际上，沿该直线 $x\to\infty$ 时 $f\to\infty$，$x\to-\infty$ 时 $f\to-\infty$，无极值。

这个例子说明：**条件极值问题不一定有解**——需要存在性论证（见 4.5 节）。

### 4.2 一个约束的 Lagrange 乘数法（A-3, A-4, A-5）

#### 4.2.1 二元函数情形（A-4）

由定理 12.14，在条件极值点处 $\nabla f = \lambda \nabla g$ 且 $g=0$。将这两个条件合并，定义 **Lagrange 函数**（Lagrangian function）：

$$\boxed{L(x,y,\lambda) = f(x,y) - \lambda g(x,y)}$$

令 $L$ 对所有自变量（包括 $\lambda$）的偏导为零，恰好给出条件极值的必要条件：

$$\begin{cases}
L_x = \dfrac{\partial L}{\partial x} = f_x(x,y) - \lambda g_x(x,y) = 0 \\[6pt]
L_y = \dfrac{\partial L}{\partial y} = f_y(x,y) - \lambda g_y(x,y) = 0 \\[6pt]
L_\lambda = \dfrac{\partial L}{\partial \lambda} = -g(x,y) = 0
\end{cases}$$

这就是 **Lagrange 方程组**（Lagrange system）。

**为什么这样构造**：
- $L_x=0$ 和 $L_y=0$ 等价于 $\nabla f = \lambda \nabla g$。
- $L_\lambda = 0$ 等价于 $g=0$（约束本身）。
- 通过引入 $\lambda$，将**条件极值问题**转化为 $L$ 的**无条件驻点问题**——尽管 $L$ 多了 $\lambda$ 这个变量，但方程组总变量数和方程数相等，可解。

**Lagrange 乘数法的特点**：
- 不需要显式解出 $y=\varphi(x)$，直接处理约束方程。
- 通过引入乘子 $\lambda$，将条件极值转化为驻点问题。
- 三方程三未知数 $(x,y,\lambda)$，**变量数 = 方程数**。

#### 4.2.2 三元函数情形（A-5）

对于 $f(x,y,z)$ 在约束 $g(x,y,z)=0$ 下的条件极值，推广完全类似：

**Lagrange 函数**：
$$L(x,y,z,\lambda) = f(x,y,z) - \lambda g(x,y,z)$$

**Lagrange 方程组**：
$$\begin{cases}
L_x = f_x - \lambda g_x = 0 \\[4pt]
L_y = f_y - \lambda g_y = 0 \\[4pt]
L_z = f_z - \lambda g_z = 0 \\[4pt]
L_\lambda = -g = 0
\end{cases}$$

四方程四未知数 $(x,y,z,\lambda)$，变量数 = 方程数。

**一般形式**（$n$ 元函数，$1$ 个约束）：
$$L(x_1,\ldots,x_n,\lambda) = f(x_1,\ldots,x_n) - \lambda g(x_1,\ldots,x_n)$$
方程组共 $n+1$ 个方程：
$$\begin{cases}
\dfrac{\partial L}{\partial x_i} = \dfrac{\partial f}{\partial x_i} - \lambda \dfrac{\partial g}{\partial x_i} = 0,\quad i=1,\ldots,n \\[6pt]
\dfrac{\partial L}{\partial \lambda} = -g = 0
\end{cases}$$

### 4.3 两个约束的 Lagrange 乘数法（三元函数情形）（A-6）

当有两个约束 $g_1(x,y,z)=0$ 和 $g_2(x,y,z)=0$ 时，由 3.2 节的推广，条件极值点处：
$$\nabla f = \lambda \nabla g_1 + \mu \nabla g_2$$

引入两个 Lagrange 乘子 $\lambda,\mu$，构造 Lagrange 函数：
$$\boxed{L(x,y,z,\lambda,\mu) = f(x,y,z) - \lambda g_1(x,y,z) - \mu g_2(x,y,z)}$$

**Lagrange 方程组**（五方程五未知数）：
$$\begin{cases}
L_x = f_x - \lambda g_{1x} - \mu g_{2x} = 0 \\[4pt]
L_y = f_y - \lambda g_{1y} - \mu g_{2y} = 0 \\[4pt]
L_z = f_z - \lambda g_{1z} - \mu g_{2z} = 0 \\[4pt]
L_\lambda = -g_1 = 0 \\[4pt]
L_\mu = -g_2 = 0
\end{cases}$$

**一般形式**（$n$ 元函数，$m$ 个约束，$m<n$）：
$$L(x_1,\ldots,x_n,\lambda_1,\ldots,\lambda_m) = f - \sum_{j=1}^m \lambda_j g_j$$
方程组共 $n+m$ 个方程（$n$ 个自变量的偏导 + $m$ 个约束条件的偏导）。

### 4.4 方程组的求解策略（A-7）

Lagrange 方程组的求解没有统一的机械算法，但有若干典型策略。

#### 策略一：消去乘子法（最常用）

从 $L_{x_i}=0$ 等方程中解出变量与乘子的关系，代入约束消去乘子。

**典型操作**：如果 $L_x=0$ 和 $L_y=0$ 分别给出 $f_x = \lambda g_x$ 和 $f_y = \lambda g_y$，在 $g_x\neq 0$ 时可得：
$$\lambda = \frac{f_x}{g_x} = \frac{f_y}{g_y}$$

这个等式给出了 $x,y$ 之间的关系，代入约束即可求解。

**例**：$f(x,y)=x^2+y^2$ 在 $x+y=1$ 下求极值。
- $L = x^2+y^2 - \lambda(x+y-1)$
- $L_x = 2x - \lambda = 0 \Rightarrow \lambda = 2x$
- $L_y = 2y - \lambda = 0 \Rightarrow \lambda = 2y$
- 因此 $x=y$，代入 $x+y=1$ 得 $x=y=\frac12$，$\lambda=1$。
- 候选点 $(\frac12,\frac12)$。

#### 策略二：相除消元法

当 $L_x=0$ 和 $L_y=0$ 都有 $\lambda$ 项时，可两式相除消去 $\lambda$，得到 $x,y$ 的关系。

**例**（三元情形）：$f(x,y,z)=x-2y+2z$ 在 $x^2+y^2+z^2=1$ 下：
- $L_x = 1-2\lambda x = 0 \Rightarrow x = \frac{1}{2\lambda}$
- $L_y = -2-2\lambda y = 0 \Rightarrow y = -\frac{1}{\lambda}$
- $L_z = 2-2\lambda z = 0 \Rightarrow z = \frac{1}{\lambda}$

将 $x,y,z$ 都用 $\lambda$ 表示后代入 $x^2+y^2+z^2=1$ 得关于 $\lambda$ 的方程，解之即可。

#### 策略三：分类讨论法

当方程组中出现因式分解时（如 $x\cdot(\text{表达式})=0$），需要分情况讨论。

**操作规范**：
1. 列出所有可能的分支情形（如 $x=0$ 或 $(\text{表达式})=0$）
2. 逐支求解，检查是否自洽
3. 排除不符合约束条件的分支

**典型情形**：从 $L_x=0$ 得 $2x(1+\lambda) - \mu = 0$，从 $L_y=0$ 得 $2y(1+\lambda) - \mu = 0$。相减得 $2(1+\lambda)(x-y)=0$，需分 $\lambda=-1$ 或 $x=y$ 两种情形讨论。

#### 策略四：利用约束的几何意义排除候选点

约束函数本身往往隐含了变量的取值范围（如 $z=x^2+y^2 \geq 0$，$x^2+y^2+z^2=R^2$ 隐含 $-R\leq x,y,z\leq R$），利用这些信息可以排除不合理的候选解。

**例**：当分类讨论得到 $\lambda=-1$ 的一个分支时，代入解得 $z=-\frac12$，但约束 $z=x^2+y^2\geq 0$ 将之排除。这个几何判断避免了无谓的计算。

#### 策略五：利用对称性简化

当目标函数和约束条件具有对称性时（如 $f$ 和 $g$ 在 $x,y$ 的置换下不变），可先假设对称关系 $x=y$ 再代入求解。

**四种策略的关系图**：
```
Lagrange 方程组
      │
      ├─ 消去乘子法 — 最通用，用 λ 表示变量
      │
      ├─ 相除消元法 — 消去 λ 得变量关系
      │
      ├─ 分类讨论法 — 因子分解时用，逐支分析
      │
      └─ 几何排除法 — 利用变量范围剔除不合理解
```

### 4.5 存在性论证（A-8）

Lagrange 乘数法只给出候选点，这些候选点不一定是真正的极值点。如何判断哪些候选点对应实际的最值？

**定理（有界闭集上的最值定理）**：设 $S \subseteq \mathbb{R}^n$ 是有界闭集（紧集），$f: S \to \mathbb{R}$ 是连续函数。则 $f$ 在 $S$ 上必能取到最大值和最小值。

**在条件极值中的应用**：
- 若约束集合 $S = \{\mathbf{x}: g_i(\mathbf{x}) = 0\}$ 是**有界闭集**，则 $f$ 在 $S$ 上既有最大值又有最小值。
- 此时，所有候选点中函数值最大者即为最大值，最小者即为最小值。
- 无需像无条件极值那样逐一判定每个候选点的极值类型——只需比较函数值。

**常见的紧约束集合**：
| 约束集合 | 紧性 | 原因 |
|----------|------|------|
| 球面 $x^2+y^2+z^2=R^2$ | 有界闭集 | 所有坐标有界 $(|x|,|y|,|z|\leq R)$，闭合曲面 |
| 椭圆 $x^2/a^2+y^2/b^2=1$ | 有界闭集 | $|x|\leq a,\ |y|\leq b$，封闭曲线 |
| 圆柱面 $x^2+y^2=R^2,\ z\in[a,b]$ | 有界闭集 | 截面有界，两端封闭 |
| 两曲面交线（若交线有界） | 有界闭集 | 光滑交线若闭合则有界 |

**非紧的情形**：
- 线性约束如 $x+y=1$（一条无限延伸的直线）——不是有界集，最值可能不存在。
- 抛物面 $z=x^2+y^2$ 本身是无界的，但与其他有界约束的交线可能是有界的。

**判定流程**：若约束集合是紧集，则最值存在 → 比较候选点函数值即可。若约束集合不是紧集，需要额外分析函数在"边界"（如无穷远处）的性态。

### 4.6 标准化求解流程——Lagrange 乘数法五步法（A-9）

综合以上内容，求解条件极值的标准流程如下。它与文件 06 中的无条件极值五步法形成对照。

**Lagrange 乘数法五步法**：

| 步骤 | 操作 | 说明 |
|:----:|:-----|:-----|
| 第 1 步：识别问题 | 确定目标函数 $f$ 和约束 $g_i=0$ | 明确哪些是变量，哪些是约束 |
| 第 2 步：构造 Lagrange 函数 | $L = f - \sum_{i=1}^m \lambda_i g_i$ | 一个约束一个乘子 |
| 第 3 步：列出方程组 | 对所有变量（含乘子）求偏导并令其为零 | $\partial L/\partial x_j=0$，$\partial L/\partial \lambda_i=0$ |
| 第 4 步：解方程组 | 使用 4.4 节的策略求解 | 得候选点 $(x_i,\lambda_i)$ |
| 第 5 步：判定极值 | 比较函数值（紧集）或进一步分析 | 若约束集为紧集，直接比较得最值 |

**五步法流程图**：
```
第 1 步：确定 f 和约束 {g_i=0}
         │
         ▼
第 2 步：构造 L = f - Σλᵢgᵢ
         │
         ▼
第 3 步：写出方程组
       ┌─┴──┐
    ∂L/∂xⱼ=0  ∂L/∂λᵢ=0 (= -gᵢ=0)
       └─┬──┘
         │
         ▼
第 4 步：解方程组（消乘子/分类讨论/几何排除）
         │
         ▼
第 5 步：比较候选点 f 值 → 得最值
         （紧集上必存在最值）
```

**与无条件极值五步法（文件 06）的对照**：

| 阶段 | 无条件极值五步法 | 条件极值 Lagrange 五步法 |
|:----:|:-----------------|:------------------------|
| 预备 | 确定 $f$ 的定义域 | 确定 $f$ **和**约束 $g_i=0$ |
| 引入 | 无需新变量 | 引入乘子 $\lambda_i$，构造 $L$ |
| 方程组 | $f_x=0,\ f_y=0$ | $L_{x_j}=0,\ L_{\lambda_i}=0$ |
| 求解 | 解偏导方程组 | 解含乘子的方程组（策略更多元） |
| 判定 | Hessian 判别法（$AC-B^2$） | 紧集上比较函数值；非紧集需分析 |
| 难点 | $H=0$ 退化情形 | 方程组分类讨论、几何排除 |

---

## 5. 例题

### 例题 1（两种方法对照——二元函数、一个约束）

求 $f(x,y) = 2x^2 + y^2$ 在约束 $x + y = 3$ 下的极值。

**解：**

**方法一：代入消元法**

由约束 $x+y=3$ 解出 $y = 3-x$。代入目标函数：
$$\Phi(x) = f(x,3-x) = 2x^2 + (3-x)^2 = 2x^2 + 9 - 6x + x^2 = 3x^2 - 6x + 9$$

求驻点（Fermat 引理）：
$$\Phi'(x) = 6x - 6 = 0 \quad\Rightarrow\quad x = 1$$
对应 $y = 3-1 = 2$。候选点 $(1,2)$。

判定极值类型：
$$\Phi''(x) = 6 > 0$$
因此 $\Phi$ 在 $x=1$ 处取极小值，故条件极小值为：
$$f(1,2) = 2\cdot1^2 + 2^2 = 2+4 = 6$$

**方法二：Lagrange 乘数法**

目标函数 $f(x,y) = 2x^2 + y^2$，约束 $g(x,y) = x+y-3 = 0$。

构造 Lagrange 函数：
$$L(x,y,\lambda) = (2x^2 + y^2) - \lambda(x+y-3)$$

写出方程组：
$$\begin{cases}
L_x = 4x - \lambda = 0 \\[4pt]
L_y = 2y - \lambda = 0 \\[4pt]
L_\lambda = -(x+y-3) = 0
\end{cases}$$

**求解**：由 $L_x=0$ 得 $\lambda = 4x$，由 $L_y=0$ 得 $\lambda = 2y$。因此 $4x = 2y$，即 $y = 2x$。代入 $L_\lambda=0$（即 $x+y=3$）：
$$x + 2x = 3 \quad\Rightarrow\quad 3x = 3 \quad\Rightarrow\quad x = 1,\ y = 2$$

与代入消元法结果一致。$\lambda = 4x = 4$。

**存在性分析**：约束 $x+y=3$ 是一条无限直线，不是有界集。但 $\Phi(x)=3x^2-6x+9$ 是开口向上的二次函数，有全局最小值（$\Phi(x) \to \infty$ 当 $x\to\pm\infty$）。因此该极小值也是全局最小值。$f_{\min} = 6$。

**小结**：代入消元法和 Lagrange 乘数法得到相同的候选点。代入消元法更直接（当约束可解时），Lagrange 乘数法更系统（适用于更复杂的情形）。

---

### 例题 2（三元函数、一个约束、紧集存在性论证）

求函数 $f(x,y,z) = 2x - y + 2z$ 在球面 $x^2 + y^2 + z^2 = 4$ 上的最大值和最小值。

**解：**

**第 1 步：识别问题**

目标函数 $f(x,y,z) = 2x - y + 2z$，约束 $g(x,y,z) = x^2+y^2+z^2 - 4 = 0$。

**第 2 步：构造 Lagrange 函数**
$$L(x,y,z,\lambda) = (2x - y + 2z) - \lambda(x^2+y^2+z^2-4)$$

**第 3 步：列出方程组**
$$\begin{cases}
L_x = 2 - 2\lambda x = 0 \quad\Rightarrow\quad 2\lambda x = 2 \quad\Rightarrow\quad x = \dfrac{1}{\lambda} \\[6pt]
L_y = -1 - 2\lambda y = 0 \quad\Rightarrow\quad 2\lambda y = -1 \quad\Rightarrow\quad y = -\dfrac{1}{2\lambda} \\[6pt]
L_z = 2 - 2\lambda z = 0 \quad\Rightarrow\quad 2\lambda z = 2 \quad\Rightarrow\quad z = \dfrac{1}{\lambda} \\[6pt]
L_\lambda = -(x^2+y^2+z^2-4) = 0
\end{cases}$$

**第 4 步：解方程组**

将 $x,y,z$ 用 $\lambda$ 的表达式代入约束方程：
$$\left(\frac{1}{\lambda}\right)^2 + \left(-\frac{1}{2\lambda}\right)^2 + \left(\frac{1}{\lambda}\right)^2 = 4$$

计算：
$$\frac{1}{\lambda^2} + \frac{1}{4\lambda^2} + \frac{1}{\lambda^2} = 4$$

通分：
$$\frac{4}{4\lambda^2} + \frac{1}{4\lambda^2} + \frac{4}{4\lambda^2} = \frac{9}{4\lambda^2} = 4$$

解得：
$$4\lambda^2 = \frac{9}{4} \quad\Rightarrow\quad \lambda^2 = \frac{9}{16} \quad\Rightarrow\quad \lambda = \pm\frac{3}{4}$$

**分支 1**：$\lambda = \dfrac{3}{4}$
$$x = \frac{1}{3/4} = \frac{4}{3},\quad y = -\frac{1}{2\cdot 3/4} = -\frac{2}{3},\quad z = \frac{1}{3/4} = \frac{4}{3}$$
候选点 $P_1 = \left(\dfrac{4}{3},\;-\dfrac{2}{3},\;\dfrac{4}{3}\right)$。

**分支 2**：$\lambda = -\dfrac{3}{4}$
$$x = \frac{1}{-3/4} = -\frac{4}{3},\quad y = -\frac{1}{2\cdot(-3/4)} = \frac{2}{3},\quad z = \frac{1}{-3/4} = -\frac{4}{3}$$
候选点 $P_2 = \left(-\dfrac{4}{3},\;\dfrac{2}{3},\;-\dfrac{4}{3}\right)$。

**第 5 步：判定极值**

**存在性论证**：约束 $x^2+y^2+z^2=4$ 是以原点为球心、半径为 2 的球面。球面是 $\mathbb{R}^3$ 中的有界闭集（紧集）。$f(x,y,z)=2x-y+2z$ 是连续函数。由有界闭集上的最值定理，$f$ 在该球面上必存在最大值和最小值。

计算候选点的函数值：
$$f(P_1) = 2\cdot\frac{4}{3} - \left(-\frac{2}{3}\right) + 2\cdot\frac{4}{3} = \frac{8}{3} + \frac{2}{3} + \frac{8}{3} = \frac{18}{3} = 6$$
$$f(P_2) = 2\cdot\left(-\frac{4}{3}\right) - \frac{2}{3} + 2\cdot\left(-\frac{4}{3}\right) = -\frac{8}{3} - \frac{2}{3} - \frac{8}{3} = -\frac{18}{3} = -6$$

因此：
$$f_{\max} = 6\quad\text{（在 }P_1\text{ 处取到）},\qquad f_{\min} = -6\quad\text{（在 }P_2\text{ 处取到）}$$

**验证**：由 Cauchy-Schwarz 不等式（或直接观察）：
$$|2x - y + 2z| \leq \sqrt{2^2+(-1)^2+2^2}\cdot\sqrt{x^2+y^2+z^2} = \sqrt{9}\cdot\sqrt{4} = 3\cdot2 = 6$$
在球面上 $x^2+y^2+z^2=4$，故 $|f| \leq 6$，与计算结果一致。

---

### 例题 3（两个约束的分类讨论）

求函数 $f(x,y,z) = x + 2y$ 在约束 $x^2 + y^2 = 5$ 和 $x + y + z = 0$ 下的最大值。

**解：**

**第 1 步：识别问题**

目标函数 $f(x,y,z) = x + 2y$，两个约束：
$$g_1(x,y,z) = x^2 + y^2 - 5 = 0,\quad g_2(x,y,z) = x + y + z = 0$$

**第 2 步：构造 Lagrange 函数**

引入两个乘子 $\lambda,\mu$：
$$L(x,y,z,\lambda,\mu) = (x+2y) - \lambda(x^2+y^2-5) - \mu(x+y+z)$$

**第 3 步：列出方程组**
$$\begin{cases}
L_x = 1 - 2\lambda x - \mu = 0 \quad\Rightarrow\quad \mu = 1 - 2\lambda x \quad①\\[4pt]
L_y = 2 - 2\lambda y - \mu = 0 \quad\Rightarrow\quad \mu = 2 - 2\lambda y \quad②\\[4pt]
L_z = 0 - 0 - \mu = 0 \quad\Rightarrow\quad \mu = 0 \quad③\\[4pt]
L_\lambda = -(x^2+y^2-5) = 0 \quad\Rightarrow\quad x^2+y^2=5 \quad④\\[4pt]
L_\mu = -(x+y+z) = 0 \quad\Rightarrow\quad z = -x-y \quad⑤
\end{cases}$$

**第 4 步：解方程组**

由 ③ 得 $\mu = 0$。代入 ① 和 ②：
$$1 - 2\lambda x = 0 \quad\Rightarrow\quad x = \frac{1}{2\lambda}$$
$$2 - 2\lambda y = 0 \quad\Rightarrow\quad y = \frac{1}{\lambda}$$

代入 ④（$x^2+y^2=5$）：
$$\left(\frac{1}{2\lambda}\right)^2 + \left(\frac{1}{\lambda}\right)^2 = 5 \quad\Rightarrow\quad \frac{1}{4\lambda^2} + \frac{1}{\lambda^2} = 5$$

通分：
$$\frac{1}{4\lambda^2} + \frac{4}{4\lambda^2} = \frac{5}{4\lambda^2} = 5 \quad\Rightarrow\quad 4\lambda^2 = 1 \quad\Rightarrow\quad \lambda^2 = \frac{1}{4} \quad\Rightarrow\quad \lambda = \pm\frac{1}{2}$$

**分支 1**：$\lambda = \dfrac{1}{2}$
$$x = \frac{1}{2\cdot 1/2} = 1,\quad y = \frac{1}{1/2} = 2,\quad z = -x-y = -1-2 = -3$$
候选点 $P_1 = (1,2,-3)$。

**分支 2**：$\lambda = -\dfrac{1}{2}$
$$x = \frac{1}{2\cdot(-1/2)} = -1,\quad y = \frac{1}{-1/2} = -2,\quad z = -x-y = 1+2 = 3$$
候选点 $P_2 = (-1,-2,3)$。

**第 5 步：判定极值**

**存在性论证**：约束 $x^2+y^2=5$ 是圆柱面（在 $xy$ 平面上的投影是半径为 $\sqrt{5}$ 的圆），$x+y+z=0$ 是平面。两者的交线是空间椭圆，是**有界闭集**（椭圆是闭合有界曲线）。$f$ 连续，因此最值存在。

计算函数值：
$$f(P_1) = 1 + 2\cdot2 = 1+4 = 5$$
$$f(P_2) = -1 + 2\cdot(-2) = -1-4 = -5$$

因此最大值为 $5$（在 $P_1$ 处取到），最小值为 $-5$（在 $P_2$ 处取到）。

**验证**：本题的特殊结构（$L_z=0$ 直接给出 $\mu=0$）使得两个约束问题本质上退化为一个约束问题（$\mu=0$ 后，第二约束仅用于确定 $z$，不影响 $x,y$ 的优化）。这符合直觉：目标函数 $x+2y$ 不含 $z$，约束 $x+y+z=0$ 只确定了 $z$ 而不限制 $x,y$。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 条件极值点处梯度 $\nabla f=0$（将无条件极值的必要条件套用到条件极值） | 条件极值点处 $\nabla f$ 不一定为零，而是与约束梯度的线性组合平行：$\nabla f = \lambda\nabla g$（一个约束）或 $\nabla f = \sum \lambda_i \nabla g_i$（多约束）|
| Lagrange 乘子 $\lambda$ 是多余的变量，与 $x,y$ 地位相同 | $\lambda$ 是**辅助变量**（Lagrange 乘子），引入它是为了将条件极值转化为 $L$ 的无条件驻点问题。其具体数值在求解过程中确定，本身不是优化目标 |
| 代入消元法总是可行的，比 Lagrange 乘数法更简单 | 代入消元法要求约束能显式解出某变量——很多约束无法全局解出（如 $x^2+y^2=1$），此时 Lagrange 乘数法是更系统的工具 |
| 找到候选点后，可以用 Hessian 矩阵（无条件极值判别法）判定条件极值类型 | 条件极值的二阶充分条件是 **bordered Hessian**（加边 Hessian），不能在原 Hessian 上直接套用无条件极值的判别法。通常用实际问题背景或紧集比较函数值来判定 |
| 约束条件一定是等式形式 $g(x,y)=0$ | 本文件仅讨论等式约束。不等式约束（$g(x,y)\leq 0$）属于更一般的优化理论（KKT 条件），不在本文件范围内 |
| 另一个约束一定意味着另一个独立方程——有时两个约束中一个是冗余的 | 两个约束需要 $\nabla g_1$ 和 $\nabla g_2$ 线性无关才构成真正有效的约束。若 $\nabla g_1 \parallel \nabla g_2$，实际只有一个独立约束 |
| 所有约束集合都是紧集，因此最值总是存在 | 只有**有界闭集**（如球面、椭圆）才是紧集。线性约束（如直线或平面）是无界的，最值不一定存在，需单独分析 |

### 检查点

- [ ] 能否给出条件极值的精确定义（定义 12.13），并区分条件极值与无条件极值？
- [ ] 能否用几何方式解释为什么条件极值点处 $\nabla f \parallel \nabla g$（一个约束）？
- [ ] 能否独立推导两个约束情形 $\nabla f = \lambda \nabla g_1 + \mu \nabla g_2$ 的几何意义？
- [ ] 能否写出一个约束下二元函数和三元函数的 Lagrange 方程组？
- [ ] 能否写出两个约束下三元函数的 Lagrange 方程组？
- [ ] 能否列举四种 Lagrange 方程组求解策略，并举例说明每种策略？
- [ ] 能否解释为什么有界闭集（紧集）上的连续函数必存在最值，并判断给定的约束集合是否为紧集？
- [ ] 能否说出 Lagrange 乘数法五步法的每一步具体操作？
- [ ] 能否区分代入消元法和 Lagrange 乘数法的适用条件和优缺点？
- [ ] 能否通过实例演示分类讨论法在解 Lagrange 方程组中的应用？

---

## 练习题

### 基础巩固

**1.** 用代入消元法求 $f(x,y) = x^2 + 2y^2$ 在约束 $2x + y = 3$ 下的条件极值，并判定极值类型。

<details><summary>参考答案</summary>

由 $2x+y=3$ 解出 $y = 3-2x$。代入目标函数：
$$\Phi(x) = f(x,3-2x) = x^2 + 2(3-2x)^2 = x^2 + 2(9 - 12x + 4x^2) = x^2 + 18 - 24x + 8x^2 = 9x^2 - 24x + 18$$

求导：$\Phi'(x) = 18x - 24 = 0 \Rightarrow x = \frac{4}{3}$。对应 $y = 3 - 2\cdot\frac{4}{3} = 3 - \frac{8}{3} = \frac{1}{3}$。

$\Phi''(x) = 18 > 0$，因此在 $\left(\dfrac{4}{3},\dfrac{1}{3}\right)$ 处为条件极小值：
$$f_{\min} = \left(\frac{4}{3}\right)^2 + 2\left(\frac{1}{3}\right)^2 = \frac{16}{9} + \frac{2}{9} = \frac{18}{9} = 2$$

注：约束 $2x+y=3$ 是直线（非紧集），但 $\Phi(x)=9x^2-24x+18$ 是开口向上的二次函数，有全局最小值。

</details>

---

**2.** 用 Lagrange 乘数法求 $f(x,y) = xy$ 在约束 $x + 2y = 4$ 下的条件极值。

<details><summary>参考答案</summary>

构造 Lagrange 函数 $L = xy - \lambda(x+2y-4)$。方程组：
$$\begin{cases}
L_x = y - \lambda = 0 \Rightarrow y = \lambda \\[4pt]
L_y = x - 2\lambda = 0 \Rightarrow x = 2\lambda \\[4pt]
L_\lambda = -(x+2y-4) = 0 \Rightarrow x+2y = 4
\end{cases}$$

代入：$2\lambda + 2\lambda = 4 \Rightarrow 4\lambda = 4 \Rightarrow \lambda = 1$。因此 $x=2$，$y=1$。

候选点 $(2,1)$，条件极值 $f(2,1)=2$。

**极值类型判定**：将 $y=2-x/2$（由约束 $x+2y=4$ 解出）代入 $f$：
$$\Phi(x) = x\left(2-\frac{x}{2}\right) = 2x - \frac{x^2}{2}$$
$\Phi'(x) = 2 - x = 0 \Rightarrow x=2$（一致）。$\Phi''(x) = -1 < 0$，因此为条件极大值 $f_{\max}=2$。

</details>

---

**3.** 求 $f(x,y,z) = x + 2y + 3z$ 在球面 $x^2 + y^2 + z^2 = 14$ 上的最大值和最小值。

<details><summary>参考答案</summary>

$L = (x+2y+3z) - \lambda(x^2+y^2+z^2-14)$

方程组：
$$\begin{cases}
L_x = 1 - 2\lambda x = 0 \Rightarrow x = \dfrac{1}{2\lambda} \\[6pt]
L_y = 2 - 2\lambda y = 0 \Rightarrow y = \dfrac{1}{\lambda} \\[6pt]
L_z = 3 - 2\lambda z = 0 \Rightarrow z = \dfrac{3}{2\lambda} \\[6pt]
L_\lambda = -(x^2+y^2+z^2-14) = 0
\end{cases}$$

代入约束：
$$\left(\frac{1}{2\lambda}\right)^2 + \left(\frac{1}{\lambda}\right)^2 + \left(\frac{3}{2\lambda}\right)^2 = 14$$

计算：
$$\frac{1}{4\lambda^2} + \frac{1}{\lambda^2} + \frac{9}{4\lambda^2} = \frac{1+4+9}{4\lambda^2} = \frac{14}{4\lambda^2} = 14$$

得 $4\lambda^2 = 1$，$\lambda = \pm\dfrac{1}{2}$。

$\lambda = \dfrac12$：$x=1,\ y=2,\ z=3$，$f=1+4+9=14$。
$\lambda = -\dfrac12$：$x=-1,\ y=-2,\ z=-3$，$f=-1-4-9=-14$。

球面 $x^2+y^2+z^2=14$ 是紧集，最值存在。$f_{\max}=14$，$f_{\min}=-14$。

</details>

---

### 迁移应用

**4.** 求函数 $f(x,y,z) = x + z$ 在约束 $x^2 + y^2 = 2$ 和 $y + z = 1$ 下的最大值和最小值。

<details><summary>参考答案</summary>

构造 $L = (x+z) - \lambda(x^2+y^2-2) - \mu(y+z-1)$

方程组：
$$\begin{cases}
L_x = 1 - 2\lambda x = 0 \quad\Rightarrow\quad x = \dfrac{1}{2\lambda} \\[6pt]
L_y = 0 - 2\lambda y - \mu = 0 \quad\Rightarrow\quad \mu = -2\lambda y \\[6pt]
L_z = 1 - 0 - \mu = 0 \quad\Rightarrow\quad \mu = 1 \\[6pt]
L_\lambda = -(x^2+y^2-2)=0 \\[4pt]
L_\mu = -(y+z-1)=0
\end{cases}$$

由 $\mu=1$ 和 $\mu=-2\lambda y$ 得 $-2\lambda y = 1$，即 $y = -\dfrac{1}{2\lambda}$。

代入 $x^2+y^2=2$：
$$\left(\frac{1}{2\lambda}\right)^2 + \left(-\frac{1}{2\lambda}\right)^2 = \frac{1}{4\lambda^2} + \frac{1}{4\lambda^2} = \frac{1}{2\lambda^2} = 2$$

得 $\lambda^2 = \dfrac14$，$\lambda = \pm\dfrac12$。

$\lambda=\dfrac12$：$x=1$，$y=-1$，由 $y+z=1$ 得 $z=2$。$f=1+2=3$。
$\lambda=-\dfrac12$：$x=-1$，$y=1$，$z=0$。$f=-1+0=-1$。

约束 $x^2+y^2=2$（圆柱面）与 $y+z=1$（平面）的交线是空间椭圆（有界闭集），最值存在。

因此 $f_{\max}=3$，$f_{\min}=-1$。

</details>

---

**5.** 用 Lagrange 乘数法推导：在周长固定的三角形中，等边三角形的面积最大。

**提示**：设三角形边长分别为 $a,b,c$，半周长为 $s = \dfrac{a+b+c}{2}$，面积为 $A = \sqrt{s(s-a)(s-b)(s-c)}$（海伦公式）。固定周长 $a+b+c = 2s$（常数），在约束 $a+b+c=2s$ 下求 $A^2 = s(s-a)(s-b)(s-c)$ 的最大值。为简化，求 $\ln A^2$ 的极值（单调变换不改变极值点）。

<details><summary>参考答案</summary>

固定半周长 $s$ 为常数，在约束 $a+b+c=2s$ 下最大化 $A$。等价于最大化 $A^2 = s(s-a)(s-b)(s-c)$，也等价于最大化：
$$F(a,b,c) = \ln(s-a) + \ln(s-b) + \ln(s-c)$$
（因为 $\ln A^2 = \ln s + \sum \ln(s-a_i)$，其中 $\ln s$ 为常数不影响极值点）

约束：$g(a,b,c) = a+b+c-2s = 0$。

构造 $L = [\ln(s-a) + \ln(s-b) + \ln(s-c)] - \lambda(a+b+c-2s)$

方程组：
$$\begin{cases}
L_a = -\dfrac{1}{s-a} - \lambda = 0 \quad\Rightarrow\quad \dfrac{1}{s-a} = -\lambda \\[6pt]
L_b = -\dfrac{1}{s-b} - \lambda = 0 \quad\Rightarrow\quad \dfrac{1}{s-b} = -\lambda \\[6pt]
L_c = -\dfrac{1}{s-c} - \lambda = 0 \quad\Rightarrow\quad \dfrac{1}{s-c} = -\lambda
\end{cases}$$

因此 $s-a = s-b = s-c$，即 $a=b=c$。代入 $a+b+c=2s$ 得 $3a=2s$，$a=2s/3$。

因此 $a=b=c=2s/3$，是等边三角形。由实际意义（所有边长大于 $0$ 的三角形组成有界闭集——实际上是闭的，但需注意退化三角形 $a+b=c$ 的边界情况），等边三角形对应最大面积。

</details>
