# 01. 函数极限的定义（ε-δ定义与Heine归并原则）

> 所属章节：第三章 函数极限与连续函数  |  文件序号：01  |  难度：基础
> 常见混淆点：ε-δ定义中的δ依赖于ε和x₀；0<|x-x₀|表明不要求函数在x₀有定义；Heine定理中要求所有子列都收敛到同一值才能保证极限存在
> 知识点关键词：函数极限，ε-δ定义，Heine定理，归并原则，单侧极限，去心邻域

## 1. 学习目标与先修前置

### 学习目标
- 理解从数列极限 ε-N 到函数极限 ε-δ 的思想迁移（离散到连续），掌握核心逻辑结构的异同
- 掌握函数极限 $\displaystyle\lim_{x\to x_0} f(x) = L$ 的严格 ε-δ 定义，能用逻辑符号完整表述
- 掌握用 ε-δ 定义证明函数极限的两段式方法（分析找 δ + 正序验证）
- 理解极限不存在的 ε-δ 否定表述，能处理无有限极限的情形
- 掌握 Heine 定理（归并原则）的陈述与证明，理解其作为函数极限与数列极限之间桥梁的作用
- 掌握单侧极限（左极限、右极限）的定义及双侧极限存在当且仅当左右极限存在且相等
- 能运用 Heine 定理证明函数极限不存在（构造两个子列收敛到不同值）
- 能分析分段函数在分段点处的极限存在性

### 先修知识
- 文件 08（第二章）：数列极限的 ε-N 定义（定义 8.2）、两段式证明法、放缩法
- 文件 09（第二章）：极限的四则运算法则（定理 9.1）
- 文件 10（第二章）：夹逼准则（迫敛性定理 10.1）、单调有界定理（定理 10.2）
- 文件 02：映射的定义 $f: A \to B$、定义域与陪域
- 读者应熟悉实数的绝对值 $|x|$、开区间 $(a,b)$、邻域的概念
- 读者应熟悉三角不等式 $|x+y| \leq |x| + |y|$

---

## 2. 背景与应用场景

在前一章中，我们详细研究了**数列极限**——当指标 $n$（离散的正整数）趋向无穷时，数列 $a_n$ 的行为。然而，数学和物理世界充满了**连续变化**的量：

- **物理中的瞬时速度**：物体在时刻 $t$ 的瞬时速度是平均速度 $\frac{s(t+\Delta t)-s(t)}{\Delta t}$ 当 $\Delta t \to 0$ 时的极限
- **几何中的切线斜率**：曲线 $y=f(x)$ 在点 $(x_0, f(x_0))$ 处的切线斜率是割线斜率 $\frac{f(x)-f(x_0)}{x-x_0}$ 当 $x \to x_0$ 时的极限
- **化学中的反应速率**：反应物浓度 $C(t)$ 在时刻 $t$ 的瞬时变化率是 $\frac{C(t+\Delta t)-C(t)}{\Delta t}$ 当 $\Delta t \to 0$ 时的极限

所有这些问题的共同核心是：**当自变量 $x$ 趋于某个值 $x_0$ 时，函数值 $f(x)$ 会趋于什么？**

数列极限 $\lim_{n\to\infty} a_n = L$ 刻画的是**离散指标**（$n$ 是正整数）趋近无穷的过程。函数极限 $\lim_{x\to x_0} f(x) = L$ 刻画的是**连续变量**（$x$ 是实数）趋近某个点 $x_0$ 的过程。

两者的根本区别在于：
- $n$ 只能取正整数，且趋于无穷时是"跳跃式"地前进——我们只关心 $n=1,2,3,\dots$ 这些离散点上的值
- $x$ 可以取 $x_0$ 附近的所有实数（除了 $x_0$ 自身），以"连续流动"的方式趋近 $x_0$

尽管有这些区别，极限的精髓是一致的——**"任意接近"的逻辑**。数列极限用 ε-N 定义，函数极限则用 ε-δ 定义。这是从离散到连续的一次自然但深刻的思想迁移。

此外，我们将学习**Heine 定理**（归并原则），它揭示了函数极限与数列极限之间的内在联系：函数极限可以通过数列极限来刻画。这不仅是理论上的桥梁，更是实际计算和证明的利器。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 示例 |
|------|------|------|
| $\lim_{x\to x_0} f(x) = L$ | $x$ 趋于 $x_0$ 时 $f(x)$ 的极限为 $L$ | $\lim_{x\to 2} (3x-1) = 5$ |
| $f(x) \to L\ (x\to x_0)$ | 与极限记号同义 | $x^2 \to 4\ (x\to 2)$ |
| $\varepsilon$ | 任意给定的正数，衡量函数值与极限的接近程度 | $\varepsilon = 0.001$ |
| $\delta$ | 存在的正数，衡量 $x$ 与 $x_0$ 的接近程度（依赖于 $\varepsilon$） | $\delta = 0.0005$ |
| $0<|x-x_0|<\delta$ | $x$ 在 $x_0$ 的去心 $\delta$ 邻域内（$x \neq x_0$） | 开区间 $(x_0-\delta, x_0+\delta)$ 去掉 $x_0$ |
| $\lim_{x\to x_0^+} f(x)$ | 右极限（$x$ 从右侧趋于 $x_0$） | $\lim_{x\to 0^+} \frac{|x|}{x} = 1$ |
| $\lim_{x\to x_0^-} f(x)$ | 左极限（$x$ 从左侧趋于 $x_0$） | $\lim_{x\to 0^-} \frac{|x|}{x} = -1$ |
| $\lim_{x\to\infty} f(x)$ | $x$ 趋于无穷时 $f(x)$ 的极限 | $\lim_{x\to\infty} \frac{1}{x} = 0$ |

### 关键概念的精确定义

#### 3.1 函数极限的 ε-δ 定义

**定义 1.1（$\displaystyle\lim_{x\to x_0} f(x) = L$）**：设函数 $f$ 在点 $x_0$ 的某个去心邻域内有定义（$x_0$ 本身可以没有定义），$L$ 是一个常数。若对任意给定的 $\varepsilon > 0$，都存在 $\delta > 0$（依赖于 $\varepsilon$ 和 $x_0$），使得当 $x$ 满足 $0 < |x - x_0| < \delta$ 时，恒有
$$|f(x) - L| < \varepsilon$$
则称当 $x$ 趋于 $x_0$ 时 $f(x)$ 的极限为 $L$，记作
$$\lim_{x\to x_0} f(x) = L \quad \text{或} \quad f(x) \to L\ (x \to x_0)$$

**用逻辑符号完整表达**：
$$\lim_{x\to x_0} f(x) = L \iff \forall\varepsilon>0,\ \exists\delta>0,\ \forall x\,(0<|x-x_0|<\delta):\ |f(x)-L|<\varepsilon$$

**定义的精读——与 ε-N 定义的逐项类比**：

| 对比项 | ε-N 定义（数列极限） | ε-δ 定义（函数极限） |
|--------|---------------------|---------------------|
| 自变量变化 | $n \to \infty$（离散，正整数） | $x \to x_0$（连续，实数） |
| "足够接近"的条件 | $n > N$（$N$ 是正整数门槛） | $0 < |x-x_0| < \delta$（$\delta > 0$ 是实数半径） |
| 门槛参数 | $N$（依赖于 $\varepsilon$） | $\delta$（依赖于 $\varepsilon$ 和 $x_0$） |
| 接近极限的条件 | $|a_n - L| < \varepsilon$ | $|f(x) - L| < \varepsilon$ |
| 逻辑结构 | $\forall\varepsilon>0,\ \exists N,\ \forall n>N:\ \cdots$ | $\forall\varepsilon>0,\ \exists\delta>0,\ \forall x(0<|x-x_0|<\delta):\ \cdots$ |
| 是否包含"自身" | $n$ 永远不等于 $\infty$（自然） | $x \neq x_0$（显式要求） |

**去心邻域 $0 < |x-x_0| < \delta$ 的几何含义**：
- $|x-x_0| < \delta$：$x$ 落在以 $x_0$ 为中心、$\delta$ 为半径的开区间 $(x_0-\delta, x_0+\delta)$ 内
- $0 < |x-x_0|$：$x \neq x_0$，即排除中心点 $x_0$ 自身
- 因此 $0 < |x-x_0| < \delta$ 表示 $x$ 是区间 $(x_0-\delta, x_0+\delta)$ 中除了中心 $x_0$ 以外的所有点

**为什么要求 $x \neq x_0$？** 极限刻画的是 $x$ "趋近" $x_0$ 时 $f(x)$ 的行为，与 $f$ 在 $x_0$ 处的函数值**无关**。函数可能在 $x_0$ 处无定义（如 $f(x)=\frac{x^2-1}{x-1}$ 在 $x=1$ 处无定义），但仍可能有极限。即使 $f$ 在 $x_0$ 处有定义，极限也可能不等于函数值——这正是"连续"概念将要讨论的内容（后续文件）。

**几何解释**：
$\varepsilon$-$\delta$ 定义有直观的几何解释。在坐标平面上，以 $(x_0, L)$ 为中心，作宽 $2\delta$、高 $2\varepsilon$ 的矩形窗口。$\lim_{x\to x_0} f(x) = L$ 意味着：无论 $\varepsilon$ 多么小，都能找到一个 $\delta$，使得 $f$ 的图像在 $(x_0-\delta, x_0+\delta)$ 上去掉 $x_0$ 点的部分完全落在水平带形区域 $(L-\varepsilon, L+\varepsilon)$ 内。

```
    y
    ↑
 L+ε ────┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈────
    |                      |
 L   ••••••••••••••••••••  |  <-- f(x) 在这段内
    |                      |
 L-ε ────┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈────
    |   |    |    |    |
    +---x₀-δ---x₀---x₀+δ---> x
```

#### 3.2 函数极限不存在的 ε-δ 否定

数列极限不存在的否定表述我们已经熟悉。对于函数极限，$\lim_{x\to x_0} f(x) \neq L$（即极限不存在或存在但不等于 $L$）的否定表述为：

$$\exists\varepsilon_0>0,\ \forall\delta>0,\ \exists x\,(0<|x-x_0|<\delta):\ |f(x)-L|\geq\varepsilon_0$$

即存在某个 $\varepsilon_0 > 0$，无论 $\delta$ 取多么小，都能在 $x_0$ 的去心 $\delta$ 邻域内找到某个 $x$，使得 $|f(x)-L| \geq \varepsilon_0$。

**无穷极限的情形**：函数也可以"发散到无穷"，如 $f(x)=\frac{1}{(x-1)^2}$ 在 $x\to 1$ 时。这种情况下我们说 $\lim_{x\to 1} \frac{1}{(x-1)^2} = +\infty$，定义如下：

$$\forall M>0,\ \exists\delta>0,\ \forall x\,(0<|x-1|<\delta):\ \frac{1}{(x-1)^2} > M$$

即无论 $M$ 取多大的正数，都能找到一个 $\delta$，使得在 $x_0$ 的去心 $\delta$ 邻域内函数值恒大于 $M$。

#### 3.3 单侧极限：左极限与右极限

**定义 1.2（右极限）**：设函数 $f$ 在 $(x_0, x_0+\eta)$（$\eta > 0$）上有定义。若对任意 $\varepsilon > 0$，存在 $\delta > 0$，使得当 $0 < x - x_0 < \delta$ 时，有 $|f(x) - L| < \varepsilon$，则称 $L$ 为 $f$ 在 $x_0$ 处的**右极限**，记作
$$\lim_{x\to x_0^+} f(x) = L \quad \text{或} \quad f(x_0^+) = L$$

**定义 1.3（左极限）**：设函数 $f$ 在 $(x_0-\eta, x_0)$（$\eta > 0$）上有定义。若对任意 $\varepsilon > 0$，存在 $\delta > 0$，使得当 $-\delta < x - x_0 < 0$（即 $x_0-\delta < x < x_0$）时，有 $|f(x) - L| < \varepsilon$，则称 $L$ 为 $f$ 在 $x_0$ 处的**左极限**，记作
$$\lim_{x\to x_0^-} f(x) = L \quad \text{或} \quad f(x_0^-) = L$$

左极限和右极限统称为**单侧极限**。双侧极限（即 $\lim_{x\to x_0} f(x)$）存在当且仅当左右极限存在且相等。我们将在第 4.3 节严格证明这一点。

#### 3.4 Heine 定理（归并原则/归结原则）

**定理 1.1（Heine 定理/归并原则）**：设函数 $f$ 在 $x_0$ 的某个去心邻域内有定义。则 $\displaystyle\lim_{x\to x_0} f(x) = L$ 的充要条件是：对任意数列 $\{x_n\}$，若 $x_n \to x_0$ 且 $x_n \neq x_0$（对所有 $n$），都有 $\displaystyle\lim_{n\to\infty} f(x_n) = L$。

Heine 定理是连接函数极限与数列极限的桥梁。它的核心思想是：**函数极限的存在性可以通过检验所有可能的"趋近路径"（即数列）来判定**。如果沿某条路径趋近时函数值不收敛到 $L$，或者沿两条不同路径趋近时函数值收敛到不同极限，则函数极限不存在。

---

## 4. 原理与方法

### 4.1 用 ε-δ 定义证明函数极限——两段式方法

证明 $\displaystyle\lim_{x\to x_0} f(x) = L$ 的标准方法称为**两段式方法**，与数列极限证明中的两段式方法（文件 08）完全类似。

#### 4.1.1 方法框架

**第一段：分析（找 δ）**
从 $|f(x) - L| < \varepsilon$ 出发，通过代数变形和放缩，将 $|f(x)-L|$ 与 $|x-x_0|$ 联系起来，反推出 $\delta$ 需要满足的条件——即找到一个 $\delta$（作为 $\varepsilon$ 的函数），使得当 $0 < |x-x_0| < \delta$ 时 $|f(x)-L| < \varepsilon$ 自动成立。

**第二段：正序证明（验证）**
按 ε-δ 定义的逻辑顺序写出完整证明：
1. 任取 $\varepsilon > 0$
2. 取 $\delta = \cdots$（从分析段得到的表达式）
3. 设 $0 < |x-x_0| < \delta$，计算 $|f(x)-L|$，证明其 $< \varepsilon$

与数列极限两段式方法的对比：

| 步骤 | ε-N（数列） | ε-δ（函数） |
|------|------------|------------|
| 分析目标 | 从 $|a_n-a|<\varepsilon$ 解出 $n > N$ 的条件 | 从 $|f(x)-L|<\varepsilon$ 解出 $|x-x_0|<\delta$ 的条件 |
| 门槛参数 | $N$（正整数，$\varepsilon$ 的函数） | $\delta$（正实数，$\varepsilon$ 的函数） |
| 验证条件 | $n > N \Rightarrow |a_n-a|<\varepsilon$ | $0<|x-x_0|<\delta \Rightarrow |f(x)-L|<\varepsilon$ |

#### 4.1.2 线性函数的 ε-δ 证明（基础模板）

**例 1**：证明 $\displaystyle\lim_{x\to 2} (5x-3) = 7$。

**分析段**：
$$|(5x-3) - 7| = |5x - 10| = 5|x-2|$$

需要 $5|x-2| < \varepsilon \iff |x-2| < \frac{\varepsilon}{5}$。

因此取 $\delta = \frac{\varepsilon}{5}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \frac{\varepsilon}{5}$。当 $0 < |x-2| < \delta$ 时，
$$|(5x-3)-7| = 5|x-2| < 5 \cdot \frac{\varepsilon}{5} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 2} (5x-3) = 7$。证毕。

**关键点**：对于线性函数 $f(x)=ax+b$，总有 $|f(x)-f(x_0)| = |a| \cdot |x-x_0|$，因此可直接取 $\delta = \varepsilon/|a|$。

#### 4.1.3 二次函数的 ε-δ 证明（含放缩技巧）

**例 2**：证明 $\displaystyle\lim_{x\to 1} x^2 = 1$。

**分析段**：
$$|x^2 - 1| = |(x-1)(x+1)| = |x-1| \cdot |x+1|$$

这里 $|x+1|$ 是随 $x$ 变化的，不能直接取 $\delta = \varepsilon/|x+1|$——因为 $\delta$ 不能依赖于 $x$。我们需要先控制 $x$ 的范围，使 $|x+1|$ 有上界。

**第 1 步：限制 $|x-1| < 1$**，即 $x$ 在区间 $(0, 2)$ 内。此时：
$$|x+1| \leq |x| + 1 < 2 + 1 = 3$$
（严格地说：$x \in (0, 2)$ 时 $x+1 \in (1, 3)$，所以 $|x+1| < 3$）

**第 2 步：推导 $\delta$ 的条件**：
$$|x^2-1| = |x-1| \cdot |x+1| < 3|x-1|$$

需要 $3|x-1| < \varepsilon \iff |x-1| < \frac{\varepsilon}{3}$。

**第 3 步：组合条件**：同时需要 $|x-1| < 1$（保证放缩有效）和 $|x-1| < \frac{\varepsilon}{3}$（保证最终不等式）。因此取：
$$\delta = \min\left\{1,\ \frac{\varepsilon}{3}\right\}$$

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \min\left\{1,\ \frac{\varepsilon}{3}\right\}$。

当 $0 < |x-1| < \delta$ 时，由于 $\delta \leq 1$，有 $|x-1| < 1$，因此 $x \in (0, 2)$，$|x+1| < 3$。又由于 $\delta \leq \frac{\varepsilon}{3}$，有 $|x-1| < \frac{\varepsilon}{3}$。于是：
$$|x^2 - 1| = |x+1| \cdot |x-1| < 3 \cdot \frac{\varepsilon}{3} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 1} x^2 = 1$。证毕。

**关键技巧总结**：当 $|f(x)-L| = |x-x_0| \cdot h(x)$ 且 $h(x)$ 在 $x_0$ 附近有界时，先限制 $|x-x_0| < 1$（或其他常数），得到 $h(x)$ 的上界 $M$，然后取 $\delta = \min\{1,\ \varepsilon/M\}$。

为什么可以这样取？因为 $\delta$ 是"取最小值"：$\delta \leq 1$ 保证了放缩前提成立，$\delta \leq \varepsilon/M$ 保证了最终不等式成立。两个条件同时满足。

#### 4.1.4 含分母的 ε-δ 证明（保证分母有正下界）

**例 3**：证明 $\displaystyle\lim_{x\to 2} \frac{1}{x} = \frac{1}{2}$。

**分析段**：
$$\left|\frac{1}{x} - \frac{1}{2}\right| = \left|\frac{2-x}{2x}\right| = \frac{|x-2|}{2|x|}$$

关键在于处理 $|x|$ 在分母上——需要保证 $x$ 不靠近 $0$，即 $|x|$ 有正下界。限制 $|x-2| < 1$，即 $x \in (1, 3)$。此时 $|x| > 1$，从而 $\frac{1}{|x|} < 1$。

于是：
$$\left|\frac{1}{x} - \frac{1}{2}\right| = \frac{|x-2|}{2|x|} < \frac{|x-2|}{2}$$

需要 $\frac{|x-2|}{2} < \varepsilon \iff |x-2| < 2\varepsilon$。同时要保证 $|x-2| < 1$。取 $\delta = \min\{1, 2\varepsilon\}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \min\{1, 2\varepsilon\}$。

当 $0 < |x-2| < \delta$ 时，$|x-2| < 1$，故 $x \in (1, 3)$，$|x| > 1$，$\frac{1}{|x|} < 1$。又 $|x-2| < 2\varepsilon$，于是：
$$\left|\frac{1}{x} - \frac{1}{2}\right| = \frac{|x-2|}{2|x|} < \frac{|x-2|}{2} < \frac{2\varepsilon}{2} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 2} \frac{1}{x} = \frac{1}{2}$。证毕。

### 4.2 函数极限不存在的判定

#### 4.2.1 用 ε-δ 否定表述证明极限不存在

要证明 $\lim_{x\to x_0} f(x)$ 不存在，或者 $\lim_{x\to x_0} f(x) \neq L$（不收敛到某个特定的 $L$），需要用到否定表述：

$$\exists\varepsilon_0>0,\ \forall\delta>0,\ \exists x\,(0<|x-x_0|<\delta):\ |f(x)-L|\geq\varepsilon_0$$

即存在一个"顽固的" $\varepsilon_0$，无论 $\delta$ 多小，都能在 $x_0$ 附近找到一点使得函数值远离 $L$。

**例 4**：证明 $\displaystyle\lim_{x\to 1} \frac{1}{(x-1)^2}$ 不存在有限极限（即发散到无穷）。

**分析**：直观上，当 $x$ 趋近 $1$ 时，分母 $(x-1)^2$ 趋近 $0$，分数值趋近无穷。我们要证明对任意有限数 $L$，$f(x)$ 都不以 $L$ 为极限。

**证明**：假设存在有限极限 $L$。取 $\varepsilon_0 = 1$。对任意 $\delta > 0$（取 $\delta$ 尽量小，如 $\delta < 1$），取 $x$ 满足 $0 < |x-1| < \delta$。

我们需要证明无论如何都能找到这样的 $x$ 使得 $|f(x)-L| \geq 1$。

由于 $\frac{1}{(x-1)^2}$ 在 $x\to 1$ 时无界增大，对任意 $M > 0$，只要 $|x-1| < \frac{1}{\sqrt{M}}$，就有 $\frac{1}{(x-1)^2} > M$。

取 $M = |L| + 2$，则当 $|x-1| < \min\left\{\delta,\ \frac{1}{\sqrt{|L|+2}}\right\}$ 时：
$$\left|\frac{1}{(x-1)^2} - L\right| \geq \left|\frac{1}{(x-1)^2}\right| - |L| > (|L|+2) - |L| = 2 > 1 = \varepsilon_0$$

因此 $\lim_{x\to 1} \frac{1}{(x-1)^2}$ 不存在有限极限。

实际上，$\lim_{x\to 1} \frac{1}{(x-1)^2} = +\infty$ 在广义极限的意义下成立。

#### 4.2.2 用 Heine 定理证明函数极限不存在

Heine 定理为证明函数极限不存在提供了一个非常实用的方法：**只需构造两个不同的数列 $\{x_n\}$、$\{y_n\}$，都收敛到 $x_0$ 且 $x_n \neq x_0$，$y_n \neq x_0$，但 $\lim f(x_n) \neq \lim f(y_n)$（或其中至少一个极限不存在）**。

**例 5**：证明 $\displaystyle\lim_{x\to 0} \sin\frac{1}{x}$ 不存在。

**构造数列**：
$$x_n = \frac{1}{2n\pi},\quad y_n = \frac{1}{2n\pi + \frac{\pi}{2}}$$

当 $n \to \infty$ 时，$x_n \to 0$，$y_n \to 0$，且 $x_n \neq 0$，$y_n \neq 0$。

计算函数值：
$$f(x_n) = \sin(2n\pi) = 0,\quad \forall n$$
$$f(y_n) = \sin\left(2n\pi + \frac{\pi}{2}\right) = 1,\quad \forall n$$

因此 $\displaystyle\lim_{n\to\infty} f(x_n) = 0$，而 $\displaystyle\lim_{n\to\infty} f(y_n) = 1$。两个数列的极限不同。

由 Heine 定理的逆否命题：存在两个趋于 $0$ 的数列使得 $f(x_n)$ 和 $f(y_n)$ 收敛到不同极限，因此 $\displaystyle\lim_{x\to 0} \sin\frac{1}{x}$ 不存在。

**直观理解**：当 $x \to 0$ 时，$\frac{1}{x}$ 的绝对值趋于无穷大，$\sin\frac{1}{x}$ 在 $[-1,1]$ 之间无限次振荡，振荡频率越来越高，不趋近于任何固定值。

### 4.3 Heine 定理的严格证明

现在我们给出定理 1.1 的完整证明。

**定理 1.1（Heine 定理/归并原则）**：$\displaystyle\lim_{x\to x_0} f(x) = L$ 当且仅当对任意数列 $\{x_n\}$，若 $x_n \to x_0$ 且 $x_n \neq x_0$（对所有 $n$），都有 $\displaystyle\lim_{n\to\infty} f(x_n) = L$。

**证明**：

**必要性（$\Rightarrow$）**：设 $\lim_{x\to x_0} f(x) = L$。任取数列 $\{x_n\}$ 满足 $x_n \to x_0$ 且 $x_n \neq x_0$。

由函数极限定义，对任意 $\varepsilon > 0$，存在 $\delta > 0$，使得当 $0 < |x-x_0| < \delta$ 时 $|f(x)-L| < \varepsilon$。

由 $x_n \to x_0$ 且 $x_n \neq x_0$，对上述 $\delta > 0$，存在 $N \in \mathbb{N}^+$，使得当 $n > N$ 时 $0 < |x_n - x_0| < \delta$。

因此当 $n > N$ 时，$|f(x_n) - L| < \varepsilon$。由数列极限定义，$\lim_{n\to\infty} f(x_n) = L$。必要性得证。

**充分性（$\Leftarrow$）**：我们证明其逆否命题——若 $\lim_{x\to x_0} f(x) \neq L$（即函数极限不存在或存在但不等于 $L$），则存在数列 $\{x_n\}$ 满足 $x_n \to x_0$ 且 $x_n \neq x_0$，使得 $\lim_{n\to\infty} f(x_n) \neq L$。

$\lim_{x\to x_0} f(x) \neq L$ 的 ε-δ 否定表述为：
$$\exists\varepsilon_0>0,\ \forall\delta>0,\ \exists x\,(0<|x-x_0|<\delta):\ |f(x)-L|\geq\varepsilon_0$$

对每个 $n \in \mathbb{N}^+$，取 $\delta = \frac{1}{n}$，则存在 $x_n$ 满足 $0 < |x_n - x_0| < \frac{1}{n}$ 且 $|f(x_n) - L| \geq \varepsilon_0$。

这样构造的 $\{x_n\}$ 满足：
- $x_n \to x_0$（因为 $|x_n - x_0| < 1/n \to 0$）
- $x_n \neq x_0$（因为 $0 < |x_n - x_0|$）
- $|f(x_n) - L| \geq \varepsilon_0$ 对所有 $n$ 成立，因此 $\lim f(x_n) \neq L$（事实上 $\{f(x_n)\}$ 不收敛到 $L$）

因此逆否命题成立，充分性得证。证毕。

**重要推论**：要证明 $\lim_{x\to x_0} f(x)$ 不存在，只需构造两个数列 $\{x_n\}$、$\{y_n\}$，都收敛到 $x_0$ 且不等于 $x_0$，而 $\{f(x_n)\}$ 和 $\{f(y_n)\}$ 收敛到不同的极限。这正是例 5 中使用的方法。

### 4.4 单侧极限与双侧极限的关系

**定理 1.2（极限存在与左右极限的关系）**：$\displaystyle\lim_{x\to x_0} f(x) = L$ 当且仅当 $\displaystyle\lim_{x\to x_0^-} f(x) = \lim_{x\to x_0^+} f(x) = L$。

**证明**：

**（$\Rightarrow$ 必要性）**：设 $\lim_{x\to x_0} f(x) = L$。由定义 1.1，对任意 $\varepsilon > 0$，存在 $\delta > 0$，使得当 $0 < |x-x_0| < \delta$ 时 $|f(x)-L| < \varepsilon$。

当 $0 < x - x_0 < \delta$（右邻域）时，自然满足 $0 < |x-x_0| < \delta$，因此 $|f(x)-L| < \varepsilon$。由定义 1.2，$\lim_{x\to x_0^+} f(x) = L$。

同理，当 $-\delta < x - x_0 < 0$（左邻域）时，也满足 $0 < |x-x_0| < \delta$，因此 $|f(x)-L| < \varepsilon$。由定义 1.3，$\lim_{x\to x_0^-} f(x) = L$。必要性得证。

**（$\Leftarrow$ 充分性）**：设 $\lim_{x\to x_0^-} f(x) = L$ 且 $\lim_{x\to x_0^+} f(x) = L$。

由左极限定义，对任意 $\varepsilon > 0$，存在 $\delta_1 > 0$，使得当 $-\delta_1 < x - x_0 < 0$ 时 $|f(x)-L| < \varepsilon$。

由右极限定义，对任意 $\varepsilon > 0$，存在 $\delta_2 > 0$，使得当 $0 < x - x_0 < \delta_2$ 时 $|f(x)-L| < \varepsilon$。

取 $\delta = \min\{\delta_1, \delta_2\}$。则当 $0 < |x-x_0| < \delta$ 时，$x$ 要么在左邻域（$-\delta < x-x_0 < 0$），要么在右邻域（$0 < x-x_0 < \delta$）。无论哪种情况，都有 $|f(x)-L| < \varepsilon$。因此 $\lim_{x\to x_0} f(x) = L$。证毕。

**定理的意义**：它告诉我们，函数在 $x_0$ 处的极限存在当且仅当左右极限都存在且相等。这个定理对判断分段函数在分段点处的极限是否存在特别有用——只需分别计算左右极限并比较即可。

### 4.5 分段函数在分段点处的极限分析

分段函数在不同区间上有不同的表达式。分析分段点处的极限时，需要分别计算左极限和右极限。

**例 6**：设 $f(x) = \begin{cases} x^2, & x < 1 \\ 2x, & x \geq 1 \end{cases}$。判断 $\displaystyle\lim_{x\to 1} f(x)$ 是否存在。

**左极限**（$x \to 1^-$，即 $x < 1$ 时 $f(x) = x^2$）：
$$\lim_{x\to 1^-} f(x) = \lim_{x\to 1} x^2 = 1^2 = 1$$

**右极限**（$x \to 1^+$，即 $x \geq 1$ 时 $f(x) = 2x$）：
$$\lim_{x\to 1^+} f(x) = \lim_{x\to 1} 2x = 2 \cdot 1 = 2$$

左极限 $= 1$，右极限 $= 2$，两者不相等。由定理 1.2，$\lim_{x\to 1} f(x)$ 不存在。

**数值验证**：
- 从左侧趋近 $1$：$f(0.9)=0.81$，$f(0.99)=0.9801$，$f(0.999)=0.998001$ $\to 1$
- 从右侧趋近 $1$：$f(1.1)=2.2$，$f(1.01)=2.02$，$f(1.001)=2.002$ $\to 2$

两侧趋近的值不同，因此极限不存在。

---

## 5. 例题

### 例题 1：线性函数的 ε-δ 证明

证明 $\displaystyle\lim_{x\to 3} (2x+1) = 7$。

**解**：

**分析段**：
$$|(2x+1) - 7| = |2x - 6| = 2|x-3|$$

需要 $2|x-3| < \varepsilon \iff |x-3| < \frac{\varepsilon}{2}$。

取 $\delta = \frac{\varepsilon}{2}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \frac{\varepsilon}{2}$。当 $0 < |x-3| < \delta$ 时，
$$|(2x+1)-7| = 2|x-3| < 2 \cdot \frac{\varepsilon}{2} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 3} (2x+1) = 7$。证毕。

---

### 例题 2：二次函数的 ε-δ 证明（min 技巧）

证明 $\displaystyle\lim_{x\to 2} x^2 = 4$。

**解**：

**分析段**：
$$|x^2 - 4| = |(x-2)(x+2)| = |x-2| \cdot |x+2|$$

先限制 $|x-2| < 1$，即 $x \in (1, 3)$。此时 $|x+2| < 5$（$x+2 \in (3, 5)$）。

于是：
$$|x^2-4| = |x-2| \cdot |x+2| < 5|x-2|$$

需要 $5|x-2| < \varepsilon \iff |x-2| < \frac{\varepsilon}{5}$。

同时要保证 $|x-2| < 1$，因此取 $\delta = \min\left\{1,\ \frac{\varepsilon}{5}\right\}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \min\left\{1,\ \frac{\varepsilon}{5}\right\}$。

当 $0 < |x-2| < \delta$ 时，由于 $\delta \leq 1$，有 $|x-2| < 1$，因此 $x \in (1, 3)$，$|x+2| < 5$。又由于 $\delta \leq \frac{\varepsilon}{5}$，有 $|x-2| < \frac{\varepsilon}{5}$。于是：
$$|x^2 - 4| = |x+2| \cdot |x-2| < 5 \cdot \frac{\varepsilon}{5} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 2} x^2 = 4$。证毕。

---

### 例题 3：Heine 定理证明函数极限不存在

证明 $\displaystyle\lim_{x\to 0} \cos\frac{1}{x}$ 不存在。

**解**：

取两个数列：
$$x_n = \frac{1}{2n\pi},\quad y_n = \frac{1}{(2n+1)\pi}$$

当 $n \to \infty$ 时，$x_n \to 0$，$y_n \to 0$，且 $x_n \neq 0$，$y_n \neq 0$。

计算函数值：
$$f(x_n) = \cos(2n\pi) = 1,\quad \forall n$$
$$f(y_n) = \cos((2n+1)\pi) = -1,\quad \forall n$$

因此 $\displaystyle\lim_{n\to\infty} f(x_n) = 1$，$\displaystyle\lim_{n\to\infty} f(y_n) = -1$。两个数列的极限不同。

由 Heine 定理（定理 1.1）的逆否命题，$\displaystyle\lim_{x\to 0} \cos\frac{1}{x}$ 不存在。

**注**：与 $\sin\frac{1}{x}$ 类似，$\cos\frac{1}{x}$ 在 $x=0$ 附近也无限振荡，极限不存在。

---

### 例题 4：分段函数在分段点处的极限

设 $f(x) = \begin{cases} 3x+1, & x < 2 \\ 5, & x = 2 \\ x^2, & x > 2 \end{cases}$。判断 $\displaystyle\lim_{x\to 2} f(x)$ 是否存在。

**解**：

**左极限**：$x \to 2^-$ 时 $f(x) = 3x+1$，$\displaystyle\lim_{x\to 2^-} f(x) = 3 \cdot 2 + 1 = 7$。

**右极限**：$x \to 2^+$ 时 $f(x) = x^2$，$\displaystyle\lim_{x\to 2^+} f(x) = 2^2 = 4$。

左极限 $= 7$，右极限 $= 4$，两者不相等。由定理 1.2，$\displaystyle\lim_{x\to 2} f(x)$ 不存在。

**注意**：$f(2) = 5$ 对极限没有影响——即使 $f(2)$ 等于左极限或右极限，只要左右极限不相等，极限就不存在。

---

### 例题 5：用 Heine 定理和单侧极限分析含绝对值函数

分析 $\displaystyle\lim_{x\to 0} \frac{|x|}{x}$ 是否存在。

**解**：

函数 $f(x) = \frac{|x|}{x}$ 在 $x=0$ 处无定义，但在 $x=0$ 的去心邻域内有定义。

当 $x > 0$ 时，$|x| = x$，所以 $\frac{|x|}{x} = 1$。
当 $x < 0$ 时，$|x| = -x$，所以 $\frac{|x|}{x} = -1$。

**方法一（单侧极限）**：

右极限：$\displaystyle\lim_{x\to 0^+} \frac{|x|}{x} = \lim_{x\to 0^+} 1 = 1$

左极限：$\displaystyle\lim_{x\to 0^-} \frac{|x|}{x} = \lim_{x\to 0^-} (-1) = -1$

左右极限不相等，由定理 1.2，$\displaystyle\lim_{x\to 0} \frac{|x|}{x}$ 不存在。

**方法二（Heine 定理）**：

取 $x_n = \frac{1}{n}$（$x_n \to 0^+$），则 $f(x_n) = 1 \to 1$。
取 $y_n = -\frac{1}{n}$（$y_n \to 0^-$），则 $f(y_n) = -1 \to -1$。

两个数列的极限不同，由 Heine 定理，原极限不存在。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| ε-δ 定义要求 $|f(x)-L| < \varepsilon$ 对所有 $x$（包括 $x = x_0$）都成立 | 定义要求 $0 < |x-x_0| < \delta$，即 $x \neq x_0$。极限与 $f(x_0)$ 的值无关——函数在 $x_0$ 处可以无定义，或有定义但 $f(x_0) \neq L$ |
| $\delta$ 必须由 $\varepsilon$ 通过一个显式公式给出，且必须是"最小的"那个 | $\delta$ 只需要存在性，不需要最小。只要找到一个 $\delta$ 使得条件成立即可。不同证法可能给出不同的 $\delta$ |
| 在 ε-δ 证明中，分析段是正式证明的一部分 | 分析段是寻找 $\delta$ 的草稿过程。正式证明只需：任取 $\varepsilon > 0$，取 $\delta = \cdots$，验证当 $0<|x-x_0|<\delta$ 时 $|f(x)-L|<\varepsilon$ |
| $\lim_{x\to x_0} f(x) = L$ 意味着 $f(x)$ 在 $x_0$ 附近单调地趋近 $L$ | 函数趋近极限的方式可以非常复杂——可以在 $L$ 附近振荡趋近（如 $f(x)=x\sin\frac{1}{x}$ 当 $x\to 0$），只要振幅越来越小即可 |
| Heine 定理只是"函数极限等于数列极限"，不需要注意 $x_n \neq x_0$ 的条件 | Heine 定理要求 $x_n \to x_0$ 且 $x_n \neq x_0$。如果允许 $x_n = x_0$，则当 $f$ 在 $x_0$ 处无定义时，$f(x_n)$ 可能无意义；即使有定义，也会引入 $f(x_0)$ 的影响，而函数极限与 $f(x_0)$ 无关 |
| 要求所有数列 $\{x_n\}$（$x_n \to x_0$，$x_n \neq x_0$）都使 $f(x_n) \to L$ 才能断定函数极限存在——这怎么可能验证"所有"数列？ | 充分性的证明用的是反证法（逆否命题），不需要验证所有数列。正面使用 Heine 定理时，通常是用它来证明极限不存在（只需找到两个特定数列即可） |
| 左右极限存在且相等时，函数在该点一定连续 | 左右极限存在且相等只保证极限存在（定理 1.2），但不保证极限等于函数值。连续需要额外要求 $\lim_{x\to x_0} f(x) = f(x_0)$ |
| ε-δ 证明中取 $\delta = \min\{1, \varepsilon/5\}$ 时，$\delta$ 依赖于 $\varepsilon$ 是合理的，但同时也依赖于"1"这个常数——是否可以取其他值？ | 完全可以。取 $\delta = \min\{0.5, \varepsilon/5\}$ 或 $\delta = \min\{2, \varepsilon/5\}$ 都可以，只要确保在 $|x-x_0|<\delta$ 的范围内 $|f(x)-L|$ 有简单的上界即可 |

### 检查点

- [ ] 能否写出 $\displaystyle\lim_{x\to x_0} f(x) = L$ 的 ε-δ 定义（用逻辑符号）？能否写出与 ε-N 定义的核心区别？
- [ ] 能否解释为什么定义中要求 $0 < |x-x_0|$（去心条件）？能否举出函数在 $x_0$ 处无定义但极限存在的例子？
- [ ] 能否完成 $\displaystyle\lim_{x\to 3} (4x-5) = 7$ 的 ε-δ 证明？
- [ ] 能否完成 $\displaystyle\lim_{x\to 1} x^2 = 1$ 的 ε-δ 证明（需要 $\delta = \min\{1, \varepsilon/3\}$ 技巧）？
- [ ] 能否写出 $\displaystyle\lim_{x\to x_0} f(x) \neq L$ 的 ε-δ 否定表述？
- [ ] 能否写出左极限和右极限的精确定义？能否陈述并证明定理 1.2（极限存在 $\iff$ 左右极限存在且相等）？
- [ ] 能否写出 Heine 定理（定理 1.1）的完整陈述？能否给出必要性的证明？
- [ ] 能否用 Heine 定理证明 $\displaystyle\lim_{x\to 0} \sin\frac{1}{x}$ 不存在？
- [ ] 能否判断分段函数在分段点处的极限是否存在？举例说明左右极限不相等的情况。
- [ ] 能否解释 ε-δ 证明中 $\delta = \min\{1, \varepsilon/M\}$ 这种取法的原理？

---

## 练习题

### 基础巩固

**1.** 用 ε-δ 定义证明下列极限：

(1) $\displaystyle\lim_{x\to 1} (3x+2) = 5$

(2) $\displaystyle\lim_{x\to -1} x^2 = 1$

(3) $\displaystyle\lim_{x\to 3} \frac{x^2-9}{x-3} = 6$

<details><summary>参考答案</summary>

**(1)** $\displaystyle\lim_{x\to 1} (3x+2) = 5$

**分析段**：
$$|(3x+2)-5| = |3x-3| = 3|x-1|$$

需要 $3|x-1| < \varepsilon \iff |x-1| < \frac{\varepsilon}{3}$。

取 $\delta = \frac{\varepsilon}{3}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \frac{\varepsilon}{3}$。当 $0 < |x-1| < \delta$ 时，
$$|(3x+2)-5| = 3|x-1| < 3 \cdot \frac{\varepsilon}{3} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 1} (3x+2) = 5$。证毕。

---

**(2)** $\displaystyle\lim_{x\to -1} x^2 = 1$

**分析段**：
$$|x^2-1| = |(x-1)(x+1)| = |x+1| \cdot |x-1|$$

先限制 $|x+1| < 1$，即 $x \in (-2, 0)$。此时 $x-1 \in (-3, -1)$，$|x-1| < 3$。

于是：
$$|x^2-1| = |x+1| \cdot |x-1| < 3|x+1|$$

需要 $3|x+1| < \varepsilon \iff |x+1| < \frac{\varepsilon}{3}$。

取 $\delta = \min\left\{1,\ \frac{\varepsilon}{3}\right\}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \min\left\{1,\ \frac{\varepsilon}{3}\right\}$。

当 $0 < |x+1| < \delta$ 时，$|x+1| < 1$，所以 $x \in (-2, 0)$，$|x-1| < 3$。又 $|x+1| < \frac{\varepsilon}{3}$，于是：
$$|x^2-1| = |x+1| \cdot |x-1| < 3 \cdot \frac{\varepsilon}{3} = \varepsilon$$

因此 $\displaystyle\lim_{x\to -1} x^2 = 1$。证毕。

---

**(3)** $\displaystyle\lim_{x\to 3} \frac{x^2-9}{x-3} = 6$

**分析段**：当 $x \neq 3$ 时，
$$\frac{x^2-9}{x-3} = \frac{(x-3)(x+3)}{x-3} = x+3$$

因此 $\left|\frac{x^2-9}{x-3} - 6\right| = |(x+3)-6| = |x-3|$。

需要 $|x-3| < \varepsilon$。取 $\delta = \varepsilon$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \varepsilon$。当 $0 < |x-3| < \delta$ 时（$x \neq 3$，化简有效）：
$$\left|\frac{x^2-9}{x-3} - 6\right| = |(x+3)-6| = |x-3| < \varepsilon$$

因此 $\displaystyle\lim_{x\to 3} \frac{x^2-9}{x-3} = 6$。证毕。

</details>

---

**2.** 判断下列分段函数在指定点的极限是否存在：

(1) $f(x) = \begin{cases} x+1, & x < 0 \\ 2x, & x \geq 0 \end{cases}$，$x_0 = 0$

(2) $g(x) = \begin{cases} x^2, & x < 2 \\ 4, & x = 2 \\ 3x-2, & x > 2 \end{cases}$，$x_0 = 2$

<details><summary>参考答案</summary>

**(1)** $f(x) = \begin{cases} x+1, & x < 0 \\ 2x, & x \geq 0 \end{cases}$，$x_0 = 0$

**左极限**：$x \to 0^-$ 时 $f(x) = x+1 \to 0+1 = 1$。

**右极限**：$x \to 0^+$ 时 $f(x) = 2x \to 2 \cdot 0 = 0$。

左极限 $= 1$，右极限 $= 0$，不相等。由定理 1.2，$\lim_{x\to 0} f(x)$ 不存在。

---

**(2)** $g(x) = \begin{cases} x^2, & x < 2 \\ 4, & x = 2 \\ 3x-2, & x > 2 \end{cases}$，$x_0 = 2$

**左极限**：$x \to 2^-$ 时 $g(x) = x^2 \to 2^2 = 4$。

**右极限**：$x \to 2^+$ 时 $g(x) = 3x-2 \to 3\cdot 2 - 2 = 4$。

左极限 $= 4$，右极限 $= 4$，相等。由定理 1.2，$\lim_{x\to 2} g(x) = 4$。

注意 $g(2) = 4$ 也等于极限值，因此函数在 $x=2$ 处连续。但如果将条件改为 $g(2) = 5$，极限仍为 $4$——极限与函数值无关。

</details>

---

### 迁移应用

**3.** 用 Heine 定理证明 $\displaystyle\lim_{x\to 0} \sin\frac{1}{x}$ 不存在（用不同的数列构造）。

<details><summary>参考答案</summary>

**证明**：

取两个数列：
$$x_n = \frac{1}{n\pi},\quad y_n = \frac{2}{(2n+1)\pi}$$

当 $n \to \infty$ 时，$x_n \to 0$，$y_n \to 0$，且 $x_n \neq 0$，$y_n \neq 0$。

计算函数值：
$$f(x_n) = \sin(n\pi) = 0,\quad \forall n$$
$$f(y_n) = \sin\left(\frac{(2n+1)\pi}{2}\right) = \begin{cases} 1, & n \text{ 为偶数} \\ -1, & n \text{ 为奇数} \end{cases}$$

即 $\{f(y_n)\}$ 在 $1$ 和 $-1$ 之间振荡，$\lim_{n\to\infty} f(y_n)$ 不存在（不收敛到任何值）。

由 Heine 定理：存在趋于 $0$ 的数列 $\{y_n\}$ 使得 $\{f(y_n)\}$ 发散，因此 $\displaystyle\lim_{x\to 0} \sin\frac{1}{x}$ 不存在。

**注**：与例题 3 中取 $f(y_n) \equiv 1$ 不同，这里我们构造了使 $f(y_n)$ 本身就不收敛的数列，同样可以说明问题。

</details>

---

**4.** 用 ε-δ 定义证明 $\displaystyle\lim_{x\to 1} (x^2 + 2x) = 3$。

<details><summary>参考答案</summary>

**分析段**：
$$|(x^2+2x)-3| = |x^2+2x-3| = |(x-1)(x+3)| = |x-1| \cdot |x+3|$$

先限制 $|x-1| < 1$，即 $x \in (0, 2)$。此时 $x+3 \in (3, 5)$，$|x+3| < 5$。

于是：
$$|x^2+2x-3| = |x-1| \cdot |x+3| < 5|x-1|$$

需要 $5|x-1| < \varepsilon \iff |x-1| < \frac{\varepsilon}{5}$。

取 $\delta = \min\left\{1,\ \frac{\varepsilon}{5}\right\}$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \min\left\{1,\ \frac{\varepsilon}{5}\right\}$。

当 $0 < |x-1| < \delta$ 时，$|x-1| < 1$，所以 $x \in (0, 2)$，$|x+3| < 5$。又 $|x-1| < \frac{\varepsilon}{5}$，于是：
$$|(x^2+2x)-3| = |x-1| \cdot |x+3| < 5 \cdot \frac{\varepsilon}{5} = \varepsilon$$

因此 $\displaystyle\lim_{x\to 1} (x^2+2x) = 3$。证毕。

</details>

---

**5.** （综合）设 $f(x) = \begin{cases} x\sin\frac{1}{x}, & x \neq 0 \\ 0, & x = 0 \end{cases}$。证明 $\displaystyle\lim_{x\to 0} f(x) = 0$。

<details><summary>参考答案</summary>

**分析段**：当 $x \neq 0$ 时，
$$|f(x)-0| = \left|x\sin\frac{1}{x}\right| = |x| \cdot \left|\sin\frac{1}{x}\right| \leq |x| \cdot 1 = |x|$$

（因为 $|\sin\theta| \leq 1$ 对任意 $\theta$ 成立）

因此 $|f(x)-0| \leq |x|$。

需要 $|x| < \varepsilon$。取 $\delta = \varepsilon$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $\delta = \varepsilon$。当 $0 < |x-0| < \delta$ 时，
$$|f(x)-0| = \left|x\sin\frac{1}{x}\right| \leq |x| < \varepsilon$$

因此 $\displaystyle\lim_{x\to 0} f(x) = 0$。证毕。

**关键观察**：这里 $\left|\sin\frac{1}{x}\right| \leq 1$ 起到了关键作用——我们不需要 $\sin\frac{1}{x}$ 收敛，只需要它有界。这正是"有界乘无穷小"性质（命题 9.1(3)）的体现：$x \to 0$ 时 $x$ 是无穷小量，$\sin\frac{1}{x}$ 有界，乘积为无穷小量。

注意这与 $\lim_{x\to 0} \sin\frac{1}{x}$ 不存在并不矛盾——多乘了一个趋于 $0$ 的因子 $x$ 后，振荡被"压制"了。

</details>
