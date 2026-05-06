# 01. Euclid 空间的基本概念——向量、内积、范数与距离

> 所属章节：第十一章 Euclid 空间上的极限和连续  |  文件序号：01  |  难度：基础
> 常见混淆点：Cauchy-Schwarz 不等式的证明中，构造二次函数 $f(\lambda) = \|\mathbf{x} - \lambda\mathbf{y}\|^2$ 的动机——为什么引入参数 $\lambda$？核心目的是将内积的绝对值上界问题转化为判别式条件；$\mathbb{R}^n$ 中点列收敛的 $\varepsilon$-$N$ 定义与一维 $\varepsilon$-$N$ 定义的结构完全一致，仅将绝对值 $|\cdot|$ 替换为范数 $\|\cdot\|$，但初学者常忽略这一替换的深层含义——$\|\cdot\|$ 已经不再是 $\mathbb{R}$ 上的绝对值，而是 $\mathbb{R}^n$ 上的度量

## 1. 学习目标与先修前置

### 学习目标
- 理解 $\mathbb{R}^n$ 作为向量空间的基本结构：向量加法、标量乘法、标准基与坐标表示
- 掌握 Euclid 内积 $\langle \mathbf{x}, \mathbf{y} \rangle$ 的定义及其三条基本性质（对称性、双线性、正定性）
- 掌握 Euclid 范数 $\|\mathbf{x}\|$ 的定义及其三条基本性质（非负性、齐次性、三角不等式）
- 理解 Euclid 距离 $d(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\|$ 的定义
- 掌握 Cauchy-Schwarz 不等式的判别式证明方法
- 理解并证明向量版的三角不等式（利用 Cauchy-Schwarz 不等式）
- 掌握 $\mathbb{R}^n$ 中点列收敛的 $\varepsilon$-$N$ 定义，并与一维数列极限的定义建立类比
- 理解并证明分量收敛与整体收敛的等价关系

### 先修知识
- 文件 08（第二章）：数列极限的 $\varepsilon$-$N$ 定义、绝对值 $|\cdot|$ 的距离含义、三角不等式 $|x+y| \leq |x| + |y|$
- 读者应熟悉实数集 $\mathbb{R}$ 的基本运算和不等式性质
- 读者应了解 $\forall$（任意）和 $\exists$（存在）逻辑量词
- 读者应了解二次函数 $f(\lambda) = a\lambda^2 + b\lambda + c$ 的判别式 $\Delta = b^2 - 4ac$ 与非负性条件的关系

---

## 2. 背景与应用场景

在前十章中，我们研究的对象始终是**一元**函数和数列：数列 $\{a_n\}$ 的每一项是实数，函数 $f(x)$ 的自变量和因变量都是实数。距离由绝对值 $|x-y|$ 度量。

然而，现实世界中的问题很少只涉及一个变量。例如：
- 一个质点在三维空间中的位置需要三个坐标 $(x, y, z)$ 确定
- 气象观测站记录的温度、气压、湿度构成一个三维数据点
- 经济学中，一种商品的价格由供需、成本、政策等多个因素共同影响

这些场景中，"数据点"不再是单个实数，而是 $n$ 个实数构成的有序组。自然地，我们需要将分析学的基本工具——距离、极限、连续性——从一维推广到 $n$ 维。

本章的目标是建立 $\mathbb{R}^n$ 上的分析学基础。第一步（即本文件）是定义 $\mathbb{R}^n$ 的基本结构：向量、内积、范数和距离，并建立最基本的分析工具——$\mathbb{R}^n$ 中点列收敛的定义及其与分量收敛的等价关系。

**从一维到 $n$ 维的关键类比**：
- 一维中的数 $x \in \mathbb{R}$ → $n$ 维中的向量 $\mathbf{x} = (x_1, x_2, \dots, x_n) \in \mathbb{R}^n$
- 一维中的绝对值 $|x-y|$ → $n$ 维中的范数 $\|\mathbf{x} - \mathbf{y}\|$
- 一维中的数列 $\{a_n\}$ → $n$ 维中的点列 $\{P_k\}$
- 一维中的 $\lim_{n\to\infty} a_n = a$ → $n$ 维中的 $\lim_{k\to\infty} P_k = P$

这些类比不仅仅是形式的相似——在下一章学习多元函数的极限时，这些概念将发挥根本性的作用。

---

## 3. 核心概念与符号约定

### 符号表

| 符号 | 含义 | 示例 |
|------|------|------|
| $\mathbb{R}^n$ | $n$ 维实向量空间 | $\mathbb{R}^2 = \{(x,y) \mid x,y \in \mathbb{R}\}$ |
| $\mathbf{x} = (x_1, \dots, x_n)$ | $n$ 维向量（通常用粗体表示） | $\mathbf{x} = (1, -2, 3)$ |
| $\mathbf{e}_i$ | 第 $i$ 个标准基向量 | $\mathbf{e}_1 = (1,0,\dots,0)$ |
| $\langle \mathbf{x}, \mathbf{y} \rangle$ | Euclid 内积 | $\langle (1,2), (3,4) \rangle = 11$ |
| $\|\mathbf{x}\|$ | Euclid 范数（长度） | $\|(3,4)\| = 5$ |
| $d(\mathbf{x}, \mathbf{y})$ | Euclid 距离 | $d((1,2),(4,6)) = 5$ |
| $\{P_k\}_{k=1}^\infty$ | $\mathbb{R}^n$ 中的点列 | $P_k = (1/k, 1-1/k)$ |
| $P_k \to P$ | 点列 $P_k$ 收敛到 $P$ | $(1/k, 1-1/k) \to (0,1)$ |
| $\mathbb{N}^+$ | 正整数集 $\{1,2,3,\dots\}$ | |
| $\forall$ | 任意（对所有的） | $\forall \varepsilon > 0$ |
| $\exists$ | 存在 | $\exists N \in \mathbb{N}^+$ |

### 3.1 $\mathbb{R}^n$ 作为向量空间

**定义 11.1（$\mathbb{R}^n$ 与向量空间结构）**：设 $n$ 是正整数。记
$$\mathbb{R}^n = \{(x_1, x_2, \dots, x_n) \mid x_i \in \mathbb{R},\ i = 1,2,\dots,n\}$$
即所有 $n$ 元实数组构成的集合。$\mathbb{R}^n$ 中的元素称为**向量**（vector），通常用粗体字母表示：$\mathbf{x} = (x_1, x_2, \dots, x_n)$，其中 $x_i$ 称为 $\mathbf{x}$ 的**第 $i$ 个分量**（component）。

在 $\mathbb{R}^n$ 上定义两种运算：

**向量加法**：对任意 $\mathbf{x} = (x_1,\dots,x_n)$，$\mathbf{y} = (y_1,\dots,y_n) \in \mathbb{R}^n$，
$$\mathbf{x} + \mathbf{y} = (x_1 + y_1,\ x_2 + y_2,\ \dots,\ x_n + y_n)$$
即对应分量分别相加。

**标量乘法**：对任意 $\mathbf{x} = (x_1,\dots,x_n) \in \mathbb{R}^n$ 和任意实数 $c \in \mathbb{R}$，
$$c\mathbf{x} = (c x_1,\ c x_2,\ \dots,\ c x_n)$$
即每个分量乘以 $c$。

在定义了这两种运算后，$\mathbb{R}^n$ 构成一个**向量空间**（vector space），满足以下八条基本性质（验证从略，读者可自行验证）：
1. $\mathbf{x} + \mathbf{y} = \mathbf{y} + \mathbf{x}$（加法交换律）
2. $(\mathbf{x} + \mathbf{y}) + \mathbf{z} = \mathbf{x} + (\mathbf{y} + \mathbf{z})$（加法结合律）
3. 存在零向量 $\mathbf{0} = (0,0,\dots,0)$，使得 $\mathbf{x} + \mathbf{0} = \mathbf{x}$
4. 对每个 $\mathbf{x}$，存在负向量 $-\mathbf{x} = (-x_1,\dots,-x_n)$，使得 $\mathbf{x} + (-\mathbf{x}) = \mathbf{0}$
5. $1\mathbf{x} = \mathbf{x}$
6. $c(d\mathbf{x}) = (cd)\mathbf{x}$（标量乘法结合律）
7. $(c + d)\mathbf{x} = c\mathbf{x} + d\mathbf{x}$（分配律 I）
8. $c(\mathbf{x} + \mathbf{y}) = c\mathbf{x} + c\mathbf{y}$（分配律 II）

**例 1**：在 $\mathbb{R}^3$ 中，计算 $\mathbf{x} + \mathbf{y}$ 和 $3\mathbf{x}$，其中 $\mathbf{x} = (1,2,-1)$，$\mathbf{y} = (2,-1,3)$。

$$\mathbf{x} + \mathbf{y} = (1+2,\ 2+(-1),\ -1+3) = (3, 1, 2)$$
$$3\mathbf{x} = (3\times 1,\ 3\times 2,\ 3\times (-1)) = (3, 6, -3)$$

### 3.2 标准基与坐标表示

**定义 11.2（标准基）**：$\mathbb{R}^n$ 中的**标准基**（standard basis）由以下 $n$ 个向量构成：对 $i = 1,2,\dots,n$，
$$\mathbf{e}_i = (\underbrace{0,\dots,0}_{i-1\text{ 个}},\ 1,\ 0,\dots,0)$$
即第 $i$ 个分量为 $1$，其余分量为 $0$。

例如，在 $\mathbb{R}^3$ 中：
$$\mathbf{e}_1 = (1,0,0),\quad \mathbf{e}_2 = (0,1,0),\quad \mathbf{e}_3 = (0,0,1)$$

**向量的坐标表示**：任何向量 $\mathbf{x} = (x_1, x_2, \dots, x_n) \in \mathbb{R}^n$ 可以唯一地表示为标准基的线性组合：
$$\mathbf{x} = x_1\mathbf{e}_1 + x_2\mathbf{e}_2 + \cdots + x_n\mathbf{e}_n = \sum_{i=1}^n x_i\mathbf{e}_i$$

这给出了向量的几何意义：分量 $x_i$ 就是向量 $\mathbf{x}$ 在第 $i$ 个坐标轴方向上的"坐标"。

**例 2**：在 $\mathbb{R}^3$ 中，向量 $\mathbf{x} = (1,2,-1)$ 可以写为：
$$\mathbf{x} = 1\cdot\mathbf{e}_1 + 2\cdot\mathbf{e}_2 + (-1)\cdot\mathbf{e}_3 = \mathbf{e}_1 + 2\mathbf{e}_2 - \mathbf{e}_3$$

### 3.3 Euclid 内积

**定义 11.3（Euclid 内积）**：对 $\mathbf{x} = (x_1,\dots,x_n)$，$\mathbf{y} = (y_1,\dots,y_n) \in \mathbb{R}^n$，定义它们之间的**Euclid 内积**（Euclidean inner product）为
$$\langle \mathbf{x}, \mathbf{y} \rangle = \sum_{i=1}^n x_i y_i = x_1 y_1 + x_2 y_2 + \cdots + x_n y_n$$

**性质**：Euclid 内积满足以下三条基本性质（对任意 $\mathbf{x},\mathbf{y},\mathbf{z} \in \mathbb{R}^n$ 和 $c \in \mathbb{R}$）：

1. **对称性**：$\langle \mathbf{x}, \mathbf{y} \rangle = \langle \mathbf{y}, \mathbf{x} \rangle$
   - 证明：$\langle \mathbf{x}, \mathbf{y} \rangle = \sum x_i y_i = \sum y_i x_i = \langle \mathbf{y}, \mathbf{x} \rangle$，其中等号成立是因为实数乘法交换。

2. **双线性**：对两个变量分别线性——
   $$\langle \mathbf{x} + \mathbf{y}, \mathbf{z} \rangle = \langle \mathbf{x}, \mathbf{z} \rangle + \langle \mathbf{y}, \mathbf{z} \rangle$$
   $$\langle c\mathbf{x}, \mathbf{y} \rangle = c\langle \mathbf{x}, \mathbf{y} \rangle$$
   $$\langle \mathbf{x}, \mathbf{y} + \mathbf{z} \rangle = \langle \mathbf{x}, \mathbf{y} \rangle + \langle \mathbf{x}, \mathbf{z} \rangle$$
   $$\langle \mathbf{x}, c\mathbf{y} \rangle = c\langle \mathbf{x}, \mathbf{y} \rangle$$
   - 证明（仅验证第一式，其余类似）：
     $$\langle \mathbf{x}+\mathbf{y}, \mathbf{z} \rangle = \sum_{i=1}^n (x_i + y_i)z_i = \sum_{i=1}^n (x_i z_i + y_i z_i) = \sum_{i=1}^n x_i z_i + \sum_{i=1}^n y_i z_i = \langle \mathbf{x}, \mathbf{z} \rangle + \langle \mathbf{y}, \mathbf{z} \rangle$$

3. **正定性**：$\langle \mathbf{x}, \mathbf{x} \rangle \geq 0$，且等号成立当且仅当 $\mathbf{x} = \mathbf{0}$。
   - 证明：$\langle \mathbf{x}, \mathbf{x} \rangle = \sum_{i=1}^n x_i^2$ 是 $n$ 个非负实数的和，因此 $\geq 0$。若 $\mathbf{x} = \mathbf{0} = (0,\dots,0)$，则 $\sum x_i^2 = 0$；反之，若 $\sum x_i^2 = 0$，则由 $x_i^2 \geq 0$ 知每个 $x_i^2 = 0$，从而每个 $x_i = 0$，即 $\mathbf{x} = \mathbf{0}$。

**例 3**：在 $\mathbb{R}^3$ 中计算 $\langle \mathbf{x}, \mathbf{y} \rangle$，其中 $\mathbf{x} = (1,2,-1)$，$\mathbf{y} = (2,-1,3)$。

$$\langle \mathbf{x}, \mathbf{y} \rangle = 1 \times 2 + 2 \times (-1) + (-1) \times 3 = 2 - 2 - 3 = -3$$

### 3.4 Euclid 范数

**定义 11.4（Euclid 范数）**：对 $\mathbf{x} = (x_1,\dots,x_n) \in \mathbb{R}^n$，定义 $\mathbf{x}$ 的**Euclid 范数**（Euclidean norm）为
$$\|\mathbf{x}\| = \sqrt{\langle \mathbf{x}, \mathbf{x} \rangle} = \sqrt{\sum_{i=1}^n x_i^2}$$

**说明**：范数 $\|\mathbf{x}\|$ 度量向量 $\mathbf{x}$ 的"长度"，是 $\mathbb{R}$ 上绝对值 $|\cdot|$ 在 $\mathbb{R}^n$ 中的自然推广：
- 当 $n = 1$ 时，$\|x\| = \sqrt{x^2} = |x|$，范数回归为绝对值。

**性质**：Euclid 范数满足以下三条基本性质（对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ 和 $c \in \mathbb{R}$）：

1. **非负性**：$\|\mathbf{x}\| \geq 0$，且 $\|\mathbf{x}\| = 0$ 当且仅当 $\mathbf{x} = \mathbf{0}$。
   - 证明：由内积的正定性，$\|\mathbf{x}\|^2 = \langle \mathbf{x}, \mathbf{x} \rangle \geq 0$，故 $\|\mathbf{x}\| \geq 0$。且 $\|\mathbf{x}\| = 0 \iff \|\mathbf{x}\|^2 = 0 \iff \langle \mathbf{x}, \mathbf{x} \rangle = 0 \iff \mathbf{x} = \mathbf{0}$。

2. **齐次性（正齐次性）**：$\|c\mathbf{x}\| = |c| \cdot \|\mathbf{x}\|$。
   - 证明：$$\|c\mathbf{x}\|^2 = \langle c\mathbf{x}, c\mathbf{x} \rangle = c^2 \langle \mathbf{x}, \mathbf{x} \rangle = c^2 \|\mathbf{x}\|^2$$ 两边开平方得 $\|c\mathbf{x}\| = |c| \cdot \|\mathbf{x}\|$（注意 $|c|$ 是 $c$ 的绝对值，因为范数始终非负）。

3. **三角不等式**（将在第 4.2 节证明）：$\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$。

**例 4**：计算 $\mathbf{x} = (1,2,-1)$ 的 Euclid 范数。
$$\|\mathbf{x}\| = \sqrt{1^2 + 2^2 + (-1)^2} = \sqrt{1 + 4 + 1} = \sqrt{6}$$

### 3.5 Euclid 距离

**定义 11.5（Euclid 距离）**：对 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，定义它们之间的**Euclid 距离**（Euclidean distance）为
$$d(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\| = \sqrt{\sum_{i=1}^n (x_i - y_i)^2}$$

**说明**：距离 $d(\mathbf{x}, \mathbf{y})$ 度量两个向量之间的"远近"，是一维绝对值 $|x-y|$ 在 $\mathbb{R}^n$ 中的推广。当 $n=1$ 时，$d(x,y) = |x-y|$。

Euclid 距离满足以下基本性质（由范数的性质直接推出）：
1. **非负性**：$d(\mathbf{x}, \mathbf{y}) \geq 0$，且 $d(\mathbf{x}, \mathbf{y}) = 0$ 当且仅当 $\mathbf{x} = \mathbf{y}$
2. **对称性**：$d(\mathbf{x}, \mathbf{y}) = d(\mathbf{y}, \mathbf{x})$
3. **三角不等式**：$d(\mathbf{x}, \mathbf{z}) \leq d(\mathbf{x}, \mathbf{y}) + d(\mathbf{y}, \mathbf{z})$

**例 5**：计算 $\mathbf{x} = (1,2,-1)$ 与 $\mathbf{y} = (2,-1,3)$ 之间的 Euclid 距离。

首先计算差值向量：
$$\mathbf{x} - \mathbf{y} = (1-2,\ 2-(-1),\ -1-3) = (-1, 3, -4)$$
然后计算距离：
$$d(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\| = \sqrt{(-1)^2 + 3^2 + (-4)^2} = \sqrt{1 + 9 + 16} = \sqrt{26}$$

**验证**：$\sqrt{26} \approx 5.099$。作为对照，一维中 $\mathbf{x}$ 和 $\mathbf{y}$ 各分量的差的绝对值分别为 $1$、$3$、$4$，它们的平方和为 $26$，开平方后约 $5.099$。

---

## 4. 原理与方法

### 4.1 Cauchy-Schwarz 不等式

Cauchy-Schwarz 不等式是分析学中最基本、最重要的不等式之一。它将内积的绝对值与范数的乘积联系起来。

**定理 11.1（Cauchy-Schwarz 不等式）**：对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，
$$|\langle \mathbf{x}, \mathbf{y} \rangle| \leq \|\mathbf{x}\| \cdot \|\mathbf{y}\|$$
等号成立当且仅当 $\mathbf{x}$ 与 $\mathbf{y}$ 线性相关（即存在 $\lambda \in \mathbb{R}$ 使得 $\mathbf{x} = \lambda \mathbf{y}$ 或 $\mathbf{y} = \lambda \mathbf{x}$）。

**证明（判别式法）**：

**第一步：构造辅助函数**

对于任意实数 $\lambda$，考虑向量 $\mathbf{x} - \lambda\mathbf{y}$。由范数的非负性，有
$$\|\mathbf{x} - \lambda\mathbf{y}\|^2 \geq 0$$

展开这个平方（利用内积的双线性）：
$$
\begin{aligned}
\|\mathbf{x} - \lambda\mathbf{y}\|^2 &= \langle \mathbf{x} - \lambda\mathbf{y},\ \mathbf{x} - \lambda\mathbf{y} \rangle \\
&= \langle \mathbf{x}, \mathbf{x} \rangle - \lambda\langle \mathbf{x}, \mathbf{y} \rangle - \lambda\langle \mathbf{y}, \mathbf{x} \rangle + \lambda^2\langle \mathbf{y}, \mathbf{y} \rangle \\
&= \|\mathbf{x}\|^2 - 2\lambda\langle \mathbf{x}, \mathbf{y} \rangle + \lambda^2\|\mathbf{y}\|^2
\end{aligned}
$$

其中用到了 $\langle \mathbf{x}, \mathbf{y} \rangle = \langle \mathbf{y}, \mathbf{x} \rangle$（对称性）。

因此，对任意 $\lambda \in \mathbb{R}$，有
$$\|\mathbf{y}\|^2 \lambda^2 - 2\langle \mathbf{x}, \mathbf{y} \rangle \lambda + \|\mathbf{x}\|^2 \geq 0$$

**第二步：应用判别式**

上式是 $\lambda$ 的二次函数 $f(\lambda) = a\lambda^2 + b\lambda + c$，其中：
$$a = \|\mathbf{y}\|^2,\quad b = -2\langle \mathbf{x}, \mathbf{y} \rangle,\quad c = \|\mathbf{x}\|^2$$

$f(\lambda) \geq 0$ 对任意 $\lambda \in \mathbb{R}$ 成立。这意味着该二次函数要么没有实根（开口向上且图像在 $x$ 轴上方），要么只有一个重根（图像与 $x$ 轴相切）。用判别式表示，即：
$$\Delta = b^2 - 4ac \leq 0$$

代入 $a, b, c$ 的具体表达式：
$$
\begin{aligned}
\Delta &= (-2\langle \mathbf{x}, \mathbf{y} \rangle)^2 - 4 \cdot \|\mathbf{y}\|^2 \cdot \|\mathbf{x}\|^2 \\
&= 4\langle \mathbf{x}, \mathbf{y} \rangle^2 - 4\|\mathbf{x}\|^2\|\mathbf{y}\|^2 \leq 0
\end{aligned}
$$

两边同时除以 $4$，得：
$$\langle \mathbf{x}, \mathbf{y} \rangle^2 \leq \|\mathbf{x}\|^2 \|\mathbf{y}\|^2$$

由于 $\langle \mathbf{x}, \mathbf{y} \rangle^2$、$\|\mathbf{x}\|^2$、$\|\mathbf{y}\|^2$ 均为非负数，两边取平方根（非负数的平方根单调递增），得：
$$|\langle \mathbf{x}, \mathbf{y} \rangle| \leq \|\mathbf{x}\| \cdot \|\mathbf{y}\|$$

**第三步：等号成立条件**

等号成立 $\iff \Delta = 0 \iff$ 二次函数有重根 $\iff$ 存在 $\lambda_0 \in \mathbb{R}$ 使得 $f(\lambda_0) = 0$。

$f(\lambda_0) = 0$ 意味着 $\|\mathbf{x} - \lambda_0\mathbf{y}\|^2 = 0$，即 $\|\mathbf{x} - \lambda_0\mathbf{y}\| = 0$。由范数的非负性（定义 11.4 性质 1），
$$\|\mathbf{x} - \lambda_0\mathbf{y}\| = 0 \iff \mathbf{x} - \lambda_0\mathbf{y} = \mathbf{0} \iff \mathbf{x} = \lambda_0\mathbf{y}$$

因此等号成立当且仅当 $\mathbf{x}$ 与 $\mathbf{y}$ 线性相关（$\mathbf{x} = \lambda_0\mathbf{y}$，或者对称地 $\mathbf{y} = \lambda_1\mathbf{x}$，只要一个向量是另一个的标量倍）。

**例 6**：验证 Cauchy-Schwarz 不等式对 $\mathbf{x} = (1,2,-1)$ 和 $\mathbf{y} = (2,-1,3)$ 成立。

计算内积绝对值：$|\langle \mathbf{x}, \mathbf{y} \rangle| = |{-3}| = 3$
计算范数乘积：$\|\mathbf{x}\| \cdot \|\mathbf{y}\| = \sqrt{6} \times \sqrt{14} = \sqrt{84} = 2\sqrt{21} \approx 9.165$

验证：$3 \leq 9.165$，不等式成立。等号不成立，因为 $\mathbf{x}$ 与 $\mathbf{y}$ 不成比例。

### 4.2 三角不等式（向量版）

**定理 11.2（三角不等式）**：对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，
$$\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$$

**证明**：

我们从 $\|\mathbf{x} + \mathbf{y}\|^2$ 出发，利用内积展开，然后应用 Cauchy-Schwarz 不等式。

**第一步**：展开平方。
$$
\begin{aligned}
\|\mathbf{x} + \mathbf{y}\|^2 &= \langle \mathbf{x} + \mathbf{y},\ \mathbf{x} + \mathbf{y} \rangle \\
&= \langle \mathbf{x}, \mathbf{x} \rangle + \langle \mathbf{x}, \mathbf{y} \rangle + \langle \mathbf{y}, \mathbf{x} \rangle + \langle \mathbf{y}, \mathbf{y} \rangle \\
&= \|\mathbf{x}\|^2 + 2\langle \mathbf{x}, \mathbf{y} \rangle + \|\mathbf{y}\|^2
\end{aligned}
$$
其中用到 $\langle \mathbf{x}, \mathbf{y} \rangle = \langle \mathbf{y}, \mathbf{x} \rangle$。

**第二步**：应用 Cauchy-Schwarz 不等式。
$$\langle \mathbf{x}, \mathbf{y} \rangle \leq |\langle \mathbf{x}, \mathbf{y} \rangle| \leq \|\mathbf{x}\| \cdot \|\mathbf{y}\|$$

这里第一个 $\leq$ 是因为任何实数不超过其绝对值（对任意 $a \in \mathbb{R}$，有 $a \leq |a|$），第二个 $\leq$ 是 Cauchy-Schwarz 不等式。

因此：
$$
\begin{aligned}
\|\mathbf{x} + \mathbf{y}\|^2 &= \|\mathbf{x}\|^2 + 2\langle \mathbf{x}, \mathbf{y} \rangle + \|\mathbf{y}\|^2 \\
&\leq \|\mathbf{x}\|^2 + 2\|\mathbf{x}\| \cdot \|\mathbf{y}\| + \|\mathbf{y}\|^2 \\
&= (\|\mathbf{x}\| + \|\mathbf{y}\|)^2
\end{aligned}
$$

**第三步**：取平方根。

$\|\mathbf{x} + \mathbf{y}\| \geq 0$ 且 $\|\mathbf{x}\| + \|\mathbf{y}\| \geq 0$，因此对平方不等式两边取非负平方根，不等号方向保持不变：
$$\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$$

证毕。

**关于等号成立条件**：三角不等式的等号成立当且仅当 $\langle \mathbf{x}, \mathbf{y} \rangle = \|\mathbf{x}\| \cdot \|\mathbf{y}\|$ 且 $\langle \mathbf{x}, \mathbf{y} \rangle \geq 0$。由 Cauchy-Schwarz 等号条件，这等价于 $\mathbf{x} = \lambda\mathbf{y}$ 且 $\lambda \geq 0$（即两个向量同方向）。读者可自行验证。

**与一维三角不等式的关系**：当 $n=1$ 时，$\|x+y\| = |x+y|$，$\|x\| + \|y\| = |x| + |y|$，于是 $\|\mathbf{x}+\mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$ 退化为一维三角不等式 $|x+y| \leq |x| + |y|$。

### 4.3 $\mathbb{R}^n$ 中点列收敛的 $\varepsilon$-$N$ 定义

我们现在将第二章中一维数列极限的 $\varepsilon$-$N$ 定义推广到 $\mathbb{R}^n$。

**定义 11.6（$\mathbb{R}^n$ 中点列的收敛）**：设 $\{P_k\}_{k=1}^\infty$ 是 $\mathbb{R}^n$ 中的一个点列（即每个 $P_k \in \mathbb{R}^n$），$P \in \mathbb{R}^n$ 是一个固定点。若对任意给定的 $\varepsilon > 0$，总存在正整数 $N$（依赖于 $\varepsilon$），使得当 $k > N$ 时恒有
$$\|P_k - P\| < \varepsilon$$
则称点列 $\{P_k\}$ **收敛**（converges）于 $P$，记作
$$\lim_{k\to\infty} P_k = P \quad \text{或} \quad P_k \to P \ (k \to \infty)$$
若 $\{P_k\}$ 不存在极限，则称该点列**发散**（diverges）。

**与一维定义的比较**：

一维数列极限定义：
$$\lim_{n\to\infty} a_n = a \iff \forall\varepsilon>0,\ \exists N\in\mathbb{N}^+,\ \forall n>N:\ |a_n - a| < \varepsilon$$

$\mathbb{R}^n$ 点列收敛定义：
$$\lim_{k\to\infty} P_k = P \iff \forall\varepsilon>0,\ \exists N\in\mathbb{N}^+,\ \forall k>N:\ \|P_k - P\| < \varepsilon$$

**比较分析**：
- 两个定义的结构完全相同——都是 $\forall\varepsilon>0,\ \exists N,\ \forall k>N$ 的三段式
- 唯一的变化是：一维中的绝对值 $|a_n - a|$ 被替换为 $n$ 维中的范数 $\|P_k - P\|$
- 这并非偶然——范数 $\|\cdot\|$ 在 $\mathbb{R}^n$ 中扮演的角色与绝对值 $|\cdot|$ 在 $\mathbb{R}$ 中扮演的角色完全相同：它们都度量"距离"（定义 11.5）

**说明**：由于 $\|P_k - P\|$ 是一个非负实数，$\varepsilon$-$N$ 定义中的所有**数值不等式操作**（如放缩、找 $N$）都与一维完全相同。区别仅在于 $P_k - P$ 现在是一个向量，其"大小"由范数而非绝对值度量。

**例 7（发散点列）**：$P_k = ((-1)^k,\ 0) \in \mathbb{R}^2$。点列在 $(1,0)$ 和 $(-1,0)$ 之间来回跳动，不收敛到任何固定点（证明思路与一维中 $a_n = (-1)^n$ 完全类似）。

### 4.4 分量收敛与整体收敛的等价关系

这是本章最重要的定理之一。它告诉我们：$\mathbb{R}^n$ 中点列的收敛等价于每个分量的数列的收敛。这使得我们可以将 $n$ 维收敛问题拆解为 $n$ 个一维收敛问题分别处理。

**定理 11.3（分量收敛等价定理）**：设 $P_k = (x_k^{(1)}, x_k^{(2)}, \dots, x_k^{(n)}) \in \mathbb{R}^n$ 是一个点列，$P = (x^{(1)}, x^{(2)}, \dots, x^{(n)}) \in \mathbb{R}^n$。则
$$P_k \to P \ \text{在} \ \mathbb{R}^n \text{中} \iff \text{对每个 } i = 1,\dots,n,\ \text{有 } x_k^{(i)} \to x^{(i)} \ \text{在} \ \mathbb{R} \text{中}$$

即点列收敛当且仅当每个分量数列都收敛。

**证明**：我们需要证明两个方向。

**（$\Rightarrow$ 方向）整体收敛 $\Rightarrow$ 每个分量收敛**：

假设 $P_k \to P$ 在 $\mathbb{R}^n$ 中成立。由定义 11.6，对任意 $\varepsilon > 0$，存在 $N \in \mathbb{N}^+$，使得当 $k > N$ 时，
$$\|P_k - P\| < \varepsilon$$

现在固定某个分量指标 $i \in \{1,2,\dots,n\}$。注意到：
$$\|P_k - P\| = \sqrt{\sum_{j=1}^n (x_k^{(j)} - x^{(j)})^2} \geq \sqrt{(x_k^{(i)} - x^{(i)})^2} = |x_k^{(i)} - x^{(i)}|$$

这里不等号成立是因为求和符号中的每一项 $(x_k^{(j)} - x^{(j)})^2$ 都是非负的，去掉除 $j=i$ 以外的所有非负项，和只可能变小或不变。

因此，当 $k > N$ 时，
$$|x_k^{(i)} - x^{(i)}| \leq \|P_k - P\| < \varepsilon$$

由于 $\varepsilon$ 的任意性，这正是 $x_k^{(i)} \to x^{(i)}$ 在 $\mathbb{R}$ 中的定义。每个分量都成立，故所有分量数列收敛。

**（$\Leftarrow$ 方向）每个分量收敛 $\Rightarrow$ 整体收敛**：

假设对每个 $i = 1,\dots,n$，有 $x_k^{(i)} \to x^{(i)}$ 在 $\mathbb{R}$ 中成立。由数列极限的定义，对任意 $\varepsilon > 0$，对每个 $i$，存在 $N_i \in \mathbb{N}^+$，使得当 $k > N_i$ 时，
$$|x_k^{(i)} - x^{(i)}| < \frac{\varepsilon}{n}$$

现在取 $N = \max\{N_1, N_2, \dots, N_n\}$。当 $k > N$ 时，对所有 $i = 1,\dots,n$ 同时有 $|x_k^{(i)} - x^{(i)}| < \varepsilon/n$。

下面估计 $\|P_k - P\|$：

$$
\begin{aligned}
\|P_k - P\| &= \sqrt{\sum_{i=1}^n (x_k^{(i)} - x^{(i)})^2} \\
&\leq \sum_{i=1}^n |x_k^{(i)} - x^{(i)}| \quad \text{（理由见下方说明）} \\
&< \sum_{i=1}^n \frac{\varepsilon}{n} = \varepsilon
\end{aligned}
$$

因此 $\|P_k - P\| < \varepsilon$ 对任意 $k > N$ 成立，即 $P_k \to P$ 在 $\mathbb{R}^n$ 中。

**关于 $\sqrt{\sum a_i^2} \leq \sum |a_i|$ 的说明**：证明中需要这一不等式来将分量估计"合并"为整体估计。验证如下：对任意实数 $a_1,\dots,a_n$，
$$\left(\sum_{i=1}^n |a_i|\right)^2 = \sum_{i=1}^n a_i^2 + 2\sum_{1 \leq i < j \leq n} |a_i||a_j| \geq \sum_{i=1}^n a_i^2$$
两边取非负平方根即得 $\sqrt{\sum a_i^2} \leq \sum |a_i|$。

证毕。

**定理 11.3 的实用价值**：该定理表明 $\mathbb{R}^n$ 中收敛性的判定可以归结为 $n$ 个一维收敛性的判定。例如，要判断 $P_k = \left(1 + \frac{1}{k},\ \frac{k}{k+1},\ \frac{1}{k^2}\right)$ 是否收敛，只需判断三个分量数列 $x_k^{(1)} = 1 + 1/k$，$x_k^{(2)} = k/(k+1)$，$x_k^{(3)} = 1/k^2$ 是否分别收敛。

---

## 5. 例题

### 例题 1：内积、范数与距离的计算

给定 $\mathbf{x} = (1, -2, 0, 3)$，$\mathbf{y} = (2, 1, -1, 1)$（均为 $\mathbb{R}^4$ 中的向量）。

(1) 计算 $\langle \mathbf{x}, \mathbf{y} \rangle$。
(2) 计算 $\|\mathbf{x}\|$ 和 $\|\mathbf{y}\|$。
(3) 计算 $d(\mathbf{x}, \mathbf{y})$。

**解**：

(1) 内积：
$$\langle \mathbf{x}, \mathbf{y} \rangle = 1\cdot 2 + (-2)\cdot 1 + 0\cdot (-1) + 3\cdot 1 = 2 - 2 + 0 + 3 = 3$$

(2) 范数：
$$\|\mathbf{x}\| = \sqrt{1^2 + (-2)^2 + 0^2 + 3^2} = \sqrt{1 + 4 + 0 + 9} = \sqrt{14}$$
$$\|\mathbf{y}\| = \sqrt{2^2 + 1^2 + (-1)^2 + 1^2} = \sqrt{4 + 1 + 1 + 1} = \sqrt{7}$$

(3) 距离。先计算差值向量：
$$\mathbf{x} - \mathbf{y} = (1-2,\ -2-1,\ 0-(-1),\ 3-1) = (-1, -3, 1, 2)$$
欧氏距离：
$$d(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\| = \sqrt{(-1)^2 + (-3)^2 + 1^2 + 2^2} = \sqrt{1 + 9 + 1 + 4} = \sqrt{15}$$

**验证**：由三角不等式，应有 $d(\mathbf{x}, \mathbf{y}) \leq \|\mathbf{x}\| + \|\mathbf{y}\|$。计算得 $\sqrt{15} \approx 3.873$，$\sqrt{14} + \sqrt{7} \approx 3.742 + 2.646 = 6.388$，不等式成立。

---

### 例题 2：Cauchy-Schwarz 不等式的数值验证与等号条件

(1) 验证 CS 不等式对 $\mathbf{x} = (1, 2, 3)$ 和 $\mathbf{y} = (1, 1, 1)$ 成立。
(2) 找出使 CS 不等式取等号的 $\mathbf{x}$ 和 $\mathbf{y}$ 的例子。

**解**：

(1) 验证：
$$\langle \mathbf{x}, \mathbf{y} \rangle = 1\cdot 1 + 2\cdot 1 + 3\cdot 1 = 6$$
$$\|\mathbf{x}\| = \sqrt{1^2 + 2^2 + 3^2} = \sqrt{14} \approx 3.742$$
$$\|\mathbf{y}\| = \sqrt{1^2 + 1^2 + 1^2} = \sqrt{3} \approx 1.732$$
$$\|\mathbf{x}\| \cdot \|\mathbf{y}\| = \sqrt{14} \cdot \sqrt{3} = \sqrt{42} \approx 6.481$$

$|\langle \mathbf{x}, \mathbf{y} \rangle| = 6 \leq 6.481$，不等式成立。

(2) 等号成立的条件是 $\mathbf{x} = \lambda \mathbf{y}$（或 $\mathbf{y} = \lambda \mathbf{x}$）。例如取 $\mathbf{x} = (2, -4, 6)$，$\mathbf{y} = (1, -2, 3)$，则 $\mathbf{x} = 2\mathbf{y}$。

验证：$\langle \mathbf{x}, \mathbf{y} \rangle = 2\cdot 1 + (-4)\cdot(-2) + 6\cdot 3 = 2 + 8 + 18 = 28$
$$\|\mathbf{x}\| = \sqrt{4 + 16 + 36} = \sqrt{56} = 2\sqrt{14}$$
$$\|\mathbf{y}\| = \sqrt{1 + 4 + 9} = \sqrt{14}$$
$$\|\mathbf{x}\| \cdot \|\mathbf{y}\| = 2\sqrt{14} \cdot \sqrt{14} = 2 \times 14 = 28$$
等号 $28 = 28$ 成立。

---

### 例题 3：用 $\varepsilon$-$N$ 定义证明点列收敛

设 $P_k = \left(\frac{1}{k},\ \frac{k-1}{k}\right) \in \mathbb{R}^2$。证明 $P_k \to (0, 1)$。

**解法一（直接证明）**：

**分析段**：先计算 $\|P_k - P\|$，其中 $P = (0,1)$。
$$P_k - P = \left(\frac{1}{k} - 0,\ \frac{k-1}{k} - 1\right) = \left(\frac{1}{k},\ -\frac{1}{k}\right)$$
$$\|P_k - P\| = \sqrt{\left(\frac{1}{k}\right)^2 + \left(-\frac{1}{k}\right)^2} = \sqrt{\frac{1}{k^2} + \frac{1}{k^2}} = \sqrt{\frac{2}{k^2}} = \frac{\sqrt{2}}{k}$$

需要 $\frac{\sqrt{2}}{k} < \varepsilon$，即 $k > \frac{\sqrt{2}}{\varepsilon}$。取 $N = \left\lceil\frac{\sqrt{2}}{\varepsilon}\right\rceil$。

**正序证明段**：

任取 $\varepsilon > 0$。取 $N = \left\lceil\frac{\sqrt{2}}{\varepsilon}\right\rceil \in \mathbb{N}^+$。当 $k > N$ 时，
$$\|P_k - (0,1)\| = \sqrt{\left(\frac{1}{k}\right)^2 + \left(-\frac{1}{k}\right)^2} = \frac{\sqrt{2}}{k} < \frac{\sqrt{2}}{N} \leq \varepsilon$$
故 $\displaystyle\lim_{k\to\infty} P_k = (0,1)$。证毕。

**解法二（利用分量收敛定理）**：

$x_k^{(1)} = \frac{1}{k} \to 0$（已知 $\lim \frac{1}{k} = 0$，第二章例 5）。
$x_k^{(2)} = \frac{k-1}{k} = 1 - \frac{1}{k} \to 1$（由极限运算法则，$\lim (1 - \frac{1}{k}) = 1$）。

两个分量数列均收敛，由定理 11.3，$P_k \to (0,1)$。

**注**：解法二利用了定理 11.3，将二维问题分解为两个一维问题，更加简便。这正体现了分量收敛等价定理的实用价值。

---

## 6. 常见误区与检查点

### 常见误区

| 错误理解 | 正确理解 |
|----------|----------|
| 认为 $\langle \mathbf{x}, \mathbf{y} \rangle = \|\mathbf{x}\| \cdot \|\mathbf{y}\|$ 恒成立 | $\langle \mathbf{x}, \mathbf{y} \rangle = \|\mathbf{x}\| \cdot \|\mathbf{y}\|$ 仅当 $\mathbf{x}$ 与 $\mathbf{y}$ 同方向时成立（CS 不等式等号条件）。一般情况下 $|\langle \mathbf{x}, \mathbf{y} \rangle| \leq \|\mathbf{x}\| \cdot \|\mathbf{y}\|$ |
| 认为 $\mathbb{R}^n$ 中的 $\varepsilon$-$N$ 定义与一维定义有本质区别 | 两个定义的结构完全一致——都是 $\forall\varepsilon>0,\ \exists N,\ \forall k>N$ 三段式，唯一的区别是将绝对值 $|\cdot|$ 替换为范数 $\|\cdot\|$ |
| 试图用一维的绝对值代替 $\mathbb{R}^n$ 中的范数来定义收敛——如 $\forall\varepsilon>0$，要求 $|x_k^{(1)} - x^{(1)}| + \cdots + |x_k^{(n)} - x^{(n)}| < \varepsilon$ | 标准定义使用的是 $\|P_k - P\| = \sqrt{\sum (x_k^{(i)} - x^{(i)})^2} < \varepsilon$。不过，由 $\sqrt{\sum a_i^2} \leq \sum |a_i|$ 和 $\sum |a_i| \leq \sqrt{n} \cdot \sqrt{\sum a_i^2}$（Cauchy-Schwarz 不等式的推论），这两种定义实际上等价 |
| 混淆 $\langle \mathbf{x}, \mathbf{y} \rangle$ 与 $\mathbf{x} \cdot \mathbf{y}$ 的记号——认为它们是不同的概念 | $\langle \mathbf{x}, \mathbf{y} \rangle$ 和 $\mathbf{x} \cdot \mathbf{y}$ 在本章中代表同一个东西：Euclid 内积。$\langle \cdot, \cdot \rangle$ 是一般内积的记号，$\cdot$ 是点积的记号，在 Euclid 空间中二者一致 |
| 认为 $\|\mathbf{x}\| = \sqrt{\sum x_i^2}$ 是 $\mathbb{R}^n$ 上唯一的范数 | 这是 Euclid 范数（2-范数），是最常用的范数，但不是唯一的。$\mathbb{R}^n$ 上还可以定义 $p$-范数 $\|\mathbf{x}\|_p = (\sum |x_i|^p)^{1/p}$，其中 $p \geq 1$。本教材仅使用 Euclid 范数 |

### 检查点

- [ ] 能否写出 $\mathbb{R}^n$ 中向量加法和标量乘法的定义？
- [ ] 能否解释标准基 $\{\mathbf{e}_1,\dots,\mathbf{e}_n\}$ 的作用？
- [ ] 能否写出 Euclid 内积的定义并验证其三条性质（对称性、双线性、正定性）？
- [ ] 能否写出 Cauchy-Schwarz 不等式的判别式证明中的关键步骤？
- [ ] 能否从 Cauchy-Schwarz 不等式推导出三角不等式？
- [ ] 能否写出 $\mathbb{R}^n$ 中点列收敛的 $\varepsilon$-$N$ 定义，并与一维定义逐项对比？
- [ ] 能否完整重述并证明分量收敛等价定理（两个方向）？
- [ ] 是否理解 $\sqrt{\sum a_i^2} \leq \sum |a_i|$ 在分量收敛证明中的作用？
- [ ] 能否举例说明一个 $\mathbb{R}^2$ 中的发散点列？

---

## 练习题

### 基础巩固

**1.** 在 $\mathbb{R}^3$ 中，给定 $\mathbf{x} = (2, -1, 3)$，$\mathbf{y} = (0, 4, -2)$。

(1) 计算 $\mathbf{x} + \mathbf{y}$ 和 $2\mathbf{x} - 3\mathbf{y}$。
(2) 计算 $\langle \mathbf{x}, \mathbf{y} \rangle$。
(3) 计算 $\|\mathbf{x}\|$，$\|\mathbf{y}\|$。
(4) 计算 $d(\mathbf{x}, \mathbf{y})$。

<details><summary>参考答案</summary>

(1) $$\mathbf{x} + \mathbf{y} = (2+0,\ -1+4,\ 3+(-2)) = (2, 3, 1)$$
$$2\mathbf{x} - 3\mathbf{y} = (4, -2, 6) - (0, 12, -6) = (4, -14, 12)$$

(2) $$\langle \mathbf{x}, \mathbf{y} \rangle = 2\cdot 0 + (-1)\cdot 4 + 3\cdot (-2) = 0 - 4 - 6 = -10$$

(3) $$\|\mathbf{x}\| = \sqrt{2^2 + (-1)^2 + 3^2} = \sqrt{4 + 1 + 9} = \sqrt{14}$$
$$\|\mathbf{y}\| = \sqrt{0^2 + 4^2 + (-2)^2} = \sqrt{0 + 16 + 4} = \sqrt{20} = 2\sqrt{5}$$

(4) $$\mathbf{x} - \mathbf{y} = (2-0,\ -1-4,\ 3-(-2)) = (2, -5, 5)$$
$$d(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\| = \sqrt{2^2 + (-5)^2 + 5^2} = \sqrt{4 + 25 + 25} = \sqrt{54} = 3\sqrt{6}$$

</details>

---

**2.** 验证 Cauchy-Schwarz 不等式对 $\mathbf{x} = (1, -1, 2, 0)$ 和 $\mathbf{y} = (2, 3, -1, 1)$（$\mathbb{R}^4$ 中）成立。

<details><summary>参考答案</summary>

$$\langle \mathbf{x}, \mathbf{y} \rangle = 1\cdot 2 + (-1)\cdot 3 + 2\cdot (-1) + 0\cdot 1 = 2 - 3 - 2 + 0 = -3$$
$$|\langle \mathbf{x}, \mathbf{y} \rangle| = 3$$

$$\|\mathbf{x}\| = \sqrt{1^2 + (-1)^2 + 2^2 + 0^2} = \sqrt{1 + 1 + 4 + 0} = \sqrt{6} \approx 2.449$$
$$\|\mathbf{y}\| = \sqrt{2^2 + 3^2 + (-1)^2 + 1^2} = \sqrt{4 + 9 + 1 + 1} = \sqrt{15} \approx 3.873$$
$$\|\mathbf{x}\| \cdot \|\mathbf{y}\| = \sqrt{6} \cdot \sqrt{15} = \sqrt{90} = 3\sqrt{10} \approx 9.487$$

$3 \leq 9.487$，不等式成立。等号不成立（$\mathbf{x}$ 与 $\mathbf{y}$ 不成比例）。

</details>

---

**3.** 用 $\varepsilon$-$N$ 定义证明 $P_k = \left(\frac{1}{k^2},\ \frac{2k+1}{k+1}\right) \in \mathbb{R}^2$ 收敛，并求出其极限。

<details><summary>参考答案</summary>

**猜测极限**：
$$x_k^{(1)} = \frac{1}{k^2} \to 0, \quad x_k^{(2)} = \frac{2k+1}{k+1} \to 2$$
由定理 11.3，猜测 $P_k \to (0, 2)$。

**直接证明**：

计算差值向量：
$$P_k - (0,2) = \left(\frac{1}{k^2},\ \frac{2k+1}{k+1} - 2\right) = \left(\frac{1}{k^2},\ \frac{2k+1 - 2(k+1)}{k+1}\right) = \left(\frac{1}{k^2},\ \frac{-1}{k+1}\right)$$

范数：
$$\|P_k - (0,2)\| = \sqrt{\frac{1}{k^4} + \frac{1}{(k+1)^2}}$$
做放缩：
$$\|P_k - (0,2)\| \leq \sqrt{\frac{1}{k^4} + \frac{1}{k^2}} = \frac{1}{k}\sqrt{\frac{1}{k^2} + 1} \leq \frac{1}{k}\sqrt{2} \quad (k \geq 1)$$

任取 $\varepsilon > 0$。需要 $\frac{\sqrt{2}}{k} < \varepsilon$，即 $k > \frac{\sqrt{2}}{\varepsilon}$。
取 $N = \left\lceil\frac{\sqrt{2}}{\varepsilon}\right\rceil$。当 $k > N$ 时，
$$\|P_k - (0,2)\| \leq \frac{\sqrt{2}}{k} < \frac{\sqrt{2}}{N} \leq \varepsilon$$
故 $P_k \to (0,2)$。

**或利用分量收敛定理**：
$\frac{1}{k^2} \to 0$（已知），$\frac{2k+1}{k+1} = 2 - \frac{1}{k+1} \to 2$（由极限运算法则）。两个分量均收敛，由定理 11.3 知 $P_k \to (0,2)$。

</details>

---

### 迁移应用

**4.** 设 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$。证明平行四边形定律：
$$\|\mathbf{x} + \mathbf{y}\|^2 + \|\mathbf{x} - \mathbf{y}\|^2 = 2(\|\mathbf{x}\|^2 + \|\mathbf{y}\|^2)$$
并给出几何解释。

<details><summary>参考答案</summary>

**证明**：
利用内积展开两个平方项：
$$
\begin{aligned}
\|\mathbf{x} + \mathbf{y}\|^2 &= \langle \mathbf{x} + \mathbf{y},\ \mathbf{x} + \mathbf{y} \rangle = \|\mathbf{x}\|^2 + 2\langle \mathbf{x}, \mathbf{y} \rangle + \|\mathbf{y}\|^2 \\
\|\mathbf{x} - \mathbf{y}\|^2 &= \langle \mathbf{x} - \mathbf{y},\ \mathbf{x} - \mathbf{y} \rangle = \|\mathbf{x}\|^2 - 2\langle \mathbf{x}, \mathbf{y} \rangle + \|\mathbf{y}\|^2
\end{aligned}
$$

两式相加：
$$\|\mathbf{x} + \mathbf{y}\|^2 + \|\mathbf{x} - \mathbf{y}\|^2 = 2\|\mathbf{x}\|^2 + 2\|\mathbf{y}\|^2 = 2(\|\mathbf{x}\|^2 + \|\mathbf{y}\|^2)$$

**几何解释**（$\mathbb{R}^2$ 或 $\mathbb{R}^3$ 中）：考虑以 $\mathbf{x}$ 和 $\mathbf{y}$ 为边的平行四边形，$\mathbf{x}+\mathbf{y}$ 和 $\mathbf{x}-\mathbf{y}$ 是其两条对角线。平行四边形定律表明：两条对角线的平方和等于四条边的平方和。

</details>

---

**5.** 设 $\{P_k\}$ 是 $\mathbb{R}^n$ 中的收敛点列，$P_k \to P$。证明 $\{P_k\}$ 有界，即存在 $M > 0$ 使得 $\|P_k\| \leq M$ 对所有 $k \in \mathbb{N}^+$ 成立。（提示：参考第二章中"收敛数列必有界"的证明思路，将绝对值替换为范数。）

<details><summary>参考答案</summary>

**证明**：

由 $P_k \to P$，取 $\varepsilon = 1$，则存在 $N \in \mathbb{N}^+$，当 $k > N$ 时
$$\|P_k - P\| < 1$$

由三角不等式（定理 11.2），对 $k > N$：
$$\|P_k\| = \|P_k - P + P\| \leq \|P_k - P\| + \|P\| < 1 + \|P\|$$

对于前 $N$ 项 $P_1, P_2, \dots, P_N$，这是有限个向量，它们的范数有最大值：
$$M_0 = \max\{\|P_1\|, \|P_2\|, \dots, \|P_N\|\}$$

取 $M = \max\{M_0,\ 1 + \|P\|\}$，则对所有 $k \in \mathbb{N}^+$，有 $\|P_k\| \leq M$。证毕。

**与一维证明的对比**：这个证明完全平行于第二章定理 8.2（收敛数列必有界）的证明，仅将绝对值 $|\cdot|$ 替换为范数 $\|\cdot\|$，三角不等式 $|x+y| \leq |x| + |y|$ 替换为 $\|\mathbf{x}+\mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$。

</details>

---

**6.** 证明反向三角不等式：对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，
$$\big|\|\mathbf{x}\| - \|\mathbf{y}\|\big| \leq \|\mathbf{x} - \mathbf{y}\|$$

<details><summary>参考答案</summary>

**证明**：

由三角不等式 $\|\mathbf{x}\| = \|\mathbf{x} - \mathbf{y} + \mathbf{y}\| \leq \|\mathbf{x} - \mathbf{y}\| + \|\mathbf{y}\|$，移项得
$$\|\mathbf{x}\| - \|\mathbf{y}\| \leq \|\mathbf{x} - \mathbf{y}\|$$

同理，交换 $\mathbf{x}$ 和 $\mathbf{y}$ 的角色：$\|\mathbf{y}\| - \|\mathbf{x}\| \leq \|\mathbf{y} - \mathbf{x}\| = \|\mathbf{x} - \mathbf{y}\|$，即
$$\|\mathbf{y}\| - \|\mathbf{x}\| \leq \|\mathbf{x} - \mathbf{y}\|$$

由于 $\big|\|\mathbf{x}\| - \|\mathbf{y}\|\big|$ 要么等于 $\|\mathbf{x}\| - \|\mathbf{y}\|$（当 $\|\mathbf{x}\| \geq \|\mathbf{y}\|$ 时），要么等于 $\|\mathbf{y}\| - \|\mathbf{x}\|$（当 $\|\mathbf{y}\| \geq \|\mathbf{x}\|$ 时），而两种情况下都有 $\big|\|\mathbf{x}\| - \|\mathbf{y}\|\big| \leq \|\mathbf{x} - \mathbf{y}\|$。证毕。

这个不等式的几何意义是：两条线段长度之差的绝对值不超过它们端点之间的距离。

</details>
