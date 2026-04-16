# 第十六章 Fourier 级数

## § 1 中 Achilles 追乌龟的路程所构成的级数 $ S_{1}(1+q+q^{2}+\cdots+q^{n-1}+q^{n}+\cdots) $，例 9.1.7 中的级数 $ 0.05\sum_{n=0}^{\infty}\left(\frac{4}{5}\right)^{n} $ 及我们熟知的 p 级数 $ \sum_{n=0}^{\infty}\frac{1}{n^{p}} $，均有一个明显的特征：它们的各个项都是正数.

### 1. (1) $ \frac{A}{\pi} + \frac{A}{2} \sin x - \frac{2A}{\pi} \sum_{k=1}^{\infty} \frac{\cos 2kx}{4k^{2} - 1} $.

### $$ \frac{2A}{\pi}-\frac{4A}{\pi}\sum_{k=1}^{\infty}\frac{\cos2kx}{4k^{2}-1}. $$

### 2. (1)

### $$ \frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\sin(2k-1)x}{2k-1}. $$

### $$ \frac{2}{\pi}-\frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\left(-1\right)^{k}}{4k^{2}-1}\cos2kx. $$

### (3)

### $$ (3)-\frac{5}{6}\pi^{2}+\sum_{n=1}^{\infty}\frac{2\left(-1\right)^{n}}{n^{2}}\cos nx. $$

### (4)

### $$ -\frac{\pi}{4}+\frac{2}{\pi}\sum_{k=0}^{\infty}\frac{\cos\left(2k+1\right)x}{\left(2k+1\right)^{2}}+\sum_{n=1}^{\infty}\frac{\left(-1\right)^{n+1}}{n}\sin nx. $$

### $$ -\frac{\left(a-b\right)\pi}{4}+\frac{2\left(a-b\right)}{\pi}\sum_{k=0}^{\infty}\frac{\cos\left(2k+1\right)x}{\left(2k+1\right)^{2}}+\left(a+b\right)\sum_{n=1}^{\infty}\frac{\left(-1\right)^{n+1}}{n}\sin nx. $$

### 3. (1) $ 2 \sum_{n=1}^{\infty} \frac{\left[1 - 2(-1)^{n}\right]}{n} \sin nx. $

### (2)

### $$ \frac{2}{\pi}\sum_{n=1}^{\infty}\frac{n\left[1-(-1)^{n}\mathrm{e}^{-2\pi}\right]}{n^{2}+4}\sin nx. $$

### (3)

### $$ \sum_{n=1}^{\infty}\left[\frac{2}{n}(-1)^{n+1}+\frac{4}{\pi n^{2}}\sin\frac{n\pi}{2}\right]\sin nx. $$

### (4)

### $$ \frac{1}{\pi}\sin\frac{\pi}{2}x+\frac{2}{\pi}\sum_{n=2}^{\infty}\frac{n-\sin\frac{n\pi}{2}}{n^{2}-1}\sin\frac{n\pi}{2}x. $$

### 4. (1)

### $$ \frac{\pi^{2}}{6}-\sum_{k=1}^{\infty}\frac{\cos2kx}{k^{2}}. $$

### (2) $ \frac{1}{\pi}(\mathrm{e}^{\pi} - 1) + \frac{2}{\pi}\sum_{n=1}^{\infty}\frac{\left[(-1)^{n}\mathrm{e}^{\pi} - 1\right]}{n^{2} + 1}\cos nx. $

### (3) $ \left(\frac{1}{\pi} + \frac{1}{2}\right) - \frac{1}{\pi} \cos 2x - \frac{2}{\pi} \sum_{n=2}^{\infty} \frac{1}{n^{2} - 1} \left(\frac{1}{n} \sin \frac{n\pi}{2} - 1\right) \cos 2nx. $

### (4) $ \frac{\pi}{4} + \frac{4}{\pi} \sum_{n=1}^{\infty} \frac{\left[(-1)^{n} - \cos \frac{n\pi}{2}\right]}{n^{2}} \cos nx. $

### 5. $ f(x) \sim \frac{a_{0}}{2} + \sum_{n=1}^{\infty} (a_{n} \cos nx + b_{n} \sin nx) $，其中

### $$ a_{n}=\frac{1}{\pi}\int_{a}^{n+2\pi}f(x)\cos n x\mathrm{d}x(n=0,1,2,\cdots), $$

### $$ b_{_{n}}=\frac{1}{\pi}\int_{_{a}}^{^{a+2\pi}}f(x)\sin n x\mathrm{d}x(n=1,2,\cdots). $$

### 6. (1) $ \sum_{n=1}^{\infty}\frac{1}{n}\sin nx. $

### (2) $ \frac{4}{3}\pi^{2}+4\sum_{n=1}^{\infty}\left(\frac{1}{n^{2}}\cos nx-\frac{\pi}{n}\sin nx\right) $

### (3) $ \frac{1}{2} - \frac{1}{\pi} \sum_{n=1}^{\infty} \frac{1}{n} \sin 2\pi nx. $

### (4)

### $$ \frac{1}{6}\left(1-\mathrm{e}^{-3}\right)+\sum_{n=1}^{\infty}\left[\frac{3\left(1-\left(-1\right)^{n}\mathrm{e}^{-3}\right)}{n^{2}\pi^{2}+9}\cos n\pi x-\frac{n\pi\left(1-\left(-1\right)^{n}\mathrm{e}^{-3}\right)}{n^{2}\pi^{2}+9}\sin n\pi x\right]. $$

### (5) $ \frac{C}{2} - \frac{2C}{\pi} \sum_{n=1}^{\infty} \frac{1}{2n - 1} \sin \frac{(2n - 1)\pi}{T} x. $

### 7. $ -\frac{5}{4\pi}(2-\sqrt{2})-\frac{5}{4\pi}\cos\omega t+\left(\frac{5}{4\pi}+\frac{35}{8}\right)\sin\omega t+ $

### $$ \frac{5}{2\pi}\sum_{n=2}^{\infty}\left[\frac{1}{n+1}\cos\frac{\left(n+1\right)\pi}{4}-\frac{1}{n-1}\cos\frac{\left(n-1\right)\pi}{4}+\frac{2}{n^{2}-1}\right]\cos n\omega t+ $$

### $$ \frac{5}{2\pi}\sum_{n=2}^{\infty}\left[\frac{1}{n+1}\sin\frac{\left(n+1\right)\pi}{4}-\frac{1}{n-1}\sin\frac{\left(n-1\right)\pi}{4}\right]\sin n\omega t. $$

### 9. (1) $ \tilde{f}(x)=\left\{\begin{aligned}-f(\pi+x),&x\in\left(-\pi,-\frac{\pi}{2}\right),\\f(-x),&x\in\left(-\frac{\pi}{2},0\right),\\f(x),&x\in\left(0,\frac{\pi}{2}\right),\\-f(\pi-x),&x\in\left(\frac{\pi}{2},\pi\right).\end{aligned}\right. $

### (2)

### $$ \begin{aligned}\tilde{f}(x)=\left\{\begin{aligned}&f(\pi+x),&x\in\left(-\pi,-\frac{\pi}{2}\right)\\&-f(-x),&x\in\left(-\frac{\pi}{2},0\right),\\&f(x),&x\in\left(0,\frac{\pi}{2}\right),\\&-f(\pi-x),&x\in\left(\frac{\pi}{2},\pi\right).\end{aligned}\right.\end{aligned} $$

### 10. (1) $ \tilde{a}_{n} = a_{n} (n = 0, 1, 2, \cdots) $, $ \tilde{b}_{n} = -b_{n} (n = 1, 2, \cdots) $.

### (2) $ \tilde{a}_{n}=a_{n}\cos nC+b_{n}\sin nC\quad(n=0,1,2,\cdots) $

### $$ \hat{b}_{n}=b_{n}\cos n C-a_{n}\sin n C\quad(n=1,2,\cdots). $$

### (3) $ \bar{a}_{0}=a_{0}^{2},\bar{a}_{n}=a_{n}^{2}-b_{n}^{2},\bar{b}_{n}=2a_{n}\bar{b}_{n}(n=1,2,\cdots) $

## § 2 Fourier 级数的收敛判别法

### 1. 提示: 因为 $ \lim_{x \to +\infty} \psi(x) = 0 $, 所以存在 N > 0, 使得当 $ x \geq N $ 时, $ \left| \psi(x) \right| < 1 $. 利用积分第二中值定理可得 $ \left| \int_N^A \psi(x) \sin px dx \right| < \frac{4}{p} (\forall A > N) $, 因此 $ \left| \int_N^{+\infty} \psi(x) \sin px dx \right| \leq \frac{4}{p} $. 而由 Riemann 引理, $ \lim_{p \to +\infty} \int_0^N \psi(x) \sin px dx = 0 $. 因此当 $ p \to +\infty $ 时, $ \int_-\infty^+\infty \psi(x) \sin px dx = \int_-\infty^N \psi(x) \sin px dx + \int_-\infty^+\infty \psi(x) \sin px dx \to 0 $.

### 2. 提示: 易知

### $$ \int_{-\pi}^{\pi}\psi\left(u\right)\frac{\cos\frac{u}{2}-\cos pu}{2\sin\frac{u}{2}}\mathrm{d}u=\int_{0}^{\pi}\left[\psi\left(u\right)-\psi\left(-u\right)\right]\frac{\cos\frac{u}{2}-\cos pu}{2\sin\frac{u}{2}}\mathrm{d}u, $$

### 于是

### $$ \int_{-\pi}^{\pi}\psi\left(u\right)\frac{\cos\frac{u}{2}-\cos pu}{2\sin\frac{u}{2}}\mathrm{d}u-\frac{1}{2}\int_{0}^{\pi}\left[\psi\left(u\right)-\psi\left(-u\right)\right]\cot\frac{u}{2}\mathrm{d}u=\frac{1}{2}\int_{0}^{\pi}\left[\psi\left(u\right)-\psi\left(-u\right)\right]\frac{\cos pu}{\sin\frac{u}{2}}\mathrm{d}u. $$

### 而

### $$ \lim_{u\to0^{+}}\frac{\psi\left(u\right)-\psi\left(-u\right)}{2\sin\frac{u}{2}}=\lim_{u\to0^{+}}\frac{\psi\left(u\right)-\psi\left(0\right)-\left\lbrack\psi\left(-u\right)-\psi\left(0\right)\right\rbrack}{u}\frac{\frac{u}{2}}{\sin\frac{u}{2}}=\psi_{+}^{\prime}\left(0\right)+\psi_{-}^{\prime}\left(0\right). $$

### 利用 Riemann 引理可得

### $$ \lim_{p\to+\infty}\frac{1}{2}\int_{0}^{\pi}\left[\psi\left(u\right)-\psi\left(-u\right)\right]\frac{\cos p u}{\sin\frac{u}{2}}\mathrm{d}u=0. $$

### 3. 提示: 由于

### $$ \int_{-\delta}^{\delta}\left\{\psi\left(u\right)-\frac{1}{2}\left[\psi\left(0+\right)\right.\right.\left.\left.+\psi\left(0-\right)\right]\right\}\frac{\sin p u}{u}\mathrm{d}u=\int_{0}^{\delta}\left|\left[\psi\left(u\right)-\psi\left(0+\right)\right.\right.\left.\left.+\left[\psi\left(-u\right)-\psi\left(0-\right)\right]\right|\frac{\sin p u}{u}\mathrm{d}u, $$

### 利用 Dirichlet 引理即得结论.

### 8. $ \frac{1}{3} $

## § 3 Fourier 级数的性质

### $$ x^{2}=\frac{\pi^{2}}{3}+4\sum_{n=1}^{\infty}\frac{\left(-1\right)^{n}}{n^{2}}\cos nx,x\in\left(-\pi,\pi\right); $$

### $$ x^{3}=2\sum_{n=1}^{\infty}\frac{\left(-1\right)^{n}\left(6-\pi^{2}n^{2}\right)}{n^{3}}\sin nx,x\in\left(-\pi,\pi\right). $$

### 5. $ \sum_{n=1}^{\infty}\frac{1}{\left(2n-1\right)^{4}}=\frac{\pi^{4}}{96}. $

### 6. $ \sum_{n=1}^{\infty}\frac{1}{n^{4}}=\frac{\pi^{4}}{90}. $

### 7. 提示: 利用分部积分法可得 $ b_{n}^{\prime\prime} = -n^{2}b_{n} $. 由于

### $$ \sqrt{\left|\begin{array}{c}b_{n}\end{array}\right|}=\frac{1}{n}\sqrt{\left|\begin{array}{c}n^{2}\\ b_{n}\end{array}\right|}\leqslant\frac{1}{2}\left(\frac{1}{n^{2}}+\left|\begin{array}{c}n^{2}\\ b_{n}\end{array}\right|\right)=\frac{1}{2}\left(\frac{1}{n^{2}}+\left|\begin{array}{c}b_{n}^{n}\\ \end{array}\right|\right)\quad(n=1,2,\cdots), $$

### 所以

### $$ \sum_{n=1}^{\infty}\sqrt{\left|b_{n}\right|}\leqslant\frac{1}{2}\Big(\sum_{n=1}^{\infty}\frac{1}{n^{2}}+\sum_{n=1}^{\infty}\left|b_{n}^{n}\right|\Big)\\ =\frac{1}{2}\Big(\frac{\pi^{2}}{6}+\sum_{n=1}^{\infty}\left|b_{n}^{n}\right|\Big)\\ <\frac{1}{2}\Big(2+\sum_{n=1}^{\infty}\left|b_{n}^{n}\right|\Big). $$

### 8. 提示: 利用 Parseval 等式可知 $ \int_{-\infty}^{\pi} f^{2}(x) \, dx = 0 $, 于是 $ f(x) \equiv 0 $.

## § 4 Fourier 变换和 Fourier 积分

### $$ \frac{A}{\mathrm{i}\omega}(1-\mathrm{e}^{-\mathrm{i}\omega\delta}). $$

### (2) $ \frac{2a}{a^{2}+\omega^{2}}. $

### (3) $ \sqrt{\frac{\pi}{a}}e^{-\frac{\omega^{2}}{4a}} $. (4) $ \frac{1}{2+\mathrm{i}\omega} $.

### (5) $ A \left[ \frac{\sin(\omega - \omega_{0}) \delta}{(\omega - \omega_{0})} + \frac{\sin(\omega + \omega_{0}) \delta}{(\omega + \omega_{0})} \right] $.

### 2. 正弦变换： $ \frac{\omega}{a^{2}+\omega^{2}} $；余弦变换： $ \frac{a}{a^{2}+\omega^{2}} $

### $$ \begin{aligned}&f_{1}*f_{2}(x)=\left\{\begin{aligned}\\ &0,&x\leqslant0,\\&\frac{1}{2}(\sin x-\cos x+\mathrm{e}^{-x}),&0<x\leqslant\frac{\pi}{2},\\&\frac{1}{2}\mathrm{e}^{-x}(1+\mathrm{e}^{\frac{\pi}{2}}),&x>\frac{\pi}{2}.\\ &\end{aligned}\right.\\ \end{aligned} $$

## § 5 快速 Fourier 变换

### 1. 提示: 先将圆频率 $ \omega $ 写成频率形式 $ 2\pi s $，再对充分大的 N，在区间 $ [-N, N] $ 以间隔 $ \Delta x $ 对被积函数抽样（参见图 16.5.2），在每个小区间内利用矩形公式近似代替积分，则

### $$ \widehat{f}\left(\omega\right)\approx\int_{-N}^{N}f(x)\mathrm{e}^{-2\pi\mathrm{s}i}\mathrm{d}x\approx\sum_{n=-M}^{M}f(n\Delta x)\mathrm{e}^{-2\pi\mathrm{s}(n\Delta x)\mathrm{i}}\Delta x, $$

### 再适当代换整理,就可以得到离散 Fourier 变换形式.

### 2. 提示: 设 $ \xi \neq 1 $ 是方程 $ x^{N} = 1 $ 的一个根, 则 $ \sum_{n=0}^{N-1} \xi^{n} = 0 $.

索引

### (名词后面所标数字分别为首次出现的章和节)

### Abel 变换 9.4

### Abel 第二定理 10.3

### Abel 判别法(含参变量反常积分) 15.2

### Abel 判别法(函数项级数) 10.2

### Abel 判别法(数项级数) 9.4

### Abel 引理 9.4

### Bernstein 多项式 10.5

### Bessel 不等式 16.3

### Beta 函数 15.3

### Bolzano-Weierstrass 定理(高维) 11.1

### Cantor 闭区域套定理 11.1

### Cauchy-Hadamard 定理 10.3

### Cauchy 乘积 9.4

### Cauchy 判别法(数项级数) 9.3

### Cauchy 收敛原理(高维) 11.1

### Cauchy 收敛原理(含参变量反常积分) 15.2

### Cauchy 收敛原理(函数项级数) 10.2

### Cauchy 收敛原理(数项级数) 9.4

### Cauchy 余项 10.4

### d'Alembert 判别法(幂级数) 10.3

### d'Alembert 判别法(数项级数) 9.3

### Dini-Lipschitz 判别法 16.2

### Dini 定理(含参变量反常积分) 15.2

### Dini 定理(函数项级数) 10.2

### Dini 条件 16.2

### Dirichlet 引理 16.2

### Dirichlet-Jordan 判别法 16.2

### Dirichlet 积分 15.2

### Dirichlet 积分(Fourier 级数) 16.2

### Dirichlet 判别法(含参变量反常积分) 15.2

### Dirichlet 判别法(函数项级数) 10.2

### Dirichlet 判别法(数项级数) 9.4

### Euclid 范数 11.1

### Euclid 空间 11.1

### Euler-Fourier 公式 16.1

### FFT 16.5

### Fourier 变换 16.4

### Fourier 变换的导数 16.4

### Fourier 变换的卷积 16.4

### Fourier 积分 16.4

### Fourier 积分的三角形式 16.4

### Fourier 级数 16.1

### Fourier 级数的部分和 16.1

### Fourier 级数的复数形式 16.4

### Fourier 级数的平方收敛性质 16.3

### Fourier 级数的逐项积分定理 16.3

### Fourier 级数的逐项微分定理 16.3

### Fourier 逆变换 16.4

### Fourier 系数 16.1

### Fourier 余弦变换 16.4

### Fourier 正弦变换 16.4

### Gamma 函数 15.3

### Gauss 系数 14.1

### Green 第二公式 14.5

### Green 第一公式 14.5

### Hamilton 算子 14.5

### Heine-Borel 定理 11.1

### Hesse 矩阵 12.6

### Hölder 条件 16.2

### Jacobi 矩阵 12.1

### Jacobi 行列式 12.1 标准区域 14.3

### Jordan 曲线 14.3 部分和数列 9.1

### Lagrange 乘数法 12.7 部分积数列 9.5

### Lagrange 函数 12.7

### Laplace 方程 14.5

### Laplace 算子 14.5

### Legendre 公式 15.3

### Leibniz 级数 9.4

### Leibniz 判别法 9.4

### Lipschitz 条件 16.2

### Mobius 带 14.2

### Parseval 等式 16.3

### Peano 曲线 13.1

### Poisson 积分 13.4

### p 级数 9.1

### Raabe 判别法 9.3

### Riemann 定理 9.4

### Riemann 引理 16.2

### Schwarz 反例 14.1

### Schwarz 不等式 11.1

### Stirling 公式 15.3

### Stirling 公式（正整数） 9.5

### Taylor 公式（多元函数） 12.3

### Taylor 级数 10.4

### Taylor 系数 10.4

### Taylor 展开 10.4

### Viete 公式 9.5

### Wallice 公式 9.5

### Weierstrass 判别法（函数项级数） 10.2

### Weierstrass 第二逼近定理 16.3

### Weierstrass 第一逼近定理 10.5

### Weierstrass 判别法（含参变量反常积分） 15.2

### B 14.5

### 场 14.5

### 重极限 11.2

### D 14.2

### 单侧曲面 14.3, 14.5

### 单连通区域 11.3

### 道路 11.3

### 道路连通 16.3

### 等周问题 15.3

### 第二类 Euler 积分 14.2

### 第二类曲面积分 14.2

### 第二类曲线积分 15.3

### 第一类 Euler 积分 14.1

### 第一类曲面积分 14.1

### 第一类曲线积分 10.1

### 点态收敛 14.2

### 定向曲线 13.1

### 多重积分 10.5

### 多项式一致逼近 11.2

### 多元函数 11.2

### 多元函数组 11.2

### E 11.2

### 二重积分 13.1

### 二重积分变量代换公式 13.3

### 二重极限 11.2

### 二次极限 11.2

### 二阶偏导数 12.1

### 二阶微分 12.1

### 二维单连通区域 14.3

### 二维复连通区域 14.3

### F 13.4

### 法平面 12.5

### 法向量 12.5

### 反常重积分 13.4

### 反常二重积分 13.4

### 方向导数 12.1

### 分段单调 16.2

### 分段可导 16.4

### 复合映射(高维) 11.2

### 复连通区域 14.3

G

### 高阶偏导数 12.1

### 高阶微分 12.1

### 更序级数 9.4

### 孤立点 11.1

### 光滑曲面 12.5

### 光滑曲线 12.5

### H

### 含参变量常义积分 15.1

### 含参变量反常积分 15.2

### 函数项级数 10.1

### 函数序列 10.1

### 和函数 10.1

### 划分 13.1

### 环量 14.5

### J

### 积分次序交换定理(含参变量常义积分) 15.1

### 积分次序交换定理(含参变量反常积分) 15.2

### 积分号下求导定理(含参变量常义积分) 15.1

### 积分号下求导定理(含参变量反常积分) 15.2

### 积分路径 14.1

### 积分判别法 9.3

### 积分曲面 14.1

### 积分余项 10.4

### 积分中值定理(多元函数) 13.2

### 基本点列(高维) 11.1

### 基波 16.4

### 基频 16.4

### 级数 9.1

### 极限点 9.2

### 极值(多元函数) 12.6

### 极值点 12.6

### 几何级数 9.1

### 加法交换律 9.4

### 简单闭曲线 14.3

### 简单曲面 14.1

### 交错级数 9.4

### 紧集 11.1

### 紧集上的连续映射 11.3

### 局部性原理 16.2

### 聚点 11.1

### 卷积 16.4

### 绝对收敛 9.4

### K

### 开覆盖 11.1

### 开集 11.1

### 开区域 11.3

### 可求面积 13.1

### 可微（多元函数） 12.1

### 空间曲线的参数方程 12.5

### 快速 Fourier 变换 16.5

### L

### 累次积分 13.2

### 累次极限 11.2

### 离散 Fourier 变换 16.5

### 离散 Fourier 逆变换 16.5

### 连通集 11.3

### 连续函数 11.2

### 连续性定理（含参变量常义积分） 15.1

### 连续性定理（含参变量反常积分） 15.2

### 连续映射（高维） 11.2

### 链式法则 12.2

### 零边界区域 13.1

### M

### 幂级数 10.3

### 幂级数展开 10.4

### 面积 13.1

### 面积元素 13.1

### 目标函数 12.7

### N

### 内闭一致收敛 10.2

### 内点 11.1

### 内积 11.1

### 拟合曲线 12.6

### 逆映射定理 12.4

### P

### 偏导函数

### 偏导数

### Q

### 奇点

### 切平面

### 切向量

### 曲顶柱体

### 曲线积分与路径无关

### 曲线坐标

### 曲线坐标网

### 全微分

### R

### 任意项级数

### S

### 三角级数

### 散度

### 散度场

### 上极限

### 势函数

### 收敛半径

### 收敛点

### 收敛因子

### 收敛域（函数项级数）

### 收敛域（含参变量反常积分）

### 数量场

### 数列卷积

### 数项级数

### 双侧曲面

### T

### 梯度

### 梯度场

### 体积元素

### 调和函数

### 调和级数

### 条件极值

### 条件收敛

### 通量

### 正项级数 9.3 驻点 12.6

### 正弦级数 16.1 最佳平方逼近元素 16.3

### 中间值定理(高维) 11.3 最小二乘法 12.6

### 逐项求导(函数项级数) 10.1 最值(多元函数) 12.6

### 逐项求积分(函数项级数) 10.1 最值定理(高维) 11.3

### 逐项求极限(函数项级数) 10.1 坐标变换 13.3

郑重声明

### 高等教育出版社依法对本书享有专有出版权。任何未经许可的复制、销售行为均违反《中华人民共和国著作权法》，其行为人将承担相应的民事责任和行政责任；构成犯罪的，将被依法追究刑事责任。为了维护市场秩序，保护读者的合法权益，避免读者误用盗版书造成不良后果，我社将配合行政执法部门和司法机关对违法犯罪的单位和个人进行严厉打击。社会各界人士如发现上述侵权行为，希望及时举报，本社将奖励举报有功人员。

### 反盗版举报电话（010）58581999 58582371 58582488

### 反盗版举报传真（010）82086060

### 反盗版举报邮箱 dd@hep.com.cn

### 通信地址 北京市西城区德外大街4号

### 高等教育出版社法律事务与版权管理部

### 邮政编码 100120

### 防伪查询说明

### 用户购书后刮开封底防伪涂层，利用手机微信等软件扫描二维码，会跳转至防伪查询网页，获得所购图书详细信息。用户也可将防伪二维码下的20位密码按从左到右、从上到下的顺序发送短信至106695881280，免费查询所购图书真伪。

### 反盗版短信举报

### 编辑短信“JB,图书名称,出版社,购买地点”发送至10669588128

### 防伪客服电话

### (010) 58582300

### <div style="text-align: center;"><img src="https://pplines-online.bj.bcebos.com/deploy/official/paddleocr/pp-ocr-vl-15//5720ee92-a7c9-4f76-b719-73f5b6ac16f7/markdown_1/imgs/img_in_image_box_91_1258_232_1393.jpg?authorization=bce-auth-v1%2FALTAKzReLNvew3ySINYJ0fuAMN%2F2026-03-30T13%3A31%3A00Z%2F-1%2F%2Fe6b80b536757c7fcada38f5aa6f5fed4801ce39ad1e0c9069498551bc53f410a" alt="Image" width="13%" /></div>

### 数字课程网站

### 网址：http://abook.hep.com.cn/1257652

### http://abook.hep.edu.cn/1257652

### 数学课程编号 使用说明 常见书内数学课程说明

### <div style="text-align: center;"><img src="https://pplines-online.bj.bcebos.com/deploy/official/paddleocr/pp-ocr-vl-15//5720ee92-a7c9-4f76-b719-73f5b6ac16f7/markdown_1/imgs/img_in_image_box_753_1222_923_1347.jpg?authorization=bce-auth-v1%2FALTAKzReLNvew3ySINYJ0fuAMN%2F2026-03-30T13%3A31%3A00Z%2F-1%2F%2Fab012081968505666db2f6ec4babf80531008f6d790c0da6bb103cb50fba16aa" alt="Image" width="16%" /></div>

### 9 "787040"516302">

### 定价 54.00 元

