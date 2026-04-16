# 第十一章 Euclid 空间上的极限和连续

## § 1 Euclid 空间上的基本定理

 $$ S^{o}=\left\{\left(x,y\right)\mid x>0,y\neq0\right\};\partial S=\left\{\left(x,y\right)\mid x=0 或 x>0,y=0\right\};\bar{S}=\left\{\left(x,y\right)\mid x\geqslant0\right\}. $$ 

 $$ \left(2\right)S^{o}=\left|\left(x,y\right)\right|0<x^{2}+y^{2}<1\left|\right;\partial S=\left|\left(x,y\right)\right|x^{2}+y^{2}=0 或 x^{2}+y^{2}=1\left|\right;\widetilde{S}=\left|\left(x,y\right)\right|x^{2}+y^{2}\leqslant1\left|\right]. $$ 

 $$ S^{\prime \prime}=\varnothing;\partial S=\left\{\left(x,y\right)\mid0<x\leqslant1,y=\sin\frac{1}{x} 或 x=0,-1\leqslant y\leqslant1\right\}; $$ 

 $$ \bar{S}=\left|\left(x,y\right)\right|0<x\leqslant1,y=\sin\frac{1}{x} 或 x=0,-1\leqslant y\leqslant1\left|\right|. $$ 

5. (1)  $ S' = \left| \pm 1 \right\rangle $.

(2)  $ S' = \varnothing $.

(3)  $ S' = \left\{(x, y) \mid y^2 - x^2 + 1 \leqslant 0\right\} $

## § 2 多元连续函数

1. (1)  $ D = \{(x, y) \mid x^{2} + y^{2} < 1, y > x\} $.

(2)  $ D = \{(x, y, z) \mid x > 0, y > 0, z > 0\} $.

(3)  $ D = \{(x, y, z) \mid r^2 \leqslant x^2 + y^2 + z^2 \leqslant R^2\} $.

(4)  $ D = \left\{(x, y, z) \mid |z| \leqslant x^2 + y^2, x^2 + y^2 \neq 0\right\} $.

2.  $ f(x)=\frac{1}{\left(1+x^{2}\right)^{\frac{3}{2}}} $

3.  $ f(x)=x^{2}+2x, z(x,y)=x+\sqrt{y}-1. $

4.（1）不存在.(2）不存在.(3）不存在.

 $$ \frac{x^{4}+y^{8}}{3}=\frac{\frac{1}{2}x^{4}+\frac{1}{2}x^{4}+y^{8}}{3}\geqslant\sqrt[3]{\frac{1}{4}x^{8}y^{8}}. $$ 

7.（1）1.（2）+ $ \infty $.（3） $ \frac{1}{2} $.（4）2.（5）1.（6）0.（7）+ $ \infty $.（8）0.

8.（1）两个二次极限存在为0，二重极限不存在.

(2) 两个二次极限存在分别为 1 和 -1，二重极限不存在.

(3) 两个二次极限不存在，二重极限存在为 0.

11. 提示: 利用 Lagrange 中值定理  $ f(x) - f(y) = f'(\xi)(x - y) $.

12. 提示: 利用  $ \left|f(x,y)-f(x_{0},y_{0})\right|\leqslant\left|f(x,y)-f(x,y_{0})\right|+\left|f(x,y_{0})-f(x_{0},y_{0})\right| $.

## § 3 连续函数的性质

3. 提示:  $ f\left(1-\frac{1}{2n},1-\frac{1}{2n}\right)-f\left(1-\frac{1}{n},1-\frac{1}{n}\right)\rightarrow+\infty $

5.（1）提示：任取一点 $ (x_{0},y_{0}) $，由 $ \lim_{x^{2}+y^{2}\to+\infty}f(x,y)=+\infty $，可知存在R>0，当 $ x^{2}+y^{2}>R^{2} $，成立 $ f(x,y)>f(x_{0},y_{0}) $。 $ f(x,y) $在紧集 $ \{(x,y)\mid x^{2}+y^{2}\leq R^{2}\} $上必定取到最小值，且此最小值就是它在 $ R^{2} $上的最小值。

(2) 提示: 任取  $ (x_{0}, y_{0}) $, 设  $ f(x_{0}, y_{0}) > 0 $, 由  $ \lim_{x^{2} + y^{2} \to +\infty} f(x, y) = 0 $, 可知存在 R > 0, 当  $ x^{2} + y^{2} > R^{2} $, 成立  $ f(x, y) < f(x_{0}, y_{0}) $, 则  $ f(x, y) $ 在紧集  $ \{(x, y) \mid x^{2} + y^{2} \leq R^{2}\} $ 上必定取到最大值, 且此最大值就是它在  $ \mathbb{R}^{2} $ 上的最大值; 若  $ f(x_{0}, y_{0}) < 0 $, 由  $ \lim_{x^{2} + y^{2} \to +\infty} f(x, y) = 0 $, 可知存在 R > 0, 当  $ x^{2} + y^{2} > R^{2} $, 成立  $ f(x, y) > f(x_{0}, y_{0}) $, 则  $ f(x, y) $ 在紧集  $ \{(x, y) \mid x^{2} + y^{2} \leq R^{2}\} $ 上必定取到最小值, 且此最小值就是它在  $ \mathbb{R}^{2} $ 上的最小值.

6. 提示: 单位球面是  $ R^{n} $ 上的紧集, 设 f 在单位球面上的最小、最大值分别为 a 和 b, 再利用  $ f(x)=\|x\|f\left(\frac{x}{\|x\|}\right) $.

8. 提示: 设  $ \zeta \in \partial D $，证明对任意点列  $ \{x_n\} (x_n \in D, x_n \to \zeta) $，点列  $ \{f(x_n)\} $ 收敛，且极限只与  $ \zeta $ 有关，而与点列  $ \{x_n\} $ 的选取无关，记该极限为  $ g(\zeta) $，令

 $$ \widetilde{f}(x)=\left\{\begin{aligned}&f(x),&x\in D,\\ &g(x),&x\in\partial D,\end{aligned}\right. $$ 

再证明  $ \tilde{f} $ 在  $ \bar{D} $ 连续.

