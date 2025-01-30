---
exports:
  - format: pdf
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
numbering:
  code: false
  equation: true
  headings: true
---

# Fourier transform

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

We are interested in Fourier transform of a function defined on a bounded interval and we can think of this function being extended periodically to $\re$.

For results on convergence of Fourier series, see

1. [@Tveito2005], Chapter 8 and 9
1. [@Davis1963], Chapter 12
1. [@Shen2011]

## Fourier series

Let $f : [-\pi,\pi] \to \re$ be a function which is extended periodically to the whole real line. The Fourier series is defined as

$$
Sf(x) = a_0 + \sum_{k=1}^\infty [ a_k \cos(k x) + b_k \sin(k x) ]
$$

where

$$
a_0 &= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ud x \\
a_k &= \frac{1}{\pi} \int_{-\pi}^\pi f(x) \cos(k x) \ud x, \qquad k=1, 2, \ldots \\
b_k &= \frac{1}{\pi} \int_{-\pi}^\pi f(x) \sin(k x) \ud x, \qquad k=1, 2, \ldots
$$

The Fourier series can be written as

$$
Sf(x) = \sum_{k=-\infty}^\infty \hat f_k \ee^{\ii k x}, \qquad \hat f_k = \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x
$$

where

$$
\hat f_0 = a_0, \qquad \hat f_k = \begin{cases}
\half(a_k + \ii b_k), & k = -1, -2, \ldots \\
\half(a_k - \ii b_k), & k = +1, +2, \ldots
\end{cases}
$$

Since $f(x)$ is real, $a_k, b_k$ are real and we have

$$
\hat f_{-k} = \conj{\hat f_k}, \qquad k=1,2,\ldots
$$

:::{prf:remark}
We can take any interval of length $2\pi$, e.g., $[0,2\pi]$. Some books use the interval $[-1,1]$ with period of 2 to define the Fourier series.
:::

## Convergence of Fourier series

When is $f = Sf$ ? To talk of such questions, let us define the truncated series

$$
S_n f(x) = a_0 + \sum_{k=1}^n [ a_k \cos(k x) + b_k \sin(k x) ] = \sum_{k=-n}^n \hat f_k \ee^{\ii k x}
$$

There are three types of convergence wrt $n \to \infty$.

1. pointwise convergence
    $$
    S_n f(x) \to f(x), \qquad \forall x
    $$

1. uniform convergence
    $$
    \norm{S_n f - f}_\infty \to 0
    $$

1. mean square convergence
    $$
    \norm{S_n f - f}_2 \to 0
    $$

where the norms are defined as

$$
\norm{f}_\infty = \max_{x} |f(x)|, \qquad \norm{f}_2 = \left( \int_{-\pi}^\pi |f(x)|^2 \ud x \right)^\half
$$

The convergence depends on the smoothness of the function, which in turn controls the rate of decay of $\hat f_k$. Clearly we need $|\hat f_k| \to 0$ as $|k| \to \infty$ for the series to converge and the rate of decay needs to sufficiently fast.

:::{prf:definition}
1. The function $f : [-\pi,\pi] \to \re$ is said to be piecewise continuous if it is continuous at all but a finite number of points where both left and right limits exist.

1. The function $f : [-\pi, \pi] \to \re$ is said to be piecewise differentiable if it is differentiable everywhere except for a finite number of points where the one-sided derivatives both exist.

1. We say that $f \in \cts_p^0[-\pi,\pi]$ if $f : (-\pi,\pi) \to \re$ is continuous and $f(-\pi) = f(\pi)$; i.e., the periodic extension to $\re$ is continuous.

1. We say that $f \in \cts_p^\nu[-\pi,\pi]$ if $f : (-\pi,\pi) \to \re$ is $\nu$-times continuously differentiable and 

    $$
    f^{(s)}(-\pi) = f^{(s)}(\pi), \qquad s = 0,1,\ldots,\nu
    $$ 

    i.e., the periodic extension to $\re$ is in $\cts^\nu(\re)$.
:::

:::{prf:theorem}
:label: thm:fdecay
(1) If $f : [-\pi, \pi] \to \re$ is piecewise continuous then

$$
|\hat f_k| = \order{\frac{1}{|k|}}, \qquad |k| \to \infty
$$

(2) If $f \in \cts_p^{s-1}[-\pi,\pi]$ and $f^{(s)}$ is piecewise continuous, then
$$
|\hat f_k| = \order{\frac{1}{|k|^{s+1}}}, \qquad |k| \to \infty
$$
:::

:::{prf:proof}
See [@Gander2018], Theorem 4.2
:::

:::{prf:example} Piecewise continuous function 

$$
u(x) = \begin{cases}
1, & \half\pi < x \le \frac{3}{2}\pi \\
0, & 0 < x \le \half \pi, \quad \frac{3}{2}\pi < x < 2\pi
\end{cases}
$$ 

$u$ is of bounded variation in $[0,2\pi]$. Its Fourier
coefficients are 

$$
\hatu_k = \begin{cases}
\pi, & k = 0 \\
0, & k \textrm{ even} \\
\frac{(-1)^{(k-1)/2}}{k}, & k \textrm{ odd}
\end{cases}
$$ 

This corresponds to Case (1) of [](#thm:fdecay).
:::

:::{prf:example} Continuous function

$$
u(x) = \begin{cases}
0, & x \le \half\pi, x \ge \frac{3}{2}\pi \\
2x/\pi - 1, & \half\pi \le x \le \pi \\
3 - 2x/\pi, & \pi \le x \le \frac{3}{2}\pi
\end{cases}
$$ 

$$
|\hatu_k| = \order{ \frac{1}{|k|^2} }
$$

This corresponds to Case (1) of [](#thm:fdecay) with $s=1$.
:::

:::{prf:example} Infinitely differentiable but not periodic 

$$
u(x) = \sin(x/2), \qquad x \in [0,2\pi]
$$ 

$u(0^+) = u((2\pi)^{-1})$ but $u'(0^+) \ne u'((2\pi)^-)$. Its Fourier coefficients are

$$
\hatu_k = \frac{2}{\pi(1 - 4 k^2)}
$$ 

This corresponds to Case (1) of [](#thm:fdecay) with $s=1$.
:::

:::{prf:example} Infinitely differentiable and periodic

$$
u(x) = \frac{3}{5 - 4\cos x}, \qquad x \in [0,2\pi]
$$ 

Its Fourier coefficients are

$$
\hatu_k = \frac{1}{2^{|k|}}
$$
:::

:::{prf:theorem}
(1) If $f : [-\pi,\pi] \to \re$ is piecewise continuous, then 

$$
\lim_{n \to \infty} \norm{S_n f - f}_2 = 0
$$

Moreover, Parseval's relation holds.

(2) If $f \in \cts_p^0[-\pi,\pi]$ and $g = f'$ is piecewise continuous, then 

$$
\lim_{n \to \infty} \norm{S_n f - f}_\infty = 0
$$

Moreover, $\hat g_k = \ii k \hat f_k$.
:::

:::{prf:proof}
(1) [@Tveito2005], Theorem 9.4 and Corollary 9.1

(2) [@Tveito2005], Theorem 9.3 and Theorem 8.1
:::
