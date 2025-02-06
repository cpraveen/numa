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

Periodic functions can be approximated by a series of sine and cosine functions.  This leads us to the Fourier transform of a function defined on a bounded interval and we can think of this function being extended periodically to $\re$.

+++

:::{seealso}
For results on convergence of Fourier series, see

1. [@Tveito2005], Chapter 8 and 9, uses elementary analysis techniques. It is recommended to read these two chapters to understand the following results.
1. [@Davis1963], Chapter 12.
1. [@Shen2011], uses Sobolev spaces.
:::

+++

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

Due to periodicity $k$ only takes integer values.

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
We can take any interval of length $2\pi$, e.g., $[0,2\pi]$. Some books like [@Tveito2005] use the interval $[-1,1]$ with period of 2 to define the Fourier series.
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

Uniform convergence is the strongest property, and it implies pointwise and mean square convergence [@Tveito2005].

Clearly we need $|\hat f_k| \to 0$ as $|k| \to \infty$ for the series to converge and the rate of decay needs to be sufficiently fast. The rate of decay of $\hat f_k$ depends on the smoothness of the function.

+++

:::{prf:definition}
1. The function $f : [-\pi,\pi] \to \re$ is said to be **piecewise continuous** if it is continuous at all but a finite number of points where both left and right limits exist.

1. The function $f : [-\pi, \pi] \to \re$ is said to be **piecewise differentiable** if it is differentiable everywhere except for a finite number of points where the one-sided derivatives exist.

1. We say that $f \in \cts_p^0[-\pi,\pi]$ if $f : (-\pi,\pi) \to \re$ is continuous and $f(-\pi) = f(\pi)$; i.e., the periodic extension to $\re$ is continuous.

1. We say that $f \in \cts_p^\nu[-\pi,\pi]$ if $f : (-\pi,\pi) \to \re$ is $\nu$-times continuously differentiable and 

    $$
    f^{(s)}(-\pi) = f^{(s)}(\pi), \qquad s = 0,1,\ldots,\nu
    $$ 

    i.e., the periodic extension to $\re$ is in $\cts^\nu(\re)$.
:::

+++

:::{prf:theorem} Decay of Fourier coefficients
:label: thm:fdecay
(1) If $f : [-\pi, \pi] \to \re$ is piecewise continuous then

$$
|\hat f_k| = \order{\frac{1}{|k|}}, \qquad |k| \to \infty
$$

(2) If $f \in \cts_p^{s-1}[-\pi,\pi]$ and $f^{(s)}$ is piecewise continuous for some $s \ge 1$, then
$$
|\hat f_k| = \order{\frac{1}{|k|^{s+1}}}, \qquad |k| \to \infty
$$
:::

The term *piecewise continuous* can be replaced with *bounded variation*.  The proof is based on [@Gander2018], Theorem 4.2.

:::{prf:proof}

(1a) First consider $f(x)$ to be a step function

$$
f(x) = y_j, \qquad x \in (x_{j-1},x_j), \qquad 1 \le j \le m
$$

for some partition $\{x_0, x_1, \ldots, x_m\}$. Then

\begin{align}
\hat f_k 
&= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x \\
&= \frac{1}{2\pi \ii k} \left[ y_1 + \ee^{-\ii k x_1} (y_2 - y_1) + \ee^{-\ii k x_2} (y_3 - y_2) + \right. \\
&  \qquad\qquad \left. \ldots + \ee^{-\ii k x_{n-1}} (y_m - y_{m-1}) - y_m \right]
\end{align}

and hence

$$
|\hat f_k| \le \frac{\TV(f) + |f(\pi) - f(-\pi)|}{2 \pi |k|} = \order{ \frac{1}{|k|} }
$$

(1b) Now let $f(x)$ be a general piecewise continuous function; it is Riemann integrable and for any $\epsilon > 0$, we can find a partition $\{ x_0, x_1, \ldots, x_m \}$ such that

$$
\left| \int_{-\pi}^\pi f(x) \ud x - \sum_{j=1}^n f(\xi_j) (x_j - x_{j-1}) \right| \le \epsilon, \qquad \xi_j \in [x_{j-1}, x_j]
$$

and we can choose $\xi_1 = x_0 = -\pi$, $\xi_m = x_m = \pi$. Consider the step function

$$
\tilde f(x) = f(\xi_j), \qquad x \in [x_{j-1}, x_j]
$$

for which

$$
\TV(\tilde f) = \sum_{j=2}^m |f(\xi_j) - f(\xi_{j-1})| \le \TV(f)
$$

Then the Fourier transform

\begin{align}
\hat f_k 
&= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x \\
&\approx \frac{1}{2\pi} \int_{-\pi}^\pi \tilde f(x) \ee^{-\ii k x} \ud x \\
&\le \frac{\TV(\tilde f) + |\tilde f(\pi) - \tilde f(-\pi)|}{2 \pi |k|} \\
&\le \frac{\TV(f) + |f(\pi) - f(-\pi)|}{2 \pi |k|}
\end{align}

(2) If $f \in \cts^0_p[-\pi,\pi]$ and piecewise differentiable, i.e., $s=1$, then using integration by parts

\begin{align}
\hat f_k
&= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x \\
&= -\frac{1}{2 \pi \ii k} \left[ f(x) \ee^{-\ii k x} \right]_{-\pi}^\pi + \frac{1}{k} \left[ \frac{1}{2 \pi \ii} \int_{-\pi}^\pi f'(x) \ee^{-\ii k x} \ud x \right] \\
& \qquad \textrm{Using the result of Part (1) in the second term above} \\
&= -\frac{\cos(k\pi)}{2 \pi \ii k} [f(\pi) - f(-\pi)] + \frac{1}{k} \order{\frac{1}{k}} \\
&= \order{\frac{1}{k^2}}
\end{align}

For any $s > 1$, we can perform several integration by parts and reach the conclusion of Part (2).
:::

:::{prf:example} Piecewise continuous function 

$$
u(x) = \begin{cases}
1, & \half\pi < x \le \frac{3}{2}\pi \\
0, & 0 < x \le \half \pi, \quad \frac{3}{2}\pi < x < 2\pi
\end{cases}
$$ 

```{code-cell}
:tags: remove-input
u = lambda x: (x > 0.5*pi)*(x < 1.5*pi)*(1) + 0.0
x = linspace(0,2*pi,500)
plot(x,u(x))
grid(True), xlabel('x'), ylabel('u(x)');
```

Its Fourier coefficients are 

$$
\hatu_k = \begin{cases}
\pi, & k = 0 \\
0, & k \textrm{ even} \\
\frac{(-1)^{(k-1)/2}}{k}, & k \textrm{ odd}
\end{cases}, 
\qquad
|\hat u_k| = \order{\frac{1}{|k|}}
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

```{code-cell}
:tags: remove-input
u = lambda x: (x > 0.5*pi)*(x < pi)*(2*x/pi - 1) \
              + (x > pi)*(x < 1.5*pi)*(3-2*x/pi) + 0.0
x = linspace(0,2*pi,500)
plot(x,u(x))
grid(True), xlabel('x'), ylabel('u(x)');
```

$$
|\hatu_k| = \order{ \frac{1}{|k|^2} }
$$

This corresponds to Case (2) of [](#thm:fdecay) with $s=1$.
:::

:::{prf:example} Infinitely differentiable but not periodic 

$$
u(x) = \sin(x/2), \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
u = lambda x: sin(x/2)
v = lambda x: cos(x/2)/2
x = linspace(0,2*pi,500)
plot(x,u(x),label='u(x)')
plot(x,v(x),label='u\'(x)')
legend(), grid(True), xlabel('x');
```

Its Fourier coefficients are

$$
\hatu_k = \frac{2}{\pi(1 - 4 k^2)}, \qquad |\hat u_k| = \order{\frac{1}{k^2}}, \qquad |k| \to \infty
$$ 

$u(0) = u(2\pi)$ but 

$$
u'(0) = \half \cos(0) = \half \qquad\ne\qquad u'(2\pi) = \half \cos(\pi) = -\half
$$

hence  $u \in \cts^0_p[0,2\pi]$ but $u \notin \cts^1_p[0,2\pi]$. This corresponds to Case (1) of [](#thm:fdecay) with $s=1$.
:::

:::{prf:example} Infinitely differentiable and periodic

$$
u(x) = \frac{3}{5 - 4\cos x}, \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
u = lambda x: 3.0/(5.0 - 4.0*cos(x))
x = linspace(0,2*pi,500)
plot(x,u(x))
grid(True), xlabel('x'), ylabel('u(x)');
```

Its Fourier coefficients are

$$
\hatu_k = \frac{1}{2^{|k|}}
$$

$u \in \cts^s_p[0,2\pi]$ for every $s \ge 0$, and the Fourier coefficients decay faster than any power of $1/|k|$.
:::

:::{prf:theorem} Convergence of Fourier series
(1) If $f : [-\pi,\pi] \to \re$ is piecewise continuous, then 

$$
\lim_{n \to \infty} \norm{S_n f - f}_2 = 0
$$

Moreover, Parseval's relation holds

$$
\half a_0^2 + \sum_{k=1}^\infty (a_k^2 + b_k^2) = \sum_{k=-\infty}^\infty |\hat f_k|^2 = \norm{f}_2^2
$$

(2) If $f \in \cts_p^0[-\pi,\pi]$ and $g = f'$ is piecewise continuous, then the Fourier series converges absolutely and uniformly,

$$
\lim_{n \to \infty} \norm{S_n f - f}_\infty = 0
$$

Moreover, $\hat g_k = \ii k \hat f_k$.
:::

:::{prf:proof}
(1) [@Tveito2005], Theorem 9.4 and Corollary 9.1

(2) [@Tveito2005], Theorem 9.3 and Theorem 8.1 or [@Davis1963], Theorem 12.1.5.
:::
