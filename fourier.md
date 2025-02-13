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

If $f(x)$ is real, $a_k, b_k$ are real and we have

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

Note that we have an expansion in terms of an orthogonal basis $\phi_k(x) = \ee^{\ii k x}$, since

$$
\int_{-\pi}^\pi \phi_j(x) \phi_k(x) \ud x = 2 \pi \delta_{jk}
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
\norm{f}_\infty = \max_{x \in [-\pi,\pi]} |f(x)|, \qquad \norm{f}_2 = \left( \int_{-\pi}^\pi |f(x)|^2 \ud x \right)^\half
$$

Uniform convergence is the strongest property, and it implies pointwise and mean square convergence [@Tveito2005].

Clearly we need $|\hat f_k| \to 0$ as $|k| \to \infty$ for the series to converge and the rate of decay needs to be sufficiently fast. The rate of decay of $\hat f_k$ depends on the smoothness of the function.

+++

:::{prf:definition}
1. The function $f : [-\pi,\pi] \to \re$ is said to be **piecewise continuous** if it is continuous at all but a finite number of points where both left and right limits exist.

1. The function $f : [-\pi, \pi] \to \re$ is said to be **piecewise differentiable** if it is differentiable everywhere except for a finite number of points where the one-sided derivatives exist.

1. We say that $f \in \cts_p^0[-\pi,\pi]$ if it is continuous and $f(-\pi) = f(\pi)$; i.e., the periodic extension of $f$ to $\re$ is continuous.

1. We say that $f \in \cts_p^\nu[-\pi,\pi]$ if it is $\nu$-times continuously differentiable and 

    $$
    f^{(s)}(-\pi) = f^{(s)}(\pi), \qquad s = 0,1,\ldots,\nu
    $$ 

    i.e., the periodic extension of $f$ to $\re$ is in $\cts^\nu(\re)$.
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

The term *piecewise continuous* can be replaced with *bounded variation*.  The proof is based on an adaptation of [@Gander2018], Theorem 4.2.

:::{prf:proof}

(1a) First consider $f(x)$ to be a step function

$$
f(x) = y_j, \qquad x \in (x_{j-1},x_j), \qquad 1 \le j \le m
$$

for some partition $\{x_0, x_1, \ldots, x_m\}$ of $[-\pi,\pi]$. Then

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

where

$$
\TV(f) = \sum_{j=2}^m |y_j - y_{j-1}| < \infty
$$

(1b) Now let $f(x)$ be a general piecewise continuous function; it is Riemann integrable and for any $\epsilon > 0$, there is a partition $P = \{ x_0, x_1, \ldots, x_m \}$ such that [@Rudin1976]

$$
U(f,P) - L(f,P) = \sum_{j=1}^m M_j (x_j - x_{j-1}) - \sum_{j=1}^m m_j (x_j - x_{j-1}) \le \epsilon
$$

where

$$
m_j = \min_{x \in [x_{j-1},x_j]} f(x) = f(\xi_j), \qquad
M_j = \max_{x \in [x_{j-1},x_j]} f(x)
$$

Consider the step function

$$
\tilde f(x) = m_j, \qquad x \in (x_{j-1}, x_j)
$$

for which

$$
\TV(\tilde f) = \sum_{j=2}^m |f(\xi_j) - f(\xi_{j-1})| \le \TV(f)
$$

Now [@Rudin1976]

\begin{align}
\left| \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x - \int_{-\pi}^\pi \tilde f(x) \ee^{-\ii k x} \ud x \right| 
&\le \int_{-\pi}^\pi |f(x) - \tilde f(x)| \ud x \\
&= \sum_{j=1}^m \int_{x_{j-1}}^{x_j} |f(x) - \tilde f(x)| \ud x \\
&= \sum_{j=1}^m \int_{x_{j-1}}^{x_j} (f(x) - m_j) \ud x \\
&\le  \sum_{j=1}^m \int_{x_{j-1}}^{x_j} (M_j - m_j) \ud x \\
&= U(f,P) - L(f,P) \\
&\le \epsilon
\end{align}

Then the Fourier transform

\begin{align}
\hat f_k 
&= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x \\
&\approx \frac{1}{2\pi} \int_{-\pi}^\pi \tilde f(x) \ee^{-\ii k x} \ud x \\
&\le \frac{\TV(\tilde f) + |\tilde f(\pi) - \tilde f(-\pi)|}{2 \pi |k|} \qquad \textrm{using Part (1a)} \\
&\lessapprox \frac{\TV(f) + |f(\pi) - f(-\pi)|}{2 \pi |k|} \\
&= \order{\frac{1}{|k|}}
\end{align}

(2) If $f \in \cts^0_p[-\pi,\pi]$ and piecewise differentiable, i.e., $s=1$, then using integration by parts

\begin{align}
\hat f_k
&= \frac{1}{2\pi} \int_{-\pi}^\pi f(x) \ee^{-\ii k x} \ud x \\
&= -\frac{1}{2 \pi \ii k} \left[ f(x) \ee^{-\ii k x} \right]_{-\pi}^\pi + \frac{1}{\ii k} \left[ \frac{1}{2 \pi} \int_{-\pi}^\pi f'(x) \ee^{-\ii k x} \ud x \right] \\
& \qquad \textrm{Using the result of Part (1) in the second term above} \\
&= -\frac{\cos(k\pi)}{2 \pi \ii k} [f(\pi) - f(-\pi)] + \frac{1}{\ii k} \order{\frac{1}{k}} \\
&= \order{\frac{1}{k^2}}, \qquad \textrm{since } f(-\pi) = f(\pi)
\end{align}

We see that every time we do an integration by parts, we get a factor of $1/k$.  For any $s > 1$, we can perform several integration by parts and reach the conclusion of Part (2).
:::

:::{prf:example} Piecewise continuous function 
:label: ex:sqhat

$$
f(x) = \begin{cases}
1, & \half\pi < x \le \frac{3}{2}\pi \\
0, & \textrm{otherwise}
\end{cases}, \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
f = lambda x: (x > 0.5*pi)*(x < 1.5*pi)*(1)
x = linspace(0,2*pi,500)
plot(x,f(x))
grid(True), xlabel('x'), ylabel('f(x)');
```

Its Fourier coefficients are 

$$
\hat f_k = \begin{cases}
\half, & k = 0 \\
0, & k \textrm{ even} \\
\frac{(-1)^{(k+1)/2}}{\pi k}, & k \textrm{ odd}
\end{cases}, 
\qquad
|\hat f_k| = \order{\frac{1}{|k|}}
$$ 

This corresponds to Part (1) of [](#thm:fdecay).
:::

:::{prf:example} Continuous function
:label: ex:trihat

$$
f(x) = \begin{cases}
2x/\pi - 1, & \half\pi \le x \le \pi \\
3 - 2x/\pi, & \pi \le x \le \frac{3}{2}\pi \\
0, & \textrm{otherwise}
\end{cases}, \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
f = lambda x: (x > 0.5*pi)*(x < pi)*(2*x/pi - 1) \
              + (x > pi)*(x < 1.5*pi)*(3-2*x/pi)
x = linspace(0,2*pi,500)
plot(x,f(x))
grid(True), xlabel('x'), ylabel('f(x)');
```

$$
|\hat f_k| = \order{ \frac{1}{|k|^2} }
$$

This corresponds to Part (2) of [](#thm:fdecay) with $s=1$.
:::

+++

:::{prf:example} Infinitely differentiable but not periodic 

$$
f(x) = \sin(x/2), \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
f = lambda x: sin(x/2)
df= lambda x: cos(x/2)/2
x = linspace(0,2*pi,500)
plot(x,f(x),label='f(x)')
plot(x,df(x),label='f\'(x)')
legend(), grid(True), xlabel('x');
```

Its Fourier coefficients are

$$
\hat f_k = \frac{2}{\pi(1 - 4 k^2)}, \qquad |\hat f_k| = \order{\frac{1}{k^2}}, \qquad |k| \to \infty
$$ 

Now $f(0) = f(2\pi)$ but 

$$
f'(0) = \half \cos(0) = \half \qquad\ne\qquad f'(2\pi) = \half \cos(\pi) = -\half
$$

and hence  $f \in \cts^0_p[0,2\pi]$ but $f \notin \cts^1_p[0,2\pi]$. This corresponds to Part (1) of [](#thm:fdecay) with $s=1$. 

We can see this with integration-by-parts;

\begin{align}
\hat f_k
&= \frac{1}{2\pi} \int_{0}^{2\pi} f(x) \ee^{-\ii k x} \ud x \\
&= \frac{1}{2 \pi \ii k} \int_{0}^{2\pi} f'(x) \ee^{-\ii k x} \ud x \qquad \textrm{since } f(0) = f(2\pi) \\
&= \order{\frac{1}{k^2}} [ f'(2\pi) - f'(0)] + \frac{1}{2 \pi (\ii k)^2} \int_{0}^{2\pi} f''(x) \ee^{-\ii k x} \ud x
\end{align}

Though we can perform more integrations giving rise to terms of order $1/|k|^3, 1/|k|^4, \ldots$, 

$$
|\hat f_k| = \order{\frac{1}{k^2}} [ f'(2\pi) - f'(0)] + \order{\frac{1}{|k|^3}} + \order{\frac{1}{|k|^4}} + \ldots
$$

the size of $|\hat f_k|$ is $\order{1/k^2}$ for $|k| \to \infty$, since $f'(0) \ne f'(2\pi)$.
:::

+++

:::{prf:remark}

We can use Sympy to compute the Fourier transform. This code is for the above function.

```{code-cell}
def compute_fk():
    from sympy import pi, sin, exp, Symbol, integrate, pprint, simplify
    k = Symbol('k', integer=True)
    x = Symbol('x', real=True)
    fk = 1/(2*pi) * integrate(sin(x/2)*exp(-1j*k*x),(x,0,2*pi))
    pprint(simplify(fk))

compute_fk()
```

:::

+++

:::{prf:example} Infinitely differentiable and periodic

$$
f(x) = \frac{3}{5 - 4\cos x}, \qquad x \in [0,2\pi]
$$ 

```{code-cell}
:tags: remove-input
f = lambda x: 3.0/(5.0 - 4.0*cos(x))
x = linspace(0,2*pi,500)
plot(x,f(x))
grid(True), xlabel('x'), ylabel('f(x)');
```

Its Fourier coefficients are

$$
\hat f_k = \frac{1}{2^{|k|}}
$$

$f \in \cts^s_p[0,2\pi]$ for every $s \ge 0$, and the Fourier coefficients decay faster than any power of $1/|k|$, which is consistent with exponential decay.
:::

+++

:::{prf:theorem} Convergence of Fourier series
:label: thm:fserconv
(1) If $f : [-\pi,\pi] \to \re$ is piecewise continuous, then 

$$
\lim_{n \to \infty} \norm{S_n f - f}_2 = 0
$$

Moreover, Parseval's relation holds

$$
\sum_{k=-\infty}^\infty |\hat f_k|^2 = \norm{f}_2^2
$$

(2) If $f \in \cts_p^0[-\pi,\pi]$ and $f'$ is piecewise continuous, then the Fourier series converges absolutely and uniformly,

$$
\lim_{n \to \infty} \norm{S_n f - f}_\infty = 0
$$

Moreover, the Fourier transform of $g=f'$ is

$$
\hat g_k = \ii k \hat f_k
$$
:::

:::{prf:proof}
(1) [@Tveito2005], Theorem 9.4 and Corollary 9.1

(2) [@Tveito2005], Theorem 9.3 and Theorem 8.1, or, [@Davis1963], Theorem 12.1.5.
:::

+++

:::{prf:example} Continuation of [](#ex:sqhat)

Here are the truncated Fourier series for the square hat function.

```{code-cell}
def fn(x,n):
    f = 0.5 * ones_like(x)
    for k in range(1,n,2):
        e1 = (-1)**(( k+1)//2) * exp( 1j*k*x) / ( k*pi)
        e2 = (-1)**((-k+1)//2) * exp(-1j*k*x) / (-k*pi)
        f += real(e1 + e2)
    return f

x = linspace(0,2*pi,1000)
f = (x > 0.5*pi)*(x < 1.5*pi)*(1) + 0.0
for i,n in enumerate([10,20,30,50]):
    subplot(2,2,i+1)
    plot(x, f, x, fn(x,n))
    text(pi,0.5,'n = ' + str(n),ha='center')
    grid(True);
```

As $n$ increases, we see that the errors do not decrease around the discontinuities and there is no convergence in maximum norm, but there is convergence in 2-norm as stated in Part (1) of [](#thm:fserconv).

We have observed this kind of Gibbs oscillations when we interpolate a  discontinuous function with polynomials.

Moreover, it looks like there is pointwise convergence away from the discontinuities, and this is the case, see next Theorem.
:::

+++

:::{prf:theorem}
Let $f$ be piecewise continuous and its periodic extension be one-sided differentiable for all $x \in \re$. Then

$$
(S_n f)(x) \to \half [f (x-) + f(x+)], \qquad x \in \re
$$

Hence, if $f$ is continuous at $x$, then $(S_n f)(x) \to f(x)$.
:::

:::{prf:proof}
See [@Tveito2005], Theorem 9.2.
:::

+++

:::{exercise}
Compute the Fourier transform of the triangular hat function in [](#ex:trihat). Verify the Fourier transform of the functions in the other examples also.

Plot the truncated Fourier series $(S_n f)(x)$ of all the functions in the above examples for a few values of $n$ as we did for the square hat function.

To test convergence as $n \to \infty$, compute $S_n f$ on a uniform grid $\{ x_j \}$ of say $m = 1000$ points and measure the error as

\begin{align}
\norm{S_n f - f}_\infty &\approx \max_{j}  |S_n f(x_j) - f(x_j)| \\
\norm{S_n f - f}_2 &\approx \left( \frac{2\pi}{m} \sum_{j}  [S_n f(x_j) - f(x_j)]^2 \right)^\half
\end{align}

Plot error versus $n$ in loglog scale or semilogy scale.
:::
