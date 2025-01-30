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

# Convergence of Chebyshev interpolation

```{include} math.md
```

```{code-cell}
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

[@Trefethen2019] provides more details and proof regarding these results.

## Chebyshev series

:::{prf:theorem}
If $f : [-1,1] \to \re$ is Lipschitz continuous, it has a unique representation as a Chebyshev series

$$
f(x) = \sum_{k=0}^\infty a_k T_k(x)
$$

which is absolutely and uniformly convergent. The coefficients are given by

$$
a_0 = \frac{1}{\pi} \int_{-1}^1 \frac{f(x)}{\sqrt{1-x^2}} \ud x, \qquad
a_k = \frac{2}{\pi} \int_{-1}^1 \frac{f(x) T_k(x)}{\sqrt{1-x^2}} \ud x, \qquad k \ge 1
$$
:::

:::{prf:theorem} Decay of Chebyshev coefficients
:label: thm:chebcoefconv
For any integer $\nu \ge 0$, let $f$ and its derivatives through $f^{(\nu-1)}$ be absolutely continuous on $[-1,1]$ and suppose the $\nu$'th derivative $f^{(\nu)}$ is of bounded variation $V$. Then for $k \ge \nu + 1$, the Chebyshev coefficients of $f$ satisfy

$$
|a_k| \le \frac{2 V}{\pi k (k-1) \ldots (k-\nu)} \le \frac{2V}{\pi(k - \nu)^{\nu+1}}
$$
:::

## Projection and interpolation

The Chebyshev projection is obtained by truncating the Chebyshev series

$$
f_n(x) = \sum_{k=0}^\infty a_k T_k(x)
$$

Let $p_n(x)$ be the interpolation of $f(x)$ at $n+1$ Chebyshev points, and let us write it as

$$
p_n(x) = \sum_{k=0}^n c_k T_k(x)
$$

### Differentiable functions

:::{prf:theorem} Convergence for differentiable functions
If $f$ satisfies the conditions of [](#thm:chebcoefconv) with $V$ again denoting the total variation of $f^{(\nu)}$ for some $\nu \ge 1$, then for any $n > \nu$, its Chebyshev projections satisfy

$$
||f - f_n|| \le \frac{2 V}{\pi \nu (n - \nu)^\nu}
$$

and its Chebyshev interpolants satisfy

$$
||f - p_n|| \le \frac{4 V}{\pi \nu (n - \nu)^\nu}
$$
:::

Thus $f^{(\nu)}$ of bounded variation implies that the convergence is at an algebraic rate of $\order{n^{-\nu}}$ for $n \to \infty$. The more derivative the function has, the faster is the convergence.

For $\nu = 0$, e.g., $f(x) = \textrm{sign}(x)$, we cannot hope for convergence since polynomials are continuous and sign function is not.

For $\nu = 1$, e.g., $f(x) = |x|$, we have $V=2$ and the above theorem predicts convergence since

$$
\norm{f - f_n} \le \frac{4}{\pi(n-1)}, \qquad \norm{f - p_n} \le \frac{8}{\pi (n-1)}
$$

+++

:::{prf:example}
Numerically compute the error in the interpolation of $f(x) = |x|$.

The degree 10 interpolant looks like this.

```{code-cell}
f = lambda x: abs(x)
n = 10
xp = linspace(-1,1,1000)
theta = linspace(0,pi,n+1)
x = -cos(theta)
y = f(x)
yp = barycentric_interpolate(x,y,xp)
plot(x,y,'o',label='Data')
plot(xp,f(xp),label='Exact')
plot(xp,yp,label='Interpolant')
legend(), xlabel('x'), ylabel('y')
```

Now compute error norm for increasing degree. The norm is only an approximation.

```{code-cell}
n = 10
nvalues,evalues = [],[]
for i in range(10):
    theta = linspace(0,pi,n+1)
    x = -cos(theta)
    y = f(x)
    yp = barycentric_interpolate(x,y,xp)
    e = abs(f(xp) - yp).max()
    nvalues.append(n), evalues.append(e)
    n = 2 * n
nvalues = array(nvalues)
loglog(nvalues,evalues,'o-',label="Error")
loglog(nvalues, 8/(pi*(nvalues-1)),label="$8/(\\pi (n-1))$")
legend(), xlabel('n'), ylabel('||f-p||')
:::

+++

### Analytic functions

For any $\rho > 1$, the Bernstein ellipse $E_\rho$ is the ellipse with foci at $x=-1$ and $x=+1$ and 

$$
\rho = \textrm{semi-minor axis + semi-major axis}
$$

It can be obtained by mapping the circle of radius $\rho$, $\xi = \rho \exp(\ii \theta)$, under the Joukowsky map

$$
z = \half (\xi + \xi^{-1})
$$

```{code-cell}
theta = linspace(0,2*pi,500)
for rho in arange(1.1,2.1,0.1):
    xi = rho * exp(1j * theta)
    z = 0.5 * (xi + 1/xi)
    plot(real(z), imag(z),'k-')
plot([-1,1],[0,0],'r-'), axis('equal')
```

:::{prf:theorem} Decay of Chebyshev coefficients: analytic functions
:label: thm:chebcoefconva
Let $f : [-1,1] \to \re$ be analytically continuable to the open Bernstein ellipse $E_\rho$ for some $\rho > 1$, where it satisfies $|f(x)| \le M$ for some $M$. Then its Chebyshev coefficients satisfy

$$
|a_0| \le M, \qquad |a_k| \le \frac{2M}{\rho^k}, \qquad k \ge 1
$$
:::

:::{prf:theorem} Convergence for analytic functions
If $f$ has the properties of [](#thm:chebcoefconva), then for each $n \ge 0$, its Chebyshev projections satisfy

$$
\norm{f - f_n} \le \frac{2 M \rho^{-n}}{\rho - 1}
$$

and its Chebyshev interpolants satisfy

$$
\norm{f - p_n} \le \frac{4 M \rho^{-n}}{\rho - 1}
$$
:::
