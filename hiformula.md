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
math:
  '\ii': '\mathfrak{i}'
  '\df': '\frac{\partial #1}{\partial #2}'
  '\ud': '\textrm{d}'
  '\norm': '\|#1\|'
  '\poly': '\mathbb{P}'
---

# Hermite integral formula

## Complex Analysis

We first recall some results from complex variables. 

:::{prf:definition}
A function $f(z) = u(x,y) + \ii v(x,y)$ is differentiable at a point $z=x+\ii y$ if

$$
\lim_{\Delta z \to 0} \frac{f(z + \Delta z) - f(z)}{\Delta z}
$$ 

exists and is unique.
:::

The necessary and sufficient condition is given by

:::{prf:theorem}
A function $f(z) = u(x,y) + \ii v(x,y)$ is differentiable at a point $z = x + \ii y$ of a region in the complex plane if and only if the partial derivatives $u_x, u_y, v_x, v_y$ are continuous and satisfy the
Cauchy-Riemann equations

$$
\df{u}{y} = -\df{v}{x}, \qquad \df{u}{x} = \df{v}{y}
$$
:::

:::{prf:definition}
1. A function $f(z)$ is analytic at $z=z_0$ if it is differentiable in a neighbourhood of $z_0$. 

2. A function $f(z)$ is analytic in a region of the complex plane if it is analytic at every point of this
region. 

3. A function $f(z)$ is an entire function if it is analytic in the finite part of the complex plane.
:::

:::{prf:theorem}
Let $f(z)$ be anaytic interior to and on a simple closed contour $C$ in
the complex plane. Then at any interior point $z$

$$
f(z) = \frac{1}{2\pi\ii} \oint_C \frac{f(\xi)}{\xi - z} \ud \xi
$$
:::

:::{prf:theorem}
Let $f(z)$ be analytic for $|z-z_0| \le R$. Then

$$
f(z) = \sum_{j=0}^\infty b_j (z-z_0)^j, \qquad b_j = \frac{1}{j!} f^{(j)}(z_0)
$$

converges uniformly in $|z-z_0| \le R_1 < R$.
:::

## Laurent expansion and residue

If $f(z)$ is not analytic at a single point $z_0$, then we can use a Laurent series. Let $f(z)$ be analytic in the region

$$
D = \{ z : 0 < |z-z_0| < \rho \}
$$ 

and let $z_0$ be an isolated singular point of $f(z)$. The Laurent expansion of $f(z)$ in $D$ is given by 

$$
f(z) = \sum_{n=-\infty}^\infty C_n (z-z_0)^n
$$ 

with

$$
C_n = \frac{1}{2\pi \ii} \oint_C \frac{f(z) \ud z}{(z-z_0)^{n+1}}
$$

where $C$ is a simple closed contour lying in $D$. The coefficient $C_{-1}$ is called the *residue* of $f(z)$ at $z_0$. When $n=-1$, we have 

$$
\oint_C f(z) \ud z = 2\pi \ii C_{-1}
$$

## Lagrange interpolation

Let $x_0,x_1,\ldots,x_n$ be a set of $n+1$ distinct interpolation points. Define 

$$
\ell(x) = (x-x_0)\ldots(x-x_n)
$$ 

The Lagrange polynomials are given by

$$
\ell_j(x) = \frac{\ell(x)}{\ell'(x_j)(x-x_j)}
$$ 

and the interpolant of degree $n$ is 

$$
p(x) = \sum_{j=0}^n f(x_j) \ell_j(x)
$$ 

Let $x$ be any point which is distinct from the interpolation points. Let $\Gamma_j$ be a contour in the complex plane that encloses $x_j$ but none of the other points, nor the point $x$. Using the Cauchy integral formula, we see that 

$$
\begin{aligned}
& \frac{1}{2\pi\ii} \oint_{\Gamma_j} \frac{\ud t}{\ell(t) (x-t)} \\
=& \frac{1}{2\pi\ii} \oint_{\Gamma_j} \frac{\phi(t)}{t-x_j} \ud t \\
& \textrm{where}  \qquad \phi(t) =
\frac{1}{(t-x_0)\ldots(t-x_{j-1})(t-x_{j+1})\ldots(t-x_n)(x-t)} \\
=& \phi(x_j) \\
=& \frac{1}{(x_j-x_0)\ldots(x_j-x_{j-1})(x_j-x_{j+1})\ldots(x_j-x_n)(x-x_j)} \\
=& \frac{1}{\ell'(x_j)(x-x_j)}
\end{aligned}
$$ 

Hence

$$
\ell_j(x) = \frac{\ell(x)}{\ell'(x_j)(x-x_j)} = \frac{1}{2\pi\ii} \oint_{\Gamma_j}
\frac{\ell(x) \ud t}{\ell(t) (x-t)}
$$ 

Now let $\Gamma'$ be a curve which encloses all of the points $\{x_j\}$ but not the point $x$ and let $f$ be analytic on and interior to $\Gamma'$. Then we can combine the $\Gamma_j$ integrals to get an expression for the interpolant

$$
\begin{aligned}
p(x) =& \frac{1}{2\pi\ii} \sum_j \oint_{\Gamma_j} f(x_j) \frac{\ell(x) \ud t}{\ell(t) (x-
t)} \\
=& \frac{1}{2\pi\ii} \sum_j \oint_{\Gamma_j} f(t) \frac{\ell(x) \ud t}{\ell(t) (x-t)} \\
=& \frac{1}{2\pi\ii} \oint_{\Gamma'} \frac{\ell(x) f(t) \ud t}{\ell(t) (x-t)}
\end{aligned}
$$ 

by using the same type of argument we used above.

Now suppose we enlarge the contour of integration to a new contour $\Gamma$ that encloses $x$ as well as $\{x_j\}$, and we assume $f$ is analytic on and inside $\Gamma$. Then 

$$
\begin{aligned}
\frac{1}{2\pi\ii} \oint_{\Gamma} \frac{\ell(x) f(t) \ud t}{\ell(t) (x-t)} =&
\frac{1}{2\pi\ii} \oint_{\Gamma'} \frac{\ell(x) f(t) \ud t}{\ell(t) (x-t)} + \frac{1}
{2\pi\ii} \oint_{\Gamma''} \frac{\ell(x) f(t) \ud t}{\ell(t) (x-t)} \\
=& p(x) - f(x)
\end{aligned}
$$ 

where we used the Cauchy integral formula to get $-f(x)$.

:::{prf:theorem}
Let $f$ be analytic in a region $\Omega$ containing distinct points $x_0,x_1,\ldots,x_n$, and let $\Gamma$ be a contour in $\Omega$ enclosing these points in the positive direction. The polynomial interpolant $p \in \poly_n$ to $f$ at these nodes is

$$
p(x) = \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{f(t)[\ell(t) - \ell(x)]}{\ell(t) (t-x)} \ud t
$$

and if $x$ is enclosed by $\Gamma$, the error in the interpolant is

$$
f(x) - p(x) = \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{\ell(x)}{\ell(t)} \frac{f(t)}{(t-x)}
\ud t
$$
:::

:::{prf:proof}
(1) If $\Gamma$ encloses $x$, then by Cauchy integral formula, $f(x)$ can be written as

$$
f(x) = \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{f(t)}{(t-x)} \ud t = \frac{1}{2\pi\ii}
\oint_{\Gamma} \frac{\ell(t) f(t)}{\ell(t) (t-x)} \ud t
$$ 

From previously derived formula for $p(x)-f(x)$, 

$$
\begin{aligned}
p(x) =& f(x) + \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{\ell(x) f(t)}{\ell(t) (x-t)} \ud t\\
=& \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{\ell(t) f(t)}{\ell(t) (t-x)} \ud t - \frac{1}
{2\pi\ii} \oint_{\Gamma} \frac{\ell(x) f(t)}{\ell(t) (t-x)} \ud t\\
=& \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{f(t)[\ell(t) - \ell(x)]}{\ell(t) (t-x)} \ud t
\end{aligned}
$$ 

In this formula, the integrand does not have a pole at $t=x$ and hence it is valid even if $\Gamma$ does not enclose $x$.

(2) The formula for $f(x)-p(x)$ has already been derived.
:::

## Effect of point distribution $\{x_j\}$

On a fixed contour $\Gamma$, the quantities $f(t)$ and $t-x$ in the error formula are independent of $\{x_j\}$. The error depends on the ratio

$$
\frac{\ell(x)}{\ell(t)} = \frac{ \prod\limits_{j=0}^n (x-x_j)}{\prod\limits_{j=0}^n(t-
x_j)}
$$ 

If $\Gamma$ is far away from $\{x_j\}$, then for each $t \in \Gamma$,

$$
\left| \frac{\ell(x)}{\ell(t)} \right| \approx \frac{|\ell(x)|}{|t|^{n+1}}, \qquad |t| \to
\infty
$$ 

this ratio will shrink exponentially as $n \to \infty$, and if this happens, we may conclude that $p(x)$ converges exponentially to $f(x)$ as $n \to \infty$. The crucial condition is that *it must be possible to continue $f$ analytically far out to $\Gamma$*.

:::{prf:example}
Let $S$ be the *stadium* in the complex plane consisting of all points lying at a distance $\le 2$ from $[-1,+1]$. Suppose $f$ is analytic in a larger region $\Omega$ that includes a curve $\Gamma$ enclosing $S$.  Since

$$
|x-x_j| \le 2, \qquad |t-x_j| > 2, \qquad \frac{|t-x_j|}{|x-x_j|} > 1
$$

there exists a $\gamma$ such that

$$
\frac{|t-x_j|}{|x-x_j|} \ge \gamma > 1
$$ 

This implies

$$
\left| \frac{\ell(x)}{\ell(t)} \right| \le \gamma^{-(n+1)}, \qquad \forall x \in [-1,+1],
\quad t \in \Gamma
$$ 

and thus by the error formula

$$
\norm{f-p} = O(\gamma^{-n})
$$ 

Note that this conclusion applies regardless of the distribution of interpolation points in $[-1,+1]$.  They could be equally spaced or random[^1].
:::

:::{note}
So convergence of polynomial interpolants to analytic functions on $[-1,+1]$ is all about how small $\ell(x)$ is on $[-1,+1]$ compared to how large it is on a contour $\Gamma$ inside which $f$ is analytic.
:::

[^1]: In practice, rounding errors may spoil the convergence.
