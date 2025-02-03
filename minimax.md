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

# Minimax approximation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

## Best approximation

Given the degree $n \ge 0$, find a $p^* \in \poly_n$ for which the error
is minimized, i.e.,

$$
\norm{f-p^*} \le \norm{f-p}, \qquad \forall p \in \poly_n
$$ 

If such a
$p^*$ exists, then

$$
E_n := \inf_{p \in \poly_n} \norm{f-p} = \norm{f-p^*}
$$ 

We will use
the maximum norm 

$$
\norm{\phi} = \max_{-1 \le x \le +1} |\phi(x)|
$$ 

This type of approximation is also called a *minimax* approximation and $E_n$ is the minimum achievable error.

## Degree 1 minimax for exp

Consider the function 

$$
f(x) = \exp(x), \qquad x \in [-1,1]
$$ 

and let us try to find the minimax approximation of the form

$$
p_1(x) = a_0 + a_1 x
$$ 

First try to interpolate at $x=-1,1$; we observe that error is small at end points and large in the middle. The error $f-p_1$ does not oscillate. Clearly, by sliding the straight line around, we see that we can reduce the error in the middle while increasing it at the end points. The error curve will then oscillate taking negative and positive values. The best choice is to make the maxima and minima of the error to be of same magnitude by adjusting the slope and position of the straight line. We observe that there are two points $x_1$ and $x_2$ where $p_1$ interpolates $f$

$$
p_1(x_1) = f(x_1), \qquad p_2(x_2) = f(x_2)
$$ 

but this is not useful to solve the problem since we do not know $x_1$ and $x_2$.

Instead by looking at the error curve, we see that the maximum error is achieved at the end points

$$
f(-1) - p_1(-1) = \ee^{-1} - (a_0 - a_1) = E_1, \qquad f(1)-p_1(1) = \ee^1 - (a_0 + a_1) = E_1
$$ 

There is an intermediate point $x_3$ where

$$
f(x_3) - p_1(x_3) = \ee^{x_3} - (a_0 + a_1 x_3) = -E_1, \qquad f'(x_3) - p_1'(x_3) = \ee^{x_3} - a_1 = 0
$$ 

We have four equations for the four unknowns $a_0, a_1, E_1, x_3$ and the solution is 

$$
\begin{aligned}
a_1 &= \half(\ee - \ee^{-1}) \approx 1.1752 \\
x_3 &= \log(a_1) \approx 0.1614 \\
E_1 &= \half \ee^{-1} + \frac{x_3}{4}(\ee - \ee^{-1}) \approx 0.2788 \\
a_0 &= E_1 + (1-x_3) a_1 \approx 1.2643
\end{aligned}
$$ 

so that 

$$
p_1^*(x) \approx 1.2643 + 1.1752 x
$$

## Cubic minimax for exp

We can find a cubic polynomial that gives the best approximation in maximum norm. The solution can be computed by the Remez algorithm and yields

$$
p_3^*(x) = 0.994579 + 0.995668 x + 0.542973 x^2 + 0.179533 x^3
$$ 

We see that the error takes maximum and minimum values of the same magnitude $E_3 \approx 0.0055$ alternately at five points. We say that the error
equioscillates at 3+2 = 5 points.

<figure>
<div class="center">
<p><img src="figs/minimax_exp1" style="width:48.0%" alt="image" /> <img
src="figs/minimax_exp2" style="width:48.0%" alt="image" /></p>
</div>
<figcaption>Cubic minimax of <span
class="math inline"><em>f</em>(<em>x</em>) = exp (<em>x</em>)</span> on
<span class="math inline">[−1, 1]</span></figcaption>
</figure>

:::{prf:remark}
Compare the cubic minimax approximation to the cubic interpolating
polynomial.
:::

:::{prf:example} Using chebfun
We can use `chebfun` to compute and visualize the best approximation.
For the function $f(x) = |x|$ in $[-1,+1]$, the code `minimax_abs.m`
shows the best approximation and error curve.
:::

:::{prf:theorem}
(1) A continuous function $f:[-1,+1] \to \re$ has a unique best approximation $p^* \in \poly_n$. 

(2) Any polynomial $p \in \poly_n$ is equal to $p^*$ iff $f-p$ equioscillates in atleast $n+2$ extreme points.
:::

:::{prf:proof}
(1) **Existence of best approximation**: The function

$$
\phi: \poly_n \to \re,  \qquad \phi(p) = \norm{f-p}
$$ 

is a continuous function[^1]. Since one candidate approximation is the zero polynomial, we know that if $p^*$ exists, then $\norm{f- p^*} \le \norm{f-0} = \norm{f}$, and hence $p^*$ lies in

$$
V_n = \{ p \in \poly_n : \norm{f-p} \le \norm{f} \} \subset \poly_n
$$

$V_n$ is a closed[^2] and bounded[^3] subset of a finite dimensional space $\poly_n$. Hence $V_n$ is compact; any continuous function on a compact set achieves its minimum and maximum. This proves the existence of $p^*$.

(2) **Equioscillation $\implies$ optimality**: Suppose $f-p$ takes equal extreme values with alternating signs at $n+2$ points $x_0 < x_1 < \ldots x_{n+1}$ and suppose there exists a $q \in \poly_n$ such that $\norm{f-q} < \norm{f-p}$. Now by our assumption of equioscillation 

$$
(f-p)(x_i) = (-1)^i E, \qquad i=0,1,\ldots,n+1
$$ 

so that 

$$
(p-q)(x_i) = (p-f)(x_i) + (f-q)(x_i) = -(-1)^i E + (f-q)(x_i)
$$

Hence 

$$
(p-q)(x_0) = -E + \underbrace{(f-q)(x_0)}_{<E} < 0
$$

$$
(p-q)(x_1) = E + \underbrace{(f-q)(x_1)}_{> -E} > 0
$$ 

etc. Hence $p-q$ takes non-zero values with alternating signs at the $n+2$ equioscillation points, implying that $p-q=0$ in atleast $n+1$ points in between. Since $p-q$ is of degree $n$, this implies that $p=q$ everywhere.

(3) **Optimality $\implies$ equioscillation at $n+2$ points**: Given that $E_n = \norm{f-p}$, suppose $f-p$ equioscillates at fewer than $n+2$ points. Without loss of generality, suppose the leftmost extremum is one where $f-p$ takes the value $- E_n$. Then there are numbers $-1 < \xi_1 < \ldots < \xi_k < +1$ with $k \le n$ where $f-p$ is zero. Then $p(x) - \epsilon (\xi_1-x)\ldots(\xi_k-x)$ is a better approximation than $p$ for small $\epsilon$. We can check this as follows. If $x \in [-1,\xi_1]$, then

$$
(f-p)(x) + \epsilon (\xi_1-x)\ldots(\xi_k-x) \ge -E_n + \epsilon (+) > -E_n
$$

If $x \in [\xi_1,\xi_2]$, then

$$
(f-p)(x) + \epsilon (\xi_1-x)\ldots(\xi_k-x) \le E_n + \epsilon (-) < E_n
$$

etc. Hence $p$ is not the best approximation. However, if $f-p$ equioscillates at atleast $n+2$, then we cannot improve it further by this argument.

(3) **Uniqueness**: Suppose $p$ and $q$ are distinct best approximations. Define $r = \half (p+q)$. Since 

$$
E_n = \norm{f-p} = \norm{f-q}
$$ 

we have

$$
\norm{f-r} \le \half \norm{f-p} + \half \norm{f-q} = E_n
$$ 

We cannot have $\norm{f-r} < E_n$, hence $\norm{f-r} = E_n$ and $r$ is also a best approximation. So $r$ must equioscillate at $n+2$ extreme points. At any of these points $x_k$, suppose

$$
(f-r)(x_k) &= E_n \\
&\Downarrow \\
\frac{f(x_k)-p(x_k)}{2} + \frac{f(x_k)-q(x_k)}{2} &= E_n
$$ 

But

$$
-E_n \le (f-p)(x_k) \le E_n, \qquad -E_n \le (f-q)(x_k) \le E_n
$$ 

We cannot have 

$$
(f-p)(x_k) < E_n \qquad\textrm{and/or}\qquad (f-q)(x_k) < E_n
$$ 

and hence

$$
E_n = f(x_k) - p(x_k) = f(x_k) - q(x_k) \qquad \Longrightarrow \qquad p(x_k) = q(x_k)
$$

A similar conclusion follows if $(f-r)(x_k) = -E_n$.  Since $p=q$ at $n+2$ distinct points but each is of degree $\le n$, we conclude that $p=q$ everywhere.
:::

:::{prf:remark}
The error curve for a best approximation may have more than $n+2$ points
of equioscillation, and this will happen if $f$ and $n$ are both even or
both odd. For the function $f(x)=|x|$, the degree $2$ polynomial
equioscillates at 5 points and the degree 4 polynomial equioscillates at
7 points.
:::

[^1]: This is true for any normed vector space and is a consequence of
    $|\norm{x} - \norm{y}| \le \norm{x-y}$.

[^2]: Show that it is closed by showing all limit points belong to
    $V_n$.

[^3]: It is bounded because
    $$\norm{p} = \norm{p-f+f} \le \norm{p-f} + \norm{f} \le \norm{f} + \norm{f} =
    2\norm{f} < \infty$$
