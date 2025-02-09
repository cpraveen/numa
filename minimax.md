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

If such a $p^*$ exists, then define

$$
E_n = E_n(f) := \min_{p \in \poly_n} \norm{f-p} = \norm{f-p^*}
$$ 

We will use the maximum norm 

$$
\norm{\phi} = \max_{x \in [-1,+1]} |\phi(x)|
$$ 

This type of approximation is also called a *minimax* approximation and $E_n$ is the minimum achievable error.

:::{prf:example} Degree 1 minimax for $\exp(x)$
Consider the function 

$$
f(x) = \exp(x), \qquad x \in [-1,1]
$$ 

and let us try to find the minimax approximation of the form

$$
p_1(x) = a_0 + a_1 x
$$ 

First try to interpolate at $x=-1,1$; 

```{code-cell}
:tags: hide-input
xmin, xmax = -1.0, 1.0
f = lambda x: exp(x)

xg = linspace(xmin,xmax,1000)
fg = f(xg)

xu = linspace(xmin,xmax,2)
fu = f(xu)
fug = barycentric_interpolate(xu,fu,xg)

print("Max error = ", abs(fg - fug).max())

figure(figsize=(10,4))
subplot(1,2,1)
plot(xu,fu,'o')
plot(xg, fg, label='exp(x)')
plot(xg, fug, label='$p_1(x)$')
legend(), xlabel('x'), grid(True)

subplot(1,2,2)
plot(xu,0*xu,'o')
plot(xg,fg-fug)
title('$f(x)-p_1(x)$')
xlabel('x'), grid(True);
```

we observe that error is small at end points and large in the middle. The error $f-p_1$ does not oscillate. Clearly, by sliding the straight line around, we see that we can reduce the error in the middle while increasing it at the end points.  Below we interpolate at $x=-0.5,0.7$.

```{code-cell}
:tags: hide-input
xu = array([-0.5,0.7])
fu = f(xu)
fug = barycentric_interpolate(xu,fu,xg)

print("Max error = ", abs(fg - fug).max())

figure(figsize=(10,4))
subplot(1,2,1)
plot(xu,fu,'o')
plot(xg, fg, label='exp(x)')
plot(xg, fug, label='$p_1(x)$')
legend(), xlabel('x'), grid(True)

subplot(1,2,2)
plot(xu,0*xu,'o')
plot(xg,fg-fug)
title('$f(x)-p_1(x)$')
xlabel('x'), grid(True);
```

The error is now smaller and the error curve will then oscillate taking negative and positive values. 

> The best choice is to make the maxima and minima of the error to be of same magnitude by adjusting the slope and position of the straight line. 

By looking at the error curve, we see that the maximum error is achieved at the end points

\begin{align}
f(-1) - p_1(-1) &= \ee^{-1} - (a_0 - a_1) = E_1 \\
f(1)-p_1(1)     &= \ee^1 - (a_0 + a_1) = E_1
\end{align}

There is an intermediate point $x_3$ where the error is $-E_1$

\begin{align}
f(x_3) - p_1(x_3) &= \ee^{x_3} - (a_0 + a_1 x_3) = -E_1 \\
f'(x_3) - p_1'(x_3) &= \ee^{x_3} - a_1 = 0
\end{align}

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

```{code-cell}
:tags: hide-input
p1 = lambda x: 1.2643 + 1.1752 * x

print("Max error = ", abs(fg - p1(xg)).max())

figure(figsize=(10,4))
subplot(1,2,1)
plot(xg, fg, label='exp(x)')
plot(xg, p1(xg), label='$p_1(x)$')
legend(), xlabel('x'), grid(True)

subplot(1,2,2)
plot(xg,fg-p1(xg))

a1 = 0.5*(exp(1)-1/exp(1))
x3 = log(a1)
E1 = 0.5*exp(-1) + 0.25*x3*(exp(1) - exp(-1))
plot([-1,-1],[0,E1]), text(-1,0.5*E1,'$+E_1$',ha='left')
plot([+1,+1],[0,E1]),  text(+1,0.5*E1,'$+E_1$',ha='left')
plot([x3,x3],[0,-E1]),  text(x3,-0.5*E1,'$-E_1$',ha='left')
plot([-1,x3,1],[0,0,0],'*')

title('$f(x)-p_1(x)$')
xlabel('x'), grid(True);
```

The error of minimax approximation takes the values $+E_1, -E_1, +E_1$ at some three points. We say that the error equioscillates at $1 + 2 = 3$ points.
:::

+++

:::{prf:definition}
We will say that a function $f : [a,b] \to \re$ equioscillates at $m$ points if there are $m$ points $x_i$ such that

$$
f(x_i) = (-1)^i E, \qquad 0 \le i \le m
$$

where

$$
E = \pm \max_{x \in [a,b]} |f(x)|
$$
:::

+++

As was shown by Chebyshev, see theorem below, equioscillation is a necessary and sufficient condition for the best approximation. 

+++

:::{prf:example} Cubic minimax for $\exp(x)$
We can find a cubic polynomial that gives the best approximation in maximum norm. The solution can be computed by the Remes algorithm and yields

$$
p_3^*(x) = 0.994579 + 0.995668 x + 0.542973 x^2 + 0.179533 x^3
$$ 

Compare the cubic minimax approximation to the cubic interpolating
polynomial.

```{code-cell}
:tags: hide-input
xmin, xmax = -1.0, 1.0
f = lambda x: exp(x)

# Minimax
p3 = lambda x: 0.994579 + 0.995668 * x + 0.542973 * x**2 + 0.179533 * x**3

# For plotting and error
xg = linspace(xmin,xmax,1000)
fg = f(xg)

n  = 3 # degree

# Uniform points
xu = linspace(xmin,xmax,n+1)
fu = f(xu)
fug = barycentric_interpolate(xu,fu,xg)

# Chebyshev points of I kind
m = n+1 # Number of points
j = linspace(0,m-1,m)
theta = (2*j+1)*pi/(2*m)
xc1 = cos(theta)
fc1 = f(xc1)
fc1g = barycentric_interpolate(xc1,fc1,xg)

# Chebyshev points of II kind
xc2 = cos(linspace(0,pi,n+1))
fc2 = f(xc2)
fc2g = barycentric_interpolate(xc2,fc2,xg)

print('Error, Minimax      = ', abs(fg-p3(xg)).max())
print('Error, Uniform      = ', abs(fg-fug).max())
print('Error, Chebyshev I  = ', abs(fg-fc1g).max())
print('Error, Chebyshev II = ', abs(fg-fc2g).max())

plot(xg, fg-p3(xg), label='Cubic minimax')
plot(xg, fg-fug, '--', label='Cubic interp: Uniform')
plot(xg, fg-fc1g, '--', label='Cubic interp: Chebyshev I')
plot(xg, fg-fc2g, '--', label='Cubic interp: Chebyshev II')
legend(), grid(True), xlabel('x');
```

We see that the error of minimax approximation takes maximum and minimum values of the same magnitude $E_3 \approx 0.0055$ alternately at five points;  the error equioscillates at $3+2=5$ points. The interpolation at Chebyshev points of first kind comes closest to minimax approximation.  
:::

+++

## Remes algorithm [@Veidinger1960]

+++

:::{prf:example} Using chebfun
We can use `chebfun` to compute and visualize the best approximation.  For the function $f(x) = |x|$ in $[-1,+1]$, the code `minimax_abs.m` shown below generates this figure.

```{figure} matlab/minimax_abs.svg
:width: 100%
```

For the minimax approximation of degree 4, the error oscillates at $4 + 3 = 7$ points.

```{literalinclude} matlab/minimax_abs.m
```

:::

+++

## Characterization of minimax

:::{prf:theorem} Chebyshev
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
