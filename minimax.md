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

Given integer $n \ge 0$ and a function $f$, we can approximate it by a degree $n$ interpolating polynomial $p$. The approximation depends on the choice of nodes; different choice lead to different error $f-p$. Is there a polynomial which gives the smallest error ?

## Best approximation in max norm

Let $f : [-1,1] \to \re$ be a continuous function. Given the degree $n \ge 0$, find a $p^* \in \poly_n$ for which the error
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

This type of approximation is also called a *minimax* approximation since

$$
E_n = \min_{p \in \poly_n} \max_{x \in [-1,1]} |f(x) - p(x)| = \norm{f-p^*}
$$

and $E_n$ is the minimum achievable approximation error by any polynomial of degree $n$.

+++

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

> The best choice is to make the maxima and minima of the error to be of same magnitude but opposite sign by adjusting the slope and position of the straight line. 

**Minimax approximation.** By looking at the error curve, we see that the maximum error is achieved at the end points

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
plot([-1,-1],[0,E1]), text(-1,0.5*E1,'$+E_1$',ha='center',backgroundcolor='white')
plot([+1,+1],[0,E1]),  text(+1,0.5*E1,'$+E_1$',ha='center',backgroundcolor='white')
plot([x3,x3],[0,-E1]),  text(x3,-0.5*E1,'$-E_1$',ha='center',backgroundcolor='white')
plot([-1,x3,1],[0,0,0],'*')

title('$f(x)-p_1(x)$')
xlabel('x'), grid(True);
```

The error of minimax approximation takes the values $+E_1, -E_1, +E_1$ at some three points. We say that the error equioscillates at $1 + 2 = 3$ points.
:::

+++

:::{prf:definition} Equioscillation
We will say that a function $\phi : [a,b] \to \re$ equioscillates at $m$ points if there are $m$ points $x_i$ such that

$$
\phi(x_i) = \pm (-1)^i E, \qquad 0 \le i \le m-1
$$

where

$$
E = \max_{x \in [a,b]} |\phi(x)|
$$
:::

+++

As we show in theorem below, equioscillation of the error $f-p$ is a necessary and sufficient condition for the best approximation. 

+++

:::{prf:example} Cubic minimax for $\exp(x)$
We can find a cubic polynomial that gives the best approximation in maximum norm. The solution can be computed by the Remez algorithm and yields

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

:::{prf:example} Using chebfun
We can use [chebfun](https://www.chebfun.org) to compute and visualize the best approximation.  For the function $f(x) = |x|$ in $[-1,+1]$, the Matlab code `minimax_abs.m` shown below generates this figure.

```{figure} matlab/minimax_abs.svg
:width: 100%
:align: center

Degree 4 minimax approximation of $f(x) = |x|$, $x\in [-1,1]$ and its error which equioscillates at 7 points.
```

For the minimax approximation of degree 4, the error oscillates at $4 + 3 = 7$ points.

```{literalinclude} matlab/minimax_abs.m
```

:::

+++

## Characterization of minimax

:::{prf:theorem}
(1) A continuous function $f:[-1,+1] \to \re$ has a unique best approximation $p^* \in \poly_n$. 

(2) Any polynomial $p \in \poly_n$ is equal to $p^*$ iff $f-p$ equioscillates in atleast $n+2$ extreme points.
:::

:::{prf:proof}
(1) **Existence of best approximation**: The function

$$
\phi: \poly_n \to \re,  \qquad \phi(p) = \norm{f-p}
$$ 

is a continuous[^1]. Since one candidate approximation is the zero polynomial, we know that if $p^*$ exists, then $\norm{f- p^*} \le \norm{f-0} = \norm{f}$, and hence $p^*$ lies in

$$
V_n = \{ p \in \poly_n : \norm{f-p} \le \norm{f} \} \subset \poly_n
$$

$V_n$ is a closed[^2] and bounded[^3] subset of a finite dimensional space $\poly_n$. Hence $V_n$ is compact; any continuous function on a compact set achieves its minimum and maximum. This proves the existence of $p^*$.

(2) **Equioscillation $\implies$ optimality**: Suppose $f-p$ takes equal extreme values $\pm E$ with alternating signs at $n+2$ points $x_0 < x_1 < \ldots x_{n+1}$; wlog we assume error $   = +E$ at $x=x_0$

$$
(f-p)(x_i) = (-1)^i E, \qquad i=0,1,\ldots,n+1
$$ 

Suppose there exists a $q \in \poly_n$ with smaller error,  $\norm{f-q} < \norm{f-p}$. Then

$$
(p-q)(x_i) = (p-f)(x_i) + (f-q)(x_i) = -(-1)^i E + (f-q)(x_i)
$$

Hence, since $-E \le (f-q)(x) \le +E$, we get

\begin{align}
(p-q)(x_0) &= -E + \underbrace{(f-q)(x_0)}_{<E} < 0 \\
(p-q)(x_1) &= +E + \underbrace{(f-q)(x_1)}_{> -E} > 0 \\
           &\textrm{etc.}
\end{align}

Hence $p-q$ takes non-zero values with alternating signs at the $n+2$ equioscillation points, implying that $p-q=0$ in atleast $n+1$ points in between. Since $p-q$ is of degree $n$, this implies that $p=q$ everywhere, and hence $\norm{f-q} = \norm{f - p}$.

(3) **Optimality $\implies$ equioscillation at $n+2$ points**: Let $E = \norm{f-p}$ and suppose $f-p$ equioscillates at fewer than $n+2$ points, say at $k+1 < n+2$ points. Without loss of generality, suppose the leftmost extremum is one where $f-p$ takes the value $-E$. Then there are points $-1 < \xi_1 < \ldots < \xi_k < +1$ with $k \le n$ where $f-p$ is zero. Then $p(x) - \epsilon (\xi_1-x)\ldots(\xi_k-x)$ is a better approximation of degree $n$ than $p$ for small $\epsilon > 0$. 

We can check this as follows. If $x \in [-1,\xi_1]$, then

$$
(f-p)(x) + \epsilon (\xi_1-x)\ldots(\xi_k-x) \ge -E + \epsilon (+) > -E
$$

If $x \in [\xi_1,\xi_2]$, then

$$
(f-p)(x) + \epsilon (\xi_1-x)\ldots(\xi_k-x) \le E + \epsilon (-) < E
$$

etc. Hence $p$ is not the best approximation. However, if $f-p$ equioscillates at atleast $n+2$ points, then we cannot improve it further by this argument since $p(x) - \epsilon (\xi_1-x)\ldots(\xi_k-x)$ will be of degree $\ge n+1$.

(4) **Uniqueness**: Suppose $p$ and $q$ are distinct best approximations. Define $r = \half (p+q)$. Since 

$$
E_n = \norm{f-p} = \norm{f-q}
$$ 

we have

$$
\norm{f-r} \le \half \norm{f-p} + \half \norm{f-q} = E_n
$$ 

We cannot have $\norm{f-r} < E_n$, hence $\norm{f-r} = E_n$ and $r$ is also a best approximation. So $r$ must equioscillate in atleast $n+2$ extreme points. At any of these points $x_k$, suppose

\begin{align}
(f-r)(x_k) &= +E_n \\
&\Downarrow \\
\frac{f(x_k)-p(x_k)}{2} + \frac{f(x_k)-q(x_k)}{2} &= +E_n
\end{align}

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

+++

:::{prf:remark}
The error curve for a best approximation may have more than $n+2$ points
of equioscillation, and this will happen if $f$ and $n$ are both even or
both odd. For the function $f(x)=|x|$, the degree $2$ polynomial
equioscillates at 5 points and the degree 4 polynomial equioscillates at
7 points.
:::

+++

## Remez algorithm

The algorithm to find the best approximation was proposed by E. Remez (1934) and is explained in [@Veidinger1960]. Let $f : [-1,1] \to \re$ and let us write $p(x)$ in terms of Chebyshev polynomials

$$
p(x) = \sum_{j=0}^n a_j T_j(x)
$$

1. Choose $n+2$ points $\{x_0, x_1, \ldots, x_{n+1}\} \subset [-1,1]$ and such that
    \begin{align}
    f(x_i) - p(x_i) &= (-1)^i E \\
    p(x_i)  + (-1)^i E &= f(x_i) \\
    \sum_{j=0}^n a_j T_j(x_i) + (-1)^i E &= f(x_i), \qquad 0 \le i \le n+1
    \end{align}
    This is a linear system of $n+2$ equations for the $n+2$ unknowns $\{a_0, a_1, \ldots, a_{n+1}, E\}$ which can be solved. The $E$ may be positive or negative.

1. If the condition $E = \pm \max_x |f(x)|$ is satisfied, then stop the iteration. Otherwise continue with the next step.

1. The error oscillates in sign at the $n+2$ points so it is zero at some $n+1$ intermediate points $\{z_0, z_1, \ldots, z_n\}$.
    $$
    \textrm{Find} \qquad z_j \in [x_j, x_{j+1}], \qquad f(z_j) - p(z_j) = 0,
    \qquad 0 \le j \le n
    $$
    by using some root finding method.

1. In the $n+2$ intervals $[-1,z_0], [z_0, z_1], \ldots, [z_n, 1]$, find the points $\{x_0^*, x_1^*, \ldots, x_{n+1}^*\}$ where the error $f-p$ takes maximum or minimum value. Do this by solving $(f'-p')(x)=0$ and also check at the end points of the intervals.

1. Set $x_j = x_j^*$ and go to Step 1.

We can observe that the algorithm is complicated and expensive. The need to find several roots and maxima/minima adds to the cost.

+++

:::{seealso}
For more discussion on Remez algorithm, see

1. [Notes and Matlab code by Sharif Tawfik](https://in.mathworks.com/matlabcentral/fileexchange/8094-remez-algorithm)
1. [Remez Method in Boost library](https://live.boost.org/doc/libs/1_69_0/libs/math/doc/html/math_toolkit/remez.html)
1. [E. D. de Groot, Master Thesis](https://fse.studenttheses.ub.rug.nl/15985/1/BSc_Math_2017_deGroot_ED.pdf), includes Matlab code
:::

+++

## Minimax versus interpolation

Given degree $n \ge 0$ let $p_n^*$ be the minimax approximation and $p_n$ be some other interpolation which can be written as

$$
p_n(x) = \sum_{j=0}^n f(x_j) \ell_j(x)
$$

where $\ell_j$ are the Lagrange polynomials. Then we can write

$$
p_n^*(x) = \sum_{j=0}^n p_n^*(x_j) \ell_j(x)
$$

and

\begin{align}
p_n^*(x) - p_n(x) &= \sum_{j=0}^n (p_n^*(x_j) - f(x_j)) \ell_j(x) \\
|p_n^*(x) - p_n(x)| &\le \sum_{j=0}^n |p_n^*(x_j) - f(x_j)| \cdot |\ell_j(x)| \\
&\le \left(\sum_{j=0}^n |\ell_j(x)| \right) \norm{f - p_n^*} \\
\norm{p_n^* - p_n} &\le \left( \max_x \sum_{j=0}^n |\ell_j(x)| \right) \norm{f - p_n^*}
\end{align}

Define the **Lebesgue constant**

$$
\Lambda_n = \max_x \sum_{j=0}^n |\ell_j(x)|
$$

whose value depends on the distribution of points. Therefore

\begin{align}
\norm{f - p_n}
&= \norm{f - p_n^* + p_n^* - p_n} \\
& \le \norm{f - p_n^*} + \norm{p_n^* - p_n} \\
&\le (1 + \Lambda_n) \norm{f - p_n^*}
\end{align}

The interpolation error is atmost a factor $1 + \Lambda_n$ larger than the error of the best approximation. If $\Lambda_n \gg 1$, then $p_n$ may be too far from the best approximation error.

+++

:::{prf:example} Lebesgue function
Define the Lebesgue function

$$
\lambda_n(x) = \sum_{j=0}^n |\ell_j(x)|
$$

```{code-cell}
:tags: remove-input
def Lagrange(j,xi,x):
    f = ones_like(x)
    n = len(xi)
    for k in range(n):
        if k != j:
            f = f * (x - xi[k]) / (xi[j] - xi[k])
    # If x close to xi[j], set to 1
    f[argwhere(abs(xi[j]-x) < 1.0e-13)] = 1.0
    return f

def LebesgueFunction(xi, x):
    f = zeros_like(x)
    n = len(xi)
    for i in range(n):
        f += abs(Lagrange(i,xi,x))
    return f

def PlotLebesgueFunction(xi,logy=False,ng=100):
    if logy:
        semilogy(xi,ones_like(xi),'o')
    else:
        plot(xi,0*xi,'o')
    for i in range(len(xi)-1):
        xg = linspace(xi[i],xi[i+1],ng)
        if logy:
            semilogy(xg, LebesgueFunction(xi,xg), 'k-')
        else:
            plot(xg, LebesgueFunction(xi,xg), 'k-')
    grid(True), xlabel('x');
```

**With uniform points**

```{code-cell}
:tags: remove-input
n  = 14
xi = linspace(-1,1,n+1)
PlotLebesgueFunction(xi)
title('Lebesgue function for '+str(n+1)+' uniform points');
```

```{code-cell}
:tags: remove-input
n  = 29
xi = linspace(-1,1,n+1)
PlotLebesgueFunction(xi,logy=True)
title('Lebesgue function for '+str(n+1)+' uniform points');
```

Note the log scale in above figure. The large values arise near the end points where we also observe Runge phenomenon.

**With Chebyshev points**

```{code-cell}
:tags: remove-input
n  = 14
xi = cos(linspace(0,pi,n+1))
PlotLebesgueFunction(xi)
title('Lebesgue function for '+str(n+1)+' Chebyshev points');
```

```{code-cell}
:tags: remove-input
n  = 29
xi = cos(linspace(0,pi,n+1))
PlotLebesgueFunction(xi)
title('Lebesgue function for '+str(n+1)+' Chebyshev points');
```

The values are small and grow slowly with $n$.
:::

+++

:::{exercise}
Write Python code to generate the plots of Lebesgue function as shown in above example.
:::

+++

:::{prf:theorem} Lebesgue constant
(1) The Lebesgue constants $\Lambda_n$ for degree $n \ge 0$ polynomial interpolation in any set of $n+1$ distinct points in $[-1,1]$ satisfy

$$
\Lambda_n \ge \frac{2}{\pi} \log(n+1) + C, \qquad C = \frac{2}{\pi}(\gamma + \log(4/\pi)) \approx 0.52125
$$

where $\gamma \approx 0.577$ is Euler's constant.

(2) For Chebyshev points, they satisfy

$$
\Lambda_n \le \frac{2}{\pi} \log(n+1) + 1
$$

(3) For equispaced points they satisfy

$$
\Lambda_n > \frac{2^{n-2}}{n^2}
$$
:::

See discussion after Theorem 15.2 in [@Trefethen2019] for a history of these results.

For $n=100$

* uniform points: $\Lambda_{100} > 3 \times 10^{25}$
* Chebyshev points: $\Lambda_{100} < 1.8$

+++

:::{prf:remark}
A large Lebesgue constant also means there is possibility of large errors in computation arising due to errors in data. If the function value is not known accurately and there is an error, then we are evaluating

$$
\tilde p(x) = \sum_{j=0}^m (f_j + \epsilon_j) \ell_j(x)
$$

instead of

$$
p(x) = \sum_{j=0}^m f_j \ell_j(x)
$$

The relative error

$$
\frac{|\tilde p(x) - p(x)|}{\max_j |\epsilon_j|} \le \sum_{j=0}^n |\ell_j(x)|
$$

and if the rhs is large, then there is possibility of large error in the evaluated polynomial.
:::

+++

[^1]: This is true for any normed vector space and is a consequence of
    $|\norm{x} - \norm{y}| \le \norm{x-y}$.

[^2]: Show that it is closed by showing all limit points belong to
    $V_n$.

[^3]: It is bounded because
    $$\norm{p} = \norm{p-f+f} \le \norm{p-f} + \norm{f} \le \norm{f} + \norm{f} =
    2\norm{f} < \infty$$
