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

# Spline functions

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

Consider a grid for the interval $[a,b]$

$$
a = x_0 < x_1 < \ldots < x_N = b
$$

:::{prf:definition}
We say that $s(x)$ is a *spline function* of order $m \ge 1$ if it
satisfies the following two properties.

1.  $s(x)$ is a polynomial of degree $< m$ on each sub-interval
    $[x_{i-1},x_i]$.

2.  $s^{(r)}(x)$ is continuous on $[a,b]$ for $0 \le r \le m-2$.
:::

The derivative of a spline of order $m$ is a spline of order $m-1$, and similarly for anti- derivatives. Note that $m=2$ corresponds to piecewise linear function which is continuous, and $m=3$ is a piecewise quadratic function which is continuous and has continuous derivatives [@Kincaid2002].

## Cubic spline interpolation

Consider the case of $m=4$ where $s(x)$ is a cubic polynomial in every
sub-interval $[x_{i-1},x_i]$ and we are given the values of some
function 

$$
y_i = f(x_i), \qquad i=0,1,\ldots,N
$$ 

Let us count the number
of unknowns and number of equations available. We can write, for
$i=1,2,\ldots,N$

$$
s(x) = a_i + b_i x + c_i x^2 + d_i x^3, \qquad x \in [x_{i-1},x_i]
$$

Hence we have $4N$ unknown coefficients $\{a_i,b_i,c_i,d_i\}$. The
spline function must satisfy some constraints; firstly, it must
interpolate the given function values

$$
s(x_i) = y_i, \qquad i=0,1,\ldots,N
$$ 

and its derivatives must be
continuous at the interior points

$$
s^{(r)}(x_i^-) = s^{(r)}(x_i^+), \qquad i=1,2,\ldots,N-1, \qquad r=0,1,2
$$

This gives us $(N+1)+3(N-1) = 4N-2$ equations, and we are short of two
equations since we have $4N$ unknown coefficients. There are different
ways in which the two missing equations can be supplied.

## Construction of cubic spline

Let 

$$
M_i = s''(x_i), \qquad i=0,1,\ldots,N
$$ 

Since $s(x)$ is cubic on $[x_i,x_{i+1}]$, then $s''(x)$ is linear and thus

$$
s''(x) = \frac{(x_{i+1}-x)M_i + (x-x_i)M_{i+1}}{h_i}, \qquad i=0,1,\ldots,N-1
$$

where $h_i = x_{i+1}-x_i$. By the above construction $s''(x)$ is continuous on $[a,b]$. Integrating twice, we get

$$
s(x) =& \frac{1}{6h_i}[ (x_{i+1}-x)^3 M_i + (x-x_i)^3 M_{i+1}] \\
      & + C_i (x_{i+1}-x) + D_i (x-x_i), \quad x \in [x_i,x_{i+1}]
$$ 

with $C_i,D_i,M_i$ to be
determined. Using the interpolation conditions $s(x_i)=y_i$ and
$s(x_{i+1}) = y_{i+1}$ gives

$$
C_i = \frac{y_i}{h_i} - \frac{h_iM_i}{6}, \qquad D_i = \frac{y_{i+1}}{h_i} - \frac{h_i
M_{i+1}}{6}, \qquad i=0,1,\ldots,N-1
$$ 

At this stage, the only unknowns
are the $M_i$ values. Note that $s(x)$ and $s''(x)$ are continuous on
$[a,b]$. We will now enforce continuity of $s'(x)$ at the interior
nodes. On $[x_i,x_{i+1}]$

$$
s'(x) = \frac{1}{2 h_i}[-(x_{i+1}-x)^2 M_i + (x-x_i)^2 M_{i+1}] + \frac{y_{i+1}-y_i}
{h_i} - \frac{(M_{i+1}-M_i)h_i}{6}
$$ 

and on $[x_{i-1},x_i]$

$$
s'(x) = \frac{1}{2 h_{i-1}}[-(x_{i}-x)^2 M_{i-1} + (x-x_{i-1})^2 M_{i}] + \frac{y_{i}-
y_{i-1}}{h_{i-1}} - \frac{(M_{i}-M_{i-1})h_{i-1}}{6}
$$ 

so that
continuity of derivative at $x=x_i$, $s'(x_i^-) = s(x_i^+)$, yields

$$
\frac{h_{i-1}}{6}M_{i-1} + \frac{h_{i-1}+h_i}{3}M_i + \frac{h_i}{6}M_{i+1} =
\frac{y_{i+1}-y_i}{h_i} - \frac{y_i - y_{i-1}}{h_{i-1}}, \quad i=1,2,\ldots,N-1
$$

We have $N-1$ equations for the $N+1$ unknowns $M_0, M_1, \ldots,M_N$.

## End-point derivative conditions

To obtain the two extra equations, we impose

$$
s'(x_0) = y_0', \qquad s'(x_N) = y_N'
$$ 

where $y_0'=f'(x_0)$, $y_N'=f'(x_N)$ are given to us. This yields the equations

$$
\begin{aligned}
\frac{h_0}{3} M_0 + \frac{h_0}{6} M_1 &= \frac{y_1-y_0}{h_0} - y_0' \\
\frac{h_{N-1}}{6} M_{N-1} + \frac{h_{N-1}}{3}M_N &= y_N' - \frac{y_N - y_{N-1}}
{h_{N-1}}
\end{aligned}
$$ 

We now have $N+1$ equations for the $N+1$ unknowns. They can be written in matrix form 

$$
AM = b
$$ 

where 

$$
b = \begin{bmatrix}
\frac{y_1-y_0}{h_0} - y_0' \\
\left\{ \frac{y_{i+1}-y_i}{h_i} - \frac{y_i - y_{i-1}}{h_{i-1}} \right\}_{i=1}^{N-1} \\
y_N' - \frac{y_N - y_{N-1}}{h_{N-1}}
\end{bmatrix}, \qquad
M = \begin{bmatrix}
M_0 \\ M_1 \\ \vdots \\ M_N
\end{bmatrix}
$$ 

$$
A = \left[
\begin{aligned}
& \frac{h_0}{3} \quad \frac{h_0}{6} \quad 0 \quad \cdots \\
& \textrm{diag}\left\{ \frac{h_{i-1}}{6}, \quad \frac{h_{i-1}+h_i}{3} \quad \frac{h_i}
{6} \right\}_{i=1}^{N-1} \\
& \qquad\quad\quad  \cdots \quad 0 \quad \frac{h_{N-1}}{6} \quad \frac{h_{N-1}}{3}
\end{aligned}
\right]
$$ 

The matrix $A$ is symmetric, positive definite and diagonally dominant. Hence the problem is uniquely solvable. The matrix is tridiagonal and can be solved rapidly using about $8N$ arithmetic operations (Thomas Algorithm). The resulting cubic spline is called the *complete cubic spline* interpolant, and we denote it by $s_c(x)$.

:::{prf:theorem}
Let $f \in \cts^4[a,b]$ and given a partition

$$
\tau_N: \quad a=x_0 < x_1 < \ldots < x_N = b
$$ 

define

$$
|\tau_N| := \max_{1 \le i \le N} (x_i - x_{i-1})
$$ 

Let $s_c(x)$ be the
complete cubic spline interpolant of $f(x)$ on the partition $\tau_N$

$$
\begin{aligned}
s_c(x_i) &= f(x_i), \qquad i=0,1,\ldots,N \\
s_c'(a) &= f'(a) \\
s_c'(b) &= f'(b)
\end{aligned}
$$ 

Then for some constants $c_j$

$$
\max_{a \le x \le b}|f^{(j)}(x) - s_c^{(j)}(x)| \le c_j |\tau_N|^{4-j} \max_{a \le x \le
b} |f^{(4)}(x)|
$$ 

for $j=0,1,2$. With the additional assumption

$$
\sup_N \frac{ |\tau_N|}{\min_{1 \le i \le N}(x_i - x_{i-1})} < \infty
$$

the above bound holds for $j=3$ also. Acceptable constants are

$$
c_0 = \frac{5}{384}, \qquad c_1 = \frac{1}{24}, \qquad c_2 = \frac{3}{8}
$$
:::

:::{prf:remark}
Taking $j=0$ we see that for a uniform grid of size $h = |\tau_N|$

$$
\max_{a \le x \le b}|f(x) - s_c(x)| \le \frac{5}{384} h^4 \max_{a \le x \le b} |f^{(4)}
(x)|
$$ 

The rate of convergence is 4; if $h$ is reduced to $h/2$, the
errors reduces by a factor of $1/16$.
:::

## Optimality property

The cubic spline gives approximations which have less oscillations
compared to other types of approximation. This is shown in the Python
example. We first prove the optimality result.

:::{prf:theorem}
Let $g \in \cts^2[a,b]$ and satisfies the interpolation conditions

$$
\begin{aligned}
g(x_i)  &= f(x_i), \quad i=0,1,\ldots,N \\
g'(x_0) &= f'(x_0) \\
g'(x_N) &= f'(x_N)
\end{aligned}
$$ 

Then

$$
\int_a^b |s_c''(x)|^2 \ud x \le \int_a^b |g''(x)|^2 \ud x
$$ 

with
equality iff $g(x) \equiv s_c(x)$. Thus $s_c$ oscillates least among all
smooth functions satisfying the interpolation conditions.
:::

:::{prf:proof}
Let 

$$
k(x) = s_c(x) - g(x)
$$ 

Lets compute the oscillations in $g$. 

$$
\begin{aligned}
\int_a^b |g''(x)|^2 \ud x 
&= \int_a^b [s_c''(x) - k''(x)]^2 \ud x \\
&= \int_a^b [s_c''(x)]^2 \ud x - 2 \int_a^b s_c''(x)  k''(x) \ud x + \int_a^b [k''(x)]^2
\ud x
\end{aligned}
$$ 

The second term is

$$
\int_a^b s_c''(x)  k''(x) \ud x = \sum_{i=0}^{N-1} \int_{x_i}^{x_{i+1}} s_c''(x) k''(x) \ud x
$$ 

For each integral on the right, perform integration by parts
twice 

$$
\begin{aligned}
\int_{x_i}^{x_{i+1}} s_c''(x) k''(x) \ud x 
&= [s_c'' k']_{x_i}^{x_{i+1}} - \int_{x_i} ^{x_{i+1}} s_c'''(x) k'(x) \ud x  \\
&= [s_c'' k']_{x_i}^{x_{i+1}} - [s_c''' k]_{x_i}^{x_{i+1}} + \int_{x_i}^{x_{i+1}} s_c''''(x) k(x) \ud x
\end{aligned}
$$ 

But $s_c''''=0$ since $s_c$ is a cubic polynomial and

$$
k(x_i) = s_c(x_i) - g(x_i) = f(x_i) - f(x_i) = 0, \qquad 0 \le i \le N
$$

Hence

$$
\int_a^b s_c''(x)  k''(x) \ud x = \sum_{i=0}^{N-1} [s_c'' k']_{x_i}^{x_{i+1}}
$$

Both $s_c''$ and $k'$ are continuous at $x_i$, $1 \le i \le N-1$. Hence
we get

$$
\int_a^b s_c''(x)  k''(x) \ud x = -s_c''(x_0) k'(x_0) + s_c''(x_N) k'(x_N)
$$

But 

$$
k'(x_0) = s_c'(x_0) - g'(x_0) = f'(x_0) - f'(x_0) = 0
$$ 

and
similarly $k'(x_N) = 0$, so that 

$$
\int_a^b s_c''(x)  k''(x) \ud x = 0
$$

This leads to

$$
\int_a^b |g''(x)|^2 \ud x = \int_a^b [s_c''(x)]^2 \ud x + \int_a^b [s_c''(x) - g''(x)]^2
\ud x \le \int_a^b [s_c''(x)]^2 \ud x
$$ 

We get equality iff

$$
s_c''(x) - g''(x) = 0, \qquad x \in [a,b]
$$ 

i.e., if $s_c(x) - g(x)$ is linear. The interpolation conditions then imply that $s_c(x) \equiv g(x)$.
:::

:::{prf:remark}
In the theorem, we can take $g=f$ and the theorem tells us that the
complete cubic spline $s_c$ that interpolates $f$ satisfies

$$
\int_a^b |s_c''(x)|^2 \ud x \le \int_a^b |f''(x)|^2 \ud x
$$
:::

:::{prf:example}
We will approximate the function
$$
f(x) = (1-x^2)^2 \sin(4\pi x) \exp(\sin(2\pi x)), \qquad x \in [-1,+1]
$$
We have taken an example for which
$$
f'(-1) = f'(+1) = 0
$$
which will be used in constructing the cubic spline below.

```{code-cell}
xmin, xmax = -1.0, +1.0
f = lambda x: (1-x**2)**2*sin(4*pi*x)*exp(sin(2*pi*x))

N = 20
x = cos(linspace(0,pi,N+1))
y = f(x)
xe = linspace(xmin,xmax,200)
ye = f(xe)
plot(xe,ye,x,y,'o')
legend(('f(x)','Data'));
```

Let us first try polynomial interpolation using Chebyshev points.


```{code-cell}
p = polyfit(x,y,N)
yp = polyval(p,xe)
plot(xe,ye,x,y,'o',xe,yp,'r-')
legend(('Exact','Data','Polynomial'));
```

We see that polynomial interpolation is rather oscillatory. Of course this will become better as we increase the number of samples. Using uniformly spaced points for polynomial interpolation would have given us even bad approximation. We will next try cubic spline interpolation. First we implement a function which constructs the matrix and right hand side.

```{code-cell}
def spline_matrix(x,y):
    N = len(x) - 1
    h = x[1:] - x[0:-1]
    A = zeros((N+1,N+1))
    A[0][0] = h[0]/3; A[0][1] = h[0]/6
    for i in range(1,N):
        A[i][i-1] = h[i-1]/6; A[i][i] = (h[i-1]+h[i])/3; A[i][i+1] = h[i]/6
    A[N][N-1] = h[N-1]/6; A[N][N] = h[N-1]/3

    r = zeros(N+1)
    r[0] = (y[1] - y[0])/h[0]
    for i in range(1,N):
        r[i] = (y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1]
    r[N] = -(y[N]-y[N-1])/h[N-1]
    return A,r
``` 

We will use uniformly spaced points for spline interpolation. We first construct the matrix and right hand side and solve the matrix problem. Finally, we evaluate the spline in each interval to make a plot.

```{code-cell}
x = linspace(xmin,xmax,N+1)
y = f(x)

A,r = spline_matrix(x,y)
m = solve(A,r)

plot(xe,ye,x,y,'o')
for i in range(0,N):
    h  = x[i+1] - x[i]
    xs = linspace(x[i],x[i+1],6)
    ys = ((x[i+1]-xs)**3*m[i] + (xs-x[i])**3*m[i+1])/(6*h) \
         + ((x[i+1]-xs)*y[i] + (xs-x[i])*y[i+1])/h \
         - ((x[i+1]-xs)*m[i] + (xs-x[i])*m[i+1])*h/6
    plot(xs,ys,'r-')

legend(('Exact','Data','Spline'));
```

We see that spline gives a better approximation and does not have too much oscillation.
:::

## Not-a-knot condition

If the end derivative conditions $f'(a)$, $f'(b)$ are not known, then we
need other conditions to complete the system of equations. We can demand
that $s'''(x)$ is continuous at $x_1$ and $x_{N-1}$. This is equivalent
to requiring that $s(x)$ be a cubic spline function with knots

$$
\{ x_0, x_2, x_3, \ldots, x_{N-2}, x_N\}
$$ 

while still interpolating at all node points 

$$
\{x_0,x_1,x_2,\ldots,x_{N-1},x_N\}
$$

This again leads to a tri-diagonal linear system of equations.

## Natural cubic spline

In the natural cubic spline approximation, the two end-point conditions

$$
s''(x_0) = s''(x_N) = 0
$$ 

are used. This leads to a tridiagonal system of equations. The cubic spline approximations converge slowly near the end-points. The natural cubic spline also has an optimality property, for proof, see the problems collection.

