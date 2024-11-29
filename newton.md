---
exports:
  - format: pdf
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
numbering:
  code: false
  math: false
  headings: true
math:
  '\re': '\mathbb{R}'
  '\df': '\frac{\partial #1}{\partial #2}'
  '\half': '\frac{1}{2}'
  '\shalf': '\tfrac{1}{2}'
  '\ud': '\textrm{d}'
  '\dd': '\frac{\ud #1}{\ud #2}'
  '\sign': '\textrm{sign}'
  '\cts': 'C'
  '\ii': '\mathfrak{i}'
---

# Newton method

+++

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
import sympy
```

The Newton-Raphson method is an iterative method to find the roots of a         function. It starts with an initial guess for the root and tries to improve it. Given an     initial guess $x_0$ for the root, we have

$$
f(x_0) \ne 0
$$

Let us try to find a correction $\Delta x_0$ which brings it closer to the root

$$
f(x_0 + \Delta x_0) = 0
$$

Using Taylor formula, and truncating it at the linear term

$$
f(x_0) + \Delta x_0 f'(x_0) \approx 0 \qquad \Longrightarrow \qquad \Delta x_0  \approx -
\frac{f(x_0)}{f'(x_0)}
$$

We obtain a {\em better} estimate for the root

$$
x_1 = x_0 + \Delta x_0 = x_0 - \frac{f(x_0)}{f'(x_0)}
$$

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/asy/newton.svg
:width: 60%
:align: center
```

Geometrically, we are approximating the function by the tangent line at $x_0$ which is given by

$$
y = f(x_0) + (x-x_0) f'(x_0)
$$

and finding the zero of the tangent line, see figure. We now repeat the same process for $x_1, x_2, \ldots$ until some convergence is achieved. The iterations are given by

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}, \qquad n=0,1,2,\ldots
$$

This algorithm requires that $f$ is differentiable and $f'$ is not zero in a    neighbourhood of the root. The algorithm for Newton-Raphson method is illustrated in Algorithm.

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/latex/p3.svg
:width: 100%
:align: center
```

:::{prf:example} Affine function
Consider an affine function $f(x) = ax + b$. Let $x_0$ be   any initial
guess. The Newton iterations are yield

$$
x_1 = x_0 - \frac{ax_0 + b}{a} = -\frac{b}{a}
$$

so that we converge to the root in one iteration.
:::

:::{prf:example}
Find the root of

$$
f(x) = \exp(x) - \frac{3}{2} - \arctan{x}
$$

which is in the interval $[\shalf,1]$.

We define the function using sympy and automatically compute its derivative. Then we make functions out of these.

```{code-cell}
x = sympy.Symbol('x')

# Define the function
f = sympy.exp(x)-3.0/2.0-sympy.atan(x)

# Compute its derivative
df = sympy.diff(f,x)
print("f = ", f)
print("df= ", df)

# Make functions
F  = sympy.lambdify(x, f, modules=['numpy'])
DF = sympy.lambdify(x, df, modules=['numpy'])

# Plot the functions
X = linspace(-1,1,100)
figure(figsize=(9,4))
subplot(1,2,1), plot(X,F(X)),  grid(True), title('Function')
subplot(1,2,2), plot(X,DF(X)), grid(True), title('Derivative');
```

Now we implement the Newton method. We have to specify the maximum number of iterations, an initial guess and a tolerance to decide when to stop.

```{code-cell}
M   = 20    # maximum iterations
x   = 0.5   # initial guess
eps = 1e-14 # relative tolerance on root

f = F(x)
for i in range(M):
    df = DF(x)
    dx = -f/df
    x  = x + dx
    e  = abs(dx)
    f  = F(x)
    print("%6d %22.14e %22.14e %22.14e" % (i,x,e,abs(f)))
    if e < eps * abs(x):
        break
```

:::

+++

:::{exercise}
Modify the above Newton method implemention to use while loop.
:::

+++

:::{prf:example} Absolute and relative tolerance
We consider the same function as in a previous example, but we have shifted it so that the root lies at a large distance from zero,

$$
f(x) = \exp(x-x_0) - \frac{3}{2} - \arctan(x-x_0)
$$

where $x_0$ is some large number, e.g., $x_0 = 10^{10}$. The root lies around $x_0$ and due to precision of floating point numbers, we cannot compute it to good absolute tolerance, but only to relative tolerance.

Define the function and its derivative.

```{code-cell}
x0 = 1e10

# Define the function
F = lambda x: exp(x-x0)-3.0/2.0-atan(x-x0)

# Compute its derivative
DF = lambda x: exp(x-x0) - 1.0/(1 + (x-x0)**2)
```

Now we implement the Newton method. We have to specify the maximum number of iterations, an initial guess and a tolerance to decide when to stop.

```{code-cell}
M   = 50    # maximum iterations
eps = 1e-14 # relative tolerance on root
```

**Using relative tolerance on the root**

```{code-cell}
x = x0 + 0.5   # initial guess
f = F(x)
for i in range(M):
    df = DF(x)
    dx = -f/df
    x  = x + dx
    e  = abs(dx)
    f  = F(x)
    print("%6d %22.14e %22.14e %22.14e" % (i,x,e,abs(f)))
    if e < eps * abs(x):
        break
```

Root converges in relative sense to machine zero in few iterations.

**Using absolute tolerance on the root**


```{code-cell}
x = x0 + 0.5   # initial guess
f = F(x)
for i in range(M):
    df = DF(x)
    dx = -f/df
    x  = x + dx
    e  = abs(dx)
    f  = F(x)
    print("%6d %22.14e %22.14e %22.14e" % (i,x,e,abs(f)))
    if e < eps:
        break
```

There is no convergence. The root cannot be computed to good absolute accuracy, and consequently, we cannot make $f(x)$ sufficiently small. We could reduce it to about $10^{-6}$ which is also achieved using relative tolerance as stopping criterion.

:::

+++

## Convergence analysis

Let $r$ be a root of $f$. From Taylor formula with remainder term

$$
0 = f(r) = f(x_n) + (r - x_n) f'(x_n) + \half (r - x_n)^2 f''(\xi_n), \qquad    \xi_n \mbox{ between $r$ and $x_n$}
$$

Solve for $r$ from the linear term

$$
r = \underbrace{x_n - \frac{f(x_n)}{f'(x_n)}}_{x_{n+1}} - \half (r - x_n)^2
\frac{f''(\xi_n)}{f'(x_n)}
$$

The error in the root is

$$
r - x_{n+1} = - (r - x_n)^2 \frac{f''(\xi_n)}{2 f'(x_n)}
$$

Atleast in a neighborhood of the root, we will assume that

$$
\left| \frac{f''(\xi_n)}{2 f'(x_n)} \right| \le C
$$

Then the error $e_n = |r - x_n|$ satisfies

$$
e_{n+1} \le C e_n^2
$$

To show convergence $x_n \to r$, we have to show that $e_n \to 0$.

+++

:::{prf:theorem} Convergence of Newton method
Assume that $f,f',f''$ are continuous in some neighbourhood of the root $r$,    and  $f'(r) \ne 0$. If $x_0$ is chosen sufficiently close to $r$, then $x_n \to r$. Moreover

$$
\lim_{n \to \infty} \frac{r - x_{n+1}}{(r - x_n)^2} = - \frac{f''(r)}{2f'(r)}
$$

:::

+++

:::{prf:proof}
Pick a sufficiently small interval $I = [r-\delta,r+\delta]$ on which $f' \ne   0$. This is possible to do since $f'$ is continuous and $f'(r) \ne 0$. Define

$$
M = \frac{\max_{x \in I} |f''(x)|}{\min_{x \in I} |f'(x)|}
$$

We have

$$
|r - x_1| \le M |r-x_0|^2 \qquad \Longrightarrow \qquad M |r - x_1| \le (M|r -  x_0|)^2
$$

Choose $x_0$ such that

$$
|r - x_0| \le \delta \qquad\mbox{and}\qquad M|r - x_0| < 1
$$

$$
\Longrightarrow \qquad M|r - x_1| < 1
$$

Moreover

$$
|r - x_1| \le \underbrace{M |r - x_0|}_{< 1} \cdot |r - x_0| < |r - x_0| \le    \delta
$$

Hence $x_1 \in I$. Applying the same argument inductively to $x_2, x_3, \ldots$ we obtain

$$
\forall n \ge 0 \qquad \begin{cases}
|r - x_n| \le \delta \qquad \Longrightarrow \qquad x_n \in I \\
M|r - x_n| < 1
\end{cases}
$$

To show convergence, we write

$$
|r - x_{n+1}| \le M |r - x_n|^2 \qquad \Longrightarrow \qquad M|r - x_{n+1}|    \le (M |r- x_n|)^2
$$

and inductively

$$
M|r - x_n| \le (M |r - x_0|)^{2^n}
$$

Now taking limit

$$
\lim_{n \to \infty} |r - x_n| \le \frac{1}{M} \lim_{n \to \infty} (M |r - x_0|)^{2^n} = 0
$$

since $M|r-x_0| < 1$, which proves the convergence. In the error formula

$$
r - x_{n+1} = - (r - x_n)^2 \frac{f''(\xi_n)}{2 f'(x_n)}
$$

$\xi_n$ is between $r$ and $x_n$ and hence $\xi_n \to r$. Thus

$$
\lim_{n \to \infty} \frac{r - x_{n+1}}{(r - x_n)^2} = -\lim_{n \to \infty}
\frac{f''(\xi_n)}{2 f'(x_n)} = - \frac{f''(r)}{2f'(r)}
$$

:::

+++

:::{prf:remark}
The above proof shows that $x_0$ must be close to $r$ such   that

$$
|r - x_0| < \frac{1}{M}
$$

In practice we may not be able to have a good estimate of $M$.
:::

+++

## Error estimate

From mean value theorem

$$
f(x_n) = f(x_n) - f(r) = f'(\xi_n) (x_n - r)
$$

where $\xi_n$ is between $x_n$ and $r$.

$$
r - x_n = - \frac{f(x_n)}{f'(\xi_n)}
$$

If $f'$ is not changing rapidly between $x_n$ and $r$ then

$$
f'(\xi_n) \approx f'(x_n)
$$

We get the {\em error estimate}

$$
r - x_n \approx - \frac{f(x_n)}{f'(x_n)} = x_{n+1} - x_n
$$

and the {\em relative error estimate}

$$
\frac{r - x_n}{r} \approx \frac{x_{n+1} - x_n}{x_{n+1}}
$$

Hence we can use the test

$$
|x_{n+1} - x_n| < \epsilon |x_{n+1}|, \qquad 0 < \epsilon \ll 1
$$

to decide on convergence. This is used in the code examples to stop the         iterations.

+++

:::{prf:theorem} Convex functions
Assume that $f \in \cts^2(\re)$, is increasing, convex and has a zero. Then the zero is unique and the Newton method will converge to it starting from any initial      point.
:::

+++

:::{prf:proof}
Let $r$ be a zero of $f$. The function cannot have more than one zero since for $x > r$

$$
f(x) = f(r) + \int_r^x f'(s) \ud s = \int_r^x f'(s) \ud s > 0
$$

unless $f' \equiv 0$.

We have $f' >0$ and $f'' > 0$.

$$
r - x_{n+1} = - (r - x_n)^2 \frac{f''(\xi_n)}{2 f'(x_n)} \le 0 \qquad           \Longrightarrow
\qquad x_{n+1} \ge r \qquad \forall n \ge 0
$$

Also $f(x_n) > 0$ (WHY) for $n=1,2,\ldots$ and hence

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} < x_n, \qquad n=1,2,\ldots
$$

$\{x_n\}_{n\ge 1}$ is a decreasing sequence bounded below by $r$. Hence $x^* =
\lim_{n \to \infty} x_n$ exists and

$$
x^* = x^* - \frac{f(x^*)}{f'(x^*)} \qquad \Longrightarrow \qquad f(x^*) = 0     \qquad
\Longrightarrow \qquad x^* = r
$$

:::

+++

:::{prf:example} Square root

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/asy/square_root.svg
:width: 60%
:align: center
```

For $a > 0$, the function

$$
f(x) = x^2 - a
$$

has roots $\pm\sqrt{a}$. Newton's method is

$$
x_{k+1} = \half \left( x_k + \frac{a}{x_k} \right)
$$

which uses only basic arithmetic operations (but division is not very basic).   Geometrically, it is easy to see that if $x_0 > 0$, then we converge to the positive root and  if $x_0 < 0$, then we converge to the negative root.

```{code-cell}
from math import sqrt
a = 2.0
r = sqrt(a)
```

Define the function and its derivative.

```{code-cell}
def f(x):
    return x**2 - a

def df(x):
    return 2.0*x
```

Now implement Newton method as a function.

```{code-cell}
def newton(f,df,x0,M=100,eps=1.0e-15):
    x = x0
    for i in range(M):
        dx = - f(x)/df(x)
        x  = x + dx
        print('%3d %22.14f %22.14e %22.14e'%(i,x,x-r,abs(f(x))))
        if abs(dx) < eps * abs(x):
            return x
    print('No convergence, current root = %e' % x)
```

We call Newton method with a complex initial guess for the root.

```{code-cell}
x0 = 2.0
x = newton(f,df,x0)
```

From the last column, we observe that in each iteration, the number of correct decimal places doubles.
:::

+++

## Implication of quadratic convergence
What does quadratic convergence imply in practical terms ? For Newton method

$$
e_{k+1} \le M e_k^2
$$

Suppose $M \approx 1.0$ and $e_0 = 10^{-1}$. Then

\begin{align*}
& e_1 \approx 10^{-2} \\
& e_2 \approx 10^{-4} \\
& e_3 \approx 10^{-8} \\
& e_4 \approx 10^{-16} \qquad \textrm{(limit of double precision)}
\end{align*}

If $\alpha = O(1)$, then

\begin{align*}
& x_1 \textrm{ accurate to 2 digits} \\
& x_2 \textrm{ accurate to 4 digits} \\
& x_3 \textrm{ accurate to 8 digits} \\
& x_4 \textrm{ accurate to 16 digits}
\end{align*}

Newton's method doubles the number of accurate figures in each iteration. Of    course this happens only after the method enters the quadratic convergence regime. In the   first few iterations, there may be no reduction of error or less than quadratic reduction.

\fbox{Show quadratic convergence in square root example}

+++

:::{prf:example} Reciprocal

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/asy/reciprocal.svg
:width: 60%
:align: center
```
Given a real number $a > 0$, we want to find its reciprocal. Since division     cannot be performed exactly on a computer, we would like to find a method which uses only addition and multiplication. We can cook up a function whose root is the reciprocal, e.g.,

$$
f(x) = \frac{1}{x} - a
$$

Applying Newton method to this function yields the iteration, see
figure.

$$
x_{k+1} = 2 x_k - a x_k^2
$$

which involves only subtraction and multiplication. If we start with any $0 <   x_0 \le \frac{1}{a}$, then we converge to the root as we can verify geometrically. If we start with $x_0 > \frac{1}{a}$, but sufficiently close to it, then we again     converge to the root.  If $a > 1$ then $x_1 \in (0, 1/a]$ iff $x_0 \in [1/a,2)$, and        subsequent iterates will converge. If $a>1$ and $x_0 > 2$, then $x_1 < 0$ after which we   have divergence to $-\infty$.

If $0 < a < 1$, then any initial guess in the interval $(0,1/a] \cup [1/a, 2/   a)$ will converge to the desired root. Note that $x_0 = a$ is a valid initial guess in   this case.

Let us try for $a = 10^{-10}$ and with an initial guess $x_0 = a$. 

```{code-cell}
a = 1.0e-10

def F(x):
    return 1/x - a

def DF(x):
    return -1/x**2
```

Now implement Newton method.

```{code-cell}
M   = 100      # maximum iterations
x   = a     # initial guess
eps = 1e-15 # relative tolerance on root

f = F(x)
for i in range(M):
    df = DF(x)
    dx = -f/df
    x  = x + dx
    e  = abs(dx)
    f = F(x)
    print("%6d %22.14e %22.14e %22.14e" % (i,x,e,abs(f)))
    if e < eps * abs(x):
        break
```
:::

+++

:::{prf:example} Slow convergence
The quadratic convergence of Newton method is obtained only when we are         sufficiently close to the root. It may happen that initially, the convergence is very slow.  Consider the reciprocal example with $a = 10^{-10}$ and let us start with $x_0 = a$ for which the method converges but this initial guess is too far from the root. The first few iterates are

$$
\begin{aligned}
x_1 &= 2 \cdot 10^{-10} - 10^{-30} = 2 \cdot 10^{-10} = 2 x_0 \\
x_2 &= 2 \cdot 2 \cdot 10^{-10} - 4 \cdot 10^{-30} = 2^2 \cdot 10^{-10} = 2 x_1
\end{aligned}
$$

Thus initially, the iterates double in each step and the change per step is     very small. This slow change will continue until we get close to the root, i.e.,

$$
2^k \cdot 10^{-10} \approx 10^{10}
$$

which needs  $k \approx 67$ iterations, after which we have rapid convergence.

:::

+++

## Inflection point

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/asy/inflection.svg
:width: 60%
:align: center
```

If we have a function with an inflection point, see figure,  then the Newton iterates may cycle indefinitely without any convergence. The Newton      method in step for is

$$
\Delta x_k = - \frac{f(x_k)}{f'(x_k)}, \qquad x_{k+1} = x_k + \Delta x_k
$$

We can use $|f(x)|$ as a metric to measure convergence towards the root:        $|f(x_k)|$ must decrease monotonically towards zero. We will accept the step $\Delta x_k$  if $| f(x_{k+1})| < |f(x_k)|$, otherwise halve the step size $\Delta x_k \leftarrow   \half \Delta x_k$. We can halve the step size repeatedly until reduction in $|f(x)|$  is obtained; such a step size always exists. Moreover, to prevent wild oscillations of the iterates, the step size can also be restricted, e.g., dont allow step size to   more than double in successive steps

$$
|\Delta x_k| \le 2 |\Delta x_{k-1}|
$$

+++

:::{prf:example} Complex roots
Newton method can converge to complex root if we start with a complex initial guess. The roots of

$$
f(x) = x^5 + 1
$$

lie on the unit circle in the complex plane. We can write

$$
-1 = \exp(\ii (2n+1)\pi), \qquad n=0,1,2,\ldots
$$

The fifth roots are

$$
(-1)^{1/5} = \exp(\ii (2n+1)\pi/5), \qquad n=0,1,2,3,4
$$

Let us compute the exact roots and plot them on a graph.

```{code-cell}
n = array([0,1,2,3,4])
r = exp(1j*(2*n+1)*pi/5) # 5'th roots of -1
plot(real(r),imag(r),'o')
grid(True), axis('equal')
for x in r:
    print(x,x**5)
```

Let us find a root using Newton method with starting guess $1 + \ii$.

```{code-cell}
f  = lambda x: x**5 + 1.0
df = lambda x: 5.0 * x**4

def newton(f,df,x0,M=100,eps=1.0e-15):
    x = x0
    for i in range(M):
        dx = - f(x)/df(x)
        x  = x + dx
        print(i,x,abs(f(x)))
        if abs(dx) < eps * abs(x):
            return x

x0 = 1.0 + 1.0j
x = newton(f,df,x0)
print('Root = ',x)
print('x**5 = ',x**5)
```

We do converge to a root, but which one will depend on the initial guess. Moreover, we can only get one root.
:::

+++

## Implicit functions

Suppose we have a function $y=y(x)$ which is defined implicitly by

$$
G(x,y) = 0
$$ 

We may not be able find an explicit expression $y(x)$, instead we find the value of $y$ for some discrete set of $x$ values.  This leads to a root finding problem for each $x$ and we can apply Newton method.

Take a value of $x$ and make an initial guess $y_0$. This may not satisfy $G(x,y_0) = 0$ so we want to find a correction $y_0 + \Delta y_0$ such that $G(x,y_0+\Delta y_0) = 0$. By Taylor expansion

$$
G(x,y_0) + G_y(x,y_0) \Delta y_0 = 0 \quad\Longrightarrow\quad \Delta y_0 = -
\frac{G(x,y_0)}{G_y(x,y_0)}
$$ 

and hence we get an improved estimate

$$
y_1 = y_0 + \Delta y_0 = y_0 - \frac{G(x,y_0)}{G_y(x,y_0)}
$$ 

This can be used to generate an iterative scheme

$$
y_{k+1} = y_k - \frac{G(x,y_k)}{G_y(x,y_k)}, \qquad k=0,1,2,\ldots
$$

If we have to solve this problem for several values of $x$, then we can use the previous solution as an initial guess for the root, which can speed up the convergence of Newton method.

## System of equations

Consider a function $f : \re^N \to \re^N$ such that for $x = [x_1, x_2, \ldots, x_N]^\top$ 

$$
f(x) = [f_1(x), \ldots, f_N(x)]^\top \in \re^N
$$ 

We want to find an $\alpha \in \re^N$ such that $f(\alpha) = 0$. This is a system of $N$ non-linear equations for the $N$ unknowns. We start with an initial guess $x^0 \in \re^N$ and in general $f(x^0) \ne 0$. We will try to find a correction $\Delta x^0$ to $x^0$ such that 

$$
f(x^0 + \Delta x^0) = 0
$$

Approximating by the first two terms of the Taylor expansion

$$
f(x^0) + f'(x^0) \Delta x^0 = 0
$$ 

we can find the correction by solving the matrix problem 

$$
f'(x^0) \Delta x^0 = - f(x^0)
$$ 

where $f' \in \re^{N \times N}$ is given by 

$$
[f']_{ij} = \df{f_i}{x_j}
$$ 

Once we solve for $\Delta x^0 = -[f'(x^0)]^{-1} f(x^0)$, we get a new estimate of the root

$$
x^1 = x^0 + \Delta x^0 = x^0 - [f'(x^0)]^{-1} f(x^0)
$$ 

The above steps are applied iteratively. Choose an initial guess $x^0$ and set $k=0$.

1.  Solve the matrix problem 
$$
f'(x^k) \Delta x^k = - f(x^k)
$$

2.  Update the root 
$$
x^{k+1} = x^k + \Delta x^k
$$

3.  $k=k+1$

+++

:::{prf:example} Newton method for system
For $x = (x_0,x_1,x_2) \in \re^3$, define $f : \re^3 \to \re^3$ by

$$
f(x) = \begin{bmatrix}
x_0 x_1 - x_2^2 - 1 \\
x_0 x_1 x_2 - x_0^2 + x_1^2 - 2 \\
\exp(x_0) - \exp(x_1) + x_2 - 3
\end{bmatrix}
$$

The gradient of $f$ is given by

$$
f'(x) = \begin{bmatrix}
x_1 & x_0 & -2 x_2 \\
x_1 x_2 - 2 x_0 & x_0 x_2 + 2 x_1 & x_0 x_1 \\
\exp(x_0) & -\exp(x_1) & 1
\end{bmatrix}
$$

We implement functions for these.

```{code-cell}
def f(x):
    y = zeros(3)
    y[0] = x[0]*x[1] - x[2]**2 - 1.0
    y[1] = x[0]*x[1]*x[2] - x[0]**2 + x[1]**2 - 2.0
    y[2] = exp(x[0]) - exp(x[1]) + x[2] - 3.0
    return y

def df(x):
    y = array([[x[1],               x[0],               -2.0*x[2]],
               [x[1]*x[2]-2.0*x[0], x[0]*x[2]+2.0*x[1], x[0]*x[1]],
               [exp(x[0]),          -exp(x[1]),               1.0]])
    return y
```

The Newton method uses matrix solver from Numpy.

```{code-cell}
def newton(fun,dfun,x0,M=100,eps=1.0e-14,debug=False):
    x = x0
    for i in range(M):
        g = fun(x)
        J = dfun(x)
        h = solve(J,-g)
        x = x + h
        if debug:
            print(i,x,norm(g))
        if norm(h) < eps * norm(x):
            return x
```

Now we solve

```{code-cell}
x0 = array([1.0,1.0,1.0])
x = newton(f,df,x0,debug=True)
print('Root = ',x)
print('f(x) = ',f(x))
```
