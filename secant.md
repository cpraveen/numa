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

# Secant Method

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## Quasi-Newton methods

Newton method requires the knowledge of the gradient. In many situations, it may not be possible to find the gradient. E.g., if $f(x)$ is available only as a computer code, then it may not be possible to compute the gradient[^1]. In such situations we can try to approximate the gradient and methods which make use of approximate derivatives are called *quasi-Newton methods*. For a scalar function $f(x)$, we can approximate the derivative using a finite difference approximation, e.g.,

$$
f'(x_k) \approx \frac{f(x_k+h) - f(x_k)}{h}, \qquad 0 \le h \ll 1
$$

Since we want to generate a sequence of approximations to the root, we can use two successive values of $x_k$ to estimate the derivative

$$
f'(x_k) \approx \frac{f(x_k) - f(x_{k-1})}{x_k - x_{k-1}}
$$ 

Using this derivative in place of the exact derivative in Newton method leads to the *secant method*.

## Secant method

Given two starting values $x_0$, $x_1$ which are sufficiently close to the root $r$, we construct a straight line approximation to $f(x)$ passing through the points $(x_0,f(x_0))$ and $(x_1, f(x_1))$. We find the zero of this straight line which will form the next approximation to the root $x_2$ and is given by

$$
x_2 = x_1 - f(x_1) \frac{x_1 - x_0}{f(x_1) - f(x_0)}
$$ 

Inductively, we obtain

$$
x_{k+1} = x_k - f(x_k) \frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}, \qquad k=1,2,\ldots
$$

Note that this is just Newton iteration formulae where the derivative has been approximated by the finite difference formula.

:::{prf:example}
Compute the zero of 

$$
f(x) = \exp(x) - \frac{3}{2} - \arctan(x)
$$ 

using secant method. The function looks like this.

```{code-cell}
f = lambda x: exp(x) - 3.0/2.0 - arctan(x)
x = linspace(0.0, 1.0, 100)
plot(x,f(x)), grid(True), xlabel('x'), ylabel('f(x)');
```

Let us implement the secant method.

```{code-cell}
M = 50       # maximum number of iterations
eps = 1.0e-6 # tolerance for stopping
x0, x1 = 0.5, 0.6

for i in range(M):
    f0, f1 = f(x0), f(x1)
    x2 = x1 - f1*(x1 - x0)/(f1 - f0)
    if abs(x2-x1) < abs(x2)*eps:
        break
    print("%3d %18.10e %18.10e %18.10e %18.10e" % (i, x0, x1, f0, f1))
    x0, x1 = x1, x2

print("Number of iterations = ", i)
print("Final solution = ", x2, f(x2))
```

:::

## Convergence analysis

From the iteration formula, we get by subtracting $r$ on both sides

$$
r - x_{k+1} = r - x_k + f(x_k) \frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}
$$

Define the Newton divided differences

$$
f[x_{k-1},x_k] := \frac{f(x_k) - f(x_{k-1})}{x_k - x_{k-1}}, \qquad f[x_{k-1},x_k,
r]:= \frac{f[x_k,r] - f[x_{k-1},x_k]}{r - x_{k-1}}
$$ 

then it is easy to show that

$$
f[x_{k-1},x_k] = f'(y_k), \qquad \textrm{for some $y_k$ between $x_{k-1}, x_k$}
$$

and

$$
f[x_{k-1},x_k,r] = \half f''(z_k) \qquad \textrm{for some $z_k$ between
$x_{k-1},x_k,r$}
$$ 

Hence we get the error formula

\begin{align}
r - x_{k+1} &= -(r - x_{k-1})(r - x_k) \frac{f[x_{k-1},x_k,r]}
{2f[x_{k-1},x_k]} \\
&= -(r - x_{k-1})(r - x_k)\frac{f''(z_k)}{2f'(y_k)}
\end{align}

If we can bound the last factor

$$
\left| \frac{f''(z_k)}{2f'(y_k)} \right| \le M < \infty
$$ 

then we get a bound on the error $e_k = |r - x_k|$

$$
e_{k+1} \le M e_{k-1} e_k
$$

:::{prf:theorem}
Assume that $f(x)$, $f'(x)$, $f''(x)$ are continuous for all values of $x$ in some interval around $r$, and that $f(r)=0$, $f'(r) \ne 0$. Then if the initial guesses $x_0$, $x_1$ are sufficiently close to $r$, then the iterates of the secant method will converge to $r$; the order of convergence is

$$
p = \half (1 + \sqrt{5}) \approx 1.62
$$
:::

Recall that Newton method converges with order $p=2$, so we can expect secant method to be bit slower than Newton.

:::{prf:proof}
(1) Pick a sufficiently small interval $I = [r-\delta, r+\delta]$ on which $f'(x) \ne 0$. This is possible to do since $f'(r) \ne 0$ and $f'$ is
continuous. Then define

$$
M = \frac{\max_{x \in I} |f''(x)|}{2\min_{x \in I} |f'(x)|} < \infty
$$

For any $x_0, x_1 \in I$, the errors in secant method satisfy

$$
\begin{aligned}
e_2 &\le e_1 \cdot M e_0 \\
M e_2 &\le Me_1 \cdot Me_0
\end{aligned}
$$ 

Further assume that $x_0$, $x_1$ are chosen such that

$$
\Delta = \max\{ Me_0 , Me_1\} < 1
$$ 

Then 

$$
Me_2 \le \Delta^2 < 1
$$

Also 

$$
M e_2 \le \Delta^2 < \Delta
$$ 

$$
\begin{aligned}
& \implies e_2 < \frac{\Delta}{M} = \max\{e_0, e_1\} \le \delta \\
& \implies x_2 \in I
\end{aligned}
$$ 

Applying the above argument inductively, we get 

$$
\left.
\begin{aligned}
x_k & \in I \\
M e_k & \le \Delta < 1
\end{aligned}
\right\} \quad \textrm{for } k = 0,1,2,3,\ldots
$$ 

To prove convergence, we apply the error inequality recursively. 

$$
\begin{aligned}
M e_0 & \le \Delta \\
M e_1 & \le \Delta \\
M e_2 & \le Me_1 \cdot Me_0 \le \Delta^2 \\
M e_3  & \le M e_2 \cdot M e_1 \le \Delta^3 \\
M e_4 & \le M e_3 \cdot M e_2 \le \Delta^5 \\
M e_5 & \le M e_4 \cdot M e_3 \le \Delta^8
\end{aligned}
$$ 

In general 

$$
M e_k \le \Delta^{q_k}
$$ 

with

$$
q_0 = q_1 = 1, \qquad q_k = q_{k-1} + q_{k-2}
$$ 

We see that $q_k$ is a Fibonacci sequence. Since $\Delta < 1$ and $q_k \to \infty$ as $k \to \infty$, we get

$$
\lim_{k \to \infty} e_k \le \frac{1}{M} \Delta^{q_k} = 0
$$ 

Hence the secant method converges to some root.

(2) Now let us estimate the convergence rate $p$. Assume that

$$
e_k \approx A e_{k-1}^p
$$ 

and we will try to estimate $p$. Now

$$
e_{k-1} \approx \left( \frac{e_k}{A} \right)^{1/p}
$$ 

Using the error relation for secant method

\begin{align}
e_{k+1} & \approx M e_k \cdot e_{k-1} \approx M e_k \cdot \frac{1}{A^{1/p}} e_k^{1/p} \\
A e_k^p & \approx M A^{-1/p} e_k^{1 + 1/p}
\end{align}

This implies that

$$p = 1 + \frac{1}{p} \qquad \textrm{and} \qquad A \approx M A^{-1/p}$$

The first equation $p^2 -p - 1 = 0$ has the roots

$$p = \frac{1 \pm \sqrt{5}}{2}$$ 

We choose the positive root since we need $p > 1$, hence the convergence rate is

$$p = \frac{1 + \sqrt{5}}{2} \approx 1.618$$ 

This is less than quadratic convergence, but still better than linear convergence. From the second relation we get

$$A \approx M^{p/(1+p)} \approx \left| \frac{ f''(r)}{2f'(r)} \right|
^{\frac{\sqrt{5} - 1}{2}}$$
:::

:::{prf:remark}
We see that the power $q_k$ in the above proof is a Fibonacci sequence whose $k$'th term can be written as

$$
q_k = \frac{1}{\sqrt{5}}[ r_0^{k+1} - r_1^{k+1}], \qquad k \ge 0
$$

where

$$
r_0 = \half(1 + \sqrt{5}) \approx 1.618, \qquad r_1 = \half (1 - \sqrt{5}) \approx -0.618
$$

are the roots of $p^2 - p - 1 = 0$.
:::

## Error estimate

We can use the function value to decide on convergence but this may not be a good criterion if the function is very flat around the root.  Additionally we can estimate the error in the root also. From mean value theorem 

$$
f(x_k) = f(x_k) - f(r) = f'(\xi_k) (x_k - r)
$$ 

where $\xi_k$ is between $x_k$ and $r$. 

$$r - x_k = - \frac{f(x_k)}{f'(\xi_k)}$$ 

If $f'$ is not changing rapidly between $x_k$ and $r$ then

$$f'(\xi_k) \approx \frac{f(x_{k+1}) - f(x_k)}{x_{k+1} - x_k}$$ 

We get the *error estimate*

$$
r - x_k \approx - \frac{f(x_k)}{\frac{f(x_{k+1}) - f(x_k)}{x_{k+1} - x_k}} = x_{k+1}
- x_k
$$ 

and the *relative error estimate*

$$
\frac{r - x_k}{r} \approx \frac{x_{k+1} - x_k}{x_{k+1}}
$$ 

Hence we can use the test

$$
|x_{k+1} - x_k| < \epsilon |x_{k+1}|, \qquad 0 < \epsilon \ll 1
$$ 

to decide on convergence. This is used in the code examples to stop the
iterations.

## Accuracy of FD approximation

:::{prf:example}

Compute derivative of

$$
f(x) = \sin(x)
$$

at $x=2\pi$ using the finite difference formula

$$
\frac{ f(x+h) - f(x) }{h}
$$

for $h = 10^{-1}, 10^{-2}, \ldots, 10^{-14}$.

```{code-cell}
f = lambda x: sin(x)

h = 10.0**arange(-1,-15,-1)
df= zeros(len(h))
x = 2.0*pi
f0= f(x)
for i in range(len(h)):
    f1 = f(x+h[i])
    df[i] = (f1 - f0)/h[i]
loglog(h,abs(df-1.0),'o-')
grid(True), xlabel('h'), ylabel('Error in derivative');
```

Initially, as $h$ decreases, the error in the FD approximation decreases, but for after some value, the error starts to increase with further decrease in $h$. Clearly, there seems to be an optimal value of $h$ for which the error is least.
:::


The secant method can be written as

$$
x_{k+1} = x_k - \frac{f(x_k)}{g_k}, \qquad g_k = \frac{f(x_k) - f(x_{k-1})}{x_k - x_{k-1}}
$$ 

As $x_k \to r$, the quantities $x_k - x_{k-1}$ and $f(x_k) - f(x_{k-1})$ can suffer from round-off error and $g_k$ can then be inaccurate. Let us analyze this using Taylor expansion for the finite difference approximation

$$
\frac{f(x+h) - f(x)}{h} = f'(x) + \half h f''(\xi), \qquad \textrm{$\xi$ between $x$ and $x+h$}
$$ 

We expect that as $h \to 0$, we get a more accurate approximation of the derivative. However, this assumes that we can evaluate $f(x)$ exactly and the arithmetic is exact, both of which are not true in a computer. What we compute on the computer is actually (there is also roundoff error in subtraction which is of similar order as unit round)

$$
\frac{\hat{f}(x+h) - \hat{f}(x)}{h}
$$ 

where, assuming $f(x) = O(1)$

$$
\hat{f}(x+h) = f(x+h) + \epsilon_2, \quad \hat{f}(x) = f(x) + \epsilon_1
$$

and

$$
\epsilon_1, \epsilon_2 \approx C \uround, \quad \textrm{the unit round of the computer}
$$

Hence 

$$
\begin{aligned}
\frac{\hat{f}(x+h) - \hat{f}(x)}{h} &=& \frac{f(x+h) - f(x)}{h} + \frac{C\uround}{h} \\
&=& f'(x) + \half h f''(\xi) + \frac{C\uround}{h}
\end{aligned}
$$ 

The error consists of 

1. discretization error: $\half h f''(\xi)$
1. round-off error: $\frac{C\uround}{h}$

As $h \to 0$, the discretization error decreases but the round-off error increases. There is an optimum value of $h$ for which the total error is minimized, given by

$$
h_{opt} = \left( \frac{2C \uround}{f''} \right)^\half
$$ 

For the above finite difference approximation, the optimum step size $h = O(\sqrt{\uround}) = \order{10^{-7}}$. Based on this analysis, we can modify the secant method as follows. If

$$|x_k - x_{k-1}| > \left| \frac{x_{k-1} + x_k}{2} \right| \sqrt{\uround}
$$

use the $g_k$ as given above. Otherwise, use the finite difference approximation

$$
g_k = \frac{f(x_k + h) - f(x_k)}{h}, \qquad h = |x_k| \sqrt{\uround}
$$

## Complex step method

We can take a complex step instead of a real step. Taylor expansion with a small complex step yields

$$
f(x+\ii h) = f(x) + \ii h f'(x) - \half h^2 f''(x) + O(\ii h^3)
$$

Hence 

$$
\imag f(x+\ii h) = h f'(x) + O(h^3)
$$ 

and

$$
\frac{\imag f(x+\ii h)}{h} = f'(x) + O(h^2)
$$ 

is a second order accurate approximation to the derivative. This formula does not involve subtraction of near equal quantities and hence it does not suffer from round-off error even if we take $h$ to be very small. In fact it is perfectly fine to take $h$ smaller than machine precision, e.g., $h = 10^{-20}$.

:::{prf:example}
We repeat the previous example using complex step method.

```{code-cell}
f = lambda x: sin(x)
h = 10.0**arange(-1,-15,-1)
df= zeros_like(h)
x = 2.0*pi
for i,dx in enumerate(h):
    df[i] = imag(f(x + 1j*dx)) / dx
loglog(h,abs(df-1.0),'o-')
grid(True), xlabel('h'), ylabel('Error in derivative')
print(df-1.0)
```

Once $h$ is below $10^{-7}$ the error in the derivative approximation is zero.
:::

[^1]: There are methods to differentiate a computer program called
    Automatic Differentiation.
