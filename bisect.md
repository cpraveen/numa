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

# Bisection method

```{include} math.md
```

+++

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

Let $f$ be a continuous function. Given an interval $[a,b]$ such that

$$
\sign f(a) \ne \sign f(b), \qquad f(a) \ne 0 \ne f(b)
$$

then $[a,b]$ is called a {\em non-trivial bracket} for a root of $f$. The basic idea of bisection method is: starting from a non-trivial bracket, find a new bracketing interval which is smaller in size. As the name suggests, we divide the interval into two equal and smaller intervals. First divide $[a,b]$ into $[a,c]$ and $[c,b]$ where $c = \half(a+b)$. Then there are three possibilities.

1. $f(c) = 0$, then we have found the root exactly.
1. $f(c) \ne 0$ and $\sign f(c) \ne \sign f(b)$. Then $[c,b]$ is a non-trivial bracket.
1. $f(c) \ne 0$ and $\sign f(a) \ne \sign f(c)$. Then $[a,c]$ is a non-trivial bracket.

We repeat the above process until the root is found or the length of the bracketing interval is reduced to a small desired size. The method may not give a unique value for the root since it can lie anywhere in the bracketing interval. We can take the mid-point of the bracketing interval as the estimate of the root, in which case the error in the root is at most equal to half the length of the bracketing interval.

The length of the bracketing interval reduces by half in each iteration, which means that the algorithm has a monotonic convergence behaviour. If

$$
L_0 = b-a = \textrm{initial length}
$$

then after $k$ iterations

$$
L_k = \frac{1}{2^k} L_0
$$

As a convergence criterion, we can put a tolerance on the length: stop if $L_k \le \epsilon$ for some $\epsilon > 0$ small. To achieve this tolerance, we require

$$
k \approx \log_2 \frac{L_0}{\epsilon}
$$

If $L_0 = 1$ and $\epsilon = 10^{-6}$ then we need about 20 iterations. If we demand more accuracy, then the number of iterations will increase, though only logarithmically. The bisection method is shown in Algorithm~(1). We have put an upper limit on the number of iterations and we also specify the stopping criterion in terms of the length of the bracketing interval and/or the function value.

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/latex/p1.svg
:width: 100%
:align: center
```

The way we have written the algorithm, we store the entire history of           computations which is not necessary. In actual practice, we only need to keep some of the latest results in memory and this is illustrated in Algorithm~(2).

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/latex/p2.svg
:width: 100%
:align: center
```

We are checking the tolerance in terms of absolute error in the length of       bracketing interval without regard to the size of the root. If $\epsilon = 10^{-6}$ and if the root is near one, then bisection method will return roughly six accurate digits. If  $\epsilon = 10^{-6}$ and the root is $\approx 10^{-7}$ then we cannot expect no figures of accuracy. On the other hand, if the root is $\gg 1$ and we use $\epsilon \ll 1$, then we may never be able to achieve this tolerance or it may require far too many iterations. It is better to use a relative error tolerance

$$
c = \half(a+b), \qquad \textrm{If } |b - a| < \epsilon |c| \quad \textrm{return } c
$$

+++

:::{prf:example}
First we define the function for which the root is required.

```{code-cell}
def fun(x):
    f = exp(x) - sin(x)
    return f
```

Let us plot and visualize the function.

```{code-cell}
x = linspace(-4,-2,100)
f = fun(x)
plot(x,f,'r-')
grid(True); xlabel('x'); ylabel('f');
```

From the figure, we see that [−4,−2] is a bracketing interval. We now implement the bisection method.

```{code-cell}
# Initial interval [a,b]
a, b = -4.0, -2.0

N     = 100    # Maximum number of iterations
eps   = 1.0e-4 # Tolerance on the interval
delta = 1.0e-4 # Tolerance on the function

fa, fb = fun(a), fun(b)
sa, sb = sign(fa), sign(fb)

for i in range(N):
    e = b-a
    c = a + 0.5*e
    if abs(e) < eps*abs(c):
        print("Interval size is below tolerance\n")
        break
    fc = fun(c)
    if abs(fc) < delta:
        print("Function value is below tolerance\n")
        break
    sc = sign(fc)
    if sa != sc:
        b  = c
        fb = fc
        sb = sc
    else:
        a  = c
        fa = fc
        sa = sc
    print("%5d %16.8e %16.8e %16.8e"%(i+1,c,abs(b-a),abs(fc)))
```

:::

+++

:::{prf:example}
Let us define three different functions for which we will find the root.

```{code-cell}
def f1(x):
    f = exp(x) - sin(x)
    return f

def f2(x):
    f = x**2 - 4.0*x*sin(x) + (2.0*sin(x))**2
    return f

def f3(x):
    f = x**2 - 4.0*x*sin(x) + (2.0*sin(x))**2 - 0.5
    return f
```

The following function implements the bisection method. Note that it takes some default arguments.

```{code-cell}
def bisect(fun,a,b,N=100,eps=1.0e-4,delta=1.0e-4,debug=False):
    fa, fb = fun(a), fun(b)
    sa, sb = sign(fa), sign(fb)

    if abs(fa) < delta:
        return (a,0)

    if abs(fb) < delta:
        return (b,0)

    # check if interval is admissible
    if fa*fb > 0.0:
        if debug:
            print("Interval is not admissible\n")
        return (0,1)

    for i in range(N):
        e = b-a
        c = a + 0.5*e
        if abs(e) < eps*abs(c):
            if debug:
                print("Interval size is below tolerance\n")
            return (c,0)
        fc = fun(c)
        if abs(fc) < delta:
            if debug:
                print("Function value is below tolerance\n")
            return (c,0)
        sc = sign(fc)
        if sa != sc:
            b = c
            fb= fc
            sb= sc
        else:
            a = c
            fa= fc
            sa= sc
        if debug:
            print("%5d %16.8e %16.8e %16.8e"%(i+1,c,abs(b-a),abs(fc)))

    # If we reached here, then there is no convergence
    print("No convergence in %d iterations !!!" % N)
    return (0,2)
```

**First function**

```{code-cell}
M     = 100        # Maximum number of iterations
eps   = 1.0e-4     # Tolerance on the interval
delta = 1.0e-4     # Tolerance on the function
a, b  = -4.0, -2.0 # Initial interval
r,status = bisect(f1,a,b,M,eps,delta,True)
```

**Second function**

```{code-cell}
M     = 100        # Maximum number of iterations
eps   = 1.0e-4     # Tolerance on the interval
delta = 1.0e-4     # Tolerance on the function
a, b  = -4.0, -2.0 # Initial interval
r,status = bisect(f2,a,b,M,eps,delta,True)
```

Lets visualize this function.

```{code-cell}
x = linspace(-3,3,500)
y = f2(x)
plot(x,y)
grid(True)
```

It looks like there are double roots; bisection method cannot compute such roots since the function value does not change around the root !!!

**Third function**

```{code-cell}
M     = 100        # Maximum number of iterations
eps   = 1.0e-4     # Tolerance on the interval
delta = 1.0e-4     # Tolerance on the function
a, b  = -3.0, +2.0 # Initial interval
r,status = bisect(f3,a,b,M,eps,delta,True)
```

We dont need to specify all arguments to the function as some of them have default values.

```{code-cell}
r,status = bisect(f3,a,b,debug=True)
```

:::

+++

:::{exercise}
We have used a `for` loop to perform the iterations. Modify the code to use `while` loop instead.
:::

+++

## Convergence rate

The bisection method always converges and the root $\alpha$ is reached in the   limit

$$
\alpha = \lim_{k \to \infty} c_k
$$

After $k$ iterations, if we take the mid-point of the bracketing interval as    the estimate of the root, the error is bounded by

$$
|\alpha - c_k| \le \frac{L_k}{2} = \frac{b-a}{2^{k+1}}
$$

:::{prf:definition}

(1) A sequence of iterates $\{x_k : k \ge 0\}$ is said to converge with **order** $p \ge 1$ to a point $\alpha$ if

$$
|\alpha - x_{k+1}| \le c |\alpha - x_k|^p, \qquad k \ge 0
$$

for some $c > 0$.

(2) If $p=1$, the sequence is said to {\em converge linearly} to $\alpha$. In   this case, we require $c < 1$; the constant $c$ is called the **rate of linear convergence** of $x_k$ to $\alpha$.
:::

:::{prf:remark}
If a sequence converges linearly, then

$$
|\alpha - x_{k+1}| \le c |\alpha - x_k|, \qquad 0 < c < 1
$$

For the bisection, we do not have a relation between $x_k$, $x_{k+1}$ so we     cannot establish an inequality as above. Instead we can take $L_k$ as the measure of   the error and

$$
L_{k+1} \le \half L_k
$$

We see that bisection method has linear convergence with rate $c = \half$.
:::

:::{prf:remark}
The bisection method

* is guaranteed to converge to some root
* provides reasonable error bound
* has reasonable stopping criteria

However, the convergence is slow, especially in the final iterations. Also, it  does not make use of the structure of the function $f$. E.g., if $f$ is linear, then we       should be able to quickly find its roots, say in just one iteration, but the bisection method     will be as slow as for some general non-linear function. It is still a useful method to    first obtain a small bracketing interval and then we can switch to more sophisticated and    fast
algorithms like Newton-Raphson method.
:::
