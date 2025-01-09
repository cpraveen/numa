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

# Root finding

```{include} math.md
```

+++

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

Given a continuous function $f : [a,b] \to \re$, find an $\alpha \in [a,b]$ such that $f(\alpha) = 0$. Such an $\alpha$ is called a root of $f$ or a zero of $f$. There could be many roots for a given function. At a root, the positive and negative parts of $f$ cancel one another. So we have to careful about roundoff errors and care must be taken to evaluate the function so as to avoid loss of significance errors. Since we work with finite precision on the computer, we cannot expect to find an $\alpha$ that makes the function exactly zero. At best, we hope to keep relative error under control. In this sense, zeros far from the origin cannot be computed with great absolute accuracy. Since the purpose is to solve some problem that involves approximations, it is not necessary to exactly locate the root, and an approximation is sufficient.

+++

````{prf:example} Minima of a function

To find the location where a given function $f : \re^N \to \re$ attains a local minimum/maximum, we have to solve the equation $f'(x) = 0$. Note that $g = f' : \re^N \to \re^N$ and we have to solve: find $x = (x_1, x_2, \ldots, x_N)$ such that

$$
g_i(x_1, x_2, \ldots, x_N) = \df{f}{x_i}(x_1, x_2, \ldots, x_N) = 0, \qquad 1 \le i \le N
$$

We have $N$ equations for the $N$ unknowns, a multi-dimensional root finding problem.

````

+++

````{prf:example} Solving ODE

Consider the system of ODE

$$
\dd{u(t)}{t} = f(u(t)), \qquad t \in [0,T]
$$

where $u \in \re^N$ and $f : \re^N \to \re^N$. We are given some initial condition $u(0) = u_0 \in \re^N$. The backward Euler scheme is given by: $u^0 = u_0$ and

$$
\frac{u^{n+1} - u^n}{\Delta t} = f(u^{n+1}), \qquad n=0,1,\ldots
$$

This is a set of $N$ non-linear equations for the unknown solution $u^{n+1}$. We can define

$$
H(u) = \frac{u - u^n}{\Delta t} - f(u)
$$

The solution $u^{n+1}$ we seek is the root of the function $H : \re^N \to \re^N$. We have to solve a multi-dimensional root problem for every time interval.

````

+++

## Bracketing the root

There are no general methods to find the roots of an arbitrary function. The user must have some knowledge of where the root might possibly lie. It is hence useful to first find some small interval in which we expect the root to lie.

```{prf:theorem} Intermediate Value Theorem
If $f \in \cts[a,b]$ and $y$ lies between $f(a)$ and $f(b)$, then there is atleast one point $x \in [a,b]$ such that $y = f(x)$.
```

Around a root, the function takes both positive and negative values. The IVT can be used to find an interval containing the root.

```{prf:lemma} Bracketing interval
If $\sign f(a) \ne \sign f(b)$, then there is alteast one root of $f$ in the interval $[a,b]$.
```
