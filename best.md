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
  title: true
  headings: true
---

# Best approximation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

Let $X = (X, \norm{\cdot})$ be a normed spaced and suppose a given $x \in X$ is to be approximated by a $y \in Y$, where $Y$ is a fixed subset of $X$. The distance between $x$ and $Y$ is

$$
\delta = \delta(x,Y) = \inf_{y \in Y} \norm{x - y}
$$

If there is $y_0 \in Y$ such that

$$
\norm{x -y_0} = \delta
$$

then $y_0$ is called the **best approximation** to $x$ out of $Y$, since

$$
\norm{x - y_0} \le \norm{x - y}, \qquad \forall y \in Y
$$

## Existence

:::{prf:theorem} Existence
:label: thm:bestexist
If $Y$ is a finite dimensional subspace of a normed space $X=(X, \norm{\cdot})$, then for each $x \in X$ there exists a best approximation to $x$ out of $Y$.
:::

Let us apply this result to some examples.

:::{prf:example} Minimax approximation
$$
X = \cts[a,b], \qquad \norm{x} = \max_{t\in[a,b]} |x(t)|, \qquad Y = \poly_n
$$

Then by [](#thm:bestexist), there is a $p_n \in Y$ such that

$$
\max_{t \in [a,b]} |x(t) - p_n(t)| \le \max_{t \in [a,b]} |x(t) - y(t)|, \qquad \forall y \in \poly_n
$$
:::

:::{prf:example} Approximation in $L^p[a,b]$
\begin{gather}
X = \cts[a,b], \qquad \norm{x} = \left( \int_a^b |x(t)|^p \ud x \right)^{1/p}, \qquad p \ge 1 \\
Y = \poly_n
\end{gather}
Then by [](#thm:bestexist), there is a best approximation in $Y$. We only need to assume that $x \in L^p[a,b]$.
:::

The finite dimensionality of $Y$ is essential for existence.

:::{prf:example} Polynomials
Let

$$
Y = \{ p :[0,\shalf] \to \re :  p \textrm{ is a polynomial} \}
$$

and dim$(Y) = \infty$. It is a subspace of $\cts[0,\shalf]$  and consider 

$$
x(t) = \frac{1}{1 - t}
$$

Then for every $\epsilon > 0$ there is an $N$ such that

$$
y_n(t) = 1 + t + \ldots + t^n \qquad\textrm{satisfies}\qquad \norm{x - y_n}_\infty < \epsilon, \quad n > N
$$

and hence $\delta(x,Y) = 0$. But since $x$ is not a polynomial, there is no $y_0\in Y$ satisfying $\delta(x,Y) = \norm{x - y_0}_\infty = 0$.
:::

We can state the existence result in another way.

:::{prf:theorem} Existence
:label: thm:bestexist2
Let $X$ be a vector space. Given any $x \in X$ and linearly independent elements $\{ x_1, x_2, \ldots, x_n \} \subset X$, the problem of finding

$$
\min_{\{a_i\}} \norm{x - (a_1 x_1 + \ldots + a_n x_n)}
$$

has a solution. We only need $\norm{\cdot}$ to be a seminorm on $X$ and a norm on $Y = \textrm{span}\{ x_1, x_2, \ldots, x_n \}$.
:::

Let us apply this result to some examples.

:::{prf:example} Discrete minimax
Let $X = \cts[a,b]$ and $t_0, t_1, \ldots, t_m$ be $m+1$ distinct points in $[a,b]$ with $m \ge n$. Then by [](#thm:bestexist2), for any $x \in X$ the problem of determining

$$
\min_{a_0,a_1,\ldots,a_n} \max_{0 \le i \le m}| x(t_i) - (a_0 + a_1 t_i + \ldots + a_n t_i^n)|
$$

has a solution. Here $Y = \poly_n$ and

$$
\norm{p} = \max_{0 \le i \le m} |p(t_i)|
$$

is a norm on $Y$.
:::

:::{prf:example} Discrete least squares
Let $X = \cts[a,b]$ and $t_0, t_1, \ldots, t_m$ be $m+1$ distinct points in $[a,b]$ with $m \ge n$. Then by [](#thm:bestexist2), for any $x \in X$ the problem of determining

$$
\min_{a_0,a_1,\ldots,a_n} \sum_{i=0}^m [ x(t_i) - (a_0 + a_1 t_i + \ldots + a_n t_i^n)]^2
$$

has a solution. Here $Y = \poly_n$ and

$$
\norm{p} = \sum_{i=0}^m [p(t_i)]^2
$$

is a norm on $Y$.
:::

## Convexity

Let $M$ be the set of all best approximations

$$
M = \{ y \in Y : \norm{x-y} \le \norm{x - z}, \ \forall z \in Y \}
$$

:::{prf:lemma} Convexity
1. In a normed space $(X, \norm{\cdot})$ the set $M$ of best approximations to a given point $x$ out of a subspace $Y$ of $X$ is convex.

2. The set $M$ contains either one element or an infinite number of elements.
:::

:::{prf:definition} Strict convexity
A **strictly convex norm** $\norm{\cdot}$ is a norm such that for all $x \ne y$ with $\norm{x} = \norm{y}  = 1$, we have

$$
\norm{x + y} < 2
$$

A normed space with such a norm is called a strictly convex normed space.
:::

:::{prf:example}
Show that $(\complex^n, \norm{\cdot})$ with

$$
\norm{x} = \left( \sum_{i=1}^n |x_i|^p \right)^{1/p}, \qquad p > 1
$$

is a strictly convex normed space.
:::

:::{prf:lemma}
1. Every Hilbert space is strictly convex.
1. The space $\cts[a,b]$ is not strictly convex.  
    (Hint: Use $x(t)=1$ and $y(t) = (t-a)/(b-a)$.)
1. The space $L^p[a,b]$, $1 < p < \infty$ is strictly convex. (skip proof)
:::


## Uniqueness

:::{prf:theorem} Uniqueness
In a strictly convex normed space $X$ there is at most one best approximation to an $x \in X$ out of a given subspace $Y$.
:::

We cannot use the results of this chapter to show uniqueness of minimax problem since maximum norm is not strictly convex. We have shown the uniqueness of minimax by different techniques in previous chapter, using equioscillation property.

:::{prf:corollary}
:label: cor:hfdbest
For every $x$ in a Hilbert space $X$ and every finite dimensional subspace $Y$, there is a unique best approximation to $x$ out of $Y$.
:::

## Hilbert spaces

We have quite general results in inner product spaces.

:::{prf:theorem}
(1) Let $X$ be a Hilbert space and $Y$ be a convex subset which is complete. Then for every given $x \in X$, there is a unique best approximation $y \in Y$.

(2) If $Y$ is a complete subspace of $X$, then the best approximation $y \in Y$ is such that $x - y$ is orthogonal to $Y$.
:::

:::{prf:proof}
See [@Kreyszig1978], Theorem 3.3-1 and 3.3-2.
:::


In the above result, $Y$ was not assumed to be finite dimensional.

:::{prf:theorem}
If $Y$ is a finite dimensional subspace of a normed linear space $X$ then $Y$ is complete.
:::

:::{prf:proof}
See [@Kreyszig1978], Theorem 2.4-2.
:::

From these two results, we can again conclude [](#cor:hfdbest).

### Solution in finite dimensional case

If $X$ is a Hilbert space, and $Y$ is a finite dimensional subspace, then we can construct an orthonormal basis for $Y$

$$
Y = \textrm{span}\{x_1, x_2, \ldots, x_n\}, \qquad \ip{x_i, x_j} = \delta_{ij}
$$

The best approximation

$$
y = a_1 x_1 + a_2 x_2 + \ldots + a_n x_n
$$

to $x \in X$ must satisfy $x - y \perp Y$, i.e.,

$$
\ip{x - y, x_j} = 0, \qquad 1 \le j \le n
$$

which yields the solution

$$
a_j = \ip{x, x_j}, \qquad 1 \le j \le n
$$

and

$$
y = \sum_{j=1}^n \ip{x, x_j} x_j
$$

## Uniform approximation: uniqueness

We now consider the best approximation problem in maximum norm. Since we dont have convecity, we need different techniques to prove uniqueness.

:::{prf:definition} Haar condition
A subspace $Y$ of $\cts[a,b]$ with dim$(y)=n$ is said to satisfy the Haar condition if every $y \in Y$, $y \ne 0$ has at most $n-1$ zeros in $[a,b]$.
:::

:::{prf:example} Polynomials
Let $Y = \poly_m$, then

$$
n = \textrm{dim}(Y) = m+1
$$

Every nonzero polynomial $y \in Y$ has atmost $m = n-1$ zeros; thus $\poly_m$  satisfies the Haar condition.
:::

:::{prf:lemma}
The Haar condition is equivalent to the condition that for every basis $\{ y_1, y_2, \ldots, y_n\}$ of $Y$ and every $n$-tuple of distinct points $t_1, \ldots, t_n$ in $[a,b]$,

$$
\det V \ne 0, \qquad V_{ij} = y_i(t_j)
$$
:::

Recall that $V$ is the Vandermonde matrix which arises in interpolation problem. Hence the Haar condition says that the interpolation problem in $Y$ using any $n$ distinct points has a unique solution.

:::{prf:definition}
An extremal point of an $x \in \cts[a,b]$ is a $t_0 \in [a,b]$ such that $|x(t_0)| = \norm{x}$.
:::

:::{prf:lemma} Extremal points
Suppose a subspace $Y$ of $\cts[a,b]$ with dim$(Y)=n$ satisfies the Haar condition. If for a given $x \in \cts[a,b]$ and a $y \in Y$ the function $x-y$ has less than $n+1$ extremal points, then $y$ is not a best approximation to $x$ out of $Y$.
:::

:::{prf:remark} Polynomials
Notice the similarity with the equioscillation property. If $Y = \poly_m$, then the lemma says that if $y \in \poly_m$ is such that the error $x - y$ has less than $m+2$ extremal points, then $y$ cannot be a best approximation.
:::

:::{prf:theorem} Haar uniqueness theorem
Let $Y$ be a finite dimensional subspace of $\cts[a,b]$. Then the best approximation out of $Y$ is unique for every $x \in \cts[a,b]$ if and only if $Y$ satisfies the Haar condition.
:::

:::{prf:corollary} Polynomials
Given any $x \in \cts[a,b]$, there is a unique best approximation $y \in \poly_m$ to $x$.
:::
