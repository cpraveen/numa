---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
math:
  '\re': '\mathbb{R}'
  '\df': '\frac{\partial #1}{\partial #2}'
  '\half': '\frac{1}{2}'
  '\ud': '\textrm{d}'
  '\dd': '\frac{\ud #1}{\ud #2}'
  '\sign': '\textrm{sign}'
  '\cts': 'C'
---

# Root finding

+++

Given a continuous function $f : [a,b] \to \re$, find an $\alpha \in [a,b]$ such that
$f(\alpha) = 0$. Such an $\alpha$ is called a root of $f$ or a zero of $f$. There could
be many roots for a given function. At a root, the positive and negative parts of $f$
cancel one another. So we have to careful about roundoff errors and care must be
taken to evaluate the function so as to avoid loss of significance errors. Since we work
with finite precision on the computer, we cannot expect to find an $\alpha$ that
makes the function exactly zero. At best, we hope to keep relative error under
control. In this sense, zeros far from the origin cannot be computed with great
absolute accuracy. Since the purpose is to solve some problem that involves
approximations, it is not necessary to exactly locate the root, and an approximation is
sufficient.

+++

````{prf:example} Minima of a function

To find the location where a given function
$f : \re^N \to \re$ attains a local minimum, we have to solve the equation $f'(x) =
0$. Note that $g = f' : \re^N \to \re^N$ and we have a multi-dimensional root finding
problem: find $x = (x_1, x_2, \ldots, x_N)$ such that

$$
g_i(x_1, x_2, \ldots, x_N) = \df{f}{x_i}(x_1, x_2, \ldots, x_N) = 0, \qquad 1 \le i \le N
$$

````

+++

````{prf:example} Solving ODE

Consider the system of ODE

$$
\dd{u(t)}{t} = f(u(t)), \qquad t \in [0,T]
$$

where $u \in \re^N$ and $f : \re^N \to \re^N$. We are given some initial condition $u(0) = u_0 \in \re^N$. The backward Euler scheme is given by

$$
\frac{u^{n+1} - u^n}{\Delta t} = f(u^{n+1})
$$

This is a set of $N$ non-linear equations for the unknown solution $u^{n+1}$. We can define

$$
H(u) = \frac{u - u^n}{\Delta t} - f(u)
$$

The solution $u^{n+1}$ we seek is the root of the function $H : \re^N \to \re^N$. We have to solve a root problem for every time interval.

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

+++

## Bisection method

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
:width: 800%
:align: center
```

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/latex/p2.svg
:width: 800%
:align: center
```
