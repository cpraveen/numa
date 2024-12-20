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

# Gauss quadrature

```{include} math.md
```

In the quadrature methods studied so far, the nodes were first chosen
and then we determine the weights via some function interpolation
process. The resulting methods have certain degree of precision in terms
of what polynomials they can exactly integrate.

Let us consider the task of computing
$$I(f) = \int_a^b w(x) f(x) \ud x$$ where $w(x)$ is some fixed positive
weight function by a formula of the type
$$I_n(f) = \sum_{j=1}^n w_j f(x_j)$$ Here $w_j$ are weights and should
not be confused with the weight function $w(x)$, though they will be
related to one another in some way.

```{note} Question
Fix some $n \ge 1$. Find the nodes $\{ x_j \}$ and weights $\{ w_j \}$
so that $I_n(f)$ has maximum degree of precision, i.e., maximize $m$
such that $$I_n(p) = I(p), \qquad \forall p \in \poly_m$$ The evaluation
of $f$ might be expensive, hence we would like to find quadrature rules
that use few number of nodes but have high degree of precision.
```

## Unit weight function

Let us take the interval to be $[-1,+1]$ and $w(x)=1$. Then we want to
compute
$$\int_{-1}^1 f(x) \ud x \approx I_n(f) = \sum_{j=1}^n w_j f(x_j)$$ The
error $$E_n(f) = \int_{-1}^1 f(x) - \sum_{j=1}^n w_j f(x_j)$$ Note thet
$E_n : \cts[-1,1] \to \re$ is a linear operator and in particular
$$E_n(a_0 + a_1 x + \ldots a_m x^m) = a_0 E_n(1) + a_1 E_n(x) + \ldots + a_m E_n(x^m)$$
Hence $$E_n(p) = 0, \qquad \forall p \in \poly_m$$ if and only if
$$E_n(x^i) = 0, \qquad 0 \le i \le m$$

### Case $n=1$:

The quadrature formula is $$I_1(f) = w_1 f(x_1)$$ and we have to find
$w_1$ and $x_1$.
$$E_1(1) = 0 \quad\implies\quad \int_{-1}^1 \ud x - w_1 = 0$$
$$E_1(x) = 0 \quad\implies\quad \int_{-1}^1 x \ud x - w_1 x_1 = 0$$
which yields $$w_1 = 2, \qquad x_1 = 0$$ so that
$$\int_{-1}^1 f(x) \ud x \approx 2 f(0)$$ This is just the
midpoint-rule. It uses only one point but is exact for linear
polynomials. Note that trapezoidal method uses two points but is also
exact only for linear polynomials. In terms of the number of nodes used,
the mid-point rule is more optimal since it needs less function
evaluation to give same degree of precision.

### Case $n=2$:

The quadrature formula is $$I_2(f) = w_1 f(x_1) + w_2 f(x_2)$$ and we
have to determine four quantities $w_1, w_2, x_1, x_2$. We can try to
make the formula exact for upto cubic polynomials
$$E_2(x^i) = \int_{-1}^1 x^i \ud x - [w_1 x_1^i + w_2 x_2^i] = 0, \qquad i=0,1,2,3$$
which gives us four equations $$\begin{aligned}
w_1 + w_2 &=& 2 \\
w_1 x_1 + w_2 x_2 &=& 0 \\
w_1 x_1^2 + w_2 x_2^2 &=& 2/3 \\
w_1 x_1^3 + w_2 x_2^3 &=& 0
\end{aligned}$$ The unique solution is
$$w_1 = w_2 = 1, \qquad x_1 = -\frac{1}{\sqrt{3}}, \qquad x_2 = \frac{1}{\sqrt{3}}$$
This rule uses 2 points but has degree of precision 3. Simpson rule also
has degree of precision 3 but uses 3 points.

### Case of general $n$:

There are $2n$ free parameters $\{w_j, x_j\}$, $j=1,2,\ldots,n$ and we
make the quadrature exact for polynomials of degree upto $2n-1$
$$E_n(x^i) = 0, \qquad i=0,1,2,\ldots,2n-1$$ which gives $2n$ equations
$$\sum_{j=1}^n w_j x_j^i = \begin{cases}
0 & i=1,3,\ldots,2n-1 \\
\frac{2}{i+1} & i=0,2,\ldots,2n-2
\end{cases}$$ We have system of non-linear equations and it is not
possible to solve them analytically. We will take a digression to study
some properties of orthogonal polynomials which will help us to derive
these Gauss quadrature rules.

:::{exercise}
Find two point quadrature rule of maximal degree of precision to compute
$$I(f) = \int_0^1 \sqrt{x} f(x) \ud x$$
:::

## Solution for general case

We can find the quadrature rule without having to solve a matrix problem.  Let $f : [-1,1] \to \re$ and for some integer $n \ge 1$, let $\{ x_j, 1 \le j \le n \}$ be the roots of Legendre polynomial $P_n$. Consider the interpolation

$$
f_n(x) = \sum_{j=1}^n f(x_j) \ell_j(x)
$$

where $\ell_j$ are the degree $n-1$ Lagrange polynomials supported at the roots. Since $f(x_j) = f_n(x_j)$ for $1 \le j \le n$, we can write

$$
f(x) = f_n(x) + P_n(x) R(x)
$$

for some function $R$. Now

$$
I(f) &= \int_{-1}^1 f(x) \ud x \\
&= \int_{-1}^1 f_n(x) \ud x + \int_{-1}^1 P_n(x) R(x) \ud x \\
&= \sum_{j=1}^n w_j f(x_j) + \int_{-1}^1 P_n(x) R(x) \ud x \\
&= I_n(f) + \int_{-1}^1 P_n(x) R(x) \ud x \\
$$

where

$$
w_j = \int_{-1}^1 \ell_j(x) \ud x
$$

If $f \in \poly_{2n-1}$ then $R \in \poly_{n-1}$ and

$$
\int_{-1}^1 P_n(x) R(x) \ud x = 0, \qquad \forall R \in \poly_{n-1}
$$

Thus the quadrature formula is exact for any $f \in \poly_{2n-1}$. Later, we will show that this is optimal and derive simpler formulae for the weights.
