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

# Least squares approximation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

## Best approximation in 2-norm

The minimax approximation is the best approximation since it minimizes the error at all points. This involves an optimization problem where we have to minimize the maximum norm of the error and the Remez algorithm is expensive. Optimization problems can be solved by gradient based methods but we cannot employ this due to the maximum norm which is not differentiable. To use gradient methods, let us change the norm to the $L^2$ norm

$$
\norm{g}_2 = \left( \int_a^b |g(x)|^2 \ud x \right)^2, \qquad g \in \cts[a,b]
$$

For a given $f \in \cts[a,b]$ and $n \ge 0$, define

$$
M_n(f) = \inf_{r \in \poly_n} \norm{f - r}_2
$$ 

1. Is there a best approximating polynomial $r_n^* \in \poly_n$, i.e.,
    $$
    \norm{f - r_n^*}_2 = M_n(f)
    $$ 

1. Is it unique ? 
1. How can we compute it ?

+++

:::{prf:example}
Consider the function 

$$
f(x) = \ee^x, \qquad x \in [-1,1]
$$ 

Let us find

$$
r_1(x) = b_0 + b_1 x
$$ 

by minimizing the $L^2$ error norm

$$
\norm{f - r_1}_2^2 = \int_{-1}^1[\ee^x - b_0 - b_1 x]^2 \ud x =: F(b_0, b_1)
$$

The optimality conditions are 

$$
\df{F}{b_0} = \df{F}{b_1} = 0
$$ 

yield the linear system of equations

\begin{align}
b_0 \int \ud x + b_1 \int x \ud x &= \int \ee^x \ud x \\
b_0 \int x \ud x + b_1 \int x^2 \ud x &= \int x \ee^x \ud x
\end{align}

whose solution is

$$
\begin{aligned}
b_0 &= \half \int_{-1}^1 \ee^x \ud x = \sinh(1) \approx 1.1752 \\
b_1 &= \frac{3}{2} \int_{-1}^1 x \ee^x \ud x = \frac{3}{\ee} \approx 1.1036
\end{aligned}
$$ 

so that

$$
r_1^*(x) = 1.1752 + 1.1036 x, \qquad \norm{f - r_1^*}_\infty \approx 0.44
$$

Compare this with

\begin{align}
\textrm{minimax:} & \qquad \norm{f - p_1^*}_\infty \approx 0.2788 \\
\textrm{Taylor:} & \qquad \norm{f - p_1}_\infty \approx 0.7182
\end{align}

We have the ordering of error in max norm

$$
\textrm{Error of minimax $\lt$ Error of least squares $\lt$ Error of Taylor}
$$
:::

+++

:::{prf:example} Cubic least squares
It is given by

$$
r_3^*(x) = 0.996294 + 0.997955 x + 0.536722 x^2 + 0.176139 x^3
$$ 

```{code-cell}
:tags: remove-input
r3 = lambda x: 0.996294 + 0.997955*x + 0.536722*x**2 + 0.176139*x**3
x  = linspace(-1,1,100)
plot(x,exp(x)-r3(x))
title('Error of cubic least squares approximation of exp(x)')
grid(True), xlabel('x');
```

and the error oscillates in sign but the error is somewhat larger at the end points, with $$\norm{f - r_3^*}_\infty \approx 0.0112$$ This error is similar to error in cubic interpolating polynomial.
:::

## General least squares problem

Let there be given a weight function 

$$
w(x) \ge 0, \qquad x \in [a,b]
$$

such that

1.  $\int_a^b|x|^n w(x) \ud x < \infty$, $\forall n \ge 0$

2.  If $\int_a^b w(x) g(x) \ud x = 0$ for some non-negative continuous
    function $g$ $\implies$ $g(x) \equiv 0$ on $(a,b)$.

Examples of commonly used weight functions are 

$$
\begin{aligned}
w(x) &= 1, \qquad\qquad a \le x \le b \\
w(x) &= \frac{1}{\sqrt{1-x^2}}, \quad -1 \le x \le +1 \\
w(x) &= \ee^{-x}, \qquad\quad 0 \le x \le \infty \\
w(x) &= \ee^{-x^2}, \quad\quad -\infty < x < \infty
\end{aligned}
$$

### Least squares problem

Given $f \in \cts[a,b]$, find the polynomial $r_n^* \in \poly_n$ which
minimizes 

$$
\int_a^b w(x)[f(x) - r(x)]^2 \ud x
$$ 

among all polynomials $r \in \poly_n$. 

1. Does it exist ? 
1. Is it unique ? 
1. How to find it ?

## Least squares using monomials

Find the coefficients in the polynomial

$$
r_n(x) = \sum_{j=0}^n a_j x^j, \qquad x \in [a,b] = [0,1]
$$ 

so that, with weight function $w(x) \equiv 1$

$$
F(a_0,a_1,\ldots,a_n) := \int_0^1 [f(x) - r_n(x)]^2 \ud x
$$ 

is minimized. The necessary optimality conditions are

$$
\df{F}{a_i} = 0 \quad\implies\quad \sum_{j=0}^n a_j \int_0^1 x^{i+j} \ud x = \int_0^1 f(x) x^i \ud x, \qquad i=0,1,\ldots,n
$$

This is a set of $n+1$ coupled linear equations

$$
Aa = b, \qquad A_{ij} = \frac{1}{i+j+1}, \qquad 0 \le i,j \le n
$$ 

$A$ is called the Hilbert matrix and it is highly ill conditioned

$$
\cond(A) = \order{ \frac{(1 + \sqrt{2})^{4n+4}}{\sqrt{n+1}} }
$$ 

The monomial basis $1,x,x^2,\ldots,x^n$ is linearly independent but for large powers $m,n$, the functions $x^m, x^n$ are nearly identical which causes the bad condition numbers. Note that we can write the matrix elements are

$$
A_{ij} = \int_0^1 \phi_i(x) \phi_j(x) \ud x, \qquad \phi_i(x) = x^i
$$

If we could make the matrix diagonal, i.e.

$$
A_{ij} = 0, \qquad i \ne j
$$ 

then its solution would become trivial and we do not to worry about ill-conditioning. We should use a set of basis functions which has the above property, instead of using the monomials.

## Orthogonal polynomials

Let $w : (a,b) \to (0,\infty)$ be a weight function and define the inner
product 

$$
\ip{f,g} = \int_a^b w(x) f(x) g(x) \ud x
$$ 

The inner product satisfies these properties.

1.  $\ip{\alpha f, g} = \ip{f,\alpha g} = \alpha \ip{f,g}$ for all
    $\alpha \in \re$

2.  $\ip{f_1 + f_2, g} = \ip{f_1,g} + \ip{f_2,g}$ and
    $\ip{f,g_1+g_2} = \ip{f,g_1} + \ip{f,g_2}$

3.  $\ip{f,g} = \ip{g,f}$

4.  $\ip{f,f} \ge 0$ for all $f \in \cts[a,b]$ and $\ip{f,f} = 0$ iff
    $f=0$

Such an inner product gives rise to a norm

$$
\norm{f} = \sqrt{ \ip{f,f} } = \sqrt{ \int_a^b w(x) |f(x)|^2 \ud x }
$$

Morever, we have the Cauchy-Schwarz inequality

$$
|\ip{f,g}| \le \norm{f} \cdot \norm{g}
$$

:::{prf:definition} Orthogonal functions
Two functions $f,g$ are said to be orthogonal with respect to the inner
product $\ip{\cdot,\cdot}$ if 

$$
\ip{f,g} = 0
$$
:::

+++

:::{prf:theorem} Orthogonal polynomials
There exists a sequence of polynomials $\{ \phi_n : n \ge 0\}$ with degree$(\phi_n) =n$ and 

$$
\ip{\phi_n, \phi_m} = 0 \qquad n \ne m
$$ 

In addition, we can construct the sequence with the additional properties

1.  $\ip{\phi_n,\phi_n} = 1$ for all $n$.

2.  The coefficient of $x^n$ in $\phi_n$ is positive.

then the sequence $\{ \phi_n \}$ is unique.
:::

:::{prf:proof}
The proof is constructive and uses the Gram-Schmidt orthogonalization process. 

(1) Start with the constant polynomial $\phi_0(x) = c =$ constant and normalize it

$$
1 = \ip{\phi_0,\phi_0} = c^2 \int_a^b w(x) \ud x
$$ 

which fixes the constant 

$$
c = \left( \int_a^b w(x) \ud x \right)^{-1}
$$ 

(2) To construct $\phi_1(x)$, start with 

$$
\psi_1(x) := x + a_{10} \phi_0(x)
$$ 

Make it orthogonal to $\phi_0$

$$
\ip{\psi_1, \phi_0} = 0 \quad\implies\quad \ip{x,\phi_0} + a_{10} \underbrace{\ip{\phi_0,\phi_0}}_{=1} = 0 \quad\implies\quad a_{10} = -\ip{x,\phi_0}
$$

Now normalize by defining the next function

$$
\phi_1(x) := \frac{1}{{\norm{\psi_1}}} \psi_1(x) \quad\implies\quad \ip{\phi_1,\phi_1} = 1
$$

(3) Suppose $\phi_0, \phi_1, \ldots, \phi_{n-1}$ have been found. For the next function which is of degree $n$, define

$$
\psi_n(x) = x^n + a_{n,n-1} \phi_{n-1}(x) + \ldots + a_{n,0} \phi_0(x)
$$

Determine the coefficients by making $\psi_n$ orthogonal to all the $\phi_0, \phi_1, \ldots, \phi_{n-1}$. For $j=0,1,\ldots,n-1$

\begin{align}
\ip{\psi_n,\phi_j} = 0 \quad\implies & \quad \ip{x^n,\phi_j} + a_{n,j} \underbrace{\ip{\phi_j,\phi_j}}_{=1} = 0 \\
\quad\implies & \quad a_{n,j} = -\ip{x^n, \phi_j} \\
\phi_n(x) = \frac{1}{\norm{\psi_n}} \psi_n(x) \quad\implies & \quad \ip{\phi_n,\phi_n} = 1
\end{align}
:::

+++

:::{prf:lemma} Orthonormal basis for polynomials
(1) The orthonoral set of polynomials

$$
V_n = \{\phi_0, \phi_1, \ldots, \phi_n\}
$$ 

is linearly independent and forms a basis for $\poly_n$. 

(2) Any polynomial $p \in \poly_n$ can be written as

$$
p = \sum_{j=0}^n a_j \phi_j, \qquad a_j = \ip{p,\phi_j}
$$
:::

:::{prf:proof}
(1) Let 

$$
\sum_{j=0}^n c_j \phi_j(x) = 0, \quad x \in [a,b]
$$

Taking inner product with $\phi_i$, $0 \le i \le n$ and using the orthogonality 

$$
c_i \ip{\phi_i, \phi_i} = 0 \quad\implies\quad c_i = 0
$$

This proves the linear independence of $V_n$.

(2) In the proof of Gram-Schmidt theorem, we had the intermediate result

$$
x^j = \norm{\psi_j} \phi_j(x) - (a_{j,j-1} \phi_{j-1}(x) + \ldots + a_{j,0} \phi_0(x)), \quad 1 \le j \le n
$$

Hence any $p \in \poly_n$ can be written as a linear combination of $V_n$.

**Alternate proof.** 
$\poly_n$ is of dimension $n+1$ and the set $V_n$ has $n+1$ linearly independent elements of $\poly_n$. Hence they form a basis for $\poly_n$.

(3) Any $p \in \poly_n$ can be written as

$$
p = \sum_{j=0}^n a_j \phi_j, \qquad a_j = \ip{p,\phi_j}
$$ 

and taking inner product with $\phi_i$ yields $a_i = \ip{p,\phi_i}$.
:::

+++

:::{prf:remark}
If the $\phi_i$ are only orthogonal, then any $p \in \poly_n$ can be
written as

$$
p = \sum_{j=0}^n a_j \phi_j, \qquad a_j = \frac{\ip{p,\phi_j}}{\ip{\phi_j,\phi_j}}
$$
:::

:::{prf:corollary}
Let $\{\phi_j, j=0,1,\ldots\}$ be an orthogonal family of polynomials with degree $\phi_j=j$. For any $p \in \poly_n$ 

$$
\ip{p,\phi_m} = 0, \qquad m=n+1,n+2,\ldots
$$
:::

:::{prf:proof}
Any $p\in\poly_n$ can be written as

$$
p = \sum_{j=0}^n a_j \phi_j
$$ 

Hence for $m > n$

$$
\ip{p,\phi_m} = \sum_{j=0}^n a_j \ip{\phi_j,\phi_m} = 0
$$
:::

## Legendre polynomials

Take $[a,b] = [-1,1]$ and weight function 

$$
w(x) \equiv 1
$$ 

The resulting orthogonal sequence are called the Legendre polynomials.

$$
\phi_0(x) = \frac{1}{\sqrt{2}}, \quad \phi_1(x) = \sqrt{\frac{3}{2}} x, \quad \phi_2(x) = \half \sqrt{\frac{5}{2}} (3x^2 - 1)
$$

They are usually defined in terms of the Rodrigues Formula

$$
P_0(x)=1, \qquad P_n(x) = \frac{(-1)^n}{2^n n!} \dd{^n}{x^n}(1-x^2)^n, \qquad n \ge 1
$$

which are polynomial solutions of Legendre's differential equation

$$
\dd{}{x}\left[(1-x^2) \dd{P_n}{x}\right] + n(n+1)P_n = 0, \qquad P_n(1) = 1
$$

The $P_n$ satisfy the orthogonality conditions

$$
\ip{P_n,P_n} = \frac{2}{2n+1}, \qquad \ip{P_n,P_m} = 0, \quad n \ne m
$$

and the recurrence relation

$$
P_{n+1}(x) = \frac{2n+1}{n+1} xP_n(x) - \frac{n}{n+1} P_{n-1}(x)
$$ 

The orthonormal functions are 

$$
\phi_n(x) = \sqrt{\frac{2n+1}{2}} P_n(x)
$$

## Chebyshev polynomials

Take $[a,b]=[-1,1]$ and weight function

$$
w(x) = \frac{1}{\sqrt{1-x^2}}
$$ 

The resulting sequence of orthogonal polynomials can be written in terms of the Chebyshev polynomials

\begin{align}
T_0(x) &= 1 \\
T_1(x) &= x \\
T_{n+1}(x) &= 2 x T_n(x) - T_{n-1}(x), \quad n\ge 1
\end{align}

or 

$$
T_n(x) = \cos(n \cos^{-1}x)
$$ 

They are orthogonal wrt to the above weight 

$$
\ip{T_n, T_m} = \begin{cases}
0 & n \ne m\\
\pi & m = n = 0 \\
\half\pi & m = n > 0
\end{cases}
$$ 

Hence the orthonormal sequence is

$$
\phi_0(x) = \frac{1}{\sqrt{\pi}}, \qquad \phi_n(x) = \sqrt{\frac{2}{\pi}} T_n(x), \quad n \ge 1
$$

## Laguerre polynomials

Take $[a,b) = [0,\infty)$ and weight function

$$
w(x) = \ee^{-x}, \qquad x \in [0,\infty)
$$ 

The resulting sequence of orthogonal polynomials are called Laguerre polynomials and are usually written as

$$
L_n(x) = \frac{\ee^x}{n!} \dd{^n}{x^n}(x^n \ee^{-x}), \qquad n \ge 0
$$

which are orthonormal and hence $\phi_n(x) = L_n(x)$. Moreover they obey the recurrence relation

$$
L_{n+1}(x) = \frac{2n + 1 -x}{n+1} L_n(x) - \frac{n}{n+1} L_{n-1}(x)
$$

## Least squares approximation

Recall that we are trying to minimize $\norm{f - r}$ wrt $r \in \poly_n$, where the norm comes from a weighted inner product. The solution of this problem will be simplified if we express $r$ in terms of the orthonormal polynomials $\{ \phi_j \}$ corresponding to the chosen inner product. So let 

$$
r(x) = \sum_{j=0}^n b_j \phi_j(x)
$$ 

Then for a given $f \in \cts[a,b]$

$$
\norm{f-r}^2 = \int_a^b w(x) \left[ f(x) - \sum_{j=0}^n b_j \phi_j(x) \right]^2 \ud x =: G(b_0, b_1, \ldots,b_n)
$$

is to be minimized wrt the parameters $b_0, b_1, \ldots, b_n$. Now

$$
\begin{aligned}
0 
&\le G(b_0, b_1, \ldots, b_n) \\
&= \ip{ f - \sum_{j=0}^n b_j\phi_j, f - \sum_{i=0}^n b_i \phi_i} \\
&= \ip{f,f} - 2 \sum_{j=0}^n b_j \ip{f,\phi_j} + \sum_{i=0}^n \sum_{j=0}^n b_i b_j \underbrace{\ip{\phi_i, \phi_j}}_{\delta_{ij}} \\
&= \norm{f}^2 - 2 \sum_{j=0}^n b_j \ip{f,\phi_j} + \sum_{j=0}^n b_j^2 \\
&= \norm{f}^2 - \sum_{j=0}^n \ip{f,\phi_j}^2 + \sum_{j=0}^n [\ip{f,\phi_j} - b_j]^2
\end{aligned}
$$ 

The first two terms do not depend on the $\{b_j\}$ and the third term is non-negative. Hence $G$ is minimized iff the third term is zero which implies that

$$
b_j = \ip{f,\phi_j}, \qquad j=0,1,\ldots,n
$$ 

is the unique solution of the minimization problem. Hence the least squares approximation exists, is unique and given by

$$
r_n^*(x) = \sum_{j=0}^n \ip{f,\phi_j} \phi_j(x)
$$

+++

:::{prf:remark}
Note that we make use of the inner product structure of the norm in this proof. This is a special version of the projection theorem in Hilbert spaces, see e.g. ([@Brezis2010], Theorem 5.2, Corollary 5.4) or ([@Kreyszig1978], Section 3.3-1 and 3.3-2).
:::

+++

:::{prf:remark} Recursive property
If $r_n^*$ is the best approximation of degree $n$ then the best approximation of degree $n+1$ is

$$
r_{n+1}^*(x) = r_n^*(x) + \ip{f,\phi_{n+1}} \phi_{n+1}(x)
$$
:::

+++

:::{prf:remark} Stability
From the above proof, we see that

$$
\norm{f - r_n^*}^2 = \norm{f}^2 - \sum_{j=0}^n \ip{f,\phi_j}^2 =  \norm{f}^2 - \norm{r_n^*}^2
$$

so that 

$$
\norm{r_n^*} \le \norm{f}
$$
:::

+++

:::{prf:remark}
We can also find the minimizer using Calculus. The necessary condition for an extremum is

$$
\df{G}{b_i} =-2 \int_a^b w(x) \left[ f(x) - \sum_{j=0}^n b_j \phi_j(x) \right] \phi_i(x) \ud x = 0
$$

Using orthonormality, we get the unique solution of the optimality conditions as

\begin{align}
\int_a^b w(x) f(x) \phi_i(x) \ud x &= \sum_{j=0}^n b_j \int_a^b w(x) \phi_i(x) \phi_j(x) \ud x = b_i \\
\ip{f, \phi_i} &= b_i
\end{align}

The Hessian

$$
\df{^2 G}{b_i \partial b_j} = 2 \int_a^b w(x) \phi_i(x) \phi_j(x) \ud x = 2 \delta_{ij}
$$

is positive definite, so the extremum is the unique minimum.
:::

+++

As $n \to \infty$ the best approximations in 2-norm converge to the function.

:::{prf:theorem}
Assume $[a,b]$ is finite, $f \in \cts[a,b]$. Then

$$
\lim_{n \to \infty} \norm{f - r_n^*} = 0
$$
:::

:::{prf:proof}
Since $r_n^*$ is the minimizer

$$
\norm{f - r_n^*} = \min_{r \in \poly_n} \norm{f - r}
$$ 

we have $\norm{f - r_n^*} \ge \norm{f - r_{n+1}^*}$ since $\poly_n \subset \poly_{n+1}$. Hence $\norm{f - r_n^*}$ is a non-increasing sequence. But does it decrease to zero ?

Let $\epsilon > 0$ be arbitrary. By Weirstrass theorem, there is a polynomial $Q_m$ of some degree $m \ge 0$ for which

$$
\max_{x \in [a,b]}|f(x) - Q_m(x)| \le \frac{\epsilon}{c}, \qquad c = \left( \int_a^b w(x) \ud x \right)^\half
$$

By optimality of $r_n^*$ 

$$
\begin{aligned}
\norm{f - r_m^*} 
&\le \norm{f - Q_m} \\
&= \left( \int_a^b w(x)|f(x) - Q_m(x)|^2 \ud x \right)^\half \\
&\le \left( \int_a^b w(x)  \frac{\epsilon^2}{c^2} \ud x \right)^\half \\
&= \epsilon
\end{aligned}
$$ 

Hence

$$
\norm{f - r_n^*} \le \epsilon, \qquad \forall n \ge m
$$ 

Since $\epsilon$ was arbitrary this proves the theorem.
:::

+++

:::{prf:remark}
We proved convergence in $L^2$ norm but this does not imply that $\norm{f - r_n^*}_\infty \to 0$. But if additional assumptions on differentiability of $f$ are made, then we can prove that $r_n^*$ converges pointwise to $f$. Also the requirement of continuity is also not needed for convergence in 2-norm, see next section.
:::

+++

## Chebyshev least squares approximation

The best approximation in this basis is given by

$$
C_n(x) = \sum_{j=0}^n{}' a_j T_j(x), \qquad a_j = \frac{2}{\pi} \int_{-1}^1 \frac{f(x) T_j(x)}{\sqrt{1-x^2}} \ud x
$$

where the prime indicates that the first term must be multiplied by $\half$. Let us make the change of variable

$$
x = \cos\theta, \qquad \theta \in [0,\pi]
$$ 

Then $T_j(x) = \cos(j\cos^{-1}x) = \cos(j\theta)$ and

$$
C_n(\cos\theta) = \sum_{j=0}^n{}' a_j \cos(j\theta), \qquad a_j = \frac{2}{\pi} \int_0^\pi \cos(j\theta) f(\cos\theta) \ud\theta
$$

---

Define

$$
F(\theta) = f(\cos\theta), \qquad \theta \in [0,\pi]
$$ 

and extend it to $[-\pi,0]$ as an even function

$$
F(\theta) = f(\cos\theta), \qquad \theta \in [-\pi,\pi]
$$

The Chebyshev series becomes a Fourier cosine series, compare to [](#eq:FourSn), [](#eq:Fourakbk),

$$
S_n F(\theta) = \sum_{j=0}^n{}' a_j \cos(j\theta), \qquad a_j = \frac{1}{\pi} \int_{-\pi}^\pi \cos(j\theta) F(\theta) \ud\theta
$$

Since $F(\theta)$ is even, the sine terms vanish, and this is also the full Fourier series. 

---

Now if $f(x)$ is piecewise continuous, so is $F(\theta)$, and we know that

$$
|a_j| = \order{\frac{1}{j}}, \qquad j \to \infty
$$

In this case, we have $\norm{f - C_n}_2 \to 0$.

---

If $f \in \cts^{\nu-1}[-1,1]$ for some $\nu \ge 1$ and $f^{(\nu)}$ is piecewise continuous, then $F \in \cts^{\nu-1}_p[-\pi,\pi]$ and $F^{(\nu)}$ is piecewise continuous. This means that

$$
|a_j| = \order{\frac{1}{j^{\nu+1}}}, \qquad j \to \infty
$$

In this case, we also get $\norm{f - C_n}_\infty \to \infty$.

---

Note that we do not require $f^{(s)}(-1) = f^{(s)}(1)$; because of the way $F(\theta)$ is constructed, it satisfies $F^{(s)}(-\pi) = F^{(s)}(\pi)$ if $f \in \cts^s[-1,1]$.

+++

## Chebyshev least squares and minimax

Suppose $f$ is very smooth, say $f^{(\nu)}$ of bounded variation for some $\nu \gg 1$, then $|c_j| = \order{1/j^{\nu+1}}$ and

$$
f(x) - C_n(x) = \sum_{j=n+1}^\infty a_j T_j(x) \approx a_{n+1} T_{n+1}(x)
$$

Now 

$$
|T_{n+1}(x)| \le 1, \qquad x \in [-1,+1]
$$ 

and at the $n+2$ points

$$
x_j = \cos\left( \frac{j\pi}{n+1} \right), \quad j=0,1,2,\ldots,n+1
$$

it takes extreme values

$$
T_{n+1}(x_j) = (-1)^j
$$ 

The error

$$
f(x_j) - C_n(x_j)  \approx a_{n+1} (-1)^j
$$ 

is approximately equi-oscillating at $n+2$ points. This indicates that the Chebyshev least squares approximation can be expected to be close to the minimax
approximation, $C_n \approx p_n^\star$.

+++

:::{prf:example}
For the function $f(x) = \ee^x, \quad x \in [-1,1]$

\begin{align}
f(x) - C_3(x) &\approx a_4 T_4(x), \qquad a_4= 0.00547 \\
\norm{f - C_3}_\infty &= 0.00607 \\
E_3(f) &= \norm{f - p_3^*}_\infty = 0.00553 \approx a_4
\end{align}
:::

+++

:::{prf:theorem}
$$
\norm{f - p_n^*}_\infty \le \norm{f - C_n}_\infty \le \left( 4 + \frac{4}{\pi^2} \log n \right) \norm{f - p_n^*}_\infty
$$
:::

+++

:::{prf:remark}
If $f \in \cts^r[a,b]$ then

$$
\norm{f - p_n^*}_\infty \le \frac{d_{r-1}}{n^r} \norm{f^{(r)}}_\infty
$$ 

Hence

$$
\norm{f - C_n}_\infty \le \frac{d_{r-1}}{n^r} \left( 4 + \frac{4}{\pi^2} \log n \right) \norm{f^{(r)}}_\infty
$$
:::

+++

## Chebyshev interpolation and minimax

If $f$ is very smooth, then 

$$
f(x) - C_n(x) \approx a_{n+1} T_{n+1}(x)
$$

The error at the $n+1$ roots of $T_{n+1}$ given by (Chebyshev points of first kind)

$$
x_j = \cos\left( \frac{2j + 1}{2n + 2} \pi \right), \qquad j=0,1,2,\ldots,n
$$

is nearly zero. 

Let us construct the polynomial $I_n(x)$ which interpolates $f(x)$ at the $n+1$ roots of $T_{n+1}(x)$. Then

$$
\begin{aligned}
I_n(x) 
&= \sum_{j=0}^n f(x_j) \ell_j(x) \\
&= \sum_{j=0}^n C_n(x_j) \ell_j(x) + \sum_{j=0}^n \underbrace{[f(x_j) - C_n(x_j)]}_{\approx a_{n+1}T_{n+1}(x_j) = 0} \ell_j(x) \\
&\approx \sum_{j=0}^n C_n(x_j) \ell_j(x) \\
&= C_n(x) \\
&\approx p_n^*(x)
\end{aligned}
$$ 

Thus the Chebyshev interpolation at first kind points can be expected to be close to the best approximation.  We observed this in [](#p3minimaxexp). Recall from [](#sec:optimalintpnodes) that these Chebyshev nodes minimize the node polynomial appearing in polynomial interpolation error formula.

+++

:::{prf:theorem}

$$
\norm{f - I_n}_\infty \le \left( 2 + \frac{2}{\pi} \log(n+1) \right) \norm{f - p_n^*}_\infty
$$
:::

:::{prf:example}
For $n=100$ 

$$
\norm{f - I_{100}} \le 3.275 \norm{f - p_n^*}_\infty
$$

The Chebyshev-I interpolation error is at worst about three times larger than the best approximation error. But for a smooth function $\norm{f - p_n^*}_\infty$ itself will be very small, in which case interpolation is also very accurate.

Clearly interpolation is very easy compared to the Remez algorithm for finding the best approximation.
:::
