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
math:
  '\ud': '\textrm{d}'
  '\dd': '\frac{\ud #1}{\ud #2}'
  '\cond': '\textrm{cond}'
  '\half': '\frac{1}{2}'
  '\norm': '\|#1\|'
---

# Interpolation

```{code-cell}
from pylab import *
```

Suppose we know the values of a function $y=f(x)$ at a set of discrete points, say in a tabular form

::: center
    x       f(x)
  ------- ----------
    x_0     f(x_0)
    x_1     f(x_1)
     ⋮        ⋮
    x_N     f(x_N)
:::

e.g., table of logarithms, sine, cosine, or perhaps the results of some experimental measurement. We may want to know the value of the function at some $x$ that is not present in the table. This is an interpolation problem.

Using the available information, we will first construct an *approximation* $\tilde{f} (x)$ to the unknown function $f(x)$. Then

1.  We can use $\tilde{f}$ to evaluate an approximation to $f$ at any
    value $x$ $$f(x) \approx \tilde{f}(x)$$

2.  We can use $\tilde{f}$ to compute approximations to derivatives or
    integrals of $f$
    $$
    \dd{}{x}f(x) \approx \dd{}{x}\tilde{f}(x), \qquad \int_a^b f(x) \ud x \approx \int_a^b \tilde{f}(x)\ud x
    $$

To simplify the formulae, let us write 

$$
y_i = f(x_i)
$$

:::{prf:example} Linear interpolation

Given $(x_0, y_0)$, $(x_1, y_1)$, the linear interpolation is given by

$$
p(x) = \frac{x - x_1}{x_0 - x_1} y_0 + \frac{x-x_0}{x_1-x_0} y_1
$$

Note that $$p(x_0) = y_0, \qquad p(x_1) = y_1$$ i.e., $p$ interpolates the given data.
:::

:::{prf:example} Quadratic interpolation
Given $(x_0, y_0)$, $(x_1, y_1)$, $(x_2,y_2)$ we can try to find a
quadratic polynomial $$p(x) = a_0 + a_1 x + a_2 x^3$$ that interpolates
the given data 

$$
\begin{aligned}
a_0 + a_1 x_0 + a_2 x_0^2 &=& y_0 \\
a_0 + a_1 x_1 + a_2 x_1^2 &=& y_1 \\
a_0 + a_1 x_2 + a_2 x_2^2 &=& y_2
\end{aligned}
$$ 

This is a linear system of 3 equations for the three unknown coefficients $a_0, a_1, a_2$ and can be written in matrix form

$$
\begin{bmatrix}
1 & x_0 & x_0^2 \\
1 & x_1 & x_1^2 \\
1 & x_2 & x_2^2
\end{bmatrix}
\begin{bmatrix}
a_0 \\ a_1 \\ a_2 \end{bmatrix}
=
\begin{bmatrix}
y_0 \\ y_1 \\ y_2 \end{bmatrix}
$$ 

The matrix $V$ in the above equation is called *Vandermonde* matrix. The solution exists if this matrix is invertible. Its determinant is 

$$
\det V = \det\begin{bmatrix}
1 & x_0 & x_0^2 \\
0 & x_1-x_0 & x_1^2 - x_0^2 \\
0 & x_2-x_0 & x_2^2 - x_0^2
\end{bmatrix} = (x_0-x_1)(x_1-x_2)(x_2-x_0)
$$ 

The determinant is non-zero iff $x_0,x_1,x_2$ are distinct points which is a reasonable restriction.

**Another solution.** We can write a quadratic polynomial in other ways, e.g., by shifting the origin to $x_0$, we can write

$$
q(x) = b_0 + b_1 (x-x_0) + b_2(x-x_0)^2
$$ 

Solving the interpolation problem with this form, we immediately see that $b_0 = y_0$ and the remaining two conditions lead to 

$$
\begin{bmatrix}
x_1 - x_0 & (x_1-x_0)^2 \\
x_2-x_0 & (x_2-x_0)^2 \end{bmatrix}
\begin{bmatrix}
b_1 \\ b_2 \end{bmatrix} =
\begin{bmatrix}
y_1 - y_0 \\ y_2 - y_0 \end{bmatrix}
$$ 

The determinant of the $2 \times 2$ matrix is $(x_0-x_1)(x_1-x_2)(x_2-x_0)$ which is same as above, and is non-zero if all points are distinct.

**Uniqueness.** We have two quadratic polynomials $p(x)$, $q(x)$ which solve the same interpolation problem. Are the two polynomials identical ? Define

$$
r(x) = p(x) - q(x)
$$ 

Then 

$$
r(x_0) = r(x_1) = r(x_2) = 0
$$ 

The quadratic polynomial vanishes at three distinct points which means that it has three distinct roots. But a quadratic polynomial can have atmost two distinct roots, which implies that $r$ must be the zero polynomial, i.e. $p(x) = q(x)$ for all $x$.
:::

In fact this is true for any degree of the interpolating polynomial.  Thus the interpolation problem leads to a unique polynomial, and does not depend on the choice of the basis function for the space of polynomials.

:::{prf:theorem}
If $p(x),q(x)$ are polynomials of degree $\le N$ which agree at $\ge N+1$ distinct points, then $p(x) = q(x)$.
:::

## Some questions

1.  What basis functions to choose ?

    -   polynomials: Easy to evaulate. But there are several choices,
        monomials, Lagrange, Newton, Chebychev, Legendre, etc.

    -   rational polynomials

    -   trigonometric functions: sines and cosines

    -   radial basis functions

2.  The basis functions may have compact support or be globally
    supported. Using globally supported functions gives high accuracy
    for smooth functions and leads to what are called spectral methods.
    Compactly supported basis functions are useful when functions are
    less regular, and are used in finite element methods.

3.  How to sample the data, uniformly or otherwise ?

4.  How to interpolate/approximate arbitrarily scattered data,
    especially in high dimensions ?

5.  Interpolation/approximation can be performed using only function
    values, or using derivative information when available, leading to
    Hermite methods.

6.  Should we exactly fit the data as in interpolation, or approximate
    it, e.g., using a least squares fit ? This is also a question of
    norms to measure the approximation error. Usually the data has some
    noise and then it may not be a good idea to interpolate.

## Interpolation using polynomials

Polynomials are the simplest functions we can use since they require only basic arithmetic operations. They can provide a good approximation as shown by Weirstrass Theorem.

### Interpolation using monomials

Given $N+1$ data $(x_i, y_i)$, $i=0,1,\ldots,N$ with $y_i = f(x_i)$, we
can try to find a polynomial of degree $N$

$$
p(x) = a_0 + a_1 x + \ldots + a_N x^N
$$ 

that satisfies the interpolation conditions $$p(x_i) = y_i, \qquad i=0,1,\ldots, N$$ This is a linear system of $N+1$ equations that can be written in matrix form

$$
\underbrace{\begin{bmatrix}
1 & x_0 & x_0^2 & \ldots & x_0^N \\
1 & x_1 & x_1^2 & \ldots & x_1^N \\
\vdots & \vdots & \ldots & \vdots \\
1 & x_N & x_N^2 & \ldots & x_N^N
\end{bmatrix}}_{V}
\begin{bmatrix}
a_0 \\ a_1 \\ \vdots \\ a_N \end{bmatrix} =
\begin{bmatrix}
y_0 \\ y_1 \\ \vdots \\ y_N \end{bmatrix}
$$ 

This has a unique solution if the Vandermonde matrix $V$ is non-singular. Since

$$
\det V = \prod_{j=0}^N \prod_{k=j+1}^N (x_k - x_j)
$$ 

we can solve the problem provided the points $\{ x_i \}$ are distinct.

**$V$ is non-singular.**
We can show this without computing the determinant. Assume that the $\{ x_i \}$ are distinct. It is enough to show that the only solution of $Va=0$ is $a=0$. Note that the set of $N+1$ equations $Va=0$ is of the form $$p(x_i) = 0, \qquad i=0,1,\ldots,N$$ which implies that $p(x)$ has $N+1$ distinct roots. But since $p$ is a polynomial of degree $N$, this implies that $p(x) \equiv 0$ and hence each $a_i = 0$.

### Condition number of $V$

If we want to solve the interpolation problem by solving the matrix problem, then the condition number of the matrix becomes important. Matrices with large condition numbers cannot solved accurately on a computer due to blow-up in round-off errors.

:::{prf:example}
Take $(x_0,x_1,x_2) = (100,101,102)$. Then 

$$
V = \begin{bmatrix}
1 & 100 & 10000 \\
1 & 101 & 10201 \\
1 & 102 & 10402 \end{bmatrix}, \qquad \cond(V) \approx 10^8
$$ 

If we scale the $x$ as $x \to \frac{x}{x_0}$, then

$$
\tilde{V} = \begin{bmatrix}
1 & 1.00 & 1.0000 \\
1 & 1.01 & 1.0201 \\
1 & 1.02 & 1.0402 \end{bmatrix}, \qquad \cond(\tilde{V}) \approx 10^5
$$

Or we can shift the origin to $x_0$ 

$$
\hat{V} = \begin{bmatrix}
1 & 0 & 0 \\
1 & 1 & 1 \\
1 & 2 & 4 \end{bmatrix}, \qquad \cond(\hat{V}) \approx 14
$$ 

The condition number can be improved by scaling and/or shifting the variables.
:::

:::{prf:remark}
It is usually better to map the data to the interval $[-1,+1]$ or $[0,1]$ and then solve the interpolation problem on the mapped interval.
Assuming that $x_0 < x_1 < \ldots
x_N$, define 

$$
\xi = \frac{x - x_0}{x_N - x_0}
$$ 

and find a polynomial of the form 

$$
p(\xi) = a_0 + a_1 \xi + \ldots + a_N \xi^N
$$ 

But still the condition number increases rapidly with $N$. For $N=20$ we have a condition number of $8 \times 10^8$ even when using the interval $[-1,+1]$.
:::

:::{prf:example}
We will use the [`cond`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.cond.html) function from `numpy.linalg` to compute the condition number. By default, it uses the 2-norm.

```{code-cell}
Nvalues, Cvalues = [], []
for N in range(1,30):
    x = linspace(-1.0,+1.0,N+1)
    V = zeros((N+1,N+1))
    for j in range(0,N+1):
        V[j][:] = x**j # This is transpose of V as defined above.
    Nvalues.append(N), Cvalues.append(cond(V.T))
semilogy(Nvalues, Cvalues, 'o-')
xlabel('N'), ylabel('cond(V)'), grid(True);
```

The condition number is large for even moderate value of $N=30$.
:::

### Lagrange interpolation

Lagrange interpolation provides the solution without having to solve a
matrix problem. Define

$$
\pi_i(x) = (x-x_0) \ldots (x-x_{i-1})(x-x_{i+1})\ldots (x-x_N)
$$ 

and

$$
\ell_i(x) = \frac{\pi_i(x)}{\pi_i(x_i)}
$$ 

Note that each $\ell_i$ is a
polynomial of degree $N$ and $\ell_i(x_j) = \delta_{ij}$, i.e.,

$$
\ell_i(x_i) = 0, \qquad \ell_i(x_j) = 0, \quad j \ne i
$$ 

Then consider
the polynomial of degree $N$ given by

$$
p_N(x) = \sum_{j=0}^N y_j \ell_j(x)
$$ 

By construction this satisfies

$$
p_N(x_i) = y_i, \qquad i=0,1,\ldots,N
$$ 

and hence is the solution to
the interpolation problem.

:::{prf:example}
Initially, let us use some python functions to compute the interpolating
polynomial. In particular, we use `polyfit` and `polyval` functions. We
will demonstrate this for the function

$$
f(x) = \cos(4\pi x), \qquad x \in [0,1]
$$

Note that polyfit actually performs a least squares fit, but we give $N+1$ points and ask for a degree $N$ polynomial, which gives an interpolating polynomial.

Define the function

```{code-cell}
xmin, xmax = 0.0, 1.0
f = lambda x: cos(4*pi*x)
```

Make a grid and evaluate function at those points.

```{code-cell}
N = 8 # degree, we need N+1 points
x = linspace(xmin, xmax, N+1)
y = f(x)
```

Construct polynomial using polyfit


```{code-cell}
p = polyfit(x,y,N)
```

Now we evaluate this on larger grid for better visualization


```{code-cell}
M = 100
xe = linspace(xmin, xmax, M)
ye = f(xe) # exact function
yp = polyval(p,xe)

plot(x,y,'o',xe,ye,'--',xe,yp,'-')
legend(('Data points','Exact function','Polynomial'))
title('Degree '+str(N)+' interpolation');
```

:::

:::{prf:example}
Interpolate the following two functions on uniformly spaced points

$$
f(x) = \cos(x), \qquad x \in [0,2\pi]
$$ 

and

$$
f(x) = \frac{1}{1+16x^2}, \qquad x \in [-5,+5]
$$

```{code-cell}
xmin, xmax = 0.0, 2.0*pi
fun = lambda x: cos(x)

xx = linspace(xmin,xmax,100);
ye = fun(xx);

figure(figsize=(9,8))
for i in range(1,7):
    N = 2*i;
    subplot(3,2,i)
    x = linspace(xmin,xmax,N+1);
    y = fun(x);
    P = polyfit(x,y,N);
    yy = polyval(P,xx);
    plot(x,y,'o',xx,ye,'--',xx,yy);
    axis([xmin, xmax, -1.1, +1.1])
    text(3.0,0.0,'N='+str(N))
```

The interpolating polynomials converge to the true function as $N$ increases.
:::

:::{prf:theorem}
Suppose $f(x)$ is a polynomial of degree $\le N$ and we interpolate this
at $N+1$ distinct points to construct a polynomial $p(x)$ of degree
$\le N$. Then $f(x) \equiv p(x)$.
:::

:::{prf:proof}
Define 

$$
r(x) = f(x) - p(x)
$$ 

which is of degree $\le N$.
$r(x)$ vanishes at $N+1$ distinct nodes and hence it must be zero
polynomial.
:::

:::{prf:example}
If $f(x) = a + b x$ and we interpolate this with $p(x) = a_0 + a_1 x  + a_2 x^2$ at three distinct points, then the solution will be 

$$
a_0 = a, \qquad a_1 = b, \qquad a_2 = 0
$$

:::

:::{exercise}
If $f(x)$ is a polynomial of degree $m$ and $p(x)$ interpolates it at $n+1$ distinct points with $n > m$, then show that $p(x)$ is actually a polynomial of degree $m$.
:::

## Error in polynomial approximation

:::{prf:theorem}
Let $p_N(x)$ be the degree $N$ polynomial which interpolates the given
data $(x_i,y_i)$, $i=0,1,\ldots,N$. Then the error

$$
e_N(x) = f(x) - p_N(x)
$$ 

is bounded by

$$
|e_N(x)| \le \frac{M_{N+1}}{(N+1)!} |(x-x_0)(x-x_1)\ldots(x-x_N)|
$$

:::

:::{prf:proof}
Assume for simplicity that $x_0 < x_1 < \ldots < x_N$. Since

$$
e_N(x_i) = 0, \qquad i=0,1,\ldots,N
$$ 

we can write the error as

$$
e_N(x) = f(x) - p_N(x) = (x-x_0)(x-x_1)\ldots(x-x_N) K(x)
$$ 

for some function $K(x)$. Choose an arbitrary $x_* \in [x_0,x_N]$, different from all the $x_i$ and define

$$
\Phi(x) = f(x) - p_N(x) - (x-x_0)(x-x_1)\ldots(x-x_N) K(x_*)
$$ 

If $f(x)$ has $(N+1)$ derivatives, then we can differetiate $N+1$ times to get 

$$
\Phi^{(N+1)}(x) = f^{(N+1)}(x) - (N+1)! K(x_*)
$$ 

Clearly, $\Phi(x)$ vanishes at $N+2$ points, $\{x_0,x_1,\ldots,x_N,x_*\}$. By mean value theorem, $\Phi'(x)$ vanishes at atleast $N+1$ points in the interval $[x_0,x_N]$, $\Phi''(x)$ vanishes at atleast $N$ points in $[x_0,x_N]$, etc. Continuing this way $\Phi^{(N+1)}(x)$ vanishes in atleast one point $\bar{x} \in [x_0,x_N]$

$$
f^{(N+1)}(\bar{x}) - (N+1)! K(x_*) = 0 \quad\Longrightarrow\quad K(x_*) = \frac{1}
{(N+1)!} f^{(N+1)}(\bar{x})
$$ 

and hence

$$
f(x_*) = p_N(x_*) + \frac{1}{(N+1)!} (x_*-x_0)(x_*-x_1)\ldots(x_*-x_N) f^{(N+1)}
(\bar{x})
$$ 

Since $x_*$ was arbitrary, we can call it $x$ and we obtain the error formula

$$
e_N(x) = f(x)  - p_N(x) = \frac{1}{(N+1)!} (x-x_0)(x-x_1)\ldots(x-x_N) f^{(N+1)}
(\bar{x})
$$ 

Defining 

$$
M_{N+1} = \max_{x \in [x_0,x_N]}|f^{(N+1)}(x)|
$$

we obtain the error bound.
:::

## Uniformly spaced points

Let us specialize the error estimate to the case of uniformly spaced
points in the interval $[a,b]$ with

$$
x_i = a + i h, \qquad 0 \le i \le N, \qquad h = \frac{b-a}{N}
$$ 

We have $x_0 = a$ and $x_N = b$.

:::{prf:theorem}
For any $x \in [a,b]$

$$
\prod_{i=0}^N |x - x_i| \le \frac{1}{4} h^{N+1} N!
$$
:::

:::{prf:proof}
Fix an $x$ and find $j$ such that $x_j \le x \le x_{j+1}$. Show
that[^1] 

$$
|(x-x_j)(x-x_{j+1})| \le \frac{1}{4}h^2
$$ 

Hence

$$
\prod_{i=0}^N |x - x_i| \le \frac{h^2}{4} \prod_{i=0}^{j-1}(x - x_i) \prod_{i=j+2}^N
(x_i - x)
$$ 

Now since $x \le x_{j+1}$ and $-x \le -x_j$

$$
\prod_{i=0}^N |x - x_i| \le \frac{h^2}{4} \prod_{i=0}^{j-1}(x_{j+1} - x_i)
\prod_{i=j+2}^N (x_i - x_j)
$$ 

Using

$$
x_{j+1} - x_i = (j-i+1)h, \qquad x_i - x_j = (i-j)h
$$ 

the inequality
becomes

$$
\prod_{i=0}^N |x - x_i| \le \frac{h^{N+1}}{4} \prod_{i=0}^{j-1}(j-i+1)
\prod_{i=j+2}^N (i - j) = \frac{h^{N+1}}{4} (j+1)! (N-j)!
$$ 

Finally show
that[^2] 

$$
(j+1)! (N-j)! \le N!, \qquad 0 \le j \le N-1
$$ 

which completes the proof of the theorem.
:::

:::{prf:theorem}
Assume that $|f^{(n)}(x)| \le M$ for all $x$ and $n$. Then the error of
polynomial interpolation on uniformly spaced points is bounded by

$$
|f(x) - p(x)| \le \frac{M h^{N+1}}{4(N+1)}
$$ 

and the error goes to zero as $N \to \infty$.
:::

:::{prf:proof}
The error bound follows from the two previous theorems. As $N$ increases, $h$ becomes less than one at some point, and beyond this, the right hand side in the error bound will converge to zero.
:::

:::{prf:remark}
Functions like $\cos x$, $\sin x$, $\exp(x)$ satisfy the conditions of the above theorem. These conditions are quite strong and can be relaxed considerably. E.g., if $f(x) = \sin(\alpha x)$ then $|f^{(n)}(x)| \le |\alpha|^n$ and if $|\alpha| > 1$ then the derivatives can increase. If the derivatives satisfy

$$
|f^{(n)}(x)| \le C |\alpha|^n
$$ 

then the error estimate gives

$$
|f(x) - p(x)| \le \frac{C (\alpha h)^{N+1}}{4(N+1)}
$$ 

As $N$ increases, we will satisfy $\alpha h < 1$ and beyond this point, the right hand side goes to zero.

(TODO) The problematic case is if the derivatives grow at a factorial rate 

$$
|f^{(n)}(x)| \le C n!
$$ 

so that error bound is

$$
|f(x) - p(x)| \le \frac{C N! (b-a)^{N+1}}{4 N^{N+1}}
$$ 

which does not go to zero.
:::

## Difficulty of polynomial interpolation

Do the polynomial approximations $p_N$ converge to the true function $f$ as $N \to \infty$ ? The error formula seems to suggest so, due to the factor $\frac{1}{(N+1)!}$ provided higher order derivatives of $f$ are bounded.  On uniformly spaced points, we have seen the interpolants of $\cos(x)$ converge but those of the rational function $\frac{1}{1+16x^2}$ do not.  This must be related to the behaviour of the derivatives of these functions.

:::{prf:example}
Consider $y=\ln(x)$. Then its derivatives are 

$$
\begin{aligned}
y' &=& \frac{1}{x} \\
y'' &=& -\frac{1}{x^2} \\
y''' &=& \frac{2!}{x^3} \\
&\vdots& \\
y^{(n)} &=& \frac{(-1)^{n-1} (n-1)!}{x^n}
\end{aligned}
$$ 

Even though the curve $y=\ln(x)$ looks smooth near any
value of $x$, as $n$ gets large, the derivatives become very large in
size, and tend to behave like $n!$ or worse.
:::

This is the general situation; for most functions, some higher order
derivatives tend to grow as $n!$. Even for a polynomial $p_N(x)$, the
derivatives grow in size until the $N$'th one, which is $a_N N!$, after
which they suddenly all become zero.

:::{prf:example} Runge phenomenon

Consider interpolating the following two functions on $[-1,1]$

$$
f_1(x) = \exp(-5x^2), \qquad f_2(x) = \frac{1}{1 + 16 x^2}
$$

We will try uniformly spaced points and Chebyshev points.

Let us first plot the two functions.


```{code-cell}
xmin, xmax = -1.0, +1.0

f1 = lambda x: exp(-5.0*x**2)
f2 = lambda x: 1.0/(1.0+16.0*x**2)

xx = linspace(xmin,xmax,100,True)
figure(figsize=(8,4))
plot(xx,f1(xx),xx,f2(xx))
legend(("$1/(1+16x^2)$", "$\\exp(-5x^2)$"));
```

The two functions look qualitatively similar and both are infinitely differentiable.

```{code-cell}
def interp(f,points):
    xx = linspace(xmin,xmax,100,True);
    ye = f(xx);

    figure(figsize=(9,8))
    for i in range(1,7):
        N = 2*i
        subplot(3,2,i)
        if points == 'uniform':
            x = linspace(xmin,xmax,N+1,True)
        else:
            theta = linspace(0,pi,N+1, True)
            x = cos(theta)
        y = f(x);
        P = polyfit(x,y,N);
        yy = polyval(P,xx);
        plot(x,y,'o',xx,ye,'--',xx,yy)
        axis([xmin, xmax, -1.0, +1.1])
        text(-0.1,0.0,'N = '+str(N))
```

Interpolate $f_1(x)$ on uniformly spaced points.

```{code-cell}
interp(f1,'uniform')
```

Interpolate $f_2(x)$ on uniformly spaced points.


```{code-cell}
interp(f2,'uniform')
```

The above results are not good. Let us try $f_2(x)$ on Chebyshev points.

```{code-cell}
interp(f2,'chebyshev')
```

What about interpolating $f_1(x)$ on Chebyshev points ?


```{code-cell}
interp(f1,'chebyshev')
```

This also seems fine. So uniform points sometimes works, Chebyshev points work all the time.
:::

## Polynomial factor in error bound

The error of polynomial interpolation is given by

$$
f(x) - p_N(x) = \frac{\omega_N(x)}{(N+1)!} f^{(N+1)}(\xi)
$$ 

where
$\xi \in I(x_0,x_1,\ldots,x_N,x)$ and

$$
\omega_N(x) = (x-x_0)\ldots(x-x_N)
$$

**Case $N=1$.**
In case of linear interpolation

$$
\omega_1(x) = (x-x_0)(x-x_1), \qquad h = x_1 - x_0
$$ 

Then

$$
\max_{x_0 \le x \le x_1} |\omega_1(x)| = \frac{h^2}{4}
$$ 

and the
interpolation error bound is

$$
\max_{x_0 \le x \le x_1} |f(x) - p_1(x)| \le \frac{h^2}{8} \max_{x_0 \le x \le x_1}|
f''(x)|
$$

**Case $N=2$.**
To bound $\omega_2(x)$ on $[x_0,x_2]$, shift the origin to $x_1$ so that

$$
\omega_2(x) = (x-h)x(x+h)
$$ 

Near the center of the data

$$
\max_{x_1-\half h \le x \le x_1 + \half h} |\omega_2(x)| = 0.375 h^3
$$

whereas on the whole interval

$$
\max_{x_0 \le x \le x_2} |\omega_2(x)| = \frac{2\sqrt{3}}{9}h^3 \approx 0.385 h^3
$$

In this case, the error is of similar magnitude whether $x$ is near the
center or near the edge.

**Case $N=3$.**
We can shift the nodes so that they are symmetric about the origin. Then

$$
\omega_3(x) = (x^2 - \frac{9}{4}h^2) (x^2 - \frac{1}{4}h^2)
$$ 

and

$$
\max_{x_1 \le x \le x_2}|\omega_3(x)| = \frac{9}{16}h^4 \approx 0.56 h^4
$$

$$
\max_{x_0 \le x \le x_3}|\omega_3(x)| = h^4
$$ 

In this case, the error
near the endpoints can be twice as large as the error near the middle.

**Case $N > 3$.**
The behaviour exhibited for $N=3$ is accentuated for larger degree. For
$N=6$, 

$$
\max_{x_2 \le x \le x_4}|\omega_3(x)| \approx 12.36 h^7, \qquad
\max_{x_0 \le x \le x_6}|\omega_3(x)| \approx 95.8 h^7
$$ 

and the error near the ends can be almost 8 times that near the center.

:::{prf:example}
Function $\omega_N(x)$ arising in interpolation error

For a given set of points $x_0, x_1, \ldots, x_N \in [-1,+1]$ we plot the function

$$
\omega_N(x) = (x-x_0)(x-x_1) \ldots (x-x_N), \qquad x \in [-1,+1]
$$

for uniformly spaced points.

```{code-cell}
def omega(x,xp):
    f = 1.0
    for z in xp:
        f = f * (x-z)
    return f

def plot_omega(x):
    M  = 1000
    xx = linspace(-1.0,1.0,M)
    f  = 0*xx
    for i in range(M):
        f[i] = omega(xx[i],x)
    plot(xx,f,'b-',x,0*x,'o')
    title("N = "+str(N));
    grid(True)
```


```{code-cell}
N = 3
x = linspace(-1.0,1.0,N+1)
plot_omega(x)
```


```{code-cell}
N = 4
x = linspace(-1.0,1.0,N+1)
plot_omega(x)
```

```{code-cell}
N = 6
x = linspace(-1.0,1.0,N+1)
plot_omega(x)
```

```{code-cell}
N = 20
x = linspace(-1.0,1.0,N+1)
plot_omega(x)
```

:::

## Distribution of data points

The error in polynomial interpolation is

$$
|f(x) - p_N(x)| \le \frac{|\omega_N(x)|}{(N+1)!} \max_{\xi} |f^{(N+1)}(\xi)|
$$

where 

$$
\omega_N(x) = (x-x_0)(x-x_1) \ldots (x-x_N)
$$ 

For a given function $f(x)$, we cannot do anything about the derivative term in the error estimate. For uniformly spaced data points, $\omega_N$ has large value near the end points which is also where the Runge phenomenon is observed. But we can try to minimize the magnitude of $\omega_N(x)$ by choosing a different set of nodes for interpolation. For the following discussion, let us assume that the $x_i$ are ordered and contained in the interval $[-1,+1]$.

Can we choose the nodes $\{x_i\}$ so that

$$
\max_{x \in [-1,+1]} |\omega_N(x)|
$$ 

is minimized ? We will show that

$$
\min_{\{x_i\}} \max_{x \in [-1,+1]} |\omega_N(x)| = 2^{-N}
$$ 

and the minimum is achieved for following set of nodes

$$
x_i = \cos\left( \frac{2(N-i)+1}{2N+2} \pi \right), \qquad i=0,1,\ldots,N
\label{eq:chebpts}
$$ 

These are called *Chebyshev points of first kind*.

## Chebyshev polynomials

The Chebyshev polynomials are defined on the interval $[-1,+1]$ and the
first few polynomials are 

$$
\begin{aligned}
T_0(x) &=& 1 \\
T_1(x) &=& x
\end{aligned}
$$ 

The remaining polynomials can be defined by recursion

$$
T_{n+1}(x) = 2 x T_n(x) - T_{n-1}(x)
$$ 

so that 

$$
\begin{aligned}
T_2(x) =& 2x^2 - 1 \\
T_3(x) =& 4x^3 - 3x \\
T_4(x) =& 8x^4 - 8x^2 + 1
\end{aligned}
$$ 

In Python, we can use function recursion to compute the Chebyshev polynomials.

```{code-cell}
def chebyshev(n, x):
    if n == 0:
        y = np.ones(len(x))
    elif n == 1:
        y = x.copy()
    else:
        y = 2*x*chebyshev(n-1,x) - chebyshev(n-2,x)
    return y
```

The first few of these are shown below.

```{code-cell}
:tag: hide-input
N = 200
x = np.linspace(-1.0,1.0,N)
figure(figsize=(8,6))
for n in range(0,6):
    y = chebyshev(n,x)
    plot(x,y)
grid(True), xlabel('x'), ylabel('$T_n(x)$');
```

We can write $x \in [-1,+1]$ in terms of an angle $\theta \in [0,\pi]$ 

$$
x=\cos\theta
$$

**Properties of Chebyshev polynomials**

-   $T_n(x)$ is a polynomial of degree $n$.

-   $T_n(\cos\theta) = \cos(n\theta)$. (Fourier cosine series)

-   $T_n(x) = \cos(n\cos^{-1}(x))$

-   $|T_n(x)| \le 1$

-   Extrema
    $$
    T_n\left(\cos\frac{j\pi}{n}\right) = (-1)^j, \qquad 0 \le j \le n
    $$

-   Roots
    $$
    T_n\left( \cos \frac{2j+1}{2n}\pi \right) = 0, \qquad 0 \le j \le n-1
    $$

:::{exercise}
Prove the above properties.
:::

:::{prf:definition} Monic polynomial
A polynomial whose term of highest degree has coefficient one is called a monic polynomial.
:::

:::{prf:remark}
In $T_n(x)$, the coefficient of $x^n$ is $2^{n-1}$ for $n \ge 1$. Hence
$2^{1-n} T_n(x)$ is monic polynomial of degree $n$.
:::

:::{prf:theorem}
If $p$ is a monic polynomial of degree $n$, then

$$
\norm{p}_\infty = \max_{-1 \le x \le +1} |p(x)| \ge 2^{1-n}
$$
:::

:::{prf:proof}
We prove by contradiction. Suppose that

$$
|p(x)| < 2^{1-n}, \qquad \forall |x| \le 1
$$ 

Define

$$
q(x) = 2^{1-n} T_n(x)
$$ 

The extrema of $T_n$ are at

$$
x_i = \cos(\frac{i\pi}{n}), \quad 0 \le i \le n
$$ 

and

$$
T_n(x_i) = (-1)^i \qquad \textrm{so that} \qquad q(x_i) = 2^{1-n} (-1)^i
$$

Now 

$$
(-1)^i p(x_i) \le |p(x_i)| < 2^{1-n} = (-1)^i q(x_i)
$$ 

so that

$$
(-1)^i [q(x_i) - p(x_i)] > 0, \qquad 0 \le i \le n
$$ 

$q-p$ changes sign $n$ times, so it must have atleast $n$ roots. But $q-p$ is of degree $n-1$, and hence we get a contradiction.
:::

## Optimal nodes

Since $\omega_N(x)$ is a monic polynomial of degree $N+1$, we know from
previous theorem that

$$
\max_{-1 \le x \le +1} |\omega_N(x)| \ge 2^{-N}
$$ 

Can we choose the
$x_i$ so that the minimum value of $2^{-N}$ is achieved ?

If $x_i$ are the $N+1$ distinct roots of $T_{N+1}(x)$, then

$$
\omega_N(x) = 2^{-N} T_{N+1}(x)
$$ 

so that

$$
\max_{-1 \le x \le +1} |\omega_N(x)| = 2^{-N} \max_{-1 \le x \le +1} |T_{N+1}(x)| =
2^{-N}
$$ 

The optimal nodes are given by

$$
x_i = \cos\left( \frac{2i+1}{2N+2} \pi \right), \qquad 0 \le i \le N
$$

In terms of $\theta$ variable, the spacing between these nodes is

$$
\theta_{i+1} - \theta_i = \frac{\pi}{N+1}
$$

:::{prf:example}
The roots of degree $n$ Chebyshev polynomial $T_n(x)$ are

$$
x_i = \cos\left( \frac{2i+1}{2n} \pi \right), \qquad i=0,1,\ldots,n-1
$$

The roots are shown below for $n=10,11,\ldots,19$.

```{code-cell}
c = 1
for n in range(10,20):
    j = np.linspace(0,n-1,n)
    theta = (2*j+1)*np.pi/(2*n)
    x = np.cos(theta)
    y = 0*x
    subplot(10,1,c)
    plot(x,y,'.')
    xticks([]); yticks([])
    ylabel(str(n))
    c += 1
```

Note that the roots are clustered at the end points.
:::

:::{prf:theorem}
If the nodes $\{x_i\}$ are the $N+1$ roots of the Chebyshev polynomial
$T_{N+1}$, then the error formula for polynomial interpolation in the
interval $[-1,+1]$ is

$$
|f(x) - p_N(x)| \le \frac{1}{2^N (N+1)!} \max_{|t| \le 1}|f^{(N+1)}(t)|
$$
:::

:::{prf:remark}
Note that

$$
x_0 = \cos\left( \frac{1}{2N+2} \pi \right) < 1, \qquad x_N = \cos\left( \frac{2N+1} {2N+2} \pi \right) > -1
$$ 

and the endpoints of $[-1,+1]$ are not nodes.  The nodes are ordered as $x_0 > x_1 > \ldots > x_N$. We can reorder them by defining the $x_i$ as in [](#eq:chebpts).
:::

:::{prf:remark}
In practice, we dont use the Chebyshev nodes as derived above. The important point is how the points are clustered near the ends of the interval. This type of clustering can be achieved by other node sets. If we want $N+1$ nodes, then divide $[0,\pi]$ into $N$ equal intervals so that 

$$
\theta_i = \frac{(N-i)\pi}{N}, \qquad i=0,1,\ldots,N
$$ 

and then the nodes are given by 

$$
x_i = \cos\theta_i
$$ 

These are called *Chebyshev points of second kind*.

FIXME

Chebyshev points are defined in the interval $[-1,+1]$. They are obtained as projections of uniformly spaced points on the unit circle onto the $x$-axis. For any $N \ge 2$ they are defined as follows

$$
\theta_i = \frac{\pi i}{N-1}, \qquad x_i = \cos(\theta_i), \qquad i=0,1,\ldots,N-1
$$

```{code-cell}
t = linspace(0,pi,1000)
xx, yy = cos(t), sin(t)
plot(xx,yy)

N = 10
theta = linspace(0,pi,N)
plot(cos(theta),sin(theta),'o')
for i in range(N):
    x1 = [cos(theta[i]), cos(theta[i])]
    y1 = [0.0, sin(theta[i])]
    plot(x1,y1,'k--',cos(theta[i]),0,'sr')
axis([-1.1, 1.1, 0.0, 1.1])
axis('equal'), title(str(N)+' Chebyshev points');
```

Below, we compare the polynomial $\omega_N(x)$ for uniform and Chebyshev points for $N=16$.

```{code-cell}
M  = 1000
xx = linspace(-1.0,1.0,M)

N = 16
xu = linspace(-1.0,1.0,N+1)    # uniform points
xc = cos(linspace(0.0,pi,N+1)) # chebyshev points
fu = 0*xx
fc = 0*xx
for i in range(M):
    fu[i] = omega(xx[i],xu)
    fc[i] = omega(xx[i],xc)
plot(xx,fu,'b-',xx,fc,'r-')
legend(("Uniform","Chebyshev"))
grid(True), title("N = "+str(N));
```

:::

## Back to Lagrange interpolation

The Lagrange interpolation polynomial of degree $N$ is

$$
p(x) = \sum_{j=0}^N f_j \ell_j(x) \qquad \textrm{where} \qquad
\ell_j(x) = \frac{\prod\limits_{i=0,i\ne j}^N (x-x_i)}{\prod\limits_{i=0,i \ne j}^N(x_j - x_i)}
$$ 

The denominators in $\ell_j$ can be pre-computed and stored.  Then for every $x$, the evaluation of $p(x)$ requires $O(N^2)$ flops.  Ideally, we desire a method which takes $O(N)$ flops. Also, if we add a new data pair $(x_{N+1},f_{N+1})$, then this requires a new computation from scratch. Moreover, the numerical evaluation is also unstable to round-off errors[^3].

### Improved Lagrange formula

Define 

$$
\ell(x) = (x-x_0)(x-x_1) \ldots (x-x_N)
$$

$$
w_j = \frac{1}{\prod\limits_{k=0, k \ne j}^N (x_j - x_k)} = \frac{1}{\ell'(x_j)}, \qquad j=0,1,\ldots,N
$$ 

The Lagrange polynomials can be written as

$$
\ell_j(x) = \ell(x) \frac{w_j}{x-x_j}
$$ 

and the interpolant is given by 

$$
p(x) = \ell(x) \sum_{j=0}^N \frac{w_j}{x-x_j} f_j
$$ 

This is called the *first form of barycentric interpolation formula* and is due to Rutishauser. For each $x$

-   calculation of $w_j$ requires $O(N^2)$ flops, and these are
    independent of $x$

-   calculation of $\ell(x)$ requires $O(N)$ flops

-   calculation of the other term requires $O(N)$ flops

so the complexity of this formula is $O(N)$. Once the weights $w_j$ are known, this permits interpolation of as many functions as desired in $O(N)$ operations. If a new data point $(x_{N+1},f_{N+1})$ is added

-   Divide each $w_j$ by $(x_j - x_{N+1})$ which needs $O(N+1)$ flops

-   compute $w_{N+1}$ at $O(N+1)$ flops

The total cost of updating the formula is $O(N)$ flops.

### Barycentric formula

Since $\ell_j$ form a partition of unity[^4]

$$
1 = \sum_{j=0}^N \ell_j(x) = \ell(x) \sum_{j=0}^N \frac{w_j}{x-x_j}
$$

the interpolant can be written as

$$
p(x) = \frac{ \sum_{j=0}^N \frac{w_j}{x-x_j} f_j }{ \sum_{j=0}^N \frac{w_j}{x-
x_j} }
$$ 

which is called the *second form of barycentric formula*. This has same complexity as the previous form. Higham has shown that barycentric formula is numerically stable. This code is general since it takes any node sets and computes the weights. For special node distributions, we can simplify the weight calculation as in next two sections.

### Equispaced points

If we have equispaced data points in $[-1,+1]$, with spacing $h = \frac{2}{N}$, the weights are given by

$$
w_j = (-1)^{N-j} \binom{N}{j} \frac{1}{h^n n!}
$$ 

Since the weights appear both in numerator and denominator, we can remove any factor independent of $j$ and use 

$$
w_j = (-1)^{j} \binom{N}{j}
$$ 

If the interval is $[a,b]$, we have to multiply the above $w_j$ by $2^N (b-a)^{-N}$, but this factor can also be dropped since it is common to numerator and denominator. The weights $w_j$ change by exponentially large factors, of order approximately $2^N$. This leads to ill-conditioning and Runge phenomenon. This is not a consequence of the barycentric foruma, but is intrinsic to interpolation using equispaced data.

### Chebyshev points of first kind

The points and weights are given by

$$
x_j = \cos\left( \frac{2j+1}{2N+2} \pi \right), \quad w_j = (-1)^j \sin\left( \frac{2j+1}
{2N+2} \pi \right), \quad j=0,1,\ldots,N
$$ 

These weights vary by factors of $O(N)$, since

$$
\frac{\max w_j}{\min w_j} \approx \frac{ \sin\left( \frac{2(N/2)+1}{2N+2} \pi \right) }
{ \sin\left( \frac{2(0)+1}{2N+2} \pi \right)} \approx \frac{\sin(\half\pi)}{\frac{\pi}
{2N+2}} = O(N) \quad \textrm{for large $N$}
$$ 

not exponentially as in case of uniformly spaced points. The weights can be computed in $O(N)$ operations[^5].

### Chebyshev points of second kind

The points are given by

$$
x_i = \cos\left( \frac{N-i}{N} \pi \right), \qquad i=0,1,2,\ldots,N
$$
and the weights are

$$
w_0 = \half, \qquad w_i = (-1)^i, \quad i=1,2,\ldots,N-1, \qquad w_N = \half (-1)^N
$$

[^1]: Assignment

[^2]: Assignment. For proof, see problems collection

[^3]: See page 51 in Powell, Approximation theory and methods, 1981.

[^4]: Assignment, See Atkinson, page 186, problem 2

[^5]: Newton interpolation always requires $O(N^2)$ operations for the
    divided differences.
