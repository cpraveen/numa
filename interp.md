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

# Interpolation

```{include} math.md
```

```{code-cell}
from pylab import *
```

Suppose we know the values of a function $y=f(x)$ at a set of discrete points, say in a tabular form

:::{table}
:align: center
|   $x$      |  $f(x)$   |
|   :---:    |   :---:   |
|   $x_0$    |  $f(x_0)$ |
|   $x_1$    |  $f(x_1)$ |
|   $\vdots$ |  $\vdots$ |
|   $x_N$    |  $f(x_N)$ |
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
