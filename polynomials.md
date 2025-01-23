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

# Roots of polynomials

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
import sympy
```

## Existence of roots

A polynomial of degree $n$ is a function of the form

$$
p(x) = a_0 + a_1 x + a_2 x^2 + \ldots + a_n x^n, \qquad a_n \ne 0
$$

where the coefficients $a_j$ may be real or complex.

:::{prf:theorem} Fundamental Theorem of Algebra
Every non-constant polynomial has at least one root in $\complex$.
:::

If $r_1$ is a root of $p(x)$, then it can be written as $p(x) = (x-r_1) q(x)$ where $q(x)$ is a polynomial of degree $n-1$. Continuing this way, we see that $p(x)$ can be written as 

$$
p(x) = a_n (x-r_1) (x - r_2) \ldots (x - r_n)
$$

Some of these roots may coincide; if $r_j$ appears $m$ times, then it is said to have multiplicity $m$.

:::{prf:theorem}
A polynomial of degree $n$ has exactly $n$ roots in $\complex$, with each root being counted a number of times equal to its multiplicity.
:::

Note that we may not have any real roots even if all the coefficients are real, for example, $x^2+1=0$ has no real solutions.

## Roots from eigenvalues

Eigenvalues of a square matrix are the roots of a certain polynomial called the characteristic polynomial. Conversely, the roots of given polynomial are the eigenvalues of a certain matrix.

The polynomial

$$
p(x) = x^n + c_{n-1} x^{n-1} + \ldots + c_1 x + c_0
$$

can be written as

$$
p(x) = (-1)^n \det \begin{bmatrix}
-x &    &     &        &     &  -c_0 \\
 1 & -x &     &        &     &  -c_1 \\
   &  1 & -x  &        &     &  -c_2 \\
   &    &  1  & \ddots &     &  \vdots \\
   &    &     & \ddots & -x  &  -c_{n-2} \\
   &    &     &        &  1  &  -x-c_{n-1} \\
\end{bmatrix}
$$

i.e.,

$$
p(x) = (-1)^n \det(A - xI), \qquad
A = \begin{bmatrix}
 0 &    &     &        &     &  -c_0 \\
 1 &  0 &     &        &     &  -c_1 \\
   &  1 &  0  &        &     &  -c_2 \\
   &    &  1  & \ddots &     &  \vdots \\
   &    &     & \ddots &  0  &  -c_{n-2} \\
   &    &     &        &  1  &  -c_{n-1} \\
 \end{bmatrix}
$$

Then

$$
p(r) = 0 \limplies \det(A - rI) = 0
$$

and every root of $p(x)$ is an eigenvalue of the $n \times n$ matrix $A$, called the **companion matrix**. We have taken the coefficient of $x^2$ to be one but if this is not the case, then the coefficients can be scaled and this does not change the roots.

+++

:::{prf:remark}
If $r$ is a root of $p(x)$, check that $u = (1,r,r^2,\ldots,r^{n-1})$ is a left eigenvector of $A$ with eigenvalue $r$, i.e., $u A = r u$.
:::

+++

The next function computes the roots by finding eigenvalues, we do not assume that $a_n = 1$, so we need to scale the coefficients.

```{code-cell}
# Input array a contains coefficient like this
#   p(x) = a[0] + a[1] * x + a[2] * x**2 + ... + a[n] * x**n
def roots1(a):
    n = len(a) - 1
    c = a[0:-1] / a[-1]
    A = zeros((n,n))
    A[:,-1] = -c
    A += diag(ones(n-1), -1)
    r = eigvals(A)
    return r
```

Make a polynomial of degree $n = m-1$ with random coefficients in $[-1,1]$ and find its roots.

```{code-cell}
m = 10
a = 2.0*rand(m) - 1.0
r = roots1(a)
print("roots = \n", r)
```

Check the roots by evaluating the polynomial at the roots. We use `polyval` which expects coefficients in opposite order to what we have assumed above.

```{code-cell}
print("|p(r)| = \n", abs(polyval(a[::-1],r)))
```

These values are close to machine precision which can be considered to be zero.

The function [numpy.roots](https://numpy.org/doc/stable/reference/generated/numpy.roots.html) computes roots by the same method.

```{code-cell}
r = roots(a[::-1])
print("roots = \n", r)
```

+++

:::{prf:example}
Consider the cubic polynomial

$$
x^3 - 21 x^2 + 120 x - 100 = 0
$$

whose exact roots are

$$
x_1 = 1, \qquad x_2 = x_3 = 10
$$

```{code-cell}
r = sympy.roots([1,-21,120,-100])
print('Roots = ',r)
for x in r:
    print(x)
```

Now let us perturb the highest coefficient by 1\%

$$
0.99 x^3 - 21 x^2 + 120 x - 100 = 0
$$

The roots are now

$$
x_1 \approx 1.000, \qquad x_2 \approx 11.17, \qquad x_3 \approx 9.041
$$

```{code-cell}
r = sympy.roots([sympy.Integer(99)/100,-21,120,-100])
for x in r:
    print(sympy.N(x))
```

While $x_1$ changes by a small amount, the other two roots change by about 10\% which is a large perturbation. Now change the coefficient to the other side by 1\%.

$$
1.01 x^3 - 21 x^2 + 120 x - 100 = 0
$$

The roots are now

$$
x_1 \approx 1.000, \qquad x_2, x_3 \approx 9.896 \pm 1.044 i
$$

```{code-cell}
r = sympy.roots([sympy.Integer(101)/100,-21,120,-100])
for x in r:
    print(sympy.N(x))
```

Again, the double roots are very sensitive to perturbation of the coefficient. 
:::

## Test for ill-conditioning

Consider a polynomial

$$
f(x) = a_n x^n + a_{n-1} x^{n-1} + \ldots + a_1 x + a_0
$$

Let $\alpha$ be a single root, i.e., $f(\alpha) = 0$ and $f'(\alpha) \ne 0$. Let us perturb the polynomial into a different polynomial

$$
p(x) = f(x) + \epsilon g(x)
$$

where

$$
g(x) = b_n x^n + b_{n-1} x^{n-1} + \ldots + b_1 x + b_0
$$

and let $\hat{\alpha} = \alpha + \delta$ be a root of $p(x)$. Then, assuming the
perturbation of the root is small, it can be estimated as

$$
|\delta| \approx \left| \frac{\epsilon g(\alpha)}{f'(\alpha)} \right|
$$

If $\alpha$ is a double root, then we can estimate the perturbation (Wilkinson)

$$
|\delta| \approx \left[ - \frac{2\epsilon g(\alpha)}{f''(\alpha)} \right]^\half
$$

+++

:::{prf:example}
The polynomial

$$
f(x) = x^3 - 21 x^2 + 120 x - 100
$$

has roots $x_1 = 1$, $x_2 = x_3 = 10$. Let us perturb the highest power term, i.e., take $g(x) = x^3$. Consider the single root $\alpha = x_1 = 1$, $\epsilon = -0.01$. Then

$$
|\delta| \approx \left| \frac{(-0.01)(1)^3}{81} \right| \approx 0.0001
$$

Hence the root $x_1$ is well-conditioned wrt perturbations in the coefficient $x^3$. Now consider the double root $\alpha = x_2 = x_3 = 10$. The perturbation in this root can be estimated as

$$
|\delta| \approx \left| \frac{(-2)(-0.01)10^3}{18}\right| \approx 1.054
$$

This seems to indicate that the double root is very sensitive to perturbation in the coefficient of $x^3$, which we have observed in the previous example.
:::

+++

## Wilkinson's polynomial

$$
p(x) = \prod_{i=1}^{20} (x - i), \qquad x \in [1,20]
$$

Find polynomial coefficients using sympy

```{code-cell} ipython3
x, p = sympy.symbols('x p')
p = 1
for i in range(1,21):
    p = p * (x - i)
a = sympy.Poly(p, x)
c = a.coeffs()
for i,coef in enumerate(c):
    print("c[%2d] = %d" % (i,coef))
```

```{code-cell} ipython3
print('Monomial form = ', a)
```

The coefficients are returned in this order

$$
p(x) = c[0] x^{20} + c[1] x^{19} + \ldots + c[19] x + c[20]
$$

i.e., `c[0]` is the coefficient of the largest degree term.

Implement the polynomial

```{code-cell} ipython3
# As a product of factors
def wpoly(x):
    p = 1.0
    for i in range(1,21):
        p = p * (x - i)
    return p

# Computing it as a monomial
def mpoly(c,x):
    p = 0.0
    for i,a in enumerate(c):
        p += a * x**(20-i)
    return p
```

Plot the polynomial in $[1,20]$ by sampling it on a uniform grid

```{code-cell} ipython3
xp = linspace(1,20,1000)
yp = polyval(c,xp)
plot(xp,yp)
plot(arange(1,21),zeros(20),'o',label='Roots')
title("Monomial form")
legend(), grid(True), ylim(-1e13,1e13);
```

```{code-cell} ipython3
xp = linspace(1,20,1000)
yp = wpoly(xp)
plot(xp,yp)
plot(arange(1,21),zeros(20),'o',label='Roots')
title("Factored form")
legend(), grid(True), ylim(-1e13,1e13);
```

Computing the polynomial as a monomial is subject to lot of rounding errors since the coefficients are very large.

Find the roots using `numpy.roots` which computes it using the eigenvalue approach

```{code-cell} ipython3
r = roots(c)
print(r)
```

```{code-cell}
rexact = arange(20,0,-1)
semilogy(rexact, abs(r-rexact)/rexact, 'o')
xticks(rexact), grid(True), xlabel('Exact root'), ylabel('Relative error');
```

Randomly perturb the monomial coefficients and find roots

$$
c_0(1 + \epsilon r_0) x^{20} + c_1 (1+ \epsilon r_1) x^{19} + \ldots + c_{20}(1+ \epsilon r_{20})
$$

where $\epsilon = 10^{-10}$ and each $r_j$ is an independent Gaussian random variable with mean zero and standard deviation 1.

```{code-cell} ipython3
eps   = 1e-10 # Relative perturbation
nsamp = 100   # 100 random perturbations
for i in range(nsamp):
    r = normal(0.0, 1.0, 21)
    cc = c * (1 + eps * r)
    root = roots(cc)
    plot(real(root), imag(root),'k.')
plot(arange(1,21),zeros(20), 'ro')
xlabel('Real'); ylabel('Imag');
```

The relative perturbation in coefficients is $O(10^{-10})$ while the change in roots is $O(1)$, which indicates a very high sensitivity wrt the coefficients.
