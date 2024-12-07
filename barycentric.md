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

# Barycentric form

```{code-cell}
from pylab import *
```

The Lagrange interpolation polynomial of degree $N$ is

$$
p(x) = \sum_{j=0}^N f_j \ell_j(x) \qquad \textrm{where} \qquad
\ell_j(x) = \frac{\prod\limits_{i=0,i\ne j}^N (x-x_i)}{\prod\limits_{i=0,i \ne j}^N(x_j - x_i)}
$$ 

The denominators in $\ell_j$ can be pre-computed and stored.  Then for every $x$, the evaluation of $p(x)$ requires $O(N^2)$ flops.  Ideally, we desire a method which takes $O(N)$ flops. Also, if we add a new data pair $(x_{N+1},f_{N+1})$, then this requires a new computation from scratch. Moreover, the numerical evaluation is also unstable to round-off errors (See page 51 in Powell, Approximation theory and methods, 1981.)

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

Since $\ell_j$ form a partition of unity (Assignment, See Atkinson, page 186, problem 2)

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

not exponentially as in case of uniformly spaced points. The weights can be computed in $O(N)$ operations. Newton interpolation always requires $O(N^2)$ operations for the divided differences.

### Chebyshev points of second kind

The points are given by

$$
x_i = \cos\left( \frac{N-i}{N} \pi \right), \qquad i=0,1,2,\ldots,N
$$
and the weights are

$$
w_0 = \half, \qquad w_i = (-1)^i, \quad i=1,2,\ldots,N-1, \qquad w_N = \half (-1)^N
$$

:::{prf:example}
Let us interpolate the Runge function on $[-1,+1]$

$$
f(x) = \frac{1}{1 + 16 x^2}
$$

using the Barycentric formula

$$
p(x) = \frac{ \sum_{i=0}^N \frac{w_i}{x - x_i} f_i }{ \sum_{i=0}^N \frac{w_i}{x - x_i} }
$$

where the weight $w_i$ is given by

$$
w_i = \frac{1}{\prod_{j=0,j\ne i}^N (x_i - x_j)}
$$

This is the function we want to interpolate.

```{code-cell}
def fun(x):
    f = 1.0/(1.0+16.0*x**2)
    return f
```

We next write a function that constructs and evaluates the Lagrange interpolation.


```{code-cell}
# X[nX], Y[nX] : Data points
# x[nx]        : points where we want to evaluate
def LagrangeInterpolation(X,Y,x):
    nx = size(x)
    nX = size(X)
    # compute the weights
    w = ones(nX)
    for i in range(nX):
        for j in range(nX):
            if i != j:
                w[i] = w[i]/(X[i]-X[j])
    # Evaluate the polynomial at x
    num= zeros(nx)
    den= zeros(nx)
    eps=1.0e-14
    for i in range(nX):
        num = num + Y[i]*w[i]/((x-X[i])+eps)
        den = den + w[i]/((x-X[i])+eps)
    f = num/den
    return f
```

```{code-cell}
xmin, xmax = -1.0, +1.0
N = 15 # Degree of polynomial
```

We first interpolate on uniformly spaced points.

```{code-cell}
X = linspace(xmin,xmax,N+1)
Y = fun(X)
x = linspace(xmin,xmax,100)
fi = LagrangeInterpolation(X,Y,x)
fe = fun(x)
plot(x,fe,'b--',x,fi,'r-',X,Y,'o')
title('Degree '+str(N)+' using uniform points')
legend(("True function","Interpolation","Data"),loc='lower center');
```

Next, we interpolate on Chebyshev points.


```{code-cell}
X = cos(linspace(0.0,pi,N+1))
Y = fun(X)
x = linspace(xmin,xmax,100)
fi = LagrangeInterpolation(X,Y,x)
fe = fun(x)
plot(x,fe,'b--',x,fi,'r-',X,Y,'o')
title('Degree '+str(N)+' using Chebyshev points')
legend(("True function","Interpolation","Data"),loc='upper right');
```

:::

:::{prf:example} Barycentric Lagrange interpolation on Chebyshev points
Let us interpolate the Runge function on $[-1,+1]$

$$
f(x) = \frac{1}{1 + 16 x^2}
$$

using the Barycentric formula

$$
p(x) =  \frac{ \sum\limits_{i=0}^N{}' \frac{(-1)^i}{x - x_i} f_i }{ \sum\limits_{i=0}^N{}' \frac{(-1)^i}{x - x_i} }
$$

where the prime on the summation means that the first and last terms must be multiplied by a factor of half.

Define the function to be interpolated.

```{code-cell}
def fun(x):
    f = 1.0/(1.0+16.0*x**2)
    return f
```

The next function evaluates the Lagrange interpolation using Chebyshev points.

```{code-cell}
def BaryInterp(X,Y,x):
    nx = size(x)
    nX = size(X)
    f  = 0*x
    # Compute weights
    w  = (-1.0)**arange(0,nX)
    w[0]    = 0.5*w[0]
    w[nX-1] = 0.5*w[nX-1]
    # Evaluate barycentric foruma at x values
    for i in range(nx):
        num, den = 0.0, 0.0
        for j in range(nX):
            if abs(x[i]-X[j]) < 1.0e-15:
                num = Y[j]
                den = 1.0
                break
            else:
                num += Y[j]*w[j]/((x[i]-X[j]))
                den += w[j]/(x[i]-X[j])
        f[i] = num/den
    return f
```

```{code-cell}
xmin, xmax = -1.0, +1.0
N = 19 # degree of polynomial
```

Let us interpolate on Chebyshev points.

```{code-cell}
X = cos(linspace(0.0,pi,N+1))
Y = fun(X)
x = linspace(xmin,xmax,100)
fi = BaryInterp(X,Y,x)
fe = fun(x)
figure(figsize=(8,5))
plot(x,fe,'b--',x,fi,'r-',X,Y,'o')
legend(("True function","Interpolation","Data"),loc='upper right')
axis([-1.0,+1.0,0.0,1.1]);
```

:::
