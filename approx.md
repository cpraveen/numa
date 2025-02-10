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

# Approximation

```{include} math.md
```

```{code-cell}
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

## Weirstrass Theorem

Weirstrass Theorem says that we can approximate continuous functions by polynomials. Take $f : [0,1] \to \re$ a bounded function and define

$$
p_n(x) = \sum_{k=0}^n \binom{n}{k} x^k (1-x)^{n-k} f(k/n), \qquad x \in [0,1]
$$

If $f$ is continuous at $x$ then 

$$
\lim_{n \to \infty} p_n(x) = f(x)
$$

If $f$ is continuous at all $x \in [0,1]$ then $p_n \to f$ uniformly in $[0,1]$, i.e.,

$$
\lim_{n \to \infty} \max_{0 \le x \le 1}|f(x) - p_n(x)| = 0
$$ 

We have written $p_n$ in terms of Bernstein polynomials

$$
B_r^n(x) = \binom{n}{r} x^r (1-x)^{n-r} \qquad \textrm{and} \qquad \binom{n}{k} = \frac{n!}{k! (n-k)!}
$$

+++

:::{prf:theorem}
Let $f:[a,b] \to \re$ be continuous. For any $\epsilon > 0$, there is a polynomial $p$ for which

$$
|f(x) - p(x)| \le \epsilon, \qquad x \in [a,b]
$$
:::

This theorem states that given any error tolerance $\epsilon > 0$, there exists a polynomial $p$ such that

$$
\max_{x\in[a,b]}|f(x) - p(x)| \le \epsilon
$$ 

i.e., polynomials can provide *uniform* approximation. This is only an existence result and does not give us a polynomial with a desired error bound. It also does not even tell us the degree of the polynomial.

:::{prf:example}

Consider 

$$
f(x) = x^2, \qquad x \in [0,1]
$$ 

Then the Bernstein polynomial satisfies 

$$
\lim_{n \to \infty} n [p_n(x) - x^2] = x (1-x) \limplies p_n(x) - x^2 \approx \frac{1}{n} x (1-x), \quad n \gg 1
$$ 

The error decreases at linear rate, which is very slow convergence. A quadratic interpolation will recovery the function exactly using just three function values.

% TODO
%$$
%1 = 1^n = (1-x + x)^n = \sum_{k=0}^n \binom{n}{k} x^k (1-x)^{n-k}
%$$
%
%$$
%p_n(x) - x^2 = \sum_{k=0}^n \binom{n}{k} ( (k/n)^2 - x^2) x^k (1-x)^{n-k}
%$$
:::

+++

## Taylor expansion

The Taylor expansion provides another way to approximate a function in a neighbourhood of some point at which the expansion is written. A Taylor expansion will exactly recover a polynomial function so we consider somewhat more non-trivial but still simple function. E.g.

$$
f(x) = \exp(x), \qquad x \in [-1,1]
$$ 

The third degree Taylor expansion around $x=0$ is

$$
p_3(x) = 1 + x + \half x^2 + \frac{1}{6} x^3
$$ 

We can estimate the error also from the remainder term

$$
\exp(x) - p_3(x) = \frac{1}{24} x^4 \exp(\xi), \qquad \textrm{$\xi$ between $0$ and $x$}
$$ 

and we get the following estimate

$$
\frac{1}{24} x^4 \le \exp(x) - p_3(x) \le \frac{\ee}{24} x^4, \qquad 0 \le x \le 1
$$

$$
\frac{1}{24\ee} x^4 \le \exp(x) - p_3(x) \le \frac{1}{24} x^4, \qquad -1 \le x \le 0
$$

By a direct calculation, the error is

$$
\max_{-1 \le x \le 1}|\exp(x) - p_3(x)| \approx 0.0516
$$ 

Plotting the error, 

```{code-cell}
x = linspace(-1,1,100)
f = lambda x: exp(x)
p3 = lambda x: 1 + x + x**2/2 + x**3/6
e  = f(x) - p3(x)
print("Max error = ", abs(e).max())
plot(x,e)
title('Error in cubic Taylor expansion')
grid(True), xlabel('x'), ylabel('$f(x)-p_3(x)$');
```

we see that it is not uniformly distributed in the interval $[-1,1]$. We have very good approximation near $x=0$ but poor approximation near the end-points of the interval.

+++

## Cubic interpolation

We can get much better cubic polynomial approximation than the Taylor polynomial. Let us interpolate a cubic polynomial at $\{-1, -\frac{1}{3}, +\frac{1}{3}, +1\}$ which yields

$$
p_3(x) = 0.99519577 + 0.99904923 x + 0.54788486 x^2 + 0.17615196 x^3
$$

The error in this approximation is shown below.

```{code-cell}
xi = linspace(-1,1,4)
yi = f(xi)
y  = barycentric_interpolate(xi,yi,x)
e  = f(x) - y
print("Max error = ",abs(e).max())
plot(x,e, xi, 0*xi, 'o')
title('Error in cubic interpolation')
grid(True), xlabel('x'), ylabel('$f(x)-p_3(x)$');
```

The maximum error is approximately $0.00998$ which is almost five times smaller than the Taylor polynomial. The error is also more uniformly distributed and oscillates in sign.
