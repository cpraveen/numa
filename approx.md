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

# Approximation of functions

```{include} math.md
```

## Weirstrass Theorem

Weirstrass Theorem says that we can approximate continuous functions by polynomials to any desired accuracy. Take $f : [0,1] \to \re$ a bounded function and define

$$
p_n(x) = \sum_{k=0}^n \binom{n}{k} f(k/n) x^k (1-x)^{n-k}, \qquad x \in [0,1]
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
B_r^n(x) = \binom{n}{r} x^r (1-x)^{n-r}
$$

:::{prf:theorem}
Let $f:[a,b] \to \re$ be continuous and let $\epsilon > 0$. Then there is a polynomial $p$ for which

$$
|f(x) - p(x)| \le \epsilon, \qquad x \in [a,b]
$$
:::

This theorem states that given any error tolerance $\epsilon$, we can find a polynomial $p$ such that

$$
\max_{x\in[a,b]}|f(x) - p(x)| \le \epsilon
$$ 

i.e., polynomials can provide *uniform* approximation. This is only an existence result and does not give us a solution with a desired error bound.

:::{prf:example}

Consider 

$$
f(x) = x^2, \qquad x \in [0,1]
$$ 

Then the Bernstein polynomial satisfies 

$$
p_n(x) - f(x) = \frac{1}{n} x (1-x)
$$ 

The error decreases at linear rate, which is very slow convergence. A quadratic interpolation will recovery the function exactly using just three function values.

$$
1 = 1^n = (1-x + x)^n = \sum_{k=0}^n \binom{n}{k} x^k (1-x)^{n-k}
$$

$$
p_n(x) - x^2 = \sum_{k=0}^n \binom{n}{k} ( (k/n)^2 - x^2) x^k (1-x)^{n-k}
$$
:::

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

Plotting the error, we see that it is not uniformly distributed in the interval $[-1,1]$. We have very good approximation near $x=0$ but poor approximation near the end-points of the interval.

FIGURE

## Cubic interpolation

We can get much better cubic polynomial approximation than the Taylor polynomial. Let us interpolate a cubic polynomial at $\{-1, -\frac{1}{3}, +\frac{1}{3}, +1\}$ which yields

$$
p_3(x) = 0.99519577 + 0.99904923 x + 0.54788486 x^2 + 0.17615196 x^3
$$

The maximum error is approximately $0.00998$ which is almost five times smaller than the Taylor polynomial. The error is also more uniformly distributed and oscillates in sign.

FIGURE
