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

# Miscellaneous methods

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

The secant method approximates the function by an affine polynomial and finds the root of the polynomial. We can improve this method by constructing a quadratic polynomial approximation.

## Muller's method

Given three initial guesses $x_0, x_1, x_2$, construct a quadratic polynomial $p(x)$ such that 

$$
p(x_i) = f(x_i), \qquad i=0,1,2
$$ 

The next approximation to the root $x_3$ is taken as the root of $p$

$$
p(x_3) = 0
$$ 

We can write the Newton form of the interpolating polynomial[^1]

$$
p(x) = f(x_2) + (x-x_2) f[x_2,x_1] + (x-x_2)(x-x_1) f[x_2,x_1,x_0]
$$

or equivalently 

$$
p(x) = f(x_2) + w(x-x_2) + f[x_2,x_1,x_0](x-x_2)^2
$$

where

$$
w = f[x_2,x_1] + (x_2 - x_1) f[x_2,x_1,x_0] = f[x_2,x_1] + f[x_2,x_0] - f[x_0,x_1]
$$

Now find the root $x_3$ that is closest to $x_2$

$$
x_3 - x_2 = \frac{-w \pm \sqrt{w^2 - 4f(x_2)f[x_2,x_1,x_0]}}{2f[x_2,x_1,x_0]}
$$

We have to choose the sign in the numerator that leads to small numerator. However this can lead to loss of significant figures. We can rewrite this as

$$
x_3 = x_2 -  \frac{2 f(x_2)}{w \pm \sqrt{w^2 - 4 f(x_2) f[x_2,x_1,x_0]}}
$$

Now we choose sign in denominator so that denominator is large in magnitude. We apply this algorithm recursively to obtain a sequence of iterates $\{ x_n \}$. If $x_n \to \alpha$ and $f'(\alpha) \ne 0$, then $\alpha$ is a root of $f$. To prove this, show that 

$$
w \to f'(\alpha)
$$ 

so that in the limit the iteration looks like

$$
\alpha = \alpha - \frac{2 f(\alpha)}{f'(\alpha) \pm \sqrt{[f'(\alpha)]^2 - 4 f(\alpha) f"(\alpha)}}
$$ 

But $f'(\alpha) \ne 0$ and so the denominator is not zero, which implies that $f(\alpha) = 0$.

The order of convergence $p \approx 1.84$ and is a root of

$$
x^3 - x^2 - x -1 = 0
$$ 

So this method converges faster than secant method but slower than Newton-Raphson method. Moreover

$$
\lim_{n \to \infty} \frac{|\alpha - x_{n+1}|}{|\alpha - x_n|^p} = \left| \frac{f^{(3)} (\alpha)}{6f'(\alpha)} \right|^{(p-1)/2}
$$ 

Muller's method can give rise to complex roots even when starting with real initial guesses.

## Inverse quadratic interpolation

We treat $x$ as a function of $f$, i.e., $x = x(f)$. Like in Muller's method, let us take $x_0, x_1,x_2$ and $f_i = f(x_i)$, $i=0,1,2$, but now construct an approximation 

$$
x = q(f; x_0,x_1,x_2)
$$ 

such that

$$
x_i = q(f_i; x_0,x_1,x_2), \qquad i=0,1,2
$$ 

Once the quadratic polynomial $q$ is found, we estimate the root as

$$
x_4 = q(0; x_0,x_1,x_2)
$$ 

This can be used to generate the iteration

$$
x_{n+1} = q(0; x_n, x_{n-1}, x_{n-2})
$$ 

Note that all the iterates will be real if we start with real initial guess.

## Linear fractional method

Suppose $f(x)$ has a vertical and horizontal asymptote, e.g., $f(x) = \frac{1}{x} - a$. Then the secant method may not converge, see figure (xxx), where $x_3$ can go to $-\infty$. Instead of linear approximation, we can use a fractional approximation.

Given $x_0, x_1, x_2$, fit a rational polynomial

$$
g(x) = \frac{x - a}{bx - c}
$$ 

using

$$
g(x_i) = f(x_i), \qquad i=0,1,2
$$ 

The solution is simplified if we shift the origin, $y = x - x_2$ and

$$
g(y) = \frac{y-a}{by-c}, \qquad g(y_i) = f(x_i), \quad i=0,1,2
$$ 

The next approximation is $x_3 = x_2 + a$.

[^1]: Problem: Show that this satisfies the interpolation conditions.
