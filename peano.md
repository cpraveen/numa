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

# Peano kernel error formulas

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

The error formulae we have derived so far assume certain derivatives exist and are continuous. This may be a restrictive assumption in many cases. Even if the function is not differentiable everywhere, the derivatives may be integrable. In this case, we can use Taylor formula with an integral remainder term to perform the error analysis.

## Trapezoidal rule

Assume that $f$ is continuously differentiable and $f''$ is integrable on $[a,b]$, i.e., 

$$
\int_a^b |f''(x)| \ud x < \infty
$$ 

**Single  interval.** Taylor's theorem gives

\begin{align}
f(x) &= \clr{red}{f(a) + (x-a) f'(a)} + \clr{blue}{ \int_a^x (x-t) f''(t) \ud t } \\
&=: \clr{red}{p_1(x)} + \clr{blue}{ R_2(x) }
\end{align}

The quadrature error $E_n$ is a linear operator

$$
E_n(F+G) = E_n(F) + E_n(G), \qquad F,G \in \cts[a,b]
$$ 

Thus

$$
E_1(f) = E_1(p_1) + E_1(R_2) = E_1(R_2)
$$ 

since Trapezoidal rule is exact for linear polynomials. Now 

$$
\begin{aligned}
E_1(R_2) 
&= \int_a^b R_2(x) \ud x - \frac{b-a}{2}[R_2(a) + R_2(b)], \qquad R_2(a) = 0 \\
&= \int_a^b \int_a^x (x-t) f''(t) \ud t \ud x - \frac{b-a}{2}\int_a^b (b-t) f''(t) \ud t
\end{aligned}
$$ 

For any integrable function $G(x,t)$, by Fubini theorem

\begin{align}
\int_a^b \int_a^x G(x,t) \ud t \ud x 
&= \int_a^b \left( \int_a^b G(x,t) \chi_{[a,x]}(t) \ud t \right) \ud x \\
&= \int_a^b \left( \int_a^b G(x,t) \chi_{[a,x]}(t) \ud x \right) \ud t \\
&= \int_a^b \left( \int_a^b G(x,t) \chi_{[t,b]}(x) \ud x \right) \ud t \\
&= \int_a^b \int_t^b G(x,t) \ud x \ud t
\end{align}

Here $\chi$ is the characteristic function and we used $\chi_{[a,x]}(t) = \chi_{[t,b]}(x)$ since $a \le t \le x \le b$. Using this result in the above error formula 

$$
\begin{aligned}
E_1(R_2) 
&= \int_a^b f''(t) \int_t^b (x-t) \ud x \ud t - \frac{b-a}{2} \int_a^b (b-t) f''(t) \ud t \\
&= \int_a^b \frac{(t-a)(t-b)}{2} f''(t) \ud t
\end{aligned}
$$ 

**Composite rule.** For the composite Trapezoid rule, the error is obtained by adding the error from each interval and can be written as

$$
E_n(f) = \int_a^b K(t) f''(t) \ud t
$$ 

where

$$
K(t) = \half (t - x_{j-1})(t - x_j), \qquad t \in [x_{j-1}, x_j], \qquad j=1,2,\ldots,n
$$

This gives the following error estimate

$$
|E_n(f)| \le \norm{K}_\infty \int_a^b |f''(t)| \ud t = \frac{h^2}{8} \int_a^b |f''(t)| \ud t
$$

Compare this to the previous error estimate

$$
|E_n(f)| = \frac{b-a}{12} h^2 |f''(\eta)| \le \frac{b-a}{12} h^2 \norm{f''}_\infty
$$

which may over-estimate the error if $f''$ has a peaky distribution, since

$$
\frac{1}{b-a} \int_a^b |f''(t)|\ud t \le \norm{f''}_\infty
$$

+++

:::{prf:example}

The function 

$$
f(x) = x^{3/2}, \qquad x \in [0,1]
$$ 

does not have finite second derivative since 

$$
f''(x) = \frac{3}{4\sqrt{x}}, \qquad \norm{f''}_\infty = \infty
$$ 

The error formula based on derivative cannot be used but since

$$
\int_0^1 |f''(x)|\ud x = \frac{3}{2}
$$ 

the error formula based on integral remainder term can be used to conclude that

$$
|E_n(f)| \le \frac{3 h^2}{16} \to 0 \qquad \textrm{as $n \to \infty$}
$$

```{code-cell}
# Performs Trapezoid quadrature
def integrate(a,b,n,f):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res1 = 0.5*y[0] + sum(y[1:n]) + 0.5*y[n]
    return h*res1
```

```{code-cell}
f = lambda x: x * sqrt(x)
a, b = 0.0, 1.0
Ie = 2.0/5.0 # Exact integral
n, N = 2, 10

e = zeros(N)
for i in range(N):
    h = (b-a)/n
    I = integrate(a,b,n,f)
    e[i] = Ie - I
    if i > 0:
        print('%6d %18.8e %14.5g %18.8e'% (n,e[i],e[i-1]/e[i],3*h**2/16))
    else:
        print('%6d %18.8e %14.5g %18.8e'%(n,e[i],0,3*h**2/16))
    n = 2*n
```

We still get $O(h^2)$ convergence and the error shown in second column is smaller than the error estimate shown on last column.
:::
