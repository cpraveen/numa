---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Special techniques

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## Change of variable

Many times, a change of variable can improve the regularity of the
integrand.

:::{prf:example} Singular integrand

$$
I = \int_0^b \frac{f(x)}{\sqrt{x}} \ud x, \qquad \textrm{$f$ is smooth}
$$

With 

\begin{gather}
x = u^2, \qquad 0 \le u \le \sqrt{b} \\
\Downarrow \\
I = 2 \int_0^{\sqrt{b}} f(u^2) \ud u
\end{gather}

:::

+++

:::{prf:example} Singular derivative

$$
\int_0^1 \sin(x) \sqrt{1-x^2} \ud x = 2\int_0^1 u^2 \sqrt{2-u^2} \sin(1-u^2) \ud u
$$

by change of variables $u = \sqrt{1-x}$. The first form has singular
derivative at $x=1$ while the second form has infinitely smooth
integrand.
:::

+++

:::{prf:example} Infinite interval of integration

$$
I = \int_1^\infty \frac{f(x)}{x^p} \ud x, \qquad p > 1
$$ 

Assume $\lim_{x \to \infty} f(x)$ exists and $f(x)$ is smooth on $[1,\infty)$. Change of variable

\begin{gather}
x = \frac{1}{u^\alpha}, \qquad \ud x = -\frac{\alpha}{u^{1+\alpha}} \ud u, \quad \alpha > 0 \\
\Downarrow \\
I = \alpha \int_0^1 u^{(p-1)\alpha-1} f(1/u^\alpha) \ud u
\end{gather}

Integrand can be singular at $u=0$, so pick $\alpha$ so that the exponent $(p-1)\alpha-1$ is large $(\ge 1)$. 

As a specific example consider

\begin{gather}
I = \int_1^\infty \frac{f(x)}{x\sqrt{x}} \ud x \\
x = \frac{1}{u^4} \quad\implies\quad I = 4 \int_0^1 u f(1/u^4) \ud u
\end{gather}

For $x \to \infty$, $f(x)$ behaves like

\begin{gather}
f(x) = c_0 + \frac{c_1}{x} + \frac{c_2}{x^2} + \ldots \\
\Downarrow \\
uf(1/u^4) = c_0 u + c_1 u^5 + c_2 u^9 + \ldots
\end{gather}

and the integrand behaves nicely as $u \to 0$.
:::

+++

:::{prf:example} Chebyshev quadrature

$$
I = \int_{-1}^1 \frac{f(x)}{\sqrt{1-x^2}} \ud x
$$ 

For this weight function, the orthogonal polynomials are Chebyshev polynomials. The integration nodes are roots of $T_n(x) = \cos(n \cos^{-1}x)$ (Chebyshev points of I kind)

$$
x_j = \cos\left( \frac{2j-1}{2n}\pi \right), \qquad j=1,2,\ldots,n
$$

and weights 

$$
w_j = \frac{\pi}{n}
$$ 

so that

$$
I \approx I_n = \frac{\pi}{n} \sum_{j=1}^n f(x_j)
$$ 

This corresponds to doing the change of variable $x = \cos\theta$ and applying the mid-point rule in the $\theta$ variable

$$
I = \int_0^\pi f(\cos\theta) \ud\theta \approx \frac{\pi}{n} \sum_{j=1}^n f(\cos\theta_j), \qquad \theta_j = \frac{(2j-1)\pi}{2n}
$$

:::

## Analytic treatment of singularity

Consider an integral with a singular term

$$
I = \int_0^b f(x) \log(x) \ud x
$$ 

Since $\log(x)$ is singular at $x=0$, we have to compute this with some care. Divide the interval into two parts

$$
I  = \int_0^\epsilon f(x) \log(x) \ud x + \int_{\epsilon}^b f(x) \log(x) \ud x =: I_1 + I_2
$$

$I_2$ can be computed by some standard technique and $I_1$ is computed semi-analytically. Assume $f(x)$ has convergent Taylor series on $[0,\epsilon]$

$$
f(x) = \sum_{j=0}^\infty a_j x^j, \qquad x \in [0,\epsilon]
$$ 

Then

$$
I_1 = \int_0^\epsilon \left( \sum_j a_j x^j \right) \log(x) \ud x = \sum_{j=0}^\infty \frac{a_j \epsilon^{j+1}}{j+1} \left[ \log(\epsilon) - \frac{1}{j+1} \right]
$$

+++

:::{prf:example} Singular integrand

$$
I = \int_0^{4\pi} \cos(x) \log(x) \ud x
$$ 

Then

$$
I_1 
&= \int_0^{\epsilon} \cos(x) \log(x) \ud x \\
&= \epsilon[\log(\epsilon)-1] - \frac{\epsilon^3}{6}[\log(\epsilon) - 1/3] + \frac{\epsilon^5}{600}[\log(\epsilon) - 1/5] + \ldots
$$

The terms in the above series become smaller, so we can truncate this
series at some term.
:::

## Product integration

Consider an integral with a singular term

$$
I(f) = \int_a^b w(x) f(x) \ud x
$$ 

where $w(x)$ is singular and $f(x)$ is a nice function. Let $\{ f_n \}$ be a sequence of approximations to $f$ such that

1.  $|w|$ is integrable.

2.  $\lim_{n\to\infty} \norm{f - f_n}_\infty = 0$

3.  The integrals 

    $$
    I_n(f) = \int_a^b w(x) f_n(x) \ud x
    $$ 

    can be analytically computed.

Then

$$
|I(f) - I_n(f)| \le \norm{f-f_n}_\infty \int_a^b |w(x)| \ud x \to 0
$$

+++

:::{prf:example} Singular integrand

$$
I(f) = \int_0^b f(x) \log(x) \ud x
$$ 

Let $f_n$ be piecewise linear approximation on a uniform partition

$$
x_j = jh, \quad j=0,1,\ldots,n \qquad\textrm{where}\qquad h = \frac{b}{n}
$$

and

$$
f_n(x) = \frac{1}{h}[ (x_j - x) f_{j-1} + (x-x_{j-1}) f_j], \qquad x_{j-1} \le x \le x_j
$$

for which we know the error estimate from [](#sec:pwlininterp)

$$
\norm{f-f_n}_\infty \le \frac{h^2}{8} \norm{f''}_\infty
$$ 

Then, with $I_n(f) = I(f_n)$

$$
|I(f) - I_n(f)| \le \frac{h^2}{8} \norm{f''}_\infty \int_0^b |\log(x)| \ud x \to 0
$$

$I_n(f)$ can be computed analytically. In particular, in the first sub-interval $[x_0,x_1]=[0,h]$ 

$$
f_n(x) = \frac{h-x}{h} f_0 + \frac{x}{h} f_1
$$

we have

$$
\int_{x_0}^{x_1} \log(x) f_n(x) \ud x = \frac{f_0}{h} \int_0^h (h-x)\log(x) \ud x + \frac{f_1}{h} \int_0^h x \log(x) \ud x
$$

and these are finite since

$$
\int_0^h \log(x) \ud x = h[\log(h)-1], \qquad \int_0^h x\log(x)\ud x = \frac{h^2}{4}[ 2\log(h) - 1 ]
$$

The integrals in other intervals can also be computed analytically.
:::
