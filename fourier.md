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

# Fourier transform

```{include} math.md
```

```{code-cell}
from pylab import *
```

## Fourier series

Let $f : [-\pi,\pi] \to \re$ be a function which is extended periodically to the whole real line. The Fourier series is defined as

$$
Sf(x) = a_0 + \sum_{k=1}^\infty [ a_k \cos(k x) + b_k \sin(k x) ]
$$

where

$$
a_0 &= \int_{-\pi}^\pi f(x) \ud x \\
a_k &= \int_{-\pi}^\pi f(x) \cos(k x) \ud x, \qquad k=1, 2, \ldots \\
b_k &= \int_{-\pi}^\pi f(x) \sin(k x) \ud x, \qquad k=1, 2, \ldots
$$


:::{prf:example}
Discontinuous function 

$$
u(x) = \begin{cases}
1, & \half\pi < x \le \frac{3}{2}\pi \\
0, & 0 < x \le \half \pi, \quad \frac{3}{2}\pi < x < 2\pi
\end{cases}
$$ 

$u$ is of bounded variation in $[0,2\pi]$. Its Fourier
coefficients are 

$$
\hatu_k = \begin{cases}
\pi, & k = 0 \\
0, & k \textrm{ even} \\
\frac{(-1)^{(k-1)/2}}{k}, & k \textrm{ odd}
\end{cases}
$$ 

$P_N u$ converges pointwise at a linear rate everywhere
except at points of discontinuity.
:::

:::{prf:example}
Infinitely differentiable but not periodic 

$$
u(x) = \sin(x/2)
$$ 

is
infinitely differentiable in $[0,2\pi]$, $u$ is periodic but
$u'(0^+) \ne u'(2\pi^-)$. Its Fourier coefficients are

$$
\hatu_k = \frac{2}{\pi(1 - 4 k^2)}
$$ 

Pointwise quadratic convergence
except at end points where it converges linearly.
:::

:::{prf:example}

$$
u(x) = \frac{3}{5 - 4\cos x}
$$ 

Infinitely differentiable and periodic
in $[0,2\pi]$. Its Fourier coefficients are

$$
\hatu_k = \frac{1}{2^{|k|}}
$$
:::

:::{prf:example}

$$
u(x) = \begin{cases}
0, & x \le \half\pi, x \ge \frac{3}{2}\pi \\
2x/\pi - 1, & \half\pi \le x \le \pi \\
3 - 2x/\pi, & \pi \le x \le \frac{3}{2}\pi
\end{cases}
$$ 

$$
|\hatu_k| = \order{ \frac{1}{|k|^2} }
$$
:::
