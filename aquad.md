---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Adaptive quadrature

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

The error in some numerical quadrature depends on the number of nodes used which the user has to select. It is desirable to have an automatic algorithm where the user specifies the error $\epsilon$ that is acceptable and the algorithm should select the number of nodes $n$ so that 

$$
|I(f) - I_n(f)| \le \epsilon
$$ 

The error clearly depends on the smoothness of the integrand $f$ also. If $f$ has a singularity, e.g., $f(x) = \sqrt{x}$ has singular derivatives at $x=0$, then the error will be large. To have an acceptable error we must use a large $n$ (small $h$). But this will increase the number of nodes. Actually, most of the error will come from the region near the singularity and we can use a small $h$ only in those regions, and not everywhere. This requires the ability to use non-uniformly spaced nodes in the numerical integration formula, and all the composite rules allow us to do this.

## Adaptive Simpson quadrature

Let $n$ be even and we use $n+1$ points. The composite Simpson quadrature is

$$
I_n(f) = \sum_{j=1}^{n/2} \left( \frac{x_{2j} - x_{2j-2}}{6} \right)( f_{2j-2} + 4 f_{2j-1} + f_{2j}), \qquad x_{2j-1} = \half(x_{2j-2} + x_{2j})
$$

with error formula

$$
I(f) - I_n(f) = -\frac{1}{2880} \sum_{j=1}^{n/2}(x_{2j} - x_{2j-2})^5 f^{(4)}(\xi_j), \qquad \xi_j \in [x_{2j-2}, x_{2j}]
$$

Note that we have a decomposition of the total error into error contributions from each of the sub-intervals. If $f^{(4)}$ is large in some sub-interval, we should choose the length of that sub-interval $x_{2j} - x_{2j-2}$ to be small, in order to reduce the error. 

> The best strategy is to equi-distribute the error among all the sub-intervals.  

Define 

$$
M_j = \max_{x \in [x_{2j-2},x_{2j}]} |f^{(4)}(x)|
$$ 

then

$$
|I(f) - I_n(f)| \le \frac{1}{2880} \sum_{j=1}^{n/2}(x_{2j} - x_{2j-2})^5 M_j
$$

We want 

$$
|I(f) - I_n(f)| \le \epsilon
$$ 

This can be achieved if

$$
\frac{1}{2880} (x_{2j} - x_{2j-2})^5 M_j \le \frac{\epsilon}{n/2} = \frac{2 \epsilon}{n}, \qquad j=1,2,\ldots,n/2
$$

If this condition is not satisfied, then we iteratively refine the sub-intervals until desired error level is reached.

### Algorithm I

We will assume that we can estimate $M_j$ and hence we can estimate the error. The algorithm will iteratively increase the number of sub-intervals until the desired error level is reached.

1.  Start with three nodes, i.e., one Simpson rule, and compute $I_2(f)$
    and estimate the error. If $|I(f) - I_2(f)| \le \epsilon$, then we
    stop and the answer is $I_2(f)$.

2.  Otherwise, insert a point mid-way between the three points and use
    two Simpson rules. Estimate the error from each rule. If both errors
    are $\le \epsilon/2$, then stop.

3.  If not, suppose first error $> \epsilon/2$ and second error
    $\le \epsilon/2$. Refine the first Simpson rule by inserting two
    points.

4.  Now the error from the first two rules must be $\le \epsilon/4$. If
    this is true, then we stop.

5. Otherwise, continue with sub-division.

### Algorithm II

We may not be able to estimate the derivatives, and the only information available is the function values. We can approximate the derivatives by first constructing a local approximation and differentiating it. But we have seen other ways to estimate the error in integral, e.g., using Richardson extrapolation.

Set $\alpha=a$ and $\beta=b$ and define

$$
I_{\alpha,\beta} = \int_\alpha^\beta f(x) \ud x
$$ 

The Simpson rule is

$$
I_{\alpha,\beta}^{(1)} = \frac{h}{3}\left[f(\alpha) + f\left(\frac{\alpha+\beta}{2}\right) + f(\beta) \right], \qquad h = \half (\beta - \alpha)
$$

To estimate the error, we compute a refined estimate as

$$
I_{\alpha,\beta}^{(2)} = I_{\alpha,\gamma}^{(1)} + I_{\gamma,\beta}^{(1)}, \qquad \gamma = \half(\alpha + \beta)
$$

and use error estimate based on Richarson extrapolation, see [](#sec:simprichexterr)

$$
E_{\alpha,\beta} = \frac{1}{15} |I_{\alpha,\beta}^{(2)} - I_{\alpha,\beta}^{(1)}| \approx |I_{\alpha,\beta} - I_{\alpha,\beta}^{(2)}|
$$ 

as the error estimate.

1.  If $E_{\alpha,\beta} \le \epsilon$
    then keep interval $[\alpha,\beta]$ and use $I_{\alpha,\beta}^{(2)}$
    as estimate of $I_{\alpha,\beta}$.

2.  Otherwise, divide $[\alpha,\beta]$ into $[\alpha,\gamma]$ and
    $[\gamma,\beta]$. Now we want to satisfy

    $$
    E_{\alpha,\gamma} \le \frac{\epsilon}{2}, \qquad E_{\gamma,\beta} \le \frac{\epsilon}{2}
    $$

3.  If both conditions are satisfied, then the integral is
    $I_{\alpha,\gamma}^{(2)} + I_{\gamma,\beta}^{(2)}$.

4.  Otherwise, suppose

    $$
    E_{\alpha,\gamma} > \frac{\epsilon}{2}, \qquad E_{\gamma,\beta} \le \frac{\epsilon}{2}
    $$

    Then we keep $I_{\gamma,\beta}^{(2)}$ and divide $[\alpha,\gamma]$
    into $[\alpha,\delta]$ and $[\delta,\gamma]$ with
    $\delta=\half(\alpha + \gamma)$. Now try to satisfy

    $$
    E_{\alpha,\delta} \le \frac{\epsilon}{4}, \qquad E_{\delta,\gamma} \le \frac{\epsilon}{4}
    $$

    If both conditions are satisfied, then integral is
    $I_{\alpha,\delta}^{(2)} + I_{\delta,\gamma}^{(2)} + I_{\gamma,\beta}^{(2)}$.

5.  Continue this process until the errors in all sub-intervals is less
    than the tolerance.


This algorithm is implemented in the following code.

* `simpson1` performs Simpson quadrature in one interval $[a,b]$.
* `asimpson` performs adaptive quadrature in a recursive manner.

```{code-cell}
from scipy.integrate import simpson

# Apply Simpson using points {a, (a+b)/2, c}
def simpson1(f,a,b):
    h = 0.5*(b-a)
    c = 0.5*(a+b)
    return (h/3.0)*(f(a) + 4*f(c) + f(b))

# Performs adaptive Simpson to make error < eps
def asimpson(f,a,b,eps):
    c = 0.5*(a+b)
    i1 = simpson1(f,a,b)
    i2 = simpson1(f,a,c) + simpson1(f,c,b)
    err = (1.0/15.0) * abs(i2 - i1)
    if err < eps:
        print('[%12.6g %12.6g] %12.6g %12.6g %12.6g' % (a,b,err,eps,i2))
        return i2
    else:
        eps1 = 0.5*eps
        return asimpson(f,a,c,eps1) + asimpson(f,c,b,eps1)
```

+++

:::{prf:example}
Let us test it for

$$
I = \int_0^1 \sqrt{x} \ud x = \frac{2}{3}
$$

Here is the result of adaptive quadrature.

```{code-cell}
f = lambda x: sqrt(x)
a, b = 0.0, 1.0
eps = 0.0001
I  = asimpson(f,a,b,eps)
Ie = 2.0/3.0
print('\nIntegral,err,eps = %14.6e %14.6e %14.6e'%(I,abs(I-Ie),eps))
```

The interval $[0, 0.5]$ is sub-divided into smaller ones which is where the integrand is singular; the other interval $[0.5, 1]$ is not divided. The final error is significantly less than the specified error level of `eps = 0.0001`.

If we use uniformly spaced points, then

```{code-cell}
# Simpson rule using uniform sample points
# Choose n to be even only
n, N = 2, 30
for i in range(N):
    x = linspace(a, b, n+1)
    y = f(x)
    I = simpson(y, x=x)
    print('%5d %16.6e %16.6e %16.6e' % (int(n/2),(b-a)/n,I,abs(I-Ie)))
    if abs(I-Ie) < eps:
        break
    n = n + 4
```

This requires more points (46) to reach the error tolerance, and the final error is still not as small as with adaptive quadrature which uses fewer points.
:::

:::{warning}
Error estimates $E_{\alpha,\beta}$ like the one used above are reliable only when $\beta-\alpha$ is sufficiently small. Starting with a single Simpson rule as we did above is perhaps not a good idea in general. 
:::
