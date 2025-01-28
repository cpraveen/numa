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

# Potential theory

```{include} math.md
```

If the function has a singularity within the stadium S, i.e., close to
the real line, then we have to do more analysis and the effect of point
distribution also has to be considered in the analysis. Define

$$
\gamma_n(x,t) = \left| \frac{\ell(t)}{\ell(x)} \right|^{1/(n+1)} =
\left( \frac{ \prod_{j=0}^n |t-x_j| }{ \prod_{j=0}^n |x- x_j| } \right)^{1/(n+1)}
$$ 

then

$$
\left| \frac{\ell(t)}{\ell(x)} \right| = [\gamma_n(x,t)]^{-n-1}
$$ 

If $\gamma_n(x,t) > 1$, then we get exponential convergence as $n \to \infty$. Define

$$
\alpha_n = \min_{x \in X, t \in \Gamma} \gamma_n(x,t), \qquad X = [-1,+1]
$$

> If $\alpha_n \ge \alpha > 1$ for all sufficiently large $n$, and if $f$ is analytic in the region bounded by $\Gamma$, then $p(x) \to f(x)$ at the rate $O(\alpha^{-n})$.

We can give a geometric interpretation of the condition $\alpha_n > 1$; the geometric mean distance of every $t \in \Gamma$ from $\{ x_j \}$ is strictly greater than the geometric mean distance of every $x \in X$ to $\{ x_j \}$.

+++

## Discrete potential

We want to study the behaviour of the function $\gamma_n(x,t)$. Taking logarithm

$$
\log \gamma_n(x,t) = \frac{1}{n+1} \sum_{j=0}^n \log|t-x_j| - \frac{1}{n+1}
\sum_{j=0}^n \log|x-x_j|
$$ 

Define the potential function

$$
u_n(s) = \frac{1}{n+1} \sum_{j=0}^n \log|s-x_j|
$$ 

$u_n$ is a harmonic function in the complex plane away from $\{x_j\}$. We may think of each $x_j$ as a point charge of strength $\frac{1}{n+1}$, like an electron, and of $u_n$ as the potential induced by all the charges, whose gradient defines an electric field.

In terms of the potential 

$$
\log \gamma_n(x,t) = u_n(t) - u_n(x)
$$ 

and hence

$$
\left| \frac{\ell(x)}{\ell(t)} \right| = \ee^{-(n+1)[ u_n(t) - u_n(x)]}
$$

We require $u_n(t)-u_n(x) > 0$; in particular if

$$
\min_{x \in X, t \in \Gamma}[ u_n(t) - u_n(x)] \ge \log\alpha > 0, \qquad \textrm{for some $\alpha > 1$}
$$ 

then

$$
\left| \frac{\ell(x)}{\ell(t)} \right| \le \ee^{-(n+1)\log\alpha} = \alpha^{-n-1} \to 0
$$

and the interpolants converge exponentially, $\norm{f-p} = O(\alpha^{-n})$. We can write the condition as

$$
\min_{t \in \Gamma} u_n(t) - \max_{x \in X} u_n(x) \ge \log\alpha > 0, \qquad
\textrm{for some $\alpha > 1$}
$$ 

Hence convergence depends on the difference of values taken by the potential function on the set of points $X$ where the interpolant is to be evaluated and on a contour $\Gamma$ inside which $f$ is analytic. If $f$ is analytic everywhere, we can take take $\Gamma$ far out in the complex plane and we easily satisfy the above condition. But if there is a singularity close to the real line, $\Gamma$ cannot be too far away, and then we have to analyze the condition more carefully.

+++

## From discrete to continuous potentials

We can write the potential $u_n$ as a Lebesque-Stieltjes integral

$$u_n(s) = \int_{-1}^1 \log|s-\tau| \ud\mu_n(\tau)
$$ 

where $\mu_n$ is a measure consisting of a sum of Dirac delta functions, each of strength $1/(n+1)$

$$
\mu_n(\tau) = \frac{1}{n+1} \sum_{j=0}^n \delta(\tau - x_j)
$$ 

What is the limiting measure as $n \to \infty$ ? This limit depends on the grid point distribution. Note that

$$
\int_{-1}^1 \mu_n(\tau)\ud\tau = 1, \qquad \int_a^b \mu_n(\tau) \ud \tau =
\frac{\textrm{number of nodes in $[a,b]$}}{n+1}
$$ 

The second property can also be written as

$$
\mu_n(\tau) \Delta\tau = \frac{\textrm{number of nodes in $[\tau-\half\Delta\tau, \tau+ \half\Delta\tau]$}}{\textrm{total number of nodes}} 
$$

+++

### Uniform point distribution

The number of nodes in an interval of length $\Delta\tau$ is proportional to $\Delta\tau$ so that

$$
\mu(\tau) \Delta\tau = C \Delta\tau
$$ 

and using the condition $\int_{-1}^1\mu(\tau)\ud\tau=1$ we get $C=\half$ and

$$
\mu(\tau) = \frac{1}{2}
$$ 

The limiting potential is

$$
u(s) = -1 + \half \real [(s+1)\log(s+1) - (s-1)\log(s-1)]
$$ 

At the end-point and middle, it takes the values

$$
u(\pm 1) = -1 + \log(2), \qquad u(0) = -1
$$ 

The equipotential curve $u(s) = -1 + \log(2)$ is shown in red in figure and it cuts the imaginary axis at $\pm 0.52552491457\ii$ and passes through the end
points $x=\pm 1$.

:::{figure} matlab/equipot_uniform.svg
:width: 80%
Equipotential curves for uniform points. The red curve corresponds to $u(s) = -1 + \log 2$.
:::

If $f$ has a singularity outside the red curve, then we can take $\Gamma$ to be an equipotential curve $u(s) = u_0$ with $u_0 > -1 + \log(2)$. Then

$$
\min_{t \in \Gamma} u(t) - \max_{x \in X} u(x) = u_0 + 1 - \log(2) > 0
$$

and we have exponential convergence. But if the singularity is inside the red curve, then we cannot choose a curve $\Gamma$ that satisfies the above condition.

+++

### Chebyshev points

The Chebyshev points are uniformly distributed with respect to the variable $\theta \in [0,\pi]$ where $x = -\cos\theta$, so that

$$
\mu(\tau)\Delta\tau = C \Delta\theta
$$ 

and using the condition $\int_{-1}^1\mu(\tau)\ud\tau=1$, we get

$$
C = \frac{1}{\pi}, \qquad \mu(\tau) = C \left[ \dd{x}{\theta} \right]^{-1} = \frac{1} {\pi \sqrt{1-\tau^2}}
$$ 

The limiting potential is

$$
u(s) = \log|s + \ii \sqrt{1-s^2}| - \log 2
$$ 

For $s \in [-1,+1]$, $u(s) = -\log 2$, i.e., a constant. Thus $X = [-1,+1]$ is an equipotential curve of the potential function $u(s)$. For $u_0 > -\log 2$, the equipotential curve $u(s) = u_0$ is the Bernstein ellipse $E_\rho$ with $\rho = 2 \ee^{u_0}$ which encloses the set $X$.

:::{figure} matlab/equipot_chebyshev.svg
:width: 90%

Equipotential curves for chebyshev points. The potential is $u(s) = \log|s + \ii \sqrt{1-s^2}| - \log 2$.
:::

No matter how close the singularity of $f$ is to $X$, we can always find an equipotential curve $u(s) = \rho$ with $\rho > -\log(2)$, which encloses $X$ and inside which $f$ is analytic. If we take $\Gamma$ to be this equipotential curve then

$$
\min_{t \in \Gamma} u(t) - \max_{x \in X} u(x) = \rho + \log(2) > 0
$$

and we obtain exponential convergence.

+++

:::{exercise}
The figures of equipotential curves were generated with the following Matlab code. Rewrite it using Python.

```{literalinclude} matlab/equipot.m
```
:::
