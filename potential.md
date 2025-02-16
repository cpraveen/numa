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

# Potential theory

```{include} math.md
```

```{code-cell}
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

If $f(x)$ is analytic inside a curve $\Gamma$ enclosing the interpolation points $\{ x_j \}$, the error in polynomial interpolant is

$$
f(x) - p(x) = \frac{1}{2\pi\ii} \oint_{\Gamma} \frac{\ell(x)}{\ell(t)}          \frac{f(t)}{(t-x)}
\ud t
$$

We have seen that if $f$ is analytic inside the stadium $S$, the error goes to zero exponentially.  If the function has a singularity within the stadium S, i.e., close to the real line, then we have to do more analysis and the effect of point distribution also has to be considered in the analysis. 

Define

$$
\gamma_n(x,t) = \left| \frac{\ell(t)}{\ell(x)} \right|^{1/(n+1)} =
\left( \prod_{j=0}^n \frac{ |t-x_j| }{ |x- x_j| } \right)^{1/(n+1)}
$$ 

then

$$
\left| \frac{\ell(x)}{\ell(t)} \right| = [\gamma_n(x,t)]^{-n-1}
$$ 

If $\gamma_n(x,t) > 1$, for all $t \in \Gamma$, then we get exponential convergence as $n \to \infty$. Define

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
u_n(s) = \frac{1}{n+1} \sum_{j=0}^n \log|s-x_j| = \frac{1}{n+1} \log|\omega_n(s)|
$$ 

Recall that $\omega_n(x)$ is the polynomial factor in the error estimate.  $u_n$ is a harmonic function in the complex plane away from $\{x_j\}$. We may think of each $x_j$ as a point charge of strength $\frac{1}{n+1}$, like an electron, and of $u_n$ as the potential induced by all the charges, whose gradient defines an electric field.

In terms of the potential 

$$
\log \gamma_n(x,t) = u_n(t) - u_n(x)
\qquad \textrm{and} \qquad
\left| \frac{\ell(x)}{\ell(t)} \right| = \ee^{-(n+1)[ u_n(t) - u_n(x)]}
$$

For convergence, we require $u_n(t)-u_n(x) > 0$; in particular if

$$
\min_{x \in X, t \in \Gamma}[ u_n(t) - u_n(x)] \ge \log\alpha > 0, \qquad \textrm{for some $\alpha > 1$}
$$ 

then

$$
\left| \frac{\ell(x)}{\ell(t)} \right| \le \ee^{-(n+1)\log\alpha} = \alpha^{-n-1} \to 0
$$

and the interpolants converge exponentially, $\norm{f-p}_\infty = O(\alpha^{-n})$. We can write the condition as

$$
\min_{t \in \Gamma} u_n(t) - \max_{x \in X} u_n(x) \ge \log\alpha > 0, \qquad
\textrm{for some $\alpha > 1$}
$$ 

Hence convergence depends on the difference of values taken by the potential function on the set of points $X$ where the interpolant is to be evaluated and on a contour $\Gamma$ inside which $f$ is analytic. If $f$ is analytic everywhere, we can take take $\Gamma$ far out in the complex plane and we easily satisfy the above condition since

$$
u_n(t) \approx \log|t|, \qquad |t| \gg 1
$$

But if there is a singularity close to the real line, $\Gamma$ cannot be taken too far away, and then we have to analyze the condition more carefully.

+++

:::{prf:example}
Plot contours of $u_n$ corresponding to uniform and Chebyshev points for $n=50$. The following function computes $u_n(s)$.

```{code-cell}
# x =[x0,x1,...,xn] are the interpolation points
def potential(x,xp):
    fp = zeros(shape(xp))
    for xi in x:
        fp += log(abs(xp - xi))
    return fp/len(x)
```

Make a grid in complex plane to draw contours of $u_n$.

```{code-cell}
xg = linspace(-1.5,1.5,200)
X, Y = meshgrid(xg,xg)
Z = X + 1j * Y
```

**Uniformly spaced points**

```{code-cell}
n = 50
x = linspace(-1,1,n+1)
contour(X,Y,potential(x,Z),levels=11)
contour(X,Y,potential(x,Z),levels=[-1+log(2)],colors='blue',linestyles='solid')
plot([-1,1],[0,0],'r-')
title('Contours of $u_{'+str(n)+'}$ for uniformly spaced points');
```

The blue color shows the contour line $u_n(s) = -1 + \log(2)$ which approximtely passes through the points $x = \pm 1$; the red line is the real interval $[-1,1]$.

**Chebyshev points**

```{code-cell}
n = 50
theta = linspace(0,pi,n+1)
x = cos(theta)
contour(X,Y,potential(x,Z),levels=11)
plot([-1,1],[0,0],'r-')
title('Contours of $u_{'+str(n)+'}$ for Chebyshev points');
```

It appears that the real interval $[-1,1]$ is a contour line.
:::

+++

## From discrete to continuous potentials

Since we are interested in the case $n \to \infty$, we can examine the potential in this limit.  We can write the potential $u_n$ as a Lebesque-Stieltjes integral

$$
u_n(s) = \int_{-1}^1 \log|s-\tau| \mu_n(\tau) \ud\tau
$$ 

where $\mu_n$ is a measure consisting of a sum of Dirac delta functions, each of strength $1/(n+1)$

$$
\mu_n(\tau) = \frac{1}{n+1} \sum_{j=0}^n \delta(\tau - x_j)
$$ 

What is the limiting measure as $n \to \infty$ ? The precise notion of convergence appropriate for this limit is known as "weak-star" convergence. This limit depends on the grid point distribution which is characterized by $\mu_n$. Note that

$$
\int_{-1}^1 \mu_n(\tau)\ud\tau = 1, \qquad \int_a^b \mu_n(\tau) \ud \tau =
\frac{\textrm{number of nodes in $[a,b]$}}{n+1}
$$ 

In the limit, we need

$$
\int_{-1}^1 \mu(\tau)\ud\tau = 1, \qquad \int_a^b \mu(\tau) \ud \tau =
\frac{\textrm{number of nodes in $[a,b]$}}{n+1}
$$ 

The second property can also be written as

$$
\mu(\tau) \ud\tau = \frac{\textrm{number of nodes in $[\tau-\half\ud\tau, \tau+ \half\ud\tau]$}}{\textrm{total number of nodes}} 
$$

+++

### Uniform point distribution

The number of nodes in an interval of length $\Delta\tau$ is proportional to $\Delta\tau$ so that

$$
\mu(\tau) \ud\tau = C \ud\tau
$$ 

and using the condition $\int_{-1}^1\mu(\tau)\ud\tau=1$ we get $C=\half$ and

$$
\mu(\tau) = \frac{1}{2}
$$ 

The limiting potential is

\begin{align}
u(s) 
&= \half \int_{-1}^1 \log|s-\tau| \ud\tau \\
&= -1 + \half \textrm{Real} [(s+1)\log(s+1) - (s-1)\log(s-1)]
\end{align}

At the end-points and middle, it takes the values

$$
u(\pm 1) = -1 + \log(2), \qquad u(0) = -1
$$ 

The equipotential curve $u(s) = -1 + \log(2)$ is shown in red in figure and it cuts the imaginary axis at $\pm 0.52552491457\ii$ and passes through the end
points $x=\pm 1$.

:::{figure} matlab/equipot_uniform.svg
:width: 80%
Equipotential curves for uniform points. The red curve corresponds to $u(s) = -1 + \log 2$.
:::

If $f$ has a singularity outside the red curve, then we can take $\Gamma$ to be an equipotential curve 

$$
\Gamma = \{ s \in \complex : u(s) = u_0 > -1 + \log(2) \}
$$ 

and

$$
u(x) \le -1 + \log(2), \qquad x \in X
$$

Then

$$
\min_{t \in \Gamma} u(t) - \max_{x \in X} u(x) \ge u_0 + 1 - \log(2) > 0
$$

and we have exponential convergence. But if the singularity is inside the red curve, then we cannot choose a curve $\Gamma$ that satisfies the above condition; interpolation at uniform nodes will not converge.

+++

:::{prf:example}
The function $f(x) = \frac{1}{1 + 16x^2}$ is not analytic inside the red curve, e.g., at $x = \pm 0.25\ii$,  and interpolating this at uniform points does not converge. The function $f(x) = \frac{1}{1+\tfrac{50}{18}x^2}$ is analytic inside the red curve, it has poles at $x = \pm 0.6\ii$.

```{code-cell}
f  = lambda x: 1.0/(1.0 + 100/36*x**2)
xe = linspace(-1,1,1000)

N = 30
for i in range(4):
    subplot(2,2,i+1)
    x = linspace(-1,1,N+1)
    y = f(x)
    ye = barycentric_interpolate(x,y,xe)
    print("N, Error = %3d %12.4e" % (N,abs(ye - f(xe)).max()))
    plot(x,y,'.'), plot(xe,ye)
    text(0,0.5,'N='+str(N),ha='center')
    N = N + 10
```

The errors initially decrease but then start to increase.  In principle, we should expect convergence, but due to round-off errors caused by the large difference in barycentric weights, we see errors starting to increase for larger degrees.
:::

+++

:::{warning}
Even for functions analytic in whole complex plane, uniform point interpolation will not succeed at very high degrees due to round-off errors in their evaluation; even barycentric form does not help.
:::

+++

### Chebyshev points

The Chebyshev points are uniformly distributed with respect to the variable $\theta \in [0,\pi]$ where $x = \cos\theta$, so that

$$
\mu(\tau)\ud\tau = C \ud\theta
$$ 

and using the condition 

$$
\int_{-1}^1\mu(\tau)\ud\tau = 1 = \int_0^\pi C \ud\theta
$$

we get

$$
C = \frac{1}{\pi}, \qquad \mu(\tau) = C \left[ \dd{\tau}{\theta} \right]^{-1} = \frac{1} {\pi \sqrt{1-\tau^2}}
$$ 

The limiting potential is

$$
u(s) = \int_{-1}^1 \frac{\log|s-\tau|}{\pi\sqrt{1-\tau^2}} \ud\tau = \log|s + \ii \sqrt{1-s^2}| - \log 2
$$ 

For 

$$
s \in [-1,+1], \qquad u(s) = \log\sqrt{s^2 + (1-s^2)} - \log 2 = -\log 2 = \textrm{const}
$$

Thus 

> $X = [-1,+1]$ is an equipotential curve of the potential function $u(s)$. 


:::{figure} matlab/equipot_chebyshev.svg
:width: 90%

Equipotential curves for chebyshev points. The potential is $u(s) = \log|s + \ii \sqrt{1-s^2}| - \log 2$.  For $u_0 > -\log 2$, the equipotential curve $u(s) = u_0$ is the Bernstein ellipse $E_\rho$ with $\rho = 2 \ee^{u_0}$ which encloses the set $X$.
:::

This property helps us.  No matter how close the singularity of $f$ is to $X$, we can always find an equipotential curve $u(s) = \rho$ with $\rho > -\log(2)$, which encloses $X$ and inside which $f$ is analytic. If we take $\Gamma$ to be this equipotential curve then

$$
\min_{t \in \Gamma} u(t) - \max_{x \in X} u(x) = \rho + \log(2) > 0
$$

and we obtain exponential convergence.

Thus, interpolation with Chebyshev points will converge exponentially for real analytic functions. Note that the form of $\mu$ indicates that what is needed is that the interpolation points $\{ x_j \}$ must be clustered near the end-points with a density of $(1-x^2)^{-1/2}$.

+++

:::{exercise}
The figures of equipotential curves were generated with the following Matlab code. Rewrite it using Python.

```{literalinclude} matlab/equipot.m
```
:::

:::{prf:remark}
So far we have considered infinitely differentiable or analytic functions. Chebyshev interpolants will converge for functions which are less nice than analytic, e.g., those with a few continuous derivatives. The proof of this requires somewhat different techniques.
:::
