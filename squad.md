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

# Special techniques

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## Singular integrands

:::{prf:example} Singular integral

$$
I = \int_0^1 \sqrt{x} \ud x
$$ 

The integrand does not have derivatives at $x=0$. 

```{code-cell}
from scipy.integrate import fixed_quad,trapezoid

f = lambda x: sqrt(x)
a,b = 0.0,1.0
qe = 2.0/3.0 # Exact integral

n,N = 2,10
e1,e2,nodes = zeros(N),zeros(N),zeros(N)
for i in range(N):
    x = linspace(a,b,n)
    val = trapezoid(f(x),dx=(b-a)/(n-1))
    e1[i] = abs(val - qe)
    val = fixed_quad(f,a,b,n=n)
    nodes[i] = n
    e2[i] = abs(val[0]-qe)
    if i>0:
        print('%5d %20.10e %12.4e %20.10e %12.4e' % 
              (n,e1[i],e1[i-1]/e1[i],e2[i],e2[i-1]/e2[i]))
    else:
        print('%5d %20.10e %12.4e %20.10e %12.4e' % 
              (n,e1[i],0,e2[i],0))
    n = 2*n

figure()
loglog(nodes,e1,'o-')
loglog(nodes,e2,'*-')
legend(('Trapezoid','Gauss'))
xlabel('n'), ylabel('Error');
```

With trapezoidal rule, we get ($2^{1.5} \approx 2.825$)

$$
|E_n(f)| \approx \frac{c}{n^{3/2}}
$$ 

while with Gauss quadrature

$$
|E_n(f)| \approx \frac{c}{n^3}
$$ 

Gauss quadrature still performs well when the singularity is at the end-point of the interval.
:::

:::{prf:remark}

For integrals of the form

$$
I = \int_0^1 x^\alpha f(x) \ud x, \qquad \alpha > -1, \textrm{ fractional}
$$

$f(x)$ smooth with $f(0) \ne 0$, the error of Gauss quadrature is

$$
E_n(f) \approx \frac{c(f,\alpha)}{n^{2(1+\alpha)}}
$$
:::

:::{prf:example}

$$
\int_0^1 \sqrt{|x - \shalf|} \ud x = \frac{\sqrt{2}}{3} \approx I_n
$$ 

```{code-cell}
f = lambda x: sqrt(abs(x-0.5))
a,b = 0.0,1.0
c = 0.5*(a+b)
qe = sqrt(2.0)/3.0 # Exact integral

n,N = 2,7
err1,err2,nodes = zeros(N),zeros(N),zeros(N)
for i in range(N):
    nodes[i] = n
    val = fixed_quad(f,a,b,n=n)
    err1[i] = abs(val[0]-qe)
    val = fixed_quad(f,a,c,n=n/2) + fixed_quad(f,c,b,n=n/2)
    err2[i] = abs(val[0]+val[2]-qe)
    if i>0:
        print('%5d %14.6e %12.4e %14.6e %12.4e' %
              (n,err1[i],err1[i-1]/err1[i],err2[i],err2[i-1]/err2[i]))
    else:
        print('%5d %14.6e %12.4e %14.6e %12.4e' % (n,err1[i],0,err2[i],0))
    n = 2*n

figure()
loglog(nodes,err1,'o-',nodes,err2,'s-')
xlabel('n'), ylabel('Error')
legend(('Gauss','Gauss+Gauss'));
```

Now Gauss quadrature is not so good because singularity is inside the interval. A trick here is to break it into two intervals and compute each with its own Gauss quadrature

$$
\int_0^\half \sqrt{|x - \shalf|} \ud x + \int_\half^1 \sqrt{x - \shalf} \ud x \approx I_{n/2} + I_{n/2}
$$

Now the singularity is at the end-points and Gauss quadrature is more
accurate.
:::

## Change of variable

Many times, a change of variable can improve the regularity of the
integrand.

:::{prf:example}

$$
I = \int_0^b \frac{f(x)}{\sqrt{x}} \ud x, \qquad \textrm{$f$ is smooth}
$$

With 

$$
x = u^2, \qquad 0 \le u \le \sqrt{b}
$$

$$
I = 2 \int_0^{\sqrt{b}} f(u^2) \ud u
$$

:::

:::{prf:example}

$$
\int_0^1 \sin(x) \sqrt{1-x^2} \ud x = 2\int_0^1 u^2 \sqrt{2-u^2} \sin(1-u^2) \ud u
$$

by change of variables $u = \sqrt{1-x}$. The first form has singular
derivative at $x=1$ while the second form has infinitely smooth
integrand.
:::

:::{prf:example} Infinite interval of integration

$$
I = \int_1^\infty \frac{f(x)}{x^p} \ud x, \qquad p > 1
$$ 

Assume $\lim_{x \to \infty} f(x)$ exists and $f(x)$ is smooth on $[1,\infty)$.

$$
x = \frac{1}{u^\alpha}, \qquad \ud x = -\frac{\alpha}{u^{1+\alpha}} \ud u, \quad \alpha > 0
$$

$$
I = \alpha \int_0^1 u^{(p-1)\alpha-1} f(1/u^\alpha) \ud u
$$ 

Integrand can be singular at $u=0$, so pick $\alpha$ so that the exponent $(p-1)\alpha-1$ is large $(\ge 1)$. As a specific example consider

$$
I = \int_1^\infty \frac{f(x)}{x\sqrt{x}} \ud x
$$

$$
x = \frac{1}{u^4} \quad\implies\quad I = 4 \int_0^1 u f(1/u^4) \ud u
$$

For $x \to \infty$, $f(x)$ behaves like

$$
f(x) = c_0 + \frac{c_1}{x} + \frac{c_2}{x^2} + \ldots
$$

$$
uf(1/u^4) = c_0 u + c_1 u^5 + c_2 u^9 + \ldots
$$ 

and the integrand behaves nicely as $u \to 0$.
:::

:::{prf:example}

$$
I = \int_{-1}^1 \frac{f(x)}{\sqrt{1-x^2}} \ud x
$$ 

For this weight function, the orthogonal polynomials are Chebyshev polynomials. The integration nodes are roots of $T_n(x) = \cos(n \cos^{-1}x)$

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
I - \int_0^\pi f(\cos\theta) \ud\theta \approx \frac{\pi}{n} \sum_{j=1}^n f(\cos\theta_j), \qquad \theta_j = \frac{(2j-1)\pi}{2n}
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
I_1 = \int_0^\epsilon (\sum_j a_j x^j) \log(x) \ud x = \sum_{j=0}^\infty \frac{a_j \epsilon^{j+1}}{j+1} \left[ \log(\epsilon) - \frac{1}{j+1} \right]
$$

:::{prf:example}

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

:::{prf:example}

$$
I(f) = \int_0^b f(x) \log(x) \ud x
$$ 

Let $f_n$ be piecewise linear approximation on a partition

$$
x_j = jh, \quad j=0,1,\ldots,n \qquad\textrm{where}\qquad h = \frac{b}{n}
$$

and

$$
f_n(x) = \frac{1}{h}[ (x_j - x) f_{j-1} + (x-x_{j-1}) f_j], \qquad x_{j-1} \le x \le x_j
$$

for which we know the error estimate

$$
\norm{f-f_n}_\infty \le \frac{h^2}{8} \norm{f''}_\infty
$$ 

Then

$$
|I(f) - I_n(f)| \le \frac{h^2}{8} \norm{f''}_\infty \int_0^b |\log(x)| \ud x \to 0
$$

$I_n$ can be computed analytically. In particular, in the first sub-interval $[x_0,x_1]=[0,h]$ we have

$$
\int_{x_0}^{x_1} \log(x) f_n(x) \ud x = \frac{f_0}{h} \int_0^h (h-x)\log(x) \ud x + \frac{f_1}{h} \int_0^h x \log(x) \ud x
$$

and these are finite since

$$
\int_0^h \log(x) \ud x = h[\log(h)-1], \qquad \int_0^h x\log(x)\ud x = \frac{h^2}{4}[ 2\log(h) - 1 ]
$$
:::
