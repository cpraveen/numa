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

# Gauss quadrature - II

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

Fix an integer $n \ge 1$. For a positive weight function $w(x)$, we want to compute an approximation to the integral

$$
I(f) = \int_a^b w(x) f(x) \ud x \approx I_n(f) = \sum_{j=1}^n w_j f(x_j)
$$

Moreover, we want the degree of precision of $I_n$ to be as large as
possible, i.e., $2n-1$. The main idea is to construct the Hermite
interpolant of $f(x)$ and take $I_n(f)$ to be the exact integral of the
Hermite interpolant.

## General case

Let $\{ \phi_n\}$ be an orthogonal family of polynomials on $[a,b]$ with
respect to the above weight function $w$. Denote the zeros of
$\phi_n(x)$ by 

$$
a < x_1 < x_2 < \ldots < x_n < b
$$ 

These zeroes will
come out as the nodes in the integration formula $I_n(f)$. We recall
some notation used previously. We write $\phi_n$ in terms of monomials

$$
\phi_n(x) = A_n x^n + B_n x^{n-1} + \ldots
$$ 

and define

$$
a_n = \frac{A_{n+1}}{A_n}, \qquad \gamma_n = \ip{\phi_n, \phi_n} = \int_a^b w(x) [\phi_n(x)]^2 \ud x
$$

+++

:::{prf:theorem}
For each $n \ge 1$, there is a numerical integration formula $I_n(f)$ of
degree of precision $2n-1$. Assuming $f \in \cts^{2n}[a,b]$

$$
\int_a^b w(x) f(x) \ud x = \sum_{j=1}^n w_j f(x_j) + \frac{\gamma_n}{A_n^2 (2n)!} f^{(2n)}(\eta), \qquad \eta \in [a,b]
$$

The nodes $\{x_j\}$ are the zeros of $\phi_n(x)$ and

$$
\label{eq:genintwj}
w_j = -\frac{a_n \gamma_n}{\phi_n'(x_j) \phi_{n+1}(x_j)}
$$
:::

:::{prf:proof}
**(a) Construction of the formula**: We will use Hermite
interpolation at the $n$ zeros of $\phi_n$ to construct the quadrature
rule. The interpolant is given by, see [](#sec:laghermint)

$$
H_n(x) = \sum_{j=1}^n f(x_j) h_j(x) + \sum_{j=1}^n f'(x_j) \tilde{h}_j(x) \in \poly_{2n-1}
$$

where

\begin{gather}
h_j(x) = [1 - 2\ell_j'(x_j)(x-x_j)] [\ell_j(x)]^2, \qquad \tilde{h}_j(x) = (x-x_j) [\ell_j(x)]^2 \\
\ell_j(x) = \prod_{i=1,i\ne j}^n \frac{x-x_i}{x_j - x_i}
\end{gather}

The Hermite interpolation error is

$$
e_n(x) 
&= f(x) - H_n(x) \\
&= [\psi_n(x)]^2 f[x_1,x_1,\ldots,x_n,x_n,x] \\
&= \frac{[\psi_n(x)]^2}{(2n)!} f^{(2n)}(\eta), \quad \eta \in [a,b]
$$

where 

$$
\psi_n(x) = (x-x_1)\ldots(x-x_n) = \frac{1}{A_n} \phi_n(x)
$$ 

Now

\begin{align}
\int_a^b w(x) f(x) \ud x 
&= \int_a^b w(x) H_n(x) \ud x + \int_a^n w(x) e_n(x) \ud x \\
&=: I_n(f) + E_n(f)
\end{align}

The degree of precision of $I_n$ is atleast $2n-1$ since $H_n$ is exact
for such polynomials. Moreover, for $f(x) = x^{2n}$,

$$
\begin{aligned}
E_n(x^{2n}) 
&= \int_a^b w(x) e_n(x) \ud x \\
&= \int_a^b w(x) \frac{[\psi_n(x)]^2}{(2n)!} f^{(2n)}(\eta) \ud x \\
&= \int_a^b w(x) \frac{[\psi_n(x)]^2}{(2n)!} (2n)! \ud x \\
&= \int_a^b w(x) [\psi_n(x)]^2 \ud x > 0
\end{aligned}
$$ 

Hence the degree of precision of $I_n$ is exactly $2n-1$.

**(b) Simpler formula for $I_n(f)$**: The integration formula is

$$
I_n(f) = \sum_{j=1}^n f(x_j) \int_a^b w(x) h_j(x) \ud x + \sum_{j=1}^n f'(x_j) \int_a^b w(x) \tilde{h}_j(x) \ud x
$$

Since, see [](#sec:implagint)

\begin{gather}
\ell_j(x) = \frac{\psi_n(x)}{(x-x_j) \psi_n'(x_j)} = \frac{\phi_n(x)}{(x-x_j) \phi_n'(x_j)} \\
\implies \tilde{h}_j(x) = (x-x_j) \ell_j(x) \cdot \ell_j(x) = \frac{\phi_n(x)}{\phi_n'(x_j)} \cdot \ell_j(x)
\end{gather}

Hence

$$
\int_a^b w(x) \tilde{h}_j(x) \ud x = \frac{1}{\phi_n'(x_j)}\int_a^b w(x) \phi_n(x) \ell_j(x) \ud x = 0, \qquad j=1,2,\ldots,n
$$

since degree$(\ell_j) = n-1$ and $\phi_n$ is orthogonal to all polynomials of degree $< n$. Hence

$$
I_n(f) = \sum_{j=1}^n w_j f(x_j), \qquad w_j = \int_a^b w(x) h_j(x) \ud x
$$

We do not need to know the derivatives $f'(x_j)$ !!!

**(c) Uniqueness of $I_n(f)$**: Suppose we have another numerical integration formula

$$
\int_a^b w(x) f(x) \ud x \approx \sum_{j=1}^n v_j f(z_j)
$$ 

that has degree of precision $\ge 2n-1$. For any $p \in \poly_{2n-1}$, the Hermite interpolant using nodes $\{ z_j \}$ is exact and

$$
p(x) = \sum_{j=1}^n p(z_j) g_j(x) + \sum_{j=1}^n p'(z_j) \tilde{g}_j(x)
$$

where

\begin{gather}
g_j(x) = [ 1 - 2 (x - z_j) \bar{\ell}_j (z_j) ] [\bar{\ell}_j(x)]^2 \\
\tilde{g}_j(x) = (x-z_j) [\bar{\ell}_j(x)]^2 = \frac{\bar{\ell}_j(x) \chi_n(x)}{\chi_n'(z_j)} \\
\bar{\ell}_j(x) = \prod_{i=1,i \ne j}^n \frac{x-z_i}{z_j - z_i}, \qquad \chi_n(x) = (x-z_1)\ldots(x-z_n)
\end{gather}

Since the quadrature formula is exact for $p \in \poly_{2n-1}$

$$
\sum_{j=1}^n v_j p(z_j) 
&= \int_a^b w(x) p(x) \ud x \\
&= \sum_{j=1}^n p(z_j) \int_a^b w(x) g_j(x) \ud x + \sum_{j=1}^n p'(z_j) \int_a^b w(x) \tilde{g}_j(x) \ud x
$$

Take $p(x) = \tilde{g}_i(x)$. Note that

$$
\tilde{g}_i(z_j) = 0, \qquad \tilde{g}_i'(z_j) = \delta_{ij}
$$ 

Hence we get 

$$
0 = 0 + \tilde{g}_i'(z_i) \int_a^b w(x) \tilde{g}_i(x) \ud x
$$

$$
\implies \int_a^b w(x) \tilde{g}_i(x) \ud x = \frac{1}{\chi_n'(z_i)} \int_a^b w(x) \chi_n(x) \bar{\ell}_i(x) \ud x = 0, \qquad i=1,2,\ldots,n
$$

Since degree$(\bar{\ell}_i)=n-1$ and they form a basis for $\poly_{n-1}$, we get

$$
\int_a^b w(x) \chi_n(x) p(x) = 0, \qquad \forall p \in \poly_{n-1}
$$

Hence $\chi_n$ is orthogonal to $\poly_{n-1} = \textrm{span}(\phi_0,\phi_1, \ldots,\phi_{n-1})$ and note that $\chi_n \in \poly_n$. But orthogonal polynomials are unique upto a multiplicative constant. Hence $\chi_n(x) = c \phi_n(x)$ and in particular their zeros are same, which imples that $z_j = x_j$.  Moreover, since degree$(h_i)=2n-1$

$$
w_i = \int_a^b w(x) h_i(x)\ud x = \sum_{j=1}^n v_j \underbrace{h_i(x_j)}_{\delta_{ij}} = v_i
$$

**(d) Weights are positive**: 

$$
\begin{aligned}
w_i &= \int_a^b w(x) h_i(x) \ud x \\
&= \int_a^b w(x) [1 - 2 \ell_i'(x_i)(x-x_i)] [\ell_i(x)]^2 \ud x \\
&= \int_a^b w(x) [\ell_i(x)]^2 \ud x - 2 \ell_i'(x_i) \underbrace{\int_a^b w(x) \tilde{h}_i(x) \ud x}_{=0} \\
&= \int_a^b w(x) [\ell_i(x)]^2 \ud x > 0
\end{aligned}
$$ 

**(e) Error formula**: Assume that $f \in \cts^{2n}[a,b]$. Then the integration error is 

$$
\begin{aligned}
E_n(f) 
&= \int_a^b w(x) e_n(x) \ud x \\
&= \int_a^b w(x) [\psi_n(x)]^2 f[x_1,x_1,\ldots,x_n,x_n,x] \ud x \\
&= f[x_1,x_1,\ldots,x_n,x_n,\xi] \int_a^b w(x) [\psi_n(x)]^2 \ud x, \quad \xi\in[a,b] \\
&= \frac{f^{(2n)}(\eta)}{(2n)!} \frac{1}{A_n^2} \int_a^b w(x) [\phi_n(x)]^2 \ud x \\
&= \frac{f^{(2n)}(\eta)}{(2n)!} \frac{1}{A_n^2} \gamma_n
\end{aligned}
$$ 

**(f) Formula for the weights**: We have already shown that the weights are positive. To obtain a simpler formula, take $f(x) = \ell_i(x)$ and note that degree($\ell_i)=n-1$, it can be exactly integrated

$$
\int_a^b w(x) \ell_i(x) \ud x = \sum_{j=1}^n w_j \ell_i(x_j) + E_n(\ell_i) = w_i + 0
$$

Hence (this was also shown in [](#sec:GLgeneral))

$$
w_i = \int_a^b w(x) \ell_i(x) \ud x = \frac{1}{\phi_n'(x_i)} \int_a^b \frac{w(x) \phi_n(x)}{(x-x_i)} \ud x
$$

Recall the Christoffel-Darboux identity from [](#sec:chrisdarb)

$$
\sum_{k=0}^n \frac{\phi_k(x) \phi_k(y)}{\gamma_k} = \frac{\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y)}{a_n \gamma_n (x-y)}, \qquad x \ne y
$$

Take $y = x_i$ and use $\phi_n(y)=\phi_n(x_i)=0$

$$
\sum_{k=0}^{n-1} \frac{\phi_k(x) \phi_k(x_i)}{\gamma_k} = -\frac{\phi_n(x) \phi_{n+1}(x_i)}{a_n \gamma_n (x-x_i)}
$$

Multiply by $w(x)$ and integrate both sides; by orthogonality of $\{ \phi_k \}$

$$
\int_a^b w(x) \phi_k(x) \ud x = 0, \qquad k \ge 1
$$ 

only the term $k=0$ in the sum survives, and since $\phi_0(x)=$ constant, we get

\begin{align}
\frac{1}{\gamma_0} \int_a^b w(x) \phi_0(x)\phi_0(x_i)\ud x = 1 &= - \frac{\phi_{n+1}(x_i)}{a_n \gamma_n} \int_a^b \frac{w(x) \phi_n(x)}{x-x_i} \ud x\\
&= - \frac{\phi_{n+1}(x_i)}{a_n \gamma_n} w_i \phi_n'(x_i)
\end{align}

which proves the formula for the weights.
:::

+++

## Gauss-Legendre quadrature

For weight function $w(x) = 1$ on $[-1,1]$ the orthogonal polynomials are the Legendre functions $P_n(x)$. The integral is approximated by

$$
I(f) = \int_{-1}^1 f(x) \ud x \approx I_n(f)  = \sum_{j=1}^n w_j f(x_j)
$$

The nodes $\{x_j\}$ are the zeros of $P_n(x)$ and the weights are obtained from
 [](#eq:genintwj)

$$
w_i = -\frac{2}{(n+1) P_n'(x_i) P_{n+1}(x_i)}, \qquad i=1,2,\ldots,n
$$

and the error is

\begin{gather}
E_n(f) = \frac{2^{2n+1} (n!)^4}{(2n+1) [(2n)!]^2} \frac{f^{(2n)}(\eta)}{(2n)!} = e_n \frac{f^{(2n)}(\eta)}{(2n)!} \\
e_n = \frac{2^{2n+1} (n!)^4}{(2n+1) [(2n)!]^2}
\end{gather}

The function [scipy.integrate.fixed_quad](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.fixed_quad.html) performs Gauss-Legendre quadrature.

:::{prf:example}

$$
I = \int_0^\pi \ee^x \cos(x) \ud x = -\half ( 1 + \ee^\pi )
$$ 

```{code-cell}
from scipy.integrate import fixed_quad

f = lambda x: exp(x)*cos(x)
a,b = 0.0,pi
qe = -0.5*(1.0 + exp(pi)) # Exact integral

N = arange(2,10)
err = zeros(len(N))
for (i,n) in enumerate(N):
    val = fixed_quad(f,a,b,n=n)
    err[i] = abs(val[0]-qe)
    print('%5d %24.14e' % (n,err[i]))

figure()
semilogy(N,err,'o-')
xlabel('n'), ylabel('Error');
```

The error decreases much faster than trapezoidal (see [](#ex:trapz1)) or Simpson rules (see [](#ex:simpexpcos)). With 8 nodes, the error has reached almost machine zero.
:::

+++

## Asymptotic behavior of error

Define

$$
M_m = \max_{x \in [-1,1]} \frac{|f^{(m)}(x)|}{m!}, \qquad m \ge 0
$$

For a large class of infinitely differentiable functions on $[-1,1]$

$$
\sup_m M_m < \infty
$$ 

In fact, for many functions

$$
\lim_{m \to \infty} M_m = 0
$$ 

e.g., analytic functions like $\ee^x$, $\cos(x)$, etc. For such functions, the Gauss quadrature has error

$$
|E_n(f)| \le e_n M_{2n}
$$ 

The speed of convergence depends on $e_n$.  Using Stirling approximation 

$$
n! \approx \ee^{-n} n^n \sqrt{2\pi n}
$$

we get 

$$
e_n \approx \frac{\pi}{4^n}
$$ 

E.g., $e_5 = 0.00293$ while $\pi/4^5 = 0.00307$. Hence we get the asymptotic error bound

$$
|E_n(f)| \le \frac{\pi}{4^n} M_{2n}
$$ 

This is exponential convergence of the error. Compare this to 

\begin{align}
\textrm{Trapezoid:} \qquad &  |E_n(f)| \le c/n^2 \\
\textrm{Simpson:}   \qquad & |E_n(f)| \le c/n^4
\end{align}

Hence for nice functions Gauss quadrature is very accurate as it converges exponentially fast.

## Relation to minimax approximation

We can see the optimal accuracy of Gauss quadrature by its relation to the best approximation.

:::{prf:theorem}
Assume $[a,b]$ is finite. Then the error of Gauss quadrature

$$
|E_n(f)| = |I(f) - I_n(f)| \le 2 \rho_{2n-1}(f) \int_a^b w(x) \ud x
$$

with $\rho_{2n-1} = \inf_{p \in \poly_{2n-1}} \norm{f - p}_\infty$ is
the minimax error.
:::

:::{prf:proof}
We have $E_n(p)=0$ for any $p \in \poly_{2n-1}$ and the error $E_n$ is a linear operator. Take $p = q_{2n-1}^* \in \poly_{2n-1}$ be the minimax approximation to $f$ on $[a,b]$. Then 

$$
\begin{aligned}
E_n(f) 
&= E_n(f) - E_n(q_{2n-1}^*) \\
&= E_n(f - q_{2n-1}^*) \\
&= \int_a^b w(x)[f(x) - q_{2n-1}^*(x)] \ud x - \sum_{j=1}^n w_j[f(x_j) - q_{2n-1}^*(x_j)]
\end{aligned}
$$ 

so that

$$
|E_n(f)| \le \norm{f - q_{2n-1}^*}_\infty\left[ \int_a^b w(x) \ud x + \sum_{j=1}^n w_j \right]
$$

since $w(x) \ge 0$ and $w_j > 0$. As the quadrature rule is exact for constant function 

$$
\int_a^b w(x) \ud x = \sum_{j=1}^n w_j
$$ 

which completes the proof.
:::

:::{exercise}
Apply [](#thm:InIdense) to Gauss-Legendre quadrature and show that the numerical integral $I_n(f)$ converges for every continuous function $f$.
:::

:::{prf:example} Integrals on $[0,\infty)$
Integrals of the form

$$
\int_0^\infty \ee^{-x} f(x) \ud x \approx \sum_{j=1}^n w_j f(x_j)
$$

can be computed using Gauss-Laguerre quadrature. Laguerre polynomials are orthogonal on $[0,\infty)$ wrt the weight $w(x) = \ee^{-x}$, see [](#sec:laguerrepoly). The quadrature points $x_j$ are the roots of $L_n(x)$ and the weights are given by

$$
w_j = \frac{x_j}{(n+1)^2 [L_{n+1}(x_j)]^2}
$$

We can use Gauss-Laguerre quadrature for integrals of the form

$$
I = \int_0^\infty g(x) \ud x
$$ 

by writing it as

$$
I = \int_0^\infty \ee^{-x}[ \ee^x g(x)] \ud x = \int_0^\infty \ee^{-x} f(x) \ud x
$$

See [@Atkinson2004], page 308-309.

$$
\int_0^\infty \frac{x}{\ee^x - 1} \ud x = \frac{\pi^2}{6}, \qquad
\int_0^\infty \frac{x}{(1 + x^2)^5} \ud x = \frac{1}{8}, \qquad
\int_0^\infty \frac{1}{1 + x^2} \ud x = \frac{\pi}{2}
$$

```{code-cell}
f1  = lambda x: x/(exp(x) - 1)
Ie1 = pi**2/6.0

f2  = lambda x: x/(1.0 + x**2)**5
Ie2 = 1.0/8.0

f3  = lambda x: 1.0/(1.0 + x**2)
Ie3 = 0.5*pi

n=2
for i in range(6):
    x,w = polynomial.laguerre.laggauss(n)
    I1 = sum(exp(x)*f1(x)*w)
    I2 = sum(exp(x)*f2(x)*w)
    I3 = sum(exp(x)*f3(x)*w)
    print('%5d %15.6e %15.6e %15.6e' % (n,I1-Ie1,I2-Ie2,I3-Ie3))
    n = 2*n
```
:::

+++

## Singular integrands

:::{prf:example} Singular integral

$$
I = \int_0^1 \sqrt{x} \ud x
$$ 

The integrand does not have derivatives at $x=0$. We compare trapezoid and Gauss-Legendre quadratures.

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

$f(x)$ smooth with $f(0) \ne 0$, the error of Gauss quadrature is (Donaldson and Elliott, (1972))

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

Now, as seen in second and third columns, Gauss quadrature is not so good because singularity is inside the interval. A trick here is to break the integral into two intervals and compute each with its own Gauss quadrature

$$
\int_0^\half \sqrt{|x - \shalf|} \ud x + \int_\half^1 \sqrt{x - \shalf} \ud x \approx I_{n/2} + I_{n/2}
$$

Now the singularity is at the end-points and Gauss quadrature is more
accurate, as seen in fourth and fifth columns.
:::

:::{prf:example}
\begin{align}
I_1 &= \int_0^5 \frac{\ud x}{1 + (x - \pi)^2} = \arctan(5-\pi) + \arctan(\pi) \\
I_2 &= \int_0^{2\pi} \ee^{-x} \sin(50 x) \ud x = \frac{100  \exp(-\pi)  \sinh(\pi)}{2501}
\end{align}

```{code-cell}
f1 = lambda x: 1.0/(1.0 + (x - pi)**2)
a1, b1 = 0.0, 5.0
I1 = arctan(5-pi) + arctan(pi)
f2 = lambda x: exp(-x) * sin(50*x)
a2, b2 = 0.0, 2.0*pi
I2 = 100 * exp(-pi) * sinh(pi) / 2501

n = 2
for i in range(7):
    val = fixed_quad(f1,a1,b1,n=n)
    err1 = abs(val[0]-I1)
    val = fixed_quad(f2,a2,b2,n=n)
    err2 = abs(val[0]-I2)
    print('%5d %14.6e %14.6e' % (n,err1,err2))
    n = 2*n
```

We get error of order machine zero for the first function for $n > 32$. For the second, the convergence is slow because the integrand is highly oscillatory; but once the function is well resolved by the sufficient quadrature points, the convergence is rapid.
:::
