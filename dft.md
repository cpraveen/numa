---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Trigonometric interpolation

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

A smooth function $u : [0,2\pi] \to \re$ can be written in terms of its Fourier series

$$
u(x) = Su(x) = \sum_{k=-\infty}^\infty \hat u_k \phi_k(x), \qquad \phi_k(x) = \ee^{\ii k x}
$$

where

$$
\label{eq:hatuk}
\hat u_k = \frac{1}{2\pi} \int_0^{2\pi} u(x) \ee^{-\ii k x} \ud x
$$

is the **continuous Fourier transform**. Recall that $\{ \phi_k \}$ is an orthogonal family of functions.

We can truncate the series $Su$ to obtain approximations. But this is not practically useful since we cannot compute the coefficients $\hat u_k$; either because we may not know $u(x)$ for all $x$ and also because the integrals cannot be computed exactly.

An approximation method which makes use of the function values at a discrete set of points only, will be more useful in computations.

## Interpolation problem

For any **even integer** $N>0$, consider the set of points in $[0,2\pi]$

$$
x_j = \frac{2\pi j}{N}, \qquad j=0,1,\ldots,N-1
$$ 

which are uniformly spaced with spacing

$$
h = \frac{2\pi}{N}
$$

```{code-cell}
:tags: remove-input
N = 9
h = 2*pi/N
x = h * arange(0,N)
figure(figsize=(8,2))
plot(x,0*x,'-o')
d = 0.005
text(0,d,"$x=0$",ha='center',va='bottom')
text(x[-1],d,"$x=2\\pi$",ha='center',va='bottom')
text(0,-d,"$x_0$",ha='center',va='top')
text(h,-d,"$x_1$",ha='center',va='top')
text(2*h,-d,"$x_2$",ha='center',va='top')
text(x[-3],-d,"$x_{N-2}$",ha='center',va='top')
text(x[-2],-d,"$x_{N-1}$",ha='center',va='top')
xticks([]), yticks([]);
```

Note that $x_0 = 0$ and the last point $x_{N-1} = 2\pi - h < 2\pi$; we do not include the point $x=2\pi$ in the grid since by periodicity, the value there is same as at $x=0$.

In Python, we can generate the grid points like this.

```python
h = 2 * pi / N
x = h * arange(0,N)
```

The degree $N/2$ trigonometric polynomial is

$$
I_N u(x) = \sumf \tilu_k \phi_k(x) = \sumf \tilu_k \ee^{\ii k x}
$$ 

We will determine the $N$ coefficients $\tilu_k$ from the interpolation conditions

$$
u_j = u(x_j) = I_N u(x_j), \qquad j=0,1,\ldots,N-1
$$ 

This is a system of $N \times N$ equations which requires solving a matrix. But we can get the solution without explicitly solving the matrix problem.

First prove the orthogonality relation 

$$
\label{eq:discorth}
\boxed{
\frac{1}{N} \sum_{j=0}^{N-1} \ee^{-\ii p x_j} = \begin{cases}
1, & p = Nm, \quad m=0, \pm 1, \pm 2, \ldots \\
0, & \textrm{otherwise}
\end{cases}
}
$$ 

Then computing

\begin{align}
\sum_{j=0}^{N-1} u(x_j) \ee^{-\ii k x_j} 
&= \sum_{j=0}^{N-1} \clr{red}{I_N u(x_j)} \ee^{-\ii k x_j} \\
&= \sum_{j=0}^{N-1} \clr{red}{\sum_{l=-N/2}^{N/2-1} \tilde u_l \ee^{\ii l x_j}} \ee^{-\ii k x_j} \\
&= \sum_{l=-N/2}^{N/2-1} \tilde u_l \clr{blue}{\sum_{j=0}^{N-1} \ee^{\ii (l-k)x_j}} \\
&= \sum_{l=-N/2}^{N/2-1} \tilde u_l \clr{blue}{N \delta_{lk}}
\end{align}

we obtain the discrete Fourier coefficients

$$
\label{eq:tildeuk}
\boxed{\tilu_k = \frac{1}{N} \sum_{j=0}^{N-1} u_j \ee^{-\ii k x_j}, \qquad k=-N/2, \ldots, N/2 - 1}
$$ 

This is known as the **discrete Fourier transform (DFT)** of $\{ u_j \}$.

The interpolation conditions

$$
\boxed{u_j = \sumf \tilu_k \ee^{\ii k x_j}, \qquad 0 \le j \le N-1}
$$

will be called the **inverse DFT (IDFT)** of $\{ \tilu_k \}$.

+++

:::{exercise}
Prove the discrete orthogonality relation [](#eq:discorth). This implies

$$
\frac{1}{N} \sum_{j=0}^{N-1} \ee^{\ii (l-k)x_j} = \delta_{lk}, \qquad - N/2 \le l, k \le N/2-1
$$
:::

+++

:::{prf:remark}
If we apply the Trapeoidal rule of integration to [](#eq:hatuk) using the uniformly spaced points, we obtain [](#eq:tildeuk).
:::

+++

### Lagrange form$^*$
We can derive a Lagrange form of the interpolation as follows.

$$
\begin{aligned}
I_N u(x) 
&= \sumf \tilu_k \ee^{\ii k x} \\
&= \frac{1}{N} \sumf \sum_{j=0}^{N-1} u(x_j) \ee^{\ii k (x-x_j)} \\
&= \sum_{j=0}^{N-1} u(x_j) \frac{1}{N} \sumf \ee^{\ii k (x-x_j)} \\
&= \sum_{j=0}^{N-1} u(x_j) \psi_j(x)
\end{aligned}
$$ 

where

$$
\psi_j(x) = \frac{1}{N} \sumf \ee^{\ii k (x-x_j)}
$$ 

The functions $\psi_j$ are trigonometric polynomials in $S_N$ that satisfy

$$
\psi_j(x_l) = \delta_{jl}, \qquad 0 \le j,l \le N-1
$$

:::{prf:remark}
Why did we not include the term corresponding to $k=N/2$ ? We already have $N$ equations for the $N$ unknowns. Moreover, the exponentials

$$
\phi_{-N/2}(x) = \ee^{-\ii Nx/2}, \qquad \phi_{N/2}(x) = \ee^{\ii Nx/2}
$$

take the same values at all the nodes $x_j$ and hence the two functions are not independent on the grid.
:::

Let us define the set of trigonometric functions

$$
S_N = \left\{ v : [0,2\pi] \to \re, \quad v = \sumf a_k \phi_k, \quad a_k \in \complex \right\}
$$

:::{prf:lemma}
If $u \in S_N$, then $I_N u = u$.
:::

:::{prf:proof}
If $u \in S_N$ then it is of the form

$$
u = \sumf c_k \phi_k
$$

for some coefficients $c_k$. Show that $\tilu_k = c_k$.
:::

Define the discrete inner product

$$
\ip{u,v}_N = \frac{2\pi}{N} \sum_{j=0}^{N-1} u(x_j) \conj{v(x_j)}
$$

and the usual continuous one

$$
\ip{u,v} =  \int_0^{2\pi} u(x) \conj{v(x)} \ud x
$$

In fact, the discrete one can be obtained if we apply Trapezoidal rule to the continuous inner product.

:::{prf:lemma}
For any $u,v \in S_N$ 

$$
\ip{u,v}_N = \ip{u,v}
$$
:::

:::{prf:proof}
Let 

$$
u = \sumf a_k \phi_k, \qquad v = \sumf b_k \phi_k
$$ 

Then

$$
\begin{aligned}
\ip{u,v} 
&= \int_0^{2\pi} u(x) \conj{v(x)} \ud x \\
&= \sum_k \sum_l a_k \conj{b_l} \underbrace{\int_0^{2\pi} \phi_k(x) \conj{\phi_l(x)}
\ud x}_{2\pi \delta_{kl}} \\
&= 2\pi \sumf a_k \conj{b_k}
\end{aligned}
$$ 

and 

$$
\begin{aligned}
\ip{u,v}_N 
&= \frac{2\pi}{N} \sum_{j=0}^{N-1} u(x_j) \conj{v(x_j)} \\
&= \frac{2\pi}{N} \sum_{j=0}^{N-1} u(x_j) \sumf \conj{b_k} \ee^{-\ii k x_j} \\
&= 2\pi \sumf \left( \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \ee^{-\ii k x_j} \right) \conj{b_k} \\
&= 2\pi \sumf a_k \conj{b_k}
\end{aligned}
$$
:::

As a consequence, $\ip{\cdot,\cdot}_N$ is an inner product on $S_N$ and

$$
\norm{u}_N := \sqrt{\ip{u,u}_N} = \sqrt{\ip{u,u}} = \norm{u}, \qquad \forall u \in S_N
$$

:::{prf:lemma}
For any continuous function $u$, the interpolant $I_N u$ satisfies

$$
\ip{I_N u, v}_N = \ip{u,v}_N, \qquad \forall v \in S_N
$$

:::

:::{prf:proof}
A direct computation gives 

$$
\begin{aligned}
\ip{I_N u, v}_N 
&= \frac{2\pi}{N} \sum_{j=0}^{N-1} I_N u(x_j) \conj{v(x_j)} \\
&= \frac{2\pi}{N} \sum_{j=0}^{N-1} u(x_j) \conj{v(x_j)} \\
&= \ip{u,v}_N
\end{aligned}
$$
:::

Hence the interpolant satisfies

$$
\ip{I_N u - u, v}_N = 0, \qquad \forall v \in S_N
$$ 

Thus 

> $I_N u$ is the orthogonal projection of $u$ onto $S_N$ with respect to the discrete inner product.
>  
> The error $u - I_N u$ is perpendicular to $S_N$ and $I_N u$ is the best approximation to $u$ in $S_N$ wrt the discrete norm $\norm{\cdot}_N$.

This is also related to the property that $\{ \phi_k(x) \}$ is an orthogonal family wrt to the discrete inner product

$$
\ip{\phi_j, \phi_k}_N = 2 \pi \delta_{jk}, \qquad 0 \le j,k \le N
$$

+++

:::{prf:lemma} Aliasing
:label: lem:falias
If the Fourier series $Su$ converges to $u$ at every node $x_j$, then

$$
\tilu_k = \hatu_k + \sum_{m=-\infty,m\ne 0}^\infty \hatu_{k + Nm}, \qquad k=-N/2, \ldots, N/2-1
$$
:::

:::{prf:proof}
Starting from the definition of $\tilu_k$ 

$$
\begin{aligned}
\tilu_k 
&= \frac{1}{N} \sum_{j=0}^{N-1} \clr{red}{u(x_j)} \ee^{-\ii k x_j} \\
&= \frac{1}{N} \sum_{j=0}^{N-1} \clr{red}{ \sum_{l=-\infty}^\infty \hatu_l \ee^{\ii l x_j} } \ee^{-\ii k x_j} \\
&= \sum_{l=-\infty}^\infty \hatu_l \left( \frac{1}{N} \sum_{j=0}^{N-1} \ee^{\ii (l- k)x_j} \right)
\end{aligned}
$$ 

But by [](#eq:discorth)

$$
\frac{1}{N} \sum_{j=0}^{N-1} \ee^{\ii (l-k)x_j}  = \begin{cases}
1, & \textrm{if } l-k = Nm, \quad m=0, \pm 1, \pm 2, \ldots \\
0, & \textrm{otherwise}
\end{cases}
$$ 

Hence 

$$
\tilu_k = \sum_{m=-\infty}^\infty \hatu_{k + Nm}
= \hatu_k + \sum_{m=-\infty,m\ne 0}^\infty \hatu_{k + Nm}
$$
:::

Define

$$
P_N u = \sumf \hatu_k \phi_k \in S_N
$$

which is a truncation of the Fourier series $Su$. Also define the summation notation

$$
\sum_{k=-\infty}^\infty = \sumf + \sumfr
$$

:::{prf:lemma}
Under the conditions of [](#lem:falias)

$$
\norm{u - P_N u}_2 \le \norm{u - I_N u}_2
$$
:::

:::{prf:proof}
Using [](#lem:falias), the interpolant can be written as

$$
\begin{aligned}
I_N u 
&= \sumf \tilu_k \phi_k \\
&= \clr{red}{\sumf \hatu_k \phi_k} + \clr{blue}{ \sumf \left( \sum_{m=-\infty,m\ne 0}^\infty \hatu_{k + Nm} \right) \phi_k } \\
&=: \clr{red}{P_N u} + \clr{blue}{R_N u}
\end{aligned}
$$ 

Hence 

$$
u - I_N u = (u - P_N u) - R_N u
$$

and

$$
u - P_N u = \sumfr \hatu_k \phi_k 
\qquad
\textrm{is orthogonal to}
\qquad
R_N u \in S_N
$$

we get

$$
\norm{u - I_N u}_2^2 = \norm{u - P_N u}_2^2 + \norm{R_N u}_2^2
$$

which proves the desired result.
:::

+++

The result of [](#lem:falias) and decay properties of Fourier transform are key to showing the error estimate of trigonometric interpolation.

:::{prf:lemma}
:label: lem:trigerr
Let $u \in \cts_p^{m-1}$ for $m \ge 1$ and $f^{(m)}$ is piecewise continuous. Then

\begin{align}
\norm{u - P_N u}_\infty \le \sumfr |\hatu_k| &= \order{\frac{1}{N^m}} \\
\norm{u - I_N u}_\infty \le 2 \sumfr |\hatu_k| &= \order{\frac{1}{N^m}}
\end{align}
:::

:::{prf:proof}
From the smoothness of the function, we have

$$
|\hatu_k| = \order{\frac{1}{|k|^{m+1}}}, \qquad m \ge 1
$$

The first inequality is easy since

$$
u - P_N u = \sumfr \hatu_k \phi_k \quad\Longrightarrow\quad \norm{u - P_N u}_\infty \le \sumfr |\hatu_k|
$$ 

For the second, we use the decomposition

$$
u - I_N u = u - P_N u - R_N u
$$ 

so that

$$
\norm{u - I_N u}_\infty \le \norm{u - P_N u}_\infty + \norm{R_N u}_\infty
$$

But

$$
\norm{R_N u}_\infty \le \sumf \sum_{m=-\infty,m\ne 0}^\infty |\hatu_{k+Nm}| =
\sumfr |\hatu_k|
$$ 

and hence

$$
\norm{u - I_N u}_\infty \le \sumfr |\hatu_k| + \sumfr |\hatu_k| = 2 \sumfr |\hatu_k|
$$

The error is given by

$$
\sumfr |\hatu_k| \le C \left[ \frac{1}{(N/2)^{m+1}} + \frac{1}{(N/2+1)^{m+1}} + \ldots \right]
$$

The sum on the right can be bounded by an integral

$$
\begin{aligned}
\frac{1}{(N/2)^{m+1}} + \frac{1}{(N/2+1)^{m+1}} + \ldots  
&< \int_0^\infty \frac{1}{(N/2 + x - 1)^{m+1}} \ud x \\
&= \frac{1}{m(N/2-1)^{m}} \\
&= \order{\frac{1}{N^{m}}} \quad \textrm{for } N \gg 1
\end{aligned}
$$
:::


:::{prf:remark}
The error of trigonometric interpolation is at most twice the error in the truncated Fourier series. This shows that trigonometric interpolation at uniformly spaced points gives very accurate approximations provided the function is sufficiently smooth.
:::

:::{prf:remark}
We can derive the error in $L^2$ norm, see Gander and Kwok, Theorem 4.3 and Theorem 4.7.
:::

## Computing the DFT

### DFT of a real function

If the function $u$ is real, then the DFT comes in complex conjugate pairs 

$$
\tilu_{-k} = \conj{\tilu_k}, \qquad 1 \le k \le N/2-1
$$ 

and $\tilu_0$, $\tilu_{-N/2}$ are real. The first property is easy. We have

$$
\tilu_0 = \frac{1}{N} \sum_{j=0}^{N-1} u(x_j)
$$ 

which is real, and

$$
\tilu_{-N/2} = \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \ee^{\ii \pi j} =  \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \cos(\pi j)
$$

is also real. When we want to evaluate $I_N u(x)$ at some arbitrary $x$, 

$$
I_N u(x) = \tilu_0 + \sum_{k=1}^{N/2-1} \left( \tilu_k \phi_k(x) + \conj{\tilu_k \phi_k(x)} \right) + \tilu_{-N/2} \phi_{-N/2}(x)
$$

the above expression is not real because of the last term. We can modify it as

\begin{align}
I_N u(x) &= \tilu_0 + \sum_{k=1}^{N/2-1} \left( \tilu_k \phi_k(x) +   \conj{\tilu_k \phi_k(x)} \right) \\
& \qquad + \half \tilu_{-N/2} \phi_{-N/2}(x) + \half \tilu_{-N/2} \phi_{N/2}(x)
\end{align}

i.e.,

$$
I_N u(x) = {\sum_{k=-N/2}^{N/2}}^\prime \hatu_k \phi_k(x), \qquad \hatu_{N/2} = \hatu_{-N/2}
$$

and the prime denotes that the first and last terms must be multiplied by half. This ensures that the interpolation conditions are still satisfied and $I_N u(x)$ is real for all $x \in [0,2\pi]$. 

However, in actual numerical computation, we may get an imaginary part of order machine precision, and then we have to take the real part. Then, we may as well compute it as before and take real part

$$
I_N u (x) = \textrm{Real} \sumf \tilu_k \phi_k(x)
$$

### Fast Fourier Transform (FFT)

Computing each $\tilu_k$ using [](#eq:tildeuk) requires $\order{N}$ floating point operations (FLOPS). Computing all $N$ coefficients thus requires $\order{N^2}$ FLOPS. This process can be written as a dense $N \times N$ matrix multiplying the vector of values

$$
[u(x_0), \ u(x_1), \ \ldots, \ u(x_{N-1})]^\top
$$

The matrix has special structure and using some tricks, the cost can be reduced to $\order{N \log N}$ using FFT. It was invented by Gauss and rediscovered 160 years later in [@Cooley1965].

### FFT using Numpy

Numpy has routines to compute the DFT using FFT, see [numpy.fft](https://numpy.org/doc/stable/reference/routines.fft.html). It takes the trigonometric polynomial in the form

$$
I_N u(x) = \frac{1}{N} \sum_{k=-N/2+1}^{N/2} \tilu_k \ee^{\ii k x}, \qquad x \in [0,2\pi]
$$ 

which is different from our convention.

```python
from numpy import pi,arange
import numpy.fft as fft
h = 2*pi/N
x = h*arange(0,N)
u = f(x)
u_hat = fft(u)
```

The output `u_hat` gives the DFT in the following order

$$
\tilu_k = \tilu_0, \ \tilu_1, \ \ldots, \ \tilu_{N/2}, \ \tilu_{-N/2+1}, \ \ldots, \ \tilu_{-2}, \ \tilu_{-1}
$$

[Scipy.fft](https://docs.scipy.org/doc/scipy/tutorial/fft.html) also provides equivalent routines for FFT.

See more examples of using Fourier interpolation for approximation and differentiation here <http://cpraveen.github.io/chebpy>. 

+++

## Examples

The following function computes the trigonometric interpolant and plots it. It also computes error norm by using `ne` uniformly spaced points.

```{code-cell}
# Wave numbers are arranged as k=[0,1,...,N/2,-N/2+1,-N/2,...,-1]
def fourier_interp(N,f,ne=1000,fig=True):
    if mod(N,2) != 0:
        print("N must be even")
        return
    h = 2*pi/N; x = h*arange(0,N);
    v = f(x);
    v_hat = fft(v)
    k = zeros(N)
    n = N//2
    k[0:n+1] = arange(0,n+1)
    k[n+1:] = arange(-n+1,0,1)

    xx = linspace(0.0,2*pi,ne)
    vf = real(exp(1j*outer(xx,k)) @ v_hat) / N
    ve = f(xx)

    # Plot interpolant and exact function
    if fig:
        plot(x,v,'o',xx,vf,xx,ve)
        legend(('Data','Fourier','Exact'))
        xlabel('x'), ylabel('y'), title("N = "+str(N))

    errori = abs(vf-ve).max()
    error2 = sqrt(h*sum((vf-ve)**2))
    print("Error (max,L2) = ",errori,error2)
    return errori, error2
```

To measure rate of error convergence, we can compute error with $N$ and $2N$ points

$$
E_N = \frac{C}{N^m}, \qquad E_{2N} = \frac{C}{(2N)^m}
$$

and estimate the rate 

$$
m = \frac{\log(E_N / E_{2N})}{\log 2}
$$

### Infinitely smooth, periodic function $(m = \infty)$

$$
u(x) = \exp(\sin x), \qquad x \in [0,2\pi]
$$

```{code-cell}
f1 = lambda x: exp(sin(x))
fourier_interp(24,f1);
```


### Infinitely smooth, derivative not periodic $(m = 1)$

$$
u(x) = \sin(x/2), \qquad x \in [0,2\pi]
$$

```{code-cell}
g = lambda x: sin(x/2)
e1i,e12 = fourier_interp(24,g)
```

```{code-cell}
e2i,e22 = fourier_interp(48,g,fig=False)
print("Rate  (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```

### Continuous function $(m = 1)$

$$
f(x) = \begin{cases}
1 - \frac{2}{\pi} |x - \pi|, & \pi/2 \le x \le 3\pi/2 \\
0, & \textrm{otherwise}
\end{cases}, \qquad x \in [0,2\pi]
$$

```{code-cell}
trihat = lambda x: (x >= 0.5*pi)*(x <= 1.5*pi)*(1 - 2*abs(x-pi)/pi)
e1i,e12 = fourier_interp(24,trihat)
```

```{code-cell}
e2i,e22 = fourier_interp(48,trihat,fig=False)
print("Rate  (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```

### Discontinuous function $(m=0)$

[](#lem:trigerr) does not apply in this case.

$$
f(x) = \begin{cases}
1, & \pi/2 < x < 3\pi/2 \\
0, & \textrm{otherwise}
\end{cases}, \qquad x \in [0,2\pi]
$$

```{code-cell}
f2 = lambda x: (abs(x-pi) < 0.5*pi)
e1i,e12 = fourier_interp(24,f2)
```

```{code-cell}
e2i,e22 = fourier_interp(48,f2)
print("Rate  (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```

We do not get convergence in maximum norm but we see convergence in 2-norm at a rate of $1/N$. There is convergence away from the discontinuities but overall, it is not a good approximation.

+++

:::{note}
Uniformly spaced points were not a good choice for polynomial interpolation, since they lead to Runge phenomenon. We needed to use points clustered at the boundaries to get good approximations.

Uniformly spaced points are the best choice for trigonometric interpolation. They lead to discrete orthogonality and an explicit expression for the DFT. Moreover, the highly efficient FFT is also possible because of this property.
:::

+++

:::{exercise}
With $N=100,500$, compute the FFT of all the functions in above examples and plot $|\tilu_k|$ versus $k$ for $0 \le k \le N/2$. Choose linear or log scale by deciding which one is most informative. Also plot it for the function

$$
u(x) = \sin(10x) + \sin(50x), \qquad x \in [0,2\pi]
$$

Comment on the nature of these curves.
:::

+++

:::{exercise}
To test convergence as $N \to \infty$, compute $I_N u$ on a uniform grid $\{ x_j \}$ of say $m = 1000$ points and measure the error as

\begin{align}
\norm{I_N u - u}_\infty &\approx \max_{j}  |I_N u(x_j) - u(x_j)| \\
\norm{I_N u - u}_2 &\approx \left( \frac{2\pi}{m} \sum_{j}  [I_N u(x_j) - u(x_j)]^2 \right)^\half
\end{align}

Plot error versus $N$ in loglog scale or semilogy scale; choose the scale appropriately based on smoothness of function. Do this for all the functions in the above examples.
:::

+++

:::{attention} Approximating derivatives
If $\hatu_k$ is the continuous Fourier transform of $u(x)$, then the transform of $u'(x)$ is $\ii k \hatu_k$. This suggests that using the DFT $\tilu_k$, we can take the inverse DFT of $\ii k \tilu_k$ to obtain an approximation to $u'(x)$. 

$$
\{ u_j \} \overset{DFT}{\longrightarrow} \{ \tilu_k \} \longrightarrow \{ \ii k \tilu_k \} \overset{IDFT}{\longrightarrow} \{ u_j' \}
$$

We can thus compute derivatives at all points with almost $\order{N}$ FLOPS. We will study more about this in future chapters.
:::

+++

:::{seealso}
Here are some nice applications of DFT from [@Brunton2019], also see the book website <https://databookuw.com>.

1. [De-noising a signal](https://github.com/dynamicslab/databook_python/blob/master/CH02/CH02_SEC02_2_Denoise.ipynb), [Video](https://www.youtube.com/watch?v=s2K1JfNR7Sc)
1. [Compressing an image](https://github.com/dynamicslab/databook_python/blob/master/CH02/CH02_SEC06_2_Compress.ipynb), [Video](https://www.youtube.com/watch?v=gGEBUdM0PVc), [Video](https://www.youtube.com/watch?v=uB3v6n8t2dQ)
1. [De-noising an image](https://github.com/dynamicslab/databook_python/blob/master/CH02/CH02_SEC06_3_Denoise.ipynb)
1. [More videos on FFT](https://databookuw.com/page-2/page-21/)
:::
