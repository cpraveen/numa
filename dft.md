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

# Trigonometric interpolation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## Interpolation problem

For any even integer $N>0$, consider the set of points

$$
x_j = \frac{2\pi j}{N}, \qquad j=0,1,\ldots,N-1
$$ 

Note that $x_0 = 0$ and the last point $x_{N-1} < 2\pi$; we do not include the point $x=2\pi$ in the grid since by periodicity, the value there is same as at $x=0$.

The degree $N/2$ trigonometric polynomial is

$$
I_N u(x) = \sumf \tilu_k \phi_k(x) = \sumf \tilu_k \ee^{\ii k x}
$$ 

We will determine the $N$ coefficients $\tilu_k$ from the interpolation conditions

$$
u(x_j) = I_N u(x_j), \qquad j=0,1,\ldots,N-1
$$ 

Using the orthogonality relation 

$$
\frac{1}{N} \sum_{j=0}^{N-1} \ee^{-\ii p x_j} = \begin{cases}
1, & p = Nm, \quad m=0, \pm 1, \pm 2, \ldots \\
0, & \textrm{otherwise}
\end{cases}
$$ 

we obtain the discrete Fourier coefficients

$$
\tilu_k = \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \ee^{-\ii k x_j}, \qquad k=-N/2, \ldots, N/2 - 1
$$ 

This is known as the **discrete Fourier transform** (DFT).

### Lagrange form
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
Why did we not include the term corresponding to $k=N/2$ ? This is because the exponentials

$$
\phi_{-N/2}(x) = \ee^{-\ii Nx/2}, \qquad \phi_{N/2}(x) = \ee^{\ii Nx/2}
$$

take the same values at all the nodes $x_j$ and hence the two functions are not independent.
:::

Let us define the set of trigonometric functions

$$
S_N = \left\{ v : [0,2\pi] \to \re, \quad v = \sumf a_k \phi_k \right\}
$$

:::{prf:lemma}
If $u \in S_N$, i.e., $u = \sumf c_k \phi_k$, then $I_N u = u$.
:::

:::{prf:proof}
Show that $\tilu_k = c_k$.
:::

Define the discrete inner product

$$
\ip{u,v}_N = \frac{2\pi}{N} \sum_{j=0}^{N-1} u(x_j) \conj{v(x_j)}
$$

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

Hence the interpolant satifies

$$
\ip{I_N u - u, v}_N = 0, \qquad \forall v \in S_N
$$ 

Thus $I_N u$ is the orthogonal projection of $u$ onto $S_N$ with respect to the discrete inner product.

:::{prf:lemma}
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
&= \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \ee^{-\ii k x_j} \\
&= \frac{1}{N} \sum_{j=0}^{N-1} \sum_{l=-\infty}^\infty \hatu_l \ee^{\ii l x_j}
\ee^{-\ii k x_j} \\
&= \sum_{l=-\infty}^\infty \hatu_l \left( \frac{1}{N} \sum_{j=0}^{N-1} \ee^{\ii (l- k)x_j} \right)
\end{aligned}
$$ 

But

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

:::{prf:lemma}
$$
\norm{u - P_N u}_2 \le \norm{u - I_N u}_2
$$
:::

:::{prf:proof}
Using previous Lemma, the interpolant can be written as

$$
\begin{aligned}
I_N u 
&= \sumf \tilu_k \phi_k \\
&= \sumf \hatu_k \phi_k + \sumf \left( \sum_{m=-\infty,m\ne 0}^\infty \hatu_{k +
Nm} \right) \phi_k \\
&= P_N u + R_N u
\end{aligned}
$$ 

Hence 

$$
u - I_N u = u - P_N u - R_N u
$$ 

Since

$$
u - P_N u = \sumfr \hatu_k \phi_k
$$ 

is orthogonal to $R_N u$, we get

$$
\norm{u - I_N u}_2^2 = \norm{u - P_N u}_2^2 + \norm{R_N u}_2^2
$$

which proves the desired result.
:::

:::{prf:remark}
If $u \in L^2(0,2\pi)$ or piecewise continuous, then $|\hatu_k| = \order{\frac{1}{|k|}}$ and

$$
\norm{u - P_N u}_2 \to 0, \qquad \norm{R_N u}_2 \to 0
$$ 

and hence $\norm{u - I_N u}_2 \to 0$ as $N \to \infty$. This is convergence in the mean, but we cannot prove point-wise convergence under the above assumption on $u$. This is because

$$
|u(x) - P_N u(x)| \le \sumfr |\hatu_k| \not\to 0
$$ 

For pointwise convergence, we require $\hatu_k$ to decay faster than $1/|k|$.
:::

:::{prf:lemma}
Let $u \in \cts_p^{m-1}$ for $m \ge 1$ and $f^{(m)}$ is piecewise continuous. Then

$$
\norm{u - P_N u}_\infty \le \sumfr |\hatu_k| &= \order{\frac{1}{N^m}} \\
\norm{u - I_N u}_\infty \le 2 \sumfr |\hatu_k| &= \order{\frac{1}{N^m}}
$$
:::

:::{prf:proof}
From the smoothness of the function, we have

$$
|\hatu_k| = \order{\frac{1}{|k|^{m+1}}}
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

If the function $u$ is real, then the DFT comes in complex conjugate pairs 

$$
\hatu_{-k} = \conj{\hatu_k}, \qquad 1 \le k \le N/2-1
$$ 

and $\hatu_0$, $\hatu_{-N/2}$ are real. The first property is easy. We have

$$
\hatu_0 = \frac{1}{N} \sum_{j=0}^{N-1} u(x_j)
$$ 

which is real, and

$$
\hatu_{-N/2} = \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \ee^{\ii \pi j} =  \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) \cos(\pi j)
$$

is also real. When we want to evaluate $I_N u(x)$ at some arbitrary $x$, we should modify it as

$$
I_N u(x) = {\sum_{k=-N/2}^{N/2}}^\prime \hatu_k \ee^{\ii k x}
$$

where $\hatu_{N/2} = \hatu_{-N/2}$ and the prime denotes that the first and last terms must be multiplied by half. This ensures that the interpolation conditions are satisfied and $I_N u(x)$ is real for all $x \in [0,2\pi]$.

Numpy has routines to compute the DFT using FFT. It takes the trigonometric polynomial in the form

$$
I_N u(x) = \sum_{k=-N/2+1}^{N/2} \hatu_k \ee^{\ii k x}
$$ 

which is different from our convention.

```python
from numpy import pi,arange
import numpy.fft as fft
h = 2*pi/N; x = h*arange(0,N);
u = f(x);
u_hat = fft(v)
```

The output `u_hat` gives the DFT in the following order

$$
\hatu_0, \ \hatu_1, \ \ldots, \ \hatu_{N/2}, \ \hatu_{-N/2+1}, \ \ldots, \ \hatu_{-2}, \ \hatu_{-1}
$$

See more examples here <http://cpraveen.github.io/chebpy>.

## Examples

The following function computes the trigonometric interpolant and plots it.

```{code-cell}
# Wave numbers are arranged as k=[0,1,...,N/2,-N/2+1,-N/2,...,-1]
def fourier_interp(N,f,ne=500,fig=True):
    if mod(N,2) != 0:
        print("N must be even")
        return
    h = 2*pi/N; x = h*arange(0,N);
    v = f(x);
    v_hat = fft(v)
    k = zeros(N)
    n = int(N/2)
    k[0:n+1] = arange(0,n+1)
    k[n+1:] = arange(-n+1,0,1)

    xx = linspace(0.0,2*pi,ne)
    vf = real(dot(exp(1j*outer(xx,k)), v_hat)/N)
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


### Infintely smooth, periodic function


```{code-cell}
f1 = lambda x: exp(sin(x))
fourier_interp(24,f1);
```


### Infinitely smooth, derivative not periodic

```{code-cell}
g = lambda x: sin(x/2)
e1i,e12 = fourier_interp(24,g)
```

```{code-cell}
e2i,e22 = fourier_interp(48,g,fig=False)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```


### Continuous function

```{code-cell}
def trihat(x):
    f = 0*x
    for i in range(len(x)):
        if x[i] < 0.5*pi or x[i] > 1.5*pi:
            f[i] = 0.0
        elif x[i] >= 0.5*pi and x[i] <= pi:
            f[i] = 2*x[i]/pi - 1
        else:
            f[i] = 3 - 2*x[i]/pi
    return f

e1i,e12 = fourier_interp(24,trihat)
```

```{code-cell}
e2i,e22 = fourier_interp(48,trihat,fig=False)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```

### Discontinuous function

```{code-cell}
f2 = lambda x: (abs(x-pi) < 0.5*pi)
e1i,e12 = fourier_interp(24,f2)
```

```{code-cell}
e2i,e22 = fourier_interp(48,f2)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```
