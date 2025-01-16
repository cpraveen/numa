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

# Piecewise polynomial interpolation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from matplotlib import animation
from IPython.display import HTML
```

Consider a partition of $[a,b]$ into a grid

$$
a = x_0 < x_1 < \ldots < x_N = b
$$ 

The points $x_j$ are called knots, breakpoints or nodes. Let $p(x)$ be a polynomial on each of the subintervals 

$$
[x_0,x_1], \quad [x_1,x_2], \quad \ldots, [x_{N-1},x_N]
$$

i.e., 

$$
p(x) \in \poly_r, \qquad x \in [x_i,x_{i+1}]
$$ 

Note that $p(x)$ need not be a polynomial in $[a, b]$. We say that $p(x)$ is a piecewise polynomial of degree $r$ if the degree of $p(x)$ is $\le r$ on each of the sub-intervals. In general, no restrictions on the continuity of $p(x)$ or its derivatives at the end points of the sub-intervals might exist.

## Piecewise linear interpolation

Let us demand the we want a continuous approximation. We are given the function values $y_i = f(x_i)$ at all the nodes of our grid. We can construct the polynomial of each sub-interval

$$
p(x) = \frac{x - x_{i+1}}{x_i - x_{i+1}} y_i + \frac{x - x_i}{x_{i+1} - x_i} y_{i+1},
\qquad x_i \le x \le x_{i+1}
$$ 

Clearly, this is linear and continuous.  On each interval we know the interpolation error estimate

$$
\max_{x \in [x_i,x_{i+1}]} |f(x) - p(x)| \le \frac{h_{i+1}^2}{2} \max_{t \in
[x_i,x_{i+1}]} |f''(t)|, \qquad h_{i+1} = x_{i+1} - x_i
$$ 

from which we get the global error

$$
\max_{x \in [x_0,x_N]} |f(x) - p(x)| \le \frac{h^2}{2} \max_{t \in [x_0,x_N]} |f''(t)|,
\qquad h = \max_i h_i
$$ 

Clearly, we have convergence as $N \to \infty$, $h \to 0$. The convergence is quadratic; if $h$ becomes $h/2$, the error reduces by a factor of $1/4$.

:::{prf:remark}
The polynomial degree is fixed and convergence is achieved because the spacing between the grid points becomes smaller. Note that we only require some condition on the second derivative of the function $f$ while in the global interpolation problem, we required conditions on higher derivatives, which can be difficult to satisfy in some situations.
:::

### Finite element form

We have written down the polynomial in a piecewise manner but we can also write a global view of the polynomial by defining some basis functions 

$$
\phi_i(x) = \begin{cases}
\frac{x - x_{i-1}}{x_i - x_{i-1}} & x_{i-1} \le x \le x_i \\
\frac{x_{i+1} - x}{x_{i+1} - x_i} & x_i \le x \le x_{i+1} \\
0 & \textrm{otherwise}
\end{cases}
$$ 

These are triangular hat functions with compact support.  Then the piecewise linear approximation is given by

$$
p(x) = \sum_{i=0}^N y_i \phi_i(x)
$$

### An adaptive algorithm

We can use the local error estimate to improve the approximation. We have to first estimate the error in the an interval $I_i = [x_{i-1},x_i]$ by some finite difference formula $H_i$. If we want the error to be $\le \epsilon$, then we must choose $h_i$ so that

$$
\frac{h_i^2}{8} H_i \le \epsilon
$$ 

If in some interval, it happens that 

$$
\frac{h_i^2}{8} H_i > \epsilon
$$ 

divide $I_i$ into two intervals

$$
\left[ x_{i-1} \half(x_{i-1}+x_i) \right] \qquad \left[ \half(x_{i-1}+x_i), x_i \right]
$$ 

We can start with a uniform grid of $N$ intervals and apply the above algorithm. The code finds the interval with the maximum error and divides it into two intervals. The error in the interval $[x_{i-1},x_i]$ is estimated by fitting a cubic polynomial to the data $x_{i-2},x_{i-1},x_i,x_{i+1}$ using `polyfit` and finding the maximum value of its second derivative at $x_{i-1},x_i$ to estimate $H_i$.  `P = polyfit(x,f,3)` returns a polynomial in the form

$$
P[0]*x^3 + P[1]*x^2 + P[2]*x + P[3]
$$ 

whose second derivative is

$$
6*P[0]*x + 2*P[1]
$$

:::{exercise}
Modify the program by changing the error estimate as follows. Estimate second derivative at $x_{i-1}$ by interpolating a quadratic to the data $x_{i-2},x_{i-1},x_i$ and taking its second derivative. Similarly, at $x_i$, find a quadratic polynomial to the data $x_{i-1},x_i,x_{i+1}$ and find its second derivative. Then estimate the maximum value of second derivative in the interval as the maximum of the derivatives at the end
points.
:::

:::{prf:example}

Let us try to approximate

$$
f(x) = \exp(-100(x-1/2)^2) \sin(4 \pi x), \qquad x \in [0,1]
$$

with piecewise linear approximation.

```{code-cell}
xmin, xmax = 0.0, 1.0
fun = lambda x: exp(-100*(x-0.5)**2) * sin(4*pi*x)
```

Here is the initial approximation.

```{code-cell}
N = 10 # number of initial points
x = linspace(xmin,xmax,N)
f = fun(x)

ne = 100
xe = linspace(xmin,xmax,ne)
fe = fun(xe)

fig, ax = subplots()
line1, = ax.plot(xe,fe,'-',linewidth=2)
line2, = ax.plot(x,f,'or--',linewidth=2)
ax.set_title('Initial approximation, N = ' + str(N));
```

The next function performs uniform and adaptive refinement.

```{code-cell}
def adapt(x,f,nadapt=100,err=1.0e-2,mode=1):
    N = len(x)

    for n in range(nadapt):
        h = x[1:] - x[0:-1]
        H = zeros(N-1)
        for j in range(1,N-2): # skip first and last element
           P = polyfit(x[j-1:j+3], f[j-1:j+3], 3)
           H[j] = max(abs(6*P[0]*x[j:j+2] + 2*P[1]))

        elem_err = (1.0/8.0) * h**2 * abs(H)
        i = argmax(elem_err); current_err = elem_err[i]
        if current_err < err:
           print('Satisfied error tolerance, N = ' + str(N))
           return x,f
        if mode == 0:
           x = linspace(xmin,xmax,2*N)
           f = fun(x)
        else:
           x = concatenate([x[0:i+1], [0.5*(x[i]+x[i+1])], x[i+1:]])
           f = concatenate([f[0:i+1], [fun(x[i+1])], f[i+1:]])
        N = len(x)
    return x,f
```

Let us first try uniform refinement

```{code-cell}
adapt(x, f, mode=0);
```

and adaptive refinement.

```{code-cell}
adapt(x, f, mode=1);
```

Here is an animation of adaptive refinement.

```{code-cell}
def fplot(frame_number,x,f):
    x1, f1 = adapt(x,f,frame_number)
    line2.set_data(x1,f1)
    ax.set_title('N = '+str(len(x1)))
    return line2,

anim = animation.FuncAnimation(fig, fplot, frames=22, fargs=[x,f], repeat=False)
HTML(anim.to_jshtml())
```

The adaptive algorithm puts new points in regions of large gradient, where the resolution is not sufficient, and does not add anything in other regions.
:::

## Piecewise quadratic interpolation

Consider a grid 

$$
x_0 < x_1 < \ldots < x_N
$$ 

and define the mid-points of the intervals

$$
x_\imh = \half(x_{i-1} + x_i), \qquad i=1,2,\ldots,N
$$ 

We are given the information 

$$
y_i = f(x_i), \qquad i=0,1,\ldots,N
$$

$$
y_\imh = f(x_\imh), \qquad i=1,2,\ldots,N
$$ 

In the interval $I_i = [x_{i-1},x_i]$ we know three values 

$$
y_{i-1}, y_\imh, y_i
$$ 

and we can determine a quadratic polynomial the interpolates these values

$$
p(x) = y_{i-1} \phi_0\left(\frac{x-x_{i-1}}{h_i}\right) + y_\imh \phi_1\left(\frac{x-
x_{i-1}}{h_i}\right) + y_i \phi_2\left(\frac{x-x_{i-1}}{h_i}\right)
$$

where $\phi_0,\phi_1,\phi_2$ are Lagrange polynomials given by

$$
\phi_0(\xi) = 2(\xi - \shalf)(\xi - 1), \qquad \phi_1(\xi) = 4\xi(\xi-1), \qquad \phi_2(\xi)
= 2\xi(\xi-\shalf)
$$ 

Note that the quadratic interpolation is continuous but may not be differentiable.

## Piecewise degree $k$ interpolation

Consider a grid 

$$
a = x_0 < x_1 < \ldots < x_N = b
$$ 

In each interval $[x_{i-1},x_i]$ we want to construct a degree $k$ polynomial approximation, so we need to know $k+1$ function values at $k+1$ distinct points 

$$
x_{i-1} = x_{i,0} < x_{i,1} < \ldots < x_{i,k} = x_i
$$

Then the degree $k$ polynomial is given by

$$
p(x) = \sum_{j=0}^k f(x_{i,j}) \phi_j\left( \frac{x-x_{i-1}}{h_i} \right)
$$

where $\phi_j : [0,1] \to \re$ are Lagrange polynomials with the property

$$
\phi_j(\xi_l) = \delta_{jl}, \quad 0 \le j,l \le k, \qquad \textrm{where} \quad \xi_l =
\frac{x_{i,l} - x_{i-1}}{h_i}
$$ 

Since we chose the internal nodes such that $x_{i,0} = x_{i-1}$ and $x_{i,k} = x_i$, the piecewise polynomial is continuous, but may not be differentiable.

## Error estimate in Sobolev spaces

The estimate we have derived required functions to be differentiable. Here we use weaker smoothness properties.

:::{prf:theorem}
Let $f \in H^2(0,1)$ and let $p(x)$ be the piecewise linear interpolant of $f$. Then 

$$
\norm{f - p}_{L^2(0,1)} \le C h^2 \norm{f''}_{L^2(0,1)}
$$

$$
\norm{f' - p'}_{L^2(0,1)} \le C h \norm{f''}_{L^2(0,1)}
$$
:::

:::{prf:proof}
Since $\cts^\infty(0,1)$ is dense in $H^2(0,1)$, it is enough to prove the error estimate for functions $f \in \cts^\infty(0,1)$.

\(1\) If $x \in [x_j, x_{j+1}]$ then 

$$
\begin{aligned}
f(x) - p(x) 
&= f(x) - \left[ f(x_j) + \frac{f(x_{j+1}) - f(x_j)}{x_{j+1} - x_j} (x - x_j)
\right] \\
&= \int_{x_j}^x f'(t) \ud t - \frac{x - x_j}{x_{j+1} - x_j} \int_{x_j}^{x_{j+1}} f'(t)
\ud t \\
&= (x-x_j) f'(x_j + \theta_x) - (x-x_j) f'(x_j + \theta_j), \quad {0 \le \theta_x \le x - x_j
\atop 0 \le \theta_j \le h_j} \\
&= (x-x_j) \int_{x_j + \theta_j}^{x_j + \theta_x} f''(t) \ud t
\end{aligned}
$$ 

Squaring both sides and using $|x-x_j| \le h_j$

$$
\begin{aligned}
|f(x) - p(x)|^2 
&\le h_j^2 \left[ \int_{x_j + \theta_j}^{x_j + \theta_x} f''(t) \ud t
\right]^2 \\
&\le h_j^2 \left[ \int_{x_j}^{x_{j+1}} |f''(t)| \ud t \right]^2 \\
&\le h_j^3 \int_{x_j}^{x_{j+1}} |f''(t)|^2 \ud t \quad \textrm{by Cauchy-Schwarz
inequality}
\end{aligned}
$$ 

Integrating on $[x_j,x_{j+1}]$

$$
\int_{x_j}^{x_{j+1}} |f(x) - p(x)|^2 \ud x \le h_j^4 \int_{x_j}^{x_{j+1}} |f''(t)|^2
\ud t
$$ 

Adding such an inequality from all the intervals

$$
\begin{aligned}
\norm{f-p}_{L^2(0,1)}^2 
&\le \sum_j h_j^4 \int_{x_j}^{x_{j+1}} |f''(t)|^2 \ud t \\
&\le h^4  \sum_j \int_{x_j}^{x_{j+1}} |f''(t)|^2 \ud t \\
&= h^4 \norm{f''}_{L^2(0,1)}^2
\end{aligned}
$$ 

we get the first error estimate.

\(2\) The error in derivative is 

$$
\begin{aligned}
f'(x) - p'(x) 
&= f'(x) - \frac{f(x_{j+1})  - f(x_j)}{h_j} \\
&= \frac{1}{h_j} \int_{x_j}^{x_{j+1}}[f'(x) - f'(t)]\ud t \\
&= \frac{1}{h_j} \int_{x_j}^{x_{j+1}} \int_t^x f''(y) \ud y \ud t
\end{aligned}
$$ 

Squaring both sides 

$$
\begin{aligned}
|f'(x) - p'(x)|^2 
&= \frac{1}{h_j^2} \left( \int_{x_j}^{x_{j+1}} \int_t^x f''(y) \ud y
\ud t \right)^2 \\
&\le \frac{1}{h_j^2} \left( \int_{x_j}^{x_{j+1}} \int_t^x |f''(y)| \ud y \ud t \right)^2
\\
&\le \frac{1}{h_j^2} \left( \int_{x_j}^{x_{j+1}} \int_{x_j}^{x_{j+1}} |f''(y)| \ud y
\ud t \right)^2 \\
&= \left( \int_{x_j}^{x_{j+1}} |f''(y)| \ud y \right)^2 \\
&\le h_j  \int_{x_j}^{x_{j+1}} |f''(y)|^2 \ud y \quad \textrm{by Cauchy-Schwarz
inequality}
\end{aligned}
$$ 

Integrating with respect to $x$

$$
\int_{x_j}^{x_{j+1}} |f'(x) - p'(x)|^2 \ud x \le h_j^2 \int_{x_j}^{x_{j+1}} |f''(y)|^2
\ud y
$$ 

Adding such inequalities from each interval, we obtain the second error estimate.
:::

