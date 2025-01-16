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

# Trapezoidal Rule

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## Single interval

Consider a function $f : [a,b] \to \re$ and let us approximate this by a polynomial of degree one, say $p_1(x)$. This can be obtained by interpolation 

$$
p_1(a) = f(a), \qquad p_1(b) = f(b)
$$ 

Then the approximate integral is

$$
I_1(f) = \int_a^b p_1(x) \ud x = \frac{b-a}{2}[f(a) + f(b)]
$$ 

Note that this is the area under the straight line curve $p_1(x)$ which has the shape of a trapezoid. To study the error in this approximation, we start with the error in interpolation

$$
f(x) - p_1(x) = (x-a)(x-b) f[a,b,x]
$$ 

Then 

$$
\begin{aligned}
E_1(f) 
&= \int_a^b [f(x) - p_1(x)] \ud x \\
&= \int_a^b (x-a) (x-b) f[a,b,x] \ud x
\end{aligned}
$$ 

Using integral mean value theorem[^1]

$$
E_1(f) = f[a,b,\xi] \int_a^b (x-a)(x-b) \ud x, \qquad \textrm{for some $\xi \in [a,b]$}
$$

Finally, using the properties of divided differences 

$$
\begin{aligned}
E_1(f) 
&= \left[ \half f''(\eta) \right] \left[ -\frac{1}{6} (b-a)^3 \right], \qquad \textrm{for some $\eta \in [a,b]$} \\
&= -\frac{1}{12} (b-a)^3 f''(\eta)
\end{aligned}
$$ 

The error will be small only if $b-a$ is small.

## Composite trapezoidal rule

Divide $[a,b]$ into $n$ equal intervals by the partition

$$
a = x_0 < x_1 < \ldots < x_n = b \qquad \textrm{with spacing} \qquad h = \frac{b-a}{n} 
$$

and 

$$
x_j = a + j h, \qquad j=0,1,2,\ldots, n
$$ 

Let us use the Trapezoidal rule in each interval $[x_{j-1},x_j]$

$$
I(f) = \int_a^b f(x) \ud x = \sum_{j=1}^n \int_{x_{j-1}}^{x_j} f(x) \ud x
$$

From the previous section, we know that

$$
\int_{x_{j-1}}^{x_j} f(x) \ud x = \frac{h}{2}[f(x_{j-1}) + f(x_j)] - \frac{h^3}{12}f''(\eta_j), \qquad \eta_j \in [x_{j-1}, x_j]
$$

Hence 

$$
\begin{aligned}
I(f) 
&= \sum_{j=1}^n \frac{h}{2}[f(x_{j-1}) + f(x_j)] - \sum_{j=1}^n \frac{h^3}{12}f''(\eta_j) \\
&= I_n(f) - \sum_{j=1}^n \frac{h^3}{12}f''(\eta_j)
\end{aligned}
$$ 

where

$$
I_n(f) = h [ \shalf f_0 + f_1 + f_2 + \ldots + f_{n-1} + \shalf f_n]
$$

is the composite Trapezoidal rule. The error is

$$
E_n(f) = I(f) - I_n(f) = - \frac{h^3}{12} \sum_{j=1}^n f''(\eta_j) = - \frac{n h^3}{12} \underbrace{\frac{1}{n} \sum_{j=1}^n f''(\eta_j)}_{M_n}
$$

where

$$
\min_{x \in [a,b]} f''(x) \le \min_{1 \le j \le n} f''(\eta_j) \le M_n \le \max_{1 \le j \le n} f''(\eta_j) \le \max_{x \in [a,b]} f''(x)
$$

If $f''(x)$ is continuous, then it attains all values between the minimum and maximum.

$$
\implies \exists \eta \in [a,b] \qquad \textrm{such that} \qquad M_n = f''(\eta)
$$

Hence, after using $nh=b-a$, the error is given by

$$
E_n(f) = - \frac{1}{12}(b-a) h^2 f''(\eta), \qquad \textrm{for some $\eta \in [a,b]$}
$$

The error goes to zero as $n \to \infty$ since $h \to 0$.

:::{prf:remark}
The above Trapezoidal rule can be interpreted as the exact integral of the piecewise linear interpolant of $f(x)$ at uniformly spaced points.  However, it is not necessary to use uniformly spaced points and we can apply the above ideas for arbitrarily spaced points.
:::

:::{prf:remark}
If $f(x)$ is a linear polynomial, then $f''=0$ and the Trapezoidal rule is exact for such functions. This is also clear from the way it was constructed, as integral of piecewise linear approximations.
:::

:::{prf:remark}
The error of Trapezoidal rule converges quadratically

$$
|E_n(f)| = \order{h^2} = \order{\frac{1}{n^2}}
$$ 

If $n$ is doubled, the error decreased by a factor of $\frac{1}{4}$.
:::

:::{prf:example}
Show the Python notebook for integral of $f(x)=x$ and $f(x) = x^2$.

```{code-cell}
# n intervals, n+1 points
def trapz1(a,b,n,f):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res = h * (sum(y[1:n]) + 0.5*(y[0] + y[n]))
    return res
```

$$
f(x) = x, \qquad x \in [0,1]
$$

Exact integral is 0.5

```{code-cell}
f = lambda x: x
print("Integral = ", trapz1(0.0,1.0,10,f))
```

$$
f(x) = x^2, \qquad x \in [0,1]
$$

Exact integral is $1/3$

```{code-cell}
f = lambda x: x**2
print("Integral = ", trapz1(0.0,1.0,10,f))
```

$$
f(x) = \exp(x)\cos(x), \qquad x \in [0,\pi]
$$

The exact integral is $-\frac{1}{2}(1+\exp(\pi))$.

```{code-cell}
f = lambda x: exp(x)*cos(x)
qe = -0.5*(1.0 + exp(pi)) # Exact integral

n,N = 4,10
e = zeros(N)
for i in range(N):
    e[i] = trapz1(0.0,pi,n,f) - qe
    if i > 0:
        print('%6d %24.14e %14.5e'%(n,e[i],e[i-1]/e[i]))
    else:
        print('%6d %24.14e'%(n,e[i]))
    n = 2*n
```

:::

:::{prf:remark}
The [scipy.integrate](https://docs.scipy.org/doc/scipy/reference/integrate.html)  module provides functions that perform the [Trapezoidal integration](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.trapezoid.html).

```{code-cell}
from scipy.integrate import trapezoid
f = lambda x: x**2
a,b,n = 0.0,1.0,11
x = linspace(a,b,n)
y = f(x)
print("Integral = ", trapezoid(y,x))
```

For uniformly space points we can just pass the spacing instead of x.

```{code-cell}
h = (b-a)/(n-1)
print("Integral = ", trapezoid(y,dx=h))
```
:::

## Asymptotic error estimate

We can estimate the error as

$$
|E_n(f)| \le \frac{1}{12}(b-a) h^2 \max_{x \in [a,b]}|f''(x)|
$$

provided we can estimate the second derivative. A more useful error estimate can be derived as follows.

Using the error of the composite rule, we can compute the limit

$$
\begin{aligned}
\lim_{n \to \infty} \frac{E_n(f)}{h^2} 
&= - \frac{1}{12} \lim_{n \to \infty} h \sum_{j=1}^n f''(\eta_j), \qquad \eta_j \in [x_{j-1}, x_j] \\
&= -\frac{1}{12} \int_a^b f''(x) \ud x \\
&= -\frac{1}{12}[f'(b) - f'(a)]
\end{aligned}
$$ 

Hence for large $n$

$$
E_n(f) \approx \tilde E_n(f) := -\frac{h^2}{12}[f'(b) - f'(a)]
$$

:::{prf:definition}
Let $E_n(f)$ be an exact error formula and let $\tilde E_n(f)$ be an estimate of the true error. We say that $\tilde E_n(f)$ is an asymptotic error estimate if 

$$
\lim_{n \to \infty} \frac{\tilde E_n(f)}{E_n(f)} = 1
$$

or equivalently

$$
\lim_{n \to \infty} \frac{E_n(f) - \tilde E_n(f)}{E_n(f)} = 0
$$
:::

:::{prf:example}
Consider the integral 

$$
I = \int_0^\pi \ee^x \cos x \ud x
$$ 

The exact value is 

$$
I = -\half (1 + \ee^\pi) \approx -12.0703463164
$$ 

The error decreases at a rate of $4$. Using the asymptotic error estimate

$$
E_n(f) = I(f) - I_n(f) \approx \tilde E_n(f) = -\frac{h^2}{12}[f'(b) - f'(a)]
$$

$$
\implies I(f) \approx I_n(f) -\frac{h^2}{12}[f'(b) - f'(a)]
$$ 

we can define a *corrected Trapezoidal rule* as

$$
C_n(f) = h [ \shalf f_0 + f_1 + f_2 + \ldots + f_{n-1} + \shalf f_n] - \frac{h^2}{12}[f'(b) - f'(a)]
$$

The error in the corrected rule converges at a rate of $16$.

```{code-cell}
# n = number of intervals
# There are n+1 points
def trapz2(a,b,n,f,df):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res = sum(y[1:n]) + 0.5*(y[0] + y[n])
    return h*res, h*res - (h**2/12)*(df(b) - df(a))
```

```{code-cell}
f  = lambda x: exp(x)*cos(x)
df = lambda x: exp(x)*(cos(x) - sin(x))
qe = -0.5*(1.0 + exp(pi)) # Exact integral

n,N = 4,10
e1,e2 = zeros(N),zeros(N)
for i in range(N):
    e1[i],e2[i] = trapz2(0.0,pi,n,f,df) - qe
    if i > 0:
        print('%6d %24.14e %10.5f %24.14e %10.5f' %
              (n,e1[i],e1[i-1]/e1[i],e2[i],e2[i-1]/e2[i]))
    else:
        print('%6d %24.14e %10.5f %24.14e %10.5f' %
              (n,e1[i],0,e2[i],0))
    n = 2*n
```
:::

[^1]: If $g$ does not change sign and $f$ is continuous, then
    $\int_a^b f(x) g(x)\ud x = f(c) \int_a^b g(x) \ud x$ for some
    $c \in [a,b]$.
