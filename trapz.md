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

```{code-cell}
:tags: remove-input
a, b = 0.0, 1.0
f = lambda x: 1.0 + x + 0.5* sin(pi*x)
x = linspace(a,b,100)
plot([a,a],[0,f(a)],'k--')
plot([b,b],[0,f(b)],'k--')
plot([a,b],[0,0],'bo-')
plot(x, f(x), 'r-', label='f(x)')
plot([a,b],[f(a),f(b)],'sk-')
d = 0.05
text(a,-d,"a",ha='center',va='top')
text(b,-d,"b",ha='center',va='top')
axis([-0.1, 1.1, -0.3, 2.2])
legend(), xlabel('x');
```

Note that this is the area under the straight line curve $p_1(x)$ which has the shape of a (rotated) trapezoid. To study the error in this approximation, we start with the error in interpolation [](#eq:newtinterr)

$$
f(x) - p_1(x) = (x-a)(x-b) f[a,b,x]
$$ 

Then, the error in the integral is 

$$
\begin{aligned}
E_1(f) 
&= \int_a^b [f(x) - p_1(x)] \ud x \\
&= \int_a^b (x-a) (x-b) f[a,b,x] \ud x
\end{aligned}
$$ 

:::{prf:theorem} Integral mean value theorem
If $g$ does not change sign and $f$ is continuous, then 

$$
\int_a^b f(x) g(x)\ud x = f(c) \int_a^b g(x) \ud x
$$ 

for some $c \in [a,b]$.
:::

Using integral mean value theorem

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

```{code-cell}
:tags: remove-input
figure(figsize=(8,5))
xp = linspace(0.0,1.0,10)
x  = linspace(0,1,100)
f  = lambda x: 2 + sin(2*pi*x)*cos(2*pi*x)
plot([0,1],[0,0],'k+-')
plot(xp,0*xp,'o')
plot(x,f(x),'r-',label='f(x)')
plot(xp,f(xp),'sk-')
for xx in xp:
    plot([xx,xx],[0,f(xx)],'b--')
axis([-0.1,1.1,-0.3,2.7])
d = 0.1
text(xp[0],-d,'$x_0$',va='top',ha='center')
text(xp[1],-d,'$x_1$',va='top',ha='center')
text(xp[-2],-d,'$x_{n-1}$',va='top',ha='center')
text(xp[-1],-d,'$x_n$',va='top',ha='center')
text(0,d,'$a$',ha='center',backgroundcolor='white')
text(1,d,'$b$',ha='center',backgroundcolor='white')
legend(), xlabel('x');
```

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

If $n$ is doubled, the error is decreased by a factor of $\frac{1}{4}$. We should observe that

$$
\frac{E_n(f)}{E_{2n}(f)} = \frac{1/n^2}{1/(2n)^2} = 2^2 = 4
$$
:::

:::{prf:example}
The next functions computes integral with trapezoidal rule.

```{code-cell}
# n intervals, n+1 points
def trapz1(a,b,n,f):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res = h * (sum(y[1:n]) + 0.5*(y[0] + y[n]))
    return res
```

\begin{align}
f(x) &= x, \qquad x \in [0,1] \\
\int_0^1 f(x) \ud x &= \half
\end{align}

```{code-cell}
f = lambda x: x
print("Integral = ", trapz1(0.0,1.0,10,f))
```

\begin{align}
f(x) &= x^2, \qquad x \in [0,1] \\
\int_0^1 f(x) \ud x &= \frac{1}{3}
\end{align}

```{code-cell}
f = lambda x: x**2
print("Integral = ", trapz1(0.0,1.0,10,f))
```

\begin{align}
f(x) &= \exp(x)\cos(x), \qquad x \in [0,\pi] \\
\int_0^1 f(x) \ud x &= -\frac{1}{2}(1+\exp(\pi))
\end{align}

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

The last column shows $E_n(f)/E_{2n}(f)$, which is close to $4$, showing that the error is $\order{h^2}$.
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
E_n(f) \approx -\frac{h^2}{12}[f'(b) - f'(a)] =: \tilde E_n(f)
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

+++

## Corrected trapezoidal rule

Using the asymptotic error estimate

\begin{gather}
E_n(f) = I(f) - I_n(f) \approx \tilde E_n(f) = -\frac{h^2}{12}[f'(b) - f'(a)] \\
\implies I(f) \approx I_n(f) -\frac{h^2}{12}[f'(b) - f'(a)]
\end{gather}

we can define a **corrected Trapezoidal rule** as

\begin{align}
C_n(f) &= I_n(f) -\frac{h^2}{12}[f'(b) - f'(a)] \\
&= h \left[ \shalf f_0 + f_1 + f_2 + \ldots + f_{n-1} + \shalf f_n \right] - \frac{h^2}{12}[f'(b) - f'(a)] 
\end{align}

+++

:::{prf:example}
Consider the integral 

$$
I = \int_0^\pi \ee^x \cos x \ud x
$$ 

The exact value is 

$$
I = -\half (1 + \ee^\pi) \approx -12.0703463164
$$ 

The next function implements both trapezoid and the corrected rule.

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

As seen in last column, the error in the corrected rule converges at a rate of $16$, i.e., as $\order{h^4}$.
:::
