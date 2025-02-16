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

# Simpson rule

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from scipy.interpolate import barycentric_interpolate
```

## Single interval

Trapezoidal method used a linear approximation to estimate the integral.
To improve the accuracy, let us use a quadratic approximation. For the
function $f : [a,b] \to \re$, we first construct a quadratic polynomial
$p_2(x)$ that interpolates at $a, b, c = \half(a+b)$

$$
p_2(a) = f(a), \qquad p_2(b) = f(b), \qquad p_2(c) = f(c)
$$ 

and the
approximation to the integral is

$$
I_2(f) = \int_a^b p_2(x) \ud x = \frac{h}{3}[f(a) + 4 f(c) + f(b)], \qquad h = \half(b-a)
$$

This is called Simpson rule. 


```{code-cell}
:tags: remove-input
a, b = 0.0, 1.0
c    = 0.5*(a + b)
f = lambda x: 1.0 + x + 0.5 * sin(pi*x)
x = linspace(a,b,100)
plot([a,a],[0,f(a)],'k--')
plot([b,b],[0,f(b)],'k--')
plot([c,c],[0,f(c)],'k--')
plot([a,c,b],[0,0,0],'bo-')
plot([a,c,b],[f(a),f(c),f(b)],'sk')
plot(x, f(x), 'r-', label='f(x)')
plot(x,barycentric_interpolate([a,c,b],[f(a),f(c),f(b)],x),'k-',label='Quadratic')
d = 0.05
text(a,-d,"a",ha='center',va='top')
text(b,-d,"b",ha='center',va='top')
text(c,-d,"c",ha='center',va='top')
axis([-0.1, 1.1, -0.3, 2.2])
legend(), xlabel('x');
```

The error in the integral follows from polynomial interpolation error [](#eq:newtinterr)

\begin{align}
E_2(f) &= I(f) - I_2(f) \\
&= \int_a^b [f(x) - p_2(x)] \ud x \\
&= \int_a^b \underbrace{(x-a)(x-c)(x-b)}_{\textrm{changes sign at $x=c$}} f[a,b,c,x] \ud x
\end{align}

We cannot apply the integral mean value theorem like we did for trapezoid rule. Define 

$$
w(x) = \int_a^x (t-a)(t-c)(t-b)\ud x
$$ 

Then

$$
w(a) = w(b) = 0, \qquad w(x) > 0, \quad \forall x \in (a,b)
$$ 

and

$$
w'(x) = (x-a)(x-c)(x-b)
$$ 

Hence the error can be written as

$$
\begin{aligned}
E_2(f) 
&= \int_a^b w'(x) f[a,b,c,x] \ud x \\
&= \left[ w(x) f[a,b,c,x] \right]_a^b - \int_a^b w(x) \dd{}{x} f[a,b,c,x] \ud x \\
&= - \int_a^b w(x) f[a,b,c,x,x] \ud x
\end{aligned}
$$ 

Now we can apply the integral mean value theorem

$$
\begin{aligned}
E_2(f) 
&= -f[a,b,c,\xi,\xi] \int_a^b w(x) \ud x, \quad \textrm{for some $\xi \in [a,b]$} \\
&= -\frac{1}{24} f^{(4)}(\eta) \left[ \frac{4}{15} h^5 \right], \quad \textrm{for some $\eta \in [a,b]$}
\end{aligned}
$$ 

The error in Simpson rule is given by

$$
E_3(f) = -\frac{1}{90} h^5 f^{(4)}(\eta), \quad \textrm{for some $\eta \in [a,b]$}
$$

and

$$
I(f) = I_n(f) -\frac{1}{90} h^5 f^{(4)}(\eta), \quad \textrm{for some $\eta \in [a,b]$}
$$

:::{attention}
By construction, Simpson rule is exact for quadratic polynomials. The error formula tells us that even for cubic polynomials, the error is zero.
:::

## Composite Simpson rule

Let $n \ge 2$ be even, the spacing between points $h = (b-a)/n$ and the
$n+1$ points 

$$
x_j = a + j h, \quad j=0,1,2,\ldots,n
$$ 

In the interval $[x_{2j-2}, x_{2j}]$, we have the three points $\{ x_{2j-2}, x_{2j-1}, x_{2j} \}$ using which we can apply Simpson rule.  Then

$$
\begin{aligned}
I(f) 
&= \int_a^b f(x) \ud x \\
&= \sum_{j=1}^{n/2} \int_{x_{2j-2}}^{x_{2j}} f(x) \ud x \\
&= \sum_{j=1}^{n/2}\left\{ \frac{h}{3}[f_{2j-2} + 4 f_{2j-1} + f_{2j}] - \frac{h^5}{90} f^{(4)}(\eta_j) \right\}, \quad \eta_j \in [x_{2j-2},x_{2j}]
\end{aligned}
$$ 

The composite Simpson rule is

$$
I_n(f) = \frac{h}{3}[ f_0 + 4 f_1 + 2 f_2 + \ldots + 2 f_{n-2} + 4 f_{n-1} + f_n]
$$

The error in this approximation is

$$
E_n(f) = I(f) - I_n(f) = -\frac{h^5}{90} \left(\frac{n}{2}\right) \clr{red}{ \left[ \frac{2}{n} \sum_{j=1}^{n/2} f^{(4)}(\eta_j) \right] }
$$

If the fourth derivative is continuous and since $n h = b-a$, we get

$$
E_n(f) = -\frac{1}{180} (b-a)h^4 \clr{red}{f^{(4)}(\eta)}, \quad \textrm{for some $\eta \in [a,b]$}
$$

We can also derive the asymptotic error formula

$$
E_n(f) \approx \tilde E_n(f) = -\frac{h^4}{180}[ f^{(3)}(b) - f^{(3)}(a)]
$$

The following function implements composite Simpson rule.

```{code-cell}
def simpson(a,b,n,f,df3):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res = 4.0*sum(y[1:n:2]) + 2.0*sum(y[2:n-1:2]) + y[0] + y[n]
    est = -(h**4/180.0)*(df3(b) - df3(a))
    return (h/3.0)*res, est
```

:::{prf:example}

$$
I = \int_0^\pi \ee^x \cos x \ud x = -\half[1 + \exp(\pi)]
$$

```{code-cell}
# Function f
f = lambda x: exp(x)*cos(x)
# Third derivative of f
df3 = lambda x: -2.0*exp(x)*(cos(x) + sin(x))

qe = -0.5*(1.0 + exp(pi)) # Exact integral

n,N = 2,10
e = zeros(N)
for i in range(N):
    integral,est = simpson(0.0,pi,n,f,df3)
    e[i] = qe - integral
    if i > 0:
        print('%6d %18.8e %14.5f %18.8e'%(n,e[i],e[i-1]/e[i],est))
    else:
        print('%6d %18.8e %14.5f %18.8e'%(n,e[i],0,est))
    n = 2*n
```

The convergence rate is $O(h^4)$ and the error estimate becomes very good as $n$ increases.
:::

## Compare Trapezoid and Simpson

The next function implements Trapezoid and Simpson methods.

```{code-cell}
# Performs Trapezoid and Simpson quadrature
def integrate(a,b,n,f):
    h = (b-a)/n
    x = linspace(a,b,n+1)
    y = f(x)
    res1 = 0.5*y[0] + sum(y[1:n]) + 0.5*y[n]
    res2 = 4.0*sum(y[1:n:2]) + 2.0*sum(y[2:n-1:2]) + y[0] + y[n]
    return h*res1, (h/3.0)*res2
```

The next function performs convergence test.

```{code-cell}
def test(a,b,f,Ie,n,N):
    e1,e2 = zeros(N), zeros(N)
    for i in range(N):
        I1,I2 = integrate(a,b,n,f)
        e1[i],e2[i] = Ie - I1, Ie - I2
        if i > 0:
            print('%6d %18.8e %14.5g %18.8e %14.5g'%
                 (n,e1[i],e1[i-1]/e1[i],e2[i],e2[i-1]/e2[i]))
        else:
            print('%6d %18.8e %14.5g %18.8e %14.5g'%(n,e1[i],0,e2[i],0))
        n = 2*n
```

:::{prf:example}

$$
I  = \int_0^1 x^3 \sqrt{x} \ud x = \frac{2}{9}
$$

```{code-cell}
f = lambda x: x**3 * np.sqrt(x)
a, b = 0.0, 1.0
Ie = 2.0/9.0 # Exact integral
n, N = 2, 10
test(a,b,f,Ie,n,N)
```

:::

:::{prf:example}

$$
I = \int_0^5 \frac{\ud x}{1 + (x-\pi)^2} = \tan^{-1}(5 - \pi) + \tan^{-1}(\pi)
$$

```{code-cell}
f = lambda x: 1.0/(1.0 + (x-np.pi)**2)
a, b = 0.0, 5.0
Ie = np.arctan(b-np.pi) + np.arctan(np.pi)
n, N = 2, 10
test(a,b,f,Ie,n,N)
```

:::

:::{prf:example}

$$
I = \int_0^1 \sqrt{x} \ud x = \frac{2}{3}
$$

```{code-cell}
f = lambda x: np.sqrt(x)
a, b = 0.0, 1.0
Ie = 2.0/3.0
n, N = 2, 10
test(a,b,f,Ie,n,N)
```

Since $f'(0)$ is not finite, we do not get the optimal convergence rates. But the errors still decrease and both methods converge at same rate.
:::

:::{prf:example}

$$
I = \int_0^{2\pi} \ee^{\cos x} \ud x = 7.95492652101284
$$

```{code-cell}
f = lambda x: np.exp(np.cos(x))
a, b = 0.0, 2*np.pi
Ie = 7.95492652101284
n, N = 2, 5
test(a,b,f,Ie,n,N)
```
The integrand is periodic and the error formula suggests that we should expect fast convergence. Trapezoid is more accurate than Simpson for small $n$.
:::
