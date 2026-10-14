---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Forward Euler method

```{code-cell}
from pylab import *
```

Integrating the ODE

$$
y(t) = y(t_0) + \int_{t_0}^t f(\tau, y(\tau)) \ud \tau
$$

We cannot find $y(t)$ for every $t$. Let us take small intervals in $t$, say of uniform length $h$

$$
& t_0, t_0+h, t_0+2h, \ldots \\
& t_0, t_1, t_2, \ldots
$$

i.e.,

$$
t_n = t_0 + n h, \qquad n=0,1,2,\ldots
$$

Then

$$
y(t_1) &= y(t_0) + \int_{t_0}^{t_1} f(\tau, y(\tau)) \ud \tau \\
&\approx y(t_0) + \int_{t_0}^{t_1} f(\tau, y(t_0)) \ud \tau \\
&= y_0 + h f (t_0, y_0)
$$

We obtain the approximation

$$
y_1 = y_0 + h f(t_0, y_0)
$$

This procedure can be continued to produce approximations at $t_2, t_2, \ldots$ by the following recursive scheme

$$
y_{n+1} = y_n + h f(t_n, y_n), \qquad n=0,1,2,\ldots
$$

## Geometric interpretation

```{image} https://raw.githubusercontent.com/cpraveen/numa/refs/heads/master/asy/forward_euler.svg
:width: 60%
:align: center
```

Through the point $(t_0, y(t_0))$, there is the curve of exact solution and we move along the tangent vector to this curve to reach $(t_1, y_1)$ where $y_1 = y_0 + h f_0$. 

Now through the point $(t_1, y_1)$ there is an exact solution curve and we move along the tangent at this point to reach $(t_2, y_2)$. 

This process is continued until we reach the final time.

:::{prf:remark}
Euler method can be applied with variable steps

$$
t_0 < t_1 < t_2 < \ldots, \qquad h_n = t_{n+1} - t_n
$$

and

$$
y_{n+1} = y_n + h_n f(t_n, y_n), \qquad n=0,1,2,\ldots
$$

:::

## Implementation

The following function implements Euler method.

$$
y_n = y_{n-1} + h f(t_{n-1},y_{n-1})
$$

```{code-cell}
def euler(f,t0,T,y0,h):
    N = int((T-t0)/h)
    y = zeros(N)
    t = zeros(N)
    y[0] = y0
    t[0] = t0
    for n in range(1,N):
        y[n] = y[n-1] + h*f(t[n-1],y[n-1])
        t[n] = t[n-1] + h
    return t, y
```

+++

:::{prf:example}
Consider the ODE problem

\begin{align}
y' &= y, \qquad t \ge 0 \\
y(0) &= 1
\end{align}

whose exact solution is

$$
y(t) = \exp(t)
$$

Right hand side function and exact solution

```{code-cell}
def f(t,y):
    return y

def yexact(t):
    return exp(t)
```

We run the Euler method for decreasing values of step size.

```{code-cell}
t0,y0,T = 0.0,1.0,5.0

te = linspace(t0,T,100)
ye = yexact(te)
plot(te,ye,'--',label='Exact')

H = [0.2,0.1,0.05]

for h in H:
    t,y = euler(f,t0,T,y0,h)
    plot(t,y,label="h="+str(h))

legend(loc=2), xlabel('t'), ylabel('y'), grid(True);
```

Even visually we observe that the numerical solution is converging towards the exact solution as $h \to 0$.
:::

+++

## Convergence

We expect that

$$
y_n \approx y(t_n)
$$

and as $h \to 0$, we expect that the numerical solution converges to the exact solution.

Suppose we are interested in the interval $[t_0, t_0 + T]$ We have the numerical approximation

$$
y_n = y_{n,h}, \qquad n=0,1,2,\ldots, \lfloor T/h \rfloor
$$

:::{prf:definition}
A method is said to be convergent if, for every ODE

$$
y' = f(t,y)
$$

with a Lipschitz function $f$ and every $T > 0$ it is true that

$$
\lim_{h \to 0} \max_{n=0,1,\ldots,\lfloor T/h \rfloor} \| y_{n,h} - y(t_n) \| = 0
$$
:::

:::{prf:theorem}
Euler's method is convergent.
:::

:::{prf:proof}
Given $h > 0$, define the error

$$
e_{n,h} = y_{n,h} - y(t_n), \qquad n=0,1,\ldots, \lfloor T/h \rfloor
$$

We want to show that

$$
\lim_{h \to 0} \max_n \| e_{n,h} \| = 0
$$

Bt Taylor's theorem and the ODE

$$
\begin{aligned}
y(t_{n+1}) &= y(t_n) + h y'(t_n) + O(h^2) \\
&= y(t_n) + h f(t_n, y(t_n)) + O(h^2)
\end{aligned}
$$

Euler method is

$$
y_{n+1} = y_n + h f(t_n, y_n)
$$

Then

$$
e_{n+1,h} &= e_{n,h} + h [ f(t_n, y(t_n) + e_{n,h}) - f(t_n, y(t_n)) ] + O(h^2) $$

Using triangle inequality and Lipschitz continuity of $f$

$$
\begin{aligned}
\| e_{n+1,h} \| &= \norm{e_{n,h}} + h \norm{ f(t_n, y(t_n)+e_{n,h}) - f(t_n, y(t_n)) } + c h^2 \\
&\le (1 + h \lambda) \norm{e_{n,h}} + c h^2, \qquad n=0,1,\ldots,\lfloor T/h \rfloor - 1
\end{aligned}
$$

We now claim that

$$
\norm{e_{n,h}} \le \frac{c h}{\lambda} [ (1 + h \lambda)^n - 1 ]
$$

The proof is by induction on $n$.

When $n=0$, we need $\norm{e_{0,h}} \le 0$ which is true since we know the initial condition. Assume that the inequality is true up to some $n$. Then

$$
\begin{aligned}
\norm{e_{n+1,h}} &\le (1 + h \lambda) \norm{e_{n,h}} + c h^2 \\
&\le (1 + h \lambda) \frac{c h}{\lambda} [ (1 + h \lambda)^n - 1] + c h^2 \\
&= \frac{c h}{\lambda} [ (1 + h \lambda)^{n+1} - 1 ]
\end{aligned}
$$

and it is true for $n+1$.

Since $h \lambda > 0$, using $1+x < \ee^x$ for any $x \in \re$, we have

$$
\begin{aligned}
1 + h \lambda &< \ee^{h \lambda} \\
(1 + h \lambda)^n &< \ee^{n h \lambda}
\end{aligned}
$$

Hence the error is bounded by

$$
\begin{aligned}
\norm{e_{n,h}} &\le \frac{c}{\lambda} ( \ee^{n h \lambda} - 1) h \\
\max_n \norm{e_{n,h}} &\le \frac{c}{\lambda} ( \ee^{T \lambda} - 1) h = O(h)
\end{aligned}
$$


:::

+++

:::{prf:example}
Consider the ODE problem

\begin{align}
y' &= -y + 2 \exp(-t) \cos(2t), \qquad t \ge 0 \\
y(0) &= 0
\end{align}

with exact solution

$$
y(t) = \exp(-t) \sin(2t)
$$

Right hand side function and exact solution

```{code-cell} 
def f(t,y):
    return -y + 2.0*exp(-t)*cos(2.0*t)

def yexact(t):
    return exp(-t)*sin(2.0*t)
```

**Solve for a given h**

```{code-cell}
t0,y0 = 0.0,0.0
T  = 10.0
h  = 1.0/20.0
t,y = euler(f,t0,T,y0,h)
te = linspace(t0,T,100)
ye = yexact(te)
plot(t,y,te,ye,'--')
legend(('Numerical','Exact'))
xlabel('t'), ylabel('y'), grid(True)
title('Step size = ' + str(h));
```

**Solve for several h**

Study the effect of decreasing step size. The error is plotted in log scale.

```{code-cell}
hh = 0.1/2.0**arange(5)
err= zeros(len(hh))
for (i,h) in enumerate(hh):
    t,y = euler(f,t0,T,y0,h)
    ye = yexact(t)
    err[i] = abs(y-ye).max()
    semilogy(t,abs(y-ye))
    
legend(hh), xlabel('t'), ylabel('log(error)')

for i in range(1,len(hh)):
    print('%e %e %e'%(hh[i],err[i],log(err[i-1]/err[i])/log(2.0)))

# plot error vs h
figure()
loglog(hh,err,'o-',label="Error")
loglog(hh,hh,'--',label="y=h")
legend(), xlabel('h'), ylabel('Error norm');
```
:::

:::{prf:remark}
The constant $c$ can be taken to be

$$
c = \max_{[t_0, t_0+T]} \norm{y''(t)}
$$

The previous theorem seems to give the error estimate

$$
\norm{e_{n,h}} \le \frac{c h}{\lambda} ( \ee^{T \lambda} - 1)
$$
:::

:::{prf:example}
Consider the problem

$$
y' = -100 y, \qquad y(0) = 1
$$

whose solution is 

$$
y(t) = \ee^{-100t}
$$

Then $\lambda = 100$, $c = \lambda^2 = 10^4$ and

$$
\norm{e_{n,h}} \le 100 h (\ee^{100 T} - 1)
$$

Taking $T=1$,

$$
|y_n - y(nh)| \le 2.69 \times 10^{45} h
$$

But $y_n = (1- 100 h)^n$ and the actual error is

$$
|y_n - y(nh)| = |(1-100h)^n - \ee^{-100 n h}|
$$

This is smaller by orders of magnitude than the error estimate above. In this case the error estimate is very pessimistic, it over-estimates the error.
:::
