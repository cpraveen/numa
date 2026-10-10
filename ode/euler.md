---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Forward Euler

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
y(t_1) &= y(t_0) + \int_{t_0}^t f(\tau, y(\tau)) \ud \tau \\
&\approx y(t_0) + \int_{t_0}^t f(\tau, y(t_0)) \ud \tau \\
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
