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

# Trapezoid method

```{code-cell}
import numpy as np
from matplotlib import pyplot as plt
```

+++

:::{prf:example}
Consider the ODE

$$
y' = -y + 2 \exp(-t) \cos(2t)
$$

with initial condition

$$
y(0) = 0
$$

The exact solution is

$$
y(t) = \exp(-t) \sin(2t)
$$

Exact solution

```{code-cell}
def yexact(t):
    return np.exp(-t)*np.sin(2.0*t)
```

This implements Trapezoidal method

$$
y_n = y_{n-1} + \frac{h}{2}[ f(t_{n-1},y_{n-1}) + f(t_n,y_n)]
$$

For the present example we get

$$
y_n = y_{n-1} + \frac{h}{2}[ -y_{n-1} + 2 \exp(-t_{n-1}) \cos(2t_{n-1})  - y_n + 2 \exp(-t_n) \cos(2t_n) ]
$$

Solving for $y_n$

$$
y_n = \frac{1}{1 + \frac{h}{2}} \left\{ (1 - \frac{h}{2}) y_{n-1} + h [\exp(-t_{n-1}) \cos(2t_{n-1}) + \exp(-t_n) \cos(2t_n) ] \right\}
$$

```{code-cell}
def trap(t0,T,y0,h):
    N = int((T-t0)/h)
    y = np.zeros(N)
    t = np.zeros(N)
    y[0] = y0
    t[0] = t0
    for n in range(1,N):
        t[n] = t[n-1] + h
        y[n] = (1.0-0.5*h)*y[n-1] + \
               h*(np.exp(-t[n-1])*np.cos(2.0*t[n-1]) + np.exp(-t[n])*np.cos(2.0*t[n]))
        y[n] = y[n]/(1.0 + 0.5*h)
    return t, y
```

**Solve for fixed h**

```{code-cell}
t0,y0 = 0.0,0.0
T  = 10.0
h  = 1.0/20.0
t,y = trap(t0,T,y0,h)
te = np.linspace(t0,T,100)
ye = yexact(te)
plt.plot(t,y,te,ye,'--')
plt.legend(('Numerical','Exact'))
plt.xlabel('t')
plt.ylabel('y')
plt.title('Step size = ' + str(h));
```

**Solve for decreasing h**

Study the effect of decreasing step size. The error is plotted in log scale.

```{code-cell}
hh = 0.1/2.0**np.arange(5)
err=np.zeros(len(hh))
for (i,h) in enumerate(hh):
    t,y = trap(t0,T,y0,h)
    ye = yexact(t)
    err[i] = np.abs(y-ye).max()
    plt.semilogy(t,np.abs(y-ye))
    
plt.legend(hh)
plt.xlabel('t')
plt.ylabel('log(error)')

for i in range(1,len(hh)):
    print('%e %e %e'%(hh[i],err[i],np.log(err[i-1]/err[i])/np.log(2.0)))

# plot error vs h
plt.figure()
plt.loglog(hh,err,'o-',label="Error")
plt.loglog(hh,hh**2/4,'--',label="$y=h^2/4$")
plt.legend(), plt.xlabel('h'), plt.ylabel('Error norm');
```
:::
