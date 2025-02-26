---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


# ODE using forward Euler


Consider the ODE problem

\begin{gather}
y' = y, \qquad t \ge 0 \\
y(0) = 1
\end{gather}

whose exact solution is

$$
y(t) = \exp(t)
$$

```{code-cell}
import numpy as np
from matplotlib import pyplot as plt
```

Right hand side function

```{code-cell}
def f(t,y):
    return y
```

Exact solution

```{code-cell}
def yexact(t):
    return np.exp(t)
```

This implements Euler method

$$
y_n = y_{n-1} + h f(t_{n-1},y_{n-1})
$$

```{code-cell} ipython3
def euler(t0,T,y0,h):
    N = int((T-t0)/h)
    y = np.zeros(N)
    t = np.zeros(N)
    y[0] = y0
    t[0] = t0
    for n in range(1,N):
        y[n] = y[n-1] + h*f(t[n-1],y[n-1])
        t[n] = t[n-1] + h
    return t, y
```

```{code-cell} ipython3
t0,y0,T = 0.0,1.0,5.0

te = np.linspace(t0,T,100)
ye = yexact(te)
plt.plot(te,ye,'--',label='Exact')

H = [0.2,0.1,0.05]

for h in H:
    t,y = euler(t0,T,y0,h)
    plt.plot(t,y,label="h="+str(h))

plt.legend(loc=2)
plt.xlabel('t')
plt.ylabel('y')
plt.grid(True);
```
