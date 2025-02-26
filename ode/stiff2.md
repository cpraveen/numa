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

# Stiff ODE: Forward Euler with variable step

+++

Consider the stiff ODE
$$
y' = -100(y - \sin(t)), \qquad t \ge 0
$$
$$
y(0) = 1
$$
The exact solution is
$$
y(t) = \frac{10101}{10001}e^{-100t} - \frac{100}{10001}(\cos t - 100 \sin t)
$$

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
```

The next two functions implement the right hand side of the ode and the exact solution.

```{code-cell} ipython3
def f(t,y):
    return -100.0*(y - np.sin(t))

def yexact(t):
    return 10101.0*np.exp(-100*t)/10001.0 - 100.0*(np.cos(t) - 100.0*np.sin(t))/10001.0
```

## Forward Euler

+++

$$
y_n = y_{n-1} - 100 h [ y_{n-1} - \sin(t_{n-1})]
$$
We start with a step size of $h=0.01$ until $t=0.1$ after which we increase it to $h=0.021$ which is slightly above the stability limit of 0.02.

```{code-cell} ipython3
def ForwardEuler(t0,y0,T,hlo,hhi):
    t, y = [], []
    t.append(t0)
    y.append(y0)
    h = hlo
    n = 1
    while t[n-1]+h <= T:
        if t[n-1] > 0.1:
            h = hhi
        y.append(y[n-1] + h*f(t[n-1],y[n-1]))
        t.append(t[n-1] + h)
        n += 1
    t, y = np.array(t),np.array(y)
    plt.plot(t,y,t,yexact(t))
    plt.title("Forward Euler")
    plt.grid(True)
    plt.legend(("Euler","Exact"));
    return t, y
```

Use a step size slightly below the stability limit of 0.02.

```{code-cell} ipython3
t0, y0 = 0.0, 1.0
T  = 2*np.pi
t,y = ForwardEuler(t0,y0,T,0.019,0.019)
```

The solution is oscillatory initially. Use a smaller step size of 0.01

```{code-cell} ipython3
t,y = ForwardEuler(t0,y0,T,0.01,0.01)
```

The exponential decay happens quite fast, so we may think of using a small step initially and then increase it once the decay is completed.

```{code-cell} ipython3
t,y = ForwardEuler(t0,y0,T,0.01,0.021)
```

The exponential mode might have decayed but we have to still use a small step size for the entire simulation, since this mode can grow back and lead to instability.
