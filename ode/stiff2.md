---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Stiff ODE: Forward Euler with variable step

```{code-cell}
from pylab import *
```

+++

:::{prf:example}
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

The next two functions implement the right hand side of the ode and the exact solution.

```{code-cell} ipython3
def f(t,y):
    return -100.0*(y - sin(t))

def yexact(t):
    return 10101.0*exp(-100*t)/10001.0 - 100.0*(cos(t) - 100.0*sin(t))/10001.0
```

**Forward Euler**

$$
y_n = y_{n-1} - 100 h [ y_{n-1} - \sin(t_{n-1})]
$$

The following function uses step size h = hlo for $t \in [0,0.1]$ and h = hhi for $t > 0.1$.

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
    t, y = array(t),array(y)
    plot(t,y,t,yexact(t))
    title("Forward Euler")
    grid(True)
    legend(("Euler","Exact"));
    return t, y
```

Use a step size slightly below the stability limit of 0.02 for the entire simulation.

```{code-cell} ipython3
t0, y0 = 0.0, 1.0
T  = 2*pi
t,y = ForwardEuler(t0,y0,T,0.019,0.019)
```

The solution is oscillatory initially; use a smaller step size of 0.01.

```{code-cell} ipython3
t,y = ForwardEuler(t0,y0,T,0.01,0.01)
```

The exponential decay happens quite fast, so we may think of using a small step initially and then increase it once the decay is completed.  We start with a step size of $h=0.01$ until $t=0.1$ after which we increase it  to $h=0.021$ which is slightly above the stability limit of 0.02.

```{code-cell} ipython3
t,y = ForwardEuler(t0,y0,T,0.01,0.021)
```

The exponential mode might have decayed but we have to still use a small step size for the entire simulation, since this mode can grow back and lead to instability.
:::

+++

:::{exercise}
Use the idea of different step size with Trapezoid method. Since Trapezoid is A-stable, this strategy will give good solutions.
:::
