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

# Stiff ODE

```{code-cell} ipython3
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

```{code-cell} ipython3
def f(t,y):
    return -100.0*(y - sin(t))

def yexact(t):
    return 10101.0*exp(-100*t)/10001.0 - 100.0*(cos(t) - 100.0*sin(t))/10001.0
```

**Forward Euler**

$$
y_n = y_{n-1} - 100 h [ y_{n-1} - \sin(t_{n-1})], \qquad n=1,2,\ldots
$$

```{code-cell} ipython3
def ForwardEuler(t0,y0,T,h):
    N = int((T-t0)/h)
    y = zeros(N)
    t = zeros(N)
    t[0] = t0
    y[0] = y0
    for n in range(1,N):
        y[n] = y[n-1] + h*f(t[n-1],y[n-1])
        t[n] = t[n-1] + h
    return t,y
```

Run with a step size just within the stability limit

$$
h = 0.95 \frac{2}{100} < \frac{2}{|\lambda|}
$$

```{code-cell} ipython3
t0 = 0
y0 = 1
T  = 2*pi

h  = 0.95*(2.0/100.0)
t,y = ForwardEuler(t0,y0,T,h)

plot(t,y,t,yexact(t))
title("Forward Euler, h ="+str(h))
legend(("Euler","Exact"))
```

```{code-cell} ipython3
h  = 1.0/100.0
t,y = ForwardEuler(t0,y0,T,h)

plot(t,y,t,yexact(t))
title("Forward Euler, h ="+str(h))
legend(("Euler","Exact"));
```

**Backward Euler scheme**

\begin{align}
y_n &= y_{n-1} + h [ -100(y_n - \sin(t_n) ] \\
&\Downarrow \\
y_n &= \frac{ y_{n-1} + 100 h \sin(t_n)}{1 + 100 h}
\end{align}

```{code-cell} ipython3
def BackwardEuler(t0,y0,T,h):
    N = int((T-t0)/h)
    y = zeros(N)
    t = zeros(N)
    t[0] = t0
    y[0] = y0
    for n in range(1,N):
        t[n] = t[n-1] + h
        y[n] = (y[n-1] + 100.0*h*sin(t[n]))/(1.0 + 100.0*h)
    return t,y
```

**Trapezoidal scheme**

\begin{align}
y_n &= y_{n-1} + \frac{h}{2}[ -100(y_{n-1}-\sin(t_{n-1})) - 100 (y_n - \sin(t_n))] \\
&\Downarrow \\
y_n &= \frac{(1-50h)y_{n-1} + 50h(\sin(t_{n-1}) + \sin(t_n))}{1 + 50 h}
\end{align}

```{code-cell} ipython3
def Trapezoidal(t0,y0,T,h):
    N = int((T-t0)/h)
    y = zeros(N)
    t = zeros(N)
    t[0] = t0
    y[0] = y0
    for n in range(1,N):
        t[n] = t[n-1] + h
        y[n] = ((1.0-50.0*h)*y[n-1] + 50.0*h*(sin(t[n-1])+sin(t[n])))/(1.0 + 50.0*h)
    return t,y
```

We now run both backward Euler and trapezoidal scheme.

```{code-cell} ipython3
h  = 10.0/100.0
t,y = BackwardEuler(t0,y0,T,h)
plot(t,y)

t,y = Trapezoidal(t0,y0,T,h)
plot(t,y)

plot(t,yexact(t))
title("Step size h ="+str(h))
legend(("Backward Euler","Trapezoid","Exact"));
```

```{code-cell} ipython3
h  = 2.0/100.0
t,y = BackwardEuler(t0,y0,T,h)
plot(t,y)

t,y = Trapezoidal(t0,y0,T,h)
plot(t,y)

plot(t,yexact(t))
title("Step size h ="+str(h))
legend(("Backward Euler","Trapezoid","Exact"));
```
:::

:::{exercise}
Implement BDF2 scheme with BDF1 in the first time step. This should give better results than Trapezoid and is second order accurate.
:::
