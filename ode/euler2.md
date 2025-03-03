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

# Instability of forward Euler scheme

```{code-cell}
from pylab import *
```

+++

:::{prf:example}
Consider the ODE problem

\begin{align}
y' &= -5 t y^2 + \frac{5}{t} - \frac{1}{t^2}, \qquad t \ge 1 \\
y(1) &= 1
\end{align}

with exact solution

$$
y(t) = \frac{1}{t}
$$

Right hand side function

```{code-cell} ipython3
def f(t,y):
    return -5.0*t*y**2 + 5.0/t - 1.0/t**2
```

Exact solution

```{code-cell} ipython3
def yexact(t):
    return 1.0/t
```

This implements Euler method
$$
y_n = y_{n-1} + h f(t_{n-1},y_{n-1})
$$

```{code-cell} ipython3
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

```{code-cell} ipython3
t0 = 1.0
y0 = 1.0
T  = 25.0
h  = 0.19
t,y = euler(f,t0,T,y0,h)
print('Number of iterations=',len(t))
te = linspace(t0,T,100)
ye = yexact(te)
plot(t,y,te,ye,'--')
legend(('Numerical','Exact'))
xlabel('t')
ylabel('y')
title('Step size = ' + str(h));
```

```{code-cell} ipython3
h  = 0.21
t,y = euler(f,t0,T,y0,h)
print('Number of iterations=',len(t))
plot(t,y,te,ye,'--')
legend(('Numerical','Exact'))
xlabel('t')
ylabel('y')
title('Step size = ' + str(h));
```
:::

## Adaptive stepping

+++

We will use a variable time step
$$
y_n = y_{n-1} + h_{n-1} f(t_{n-1}, y_{n-1})
$$
where
$$
h_{n-1} = \frac{1}{|f_y(t_{n-1},y_{n-1})|} = \frac{1}{10 t_{n-1} |y_{n-1}|}
$$

```{code-cell} ipython3
def aeuler(f,t0,T,y0):
    t, y = [], []
    y.append(y0)
    t.append(t0)
    time = t0; n = 1
    while time < T:
        h    = 1.0/abs(10*t[n-1]*y[n-1])
        y.append(y[n-1] + h*f(t[n-1],y[n-1]))
        time = time + h
        t.append(time)
        n = n + 1
    return array(t), array(y)
```

```{code-cell} ipython3
t,y = aeuler(f,t0,T,y0)
print('Number of iterations=',len(t))

plot(t,y,te,ye,'--')
legend(('Numerical','Exact'))
xlabel('t')
ylabel('y');

figure()
plot(range(len(t)-1),t[1:]-t[0:-1])
ylabel('Step size, h')
xlabel('Iteration');
```
