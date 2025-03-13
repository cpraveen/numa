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

# ODE System

```{code-cell} ipython3
from pylab import *
```

+++

:::{prf:example}
Consider a system of two coupled ODE

$$
y' = Ay, \qquad A = \begin{bmatrix}
-100 & 1 \\
0 & -1/10 \end{bmatrix}
$$

with initial condition

$$
y(0) = y_0
$$

The $\theta$ scheme is given by

\begin{align}
y^{n+1} &= y^n + h [(1-\theta) Ay^n + \theta Ay^{n+1}] \\
&\Downarrow \\
[I-\theta hA] y^{n+1} &= y^n + (1-\theta)h A y^n
\end{align}

The following function implements the $\theta$ scheme.

```{code-cell} ipython3
def solve(y0,h,T,theta):
    A = matrix([[-100.0, 1.0], [0.0, -0.1]])
    A1 = eye(2) - h*theta*A
    A1 = linalg.inv(A1)
    A2 = eye(2) + (1.0-theta)*h*A
    M  = A1.dot(A2)
    N = int(T/h)
    y = ndarray((N,2))
    t = zeros(N)
    y[0,:]=y0
    t[0] = 0.0
    for i in range(N-1):
        y[i+1,:] = M.dot(y[i,:])
        t[i+1] = t[i] + h
    figure(figsize=(9,4))
    subplot(1,2,1), plot(t,y[:,0])
    xlabel('t'), title('First component')
    subplot(1,2,2), plot(t,y[:,1])
    xlabel('t'), title('Second component')
    
# Initial condition
y0 = array([10.0/999.0,1.0])

# Final time
T = 25.0
```

**Forward Euler, $\theta=0$**

The stability condition is

$$
|1-100h| < 1, \qquad |1-h/10| < 1
$$

which means that we need to choose $h < 1/50$. We first try with $h=0.01$

```{code-cell} ipython3
h = 0.01
solve(y0,h,T,0.0)
```

We now try with $h=0.021$ which is slightly larger than the stability limit.

```{code-cell} ipython3
h = 0.021
solve(y0,h,T,0.0)
```

**Trapezoidal scheme, $\theta = 0.5$**

```{code-cell} ipython3
h = 0.021
solve(y0,h,T,0.5)
```
:::
