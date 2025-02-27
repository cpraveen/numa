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

# Multistep methods

```{code-cell}
import numpy as np
from matplotlib import pyplot as plt
```

+++

:::{exercise}
Show that the implicit multi-step method

$$
y_{n+2} - 3 y_{n+1} + 2 y_n = h \left[ \frac{13}{12} f_{n+2} - \frac{5}{3} f_{n+1} - \frac{5}{12} f_n \right]
$$

is a consistent method of order 2.
:::

+++

:::{prf:example}
Apply the method of previous exercise to

$$
y'(t) = 0, \qquad t \ge 0, \qquad y(0) = 1/3
$$

which gives

\begin{align}
y_{n+2} - 3 y_{n+1} + 2 y_n &= 0 \\
y_{n+2} &= 3 y_{n+1} - 2 y_n \\
\textrm{Use starting values} \qquad y_0 &= y_1 = \frac{1}{3}
\end{align}

The exact solution of the numerical scheme is $y_n = 1/3$ for all $n$.

```{code-cell}
N = 100
y = np.zeros(N)
y[0] = 1.0/3.0
y[1] = 1.0/3.0
for n in range(2,N):
    y[n] = 3.0 * y[n-1] - 2.0 * y[n-2]
plt.semilogy(abs(y))
plt.grid(True), plt.xlabel("n"), plt.ylabel("$abs(y_n)$");
```

On the computer the solutions blow up due to amplification of rounding errors, indicating some instability of the method. Consistency of the method is not sufficient to guarantee good solutions on the computer, the methods need to have some stability.
:::

+++

:::{exercise}
Show that the 2-step method

$$
y_{n+2} - 2.01 y_{n+1} + 1.01 y_n = h[0.995 f_{n+1} - 1.005 f_n]
$$

is second order accurate.
:::

+++

:::{prf:example}
Consider the ODE problem

\begin{gather}
y' = -y, \quad t \ge 0 \\
y(0) = 1
\end{gather}

whose exact solution is

$$
y(t) = \exp(-t)
$$

Let us use a 2-step, second order method

$$
y_{n+2} - 2.01 y_{n+1} + 1.01 y_n = h[0.995 f_{n+1} - 1.005 f_n]
$$

We will use the initial conditions

$$
y_0 = 1, \qquad y_1 = y(-h) = \exp(-h)
$$

where we used the exact solution for $y_1$. Then solving for $y_{n+2}$

$$
y_{n+2} =  2.01 y_{n+1} - 1.01 y_n + h[0.995 f_{n+1} - 1.005 f_n]
$$

The following code implements the trapezoid method.
Exact solution

The next functions implement this method.

```{code-cell}
def yexact(t):
    return np.exp(-t)

def f(t,y):
    return -y

def solve(t0,T,y0,h):
    N = int((T-t0)/h)
    y = np.zeros(N)
    t = np.zeros(N)
    t[0], y[0] = t0, y0
    t[1], y[1] = t0+h, np.exp(-(t0+h))
    for n in range(2,N):
        t[n] = t[n-1] + h
        y[n] = 2.01*y[n-1] - 1.01*y[n-2] \
               + h*(0.995*f(t[n-1],y[n-1]) - 1.005*f(t[n-2],y[n-2]))
    return t, y
```

**Solve for decreasing h**

```{code-cell}
T  = 15.0
t0, y0 = 0.0, 1.0
H  = [1.0/10.0,1.0/20.0,1.0/40]

te = np.linspace(t0,T,100)
ye = yexact(te)
plt.figure(figsize=(5,3))
plt.plot(te,ye,'--')

for h in H:
    t,y = solve(t0,T,y0,h)
    plt.plot(t,y)
plt.legend(('Exact','h=1/10','h=1/20','h=1/40'))
plt.xlabel('t')
plt.ylabel('y')
plt.axis([0,15,-0.1,1])
plt.grid(True);
```

The solution becomes worse with decreasing values of $h$. The scheme has the polynomials

$$
\rho(w) = w^2 - 2.01 w + 1.01, \qquad \sigma = 0.995 w - 1.005
$$

The roots of $\rho$ are $1$ and $1.01$. For the homogeneous problem, the 2-step scheme has the solution

$$
y_n = A + B (1.01)^n
$$

The second root is responsible for existence of growing solutions in the above scheme. As we reduce $h$ we take more steps to reach the final time, which leads to larger growth of the spurious solution.
:::
