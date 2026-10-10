---
exports:
  - format: pdf
    template: arxiv_nips
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# ODE

```{code-cell}
from pylab import *
```

Given an initial condition $y_0 \in \re^d$ and a right hand side function

$$
f : (t,y) \to f(t,y) \in \re^d
$$

consider the initial value problem

$$
y' &= f(t,y), \quad t \ge t_0 \\
y(t_0) &= y_0
$$

$f$ is continuous in $t$ and $y$ and Lipschitz continuous in $y$

$$
\| f(t,x) - f(t,y) \| \le \lambda \| x - y \|, \qquad x,y \in \re^d, \quad t \ge t_0
$$

* For any $y_0 \in \re^d$, there exists a unique solution. This solution is differentiable.

* The solution depends continuously on the initial data; if $\hat y$ also satisfies the ODE but with initial condition $\hat y_0$, then

$$
\|y(t) - \hat y(t) \| \le \ee^{\lambda(t-t_0)} \| y_0 - \hat y_0 \|
$$

* If $\hat y$ satisfies a perturbed ODE

$$
\hat y' = f(t,\hat y) + r(t,\hat y)
$$

where $r$ is bounded

$$
\| r \| \le M
$$

then

$$
\|y(t) - \hat y(t) \| \le \ee^{\lambda(t-t_0)} \| y_0 - \hat y_0 \| + \frac{M}{\lambda} \left( \ee^{\lambda(t-t_0)} - 1 \right)
$$
