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

# Shooting method for BVP

```{code-cell}
from pylab import *
```

+++

:::{prf:example}
Consider the BVP

$$
y'' = -y + \frac{2 (y')^2}{y}, \qquad -1 < x < +1
$$

with boundary conditions

$$
y(-1) = y(+1) = (e + e^{-1})^{-1}
$$

The exact solution is

$$
y(x) = (e^x + e^{-x})^{-1}
$$

**Formulate as an initial value problem**

$$
y'' = -y + \frac{2 (y')^2}{y}, \qquad -1 < x < +1
$$

with initial conditions

$$
y(-1) = (e + e^{-1})^{-1}, \qquad y'(-1) = s
$$

Let the solution of this IVP be denoted as $y(x;s)$. We have to determine $s$ as the root of the function

$$
\phi(s) = y(1;s) - (e + e^{-1})^{-1} = 0
$$

The exact root is

$$
s = y'(-1) = \frac{e - e^{-1}}{(e + e^{-1})^2} \approx 0.24677717
$$

**Newton method**

To find the root of $\phi(s)$, we use Newton method with an initial guess $s_0$. Then the Newton method updates the guess by

$$
s_{m+1} = s_m - \frac{\phi(s_m)}{\frac{d}{ds}\phi(s_m)}, \qquad m=0,1,2,\ldots
$$

Define

$$
z_s(x) = \frac{\partial}{\partial s} y(x;s)
$$

Then

$$
z_s'(x) = \frac{\partial}{\partial s} y'(x;s)
$$

The derivative of $\phi(s)$ is given by

$$
\frac{d}{ds}\phi(s) = \frac{\partial}{\partial s} y(1;s) = z_s(1)
$$

We can write an equation for $z_s(x)$

$$
z_s'' = \left[ -1 - 2 \left( \frac{y'}{y} \right)^2 \right] z_s + 4 \frac{y'}{y} z_s'
$$

with initial conditions

$$
z_s(-1) = 0, \qquad z_s'(-1) = 1
$$

**First order ODE system**

Define the vector

$$
u = \begin{bmatrix}
u_1 \\ u_2 \\ u_3 \\ u_4 \end{bmatrix} = \begin{bmatrix}
y \\ y' \\ z_s \\ z_s' \end{bmatrix}
$$

Then

\begin{align}
u_1' &= u_2 \\ 
u_2' &= -u_1 + \frac{2 u_2^2}{u_1} \\
u_3' &= u_4 \\
u_4' &= \left[ -1 - 2 \left( \frac{u_2}{u_1} \right)^2 \right] u_3 + 4 \frac{u_2}{u_1} u_4
\end{align}

Hence we get the first order ODE system

$$
u' = f(u) = \begin{bmatrix}
u_2 \\
-u_1 + \frac{2 u_2^2}{u_1} \\
u_4 \\
\left[ -1 - 2 \left( \frac{u_2}{u_1} \right)^2 \right] u_3 + 4 \frac{u_2}{u_1} u_4
\end{bmatrix}
$$

with initial condition

$$
u(-1) = \begin{bmatrix}
(e+e^{-1})^{-1} \\
s \\
0 \\
1
\end{bmatrix}
$$

Once we solve this IVP, we get

$$
\phi(s) = u_1(1) - (e + e^{-1})^{-1}, \qquad \frac{d}{ds}\phi(s) = z_s(1) = u_3(1)
$$

We now code the problem specific data.

```{code-cell} 
def f(u):
    rhs = zeros(4)
    rhs[0] = u[1]
    rhs[1] = -u[0] + 2.0*u[1]**2/u[0]
    rhs[2] = u[3]
    rhs[3] = (-1.0-2.0*(u[1]/u[0])**2)*u[2] + 4.0*u[1]*u[3]/u[0]
    return rhs

def initialcondition(s):
    u = zeros(4)
    u[0] = 1.0/(exp(1)+exp(-1))
    u[1] = s
    u[2] = 0.0
    u[3] = 1.0
    return u

def yexact(x):
    return 1.0/(exp(x) + exp(-x))
```

The function below implements the two step RK scheme.

```{code-cell} 
def solvephi(n,s):
    h = 2.0/n
    u = zeros((n+1,4))
    u[0] = initialcondition(s)
    for i in range(1,n+1):
        u1 = u[i-1] + 0.5*h*f(u[i-1])
        u[i] = u[i-1] + h*f(u1)
    phi = u[n][0] - 1.0/(exp(1)+exp(-1))
    dphi= u[n][2]
    return phi,dphi,u
```

Let us see the solution corresponding to the initial guess.

```{code-cell} 
n = 100
s = 0.2
p, dp, u = solvephi(n,s)
x = linspace(-1.0,1.0,n+1)
xe = linspace(-1.0,1.0,1000)
ye = yexact(xe)
plot(x,u[:,0],xe,ye)
grid(True), title("s = " + str(s))
legend(("Numerical","Exact"));
```

The initial guess is far from the solution. We will not start the shooting method.

```{code-cell} 
def shoot(s,n,plotsol=False):
    maxiter, TOL = 100, 1.0e-8
    x  = linspace(-1.0,1.0,n+1)
    if plotsol:
        figure(figsize=(8,5))
        xe = linspace(-1.0,1.0,1000)
        ye = yexact(xe)
        plot(xe,ye,label="Exact")
    it = 0
    while it < maxiter:
        p, dp, u = solvephi(n,s)
        if plotsol:
            plot(x,u[:,0],label=str(it))
        if abs(p) < TOL:
            break
        s = s - p/dp
        it += 1
        print('%5d %16.6e %16.6e'%(it,s,abs(p)))
    if plotsol:
        legend(); grid(True)
    return s,p,x,u
```

Solve with initial guess of $s=0.2$ and step size in RK method of $h = 2/n$ where $n=64$.

```{code-cell}
s, n = 0.2, 64
s,p,x,u = shoot(s,n,True)
print("Root s        = %e" % s)
print("phi(s)        = %e" % p)
sexact = (exp(1) - exp(-1))/(exp(1) + exp(-1))**2
print("Error in root = ",  sexact-s)
figure(figsize=(8,5))
plot(x,u[:,0],xe,ye)
grid(True)
title("Final solution, s = " + str(s))
legend(("Numerical","Exact"));
```
:::
