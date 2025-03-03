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
from pylab import *
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
    return exp(-t)*sin(2.0*t)
```

The trapezoid method is

\begin{align}
y_n 
&= y_{n-1} + \frac{h}{2}[ f(t_{n-1},y_{n-1}) + f(t_n,y_n)] \\
&= y_{n-1} + \frac{h}{2}[ -y_{n-1} + 2 \exp(-t_{n-1}) \cos(2t_{n-1})  - y_n + 2 \exp(-t_n) \cos(2t_n) ]
\end{align}

Solving for $y_n$

$$
y_n = \frac{1}{1 + \frac{h}{2}} \left\{ (1 - \frac{h}{2}) y_{n-1} + h [\exp(-t_{n-1}) \cos(2t_{n-1}) + \exp(-t_n) \cos(2t_n) ] \right\}
$$

This implements Trapezoidal method

```{code-cell}
def trap(t0,T,y0,h):
    N = int((T-t0)/h)
    y = zeros(N)
    t = zeros(N)
    y[0] = y0
    t[0] = t0
    for n in range(1,N):
        t[n] = t[n-1] + h
        y[n] = (1.0-0.5*h)*y[n-1] + \
               h*(exp(-t[n-1])*cos(2.0*t[n-1]) + exp(-t[n])*cos(2.0*t[n]))
        y[n] = y[n]/(1.0 + 0.5*h)
    return t, y
```

**Solve for fixed h**

```{code-cell}
t0,y0 = 0.0,0.0
T  = 10.0
h  = 1.0/20.0
t,y = trap(t0,T,y0,h)
te = linspace(t0,T,100)
ye = yexact(te)
plot(t,y,te,ye,'--')
legend(('Numerical','Exact'))
xlabel('t')
ylabel('y')
title('Step size = ' + str(h));
```

**Solve for decreasing h**

Study the effect of decreasing step size. The error is plotted in log scale.

```{code-cell}
hh = 0.1/2.0**arange(5)
err=zeros(len(hh))
for (i,h) in enumerate(hh):
    t,y = trap(t0,T,y0,h)
    ye = yexact(t)
    err[i] = abs(y-ye).max()
    semilogy(t,abs(y-ye))
    
legend(hh)
xlabel('t')
ylabel('log(error)')

for i in range(1,len(hh)):
    print('%e %e %e'%(hh[i],err[i],log(err[i-1]/err[i])/log(2.0)))

# plot error vs h
figure()
loglog(hh,err,'o-',label="Error")
loglog(hh,hh**2/4,'--',label="$y=h^2/4$")
legend(), xlabel('h'), ylabel('Error norm');
```
:::

+++

:::{prf:example} Harmonic oscillator

```{code-cell}
y0    = 1.0
omega = 4.0
h     = 0.1

A = array([[0.0,       1.0], 
           [-omega**2, 0.0]])
Afe = eye(2) + h * A
Abe = inv(eye(2) - h * A)
At1 = eye(2) - 0.5 * h * A
At2 = eye(2) + 0.5 * h * A
Atr = inv(At1) @ At2

n = int(2*pi/h)
yfe, ybe, ytr = zeros((n,2)), zeros((n,2)), zeros((n,2))
yfe[0,0] = y0
ybe[0,0] = y0 
ytr[0,0] = y0
t = zeros(n)
for i in range(1,n):
    yfe[i,:] = Afe @ yfe[i-1,:]
    ybe[i,:] = Abe @ ybe[i-1,:]
    ytr[i,:] = Atr @ ytr[i-1,:]
    t[i]     = t[i-1] + h

plot(t, cos(omega*t),'--',label="Exact")
plot(t,yfe[:,0],label="FE")
plot(t,ybe[:,0],label="BE")
plot(t,ytr[:,0],label="TR")
axis([0, 2*pi,-2,2])
legend(), xlabel("t"), ylabel("$\\theta$");
```

:::
