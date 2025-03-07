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

# Runge-Kutta method

```{include} ../math.md
```

```{code-cell}
from pylab import *
```

+++

:::{prf:example} Harmonic oscillator
We repeat the previous example with RK method.

```{code-cell}
y0    = 1.0
omega = 4.0
h     = 0.1

A = matrix([[0.0,       1.0], 
            [-omega**2, 0.0]])
hA = h * A
A2 = eye(2) + hA + (1.0/2.0)*hA**2
A4 = eye(2) + hA + (1.0/2.0)*hA**2 + (1.0/6.0)*hA**3 + (1.0/24.0)*hA**4

n = int(2*pi/h)
y2, y4 = zeros((n,2)), zeros((n,2))
y2[0,0] = y0
y4[0,0] = y0 
t = zeros(n)
for i in range(1,n):
    y2[i,:] = A2 @ y2[i-1,:]
    y4[i,:] = A4 @ y4[i-1,:]
    t[i]     = t[i-1] + h

plot(t, cos(omega*t),'--',label="Exact")
plot(t,y2[:,0],label="RK2")
plot(t,y4[:,0],label="RK4")
axis([0, 2*pi,-2,2])
legend(), xlabel("t"), ylabel("$\\theta$");
```

The solution from RK2 is growing with time, no matter how small a time step we take. RK4 will be conditionally stable, if $h < 2\sqrt{2}/\omega \approx 0.707$.
:::

:::{prf:example} Stability regions of RK

```{code-cell}
r2 = lambda z: 1 + z + z**2/2
r3 = lambda z: 1 + z + z**2/2 + z**3/6
r4 = lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24

n = 100
x = linspace(-3,0.5,n)
y = linspace(-3,3,n)
X,Y = meshgrid(x,y)
Z2 = r2(X + 1j*Y)
Z3 = r3(X + 1j*Y)
Z4 = r4(X + 1j*Y)
contour(X,Y,abs(Z2),levels=[1],colors="black")
contour(X,Y,abs(Z3),levels=[1],colors="blue")
contour(X,Y,abs(Z4),levels=[1],colors="red")
xlabel('real(z)'), ylabel('imag(z)')
title('Stability domain of RK2, RK3, RK4')
grid(True), axis('equal');
```

Zoomed view around the imaginary axis

```{code-cell}
contour(X,Y,abs(Z2),levels=[1],colors="black")
contour(X,Y,abs(Z3),levels=[1],colors="blue")
contour(X,Y,abs(Z4),levels=[1],colors="red")
axis([-0.3, 0.3, -3, 3])
xlabel('real(z)'), ylabel('imag(z)')
title('Stability domain of RK2, RK3, RK4')
grid(True);
```

RK2 does not enclose any portion of the imaginary axis, it is tangential to it, while RK3 and RK4 enclose a part of the imaginary axis. In case of RK4, the  interval $[ -\ii 2\sqrt{2}, \ii 2\sqrt{2}]$ is inside the stability domain.
:::
