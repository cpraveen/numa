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
  headings: true
---

# Trigonometric interpolation

```{include} math.md
```

```{code-cell}
from pylab import *
```

The following function computes the trigonometric interpolant and plots it.

```{code-cell}
# Wave numbers are arranged as k=[0,1,...,N/2,-N/2+1,-N/2,...,-1]
def fourier_interp(N,f,ne=500,fig=True):
    if mod(N,2) != 0:
        print("N must be even")
        return
    h = 2*pi/N; x = h*arange(0,N);
    v = f(x);
    v_hat = fft(v)
    k = zeros(N)
    n = int(N/2)
    k[0:n+1] = arange(0,n+1)
    k[n+1:] = arange(-n+1,0,1)

    xx = linspace(0.0,2*pi,ne)
    vf = real(dot(exp(1j*outer(xx,k)), v_hat)/N)
    ve = f(xx)

    # Plot interpolant and exact function
    if fig:
        plt.plot(x,v,'o',xx,vf,xx,ve)
        plt.legend(('Data','Fourier','Exact'))
        plt.title("N = "+str(N))

    errori = abs(vf-ve).max()
    error2 = sqrt(h*sum((vf-ve)**2))
    print("Error (max,L2) = ",errori,error2)
    return errori, error2
```


## Infintely smooth, periodic function


```{code-cell}
f1 = lambda x: exp(sin(x))
fourier_interp(24,f1);
```


## Infinitely smooth, derivative not periodic

```{code-cell}
g = lambda x: sin(x/2)
e1i,e12 = fourier_interp(24,g)
```

```{code-cell}
e2i,e22 = fourier_interp(48,g,fig=False)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```


## Continuous function

```{code-cell}
def trihat(x):
    f = 0*x
    for i in range(len(x)):
        if x[i] < 0.5*pi or x[i] > 1.5*pi:
            f[i] = 0.0
        elif x[i] >= 0.5*pi and x[i] <= pi:
            f[i] = 2*x[i]/pi - 1
        else:
            f[i] = 3 - 2*x[i]/pi
    return f

e1i,e12 = fourier_interp(24,trihat)
```

```{code-cell}
e2i,e22 = fourier_interp(48,trihat,fig=False)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```

## Discontinuous function

```{code-cell}
f2 = lambda x: (abs(x-pi) < 0.5*pi)
e1i,e12 = fourier_interp(24,f2)
```

```{code-cell}
e2i,e22 = fourier_interp(48,f2)
print("Rate (max,L2) = ", log(e1i/e2i)/log(2), log(e12/e22)/log(2))
```
