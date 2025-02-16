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

# Mid-point rule

```{include} math.md
```

## Single interval

The trapezoidal rule approximates the integral as the area below the
straight line approximation

$$
\int_a^b f(x) \ud x \approx (b-a) \frac{f(a) + f(b)}{2}
$$ 

Another
approximation to the area is to use the function value at the mid-point

$$
\int_a^b f(x) \ud x \approx (b-a) f\left( \frac{a+b}{2} \right)
$$ 

Let $c=(a+b)/2$. Using Taylor formula

$$
f(x) = f(c) + (x-c) f'(c) + \half (x-c)^2 f''(\xi_x), \qquad \textrm{$\xi_x$ between $x$ and $c$}
$$

Integrating both sides

$$
\int_a^b f(x) \ud x - (b-a) f\left( \frac{a+b}{2} \right) = \half \int_a^b (x-c)^2 f''(\xi_x) \ud x
$$

Using integral mean value theorem, the error is

$$
\int_a^b f(x) \ud x - (b-a) f\left( \frac{a+b}{2} \right) 
&= \half f''(\eta) \int_a^b (x-c)^2  \ud x \\
&= \frac{(b-a)^3}{24} f''(\eta), \qquad \eta \in [a,b]
$$

## Composite rule

Let us partition $[a,b]$ into $n$ intervals each of width $h = (b-a)/n$.  The mid-points of these intervals are

$$
x_j = a + (j-\shalf)h, \qquad j=0, 1, \ldots, n-1
$$ 

The integral is given by 

$$
\begin{aligned}
\int_a^b f(x) \ud x 
&= \sum_{j=1}^{n-1} \int_{x_{j-1}}^{x_j} f(x) \ud x \\
&= \sum_{j=1}^{n-1} h f(x_j) + \frac{h^3}{24} \sum_{j=1}^n f''(\eta_j), \qquad \eta_j \in [x_{j-1}, x_j] \\
&= I_n(f) + E_n(f)
\end{aligned}
$$ 

where 

$$
I_n(x) = h[f_0 + f_1 + \ldots + f_{n-1}]
$$ 

and the error is 

$$
\begin{aligned}
E_n(f) 
&= \frac{h^3}{24} \sum_{j=1}^n f''(\eta_j) = \frac{n h^3}{24} \left[ \frac{1}{n} \sum_{j=1}^n f''(\eta_j) \right] \\
&= \frac{(b-a)h^2}{24} f''(\eta), \qquad \eta \in [a,b]
\end{aligned}
$$

:::{prf:remark}
The error of composite trapezoidal rule is

$$
E_n(f) = \frac{(b-a)h^2}{12} f''(\eta), \qquad \eta \in [a,b]
$$ 

Thus the error of the mid-point rule can be half the error of the trapezoidal rule. Anyway, both methods require similar amount of information and computational work.
:::

## Integral error estimate

Using Taylor formula with integral remainder term

$$
f(x) = f(c) + (x-c)f'(c) + \int_c^x (x-t) f''(t) \ud t
$$ 

The error in mid-point rule is

$$
\int_a^b f(x) \ud x - (b-a) f(c) =  \int_a^b \int_c^x (x-t) f''(t) \ud t \ud x
$$

By change of order of integration 

$$
\begin{aligned}
\int_a^b \int_c^x (x-t) f''(t) \ud t \ud x 
&= \int_c^b \int_t^b (x-t) f''(t) \ud x \ud t - \int_a^c \int_a^t (x-t) f''(t) \ud x \ud t \\
&= \int_c^b \frac{(b-t)^2}{2} f''(t) \ud t + \int_a^c \frac{(a-t)^2}{2} f''(t) \ud t
\end{aligned}
$$ 

Using this estimate, the error of the composite rule is

$$
|E_n(f)| \le \frac{h^4}{4} \int_a^b |f''(t)| \ud t
$$ 

This error is more than that of trapezoidal rule, which seems strange.
