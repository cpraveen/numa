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

# Simpson rule

```{include} math.md
```

## Single interval

Trapezoidal method used a linear approximation to estimate the integral.
To improve the accuracy, let us use a quadratic approximation. For the
function $f : [a,b] \to \re$, we first construct a quadratic polynomial
$p_2(x)$ that interpolates at $a, b, c = \half(a+b)$

$$
p_2(a) = f(a), \qquad p_2(b) = f(b), \qquad p_2(c) = f(c)
$$ 

and the
approximation to the integral is

$$
I_2(f) = \int_a^b p_2(x) \ud x = \frac{h}{3}[f(a) + 4 f(c) + f(b)], \qquad h = \half(b-a)
$$

This is called Simpson rule. The error is

$$
E_2(f) = I(f) - I_2(f) = \int_a^b \underbrace{(x-a)(x-c)(x-b)}_{\textrm{changes sign at $x=c$}} f[a,b,c,x] \ud x
$$

Define 

$$
w(x) = \int_a^x (t-a)(t-c)(t-b)\ud x
$$ 

Then

$$
w(a) = w(b) = 0, \qquad w(x) > 0 \quad \forall x \in (a,b)
$$ 

and

$$
w'(x) = (x-a)(x-c)(x-b)
$$ 

Hence the error can be written as

$$
\begin{aligned}
E_2(f) 
&= \int_a^b w'(x) f[a,b,c,x] \ud x \\
&= [w(x) f[a,b,c,x]]_a^b - \int_a^b w(x) \dd{}{x} f[a,b,c,x] \ud x \\
&= - \int_a^b w(x) f[a,b,c,x,x] \ud x
\end{aligned}
$$ 

Now we can apply the integral mean value theorem

$$
\begin{aligned}
E_2(f) 
&= -f[a,b,c,\xi,\xi] \int_a^b w(x) \ud x, \quad \textrm{for some $\xi \in [a,b]$} \\
&= -\frac{1}{24} f^{(4)}(\eta) \left[ \frac{4}{15} h^5 \right], \quad \textrm{for some $\eta \in [a,b]$}
\end{aligned}
$$ 

The error in Simpson rule is given by

$$
E_3(f) = -\frac{1}{90} h^5 f^{(4)}(\eta), \quad \textrm{for some $\eta \in [a,b]$}
$$

By construction, Simpson rule is exact for quadratic polynomials. The
error formula tells us that even for cubic polynomials, the error is
zero.

## Composite Simpson rule

Let $n \ge 2$ be even, the spacing between points $h = (b-a)/n$ and the
$n+1$ points 

$$
x_j = a + j h, \quad j=0,1,2,\ldots,n
$$ 

Then

$$
\begin{aligned}
I(f) 
&= \int_a^b f(x) \ud x \\
&= \sum_{j=1}^{n/2} \int_{x_{2j-2}}^{x_{2j}} f(x) \ud x \\
&= \sum_{j=1}^{n/2}\left\{ \frac{h}{3}[f_{2j-2} + 4 f_{2j-1} + f_{2j}] - \frac{h^5}{90} f^{(4)}(\eta_j) \right\}, \quad \eta_j \in [x_{2j-2},x_{2j}]
\end{aligned}
$$ 

The composite Simpson rule is

$$
I_n(f) = \frac{h}{3}[ f_0 + 4 f_1 + 2 f_2 + \ldots + 2 f_{n-2} + 4 f_{n-1} + f_n]
$$

The error in this approximation is

$$
E_n(f) = I(f) - I_n(f) = -\frac{h^5}{90} \left(\frac{n}{2}\right) \left[ \frac{2}{n} \sum_{j=1}^{n/2} f^{(4)}(\eta_j) \right]
$$

If the fourth derivative is continuous, then

$$
E_n(f) = -\frac{1}{180} (b-a)h^4 f^{(4)}(\eta), \quad \textrm{for some $\eta \in [a,b]$}
$$

We can also derive the asymptotic error formula

$$
E_n(f) \approx \tilde E_n(f) = -\frac{h^4}{180}[ f^{(3)}(b) - f^{(3)}(a)]
$$

:::{prf:example}

$$
I = \int_0^\pi \ee^x \cos x \ud x
$$

:::

:::{prf:example}

$$
I  = \int_0^1 x^3 \sqrt{x} \ud x = \frac{2}{9}
$$

:::

:::{prf:example}

$$
I = \int_0^5 \frac{\ud x}{1 + (x-\pi)^2} = \tan^{-1}(5 - \pi) + \tan^{-1}(\pi)
$$

:::

:::{prf:example}

$$
I = \int_0^1 \sqrt{x} \ud x = \frac{2}{3}
$$

:::

:::{prf:example}

$$
I = \int_0^{2\pi} \ee^{\cos x} \ud x = 7.95492652101284
$$

:::
