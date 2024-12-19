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

# Newton-Cotes integration

```{include} math.md
```

The Trapezoidal and Simpson rules were based on integrating a polynomial approximation. Such methods are called Newton-Cotes integration methods.  If we choose $n+1$ distinct points in $[a,b]$ and approximate the function $f : [a,b] \to \re$ by interpolation

$$
p_n(x) = \sum_{j=0}^n f(x_j) \ell_j(x), \qquad x \in [a,b]
$$ 

the integral can be approximated by

$$
I_n(f) = \int_a^b p_n(x) \ud x = \sum_{j=0}^n f(x_j) \int_a^b \ell_j(x) \ud x = \sum_{j=0}^n w_j f(x_j)
$$

where the weights are given by 

$$
w_j = \int_a^b \ell_j(x) \ud x
$$ 

As an example for $n=3$ and using uniformly spaced points, we interpolate by a cubic polynomial and the integral is approximated as

$$
I_3(f) = \frac{3h}{8}[f_0 + 3 f_1 + 3 f_2 + f_3]
$$ 

For general $n$ we can state the following error estimates.

:::{prf:theorem}
\(a\) For $n$ even, $f \in \cts^{n+2}[a,b]$

$$
I(f) - I_n(f) = C_n h^{n+3} f^{(n+2)}(\eta), \qquad \eta \in [a,b]
$$

$$
C_n = \frac{1}{(n+2)!} \int_0^n \mu^2 (\mu-1) \ldots (\mu-n) \ud \mu
$$

(b) For $n$ odd and $f \in \cts^{n+1}[a,b]$

$$
I(f) - I_n(f) = C_n h^{n+2} f^{(n+1)}(\eta), \qquad \eta \in [a,b]
$$

$$
C_n = \frac{1}{(n+1)!} \int_0^n \mu (\mu-1) \ldots (\mu-n) \ud \mu
$$
:::

:::{prf:proof}
The proof is based on using the interpolation error estimate and integral mean value theorem as was done for trapezoidal and Simpson rules.
:::

:::{prf:definition}
A numerical integration formula $\tilde{I}(f)$ that approximates $I(f)$ is said to have degree of precision $m$ if

1.  $\tilde{I}(f) = I(f)$ for all $f \in \poly_m$

2.  $\tilde{I}(f) \ne I(f)$ for atleast one $f \in \poly_{m+1}$
:::

:::{prf:remark}
Trapezoidal rule ($n=1$) has degree of precision one and Simpson rule ($n=2$) has degree of precision three. For $n$ odd, Newton-Cotes has degree of precision of $n$ while for $n$ even, it has degree of precision $n+1$.
:::

:::{prf:example}
On uniformly spaced nodes, the interpolating polynomial $p_n(x)$ does not always converge. Hence the Newton-Cotes formulae may not converge also. Consider

$$
I = \int_{-4}^4 \frac{\ud x}{1 + x^2} = 2 \tan^{-1}(4) \approx 2.6516
$$

Hence it is not advisable to use Newton-Cotes rule for large degree $n$.  However, we can use moderate but fixed degree of the interpolating polynomial and build a composite rule. The convergence here is obtained as the length of the sub-intervals tends to zero but the polynomial degree in each sub-interval remains fixed.
:::

:::{prf:definition}
Let $\mathcal{F}$ be a set of continuous functions on a given interval $[a,b]$. We say that $\mathcal{F}$ is dense in $[a,b]$ if for every $f \in \cts[a,b]$ and every $\epsilon > 0$, $\exists f_\epsilon \in \mathcal{F}$ such that

$$
\norm{f - f_\epsilon}_\infty \le \epsilon
$$
:::

:::{prf:example}

(1) By Weirstrass theorem, the set of all polynomials is dense in $\cts[a,b]$. 

(2) Define

$$
\mathcal{F} = \{ f \in \cts[a,b] : \exists \textrm{ a partition $(x_j)_j$ of $[a,b]$ s.t. } f|_{[x_{j-1},x_j]} \in \poly_1 \}
$$

then $\mathcal{F}$ is dense in $\cts[a,b]$.
:::

:::{prf:theorem}
Let 

$$
I_n(f) = \sum_{j=0}^n w_{j,n} f(x_{j,n}), \qquad n \ge 1
$$ 

be a sequence of numerical integration formulae that approximate

$$
I(f) = \int_a^b f(x) \ud x
$$ 

Let $\mathcal{F}$ be a family dense in $\cts[a,b]$. Then 

$$
I_n(f) \to I(f) \qquad \forall f \in \cts[a,b]
$$ 

if
and only if

-   $I_n(f) \to I(f)$ for all $f \in \mathcal{F}$

-   $B = \sup_n \sum_{j=0}^n |w_{j,n}| < \infty$
:::

:::{prf:proof}
(a) The forward implication is an application of the uniform boundedness principle. Consider 

$$
I_n : \mathcal{F} \to \re
$$ 

But why is $\mathcal{F}$ a Banach space ? The set of polynomials is not a Banach space in the maximum norm.

(b) Now we prove the reverse implication. We have to use the density and
the given two conditions.

$$
I(f) - I_n(f) = I(f) - I(f_\epsilon) + I(f_\epsilon)  - I_n(f_\epsilon) + I_n(f_\epsilon) - I_n(f)
$$

$$
\begin{aligned}
|I(f) - I_n(f)| 
&\le |I(f) - I(f_\epsilon)| + |I(f_\epsilon)  - I_n(f_\epsilon)| + |I_n(f_\epsilon) - I_n(f)| \\
&\le (b-a)\norm{f - f_\epsilon}_\infty + |I(f_\epsilon)  - I_n(f_\epsilon)| + \norm{f - f_\epsilon}_\infty \underbrace{\sum_{j=0}^n |w_{j,n}|}_{\le B}
\end{aligned}
$$ 

Using density, choose $f_\epsilon$ such that

$$
\norm{f - f_\epsilon}_\infty \le \frac{\epsilon}{2(b-a + B)}
$$ 

Using
property (1), we can choose $n$ large enough so that

$$
|I(f_\epsilon)  - I_n(f_\epsilon)| \le \frac{\epsilon}{2}
$$ 

Then

$$
|I(f) - I_n(f)| \le \epsilon
$$ 

Since $\epsilon$ was arbitrary, we conclude that $I_n(f) \to I(f)$.
:::

:::{prf:example}
For the function $f(x) = 1/(1+x^2)$, $x \in [-4,4]$, the Newton-Cotes on uniform points does not converge. This is because while property (1) in above theorem is satisfied, property (2) does not hold for uniform points.
:::
