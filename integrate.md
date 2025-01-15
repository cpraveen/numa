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

# Numerical Integration

```{include} math.md
```

The problem is to compute the numerical value of the definite integral

$$
I(f) = \int_a^b f(x) \ud x
$$ 

where $f:[a,b] \to \re$ is a given function. The solution to this problem usually takes the following form: Form a partition or grid 

$$
a \le x_0 \le x_1 \le \ldots \le x_n \le b
$$

and approximate the integral by a formula of the type

$$
I_n(f) = \sum_{j=0}^n w_j f(x_j)
$$ 

where the $\{ w_j \}$ are some weights. We may be interested in the following type of questions.

1.  Given the nodes $\{ x_j \}$ how to find the weights to obtain a good
    approximation ?

2.  Can we find both the nodes and weights so that the formula is as
    accurate as possible ?

To measure the accuracy, we may establish a result like

$$
|I(f) - I_n(f)| = \order{\frac{1}{n^p}}, \qquad \textrm{for some $p > 0$}
$$

Then the approximations converge fast if $p$ is very large and this is algebraic convergence. A still better approximation is one which would converge exponentially

$$
|I(f) - I_n(f)| = \order{\ee^{-\alpha n}}, \qquad \textrm{for some $\alpha > 0$}
$$

## Integration via function approximation

Let $\{ f_n \}$ be a family of approximations to $f$ such that

$$
\lim_{n \to \infty} \norm{f - f_n}_\infty = 0
$$ 

i.e., they converge uniformly and also in a pointwise sense. Then define 

$$
I_n(f) = \int_a^b f_n(x) \ud x = I(f_n)
$$

Then the error in the integral is

$$
E_n(f) = I(f) - I_n(f) = \int_a^b [f(x) - f_n(x)] \ud x
$$ 

so that

$$
|E_n(f)| \le \int_a^b |f(x) - f_n(x)| \ud x \le (b-a) \norm{f - f_n}_\infty
$$

We thus automatically obtain convergence of the integral approximations.

:::{prf:remark}

Usually, we take the approximations $f_n$ to be polynomials. Then clearly 

$$
I_n(p) = I(p), \qquad \forall p \in \poly_n
$$

:::

:::{prf:remark}
The approximation $f_n$ can be either

-   global: a single approximation in the whole interval $[a,b]$

-   local: partition $[a,b]$ into many sub-intervals and construct $f_n$
    in a piecewise manner in each sub-interval
:::
