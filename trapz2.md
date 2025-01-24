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

# Trapezoidal Rule - II

```{include} math.md
```

## Bernoulli polynomials

The Bernoulli polynomial $B_n(x)$, $n \ge 0$ is defined implicitly by the generating function

$$
\frac{t(\ee^{xt}-1)}{\ee^t - 1} = \sum_{j=1}^\infty B_j(x) \frac{t^j}{j!}
$$

The first few polynomials are[^1]

$$
B_0(x)=1, \qquad B_1(x)=x, \qquad B_2(x) = x^2 - x
$$

$$
B_3(x) = x^3 - \frac{3}{2} x^2 + \half x, \qquad B_4(x) = x^2(1-x)^2
$$

Note that 

$$
B_k(0)=0, \qquad k \ge 1
$$ 

It can be shown that for all $j \ge 0$ 

$$
B_{2j}(x) > 0, \qquad 0 < x < 1
$$

## Bernoulli numbers

These are defined implicitly by

$$
\frac{t}{\ee^t - 1} = \sum_{j=0}^\infty B_j \frac{t^j}{j!}
$$ 

The first few numbers are

$$
B_0 = 1, \qquad B_1 = -\half, \qquad B_2 = \frac{1}{6}
$$

$$
B_4 = -\frac{1}{30}, \qquad B_6 = \frac{1}{42}, \qquad B_8 = -\frac{1}{30}
$$

$$
B_j = 0, \qquad \textrm{$j$ is odd and $j \ge 3$}
$$

## Connection

To obtain the relation between $B_j(x)$ and $B_j$, integrate the generating function of Bernoulli polynomial on $[0,1]$ to obtain

$$
1 - \frac{t}{\ee^t - 1} = \sum_{j=1}^\infty \frac{t^j}{j!} \int_0^1 B_j(x) \ud x
$$

Hence 

$$
B_0 = 1, \qquad B_j = -\int_0^1 B_j(x) \ud x, \quad j \ge 1
$$

Next, define the periodic extension of $B_j(x)$ by

$$
\bar{B}_j(x) = \begin{cases}
B_j(x) & 0 \le x \le 1 \\
\bar{B}_j(x-1) & x > 1
\end{cases}
$$

:::{prf:theorem}
Let $m \ge 0$, $n \ge 1$ and define $h = (b-a)/n$, $x_j = a + jh$, $j=0,1,\ldots,n$. Further assume that $f \in \cts^{2m+2}[a,b]$ for some $m \ge 0$. Then the error in trapezoidal method is 

$$
\begin{aligned}
E_n(f) 
&= \int_a^b f(x) \ud x - h {\sum_{j=0}^n}' f(x_j) \\
&= -\sum_{i=1}^m \frac{B_{2i} h^{2i}}{(2i)!}[ f^{(2i-1)}(b) - f^{(2i-1)}(a)] \\
& \quad + \frac{h^{2m+2}}{(2m+2)!} \int_a^b \bar{B}_{2m+2}((x-a)/h) f^{(2m+2)}(x) \ud x
\end{aligned}
$$ 

where the prime on the summation indicates that the first and last terms must be halved.
:::

:::{prf:corollary}
Suppose $f \in \cts^\infty[a,b]$ and all its odd ordered derivatives are periodic with $b-a$ an integer multiple of the period. Then the order of convergence of the trapezoidal method is greater than any power of $h$.
:::

:::{prf:proof}
The periodicity of $f$ means

$$
f^{(2j-1)}(a) = f^{(2j-1)}(b), \qquad j \ge 1
$$ 

Hence for any $m$ the error formula of trapezoidal method yields 

$$
\begin{aligned}
E_n(f) &=& \frac{h^{2m+2}}{(2m+2)!} \int_a^b \bar{B}_{2m+2}((x-a)/h) f^{(2m+2)}(x) \ud x \\
&=& - \frac{h^{2m+2}(b-a)B_{2m+2}}{(2m+2)!} f^{(2m+2)}(\xi), \qquad \xi \in [a,b]
\end{aligned}
$$ 

where we used the positivity of $B_{2j}(x)$ and integral mean value theorem. Since $m$ was arbitrary, we obtain the desired
result.
:::

:::{prf:remark}
For periodic functions, trapezoidal rule is the most accurate, even better than the Gauss rules.
:::

:::{seealso}
[@Fornberg2021]
:::

[^1]: $B_0(x)$ is defined to be one, it is not present in the definition
    of the generating function.
