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

(chap:hermint)=
# Hermite interpolation

```{include} math.md
```

## Using derivatives

So far we have developed approximation methods which make use of function values only. Sometimes, the derivative information may also be available. Methods which make use of derivatives also to construct the approximation are known are Hermite approximations. Here we will look at interpolating given function and derivative information.

:::{prf:example} Cubic hermite interpolation
Given $f(a)$, $f'(a)$, $f(b)$, $f'(b)$, we can try to construct a cubic polynomial $p(x)$ so that 

$$
\begin{aligned}
p(a) = f(a), \quad p'(a) = f'(a) \\
p(b) = f(b), \quad p'(b) = f'(b)
\end{aligned}
$$ 

The cubic polynomial satisfying the above conditions is

$$
\begin{aligned}
H_2(x) &= \left[1 + 2\frac{x-a}{b-a}\right] \left[ \frac{b-x}{b-a} \right]^2 f(a) +  \left[ 1 + 2\frac{b-x}{b-a}\right] \left[\frac{x-a}{b-a}\right]^2 f(b) \\
& \quad +  \frac{(x-a)(x-b)^2}{(b-a)^2} f'(a) - \frac{(x-a)^2 (b-x)}{(b-a)^2} f'(b) \\
&= \phi_1(x) f(a) + \phi_2(x) f(b) + \phi_3(x) f'(a) + \phi_4(x) f'(b)
\end{aligned}
$$ 

This is analogous to a Lagrange interpolation with some set of basis functions which are multiplying the given data.

We can write the Hermite interpolation in terms of Newton form using the points

$$
\{x_0, x_1, x_2, x_3 \} = \{a, a, b, b\}
$$

which is

$$
H_2(x) &= f(a) + (x-a)f[a,a] + (x-a)^2 f[a,a,b] \\
& \quad + (x-a)^2 (x-b) f[a,a,b,b]
$$

where the divided differences are given by

\begin{align}
f[a,a] &= \lim_{t \to a} f[t,a] = \lim_{t \to a} \frac{f(t) - f(a)}{t-a} = f'(a) \\
f[a,a,b] &= \frac{f[a,b] - f[a,a]}{b-a} = \frac{f[a,b] - f'(a)}{b-a} \\
f[a,a,b,b] &= \frac{f[a,b,b] - f[a,a,b]}{b-a} = \frac{\frac{f[a,b] - f'(b)}{a-b} - \frac{f[a,b] - f'(a)}{b-a}}{b-a} \\
&= \frac{f'(b) - 2f[a,b] + f'(a)}{(b-a)^2}
\end{align}

It is easy to see that $H_2(a) = f(a)$ and

$$
H_2(b) &= f(a) + (b-a)f'(a) + (b-a)^2 \frac{f[a,b] - f'(a)}{b-a} \\
&= f(a) + (b-a)f[a,b] \\
&= f(b)
$$

Moreover

$$
H_2'(x) = f[a,a] + 2(x-a)f[a,a,b] + \{ 2(x-a)(x-b) + (x-a)^2 \} f[a,a,b,b,]
$$

so that $H_2'(a)  = f[a,a] = f'(a)$ and 

$$
\begin{aligned}
H_2'(b) &= f[a,a] + 2(b-a)f[a,a,b] + (b-a)^2 f[a,a,b,b] \\
&= f'(a) + 2\{ f[a,b] - f'(a) \} + \{ f'(b) - 2f[a,b] + f'(a) \} \\
&= f'(b)
\end{aligned}
$$

This form cannot be used for numerical computation, since some of the divided differences require taking limits.
:::

## Lagrange form of interpolation

Given $(x_i,y_i,y_i')$, $i=0,1,2,\ldots,N$, find a polynomial $p(x)$
such that

$$
p(x_i) = y_i = f(x_i), \qquad p'(x_i) = y_i' = f'(x_i), \qquad i=0,1,2,\ldots,N
$$

There are $2N+2$ conditions given so we can try to find a polynomial of
degree $2N+1$ that interpolates the given data. 

Define

$$
\ell(x) = (x-x_0)\ldots(x-x_N)
$$ 

then the Lagrange polynomials are

$$
\ell_i(x) = \prod_{j=0, j \ne i}^N \frac{x-x_j}{x_i - x_j} = \frac{\ell(x)}{(x-x_i) \ell'(x_i)}
$$ 

Define 

$$
\begin{aligned}
r_i(x) &= [1 - 2\ell_i'(x_i)(x-x_i)][\ell_i(x)]^2\\
s_i(x) &= (x-x_i) [\ell_i(x)]^2
\end{aligned}
$$ 

Since $\ell_i \in \poly_N$ both $r_i$, $s_i$ are in $\poly_{2N+1}$. Since $\ell_i(x_j) = \delta_{ij}$, we have 

\begin{align}
r_i(x_j)  &= \delta_{ij}, \qquad & s_i(x_j) &= 0 \\
r_i'(x_j) &= 0, \qquad & s_i'(x_j) &= \delta_{ij}
\end{align}

Then the interpolating polynomial is given by

$$
H_N(x) = \sum_{i=0}^N y_i r_i(x) + \sum_{i=0}^N y_i' s_i(x)
$$

##### Uniqueness

Suppose there are two polynomials $H(x), G(x) \in \poly_{2N+1}$ that interpolate the same set of data. Then define

$$
R(x) = H(x) - G(x) \in \poly_{2N+1}
$$ 

and

$$
R(x_i) = R'(x_i) = 0, \qquad i=0,1,\ldots,N
$$ 

Thus $R(x)$ has $N+1$ double roots, hence atleast $2N+2$ roots, which implies that $R(x) \equiv 0$.

## Newton form of interpolation

Even though the Newton form is not useful for computations, we show it now for general degree. The polynomial interpolating $f(x)$ at nodes $z_0,\ldots,z_{2N+1}$ is

\begin{align}
p_{2N+1}(x) &= f(z_0) \\
& \quad + (x-z_0)f[z_0,z_1] \\
& \quad + (x-z_0)(x - z_1) f[z_0,z_1,z_2] \\
& \quad + \ldots \\
& \quad + (x-z_0)\ldots(x- z_{2N})f[z_0,\ldots,z_{2N+1}]
\end{align}

with error formula

$$
f(x) - p_{2N+1}(x) = (x-z_0)\ldots(x-z_{2N})(x-z_{2N+1})f[z_0,\ldots,z_{2N+1},x]
$$

We can let the nodes coincide and Newtons formula will still make sense.

$$
\begin{aligned}
z_0, z_1 \to x_0 \\
z_2, z_3 \to x_1 \\
\vdots \qquad \\
z_{2N}, z_{2N+1} \to x_N
\end{aligned}
$$ 

which gives 

\begin{align}
p_{2N+1}(x) &= f(x_0) \\
& \quad + (x-x_0)f[x_0,x_0] \\
& \quad + (x-x_0)^2 f[x_0,x_0,x_1] \\
& \quad + \ldots \\
& \quad + (x-x_0)^2 \ldots(x-x_{N-1})^2 (x-x_N) f[x_0,x_0,\ldots,x_N,x_N]
\end{align}

This is a polynomial of degree $\le 2N+1$ with error
formula

$$
f(x) - p_{2N+1}(x) = (x-x_0)^2 \ldots (x-x_N)^2 f[x_0,x_0,\ldots,x_N,x_N,x]
$$

##### Claim: $p_{2N+1} = H_N$

Assume that $f \in \cts^{2N+3}$. Due to the Newton form, we have

$$
p_{2N+1}(x_i) = f(x_i), \qquad i=0,1,\ldots,N
$$ 

Differentiating the
error formula 

$$
\begin{aligned}
f'(x) - p_{2N+1}'(x) = \ & (x-x_0)^2 \ldots (x-x_N)^2 \dd{}{x}
f[x_0,x_0,\ldots,x_N,x_N,x] \\
& +  f[x_0,x_0,\ldots,x_N,x_N,x] \dd{}{x} [(x-x_0)^2 \ldots (x-x_N)^2] \\
= \ & (x-x_0)^2 \ldots (x-x_N)^2  f[x_0,x_0,\ldots,x_N,x_N,x,x] \\
& +  2 f[x_0,x_0,\ldots,x_N,x_N,x] \sum_{i=0}^N (x-x_i) \prod_{j=0,j \ne i}^N(x-
x_j)^2
\end{aligned}
$$ 

we see that

$$
f'(x_i) - p_{2N+1}'(x_i) = 0, \qquad i=0,1,\ldots,N
$$ 

We have degree
$p_{2N+1} \le 2N+1$ and it satisfies the given data $\{x_i,y_i,y_i'\}$,
$i=0,1,\ldots,N$. By uniqueness of the Hermite interpolating polynomial,
we have $p_{2N+1} = H_N$.
