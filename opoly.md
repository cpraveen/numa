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

# Orthogonal polynomials

```{include} math.md
```

We go back to orthogonal polynomials and establish some results that are useful for further topics on numerical integration. We recall two results. If $\{\phi_n\}$ are orthogonal polynomials wrt some inner product and degree $\phi_n = n$, and if $f \in \poly_n$ we can write it as

$$
f(x) = \sum_{j=0}^n c_j \phi_j(x), \qquad c_j = \frac{\ip{f,\phi_j}}{\ip{\phi_j,\phi_j}}
$$

## Roots of orthogonal polynomials

:::{prf:theorem}
Let $\{ \phi_n : n \ge 0\}$ be an orthogonal family of polynomials on $[a,b]$ with weight function $w(x)$. For $n \ge 1$, the polynomial $\phi_n$ has exactly $n$ distinct real roots in the open interval $(a,b)$.
:::

:::{prf:proof}
Let $x_1, x_2, \ldots, x_m$ be all the distinct zeros of $\phi_n$ for which

1.  $a < x_i < b$

2.  $\phi_n$ changes sign at $x_i$

i.e., there are $m$ roots each with odd multiplicity in $(a,b)$. This includes all the simple roots also. There might be other roots outside the interval $(a,b)$ and there could be roots with even multiplicity in $(a,b)$. Since degree $\phi_n = n$, we have $m \le n$. Now since $\phi_n$ is orthogonal to $\phi_0 =$ constant,

$$
\int_a^b w(x) \phi_n(x) \ud x = 0, \qquad n \ge 1
$$ 

$\implies$ each $\phi_n$ has atleast one root in $(a,b)$, i.e., $m \ge 1$.

Let us assume that $m < n$ and derive a contradiction. We can write

$$
\phi_n(x) = (x-x_1)^{r_1} \ldots (x-x_m)^{r_m} h(x)
$$ 

Each $r_i$ is odd and other roots are roots of $h(x)$. Define

$$
B(x) = (x-x_1) \ldots (x-x_m)
$$ 

so that

$$
B(x) \phi_n(x) = (x-x_1)^{r_1+1} \ldots (x-x_m)^{r_m+1} h(x)
$$ 

$h(x)$ does not change sign in $(a,b)$ and since each $r_i+1$ is even, $B(x)\phi_n(x)$ does not change sign in $(a,b)$. Hence

$$
\int_a^b w(x) \phi_n(x) B(x) \ud x = \ip{B, \phi_n} \ne 0
$$ 

But degree of $B = m < n$ so that $\ip{B,\phi_n} = 0$ which leads to a contradiction. Hence $m=n$. Since $\phi_n$ has atmost $n$ roots, the multiplicity of each $x_i$ is one. Hence $\phi_n$ has $n$ distinct roots.
:::

## Triple recursion relation

Let $\{ \phi_n \}$ be an orthogonal sequence of polynomials on $[a,b]$ with weight function $w$ and degree $\phi_n = n$ for all $n$. Define $A_n$ and $B_n$ by 

$$
\phi_n(x) = A_n x^n + B_n x^{n-1} + \ldots
$$ 

Also write 

$$
\phi_n(x) = A_n (x- x_{n,1}) \ldots (x-x_{n,n})
$$ 

Let

$$
a_n = \frac{A_{n+1}}{A_n}, \qquad \gamma_n = \ip{\phi_n, \phi_n} > 0
$$

:::{prf:theorem}
Let $\{ \phi_n \}$ be an orthonal family of polynomials on $[a,b]$ with weight function $w$. Then for $n \ge 1$

$$
\phi_{n+1}(x) = (a_n x + b_n) \phi_n(x) - c_n \phi_{n-1}(x)
$$ 

where

$$
b_n = a_n \left[ \frac{B_{n+1}}{A_{n+1}} - \frac{B_n}{A_n} \right], \qquad c_n = \frac{A_{n-1} A_{n+1}}{A_n^2} \frac{\gamma_n}{\gamma_{n-1}}
$$
:::

:::{prf:proof}
(a) Consider the polynomial 

$$
\begin{aligned}
G(x) 
&= \phi_{n+1}(x) - a_n x \phi_n(x) \\
&= [A_{n+1} x^{n+1} + B_{n+1} x^n + \ldots] - \frac{A_{n+1}}{A_n} x [A_n x^n + B_n x^{n-1} + \ldots ] \\
&= \left[ B_{n+1} - \frac{A_{n+1} B_n}{A_n} \right] x^n + \ldots
\end{aligned}
$$ 

and hence degree $G \le n$. We can write $G$ as

$$
G(x) = d_0 \phi_0 + \ldots + d_n \phi_n(x)
$$ 

where

$$
d_i = \frac{\ip{G,\phi_i}}{\ip{\phi_i,\phi_i}} = \frac{1}{\gamma_i}[ \ip{\phi_{n+1}, \phi_i} - a_n \ip{x \phi_n, \phi_i} ], \qquad 0 \le i \le n
$$

But

$$
\ip{\phi_{n+1},\phi_i} &= 0, \qquad 0 \le i \le n \\
\ip{x\phi_n, \phi_i} &= \ip{\phi_n, x \phi_i} = 0, \qquad i+1 < n \quad i.e., \quad 0 \le i \le n-2 \\
\implies d_i &= 0, \qquad 0 \le i \le n-2
$$ 

so that

$$
G(x) = d_n \phi_n(x) + d_{n-1} \phi_{n-1}(x)
$$ 

or

$$
\phi_{n+1}(x) = (a_n x + d_n) \phi_n(x) + d_{n-1} \phi_{n-1}(x)
$$ 

We have a triple recursion formula and it remains to show that $d_{n-1} = -c_n$ and $d_n = b_n$ where $b_n, c_n$ are as defined in the theorem.

(b) $d_{n-1}$ is given by

$$
d_{n-1} = - \frac{a_n}{\gamma_{n-1}} \ip{\phi_n, x\phi_{n-1}}
$$ 

Now

$$
\begin{aligned}
x\phi_{n-1}(x) 
&= A_{n-1} x^n + B_{n-1} x^{n-1} + \ldots \\
&= A_{n-1} \frac{1}{A_n}[\phi_n(x) - B_n x^{n-1} - \ldots] + B_{n-1} x^{n-1} + \ldots
\end{aligned}
$$ 

so that

$$
\ip{\phi_n, x\phi_{n-1}} = \frac{A_{n-1}}{A_n} \ip{\phi_n, \phi_n} = \frac{A_{n-1}\gamma_n}{A_n}
$$

and hence

$$
d_{n-1} = -\frac{a_n}{\gamma_{n-1}} \frac{A_{n-1}\gamma_n}{A_n} = -\frac{A_{n-1}A_{n+1}}{A_n^2} \frac{\gamma_n}{\gamma_{n-1}} = -c_n
$$

(c) $d_n$ is given by

$$
d_{n} = - \frac{a_n}{\gamma_{n}} \ip{\phi_n, x\phi_{n}}
$$ 

Now

$$
\begin{aligned}
x\phi_n(x) 
&=  A_n x^{n+1} + B_n x^n + \ldots \\
&= A_n \frac{1}{A_{n+1}}[ \phi_{n+1}(x) - B_{n+1} x^n - \ldots ] + B_n \frac{1}{A_n}[ \phi_n(x) - B_n x^{n-1} - \ldots ] \\
& \quad + ()x^{n-1} + \ldots
\end{aligned}
$$ 

so that

$$\ip{\phi_n ,x\phi_n} = - \frac{A_n B_{n+1}}{A_{n+1}} \ip{\phi_n, x^n} + \frac{B_n}{A_n} \underbrace{\ip{\phi_n, \phi_n}}_{=\gamma_n}
$$

But

$$
\ip{\phi_n, x^n} &= \ip{\phi_n, \frac{1}{A_n}[\phi_n - B_n x^{n-1} - \ldots]} = \frac{1}{A_n} \ip{\phi_n,\phi_n} = \frac{\gamma_n}{A_n} \\
\implies \ip{\phi_n ,x\phi_n} &= -\frac{B_{n+1} \gamma_n}{A_{n+1}} + \frac{B_n \gamma_n}{A_n} \\
\implies d_n &= a_n \left[ \frac{B_{n+1}}{A_{n+1}} - \frac{B_n}{A_n} \right] = b_n
$$
:::

(sec:chrisdarb)=
## Christoffel-Darboux identity

We will use the next identity to derive simpler formula for the quadrature weights.

:::{prf:theorem} Christoffel-Darboux identity
For $\{ \phi_n \}$ an orthogonal family of polynomials with weight
function $w$

$$
\sum_{k=0}^n \frac{\phi_k(x) \phi_k(y)}{\gamma_k} = \frac{\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y)}{a_n \gamma_n (x-y)}, \qquad x \ne y
$$
:::

:::{prf:proof}
Using the triple recursion relation, we can write

$$
\begin{aligned}
\phi_n(y) \cdot [\phi_{n+1}(x) &= (a_n x + b_n) \phi_n(x) - c_n \phi_{n-1}(x) ] \\
\phi_n(x) \cdot [\phi_{n+1}(y) &= (a_n y + b_n) \phi_n(y) - c_n \phi_{n-1}(y) ]
\end{aligned}
$$ 

Subtracting

$$
\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y) &= a_n (x-y) \phi_n(x) \phi_n(y) \\
& \quad + c_n[ \phi_n(x) \phi_{n-1}(y) - \phi_{n-1}(x) \phi_n(y) ]
$$

and

$$
\begin{aligned}
& \frac{\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y)}{a_n \gamma_n (x-y)} \\
&= \frac{\phi_n(x) \phi_n(x)}{\gamma_n} + \frac{c_n}{a_n \gamma_n} \frac{ \phi_n(x) \phi_{n-1}(y) - \phi_{n-1}(x) \phi_n(y)  }{x-y} \\
&= \frac{\phi_n(x) \phi_n(x)}{\gamma_n} + \frac{ \phi_n(x) \phi_{n-1}(y) - \phi_{n-1}(x) \phi_n(y)  }{a_{n-1} \gamma_{n-1} (x-y)}
\end{aligned}
$$ 

since

$$
\frac{c_n}{a_n \gamma_n} = \frac{1}{a_{n-1} \gamma_{n-1}}
$$ 

Repeating this process, we get

$$
\frac{\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y)}{a_n \gamma_n (x-y)} = \sum_{k=1}^n \frac{\phi_k(x) \phi_k(y)}{\gamma_k} + \frac{ \phi_1(x) \phi_{0}(y) - \phi_{0}(x) \phi_1(y)  }{a_{0} \gamma_{0} (x-y)}
$$

But

$$
\frac{ \phi_1(x) \phi_{0}(y) - \phi_{0}(x) \phi_1(y)  }{a_{0} \gamma_{0} (x-y)} 
&= \frac{(A_1 x + B_1) A_0 - A_0 (A_1 y + B_1)}{(A_1/A_0) \gamma_0 (x-y)} \\
&= \frac{A_0^2}{\gamma_0} \\
&= \frac{\phi_0(x) \phi_0(y)}{\gamma_0}
$$

which completes the proof.
:::
