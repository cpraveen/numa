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


# Newton interpolation

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

## The method

Consider $N+1$ distinct points $x_0,x_1,\ldots,x_N$ at which data about a function is given. The Newton form of interpolation makes use of the basis functions

$$
1, \quad (x-x_0), \quad (x-x_0)(x-x_1), \quad \ldots, \quad (x-x_0)(x-x_1)\ldots(x- x_{N-1})
$$ 

so that the degree $N$ polynomial has the form

$$
p(x) = c_0 + c_1(x-x_0) + c_2 (x-x_0)(x-x_1) + \ldots + c_N (x-x_0) \ldots(x-x_{N-1})
$$

The two practical questions are:

1.  How to determine the coefficients $\{ c_i \}$ ?

2.  How to evaluate the resulting polynomial ?

The interpolation conditions are given by 

$$
\begin{aligned}
f_0 &= c_0 \\
f_1 &= c_0 + c_1 (x_1 - x_0) \\
f_2 &= c_0 + c_1(x_2-x_0) + c_2(x_2-x_0)(x_2-x_1) \\
 & \vdots  \\
 f_N &= c_0 + c_1 (x_N-x_0) + \ldots + c_N(x_N-x_0)\ldots(x_N-x_{N-1})
\end{aligned}
$$ 

We have a matrix equation with a lower triangular matrix which can be solved recursively starting from the first equation. There is no need to invert the matrix. The non-singularity of the matrix is clear since the determinant is the product of the diagonal terms and this is non-zero since the points $\{ x_i \}$ are distinct. The coefficients $\{ c_i \}$ can be computed in $O(N^2)$ operations but there is danger of overflow/ underflow errors if we do not do this carefully. The triangular structure means that if we add a new data point $(x_{N+1},f_{N+1})$, then the coefficients $c_0, c_1, \cdots c_N$ do not change; we have to compute the additional coefficient $c_{N+1}$.

## Divided differences

Define zeroth divided difference 

$$
f[x_i] = f_i, \qquad i=0,1,\ldots,N
$$

and the first divided difference

$$
f[x_i,x_j] = \frac{f[x_j] - f[x_i]}{x_j - x_i}, \qquad i \ne j
$$ 

The $k$'th divided difference is defined recursively by

$$
f[x_0,x_1,\ldots,x_k] = \frac{ f[x_1,x_2,\ldots,x_k] - f[x_0,x_1,\ldots,x_{k-1}]}{x_k - x_0}
$$ 

We claim that 

$$
c_k = f[x_0,x_1,\ldots,x_k]
$$ 

We can compute the divided differences by constructing the following table

$$
\begin{array}{llll}
f[x_0] &  &   & \\
f[x_1] & f[x_0,x_1] &  & \\
f[x_2] & f[x_1,x_2] & f[x_0,x_1,x_2] & \\
f[x_3] & f[x_2,x_3] & f[x_1,x_2,x_3] & f[x_0,x_1,x_2,x_3] \\
\vdots & \vdots & \vdots & \vdots
\end{array}
$$ 

The diagonal terms are the coefficients $\{c_i\}$ of Newton interpolation. To compute the the diagonal terms, it is not necessary to store all the other terms in memory.

```c
for(i=0; i<=N; ++i) c[i] = f[i];
for(j=1; j<=N; ++j)
   for(i=N; i>=j; --i)
      c[i] = (c[i] - c[i-1])/(x[i] - x[i-j]);
```

This requires $N^2$ additions and $\frac{1}{2}N^2$ divisions.

## Evaluation of polynomials

Suppose given $x$ we have to compute

$$
p(x) = a_N x^N + a_{N-1} x^{N-1} + \ldots + a_1 x + a_0
$$ 

A direct term-by-term evaluation as in the above formula can suffer from overflow/ underflow error and the cost is also high. To see how to evaluate more efficiently, let us consider some simpler examples for quadratic 

$$
p = a_2 x^2 + a_1 x + a_0 = (a_2 x + a_1) x + a_0
$$

and cubic polynomials,

$$
p = a_3 x^3 + a_2 x^2 + a_1 x + a_0 = ((a_3 x + a_2) x + a_1) x + a_0
$$

This suggests the following recursive algorithm called *Horner's method*

$$
\begin{aligned}
p &= a_N \\
p &= p \cdot x + a_{N-1} \\
p &= p \cdot x + a_{N-2} \\
   \vdots & \\
p &= p \cdot x + a_0
\end{aligned}
$$ 

which gives $p(x)$ at the end. This approach can be used for the Newton form also. 

$$
\begin{aligned}
p &= c_N \\
p &= p \cdot (x-x_{N-1}) + c_{N-1} \\
p &= p \cdot (x-x_{N-2}) + c_{N-2} \\
  & \vdots \\
p &= p \cdot (x-x_0) + c_0
\end{aligned}
$$ 

This requires $2N$ additions and $N$ multiplications so the cost is $O(N)$. The direct term-by-term computation

```c
double x[N], c[N];
double p = 0.0;
for(int i=0; i<=N; ++i)
{
   double tmp = 1.0;
   for(int j=0; j<i; ++j) tmp *= (x - x[j]);
   p += c[i] * tmp;
}
```

would require $N$ additions and

$$
1 + 2 + \ldots + N = \frac{1}{2}N(N+1)
$$ 

multiplications for a total cost of $O(N^2)$ for each evaluation of $p(x)$. This can be improved by a simple modification

```c
double x[N], c[N];
double tmp = 1.0, p = 0.0;
for(int i=0; i<=N; ++i)
{
   p += c[i] * tmp;
   tmp *= (x-x[i]);
}
```

:::{exercise}
Give the data $(x_i,f_i)$, $i=0,1,\ldots,N$ write a function to compute the coefficients of Newton interpolation. Write a second function which evaluates $p(x)$ by Horner's method.
:::

## Proof of $c_k = f[x_0,x_1,\ldots,x_k]$

The polynomial interpolating the data at $\{ x_0,x_1,\ldots,x_k \}$ is

$$
p_k(x) = c_0 + c_1 (x-x_0) + c_2 (x-x_0)(x-x_1) + \ldots c_k (x-x_0)\ldots(x-x_{k-1})
$$

Note that $c_k$ depends on $x_0,x_1,\ldots,x_k$ so we indicate this
dependence symbolically as 

$$
c_k = f[x_0,x_1,\ldots,x_k]
$$ 

We will show
that this *symbol satisfies the recursion relation for divided
differences*. Note that $c_k$ is the coefficient of $x^k$ in the
polynomial $p_k(x)$.

Let us also write the Lagrange polynomial interpolating the same data as
above. Define 

$$
\psi_k(x) = (x-x_0)(x-x_1) \ldots (x-x_k)
$$ 

so that

$$
\psi_k'(x_i) = (x_i-x_0) \ldots (x_i - x_{i-1})(x_i-x_{i+1}) \ldots (x_i-x_k)
$$

Then the Lagrange form can be written as

$$
p_k(x) = \sum_{j=0}^k \frac{\psi_k(x)}{(x-x_j)\psi_k'(x_j)} f(x_j)
$$

The coefficient of $x^k$ in this form is

$$
c_k = f[x_0,x_1,\ldots,x_k] = \sum_{j=0}^k \frac{f(x_j)}{\psi_k'(x_j)}
$$

since the interpolating polynomial is unique.

Now, using above result, we can write a formula for
$f[x_0,x_1,\ldots,x_{k-1}]$ but we have to knock off the factor
$(x_j-x_k)$ from $\psi_k'(x_j)$, so that

$$
f[x_0,x_1,\ldots,x_{k-1}] = \sum_{j=0}^{k-1} \frac{f(x_j)}{\psi_k'(x_j)} (x_j-x_k) =
\sum_{j=0}^k \frac{f(x_j)}{\psi_k'(x_j)} (x_j - x_k)
$$ 

Similarly

$$
f[x_1,x_2,\ldots,x_{k}] = \sum_{j=1}^{k} \frac{f(x_j)}{\psi_k'(x_j)} (x_j-x_0) =
\sum_{j=0}^k \frac{f(x_j)}{\psi_k'(x_j)} (x_j - x_0)
$$ 

Hence

$$
f[x_1,x_2,\ldots,x_{k}] - f[x_0,x_1,\ldots,x_{k-1}] 
&= \sum_{j=0}^k \frac{f(x_j)} {\psi_k'(x_j)} (x_k - x_0) \\
&= (x_k - x_0) f[x_0,x_1,\ldots,x_k]
$$ 

which is the recursion relation for divided differences.

:::{prf:remark}
Let $\{i_0, i_1, \ldots, i_k\}$ be an arbitrary permutation of
$\{0,1,\ldots,k\}$. Then

$$
f[x_0,x_1,\ldots,x_k] = \sum_{j=0}^k \frac{f(x_j)}{\psi_k'(x_j)} = \sum_{j=0}^k
\frac{f(x_{i_j})}{\psi_k'(x_{i_j})} = f[x_{i_0}, x_{i_1}, \ldots, x_{i_k}]
$$

Hence the value of $f[x_0,x_1,\ldots,x_N]$ is independent of the order
of the arguments.
:::

## Interpolation error

The Newton form of interpolation to the data $(x_i,f_i)$,
$i=0,1,\ldots,N$ is 

$$
\begin{aligned}
p_N(x) = f_0 & + (x-x_0)f[x_0,x_1] \\
& + (x-x_0)(x-x_1) f[x_0,x_1, x_2] \\
& + \ldots \\
& + (x-x_0) \ldots (x-x_{N-1}) f[x_0,\ldots,x_N]
\end{aligned}
$$ 

Let $t \in \re$ be distinct from the nodes $x_0,x_1,\ldots,x_N$ but otherwise arbitrary. The polynomial interpolating $f(x)$ at $x_0,x_1,\ldots,x_N,t$ is

$$
p_{N+1}(x) = p_N(x) + (x-x_0)\ldots (x-x_N) f[x_0,\ldots,x_N,t]
$$ 

But $p_{N+1}(t) = f(t)$ so that

$$
f(t) = p_N(t) + (t-x_0)\ldots (t-x_N) f[x_0,\ldots,x_N,t]
$$ 

Since $t$ was arbitrary, we could call it $x$. Then the error is

$$
f(x) - p_N(x)  = (x-x_0)\ldots (x-x_N) f[x_0,\ldots,x_N,x]
$$ 

Comparing this with the formula we have derived earlier, we conclude that

$$
f[x_0,\ldots,x_N,x] = \frac{1}{(N+1)!} f^{(N+1)}(\xi), \qquad \xi \in I(x_0,\ldots,x_N,x)
$$

## Another proof of independence of order

We will show that $f[x_0,x_1,\ldots,x_N]$ is independent of the order of the arguments. Let $p_{N-1}(x)$ interpolate the data $\{ x_0, x_1, \ldots, x_{N-1}\}$. In Newton form 

$$
\begin{aligned}
p_{N-1}(x) = f_0 & + (x-x_0)f[x_0,x_1] \\
& + (x-x_0)(x-x_1) f[x_0,x_1, x_2] \\
& + \ldots \\
& + (x-x_0) \ldots (x-x_{N-2}) f[x_0,\ldots,x_{N-1}]
\end{aligned}
$$ 

Consider a random permutation

$$
\{ x_{i_0}, x_{i_1}, \ldots, x_{i_{N-1}} \}
$$ 

and let $Q_{N-1}(x)$ interpolate this data 

$$
\begin{aligned}
Q_{N-1}(x) = f_{i_0} & + (x-x_{i_0})f[x_{i_0},x_{i_1}] \\
& + (x-x_{i_0})(x-x_{i_1}) f[x_{i_0},x_{i_1}, x_{i_2}] \\
& + \ldots \\
& + (x-x_{i_0}) \ldots (x-x_{i_{N-2}}) f[x_{i_0},\ldots,x_{i_{N-1}}]
\end{aligned}
$$ 

But these two polynomials are same $Q_{N-1} = P_{N-1}$ since they interpolate same data, and in particular the error at $x=x_N$ is same 

$$
(x_N-x_0) \ldots(x_N - x_{N-1}) & f[x_0,x_1,\ldots,x_N] = \\
& (x_N-x_{i_0}) \ldots(x_N - x_{i_{N-1}}) f[x_{i_0},x_{i_1},\ldots,x_N]
$$

Since the polynomial factors are same, we get

$$
f[x_0,x_1,\ldots,x_N] = f[x_{i_0},x_{i_1},\ldots,x_N]
$$ 

We can repeat this proof by leaving out $x_{N-1}, x_{N-2}, \ldots, x_0$ and then we obtain the complete proof.

## Forward difference operator

For $h > 0$, define the forward difference

$$
\Delta f(x) = f(x+h) - f(x)
$$ 

We will use uniformly spaced data at

$$
x_i = x_0 + ih, \qquad i=0,1,2,\ldots
$$ 

Then

$$
\Delta f(x_i) = f(x_i+h) - f(x_i) = f(x_{i+1}) - f(x_i) = f_{i+1} - f_i
$$

or, in compact form 

$$
\Delta f_i = f_{i+1} - f_i
$$ 

For $r \ge 0$, define

$$
\Delta^{r+1} f(x) = \Delta^r f(x+h) - \Delta^r f(x)
$$

$$
\Delta^0 f(x) = f(x)
$$ 

$\Delta^{r+1} f(x)$ is called the $r$-th order
forward difference of $f$ at $x$.

:::{prf:lemma}
For $k \ge 0$ 

$$
f[x_0,x_1,\ldots,x_k] = \frac{1}{k! h^k} \Delta^k f_0
$$

:::

:::{prf:proof}
For $k=0$, the result is trivially true since

$$
\frac{1}{0! h^0} \Delta^0 f_0 = f_0 = f[x_0]
$$ 

For $k=1$,

$$
f[x_0,x_1] = \frac{f_1 - f_0}{x_1 - x_0} = \frac{1}{h} \Delta f_0
$$

Assume that the result is true for all forward differences of order
$k \le r$. Then for $k=r+1$ 

$$
\begin{aligned}
f[x_0,x_1,\ldots,x_{r+1}] 
&= \frac{f[x_1,\ldots,x_{r+1}] - f[x_0,\ldots,x_r]}{x_{r+1} - x_0} \\
&= \frac{1}{(r+1)h} \left[ \frac{1}{r! h^r}\Delta^r f_1 - \frac{1}{r! h^r}\Delta^r f_0 \right] \\
&= \frac{1}{(r+1)! h^{r+1}} \Delta^{r+1}f_0
\end{aligned}
$$

:::

## Newton form on uniform data

For uniformly spaced data, we can avoid the use of divided differences
and write in terms of forward differences, which is more efficient since
we avoid division, which is an expensive operation. Define

$$
\mu = \frac{x - x_0}{h}
$$ 

Newton interpolation formula has factors of
the form $(x-x_0)\ldots(x-x_k)$. Since

$$
x-x_j = (x_0 + \mu h) - (x_0 + jh) = (\mu-j) h
$$ 

hence

$$
(x-x_0)\ldots(x-x_k) = \mu(\mu-1)\ldots(\mu-k) h^{k+1}
$$ 

Using
previous lemma, we write divided difference in terms of forward
differences

$$
p_n(x) = f_0 + \mu h \frac{\Delta f_0}{h} + \mu(\mu-1)h^2 \frac{\Delta^2 f_0}{2!
h^2} + \ldots + \mu(\mu-1)\ldots(\mu-n+1)h^n \frac{\Delta^n f_0}{n! h^n}
$$

Define the binomial coefficients

$$
\binom{\mu}{k} = \frac{\mu(\mu-1) \ldots (\mu-k+1)}{k!}, \qquad k > 0
$$

$$
\binom{\mu}{0} = 1
$$ 

Then

$$
p_n(x) = \sum_{j=0}^n \binom{\mu}{j} \Delta^j f_0, \qquad \mu = \frac{x-x_0}{h}
$$

This is the Newton forward difference form of the interpolating polynomial. A schematic of the computation is shown in [](#tab:fd).

:::{table} Forward differences
:label: tab:fd

|  $x_i$ | $f_i$ | $\Delta f_i$ | $\Delta^2 f_i$ | $\Delta^3 f_i$ | $\Delta^4 f_i$ | $\Delta^5 f_i$ |
|  :---: | :---: |  :---:       |  :---:         |    :---:       |    :---:       |    :---:       |
|        |       |              |                |                |                |                |
|  $x_0$ | $f_0$ |              |                |                |                |                |
|        |       | $\Delta f_0$ |                |                |                |                |
|  $x_1$ | $f_1$ |              | $\Delta^2 f_0$ |                |                |                |
|        |       | $\Delta f_1$ |                | $\Delta^3 f_0$ |                |                |
|  $x_2$ | $f_2$ |              | $\Delta^2 f_1$ |                | $\Delta^4 f_0$ |                |
|        |       | $\Delta f_2$ |                | $\Delta^3 f_1$ |                | $\Delta^5 f_0$ |
|  $x_3$ | $f_3$ |              | $\Delta^2 f_2$ |                | $\Delta^4 f_1$ |                |
|        |       | $\Delta f_3$ |                | $\Delta^3 f_2$ |                |                |
|  $x_4$ | $f_4$ |              | $\Delta^2 f_3$ |                |                |                |
|        |       | $\Delta f_4$ |                |                |                |                |
|  $x_5$ | $f_5$ |              |                |                |                |                |

:::

:::{prf:example}
For $n=1$, i.e., given $(x_0,f_0)$, $(x_1,f_1)$

$$
p_1(x) = f_0 + \mu \Delta f_0
$$ 

For $n=2$, i.e., given $(x_0,f_0)$,
$(x_1,f_1)$, $(x_2,f_2)$

$$
p_2(x) = p_0 + \mu \Delta f_0 + \frac{1}{2}\mu(\mu-1)\Delta^2 f_0
$$
:::

::::{prf:example}
Take the function $f(x) = \sqrt{x}$. Given data on $f(x)$ rounded to seven digits. Hence the data is accurate only to seven digits.  Approximate $f(2.15) = \sqrt{2.15} = 1.4662878$

The data and forward differences are shown in [](#tab:fd2). The sixth decimal place may not be correct. The last column has no accuracy since the non-zero digit is at sixth decimal place.

:::{table}  Forward differences for $f(x) = \sqrt{x}$. Data rounded to 7 digits. The 6'th decimal place could be erroneous.
:label: tab:fd2
:align: center

|  $x_i$ |   $f_i$    | $\Delta f_i$ | $\Delta^2 f_i$ | $\Delta^3 f_i$ | $\Delta^4 f_i$ |
|  :---: |   :---:    |   :---:      |   :---:        |   :---:        |   :---:        |
|  $2.0$ | $1.414214$ |              |                |                |                |
|        |            |  $0.034924$  |                |                |                |
|  $2.1$ | $1.449138$ |              |  $-0.000822$   |                |                |
|        |            |  $0.034102$  |                |   $0.000055$   |                |
|  $2.2$ | $1.483240$ |              |  $-0.000767$   |                |  $-0.000005$   |
|        |            |  $0.033335$  |                |   $0.000050$   |                |
|  $2.3$ | $1.516575$ |              |  $-0.000717$   |                |                |
|        |            |  $0.032618$  |                |                |                |
|  $2.4$ | $1.549193$ |              |                |                |                |

:::

::::

## Errors in data and forward differences

We can use a forward difference table to detect noise in physical data, as long as the noise is larger relative to the experimental error.

:::{prf:lemma}

$$
\Delta^r f(x_i) = h^r f^{(r)}(\xi), \qquad x_i \le \xi_i \le x_{i+r}
$$

:::

:::{prf:proof}
We can relate forward difference to divided differences

$$
\begin{aligned}
\Delta^r f_i 
&= h^r r! f[x_i, \ldots,x_{i+r}] \\
&= h^r r! \frac{f^{(r)}(\xi_i)}{r!} \\
&= h^r f^{(r)}(\xi_i)
\end{aligned}
$$
:::

:::{prf:lemma}
For any two functions $f$ and $g$, and constants $\alpha$ and $\beta$

$$
\Delta^r[\alpha f(x) + \beta g(x)] = \alpha \Delta^r f(x) + \beta \Delta^r g(x), \qquad r
\ge 0
$$

:::

:::{prf:proof}
The result is true if $r=0$ or $r=1$. Assume the result is true for all $r \le n$ and prove it for $r=n+1$. 

$$
\begin{aligned}
& \Delta^{n+1}[\alpha f(x) + \beta g(x)]  \\
=& \Delta^n[\alpha f(x+h) + \beta g(x+h)] - \Delta^n[\alpha f(x) + \beta g(x)] \\
=& [\alpha \Delta^n f(x+h) + \beta \Delta^n g(x+h)] - [\alpha \Delta^n f(x) + \beta \Delta^n g(x)] \\
=& \alpha [\Delta^n f(x+h) - \Delta^n f(x)] + \beta [\Delta^n g(x+h) - \Delta^n g(x)] \\ 
=& \alpha \Delta^{n+1} f(x) + \beta \Delta^{n+1} g(x)
\end{aligned}
$$

:::

::::{prf:example}
If the derivatives of $f(x)$ are uniformly bounded, or if $h^n f^{(n)} \to 0$ as $n \to \infty$, then the forward differences $\Delta^n f(x)$ should become smaller as $n$ increases. If the given data $\tilde{f}_i$ has errors so that 

$$
f(x_i) = \tilde{f}_i + e(x_i)
$$ 

then

$$
\Delta^r \tilde{f}_i = \Delta^r f(x_i) - \Delta^r e(x_i)
$$ 

The first term on the right decreases with $r$ but the error term can behave badly. To understand this behaviour, consider a simple case where only one data point has error 

$$
e(x_i) = \begin{cases}
0 & i \ne k \\
\epsilon & i = k
\end{cases}
$$

:::{table} Forward differences
:label: tab:fderr

|    $x_i$   |  $e(x_i)$  | $\Delta e(x_i)$ | $\Delta^2 e(x_i)$ | $\Delta^3 e(x_i)$ | $\Delta^4 f_i$ | $\Delta^5 f_i$ |
|   :---:    |  :---:     |  :---:          |  :---:            |       :---:       |       :---:    |      :---:     |
|   $\vdots$ |            |                 |                   |                   |                |                |
|  $x_{k-2}$ |    $0$     |                 |        $0$        |                   |  $\epsilon$    |                |
|            |            |       $0$       |                   |    $\epsilon$     |                |  $-5\epsilon$  |
|  $x_{k-1}$ |    $0$     |                 |    $\epsilon$     |                   | $-4\epsilon$   |                |
|            |            |   $\epsilon$    |                   |   $-3\epsilon$    |                |  $10\epsilon$  |
|    $x_k$   | $\epsilon$ |                 |   $-2\epsilon$    |                   | $6\epsilon$    |                |
|            |            |   $-\epsilon$   |                   |    $3\epsilon$    |                | $-10\epsilon$  |
|  $x_{k+1}$ |    $0$     |                 |    $\epsilon$     |                   | $-4\epsilon$   |                |
|            |            |       $0$       |                   |    $-\epsilon$    |                |  $5\epsilon$   |
|  $x_{k+2}$ |    $0$     |                 |        $0$        |                   |  $\epsilon$    |                |
|            |            |       $0$       |                   |        $0$        |                |  $-\epsilon$   |
|    $x_5$   |    $0$     |                 |        $0$        |                   |     $0$        |                |

:::


The differences $\Delta^r e(x_i)$ increase with $r$ as shown in [](#tab:fderr).  Hence the differences $\Delta^r \tilde{f}_i$ will also increase with $r$. This is an indication of error in the data. If the $\Delta^r \tilde{f}_i$ start increasing, then the Newton polynomial must be truncated.
::::


:::{prf:theorem}
Let $x_0,x_1,\ldots,x_n$ be distinct and let $f$ be $n$ times continuously differentiable in the interval of these points. Then

$$
f[x_0,x_1,\ldots,x_n] = \int_{\tau_n} f^{(n)}(t_0 x_0 + t_1 x_1 + \ldots + t_n x_n) \ud
t_1 \ldots \ud t_n
$$ 

where

$$
\tau_n = \{(t_1,\ldots,t_n) : t_i \ge 0, \ \sum_{j=1}^n t_i \le 0 \}, \qquad t_0 = 1 -
\sum_{j=1}^n t_j
$$ 

(Note that $0 \le t_0 \le 1$ and $\sum_{j=0}^n t_j = 1$.)
:::

:::{prf:proof}
See [@Atkinson2004], Theorem 3.3.
:::

## Some properties of divided differences

1.  Recall that the divided difference can be written as
    $$
    f[x_0,\ldots,x_n] = \frac{f^{(n)}(\xi)}{n!}
    $$ 
    Prove this result by induction. The Hermite-Gennochi formula relates the 
    sames quantities in terms of an integral.

2.  If $f^{(n)}$ is continuous in the interval $[\min_j x_j, \max_j x_j]$,
    then $f[x_0,\ldots,x_n]$ is a continuous function of
    $x_0, \ldots,x_n$.

3.  From the Hermite-Gennochi formula, we see that the right hand side
    makes sense even if all the $x_j$ are not distinct. If we take all
    arguments to be equal to $x$, then 
    $$
    \begin{aligned}
    f[\underbrace{x,x,\ldots,x}_{n+1 \textrm{ arguments}}] 
    =& \int_{\tau_n} f^{(n)}(x) \ud t_1 \ldots \ud t_n \\
    =& f^n(x) \textrm{ volume}(\tau_n) \\
    =& \frac{f^{(n)}(x)}{n!}
    \end{aligned}
    $$

4.  Consider $f[x_0,\ldots,x_n,x]$ as a function of $x$ and
    differentiate this. By definition of the derivative
    $$
    \begin{aligned}
    \dd{}{x}f[x_0,\ldots,x_n,x] 
    =& \lim_{h \to 0} \frac{f[x_0,\ldots,x_n,x+h] - f[x_0,\ldots,x_n,x]}{h} \\
    =& \lim_{h \to 0} \frac{f[x_0,\ldots,x_n,x+h] - f[x,x_0,\ldots,x_n]}{h} \\
    =& \lim_{h \to 0} f[x,x_0,\ldots,x_n,x+h] \\
    =& f[x,x_0,\ldots,x_n,x] \\
    =& f[x_0,\ldots,x_n,x,x]
    \end{aligned}$$ The quantity on the right has $n+3$ arguments and it
    is well defined if $f \in
    \cts^{n+2}$.
