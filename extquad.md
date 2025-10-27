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

# Extrapolation methods

```{include} math.md
```

Extrapolation methods exploit the error estimate to compute better approximations. Such methods can be used for any approximation methods but here we illustrate it for numerical integration.

## Aitken extrapolation

Suppose we have an asymptotic error formula

$$
I - I_n \approx \frac{c}{n^p}, \qquad p > 0
$$ 

If we know $I$ then we can estimate $p$ from 

$$
\frac{I - I_n}{I - I_{2n}} \approx 2^p
$$ 

How can we estimate the convergence rate $p$ when we dont know $I$ ? First compute $I_n$, $I_{2n}$, $I_{4n}$. Then 

$$
\begin{aligned}
\frac{I_{2n} - I_n}{I_{4n} - I_{2n}} 
&= \frac{ (I-I_n) - (I - I_{2n})}{(I - I_{2n}) - (I - I_{4n})} \\
&\approx \frac{ c/n^p - c/(2n)^p}{c/(2n)^p - c/(4n)^p} \\
&= 2^p
\end{aligned}
$$ 

Hence, the convergence rate is

$$
p \approx \frac{ \log\frac{I_{2n} - I_n}{I_{4n} - I_{2n}}}{\log 2}
$$ 

We can obtain a more accurate estimate of $I$ as follows. Again using the error estimate

$$
\frac{I - I_n}{I - I_{2n}} \approx 2^p \approx \frac{I - I_{2n}}{I - I_{4n}}
$$ 

or

$$
(I - I_n)(I - I_{4n}) = (I - I_{2n})^2
$$ 

Solving this

$$
I \approx I_{4n} - \frac{(I_{4n} - I_{2n})^2}{(I_{4n} - I_{2n}) - (I_{2n} - I_n)} =: \tilde I_{4n}
$$

Then $\tilde I_{4n}$ is more accurate than $I_{4n}$.

:::{prf:example}
Approximate 

$$
I = \int_0^1 x \sqrt{x} \ud x
$$ 

using Simpson rule with $n=16,32,64$ intervals (we have $n+1$ points). Then

$$
I - I_{64} \approx -4.29 \times 10^{-7}, \qquad I - \tilde I_{64} \approx 6.13 \times 10^{-10}
$$

Clearly $\tilde I_{64}$ is more accurate than $I_{64}$. Moreover

\begin{align}
\tilde I_{64} - I_{64} 
&\approx 0.399999999387 - 0.400000429413 \\
&= -4.30026 \times 10^{-7}  \\
&\approx I - I_{64}
\end{align}
:::

:::{prf:remark}
1. If we have computed $I_{4n}$ then we can compute $I_{2n}$ and $I_n$ with very little effort since the same function values are used[^1]. 
1. We can then also obtain a more accurate estimate $\tilde I_{4n}$ very easily. 
1. We also have a **computable error estimate**

    $$
    \tilde I_{4n} - I_{4n} \approx I - I_{4n}
    $$ 

    If $|\tilde I_{4n} - I_{4n}| \le \epsilon$, then we can assume that $|I - I_{4n}| \lesssim \epsilon$ and use $\tilde I_{4n}$ as the estimate of the integral.
:::

## Richardson extrapolation

This is a more sophisticated way of performing the extrapolation. Let us illustrate this technique using the trapezoidal method for which we know the error estimate of the form (see [](#sec:trapzerror))

$$
I - I_n = \frac{d_2^{(0)}}{n^2} + \frac{d_4^{(0)}}{n^4} + \ldots + \frac{d_{2m}^{(0)}}{n^{2m}} + F_{n,m}
$$

when $f \in \cts^{2m+2}[a,b]$. For $n$ even

$$
I  - I_{n/2} = \frac{4 d_2^{(0)}}{n^2} + \frac{16 d_4^{(0)}}{n^4} + \frac{64 d_6^{(0)}}{n^6} + \ldots
$$

Combining the last two equations

$$
4(I - I_n) - (I - I_{n/2}) = -\frac{12 d_4^{(0)}}{n^4} - \frac{60 d_6^{(0)}}{n^6} + \ldots
$$

so that

$$
\label{eq:ImIn1}
I = \clr{red}{\frac{1}{3}(4 I_n - I_{n/2})}  -\frac{4 d_4^{(0)}}{n^4} - \frac{20 d_6^{(0)}}{n^6} + \ldots
$$

Define

$$
\boxed{I_n^{(1)} = \frac{4}{3} I_n^{(0)} - \frac{1}{3} I_{n/2}^{(0)}, \qquad n \ge 2 \textrm{ and $n$ is even}}
$$

with 

$$
I_n^{(0)} = I_n
$$ 

$I_n^{(1)}$ is called the Richardson extrapolation of $I_n$.  Note that $I_n^{(1)}$ uses the same data as $I_n$ but has a convergence rate of $4$ while $I_n$ converges at rate $2$. 

But what is $I_n^{(1)}$ ? 

$$
\begin{aligned}
I_n^{(1)} 
&= \frac{4}{3} h [\shalf f_0 + f_1 + f_2 + \ldots + f_{n-1} + \shalf f_n] \\
& \quad - \frac{1}{3} (2h)[\shalf f_0 + f_2 + f_4 + \ldots + f_{n-2} + \shalf f_n] \\
&= \frac{h}{3}[f_0 + 4 f_1 + 2 f_2 + \ldots + 4 f_{n-2} + 2 f_{n-2} + f_n ]
\end{aligned}
$$ 

which is just Simpson rule (integral of quadratic interpolant) !!!

This procedure can be continued further. Let $n$ be a multiple of 4. We
have already shown in [](#eq:ImIn1) that

$$
I - I_n^{(1)} = \frac{d_4^{(1)}}{n^4} + \frac{d_6^{(1)}}{n^6} + \ldots
$$

where

$$
d_4^{(1)} = -4 d_4^{(0)}, \qquad d_6^{(1)} = -20 d_6^{(0)}, \ldots
$$

and so

$$
I - I_{n/2}^{(1)} = \frac{16 d_4^{(1)}}{n^4} + \frac{64 d_6^{(1)}}{n^6} + \ldots
$$

Hence

$$
16(I - I_n^{(1)}) - (I - I_{n/2}^{(1)}) = -\frac{48 d_6^{(1)}}{n^6} + \ldots
$$

and

$$
I = \clr{red}{ \frac{16 I_n^{(1)} - I_{n/2}^{(1)}}{15} } - \frac{48 d_6^{(1)}}{15 n^6} + \ldots
$$

Then

$$
\boxed{I_n^{(2)} = \frac{16 I_n^{(1)} - I_{n/2}^{(1)}}{15}, \qquad n \ge 4, \quad \textrm{$n$ is multiple of 4}}
$$

converges at the rate 6 and is same as Boole's rule (integral of cubic interpolant).

(sec:simprichexterr)=
### Richardson error estimate for $I_n^{(1)}$ (Simpson rule)

From

$$
I = \frac{16 I_n^{(1)} - I_{n/2}^{(1)}}{15} - \frac{48 d_6^{(1)}}{15 n^6} + \ldots
$$

we get 

$$
\begin{aligned}
I - I_n^{(1)} &= \frac{16 I_n^{(1)} - I_{n/2}^{(1)}}{15} - I_n^{(1)} - \frac{48 d_6^{(1)}}{15 n^6} + \ldots \\
&= \frac{1}{15}[ I_n^{(1)} - I_{n/2}^{(1)} ] + \order{\frac{1}{n^6}}
\end{aligned}
$$ 

Since $I_n^{(1)} - I_{n/2}^{(1)}  = \order{1/n^4}$, we
can use the error estimate

$$
I - I_n^{(1)} \approx \frac{1}{15}[ I_n^{(1)} - I_{n/2}^{(1)} ]
$$ 

The
right hand side can be easily computed.

## Romberg integration

Using the error formula of trapezoidal rule, we can generalize Richardson extrapolation to arbitrary orders. For $n = 2^k$, $k \ge 1$ some integer

$$
\label{eq:romiter}
I_n^{(k)} = \frac{4^k I_n^{(k-1)} - I_{n/2}^{(k-1)}}{4^k - 1}
$$ 

and 

$$
I - I_n^{(k)} = \order{\frac{1}{n^{2k+2}}}
$$

Romberg integration is defined by

$$
J_k(f) = I_{2^k}^{(k)}, \qquad k=0,1,2,\ldots
$$ 

For example, if $k=3$, compute 

$$
\begin{array}{cccc}
I_1^{(0)} & & & \\
I_2^{(0)} & I_2^{(1)} & & \\
I_4^{(0)} & I_4^{(1)} & I_4^{(2)} & \\
I_8^{(0)} & I_8^{(1)} & I_8^{(2)} & I_8^{(3)}\\
\end{array}
$$ 

The first column is trapezoid rule; $I_8^{(0)}$ is obtained using all $n = 2^3 = 8$ intervals, $I_4^{(0)}$ is obtained using the half the intervals, etc. 

The second column is generated from the first column using [](#eq:romiter) with $k=1$. The third column is generated from second column and the fourth from the third. Then $I_8^{(3)}=J_3(f)$ is
the final answer.

[^1]: This is true for methods like trapezoid, mid-point, Simpson but is
    not true for Gauss rules.
