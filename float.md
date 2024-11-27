---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
math:
  '\re': '\mathbb{R}'
  '\fl': '\textrm{fl}(#1)'
  '\half': '\frac{1}{2}'
  '\uround': '\mathfrak{u}'
  '\op': '\otimes'
  '\ap': '\hat{\op}'
  '\hatx': '\hat{x}'
  '\haty': '\hat{y}'
  '\hats': '\hat{s}'
  '\hatf': '\hat{f}'
  '\hatu': '\hat{u}'
  '\rel': '\textrm{RE}\left(#1\right)'
  '\err': '\textrm{AE}\left(#1\right)'
numbering:
  code: false
  math: false
  headings: true
---

# Computer representation of numbers

```{code-cell} ipython3
#%config InlineBackend.figure_format = 'svg'
from pylab import *
```

We are used to decimal system for counting. In general we can use any base $\beta$ for the number system. On most modern computers, the base $\beta = 2$, i.e., we use the binary number system. Any non-zero number $x$ is written as

$$
x = s \cdot (.a_1 a_2 \ldots a_t)_\beta \cdot \beta^e = \textrm{(sign)} \cdot
\textrm{(fractional part)} \cdot \textrm{(exponent)}
$$

where

$$
\begin{aligned}
s = -1 \textrm{ or } +1 \\
0 \le a_i \le \beta - 1 \\
L \le e \le U \textrm{ is an integer}
\end{aligned}
$$

In base 10, the fractional part has the value

$$
(.a_1 a_2 \ldots a_t)_\beta = \frac{a_1}{\beta} + \frac{a_2}{\beta^2} + \ldots +
\frac{a_t}{\beta^t}
$$

We will assume that $a_1 \ne 0$ which is called a *normalized floating point system* and then

$$
1 \le a _i \le \beta-1, \qquad i=1,2,\ldots,t
$$

$(\beta, t, L,U)$ specifies the arithmetic characteristics of the floating point system.

+++

## Chopping and rounding

Most real numbers cannot be exactly represented on a computer since we have finite precision. Hence they must be represented by a nearby number in the floating point system. This requires us to truncate a real number to fit into the floating point system, which can be done in two ways: chopping and rounding.

Consider a real number

$$
x = s \cdot (.a_1 a_2 \ldots a_t a_{t+1} \ldots )_\beta \cdot \beta^e, \qquad a_1 \ne 0
$$

Let

$$
\fl{x} = \textrm{approximation of $x$ in the floating point system}
$$

**chopped machine representation**

$$
\fl{x} = s \cdot (.a_1 a_2 \ldots a_t)_\beta \cdot \beta^e
$$

We throw away all digits in the fractional part that do not fit into the floating-point system.

**rounded machine representation**

$$
\fl{x} = \begin{cases}
s \cdot (.a_1 a_2 \ldots a_t)_\beta \cdot \beta^e & 0 \le a_{t+1} < \half\beta \\
s \cdot [(.a_1 a_2 \ldots a_t)_\beta + (.0\ldots01)_\beta ]\cdot \beta^e & \half\beta \le
a_{t+1} < \beta
\end{cases}
$$

We approximate by the nearest floating point number.

+++

## Error in computer representation

For most real numbers $x \ne \fl{x}$. Define the relative error

$$
\epsilon = \epsilon(x) = \frac{\fl{x} - x}{x}
$$

We have following bounds on the relative error

* Chopping

$$
-\beta^{-t + 1} \le \epsilon \le 0
$$

* Rounding

$$
-\half \beta^{-t + 1} \le \epsilon \le \half \beta^{-t+1}
$$

Thus the floating point approximation $\fl{x}$ of $x \in \re$ satisfies

$$
\fl{x} = x (1 + \epsilon)
$$

This type of analysis was introduced by Wilkinson (1963). It allows us to deal precisely with the effects of rounding/chopping in computer arithmetic operations.

+++

## Proof for chopping

Take sign $s=+1$. Then $\fl{x} \le x$ and

$$
x - \fl{x} = (.00\ldots 0a_{t+1}\ldots)_\beta \cdot \beta^e
$$

Let $\gamma=\beta-1$ (largest digit in base $\beta$). Then

$$
0 & \le x - \fl{x} \\
& \le (.00\ldots 0 \gamma\gamma\ldots)_\beta \cdot \beta^e \\
& = \gamma \left[ \frac{1}{\beta^{t+1}} + \frac{1}{\beta^{t+2}} + \ldots \right]
\beta^e \\
& = \gamma \frac{\beta^{-t-1}}{1- \beta^{-1}} \beta^e \\
& = \beta^{-t + e}
$$

The relative error is

$$
0 &\le \frac{x - \fl{x}}{x} \\
&\le \frac{\beta^{-t+e}}{(.a_1a_2\ldots)_\beta \cdot \beta^e} \\
&\le \frac{\beta^{-t}}{(.100\ldots)_\beta} \\
&= \beta^{-t+1}
$$

+++

## Accuracy of floating numbers: Unit Round

The **unit round** $\uround$

* It is a positive floating point number.
* It is the smallest such number for which
$$
\fl{1 + \uround} > 1
$$

This means that for any $\delta < \uround$, then

$$
\fl{1 + \delta} = 1
$$

1 and $1 + \delta$ are identical within the computer arithmetic. Thus the *unit round characterizes the accuracy of arithmetic computations*. We have

$$
\uround = \begin{cases}
\beta^{-t+1} & \textrm{for chopping} \\
\half\beta^{-t+1} & \textrm{for rounding}
\end{cases}
$$

Hence we get the important relation for floating point numbers

$$
\boxed{\fl{x} = x(1 + \epsilon), \qquad |\epsilon| \le \uround}
$$

+++

**With rounding on a binary computer, show that $\uround = 2^{-t}$**

We must show that

$$
\fl{1 + 2^{-t}} > 1
$$

Firstly $2^{-t} = (.10\ldots 0)_2 \cdot 2^{-t+1}$ is a floating point number.

$$
1 + 2^{-t} &= [ (.10\ldots 0)_2 + (.00\ldots 010\ldots)_2] \cdot 2^1 \\
&= (.10\ldots 01)_2 \cdot 2^1 \qquad \textrm{(1 at  $(t+1)$'th position})
$$

Hence rounding gives

$$
\fl{1 + 2^{-t}} =& (.10\ldots 10)_2 \cdot 2^1 \\
& \textrm{rounding gives 1 at $t$'th position} \\
=& 1 + 2^{-t+1} \\
>& 1
$$

If $\delta < \uround$, then

$$
1 + \delta =& [ (.10\ldots 0)_2 + (.00\ldots 00\ldots)_2] \cdot 2^1 \\
& \qquad\qquad \textrm{zero in $(t+1)$'th position} \\
=& (.10\ldots 00\ldots)_2 \cdot 2^1
$$

Hence, rounding to $t$ places gives

$$
\fl{1+\delta} = (.10\ldots 0)_2 \cdot 2^1 = 1
$$

The following code finds the unit round by finding the smallest $i$ such that $1 + 2^{-i} = 1$.

```{code-cell}
i, x = 0, 1.0
while 1.0 + x > 1.0:
    x = x / 2
    i = i + 1
print(i)
```

This shows that in double precision, $1 + 2^{-53} = 1$, and thus the unit round is $\uround = 2^{-53}$.

+++

## Underflow and overflow

The smallest positive floating point number in the system $(\beta,t,L,U)$ is

$$
x_L = (.10\ldots.0)_\beta \cdot \beta^L = \beta^{L-1}
$$

and the largest positive floating point number is

$$
x_U = (.\gamma\gamma\ldots\gamma)_\beta \cdot \beta^U = (1-\beta^{-t}) \beta^U
$$

The range of numbers that are representable in the floating-point system is

$$
x_L \le |x| \le x_U
$$

Obviously, we cannot represent a number which is outside this range.

* Overflow error: $|x| > x_U$
  * leads to fatal error in computation
  * program may terminate
* Underflow error: $|x| < x_L$
  * may not be a fatal error
  * $\fl{x}$ is set to zero and computations can continue
  * However, in some cases, accuracy may be severly lost.

+++

## IEEE floating point system

The IEEE floating point system was developed under the leadership of Prof. William Kahan of Univ. of California, Berkely.  Project 754 was formed in 1970 to produce the best possible definition of floating-point arithmetic. The IEEE 754 Standard gives

* definition of floating-point numbers
* rounding operations
* format conversion
* exception rules (invalid operation, division by zero, overflow, underflow)

### IEEE single precision

Fortran77: `real*4`, Fortran90: `real`, C/C++: `float`

* 32 bit word
* 1 bit for sign
* 8 bits for biased exponent $(L=-126, U=127)$; $2^8 = 256 = 126 + 1 + 127 + 1 + 1$, the last two are for Inf and NaN.
* 23 bits for fractional part (mantissa)
* Largest number = $(2 - 2^{-23}) \cdot 2^{127} \approx 3.4028 \times 10^{38}$
* Smallest positive normalized number = $2^{-126} \approx 1.1755 \times 10^{-38}$
* Unit roundoff, $\uround = 2^{-24} \approx 5.96 \times 10^{-8}$, corresponds to about seven decimal places.

### IEEE double precision

Fortran77: `real*8`, Fortran90: `double precision`, C/C++: `double`

* Based on two 32 bit words or one 64 bit word
* 1 bit for the sign
* 11 bits for the biased exponent $(L=-1022, U=1023)$
* 52 bits for the fractional part
* Largest number = $(2 - 2^{-52}) \cdot 2^{1023} \approx 1.7977 \times 10^{308}$
* Smallest positive normalized number = $2^{-1022} \approx 2.2251 \times 10^{-308}$
* Unit round $\uround = 2^{-53} \approx 1.11 \times 10^{-16}$

+++

## Precision of computer arithmetic and physics

With IEEE double precision, the interval $[1,2]$ is approximated by about $10^{16}$ floating point numbers. In a handful of solid or liquid, or a balloon of gas, the number of atoms/molecules in a line from one point to another is of the order of $10^8$ (cube root of Avogadro number). Such a system behaves like a continuum, enabling us to define physical quantities like density, pressure, stress, strain and temperature.  Computer arithmetic is more than a million times finer tham this.

The fundamental constants of physics are known to only a few decimal places.

* gravitational constant $G$ - 4 digits
* Planck's constant $h$ - 7 digits
* elementary charge $e$ - 7 digits
* ratio of magnetic moment of electron to the Bohr magneton $\mu_e/\mu_B$ - 12 digits

Nothing in physics is known to more than 12 or 13 digits of accuracy. IEEE numbers are orders of magnitude more precise than any number in science. Mathematical quantities like $\pi$ are of course known to more accuracy.

In the early days of numerical analysis, computers were not very precise, numerical algorithms, especially conditioning and stability were not well understood. Hence there was too much emphasis on studying rounding errors. Today, computer arithmetic is more precise and well understood, and we can control rounding errors through stability of numerical algorithms.

> The main business of numerical analysis is designing algorithms that converge quickly; rounding error analysis is rarely the central issue. If rounding errors vanished, 99\% of numerical analysis would remain. - Lloyd N. Trefethen

+++

## Propagation of errors

Let us consider the basic arithmetic operations: $+, -, *, /$. The computer version of these operations is an approximation due to rounding.

$$
\op: \textrm{exact operation} \qquad \ap: \textrm{computer approximation}
$$

Let

$$
x, y : \textrm{real numbers}, \qquad \hatx, \haty : \textrm{computer representation}
$$

where

$$
\hatx = \fl{x} = x(1+\epsilon_x), \qquad \haty = \fl{y} = y (1 + \epsilon_y)
$$

Let us write this as

$$
x = \hatx + \epsilon, \qquad y = \haty + \eta
$$

Also

$$
x \op y : \textrm{exact value}, \qquad \hatx \ap \haty : \textrm{computed value}
$$

where

$$
\hatx \ap \haty = \fl{\hatx \op \haty}
$$

$\ap$ comes about because $\hatx \op \haty$ may not be representable in the floating point system and has to be rounded !!! The error is

$$
x \op y - \hatx \ap \haty = \underbrace{[x \op y - \hatx \op \hat y]}
_{\textrm{propagated error}} + \underbrace{[\hatx \op \haty - \hatx \ap \haty]}
_{\textrm{rounding error}}
$$

The rounding error is bounded by the unit round

$$
|\hatx \op \haty - \hatx \ap \haty| \le |\hatx \op \haty| \half \beta^{-t+1}
$$

and the only way to control this is to use higher precision. Let us study the propagated error.

+++

### Multiplication

$$
xy - \hatx\haty = xy - (x-\epsilon)(y-\eta) = x\eta + y\epsilon - \epsilon\eta
$$

The relative error is

$$
\rel{\hatx\haty} := \frac{xy - \hatx\haty}{xy} = \frac{\epsilon}{x} + \frac{\eta}{y} -
\frac{\epsilon}{x} \frac{\eta}{y} = \rel{\hatx} + \rel{\haty} - \rel{\hatx}
\rel{\haty}
$$

Hence

$$
\rel{\hatx\haty} \approx \rel{\hatx} + \rel{\haty}
$$

+++

### Division

Consider division of two numbers

$$
\rel{\frac{\hatx}{\haty}} = \frac{\rel{x} - \rel{y}}{1 - \rel{\haty}}
$$

If $\rel{\haty} \ll 1$ then

$$
\rel{\frac{\hatx}{\haty}} \approx \rel{\hatx} - \rel{\haty}
$$

For multiplication and division, relative errors do not propagate rapidly.

### Addition/subtraction

Consider adding or subtracting two numbers

$$
(x \pm y) - (\hatx \pm \haty) = (x - \hatx) \pm (y - \haty) = \epsilon + \eta
$$

The absolute error is

$$
\err{\hatx \pm \haty} = \err{\hatx} \pm \err{\haty}
$$

The absolute error is well-behaved, but relative errors can be large if $x$ and $y$ are nearly equal and we are doing subtraction.

+++

:::{prf:example} Round-off errors

Evaluate

$$
f(x) = 1 - \cos x
$$

for $x = 10^{-8}$. This yields $f = 0.0$ even in double precision.

```{code-cell} ipython3
from math import cos
x = 1.0e-8
y = 1.0 - cos(x)
print("%24.14e" % y)
```

An equivalent expression is

$$
f(x) = 2 \sin^2(x/2)
$$

which yields $f = 5 \times 10^{-17}$. 

```{code-cell} ipython3
from math import sin
z = 2.0*sin(0.5*x)**2
print("%24.14e" % z)
```

The first form suffers from roundoff errors while the second one is more accurate. This may or may not be a problem. A situation where this is problematic is if you want to compute

$$
f(x) = \frac{1 - \cos(x)}{x^2}
$$

```{code-cell} ipython3
y = (1.0 - cos(x))/x**2
print("%24.14e" % y)
```

We will get an answer of zero, but the correct value is closer to $\half$.

```{code-cell} ipython3
z = 2.0*sin(0.5*x)**2 / x**2
print("%24.14e" % z)
```

:::

+++

:::{prf:example} Round-off error

Consider computing

$$
x = 10^{-8}, \qquad y = \sqrt{1+x^2} - 1
$$

```{code-cell} ipython3
from math import sqrt
x = 1.0e-8
y = sqrt(1.0 + x**2) - 1.0
print("%20.10e" % y)
```

The result is zero due to round-off error. An equivalent expression is

$$
y = \frac{x^2}{\sqrt{1 + x^2} + 1}
$$

```{code-cell} ipython3
y = x**2/(sqrt(1.0+x**2) + 1.0)
print("%20.10e" % y)
```

This is more accurate

$$
y = \sqrt{1 + x^2} - 1 \approx (1 + x^2/2)-1 = \frac{1}{2}x^2 = 5.0 \times 10^{-17}
$$

:::

+++

:::{prf:example} polynomial evaluation

Let us evaluate

$$
y = x^3 - 3 x^2 + 3 x - 1
$$

in single precision. Note that it can also be written as 

$$
y = (x-1)^3
$$

```{code-cell}
x = linspace(0.99,1.01,50,dtype=float32)
y = x**3 - 3.0*x**2 + 3.0*x - 1.0
plot(x,y,'-o',x,(x-1)**3)
xlabel('x'), ylabel('y')
legend(('Single precision','Exact'));
```

Now, let us do it in double precision.

```{code-cell} ipython3
x = linspace(0.99,1.01,50,dtype=float64)
y = x**3 - 3.0*x**2 + 3.0*x - 1.0
plot(x,y,'-o',x,(x-1)**3)
xlabel('x'), ylabel('y')
legend(('Double precision','Exact'));
```

:::

+++

## Errors in summation

Suppose we want to compute the sum

$$
s_n = \sum_{i=1}^n x_i
$$

In a computer program, we would perform a loop and accumulate the sum sequentially.  E.g., in fortran, the code may look like this:

```fortran
s = 0.0
do i=1,n
   s = s + x(i)
enddo
```

Let us assume that all the $x_i$ are floating point numbers. Let $s_1 = x_1$ and compute

$$
s_{k+1} = s_k + x_{k+1}, \qquad k=1,2,\ldots,n-1
$$

But a computer will perform roundoff in each step, so we actually compute

$$
\hats_1 = x_1, \qquad \hats_{k+1} = \fl{\hats_k + x_{k+1}} = (\hats_k + x_{k+1})(1 + \epsilon_{k+1})
$$

where each $|\epsilon_k| \le \uround$. Define the relative error

$$
\rho_k = \frac{\hats_k - s_k}{s_k}
$$

Then

$$
\rho_{k+1} =& \frac{\hats_{k+1} - s_{k+1}}{s_{k+1}} \\
=& \frac{(\hats_k + x_{k+1})(1+\epsilon_{k+1}) - (s_k + x_{k+1})}{s_{k+1}} \\
=&  \frac{[s_k(1+\rho_k) + x_{k+1}](1+\epsilon_{k+1}) - (s_k + x_{k+1})}
{s_{k+1}} \\
=& \epsilon_{k+1} + \rho_k \frac{s_k}{s_{k+1}}(1 + \epsilon_{k+1})
$$

Let us assume that each $x_i \ge 0$ so that $s_k \le s_{k+1}$. Hence

$$
|\rho_{k+1}| \le \uround + |\rho_k| (1 + \uround) = \uround + |\rho_k| \theta, \qquad \theta := 1 + \uround
$$

and so

$$
|\rho_1| =& 0 \\
|\rho_2| \le& \uround \\
|\rho_3| \le& \uround + \uround\theta = (1 + \theta) \uround\\
|\rho_4| \le& \uround + (\uround + \uround\theta)\theta = (1 + \theta + \theta^2) \uround \\
\vdots&
$$

The relative error in the sum can be bounded as

$$
|\rho_n| \le& (1 + \theta + \ldots + \theta^{n-2}) \uround \\
=& \left[ \frac{\theta^{n-1} - 1}{\theta-1} \right] \uround \\
=& (1 + \uround)^{n-1} - 1 \\
\approx& (n-1) \uround
$$

We have $(n-1)$ additions and the above estimate says that the relative error is propertional to the number of additions. Usually this is a pessimistic estimate since we took absolute values, and in actual computations, the succesive errors can cancel, leading to much smaller total error in the sum.

### A closer look

Let us estimate some of the partial sums.

:::{math}
\hats_1 =& x_1 \\
\hats_2 =& (\hats_1 + x_2)(1 + \epsilon_2) = s_2 + s_2 \epsilon_2 \\
\hats_3 =& (\hats_2 + x_3)(1 + \epsilon_3) = (s_3 + s_2 \epsilon_2)(1+\epsilon_3) =
s_3 + s_2 \epsilon_2 + s_3 \epsilon_3 + O(\uround^2) \\
\hats_4 =& (\hats_3 + x_4)(1+\epsilon_4) = (s_4 + s_2 \epsilon_2 + s_3 \epsilon_3 +
O(\uround^2))(1 + \epsilon_4) \\
=& s_4 + s_2 \epsilon_2 + s_3 \epsilon_3 + s_4 \epsilon_4 + O(\uround^2) \\
\vdots &
:::

Hence the error in the sum is

:::{math}
\hats_n - s_n =& s_2 \epsilon_2 + \ldots + s_n \epsilon_n + O(\uround^2) \\
=& x_1 (\epsilon_2 + \ldots + \epsilon_n) + x_2 (\epsilon_2 + \ldots + \epsilon_n) + x_3 (\epsilon_3 + \ldots + \epsilon_n) + \ldots + x_n \epsilon_n
:::

Note that the coefficients of $x_i$ decrease with increasing $i$. If we want to minimize the error, we can sort the numbers in increasing order $x_1 \le x_2 \le \ldots \le x_n$ and then sum them. In practice, the $\epsilon_i$ can be positive or negative, so the precise behaviour can be case dependent.

+++

:::{exercise} Error in product

Consider computing the product of $n$ floating point numbers

$$
p_n = x_1 x_2 \ldots x_n
$$

Here is a fortran code

```fortran
p = x(1)
do i=2,n
   p = p * x(i)
enddo
```

Show that the rounding error is bounded by

$$
\frac{|\hat{p}_n - p_n|}{|p_n|} \le (n-1) \uround + O(\uround^2)
$$

:::

+++

## Function evaluation

Suppose we want to compute a function value

$$
y = f(x), \qquad x \in \re
$$

There are two types of errors we have to deal with on a computer. Firstly, the real number has to be approximated by a floating point number

$$
\hatx = \fl{x}
$$

Secondly, the function $f$ may be computed only approximately as $\hatf$. The only computations a computer can do exactly are addition and multiplication, and even here we have roundoff errors. Any other type of computation, even division, can only be performed approximately. E.g., if $f(x) = \sin(x)$, a computer will evaluate it approximately by some truncated Taylor approximation of $\sin(x)$. So what we actually evaluate on a computer is

$$
\haty = \hatf(\hatx)
$$

The absolute error is

$$
\err{\haty} =& |y - \haty| \\
=& |f(x) - \hatf(x) + \hatf(x) - \hatf(\hatx)| \\
\le& |f(x) - \hatf(x)| + |\hatf(x) - \hatf(\hatx)|
$$

We have two sources of error.

1. $|f(x) - \hatf(x)|$: numerical discretization error, can be controlled by using an **accurate** algorithm
1. $|\hatf(x) - \hatf(\hatx)|$: transmission of round-off error by the numerical scheme, can be controlled by **stability** of numerical scheme and using enough precision so that $x - \hatx$ is small.

+++

:::{exercise}

The function

$$
f(x) = coth(x) - 1/x
$$

is finite at $x=0$ but the two terms are not. Note that

$$
f(x) = \frac{x}{3} - \frac{x^3}{45} + \frac{2 x^5}{945} + \ldots
$$

Approximate this in two parts $[0,\delta]$ and $[\delta,\infty)$. How to choose $\delta$ ?

:::
