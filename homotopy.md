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

# Homotopy and continuation methods

```{include} math.md
```

```{code-cell}
#%config InlineBackend.figure_format = 'svg'
from pylab import *
from scipy.integrate import odeint
from scipy.linalg import solve, norm
```

The material of this Chapter is from [@Kincaid2002].

## Main idea

Let $f : X \to Y$ and consider the problem: find $x \in X$ such that $f(x) = 0$. This problem may be difficult to solve; for Newton method, we need a good initial guess, which we may not know.

Suppose we can easily solve the problem: find $x \in X$ such that $g(x) = 0$. Define 

$$
h(t,x) = t f(x) + (1-t) g(x)
$$ 

Note that

$$
h(0,x) = g(x), \qquad h(1,x) = f(x)
$$ 

Partition the interval $[0,1]$ into 

$$
0 = t_0 < t_1 < t_2 < \ldots < t_m = 1
$$ 

We can easily solve the problem $h(t_0,x) = g(x) = 0$, let $x_0$ be the solution. Now to solve the problem $h(t_1,x)=0$, say using Newton method, we can use $x_0$ as the initial guess. Since the two problems are close, $x_0$ can be expected to be a good initial guess for the next problem. This process can be continued until we solve the problem $h(t_m,x) = f(x) = 0$ using the initial guess $x_{m-1}$.

:::{prf:example}
To find a solution of 

$$
f(x) = x + \log x - 2 = 0
$$ 

consider the related problem 

$$
h(t,x) = x + t \log x - 2 = 0
$$ 

For $t=0$ the above equation has the solution $x=2$. Solve the problem for increasing values of $t$ using previous solution as initial guess.
:::

:::{prf:definition}
Let $f : X \to Y$ and $g : X \to Y$. Then $f$ and $g$ are said to be homotopic if there exists a continuous map $h : [0,1] \times X \to Y$ such that 

$$
h(0,x) = g(x), \qquad h(1,x) = f(x)
$$ 

$h$ is called a homotopy that connects the two functions $f$ and $g$.
:::

The simplest example is the linear interpolation between the two problems, which we discussed above.

:::{prf:example}
Choose an $x_0 \in X$ define

$$
h(t,x) = t f(x) + (1-t)[\underbrace{f(x) - f(x_0)}_{g(x)}] = f(x) + (t-1) f(x_0)
$$

which is a homotopy connecting $g$ and $f$. The problem $h(0,x) = f(x)-f(x_0) = 0$ has the solution $x = x_0$.
:::

## An ODE

Instead of solving the sequence of problems, we will convert it into an ODE problem. The solution of $h(t,x)=0$ depends on the parameter $t$, so $x(t)$ satisfies 

$$
h(t, x(t)) = 0
$$ 

Differentiate wrt $t$

$$
h_t (t, x(t)) + h_x(t,x(t)) x'(t) = 0
$$ 

so that

$$
x'(t) = - [h_x(t,x(t))]^{-1} h_t(t,x(t)) =: R(t, x(t))
$$ 

For this ODE to be valid, we need $h$ to be differentiable wrt $t$ and $x$ and $h_x$ must be non-singular. Suppose we know that $x_0$ solves $h(0,x)=0$. We solve the ODE approximately upto $t=1$ with initial condition $x(0) = x_0$.

:::{prf:example}
Let $X = Y = \re^2$ and $x = (\xi_1, \xi_2) \in \re^2$. Consider the problem

$$
f(x) = 
\begin{bmatrix}
\xi_1^2 - 3 \xi_2^2 + 3 \\
\xi_1 \xi_2 + 6 
\end{bmatrix} = 0
$$ 

Let $x_0 = (1,1)$ and define the homotopy 

$$
h(t,x) = t f(x) + (1-t)(f(x) - f(x_0))
$$ 

Then

$$
h_x = f'(x) = \left[ \df{f_i}{\xi_j} \right]_{ij} = \begin{bmatrix}
2\xi_1 & - 6 \xi_2 \\
\xi_2 & \xi_1 \end{bmatrix}, \qquad h_t = f(x_0) = \begin{bmatrix}
1 \\
7 \end{bmatrix}
$$ 

$$
h_x^{-1} = \frac{1}{\Delta} \begin{bmatrix}
\xi_1 & 6 \xi_2 \\
-\xi_2 & 2\xi_1 \end{bmatrix}, \qquad \Delta = 2\xi_1^2 + 6 \xi_2^2
$$

The ODE is given by 

$$
\dd{}{t}\begin{bmatrix}
\xi_1 \\ \xi_2 \end{bmatrix} = -h_x^{-1} h_t =
-\frac{1}{\Delta} \begin{bmatrix}
\xi_1 + 42 \xi_2 \\
-\xi_2 + 14 \xi_1 \end{bmatrix}
$$ 

Integrating the ODE from $t=0$ to $t=1$ using some numerical ODE solver, we get

$$
x(1) \approx \begin{bmatrix}
-2.990 \\
2.006 \end{bmatrix}
$$ 

while the exact root is $\begin{bmatrix} -3 \\ 2 \end{bmatrix}$. We do not get the exact root because the ODE solver is approximate. Starting with the approximate solution of the ODE, perform a few steps of Newton-Raphson iterations to get a more accurate estimate of the root.
:::

:::{prf:theorem}
If $f : \re^n \to \re^n$ is continuously differentiable and if $\norm{[f'(x)]^{-1}} \le M$ on $\re^n$, then for any $x_0 \in \re^n$ there is a curve $\{ x(t) \in \re^n : 0 \le t \le 1 \}$ such that

$$
h(t,x(t)) = f(x(t)) + (t-1) f(x_0) = 0, \qquad 0 \le t \le 1
$$ 

The function $t \to x(t)$ is continuously differentiable solution of the IVP

$$
x'(t) = -[f'(x(t))]^{-1} f(x_0), \qquad x(0) = x_0
$$

:::

:::{prf:example} continue previous example
We now implement code to apply this idea.  This is the function whose roots are required and its gradient.

```{code-cell}
def f(x):
    y = zeros(2)
    y[0] = x[0]**2 - 3*x[1]**2 + 3
    y[1] = x[0]*x[1] + 6
    return y

def df(x):
    y = array([[2*x[0], -6*x[1]],
               [x[1],    x[0]]])
    return y
```

This is the rhs of the ODE obtained from homotopy method.

```{code-cell}
def F(x,t):
    y = zeros(2)
    delta = 2*x[0]**2 + 6*x[1]**2
    y[0] = -(x[0] + 42*x[1])/delta
    y[1] = -(-x[1] + 14*x[0])/delta
    return y
```

We solve the ode with a relaxed error tolerance using [scipy.odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html).

```{code-cell}
x0 = array([1.0,1.0])
t = linspace(0,1,100)
x = odeint(F,x0,t,rtol=0.1)

# plot results
plot(t,x), xlabel('t'), ylabel('x(t)'), grid(True)
legend(('$\\xi_1$', '$\\xi_2$'))

# Final solution
xf = array([x[-1,0],x[-1,1]])
print('x =',xf)
print('f(x) =',f(xf))
```

The graph shows how $(\xi_1,\xi_2)$ change with $t$; at the final solution of ODE, $f(x)$ is not close to zero. Now we can improve the final solution by applying Newton-Raphson method.

```{code-cell}
N = 10
eps = 1.0e-13

print('%18.10e %18.10e %18.10e' % (xf[0], xf[1], norm(f(xf))))
for i in range(N):
    J = df(xf)
    dx = solve(J,-f(xf))
    xf = xf + dx
    print('%18.10e %18.10e %18.10e' % (xf[0], xf[1], norm(f(xf))))
    if norm(dx) < eps*norm(xf):
        break

print('f(x) =',f(xf))
```

We converge to many digits in about three iterations.
:::

+++

## Relation to Newton method

Choose an $x_0$ and consider the homotopy 

$$
h(t,x) = f(x) - \ee^{-t} f(x_0)
$$ 

so that

$$
h(0,x) = f(x) - f(x_0), \qquad h(\infty,x) = f(x)
$$ 

The curve $x(t)$ such that 

$$
h(t,x(t)) = 0 = f(x(t)) - \ee^{-t} f(x_0)
$$ 

satisfies the ODE 

$$
0 = f'(x(t)) x'(t) + \ee^{-t} f(x_0) = f'(x(t)) x'(t) + f(x(t))
$$

so that 

$$
x'(t) = -[f'(x(t))]^{-1} f(x(t))
$$ 

Applying forward Euler scheme with $\Delta t = 1$ we get

$$
x_{n+1} = x_n - [f'(x_n)]^{-1} f(x_n)
$$ 

we get the Newton-Raphson method.

## Linear programming

Let $c,x \in \re^n$, $A \in \re^{m \times n}$, $m < n$ and consider the constrained optimization problem 

$$
\begin{cases}
\max\limits_{x \in \re^n} c^\top x \\
\textrm{subject to } A x = b, \qquad x \ge 0
\end{cases}
$$ 

Choose a starting point $x^0$ which satisfies the constraints, i.e.

$$
x^0 \in \mathcal{F} = \{ x \in \re^n : A x = b, \quad x \ge 0 \}
$$ 

Our goal is to find $x^1$ such that $c^\top x^1 > c^\top x^0$ and $x^1 \in
\mathcal{F}$. Equivalently, find a curve $x(t)$ such that

1.  $x(0) = x^0$

2.  $x(t) \ge 0$ $\forall t \ge 0$

3.  $A x(t) = b$

4.  $c^\top x(t)$ is increasing function of $t$, $t \ge 0$

Let us try to find an equation satisfied by such a curve $x(t)$

$$
x'(t) = f(x), \qquad x(0) = x^0
$$ 

To satisfy (2), choose $f$ such that $f_i(x) \to 0$ as $x_i \to 0$. One such choice is

$$
D(x) = \diag{[x_1, x_2, \ldots, x_n]}
$$

$$
f(x) = D(x) G(x), \qquad x_i' = x_i G_i(x)
$$ 

To satisfy (3) we must have $A x'(t) = 0$ or $A D(x) G(x) = 0$. Choose 

$$
G(x) = P(x) H(x)
$$

where $P : \re^n \to$ null space of $AD(x)$. So we can take

$$
P(x) = \textrm{orthogonal projection of $x$ into the null space of $AD(x)$}
$$

The null space of $AD(x)$ is a subspace of $\re^n$. By Projection Theorem in Hilbert spaces, we have

$$
\forall v \in \re^n, \qquad Pv \in \textrm{null space of } AD(x)
$$

$$\textrm{ iff } \ip{v-Pv,w} = 0 \quad \forall w \in \textrm{null space of } AD(x)
$$

In particular 

$$
\ip{v - Pv, Pv} = 0, \qquad \forall v \in \re^n
$$ 

Such a $P$ is given by 

$$
P = I - (AD)^\top [ (AD)(AD)^\top]^{-1} (AD)
$$ 

The existence of $P$ requires that $(AD)(AD)^\top = (m\times n)(n\times m) = (m \times m)$ is invertible, and hence rank $AD = m$. This is true if each
$x_i > 0$ and rank $A=m$.

Now (4) will be satisfied by the solution of the ODE if

$$
\begin{aligned}
0 \le& \dd{}{t} c^\top x(t) \\
=& c^\top f(x) \\
=& c^\top D P(x) H(x) \\
=& (Dc)^\top P(x) H(x)
\end{aligned}
$$ 

Let us choose $H(x) = D(x) c$, then 

$$
\begin{aligned}
(Dc)^\top P (Dc) 
=& v^\top P v, \qquad v = Dc \\
=& \ip{v, Pv} \\
=& \ip{v-Pv+Pv,Pv} \\
=& \ip{v-Pv,Pv} + \ip{Pv,Pv} \\
=& \ip{Pv,Pv} \\
\ge& 0
\end{aligned}
$$ 

We have constructed an ODE model whose solution has all the required properties. How do we compute $Pv$ ? Using the definition of $P$ requires computing an inverse matrix which can be expensive if the sizes involved are large. Instead, define $B=AD$ and given $v \in \re^n$, compute $Pv$ by the following steps.

1.  Solve for $z$ 
$$
B B^\top z = Bv
$$

2.  $Pv = v - B^\top z$

Now we are in a position to solve the ODE

$$
x'(t) = D(x) P(x) D(x) c, \qquad x(0) = x^0 \in \mathcal{F}
$$ 

The simplest approach is the Euler method which finds $x^{k+1} \approx x(t_k+\Delta t_k)$ from

$$
\frac{x^{k+1} - x^k}{\Delta t_k} = f(x^k), \qquad k=0,1,2,\ldots
$$

$$
x^{k+1} = x^k + \Delta t_k f(x^k)
$$ 

It is easy to verify that $Ax^{k+1}=b$ since $Af(x^k) = 0$. It remains to ensure that $x^{k+1} > 0$ which may not be automatically satisfied by the Euler method. We can satisfy this condition by choosing $\Delta t_k$ small enough, 

$$
x_i^{k+1} = x_i^k + \Delta t_k f_i(x^k) \ge 0
$$ 

e.g.,

$$
\Delta t_k = \frac{9}{10} \min_{i} \left[ \frac{-x_i^k}{f_i(x^k)} \right]
$$

The factor $9/10$ ensures the strict inequality $x^{k+1} > 0$.
