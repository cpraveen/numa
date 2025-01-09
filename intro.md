---
exports:
  - format: pdf
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Introduction

```{include} math.md
```

Numerical Analysis is concerned with the construction and analysis of algorithms to solve a mathematically posed problem on a computer so as to generate a numerical answer.  This is necessary because most mathematical equations do not have closed form solutions that can be easily evaluated. Even in cases where a closed form solution exists, it may be necessary to resort to a computer to obtain a numerical answer. The mathematical problems that we want to solve numerically arise from science and engineering, where a physical system is modeled mathematically in terms of equations. These equations may be algebraic, ordinary or partial differential equations, or a combination of these. Their solutions are functions in one or many variables, e.g., space and/or time. We assume that space and time are a continuum and use real numbers to represent them. However, the real numbers form an infinite set and we cannot hope to represent them in a computer which has finite memory and resources. On a computer, real numbers are approximated by a finite set of **floating point numbers**.

+++

## Example: Approximation of functions

A general function $f : [a,b] \to \re$ has an infinite amount of information which we cannot represent on a computer. For most applications it is enough to obtain an approximation to the function.

Given a function $f \in V$ we want to find a function $f_n \in V_n$ where $V_n$ is finite dimensional, such that

$$
\| f - f_n \| = \min_{g \in V_n} \| f - g \|
$$

The norm is usually maximum norm or $L^2$ norm, leading to **uniform approximation** and **least squares approximation**.

The above kind of approximation requires knowledge of the full function which may not be available. We can sample the function at a set of $n$ discrete points 

$$
a \le x_0 < x_1 < \ldots < x_{n-1} \le b
$$

and construct an approximation $f_n(x)$ to $f(x)$, for example by **interpolation**

$$
f_n(x_i) = f(x_i), \qquad i=0,1,\ldots,n-1
$$

The approximation $f_n$ is something we can represent in a computer because it has finite amount of information, and we can evaluate it at some unknown $x$ and hope that $f_n(x) \approx f(x)$. As we let $n \to \infty$, we hope that $f_n \to f$ in some suitable norm.

Interpolation tries to exactly match the approximation to the given value, which may not be a good idea if the data contains noise. Instead of interpolation, we can perform a least squares fit

$$
\|f - f_n\|_n = \min_{g \in V_n} \| f - g \|_n, \quad \| f - f_n \|_n^2 = \sum_{i=0}
^{n-1}[ f(x_i) - f_n(x_i)]^2
$$

:::{note} Some typical questions to ask

1. How to choose the spaces $V_n$ ? For interpolation this includes the choice of the nodes. Polynomials, rational polynomials and trigonometric functions are the most common choices for the function spaces.
1. What is the error $\| f - f_n \|$ ?
:::

+++

## Example: Differential equations

In many problems, the function that we seek may be unknown and it may be the solution of some differential equation: find $f \in V$ such that

$$
L(f) = 0
$$

where $L$ is a differential operator and $V$ is some infinite dimensional function space.  Such a problem is replaced by a finite dimensional problem: find $f_n \in V_n$ such that

$$
L_n(f_n) = 0
$$

where $V_n$ is a finite dimensional space of functions and $L_n$ is an approximation of the differential operator $L$. The space $V_n$ can be constructed by an interpolation process as in the previous example. We would like to show that $f_n \to f$ as $n \to \infty$.

As an example consider

$$
-u''(x) = f(x), \qquad x \in (0,1)
$$

with boundary conditions

$$
u(0) = u(1) = 0
$$

In the finite difference method, we divide the domain with grid points $x_i = ih$, $i=0,1,\ldots,n-1$, $h = 1/(n-1)$ and replace derivatives with finite differences

$$
u_0 &= 0 \\
- \frac{u_{i-1} - 2u_i + u_{i+1}}{h^2} &= f_i, \qquad 1 \le i \le n-2 \\
u_{n-1} &= 0
$$

This is a system of linear equations which can be put in matrix form

$$
AU = b
$$

where $A$ is a tri-diagonal matrix. The numerical solution of boundary value problems leads to matrix equations.

:::{note} Some questions
1. Does the discrete solution $U$ exists ?
1. How to find it ?
1. What is the error $\|U - u\|$ ?
:::

+++

## Example: Roots, matrix equations

Suppose we want to 

* find the roots of a function $f : \re^n \to \re^n$, i.e., find $r \in \re^n$ such that $f(r) = 0$. 
* find the solution of $Ax=b$ where $A$ is a $n \times n$ matrix and $b$ is an $n$-vector. 

Usually such problems arise as a consequence of trying to numerically solve a differential equation, e.g., for a linear PDE we may end up with

$$
L_n(F) = A F - b = 0
$$

where $F \in \re^n$ is the set of function values at the $n$ distinct nodes. To solve such problems we construct algorithms which are iterative in nature, i.e., given an initial guess $F^0$, the algorithm updates it in a series of steps,

$$
F^{k+1} = H(F^k), \qquad k=0,1,2,\ldots
$$

and we hope that as $k \to \infty$ we approach closer to the solution, i.e., $F^k \to F = A^{-1}b$. Note that the solution is a fixed point of $H$.

For solving matrix equations, we can also employ direct methods like Gaussian elimination which give the answer in a finite number of steps. However such methods may be limited to small problem sizes.

:::{note} Some questions
1. Do the iterations converge ?
1. How fast do they converge ?
:::

## Example: Integration

Given a function $f : [a,b] \to \re$, we want to compute a numerical value of the integral

$$
I(f) = \int_a^b f(x) dx
$$

We will construct an approximation of the form

$$
I_n(f) = \sum_{j=0}^{n-1} w_j f(x_j) \approx I(f)
$$

where $n$ is an integer, $w_j$ are some quadrature weights and $x_j \in [a,b]$ are corresponding nodes. 

:::{note} Some questions

1. What is the error in the approximation $I(f) - I_n(f)$ ?
1. Given $n$, how to choose the weights and nodes to get the best possible approximation ?
:::

+++

## Aims of Numerical Analysis

Numerical Analysis is concerned with constructing algorithms to solve such problems and further we would like to analyze these algorithms in terms of the following questions.

**Accuracy**: How well does the numerical solution approximate the exact solution ?

**Convergence**: As we refine the discretization, does the approximation converge to the true solution ? In case of iterative methods, do the iterations converge ?

**Convergence rate, efficiency**: How fast does it converge ? E.g., we may have

$$
\| f_n - f \| \le C n^{-\alpha}, \qquad \alpha > 0
$$

or better still

$$
\| f_n - f \| \le C \exp(-\alpha n), \qquad \alpha > 0
$$

For iterative scheme, we would like

$$
\| F^{k+1} - F \| \le \alpha \| F^k - F \|, \qquad 0 < \alpha < 1
$$

with $\alpha$ as small as possible, or better still

$$
\| F^{k+1} - F \| \le \alpha \| F^k - F \|^p, \qquad \alpha > 0, \quad p > 1
$$

with $p$ as large as possible. Both of these properties will imply convergence of the iterations.

**Stability**: Is the numerical algorithm stable ? If we change the data by a small amount, we hope that the answer given by the algorithm changes by a small amount.

**Backward stability**: Show that the approximate solution to some problem, is the exact solution of a nearby problem. E.g., if $x^*$ is an approximate solution to $Ax=b$, show that it is the exact solution to $(A+E)x=b$ where $E$ is small.

**Algorithms**: Efficient implementation of numerical methods in a working code is very important. E.g., trigonometic interpolation would be costly if implemented in a naive way, while FFT provides a fast implementation.
