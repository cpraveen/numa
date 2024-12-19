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

# Gauss quadrature (contd.)

```{include} math.md
```

For a positive weight function, we want to compute
$$I(f) = \int_a^b w(x) f(x) \ud x \approx I_n(f) = \sum_{j=1}^n w_j f(x_j)$$
Moreover, we want the degree of precision of $I_n$ to be as large as
possible, i.e., $2n-1$. The main idea is to construct the Hermite
interpolant of $f(x)$ and take $I_n(f)$ to be the exact integral of the
Hermite interpolant.

Let $\{ \phi_n\}$ be an orthogonal family of polynomials on $[a,b]$ with
respect to the above weight function $w$. Denote the zeros of
$\phi_n(x)$ by $$a < x_1 < x_2 < \ldots < x_n < b$$ These zeroes will
come out as the nodes in the integration formula $I_n(f)$. We recall
some notation used previously. We write $\phi_n$ in terms of monomials
$$\phi_n(x) = A_n x^n + B_n x^{n-1} + \ldots$$ and define
$$a_n = \frac{A_{n+1}}{A_n}, \qquad \gamma_n = \ip{\phi_n, \phi_n} = \int_a^b w(x) [\phi_n(x)]^2 \ud x$$

::: theorem
For each $n \ge 1$, there is a numerical integration formula $I_n(f)$ of
degree of precision $2n-1$. Assuming $f \in \cts^{2n}[a,b]$
$$\int_a^b w(x) f(x) \ud x = \sum_{j=1}^n w_j f(x_j) + \frac{\gamma_n}{A_n^2 (2n)!} f^{(2n)}(\eta), \qquad \eta \in [a,b]$$
The nodes $\{x_j\}$ are the zeros of $\phi_n(x)$ and
$$w_j = -\frac{a_n \gamma_n}{\phi_n'(x_j) \phi_{n+1}'(x_j)}$$
:::

::: proof
*Proof.* **(a) Construction of the formula**: We will use Hermite
interpolation at the $n$ zeros of $\phi_n$ to construct the quadrature
rule. The interpolant is given by[^1]
$$H_n(x) = \sum_{j=1}^n f(x_j) h_j(x) + \sum_{j=1}^n f'(x_j) \tilde{h}_j(x)$$
where
$$h_j(x) = [1 - 2\ell_j'(x_j)(x-x_j)] [\ell_j(x)]^2, \qquad \tilde{h}_j(x) = (x-x_j) [\ell_j(x)]^2$$
and $$\ell_j(x) = \prod_{i=1,i\ne j}^n \frac{x-x_i}{x_j - x_i}$$ The
Hermite interpolation error is
$$e_n(x) = f(x) - H_n(x) = [\psi_n(x)]^2 f[x_1,x_1,\ldots,x_n,x_n,x] = \frac{[\psi_n(x)]^2}{(2n)!} f^{(2n)}(\eta), \quad \eta \in [a,b]$$
where $$\psi_n(x) = (x-x_1)\ldots(x-x_n) = \frac{1}{A_n} \phi_n(x)$$ Now
$$\int_a^b w(x) f(x) \ud x = \int_a^b w(x) H_n(x) \ud x + \int_a^n w(x) e_n(x) \ud x = I_n(f) + E_n(f)$$
The degree of precision of $I_n$ is atleast $2n-1$ since $H_n$ is exact
for such polynomials. Moreover $$\begin{aligned}
E_n(x^{2n}) &=& \int_a^b w(x) e_n(x) \ud x \\
&=& \int_a^b w(x) \frac{[\psi_n(x)]^2}{(2n)!} (2n)! \ud x \\
&=& \int_a^b w(x) [\psi_n(x)]^2 \ud x > 0
\end{aligned}$$ Hence the degree of precision of $I_n$ is exactly
$2n-1$.\
**(b) Simpler formula for $I_n(f)$**: The integration formula is
$$I_n(f) = \sum_{j=1}^n f(x_j) \int_a^b w(x) h_j(x) \ud x + \sum_{j=1}^n f'(x_j) \int_a^b w(x) \tilde{h}_j(x) \ud x$$
Since
$$\ell_j(x) = \frac{\psi_n(x)}{(x-x_j) \psi_n'(x_j)} = \frac{\phi_n(x)}{(x-x_j) \phi_n'(x_j)}$$
$$\implies \tilde{h}_j(x) = (x-x_j) \ell_j(x) \cdot \ell_j(x) = \frac{\phi_n(x)}{\phi_n'(x_j)} \cdot \ell_j(x)$$
Hence
$$\int_a^b w(x) \tilde{h}_j(x) \ud x = \frac{1}{\phi_n'(x_j)}\int_a^b w(x) \phi_n(x) \ell_j(x) \ud x = 0, \qquad j=1,2,\ldots,n$$
since degree$(\ell_j) = n-1$ and $\phi_n$ is orthogonal to all
polynomials of degree $< n$. Hence
$$I_n(f) = \sum_{j=1}^n w_j f(x_j), \qquad w_j = \int_a^b w(x) h_j(x) \ud x$$
**(c) Uniqueness of $I_n(f)$**: Suppose we have another numerical
integration formula
$$\int_a^b w(x) f(x) \ud x \approx \sum_{j=1}^n v_j f(z_j)$$ that has
degree of precision $\ge 2n-1$. For any $p \in \poly_{2n-1}$, the
Hermite interpolant using nodes $\{ z_j \}$ is exact and
$$p(x) = \sum_{j=1}^n p(z_j) g_j(x) + \sum_{j=1}^n p'(z_j) \tilde{g}_j(x)$$
where
$$g_j(x) = [1 - 2(x-z_j) \bar{\ell}_j(z_j] [\bar{\ell}_j(x)]^2, \qquad \tilde{g}_j(x) = (x-z_j) [\bar{\ell}_j(x)]^2 = \frac{\bar{\ell}_j(x) \chi_n(x)}{\chi_n'(z_j)}$$
$$\bar{\ell}_j(x) = \prod_{i=1,i\ne j}^n \frac{x-z_i}{z_j - z_i}, \qquad \chi_n(x) = (x-z_1)\ldots(x-z_n)$$
Since the quadrature formula is exact for $p \in \poly_{2n-1}$
$$\sum_{j=1}^n v_j p(z_j) = \int_a^b w(x) p(x) \ud x = \sum_{j=1}^n p(z_j) \int_a^b w(x) g_j(x) \ud x + \sum_{j=1}^n p'(z_j) \int_a^b w(x) \tilde{g}_j(x) \ud x$$
Take $p(x) = \tilde{g}_i(x)$. Note that
$$\tilde{g}_i(z_j) = 0, \qquad \tilde{g}_i'(z_j) = \delta_{ij}$$ Hence
we get $$0 = 0 + \tilde{g}_i'(z_i) \int_a^b w(x) \tilde{g}_i(x) \ud x$$
$$\implies \int_a^b w(x) \tilde{g}_i(x) \ud x = \frac{1}{\chi_n'(z_i)} \int_a^b w(x) \chi_n(x) \bar{\ell}_i(x) \ud x = 0, \qquad i=1,2,\ldots,n$$
Since degree$(\bar{\ell}_i)=n-1$ and they form a basis for
$\poly_{n-1}$, we get
$$\int_a^b w(x) \chi_n(x) p(x) = 0, \qquad \forall p \in \poly_{n-1}$$
Hence $\chi_n$ is orthogonal to
$\poly_{n-1} = \textrm{span}(\phi_0,\phi_1, \ldots,\phi_{n-1})$ and note
that $\chi_n \in \poly_n$. But orthogonal polynomials are unique upto a
multiplicative constant. Hence $\chi_n(x) = c \phi_n(x)$ and in
particular their zeros are same, which imples that $z_j = x_j$.
Moreover, since degree$(h_i)=2n-1$
$$w_i = \int_a^b w(x) h_i(x)\ud x = \sum_{j=1}^n v_j \underbrace{h_i(x_j)}_{\delta_{ij}} = v_i$$
**(d) Weights are positive**: $$\begin{aligned}
w_i &=& \int_a^b w(x) h_i(x) \ud x \\
&=& \int_a^b w(x) [1 - 2 \ell_i'(x_i)(x-x_i)] [\ell_i(x)]^2 \ud x \\
&=& \int_a^b w(x) [\ell_i(x)]^2 \ud x - 2 \ell_i'(x_i) \underbrace{\int_a^b w(x) \tilde{h}_i(x) \ud x}_{=0} \\
&=& \int_a^b w(x) [\ell_i(x)]^2 \ud x > 0
\end{aligned}$$ **(e) Error formula**: Assume that
$f \in \cts^{2n}[a,b]$. Then the integration error is $$\begin{aligned}
E_n(f) &=& \int_a^b w(x) e_n(x) \ud x = \int_a^b w(x) [\psi_n(x)]^2 f[x_1,x_1,\ldots,x_n,x_n,x] \ud x \\
&=& f[x_1,x_1,\ldots,x_n,x_n,\xi] \int_a^b w(x) [\psi_n(x)]^2 \ud x, \quad \xi\in[a,b] \\
&=& \frac{f^{(2n)}(\eta)}{(2n)!} \frac{1}{A_n^2} \int_a^b w(x) [\psi_n(x)]^2 \ud x = \frac{f^{(2n)}(\eta)}{(2n)!} \frac{1}{A_n^2} \gamma_n
\end{aligned}$$ **(f) Formula for the weights**: We have already shown
that the weights are positive. To obtain a simpler formula, take
$f(x) = \ell_i(x)$ and note that degree($\ell_i)=n-1$
$$\int_a^b w(x) \ell_i(x) \ud x = \sum_{j=1}^n w_j \ell_j(x_j) + E_n(\ell_i) = w_i + 0$$
Hence
$$w_i = \int_a^b w(x) \ell_i(x) \ud x = \frac{1}{\phi_n'(x_i)} \int_a^b \frac{w(x) \phi_n(x)}{(x-x_i)} \ud x$$
Recall the Christoffel-Darboux identity
$$\sum_{k=0}^n \frac{\phi_k(x) \phi_k(y)}{\gamma_k} = \frac{\phi_{n+1}(x) \phi_n(y) - \phi_n(x) \phi_{n+1}(y)}{a_n \gamma_n (x-y)}, \qquad x \ne y$$
Take $y = x_i$ and use $\phi_n(y)=\phi_n(x_i)=0$
$$\sum_{k=0}^{n-1} \frac{\phi_k(x) \phi_k(x_i)}{\gamma_k} = -\frac{\phi_n(x) \phi_{n+1}(x_i)}{a_n \gamma_n (x-x_i)}$$
Multiply by $w(x)$ and integrate both side, and note that
$$\int_a^b w(x) \phi_k(x) \ud x = 0, \qquad k \ge 1$$ and since
$\phi_0(x)=$ constant, we get
$$\frac{1}{\gamma_0} \int_a^b w(x) \phi_0(x)\phi_0(x_i)\ud x = 1 = - \frac{\phi_{n+1}(x_i)}{a_n \gamma_n} \int_a^b \frac{w(x) \phi_n(x)}{x-x_i} \ud x$$
which proves the formula for the weights. ◻
:::

## Gauss-Legendre quadrature

For weight function $w(x) = 1$ on $[-1,1]$ the orthogonal polynomials
are the Legendre functions $P_n(x)$. The integral is approximated by
$$I(f) = \int_{-1}^1 f(x) \ud x \approx I_n(f)  = \sum_{j=1}^n w_j f(x_j)$$
The nodes $\{x_j\}$ are the zeros of $P_n(x)$ and the weights are
$$w_i = -\frac{2}{(n+1) P_n'(x_i) P_{n+1}(x_i)}, \qquad i=1,2,\ldots,n$$
and the error is
$$E_n(f) = \frac{2^{2n+1} (n!)^4}{(2n+1) [(2n)!]^2} \frac{f^{(2n)}(\eta)}{(2n)!} = e_n \frac{f^{(2n)}(\eta)}{(2n)!}$$

##### Example

$$I = \int_0^\pi \ee^x \cos(x) \ud x$$ The error decreases much faster
than trapezoidal or Simpson rules. With 8 nodes, the error has reached
almost machine precision.

## Asymptotic behavior of error

Define
$$M_m = \max_{x \in [-1,1]} \frac{|f^{(m)}(x)|}{m!}, \qquad m \ge 0$$
For a large class of infinitely differentiable functions
$$\sum_m M_m \le M < \infty$$ In fact, for many functions
$$\lim_{m \to \infty} M_m = 0$$ e.g., analytic functions like $\ee^x$,
$\cos(x)$, etc. For such functions, the Gauss quadrature has error
$$|E_n(f)| \le e_n M_{2n}$$ The speed of convergence depends on $e_n$.
Using Stirling approximation $$n! \approx \ee^{-n} n^n \sqrt{2\pi n}$$
we get $$e_n \approx \frac{\pi}{4^n}$$ E.g., $e_5 = 0.00293$ while
$\pi/4^5 = 0.00307$. Hence we get the asymptotic error bound
$$|E_n(f)| \le \frac{\pi}{4^n} M_{2n}$$ This is exponential convergence
of the error. Compare this to trapezoidal which has $|E_n(f)| \le c/n^2$
and Simpson which has $|E_n(f)| \le c/n^4$. Hence for nice functions
Gauss quadrature is very accurate.

## Relation to minimax approximation

We can see the optimal accuracy of Gauss quadrature by its relation to
the best approximation.

::: theorem
Assume $[a,b]$ is finite. Then the error of Gauss quadrature
$$|E_n(f)| = |I(f) - I_n(f)| \le 2 \rho_{2n-1}(f) \int_a^b w(x) \ud x$$
with $\rho_{2n-1} = \inf_{p \in \poly_{2n-1}} \norm{f - p}_\infty$ is
the minimax error.
:::

::: proof
*Proof.* We have $E_n(p)=0$ for any $p \in \poly_{2n-1}$ and the error
$E_n$ is a linear operator. Take $p = q_{2n-1}^* \in \poly_{2n-1}$ be
the minimax approximation to $f$ on $[a,b]$. Then $$\begin{aligned}
E_n(f) &=& E_n(f) - E_n(q_{2n-1}^*) \\
&=& E_n(f - q_{2n-1}^*) \\
&=& \int_a^b w(x)[f(x) - q_{2n-1}^*(x)] \ud x - \sum_{j=1}^n w_j[f(x_j) - q_{2n-1}^*(x_j)]
\end{aligned}$$ so that
$$|E_n(f)| \le \norm{f - q_{2n-1}^*}_\infty\left[ \int_a^b w(x) \ud x + \sum_{j=1}^n w_j \right]$$
since $w(x) \ge 0$ and $w_j > 0$. As the quadrature rule is exact for
constant function $$\int_a^b w(x) \ud x = \sum_{j=1}^n w_j$$ which
completes the proof.

##### Example

We can use Gauss quadrature for integrals of the form
$$I = \int_0^\infty g(x) \ud x$$ by write it as
$$I = \int_0^\infty \ee^{-x}[ \ee^x g(x)] \ud x = \int_0^\infty \ee^{-x} f(x) \ud x$$
We can use quadrature formula based on Laguerre polynomials. See
Atkinson, page 308 and 309. ◻
:::

[^1]: See Chapter [\[chap:hermint\]](#chap:hermint){reference-type="ref"
    reference="chap:hermint"}