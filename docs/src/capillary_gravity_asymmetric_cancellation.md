# Proof: Asymmetric cancellation at steady state for $\alpha > 0$

## Upstream-downstream asymmetry as $t\rightarrow\infty$

The time-dependent response $\eta_{\mathrm{tr}}(x,t)$ is given in eqns. (4.5c) and (4.5d) of the manuscript. We may observe that $\mathbb{I}_4(x,t)$ is a decaying non-singular term which can be integrated numerically. This term tends to zero as $t\to\infty$ in the form $\mathbb{I}_4\approx 1/\sqrt{t}$, thereby implying that it will not contribute to the steady state (see fig.7 of the manuscript). On the other hand, $\mathbb{I}_3(x,t)$ is a singular integral which requires to be solved by the PV method and may contain a term contributing to the steady state as $t\to\infty$.

Let us define,

```math
g(k)=k-\chi(k).
\tag{C.1}
```

and

```math
z(k)=\alpha k^2+1-\rho_r.
\tag{C.2}
```

Substituting $g(k)$ and $z(k)$ into $\mathbb{I}_3$ (eqn. (4.5d) of the manuscript) and writing the cosine function in terms of exponential functions yields,

```math
\begin{aligned}
\mathbb{I}_3(x,t)
={}&
-\frac{1+\rho_r}{2\alpha}
\int_0^\infty dk\;
\frac{
\left(k+\chi(k)\right)
}{
z(k)(k-k_l)(k-k_s)
}
\\
&\qquad\times
\Big\{
\exp\left[i\left(tg(k)-kx\right)\right]
+
\exp\left[-i\left(tg(k)-kx\right)\right]
\Big\}.
\end{aligned}
\tag{C.3}
```

Equation (C.3) is valid for $-\infty<x<\infty$.

We break $\mathbb{I}_3$ into four parts corresponding to the two poles and the two exponential terms in the numerator.

```math
\begin{aligned}
\mathbb{I}_3(x,t)
={}&
\mathbb{I}'_{3A}(x,t)
+\mathbb{I}''_{3A}(x,t)
\\
&+
\mathbb{I}'_{3B}(x,t)
+\mathbb{I}''_{3B}(x,t).
\end{aligned}
\tag{C.4}
```

where,

```math
\begin{aligned}
\mathbb{I}'_{3A}(x,t)
={}&
\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_0^\infty dk\;
\frac{
\left(k+\chi(k)\right)
\left\{\exp\left[itg(k)-ikx\right]\right\}
}{
z(k)
}
\left(\frac{1}{k-k_s}\right).
\end{aligned}
\tag{C.5}
```

```math
\begin{aligned}
\mathbb{I}''_{3A}(x,t)
={}&
\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_0^\infty dk\;
\frac{
\left(k+\chi(k)\right)
\left\{\exp\left[-itg(k)+ikx\right]\right\}
}{
z(k)
}
\left(\frac{1}{k-k_s}\right).
\end{aligned}
\tag{C.6}
```

```math
\begin{aligned}
\mathbb{I}'_{3B}(x,t)
={}&
-\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_0^\infty dk\;
\frac{
\left(k+\chi(k)\right)
\left\{\exp\left[itg(k)-ikx\right]\right\}
}{
z(k)
}
\left(\frac{1}{k-k_l}\right).
\end{aligned}
\tag{C.7}
```

```math
\begin{aligned}
\mathbb{I}''_{3B}(x,t)
={}&
-\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_0^\infty dk\;
\frac{
\left(k+\chi(k)\right)
\left\{\exp\left[-itg(k)+ikx\right]\right\}
}{
z(k)
}
\left(\frac{1}{k-k_l}\right).
\end{aligned}
\tag{C.8}
```

To evaluate these integrals in the large-time limit, we note that the function $g(k)$, namely,

```math
g(k)
=
k-(1+\rho_r)^{-1/2}k^{1/2}[z(k)]^{1/2},
\tag{C.9}
```

satisfies the following relations:

```math
g(k_s)=0=g(k_l),
\tag{C.10}
```

```math
g'(k_s)=\frac{\Delta}{2},
\tag{C.11}
```

```math
g'(k_l)=-\frac{\Delta}{2},
\tag{C.12}
```

where

```math
\Delta
=
\sqrt{1-\frac{4\alpha\beta}{1+\rho_r}}
=
\sqrt{1-\frac{\alpha}{\alpha_{\max}}}
>0,
\qquad
\alpha\leq\alpha_{\max}
=
\frac{1}{4}
\frac{(1+\rho_r)^2}{1-\rho_r}.
\tag{C.13}
```

Now $\mathbb{I}'_{3A}(x,t)$, $\mathbb{I}''_{3A}(x,t)$, $\mathbb{I}'_{3B}(x,t)$, and $\mathbb{I}''_{3B}(x,t)$ may be solved separately and added to yield $\mathbb{I}_3(x,t)$.

### Evaluation of $\mathbb{I}'_{3A}(x,t)$ and $\mathbb{I}''_{3A}(x,t)$ as $t\to\infty$

To evaluate $\mathbb{I}'_{3A}$, we Taylor expand the function $g(k)$ around $k_s$ using a new variable

```math
u=k-k_s.
\tag{C.14}
```

Then,

```math
g(k)
=
\frac{u\Delta}{2}
+
\frac{u^2}{2}g''(k_s)
+\ldots .
\tag{C.15}
```

Also,

```math
z(k)
=
(1+\rho_r)k_s
+
(1+\rho_r)(1-\Delta)u
+
\alpha u^2.
\tag{C.16}
```

and

```math
k+\chi(k)
=
2u+2k_s
-\frac{u\Delta}{2}
-\frac{u^2}{2}g''(k_s).
\tag{C.17}
```

Upon substituting eqns. (C.14)–(C.17) into eqn. (C.5), one obtains

```math
\begin{aligned}
\mathbb{I}'_{3A}(x,t)
={}&
\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_{-k_s}^{\infty}
\frac{
\left[
2u+2k_s-\frac{u\Delta}{2}
-\frac{u^2}{2}g''(k_s)
\right]
}{
(1+\rho_r)k_s
+(1+\rho_r)(1-\Delta)u
+\alpha u^2
}
\\
&\qquad\times
\exp\left\{
i\left[
\frac{tu\Delta}{2}
+\frac{tu^2}{2}g''(k_s)
-k_sx-ux
\right]
\right\}
\frac{du}{u}.
\end{aligned}
\tag{C.18}
```

With the change of variable to $w$, where

```math
w=ut,
\tag{C.19}
```

one obtains

```math
\begin{aligned}
\mathbb{I}'_{3A}
={}&
\frac{1+\rho_r}
{2\alpha(k_l-k_s)}
\int_{-k_st}^{\infty}
\frac{
\left[
\frac{2w}{t}
+2k_s
-\frac{w\Delta}{2t}
-\frac{w^2}{2t^2}g''(k_s)
\right]
}{
(1+\rho_r)k_s
+(1+\rho_r)(1-\Delta)\frac{w}{t}
+\alpha\frac{w^2}{t^2}
}
\\
&\qquad\times
\exp\left\{
i\left[
\frac{w\Delta}{2}
+\frac{w^2}{2t}g''(k_s)
-k_sx
-\frac{w}{t}x
\right]
\right\}
\frac{dw}{w}.
\end{aligned}
\tag{C.20}
```

In the limit $t\to\infty$, the lower limit of integration of the above equation, namely $-k_st$, may be extended to $-\infty$. Furthermore, one may neglect terms containing $1/t$ and higher powers of $1/t$ both in the prefactors and in the exponential. This yields

```math
\begin{aligned}
\lim_{t\rightarrow\infty}
\left[\mathbb{I}'_{3A}(x,t)\right]
={}&
\frac{
\exp\left\{-ik_sx\right\}
}{
\alpha(k_l-k_s)
}
\int_{-\infty}^{\infty}
\exp\left\{
\frac{iw\Delta}{2}
\right\}
\frac{dw}{w},
\\
&\qquad -\infty<x<\infty.
\end{aligned}
\tag{C.21}
```

Equation (C.21) is still a singular integral due to the pole at $w=0$. Following the standard procedure, the PV value of the integral may be evaluated as

```math
\operatorname{PV}\!
\left[
\int_{-\infty}^{\infty}
\exp\left\{
\pm iw\frac{\Delta}{2}
\right\}
\frac{dw}{w}
\right]
=
\pm i\pi,
\qquad
\text{if}\quad \Delta>0.
\tag{C.22}
```

With this one obtains

```math
\lim_{t\rightarrow\infty}
\left[\mathbb{I}'_{3A}(x,t)\right]
=
\frac{
\exp\left\{-ik_sx\right\}
}{
\alpha(k_l-k_s)
}
\left(i\pi\right).
\tag{C.23}
```

Proceeding exactly in the same manner for $\mathbb{I}''_{3A}$, one obtains

```math
\begin{aligned}
\lim_{t\rightarrow\infty}
\left[\mathbb{I}''_{3A}(x,t)\right]
={}&
\frac{
\exp\left\{ik_sx\right\}
}{
\alpha(k_l-k_s)
}
\int_{-\infty}^{\infty}
\exp\left\{
\frac{-iw\Delta}{2}
\right\}
\frac{dw}{w}
\\
={}&
\frac{
\exp\left\{ik_sx\right\}
}{
\alpha(k_l-k_s)
}
\left(-i\pi\right).
\end{aligned}
\tag{C.24}
```

Upon combining eqns. (C.23) and (C.24), we obtain the steady-state contribution in the asymptotic limit $t\to\infty$,

```math
\begin{aligned}
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}'_{3A}(x,t)
+
\mathbb{I}''_{3A}(x,t)
\right]
=
\frac{2\pi}
{\alpha(k_l-k_s)}
\sin(k_sx).
\end{aligned}
\tag{C.25}
```

One may perform similar evaluations of the integrals $\mathbb{I}'_{3B}$ and $\mathbb{I}''_{3B}$. For this case, define a variable $u=k-k_l$ and note that the Taylor expansion of $g(k)$ is

```math
g(k)
=
-\frac{u\Delta}{2}
+
\frac{u^2}{2}g''(k_l)
+\ldots .
\tag{C.26}
```

It must be remarked that since $\Delta$ is greater than zero, the coefficient of $u$ is negative, in contrast to the previous case (i.e., at $k=k_s$). This is the fundamental difference which causes the asymmetric solution to the problem. Upon proceeding as earlier, one obtains

```math
z(k)
=
(1+\rho_r)k_l
+
(1+\rho_r)(1+\Delta)u
+
\alpha u^2.
\tag{C.27}
```

```math
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}'_{3B}(x,t)
\right]
=
-\frac{
\exp\left\{-ik_lx\right\}
}{
\alpha(k_l-k_s)
}
\left(-i\pi\right).
\tag{C.28}
```

and

```math
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}''_{3B}(x,t)
\right]
=
-\frac{
\exp\left\{ik_lx\right\}
}{
\alpha(k_l-k_s)
}
\left(i\pi\right).
\tag{C.29}
```

Hence, the solution around the $k_l$ pole is

```math
\begin{aligned}
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}'_{3B}(x,t)
+
\mathbb{I}''_{3B}(x,t)
\right]
=
\frac{2\pi}
{\alpha(k_l-k_s)}
\sin(k_lx).
\end{aligned}
\tag{C.30}
```

Substituting eqns. (C.25) and (C.30) into eqn. (C.4) gives

```math
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}_3(x,t)
\right]
=
\frac{2\pi}
{\alpha(k_l-k_s)}
\left[
\sin(k_sx)+\sin(k_lx)
\right].
\tag{C.31}
```

The limit $t\rightarrow\infty$ of $\eta_{\mathrm{tr}}(x,t)$, given by eqn. (4.5c), is

```math
\begin{aligned}
\lim_{t\rightarrow\infty}
\left[
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
\right]
={}&
-\frac{1}{2\pi}
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}_3(x,t)
\right]
\\
&-
\frac{1}{2\pi}
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}_4(x,t)
\right].
\end{aligned}
\tag{C.32}
```

The integral $\mathbb{I}_4(x,t)$ is a decaying integral (Shown in fig.7 of the manuscript), as

```math
\lim_{t\rightarrow\infty}
\left[
\mathbb{I}_4(x,t)
\right]
=0.
\tag{C.33}
```

Substituting eqn. (C.31) into eqn. (C.32) yields

```math
\lim_{t\rightarrow\infty}
\left[
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
\right]
=
\frac{1}
{\alpha(k_l-k_s)}
\left[
-\sin(k_sx)-\sin(k_lx)
\right].
\tag{C.34}
```

Upon adding eqn. (C.34) to the time-independent solution, eqn. (4.6) of the manuscript, the final expression for the steady-state interfacial displacement $\eta(x)$ is

```math
\frac{\eta(x)}{F_0}
=
\frac{\eta_s(x)}{F_0}
+
\lim_{t\rightarrow\infty}
\left[
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
\right].
\tag{C.35}
```

i.e.,

```math
\begin{aligned}
\frac{\eta(x)}{F_0}
={}&
\frac{1}{\alpha(k_l-k_s)}
\Big\{
-\sin(k_s|x|)
-\sin(k_sx)
\\
&\qquad\qquad
+\sin(k_l|x|)
-\sin(k_lx)
\Big\}
\\
&+
\frac{k_l+k_s}{\pi\alpha}
\int_0^\infty dy\,
\frac{
y\exp\left(-|x|y\right)
}{
\left(y^2+k_l^2\right)
\left(y^2+k_s^2\right)
}.
\end{aligned}
\tag{C.36}
```

We see that the symmetry that existed in $\eta_s(x)$ is now broken by the contribution from the time-dependent term $\eta_{\mathrm{tr}}(x,t)$ in the asymptotic limit $t\to\infty$. Explicitly, the solutions for $x>0$ and $x<0$ may be written as

```math
\begin{aligned}
\frac{\eta(x)}{F_0}
={}&
-\frac{2}
{\alpha(k_l-k_s)}
\sin(k_sx)
\\
&+
\frac{k_l+k_s}{\pi\alpha}
\int_0^\infty dy\,
\frac{
y\exp\left(-xy\right)
}{
\left(y^2+k_l^2\right)
\left(y^2+k_s^2\right)
},
\qquad x>0.
\end{aligned}
\tag{C.37}
```

and

```math
\begin{aligned}
\frac{\eta(x)}{F_0}
={}&
-\frac{2}
{\alpha(k_l-k_s)}
\sin(k_lx)
\\
&+
\frac{k_l+k_s}{\pi\alpha}
\int_0^\infty dy\,
\frac{
y\exp\left(xy\right)
}{
\left(y^2+k_l^2\right)
\left(y^2+k_s^2\right)
},
\qquad x<0.
\end{aligned}
\tag{C.38}
```

This is consistent with the observation of a gravity wave (long wavelength) for $x>0$ and a capillary wave (short wavelength) for $x<0$.

