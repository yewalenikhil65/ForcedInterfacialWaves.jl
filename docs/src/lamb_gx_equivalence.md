# Proof: Equivalence between the local contribution and Lamb's local function $G(x)$

Lamb's local function $G(x)$ is given by the classical steady-state result. This function may be shown to be exactly the same as the exponential-decay (local) term in the time-independent solution of the present Initial Value Problem (IVP), given in eqn. (S.12) of the [Steady-state proof](steady_proof.md). This is shown here by applying the Cauchy-residue theorem to $G(x)$.

Combine the two cosine integrals in $G(x)$ into a single integral of $\cos(kx)$:

**(L.1)**
```math
G(x) = \int_{0}^{\infty} dk\, \frac{\cos(kx)}{(k+k_s)(k+k_l)}.
```

Writing the cosine as its exponential form and splitting the integral gives

**(L.2)**
```math
G(x) = \frac{1}{2}\left[\mathbb{I}_{11}(x) + \mathbb{I}_{12}(x)\right],
```

where

**(L.3)**
```math
\mathbb{I}_{11}(x) = \int_{0}^{\infty} dk\, \frac{\exp\left(ikx\right)}{(k+k_s)(k+k_l)},
```

**(L.4)**
```math
\mathbb{I}_{12}(x) = \int_{0}^{\infty} dk\, \frac{\exp\left(-ikx\right)}{(k+k_s)(k+k_l)}.
```

The integrals $\mathbb{I}_{11}(x)$ and $\mathbb{I}_{12}(x)$ may be evaluated by contour integration using the Cauchy-residue theorem. The corresponding closed contours (Figure L1(a) for $x>0$ and Figure L1(b) for $x<0$) enclose no poles, since $-k_s$ and $-k_l$ lie outside the first/fourth quadrant quarter-circle contours used here.

![Contours for evaluating the Lamb G(x) integrals](assets/supp5.png)

---

Consider the integral $\mathbb{I}_{11}^c(x)$ along the closed contour (Figure L1(a)) for $x>0$,

**(L.5)**
```math
\mathbb{I}_{11}^{\,c}(x) = \oint dz\, \frac{\exp\left(izx\right)}{(z+k_s)(z+k_l)}.
```

$\mathbb{I}_{11}^c(x)$ will be zero as the closed contour does not enclose any poles. Breaking the integral along the individual contour segments,

**(L.6)**
```math
\mathbb{I}_{11}^c(x) = 0 = \int_{\Gamma_1} dz\, \frac{\exp\left(izx\right)}{(z+k_s)(z+k_l)}
	+\int_{\Gamma_2} dz\,\frac{\exp\left(izx\right)}{(z+k_s)(z+k_l)}
	+\int_{\Gamma_3} dz\, \frac{\exp\left(izx\right)}{(z+k_s)(z+k_l)}.
```

The integral on the large quarter circle $\Gamma_2$ tends to zero as $R \to \infty$ for $x>0$, as argued below,

**(L.7)**
```math
\lim_{R\rightarrow\infty}\left[\mathbb{I}_{11}^{\,c}(x)\right]
	=
	\lim_{R\rightarrow\infty}
	\left[\oint d\theta\,
	iR\exp(i\theta)\,
	\frac{\exp\left(ixR\cos\theta\right)\,\exp\left(-xR\sin\theta\right)}
	{\left(R\exp(i\theta)+k_s\right)\left(R\exp(i\theta)+k_l\right)}\right].
```

The value of the integrand above is governed by the factor $\exp(-xR\sin\theta)$, which tends to zero as $R \to \infty$ for $x>0$ by Jordan's lemma, since $\sin(\theta)$ is always positive in the first quadrant. In view of this, eqn. (L.6) reduces to,

**(L.8)**
```math
\lim_{R\rightarrow\infty}\left[\int_{0}^{R}  dk\,\frac{\exp\left(ikx\right)}{(k+k_s)(k+k_l)}-i\int_{0}^{R} dy\,
	\frac{\exp\left(-yx\right)}{(iy+k_s)(iy+k_l)}\right]=0.
```

Upon completing the limiting process, the above equation can be rewritten to obtain

**(L.9)**
```math
\mathbb{I}_{11}(x)=\int_{0}^{\infty} dk\,\frac{\exp\left(ikx\right)}{(k+k_s)(k+k_l)}
	=
	i\int_{0}^{\infty} dy\,
	\frac{\exp\left(-yx\right)}{(iy+k_s)(iy+k_l)}, \qquad x>0.
```

Similarly, for $x<0$ one performs similar steps of contour integration using the contour in Figure L1(b) and obtains,

**(L.10)**
```math
\mathbb{I}_{11}(x)=\int_{0}^{\infty} dk\,\frac{\exp\left(ikx\right)}{(k+k_s)(k+k_l)}
	=
	-i\int_{0}^{\infty} dy\,
	\frac{\exp\left(yx\right)}{(iy-k_s)(iy-k_l)}, \qquad x<0.
```

The integral $\mathbb{I}_{12}(x)$ in eqn. (L.4) is the same as $\mathbb{I}_{11}(x)$ with $x$ replaced with $-x$, accordingly one may write

**(L.11)**
```math
\mathbb{I}_{12}(x)=\int_{0}^{\infty} dk\,\frac{\exp\left(-ikx\right)}{(k+k_s)(k+k_l)}
	=-i\int_{0}^{\infty} dy\,
	\frac{\exp\left(-yx\right)}{(iy-k_s)(iy-k_l)}
	, \qquad x>0,
```

and

**(L.12)**
```math
\mathbb{I}_{12}(x)=\int_{0}^{\infty} dk\,\frac{\exp\left(-ikx\right)}{(k+k_s)(k+k_l)}
	=
	i\int_{0}^{\infty} dy\,
	\frac{\exp\left(yx\right)}{(iy+k_s)(iy+k_l)}, \qquad x<0.
```

Upon plugging eqns. (L.9)–(L.12) into eqn. (L.2) and separating the real and imaginary parts (the latter vanishes), one obtains

**(L.13)**
```math
G(x)=
	(k_l+k_s)\int_{0}^{\infty} dy\,
	\frac{y\,\exp\left(-y|x|\right)}{(y^2+k_s^2)(y^2+k_l^2)}\,, \quad -\infty<x<\infty.
```

Eqn. (L.13), obtained from Lamb's local function $G(x)$, is exactly the same as the local contribution in the present steady solution — compare with the exponential term in eqn. (S.12) of the [Steady-state proof](steady_proof.md). Furthermore, it has already been shown analytically that the far-field component of the present steady solution is exactly the same as Lamb's result. Hence it is proved analytically that the long-time limit of the present IVP solution is exactly the same as Lamb's steady solution.
