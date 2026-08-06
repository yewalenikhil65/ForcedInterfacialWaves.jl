# Proof: Evaluation of the time-dependent term $\eta_{tr}(x,t)$ for pure gravity ($\alpha = 0$)

We now focus on the time-dependent contribution $\eta_{tr}(x,t)$, which can be rewritten as

**(T.1)**
```math
\frac{\eta_{tr}(x,t)}{F_0}	= -\dfrac{1}{2\pi(1-\rho_r)}\left\{\mathbb{I}_7(x,t)+\mathbb{I}_8(x,t)\right\},
```

where, we define,

**(T.2)**
```math
\mathbb{I}_7(x,t)	= \int_{0}^{\infty}dk\;\dfrac{k\;\cos\left\{\left(k(t-x) - t\sqrt{k \beta}\right)\right\}}{k- \sqrt{k \beta}}
```

and

**(T.3)**
```math
\mathbb{I}_8(x,t)	= \int_{0}^{\infty}dk\;\dfrac{k\;\cos\left\{\left(k(t-x) + t\sqrt{k \beta}\right)\right\}}{k+ \sqrt{k \beta}}.
```

## Evaluation of the non-singular integral $\mathbb{I}_8(x,t)$

The integral $\mathbb{I}_8(x,t)$ is non-singular and may be simplified to reduce the oscillatory behaviour of the integrand. Introducing the transformation $k=v^2$, eqn. (T.3) becomes,

**(T.4)**
```math
\mathbb{I}_8(x,t)	= 2 \int_{0}^{\infty}dv\;\dfrac{v^2\; \cos\left\{\left(v^2 (t-x) + vt\sqrt{\beta}\right)\right\}}{v+ \sqrt{\beta}}.
```

Using the identity

**(T.5)**
```math
\frac{v^2}{v+\sqrt{\beta}}=v-\sqrt{\beta}+\frac{\beta}{v+\sqrt{\beta}},
```

we write

**(T.6)**
```math
\mathbb{I}_8(x,t)	= \mathbb{I}_8'+\mathbb{I}_8'',
```

where,

**(T.7)**
```math
\mathbb{I}_8'(x,t)= 2 \int_{0}^{\infty}dv\; \left(v-\sqrt{\beta}\right) \cos\left\{\left(v^2 (t-x) + vt\sqrt{\beta}\right)\right\}
```

and

**(T.8)**
```math
\mathbb{I}_8''(x,t)	= 2 \int_{0}^{\infty}dv\;\dfrac{\beta\; \cos\left\{\left(v^2 (t-x) + vt\sqrt{\beta}\right)\right\}}{v+ \sqrt{\beta}}.
```

The integral $\mathbb{I}_8''(x,t)$ is numerically integrable. We now simplify $\mathbb{I}_8'(x,t)$ to remove oscillations. Let $a=t-x$ and $b=\frac{t \sqrt{\beta}}{2}$; then eqn. (T.7) can be rewritten as

**(T.9)**
```math
\mathbb{I}_8'(x,t)= 2 \int_{0}^{\infty}dv\; \left(v+\frac{b}{a}- \left(\sqrt{\beta}+\frac{b}{a}\right)\right) \cos\left\{a v^2  + 2bv\right\},
```

after rearrangement,

**(T.10)**
```math
\mathbb{I}_8'(x,t)= 2 \int_{0}^{\infty}dv\; \left(v+\frac{b}{a}\right) \cos\left\{a v^2  + 2bv\right\}- 2 \int_{0}^{\infty}dv\; \left(\sqrt{\beta}+\frac{b}{a}\right) \cos\left\{a v^2  + 2bv\right\}.
```

Let $z=av^2+2bv$; substituting this into the first integrand on the R.H.S. of the above equation reduces it to an oscillatory cosine integral (i.e. $\frac{1}{a} \int_{0}^{\infty}dz\;  \cos\left\{z\right\}$) whose value is zero, yielding

**(T.11)**
```math
\mathbb{I}_8'(x,t)= - 2 \left(\sqrt{\beta}+\frac{b}{a}\right) \int_{0}^{\infty}dv\;  \cos\left\{a v^2  + 2bv\right\},
```

substituting $a=t-x$ and $b=\frac{t \sqrt{\beta}}{2}$ into the above equation gives,

**(T.12)**
```math
\mathbb{I}_8'(x,t)= - 2 \sqrt{\beta} \left(1+\frac{t}{2(t-x)}\right) \int_{0}^{\infty}dv\;  \cos\left\{a v^2  + 2bv\right\}.
```

Hence,

**(T.13)**
```math
\mathbb{I}_8(x,t)= - 2 \sqrt{\beta} \left(1+\frac{t}{2(t-x)}\right) \int_{0}^{\infty}dv\;  \cos\left\{a v^2  + 2bv\right\}+2\beta \int_{0}^{\infty}dv\;\dfrac{ \cos\left\{\left(v^2 (t-x) + vt\sqrt{\beta}\right)\right\}}{v+ \sqrt{\beta}},
```

the integral in the first term of the R.H.S of the above equation is still an oscillatory integral. We may reduce this oscillatory behaviour by writing it in Fresnel form, i.e.,

**(T.14)**
```math
\int_{0}^{\infty}dv\;  \cos\left\{a v^2  + 2bv\right\}	=\sqrt{\frac{\pi}{2a}}\left[\cos\left(\frac{b^2}{a}\right)\left\{\frac{1}{2}-\operatorname{C}\left(b\sqrt{\frac{2}{\pi a}}\right)\right\}+\sin\left(\frac{b^2}{a}\right)\left\{\frac{1}{2}-\operatorname{S}\left(b\sqrt{\frac{2}{\pi a}}\right)\right\}\right],
```

where,

**(T.15)**
```math
\operatorname{C}\left(b\sqrt{\frac{2}{\pi a}}\right) = \int_{0}^{b\sqrt{\frac{2}{\pi a}}}dt\;  \cos\left\{\frac{\pi t^2}{2}\right\},
```

and

**(T.16)**
```math
\operatorname{S}\left(b\sqrt{\frac{2}{\pi a}}\right) = \int_{0}^{b\sqrt{\frac{2}{\pi a}}}dt\;  \sin\left\{\frac{\pi t^2}{2}\right\}.
```

The Fresnel function definitions in eqns. (T.15)–(T.16) are taken from Abramowitz & Stegun: *Handbook of Mathematical Functions*. The derivation above is valid for $a>0$. Detailed calculation shows that it continues to remain valid for $a<0$, provided we replace $\int_{0}^{\infty}dv\;  \cos\left\{a v^2  + 2bv\right\}$ with $\int_{0}^{\infty}dv\;  \cos\left\{|a| v^2  - 2bv\right\}$.

## Evaluation of the singular integral $\mathbb{I}_7(x,t)$

Now let us focus on the time-dependent singular integral $\mathbb{I}_7(x,t)$. Eqn. (T.2) can be rewritten as,

**(T.17)**
```math
\mathbb{I}_7(x,t)	=  \mathbb{I}_7'(x,t)+\mathbb{I}_7''(x,t),
```

where, we define,

**(T.18)**
```math
\mathbb{I}_7'(x,t)	= \frac{1}{2} \int_{0}^{\infty}dk\;\dfrac{k\;\exp\left\{i\left(k(t-x) - t\sqrt{k \beta}\right)\right\}}{k- \sqrt{k \beta}},
```

and

**(T.19)**
```math
\mathbb{I}_7''(x,t)	= \frac{1}{2} \int_{0}^{\infty}dk\;\dfrac{k\;\exp\left\{-i\left(k(t-x) - t\sqrt{k \beta}\right)\right\}}{k- \sqrt{k \beta}}.
```

The integrals $\mathbb{I}_7'(x,t)$ and $\mathbb{I}_7''(x,t)$ are singular at $k=\beta$, lying along the line of integration, as shown in the contour of Figure 2 (see [Pure-gravity steady-state proof](pure_gravity_steady_proof.md)).
Hence we interpret them to exist in a Principal Value (PV) sense. They may be evaluated by contour integration using the PV method, with the same contours.

Consider the integral $\left[\mathbb{I}_7'(x,t)\right]^c$ along the closed contour (Figure 2(a)) for $x<t$,

**(T.20)**
```math
\left[\mathbb{I}_7'(x,t)\right]^c = \frac{1}{2} \oint dz \;  \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}}.
```

$\left[\mathbb{I}_7'(x,t)\right]^c$ will be zero as the closed contour does not enclose any singularity. Further, upon breaking this integral along the individual segments of the contour, one may write as,

**(T.21)**
```math
\begin{aligned}
	\left[\mathbb{I}_7'(x,t)\right]^c = 0 &= \frac{1}{2}	\int_{\Gamma_1} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}} 
	+
	\frac{1}{2}\int_{\Gamma_2} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}} \\
	&+\frac{1}{2}\int_{\Gamma_3} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}} 
	+
	\frac{1}{2}\int_{\Gamma_4} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}} \\
	&+\frac{1}{2}\int_{\Gamma_5} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}}.
\end{aligned}
```

The integral on the large quarter circle $\Gamma_4$ tends to zero as $R \to \infty$ for $x<t$, as argued below,

**(T.22)**
```math
\lim_{R \rightarrow \infty} \left[\frac{1}{2}\int_{\Gamma_4} dz\, \dfrac{z\;\exp\left\{i\left(z(t-x) - t\sqrt{\beta z}\right)\right\}}{z- \sqrt{\beta z}}\right]
	=\lim_{R \rightarrow \infty}\left[\frac{1}{2} \oint d \theta \; i R \exp(i \theta)  \dfrac{R \exp(i \theta)\;\exp\left\{i\left(R \exp(i \theta)(t-x) - t\sqrt{\beta R \exp(i \theta)}\right)\right\}}{R \exp(i \theta)- \sqrt{\beta R \exp(i \theta)}}  \right].
```

The value of the integrand above is governed by the factor $\exp\left(-R \sin(\theta)(t-x)\right)$ which tends to zero as $R \to \infty$ for $x<t$ by Jordan's lemma, since $\sin(\theta)$ is always positive in the first quadrant. In view of this, eqn. (T.21) reduces to,

**(T.23)**
```math
\begin{aligned}
	&\lim_{\substack{\epsilon \to 0 \\ R \to \infty}} \left[\frac{1}{2}\int_{0}^{\beta-\epsilon}dk\, 
	\dfrac{k\;\exp\left\{i\left(k(t-x) - t\sqrt{\beta k}\right)\right\}}{k- \sqrt{\beta k}} 
	+
	\frac{1}{2}\int_{\beta+\epsilon}^{R}dk\,
	\dfrac{k\;\exp\left\{i\left(k(t-x) - t\sqrt{\beta k}\right)\right\}}{k- \sqrt{\beta k}} \right] \\ 
	&+
	\lim_{\epsilon \to 0}\left[\frac{1}{2} \int_{\pi}^{0}
	d\theta_s \, \frac{ \left\{\beta\right\} \left\{\exp \left(i \left[\beta+\epsilon \exp \left(i\theta_s\right)\right](t-x)-it\sqrt{\beta \left[\beta+\epsilon \exp \left(i\theta_s\right)\right]}\right)\right\}\, \left\{i\epsilon \exp \left(i\theta_s\right)\right\} }
	{\left\{\beta+\epsilon \exp \left(i\theta_s\right)-\sqrt{\beta \left[\beta+\epsilon \exp \left(i\theta_s\right)\right]}\right\}}\right] \\
	&+
	\frac{1}{2} \int_{\infty}^{0} dy\,\exp \left(i\frac{\pi}{2}\right)
	\frac{ \left\{\exp \left(i\frac{\pi}{2}\right)y\right\} \left\{\exp \left(i \left[\exp \left(i\frac{\pi}{2}\right)y(t-x)-t\sqrt{\beta \exp \left(i\frac{\pi}{2}\right)y}\right] \right)\right\} }{\left\{\exp \left(i\frac{\pi}{2}\right)y-\sqrt{\beta \exp \left(i\frac{\pi}{2}\right)y}\right\}}\, =0.
\end{aligned}
```

Upon completing the limiting process and identifying the first two terms with the PV of $\mathbb{I}_7'(x,t)$ and evaluating the remaining terms, one obtains,

**(T.24)**
```math
\operatorname{PV}\! \left[\mathbb{I}_7'(x,t)\right] =  i \pi \beta\, \exp \left(-i\beta x\right)
	+ \frac{1}{2} \int_{0}^{\infty} dy\, \frac{y \left[\exp \left\{-y(t-x)+t \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right] \left[\exp \left\{-it \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right]}{\sqrt{\frac{y}{2}} \sqrt{\beta} -i \left(y-\sqrt{\frac{y}{2}} \sqrt{\beta}\right)}, 
	\quad x<t.
```

Similarly, for $x>t$, one performs similar steps of contour integration using the contour of Figure 2(b) and obtains,

**(T.25)**
```math
\operatorname{PV}\! \left[\mathbb{I}_7'(x,t)\right] =  -i \pi \beta\, \exp \left(-i\beta x\right)
	+ \frac{1}{2} \int_{0}^{\infty} dy\, \frac{y \left[\exp \left\{y(t-x)-t \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right] \left[\exp \left\{-it \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right]}{\sqrt{\frac{y}{2}} \sqrt{\beta} +i \left(y-\sqrt{\frac{y}{2}} \sqrt{\beta}\right)}, \quad x>t.
```

The integral $\mathbb{I}_7''(x,t)$ in eqn. (T.19) is the converse of $\mathbb{I}_7'(x,t)$ with $k(t-x) - t\sqrt{k \beta}$ replaced with $-\left[k(t-x) - t\sqrt{k \beta}\right]$, accordingly one may write

**(T.26)**
```math
\operatorname{PV}\! \left[\mathbb{I}_7''(x,t)\right] =  -i \pi \beta\, \exp \left(i\beta x\right)
	+ \frac{1}{2} \int_{0}^{\infty} dy\, \frac{y \left[\exp \left\{-y(t-x)+t \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right] \left[\exp \left\{it \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right]}{\sqrt{\frac{y}{2}} \sqrt{\beta} + i \left(y-\sqrt{\frac{y}{2}} \sqrt{\beta}\right)},
	\quad x<t,
```

and

**(T.27)**
```math
\operatorname{PV}\! \left[\mathbb{I}_7''(x,t)\right] =  i \pi \beta\, \exp \left(i\beta x\right)
	+ \frac{1}{2} \int_{0}^{\infty} dy\, \frac{y \left[\exp \left\{y(t-x)-t \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right] \left[\exp \left\{it \sqrt{\frac{y}{2}} \sqrt{\beta}\right\}\right]}{\sqrt{\frac{y}{2}} \sqrt{\beta} - i \left(y-\sqrt{\frac{y}{2}} \sqrt{\beta}\right)},
	\quad x>t.
```

Changing the variable of integration to $v=\sqrt{\frac{y}{2}}$ in $\mathbb{I}_7'(x,t)$ (eqn. T.24) and $\mathbb{I}_7''(x,t)$ (eqn. T.26), then substituting them into eqn. (T.17), gives (for $x<t$)

**(T.28)**
```math
\mathbb{I}_7(x,t) =  2 \pi \beta\, \sin(\beta x) + 8 \int_{0}^{\infty} dv\; v^2 \left\{\exp \left(-2v^2 a+vt \sqrt{\beta}\right)\right\} \left\{\frac{\sqrt{\beta}\, \cos(vt \sqrt{\beta})+(2v-\sqrt{\beta})\, \sin(vt \sqrt{\beta})}{\beta+(2v-\sqrt{\beta})^2}\right\}.
```

Changing the variable of integration to $v=\sqrt{\frac{y}{2}}$ in $\mathbb{I}_7'(x,t)$ (eqn. T.25) and $\mathbb{I}_7''(x,t)$ (eqn. T.27), then substituting them into eqn. (T.17), gives (for $x>t$)

**(T.29)**
```math
\mathbb{I}_7(x,t) =  -2 \pi \beta\, \sin(\beta x) + 8 \int_{0}^{\infty} dv\; v^2 \left\{\exp \left(2v^2 a-vt \sqrt{\beta}\right)\right\} \left\{\frac{\sqrt{\beta}\, \cos(vt \sqrt{\beta})-(2v-\sqrt{\beta})\, \sin(vt \sqrt{\beta})}{\beta+(2v-\sqrt{\beta})^2}\right\}.
```

Substituting $\mathbb{I}_7(x,t)$ from eqn. (T.28) and $\mathbb{I}_8(x,t)$ from eqn. (T.13) into eqn. (T.1) yields $\eta_{tr}(x,t)$ for $x<t$. Similarly, substituting $\mathbb{I}_7(x,t)$ from eqn. (T.29) and $\mathbb{I}_8(x,t)$ from eqn. (T.13) into eqn. (T.1) gives $\eta_{tr}(x,t)$ for $x>t$.

It is noteworthy that the time-dependent integral $\mathbb{I}_7(x,t)$ (the first term on the right-hand side of eqns. T.28 and T.29) contains a time-independent contribution. This time-independent term, which emerges from the evaluation of the time-dependent integral, combines with the steady component $\eta_s(x)$ in such a way that, in the long-time limit, the interface becomes flat at $z=0$ for $x<0$, while gravity waves persist for $x>0$.

By contrast, solving only the corresponding time-independent boundary-value problem yields a symmetric wave pattern about $x=0$ (see eqn. P.13 of the [Pure-gravity steady-state proof](pure_gravity_steady_proof.md)), which is not physically admissible. It is therefore more appropriate to formulate the problem as an Initial-Value Problem (IVP), which naturally selects the physically relevant solution.
