# Large-time asymptotics of the two-fluid transient integrals $I_7(x,t)$ and $I_8(x,t)$

*Pure-gravity limit, $\alpha=0$*

## Transient contribution

For the two-fluid problem in the pure-gravity limit $\alpha=0$, the time-dependent contribution to the interface deformation (eqns. 4.2b and 4.2c of the manuscript) can be combined and rewritten as,

```math
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
        =
        -\frac{1}{2\pi(1-\rho_r)}
        \left[
        I_7(x,t)+I_8(x,t)
        \right],
```

where, we define,

```math
I_7(x,t)
        =
        \int_0^\infty
        \frac{k}{k-\sqrt{\beta k}}
        \cos\left[
        k(t-x)-t\sqrt{\beta k}
        \right]\,dk ,
```

and

```math
I_8(x,t)
        =
        \int_0^\infty
        \frac{k}{k+\sqrt{\beta k}}
        \cos\left[
        k(t-x)+t\sqrt{\beta k}
        \right]\,dk ,
```

```math
\rho_r=\frac{\rho_u}{\rho_l},
        \qquad
        \beta=\frac{1-\rho_r}{1+\rho_r}.
```

The objective is to determine the large-time behaviour of $I_7(x,t)$ and $I_8(x,t)$ for

```math
t\rightarrow\infty,
        \qquad
        x=\text{fixed}.
```

## Large-time asymptotics of the singular integral $I_7(x,t)$

The defining expression for $I_7(x,t)$ may be rewritten as

```math
I_7(x,t)
        =
        \int_0^\infty
        \frac{k}{k-\sqrt{\beta k}}
        \cos\left\{
        t\left(k-\sqrt{\beta k}\right)-kx
        \right\}\,dk .
```

### Transformation of $I_7$

Introduce

```math
k=\beta z^2.
```

Then

```math
\sqrt{\beta k}=\beta z,
        \qquad
        dk=2\beta z\,dz,
```

and

```math
k-\sqrt{\beta k}
        =
        \beta z(z-1).
```

Hence,

```math
\frac{k}{k-\sqrt{\beta k}}\,dk
        =
        \frac{2\beta z^2}{z-1}\,dz .
```

The phase becomes

```math
\begin{aligned}
        t\left(k-\sqrt{\beta k}\right)-kx
        &=
        t\left(\beta z^2-\beta z\right)
        -\beta xz^2\\
        &=
        \beta
        \left[
        t(z^2-z)-xz^2
        \right].

\end{aligned}
```

Therefore,

```math
\boxed{
            I_7(x,t)
            =
            \int_0^\infty
            \frac{2\beta z^2}{z-1}
            \cos\left\{
            \beta
            \left[
            t(z^2-z)-xz^2
            \right]
            \right\}\,dz
        }
```

The transformed integral has a simple pole at

```math
z=z_p=1,
```

corresponding to

```math
k=k_p=\beta.
```

In addition, the rapidly oscillating phase has a stationary point. The pole and stationary point are distinct, so their leading contributions may be evaluated separately.

### Stationary-phase contribution to $I_7$

Define the leading $t$-dependent part of the phase as

```math
\phi(z)
        =
        \beta(z^2-z).
```

The remaining factor $e^{-i\beta xz^2}$ will be retained in the amplitude.

A stationary point satisfies

```math
\phi'(z_s)=0.
```

Since

```math
\phi'(z)
        =
        \beta(2z-1),
```

we obtain

```math
\boxed{
            z_s=\frac{1}{2}.
        }
```

Therefore,

```math
k_s
        =
        \beta z_s^2
        =
        \frac{\beta}{4}.
```

Thus,

```math
k_s=\frac{\beta}{4},
        \qquad
        k_p=\beta,
```

and the stationary point is separated from the pole.

#### Complex representation

Using

```math
\cos\theta
        =
        \Re\left(e^{i\theta}\right).
```

define

```math
J_7
        =
        \int
        a(z)
        e^{\,i\beta[t(z^2-z)-xz^2]}
        \,dz,
```

where

```math
a(z)
        =
        \frac{2\beta z^2}{z-1}.
```

The complex representation of $J_7$ may be written as

```math
J_7
        =
        \int
        a(z)e^{-i\beta xz^2}
        e^{it\phi(z)}
        \,dz.
```

Define

```math
g(z)
        =
        a(z)e^{-i\beta xz^2},
```

so that

```math
J_7
        =
        \int
        g(z)e^{it\phi(z)}
        \,dz .
```

#### Expansion about the stationary point

Taylor expanding $\phi(z)$ about $z=z_s$,

```math
\phi(z)
        =
        \phi(z_s)
        +
        \phi'(z_s)(z-z_s)
        +
        \frac{\phi''(z_s)}{2}(z-z_s)^2
        +\cdots .
```

Since

```math
\phi'(z_s)=0,
```

we have

```math
\phi(z)
        \simeq
        \phi(z_s)
        +
        \frac{\phi''(z_s)}{2}(z-z_s)^2.
```

In the present problem $\phi(z)$ is quadratic, so this expression is exact for $\phi$.

Substituting into this standard form,

```math
J_{7,\mathrm{SP}}
        \simeq
        g(z_s)e^{it\phi(z_s)}
        \int
        \exp\left[
        \frac{it\phi''(z_s)}{2}(z-z_s)^2
        \right]\,dz .
```

Introduce

```math
u
        =
        \sqrt{
            \frac{t\phi''(z_s)}{2}
        }
        (z-z_s).
```

Then

```math
dz
        =
        \sqrt{
            \frac{2}
            {t\phi''(z_s)}
        }\,du.
```

Since $z_s=1/2$ is an interior stationary point, the local limits may be extended asymptotically to $-\infty<u<\infty$. Hence,

```math
J_{7,\mathrm{SP}}
        \simeq
        g(z_s)e^{it\phi(z_s)}
        \sqrt{
            \frac{2}
            {t\phi''(z_s)}
        }
        \int_{-\infty}^{\infty}
        e^{iu^2}\,du .
```

Using

```math
\int_{-\infty}^{\infty}
        e^{iu^2}\,du
        =
        \sqrt{\pi}\,e^{i\pi/4},
```

we obtain

```math
J_{7,\mathrm{SP}}
        \simeq
        g(z_s)
        e^{it\phi(z_s)}
        e^{i\pi/4}
        \sqrt{
            \frac{2\pi}
            {t\phi''(z_s)}
        }.
```

#### Evaluation at $z_s=1/2$

At the stationary point,

```math
a(z_s)
        =
        \frac{2\beta(1/4)}{1/2-1}
        =
        -\beta.
```

Furthermore,

```math
\phi(z_s)
        =
        \beta
        \left(
        \frac14-\frac12
        \right)
        =
        -\frac{\beta}{4},
```

and

```math
\phi''(z_s)
        =
        2\beta.
```

Also,

```math
g(z_s)
        =
        -\beta e^{-i\beta x/4}.
```

Substitution into the stationary-phase formula gives

```math
\begin{aligned}
        J_{7,\mathrm{SP}}
        &\simeq
        -\beta
        e^{-i\beta x/4}
        e^{-i\beta t/4}
        e^{i\pi/4}
        \sqrt{\frac{\pi}{\beta t}}\\
        &=
        -\sqrt{\frac{\pi\beta}{t}}
        e^{-i\beta(t+x)/4+i\pi/4}.

\end{aligned}
```

Taking the real part,

```math
\boxed{
            I_{7,\mathrm{SP}}(x,t)
            \simeq
            -\sqrt{\frac{\pi\beta}{t}}
            \cos\left[
            \frac{\beta(t+x)}{4}
            -\frac{\pi}{4}
            \right]
        }
```

as $t\rightarrow\infty$.

Thus,

```math
I_{7,\mathrm{SP}}
        =
        O(t^{-1/2}).
```

### Pole contribution to $I_7$

The pole in the original $k$-variable occurs at

```math
k=\beta.
```

Introduce

```math
k=\beta+h,
        \qquad
        h\rightarrow0.
```

Define

```math
D(k)
        =
        k-\sqrt{\beta k}.
```

Since

```math
D(\beta)=0,
```

Taylor expansion about $k=\beta$ gives

```math
D(k)
        =
        D'(\beta)h
        +
        \frac{D''(\beta)}{2}h^2
        +\cdots .
```

Now,

```math
D'(k)
        =
        1-\frac{\sqrt{\beta}}{2\sqrt{k}},
```

so that

```math
D'(\beta)
        =
        \frac12.
```

Also,

```math
D''(k)
        =
        \frac{\sqrt{\beta}}{4k^{3/2}},
```

and therefore

```math
D''(\beta)
        =
        \frac{1}{4\beta}.
```

Hence,

```math
D(k)
        =
        \frac{h}{2}
        +
        \frac{h^2}{8\beta}
        +\cdots ,
```

and to leading order,

```math
\boxed{
            D(k)\simeq\frac{h}{2}.
        }
```

Therefore,

```math
\begin{aligned}
        \frac{k}{k-\sqrt{\beta k}}
        &=
        \frac{\beta+h}{D(k)}\\
        &\simeq
        \frac{\beta+h}{h/2}\\
        &=
        \frac{2\beta}{h}+2.

\end{aligned}
```

The singular part is consequently

```math
\boxed{
            \frac{k}{k-\sqrt{\beta k}}
            \simeq
            \frac{2\beta}{h}.
        }
```

#### Phase near the pole

The phase is

```math
tD(k)-kx.
```

Using $D(k)\simeq h/2$ and $k=\beta+h$, the phase near the pole $k=\beta$ becomes

```math
\begin{aligned}
        tD(k)-kx
        &\simeq
        \frac{th}{2}
        -(\beta+h)x\\
        &=
        -\beta x
        +
        \left(
        \frac{t}{2}-x
        \right)h.

\end{aligned}
```

To evaluate the contribution from the neighbourhood of the pole $k=\beta$, consider a small symmetric interval

```math
\beta-\delta<k<\beta+\delta,
        \qquad
        \delta>0,
```

where $\delta$ is sufficiently small for the above local approximations to remain valid. Since

```math
k=\beta+h,
        \qquad
        h=k-\beta,
```

the limits of this local interval transform according to

```math
k=\beta-\delta
        \quad\Longrightarrow\quad
        h=-\delta,
```

and

```math
k=\beta+\delta
        \quad\Longrightarrow\quad
        h=\delta.
```

Thus, the neighbourhood $\beta-\delta<k<\beta+\delta$ of the pole $k=\beta$ is mapped to the symmetric interval

```math
-\delta<h<\delta,
```

with the pole itself located at $h=0$.

Therefore, using the leading-order approximations for the amplitude and phase near $k=\beta$, the local pole contribution is

```math
I_{7,\mathrm{pole}}
        \simeq
        2\beta\,
        \operatorname{PV}
        \int_{-\delta}^{\delta}
        \frac{
            \cos\left[
            -\beta x+
            \left(
            \dfrac{t}{2}-x
            \right)h
            \right]
        }{h}\,dh .
```

The Cauchy principal value is required because the integrand is singular at $h=0$. For an integral over this symmetric neighbourhood, it is defined as

```math
\operatorname{PV}
        \int_{-\delta}^{\delta}
        f(h)\,dh
        =
        \lim_{\epsilon\to0^+}
        \left[
        \int_{-\delta}^{-\epsilon}
        f(h)\,dh
        +
        \int_{\epsilon}^{\delta}
        f(h)\,dh
        \right].
```

Here, $\delta$ specifies the size of the local neighbourhood around the pole, whereas $\epsilon$ denotes the small symmetric interval excluded about $h=0$ in defining the Cauchy principal value.

Define

```math
A=-\beta x,
        \qquad
        B=\frac{t}{2}-x.
```

Then

```math
I_{7,\mathrm{pole}}
        \simeq
        2\beta \operatorname{PV}
        \int_{-\delta}^{\delta}
        \frac{\cos(A+Bh)}{h}\,dh.
```

Using

```math
\cos(A+Bh)
        =
        \cos A\cos(Bh)
        -
        \sin A\sin(Bh),
```

we obtain

```math
\begin{aligned}
        I_{7,\mathrm{pole}}
        \simeq
        2\beta
        \Bigg[
        &
        \cos A
        \operatorname{PV}
        \int_{-\delta}^{\delta}
        \frac{\cos(Bh)}{h}\,dh\\
        &
        -
        \sin A \operatorname{PV}
        \int_{-\delta}^{\delta}
        \frac{\sin(Bh)}{h}\,dh
        \Bigg].

\end{aligned}
```

The first integral vanishes because $\cos(Bh)/h$ is odd:

```math
\operatorname{PV}
        \int_{-\delta}^{\delta}
        \frac{\cos(Bh)}{h}\,dh
        =
        0.
```

For the second integral, set

```math
u=Bh.
```

Then

```math
\int_{-\delta}^{\delta}
        \frac{\sin(Bh)}{h}\,dh
        =
        \int_{-B\delta}^{B\delta}
        \frac{\sin u}{u}\,du .
```

For fixed $x$,

```math
B=\frac{t}{2}-x
        \rightarrow\infty
        \qquad
        (t\rightarrow\infty).
```

Therefore,

```math
\int_{-B\delta}^{B\delta}
        \frac{\sin u}{u}\,du
        \longrightarrow
        \int_{-\infty}^{\infty}
        \frac{\sin u}{u}\,du
        =
        \pi.
```

Thus,

```math
I_{7,\mathrm{pole}}
        \simeq
        -2\pi\beta\sin A.
```

Since $A=-\beta x$,

```math
\boxed{
            I_{7,\mathrm{pole}}(x,t)
            \simeq
            2\pi\beta\sin(\beta x).
        }
```

This contribution is independent of $t$.

### Large-time result for $I_7$

Combining the stationary-phase and pole contributions,

```math
\begin{aligned}
                I_7(x,t)
                \sim{}&
                2\pi\beta\sin(\beta x)
                \\
                &-
                \sqrt{\frac{\pi\beta}{t}}
                \cos\left[
                \frac{\beta(t+x)}{4}
                -\frac{\pi}{4}
                \right],
                \qquad
                t\rightarrow\infty .
            \end{aligned}
```

Hence,

```math
I_{7,\mathrm{pole}}=O(1),
        \qquad
        I_{7,\mathrm{SP}}=O(t^{-1/2}),
```

and

```math
\boxed{
            \lim_{t\rightarrow\infty}
            I_7(x,t)
            =
            2\pi\beta\sin(\beta x).
        }
```

### Numerical verification of $I_7$

For the numerical comparison we take

```math
\rho_r=0.001,
        \qquad
        \beta=\frac{1-\rho_r}{1+\rho_r},
        \qquad
        x=1.
```

The original $I_7$ integral is evaluated numerically in the Cauchy principal-value sense by excluding a small symmetric neighbourhood of the pole $k=\beta$.

The resulting asymptotic form gives

```math
\boxed{
            \sqrt{\frac{t}{\pi\beta}}
            \left[
            I_7(x,t)-2\pi\beta\sin(\beta x)
            \right]
            \sim
            -\cos\left[
            \frac{\beta(t+x)}{4}
            -\frac{\pi}{4}
            \right].
        }
```

Thus, the compensated comparison verifies the predicted $t^{-1/2}$ amplitude together with its phase and frequency.

## Large-time asymptotics of the non-singular integral $I_8(x,t)$

The defining expression for $I_8(x,t)$ may be rewritten as

```math
I_8(x,t)
        =
        \int_0^\infty
        \frac{k}{k+\sqrt{\beta k}}
        \cos\left\{
        t\left(k+\sqrt{\beta k}\right)-kx
        \right\}\,dk .
```

### Transformation of $I_8$

Introduce

```math
k=\beta z^2.
```

Then

```math
\sqrt{\beta k}=\beta z,
        \qquad
        dk=2\beta z\,dz,
```

and

```math
k+\sqrt{\beta k}
        =
        \beta z(z+1).
```

Hence,

```math
\frac{k}{k+\sqrt{\beta k}}\,dk
        =
        \frac{2\beta z^2}{z+1}\,dz.
```

The phase becomes

```math
\begin{aligned}
        t\left(k+\sqrt{\beta k}\right)-kx
        &=
        t\left(\beta z^2+\beta z\right)
        -\beta xz^2\\
        &=
        \beta
        \left[
        t(z^2+z)-xz^2
        \right].

\end{aligned}
```

Therefore,

```math
\boxed{
            I_8(x,t)
            =
            2\beta
            \int_0^\infty
            \frac{z^2}{z+1}
            \cos\left\{
            \beta
            \left[
            t(z^2+z)-xz^2
            \right]
            \right\}\,dz.
        }
```

### Absence of poles and stationary points

Since

```math
z+1>0
        \qquad
        (z\ge0),
```

there is no pole on the integration interval.

Define the leading $t$-dependent phase as

```math
\phi_8(z)
        =
        \beta(z^2+z).
```

A stationary point satisfies

```math
\phi_8'(z_s)=0.
```

Since

```math
\phi_8'(z)
        =
        \beta(2z+1),
```

the only stationary point is

```math
z_s=-\frac12,
```

which lies outside the integration interval $z\ge0$.

Thus, $I_8$ has neither a pole nor a stationary point on the path of integration. Its leading large-time behaviour is therefore controlled by the endpoint z=0.

### Endpoint scaling

The phase function is

```math
\phi_8(z)
        =
        \beta z+\beta z^2.
```

Near $z=0$,

```math
\phi_8(z)
        =
        \beta z+O(z^2),
```

so that

```math
t\phi_8(z)
        =
        \beta t z+\beta t z^2
        \sim
        \beta t z.
```

Thus, the leading variation of the cosine near the endpoint is controlled by the combination $\beta t z$. In order to resolve the region that contributes to the integral as $t\to\infty$, this combination must remain of order unity, i.e.,

```math
\beta t z=O(1).
```

Consequently, the relevant endpoint region has the scaling

```math
z
        =
        O\left(\frac{1}{\beta t}\right).
```

This motivates the introduction of the stretched endpoint variable

```math
\boxed{
            y=\beta t z
        }.
```

Equivalently,

```math
z=\frac{y}{\beta t},
        \qquad
        dz=\frac{dy}{\beta t}.
```

The amplitude transforms as

```math
\begin{aligned}
        2\beta\frac{z^2}{1+z}\,dz
        &=
        \frac{2}{\beta^2t^3}
        y^2
        \left(
        1+\frac{y}{\beta t}
        \right)^{-1}
        dy .

\end{aligned}
```

For $t\rightarrow\infty$,

```math
\left(
        1+\frac{y}{\beta t}
        \right)^{-1}
        =
        1-\frac{y}{\beta t}
        +
        O(t^{-2}).
```

The phase becomes

```math
\begin{aligned}
        \beta
        \left[
        t(z^2+z)-xz^2
        \right]
        &=
        y+
        \frac{y^2}{\beta t}
        \left(
        1-\frac{x}{t}
        \right).

\end{aligned}
```

Hence,

```math
\begin{aligned}
        I_8(x,t)
        =
        \frac{2}{\beta^2t^3}
        \int_0^\infty
        y^2
        \left(
        1+\frac{y}{\beta t}
        \right)^{-1}
        \cos\left[
        y+
        \frac{y^2}{\beta t}
        \left(
        1-\frac{x}{t}
        \right)
        \right]dy .

\end{aligned}
```

### Complex representation and asymptotic expansion

Using

```math
\cos\theta
        =
        \Re \left(e^{i\theta}\right),
```

The stretched-variable form becomes

```math
\begin{aligned}
        I_8(x,t)
        =
        \frac{2}{\beta^2t^3}
        \Re
        \int_0^\infty
        y^2e^{iy}
        \left(
        1+\frac{y}{\beta t}
        \right)^{-1}
        \exp\left[
        \frac{iy^2}{\beta t}
        \left(
        1-\frac{x}{t}
        \right)
        \right]dy.

\end{aligned}
```

For fixed $x$,

```math
\frac{x}{t}
        =
        O(t^{-1}),
```

and therefore

```math
\exp\left[
        \frac{iy^2}{\beta t}
        \left(
        1-\frac{x}{t}
        \right)
        \right]
        =
        1+
        \frac{iy^2}{\beta t}
        +
        O(t^{-2}).
```

Multiplying the amplitude and phase expansions,

```math
\left(
        1+\frac{y}{\beta t}
        \right)^{-1}
        \exp\left[
        \frac{iy^2}{\beta t}
        \left(
        1-\frac{x}{t}
        \right)
        \right]
        =
        1
        -\frac{y}{\beta t}
        +\frac{iy^2}{\beta t}
        +
        O(t^{-2}).
```

Substituting into the complex representation,

```math
\begin{aligned}
        I_8(x,t)
        \sim
        \frac{2}{\beta^2t^3}
        \Bigg[
        &
        \Re
        \int_0^\infty
        y^2e^{iy}\,dy\\
        &
        -
        \frac{1}{\beta t}
        \Re
        \int_0^\infty
        y^3e^{iy}\,dy\\
        &
        +
        \frac{1}{\beta t}
        \Re
        \left(
        i\int_0^\infty
        y^4e^{iy}\,dy
        \right)
        \Bigg].

\end{aligned}
```

### Evaluation of the oscillatory integrals

The resulting oscillatory integrals are understood in the Abel-regularized sense:

```math
\int_0^\infty
        y^n e^{iy}\,dy
        \equiv
        \lim_{\epsilon\rightarrow0^+}
        \int_0^\infty
        y^n e^{-(\epsilon-i)y}\,dy.
```

For $\epsilon>0$,

```math
\int_0^\infty
        y^n e^{-(\epsilon-i)y}\,dy
        =
        \frac{n!}{(\epsilon-i)^{n+1}},
```

and therefore

```math
\boxed{
            \int_0^\infty
            y^n e^{iy}\,dy
            =
            \frac{n!}{(-i)^{n+1}}.
        }
```

For $n=2$,

```math
\int_0^\infty
        y^2e^{iy}\,dy
        =
        -2i,
```

so that

```math
\Re
        \int_0^\infty
        y^2e^{iy}\,dy
        =
        0.
```

For $n=3$,

```math
\int_0^\infty
        y^3e^{iy}\,dy
        =
        6,
```

and hence

```math
\Re
        \int_0^\infty
        y^3e^{iy}\,dy
        =
        6.
```

For $n=4$,

```math
\int_0^\infty
        y^4e^{iy}\,dy
        =
        24i.
```

Therefore,

```math
\Re
        \left[
        i\int_0^\infty
        y^4e^{iy}\,dy
        \right]
        =
        -24.
```

### Large-time result for $I_8$

Substituting these evaluated integrals,

```math
\begin{aligned}
        I_8(x,t)
        &\sim
        \frac{2}{\beta^2t^3}
        \left[
        0
        -\frac{6}{\beta t}
        -\frac{24}{\beta t}
        \right]\\
        &=
        -\frac{60}{\beta^3t^4}.

\end{aligned}
```

Thus,

```math
I_8(x,t)
            \sim
            -\frac{60}{\beta^3t^4},
            \qquad
            t\rightarrow\infty,
            \quad
            x=\text{fixed}.
```

Therefore,

```math
I_8(x,t)=O(t^{-4}).
```

The leading contribution is independent of fixed $x$; the $x$-dependence enters only at higher order.

### Contour-rotated representation and numerical verification of $I_8$

To verify the large-time asymptotic result derived above, the integral $I_8(x,t)$ is evaluated numerically. Rather than directly evaluating the original oscillatory integral along the positive real axis, it is convenient to first transform it into an exponentially damped integral by rotating the contour into the positive imaginary axis.

The transformed integral is

```math
I_8(x,t)
        =
        2\beta
        \int_0^\infty
        \frac{z^2}{z+1}
        \cos\left\{
        \beta
        \left[
        (t-x)z^2+tz
        \right]
        \right\}
        dz .
```

Using

```math
\cos\theta
        =
        \Re\left(e^{i\theta}\right),
```

This integral can be written as

```math
I_8(x,t)
        =
        \Re\left[J_8(x,t)\right],
```

where

```math
J_8(x,t)
        =
        2\beta
        \int_0^\infty
        \frac{z^2}{z+1}
        \exp\left\{
        i\beta
        \left[
        (t-x)z^2+tz
        \right]
        \right\}
        dz .
```

#### Rotation of the integration contour

For $t>x$, consider the quarter-circle contour in the first quadrant shown schematically in the omitted contour figure. The contour consists of the positive real axis, a quarter-circle of radius $R$, and the positive imaginary axis.

Writing

```math
z=Re^{i\theta},
        \qquad
        0\leq\theta\leq\frac{\pi}{2},
```

the modulus of the quadratic exponential factor is

```math
\begin{aligned}
        \left|
        e^{\,i\beta(t-x)z^2}
        \right|
        &=
        \left|
        e^{\,i\beta(t-x)R^2
            (\cos 2\theta+i\sin 2\theta)}
        \right|\\
        &=
        e^{-\beta(t-x)R^2\sin 2\theta}.

\end{aligned}
```

Similarly,

```math
\begin{aligned}
        \left|
        e^{\,i\beta tz}
        \right|
        &=
        \left|
        e^{\,i\beta tR(\cos\theta+i\sin\theta)}
        \right|\\
        &=
        e^{-\beta tR\sin\theta}.

\end{aligned}
```

Since

```math
\sin\theta\geq0,
        \qquad
        \sin2\theta\geq0,
        \qquad
        0\leq\theta\leq\frac{\pi}{2},
```

and $t>x$, both exponential factors are non-growing on the first-quadrant arc. The linear exponential provides exponential decay away from the positive real axis, while the quadratic exponential also decays for $0<\theta<\pi/2$. Consequently, the contribution from the quarter-circle vanishes as $R\to\infty$.

The integrand has a pole at $z=-1$, which lies outside the first-quadrant contour. Therefore, no residue is enclosed.

Let

```math
f(z)
        =
        2\beta
        \frac{z^2}{z+1}
        e^{\,i\beta(t-x)z^2}
        e^{\,i\beta tz}.
```

Since $f(z)$ is analytic inside the first-quadrant contour, Cauchy’s theorem gives

```math
\oint_C f(z)\,dz=0.
```

Decomposing the contour into the positive real axis $\Gamma_1$, the quarter-circle $\Gamma_2$, and the imaginary-axis segment $\Gamma_3$,

```math
\int_{\Gamma_1}f(z)\,dz
        +
        \int_{\Gamma_2}f(z)\,dz
        +
        \int_{\Gamma_3}f(z)\,dz
        =
        0.
```

As $R\to\infty$,

```math
\int_{\Gamma_2}f(z)\,dz\longrightarrow0.
```

The segment $\Gamma_3$ is traversed downwards, from $i\infty$ to the origin. Parameterizing this segment by

```math
z=is,
        \qquad
        dz=i\,ds,
```

its orientation corresponds to $s:\infty\rightarrow0$. Hence,

```math
\int_{\Gamma_3}f(z)\,dz
        =
        \int_{\infty}^{0}f(is)\,i\,ds
        =
        -i\int_0^\infty f(is)\,ds.
```

Cauchy's theorem then gives

```math
\int_0^\infty f(z)\,dz
        =
        i\int_0^\infty f(is)\,ds.
```

Under the transformation $z=is$,

```math
z^2=-s^2,
        \qquad
        z+1=1+is,
```

and the exponential factors become

```math
e^{\,i\beta(t-x)z^2}
        =
        e^{-i\beta(t-x)s^2},
```

and

```math
e^{\,i\beta tz}
        =
        e^{-\beta ts}.
```

Therefore,

```math
\begin{aligned}
        f(is)
        &=
        2\beta
        \frac{(is)^2}{1+is}
        e^{-i\beta(t-x)s^2}
        e^{-\beta ts}\\
        &=
        \frac{-2\beta s^2}{1+is}
        e^{-\beta ts}
        e^{-i\beta(t-x)s^2}.

\end{aligned}
```

Thus, the complex integral can be written as

```math
\boxed{
            J_8(x,t)
            =
            \int_0^\infty
            \frac{-2i\beta s^2}{1+is}
            e^{-\beta ts}
            e^{-i\beta(t-x)s^2}
            \,ds
        }.
```

Taking the real part gives the required representation of $I_8$,

```math
\boxed{
            I_8(x,t)
            =
            \Re
            \left[
            \int_0^\infty
            \frac{-2i\beta s^2}{1+is}
            e^{-\beta ts}
            e^{-i\beta(t-x)s^2}
            \,ds
            \right]
        }.
```

This contour-rotated representation is the form used for the numerical evaluation in MATLAB. In contrast to the original real-axis representation, the factor

```math
e^{-\beta ts}
```

provides exponential damping as $s\to\infty$, making the contour-rotated form particularly convenient for numerical quadrature at large $t$.

#### Comparison with the large-time asymptotic result

For the numerical comparison, the parameters are chosen as

```math
\rho_r=0.001,
        \qquad
        \beta=\frac{1-\rho_r}{1+\rho_r},
        \qquad
        x=1.
```

The contour-rotated integral is evaluated numerically and compared with the large-time asymptotic result derived above,

```math
I_8(x,t)
        \sim
        -\frac{60}{\beta^3t^4},
        \qquad
        t\to\infty.
```

This asymptotic comparison further implies

```math
\boxed{
            \beta^3t^4 I_8(x,t)
            \longrightarrow
            -60,
            \qquad
            t\to\infty
        }.
```

Accordingly, the log--log comparison verifies the predicted $t^{-4}$ decay, while the compensated comparison provides a more stringent test of the leading-order coefficient. The convergence of the compensated numerical result to $-60$ confirms both the decay exponent and the coefficient in the leading large-time asymptotic expression.

## Combined large-time behaviour

The two transient integrals have the leading asymptotic forms

```math
\begin{aligned}
                I_7(x,t)
                \sim{}&
                2\pi\beta\sin(\beta x)
                -
                \sqrt{\frac{\pi\beta}{t}}
                \cos\left[
                \frac{\beta(t+x)}{4}
                -\frac{\pi}{4}
                \right],
                \\[4pt]
                I_8(x,t)
                \sim{}&
                -\frac{60}{\beta^3t^4}.
            \end{aligned}
```

Therefore,

```math
I_7(x,t)+I_8(x,t)
        \sim
        2\pi\beta\sin(\beta x)
        -
        \sqrt{\frac{\pi\beta}{t}}
        \cos\left[
        \frac{\beta(t+x)}{4}
        -\frac{\pi}{4}
        \right]
        -
        \frac{60}{\beta^3t^4}.
```

Substitution into the transient decomposition gives

```math
\begin{aligned}
                \frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
                \sim
                -\frac{1}{2\pi(1-\rho_r)}
                \Bigg[
                &
                2\pi\beta\sin(\beta x)
                \\
                &-
                \sqrt{\frac{\pi\beta}{t}}
                \cos\left[
                \frac{\beta(t+x)}{4}
                -\frac{\pi}{4}
                \right]
                \\
                &-
                \frac{60}{\beta^3t^4}
                \Bigg].
            \end{aligned}
```

The large-time asymptotic structure may therefore be summarized as

```math
I_{7,\mathrm{pole}}=O(1),
        \qquad
        I_{7,\mathrm{SP}}=O(t^{-1/2}),
        \qquad
        I_8=O(t^{-4}).
```

Thus, the pole of $I_7$ produces the time-independent limiting contribution, while the dominant decaying transient is generated by the stationary-phase contribution of $I_7$. The contribution from $I_8$ decays much more rapidly.

In particular,

```math
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
            =
            -\frac{\beta}{1-\rho_r}
            \sin(\beta x)
            +
            O(t^{-1/2}),
            \qquad
            t\rightarrow\infty.
```

Since

```math
\beta
        =
        \frac{1-\rho_r}{1+\rho_r},
```

the limiting time-independent contribution may also be written as

```math
\boxed{
            \frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
            =
            -\frac{1}{1+\rho_r}
            \sin(\beta x)
            +
            O(t^{-1/2}).
        }
```

## Summary

For fixed $x$, the large-time asymptotic behaviour of the two-fluid transient integrals in the pure-gravity limit is

```math
\boxed{
            I_7(x,t)
            \sim
            2\pi\beta\sin(\beta x)
            -
            \sqrt{\frac{\pi\beta}{t}}
            \cos\left[
            \frac{\beta(t+x)}{4}
            -\frac{\pi}{4}
            \right],
        }
```

and

```math
\boxed{
            I_8(x,t)
            \sim
            -\frac{60}{\beta^3t^4}.
        }
```

The singular integral $I_7$ contains both a non-decaying pole contribution and an oscillatory stationary-phase contribution that decays as $t^{-1/2}$. In contrast, $I_8$ has neither a pole nor a stationary point on the integration interval, and its asymptotic behaviour is controlled by the endpoint $z=0$, resulting in the much faster $t^{-4}$ decay.
