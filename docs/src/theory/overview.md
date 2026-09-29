# Theory

This section contains the mathematical formulation and numerical implementation of the forced interfacial wave IVP, following the notation of [Kadari et al. (2026)](https://arxiv.org/abs/2605.12254).

## Problem description and nondimensionalisation

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/Fig3.png" alt="Pressure forcing at a two-fluid interface" style="max-width:100%;">
  <figcaption style="text-align:center;"><strong>Figure 2.</strong> A point force (red arrow) of strength <em>F̃<sub>0</sub></em> acts at <em>t̃</em> &gt; 0 at the interface of two fluids of density <em>&rho;<sub>u</sub></em> and <em>&rho;<sub>l</sub></em>, both moving with speed <em>U</em> (as seen in the co-moving frame). The linearised IVP predicts how waves develop at the interface, in time.</figcaption>
</figure>
```

A localised pressure forcing $\tilde{p}_e = \tilde{F}_0\,\delta(\tilde{x})$ acts at the interface of two inviscid, incompressible, irrotational fluids of infinite depth. Both streams move at uniform speed $U$ rightwards. In the absence of forcing, the interface is flat at $\tilde{z} = 0$.

After nondimensionalisation using

```math
l_c = \frac{U^2}{\tilde{g}}, \quad t_c = \frac{U}{\tilde{g}}, \quad p_c = \rho_l U^2,
```

the key dimensionless groups are

```math
\alpha = \frac{\tilde{g} T}{\rho_l U^4}, \quad
\rho_r = \frac{\rho_u}{\rho_l}, \quad
\beta = \frac{1-\rho_r}{1+\rho_r}.
```

The parameter $\alpha$ measures the relative importance of surface tension. 
