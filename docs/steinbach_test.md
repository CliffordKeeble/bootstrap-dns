# The Steinbach Line Test

## Background

In incompressible turbulence, the velocity gradient tensor A_ij = ∂u_i/∂x_j encodes the local state of the flow. Its symmetric part S = (A + Aᵀ)/2 is the strain-rate tensor, whose three eigenvalues λ₁ ≥ λ₂ ≥ λ₃ (with λ₁ + λ₂ + λ₃ = 0 by incompressibility) characterise the local stretching geometry.

The **Q-R invariants** of the velocity gradient tensor are standard diagnostics in turbulence research:

    Q = -½ tr(A²)     (second invariant)
    R = -⅓ tr(A³)     (third invariant)

Their joint PDF has a characteristic "teardrop" shape first analysed by Vieillefosse (1982, 1984). The discriminant line D = 27R²/4 + Q³ = 0 separates real and complex eigenvalue regions.

## The Bootstrap Prediction

Paper 165 ("The Heptagonal Constraint") conjectures that the strain-rate tensor in developed turbulence is constrained by the **Steinbach polynomial**:

    t³ − (7/3)t + 7/27 = 0

This cubic arises naturally from the icosahedral geometry of S³/2I (the Poincaré homology sphere). Its coefficients contain only primes {3, 7}, and its discriminant is Δ = 49 = 7² — the smallest prime-square discriminant for a cyclic cubic field.

### The Steinbach Roots

    τ₁ ≈ +1.4686  (extensional)
    τ₂ ≈ +0.1117  (intermediate)
    τ₃ ≈ −1.5803  (compressive)

### Predicted Eigenvalue Ratios

    λ₂/λ₁ ≈ 0.076  (intermediate/extensional)
    λ₃/λ₁ ≈ −1.076 (compressive/extensional)

### The Steinbach Line in (Q_S, R_S) Space

For the strain-rate tensor (symmetric part only), the Steinbach polynomial defines a specific curve in (Q_S, R_S) space:

    R_S / (-Q_S)^{3/2} = constant = (7/27) / (7/3)^{3/2}

If the Steinbach constraint holds, the joint PDF of (Q_S, R_S) should cluster along this line, particularly in the high-strain tail.

## What To Look For

1. **Eigenvalue ratios:** The DNS distributions of λ₂/λ₁ and λ₃/λ₁ should peak near the Steinbach values. At low Re, the distributions are broad. At higher Re (≥ 1600), the peaks should sharpen toward the Steinbach prediction.

2. **Q_S-R_S clustering:** The high-strain tail of the joint PDF should follow the Steinbach line, not just the generic discriminant Δ = 0.

3. **Reynolds number dependence:** The Steinbach constraint is geometric — it should become more pronounced as Re increases and the turbulence becomes more developed.

## Known Limitations

- At Re ~ 125 (32³ grid), the Taylor-Green vortex hasn't fully transitioned to developed turbulence. The eigenvalue distributions are broad and the Steinbach signal, if present, is obscured.
- At Re ~ 400 (64³), partial turbulence develops. The Q-R teardrop is well-formed but the inertial range is short.
- Re ≥ 1600 (128³+) is needed for a meaningful test with a developed inertial range.

## Connection to Paper 166

Paper 166 ("The Cascade Exponent") derives Kolmogorov's −5/3 from the Newton-Girard power sums of the same Steinbach polynomial: the prime 5 emerges at exactly p₅ = −5 × 7²/3⁴. The DNS energy spectrum E(k) ∝ k^{−5/3} is a sanity check — well established experimentally, but here connected to the Steinbach algebra.

## References

- Vieillefosse, P. (1982). "Local interaction between vorticity and shear in a perfect incompressible fluid." *J. Phys. France*, 43(6), 837-842.
- Cantwell, B. J. (1992). "Exact solution of a restricted Euler equation for the velocity gradient tensor." *Physics of Fluids A*, 4(4), 782-793.
- Steinbach, P. (2015). "Field extensions and minimal discriminants." Tables of number fields with prescribed ramification.

---

*Bootstrap Universe Programme — Woodbridge, Suffolk, UK*
