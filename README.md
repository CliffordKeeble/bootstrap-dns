# Bootstrap DNS

**Pseudospectral Direct Numerical Simulation for testing predictions from the Bootstrap Universe Programme**

A lightweight Python DNS solver for incompressible Navier-Stokes, designed to extract strain-rate diagnostics that test specific predictions from Papers 163–166 of the Bootstrap Programme.

## What This Tests

The Bootstrap Programme derives turbulence properties from the geometry of the Poincaré homology sphere S³/2I and the Steinbach polynomial t³ − 7t/3 + 7/27 = 0. This solver provides computational evidence for or against those predictions.

### Paper 165 — The Heptagonal Constraint
**Prediction:** The strain-rate tensor in developed turbulence is constrained by the Steinbach polynomial. The joint PDF of the Q-R invariants should cluster along the Steinbach discriminant curve Δ = 49A⁶.

**Test:** Extract the velocity gradient tensor ∂uᵢ/∂xⱼ at every grid point, compute strain-rate eigenvalues and Q-R invariants, compare with the Steinbach constraint.

### Paper 166 — The Cascade Exponent
**Prediction:** Kolmogorov's −5/3 exponent derives from the Newton-Girard power sums of the Steinbach polynomial: the prime 5 emerges at exactly p₅, giving −(quintic emergence order)/(eigenvalue count) = −5/3.

**Test:** Verify k⁻⁵/³ energy spectrum develops in the inertial range (sanity check — well established, but confirms the solver is correct).

### Specific Diagnostics
| Test | What We Measure | What Paper 165/166 Predicts |
|------|----------------|---------------------------|
| Q-R joint PDF | Shape of (Q,R) scatter | Clustering along Steinbach discriminant |
| Eigenvalue ratios | λ₂/λ₁ and λ₃/λ₁ | Steinbach: λ₂/λ₁ ≈ 0.076, λ₃/λ₁ ≈ −1.076 |
| Energy spectrum | E(k) vs k | k⁻⁵/³ inertial range |
| Vieillefosse tail | Q-R teardrop shape | Standard benchmark (validates solver) |

## Quick Start

```bash
# Clone
git clone https://github.com/DrCliffKeeble/bootstrap-dns.git
cd bootstrap-dns

# Install dependencies
pip install -r requirements.txt

# Run with defaults (32³, Re~125 — runs in ~1 min)
python scripts/dns_solver.py

# Run production (128³, Re~1600 — needs decent hardware, ~hours)
python scripts/dns_solver.py --N 128 --Re 1600 --T 20.0 --dt 0.001

# Generate diagnostic plots from saved data
python scripts/plot_diagnostics.py results/
```

## Configuration

All parameters are set via command-line arguments:

| Argument | Default | Description |
|----------|---------|-------------|
| `--N` | 32 | Grid points per dimension (32, 64, 128, 256) |
| `--Re` | 125 | Reynolds number |
| `--T` | 5.0 | Final time |
| `--dt` | 0.004 | Time step |
| `--outdir` | `results/` | Output directory |
| `--print-every` | 50 | Diagnostic print interval (steps) |

### Recommended Configurations

| Config | Grid | Re | Time | Hardware | Purpose |
|--------|------|-----|------|----------|---------|
| Quick test | 32³ | 125 | 5.0 | Any laptop | Verify solver works |
| Development | 64³ | 400 | 10.0 | Laptop | Strain statistics, partial spectrum |
| Production | 128³ | 1600 | 20.0 | Desktop/workstation | Full inertial range, Steinbach test |
| Research | 256³ | 6000 | 30.0 | HPC cluster | Publication-quality statistics |

## Output

The solver produces:

```
results/
├── dns_state.npz          # Final velocity field (Fourier coefficients)
├── diagnostics.json       # Time series: energy, enstrophy, max|u|
├── spectrum.npz           # Energy spectra at saved times
├── strain_data.npz        # Q, R, Q_S, R_S, eigenvalues (at peak dissipation)
├── fig1_energy.png        # Energy and enstrophy evolution
├── fig2_spectrum.png      # Energy spectrum with k⁻⁵/³ reference
├── fig3_QR.png            # Q-R joint PDF with Steinbach line
└── fig4_eigenvalues.png   # Strain eigenvalue ratio distributions
```

## Method

- **Spatial discretisation:** Fourier pseudospectral with 2/3 dealiasing
- **Temporal integration:** Classical RK4
- **Nonlinear term:** Rotational (Lamb) form: −(ω × u) with Leray projection
- **Initial condition:** Taylor-Green vortex: u = sin(x)cos(y)cos(z), v = −cos(x)sin(y)cos(z), w = 0
- **Boundary conditions:** Triply periodic on [0, 2π]³

## Requirements

- Python ≥ 3.9
- NumPy
- Matplotlib
- SciPy (optional, for advanced post-processing)

## Context

This code is part of the [Bootstrap Universe Programme](https://zenodo.org/communities/bootstrap-universe/), a mathematical physics research programme deriving physical constants and structures from the geometry of S³/2I.

**Papers tested by this code:**
- Paper 165: "The Heptagonal Constraint" — Steinbach polynomial constrains strain tensor
- Paper 166: "The Cascade Exponent" — Kolmogorov −5/3 from Newton-Girard power sums

**Author:** Dr. Clifford Keeble, PhD  
**ORCID:** [0009-0003-6828-2155](https://orcid.org/0009-0003-6828-2155)  
**Location:** Woodbridge, Suffolk, UK

## Licence

MIT — see [LICENCE](LICENCE)

## Status

🟡 **Active development.** The solver is validated against standard Taylor-Green benchmarks. The Steinbach diagnostic is preliminary — higher Reynolds numbers (Re ≥ 1600) are needed for meaningful strain statistics. Contributions and independent verification welcome.

---

*"Bootstrap is telling us it's about integers."*

🐕☕⬡
