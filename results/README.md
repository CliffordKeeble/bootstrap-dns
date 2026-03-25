# Results

This directory is populated by running `scripts/dns_solver.py`.

Output files:
- `diagnostics.json` — Time series of energy, enstrophy, parameters
- `spectrum.npz` — Energy spectra at saved times
- `strain_data.npz` — Q-R invariants and eigenvalues at peak dissipation
- `dns_state.npz` — Final velocity field (Fourier coefficients, for restart)
- `fig1_energy.png` — Energy and enstrophy evolution
- `fig2_spectrum.png` — Energy spectrum with k⁻⁵/³ reference
- `fig3_QR.png` — Q-R joint PDF with Steinbach line
- `fig4_eigenvalues.png` — Strain eigenvalue ratio distributions
