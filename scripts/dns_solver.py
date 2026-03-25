#!/usr/bin/env python3
"""
Bootstrap DNS Solver
====================
Pseudospectral solver for incompressible Navier-Stokes equations.
Taylor-Green vortex initial condition.
Strain-rate diagnostics for testing Bootstrap Programme Papers 165-166.

Dr. Clifford Keeble & Claude
Bootstrap Universe Programme, March 2026
"""

import argparse
import json
import os
import sys
from time import time

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq


# ============================================================
# CONFIGURATION
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Bootstrap DNS — pseudospectral Navier-Stokes solver",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dns_solver.py                          # Quick test (32³, Re~125)
  python dns_solver.py --N 64 --Re 400         # Development run
  python dns_solver.py --N 128 --Re 1600 --T 20 --dt 0.001  # Production
        """)
    p.add_argument("--N", type=int, default=32,
                   help="Grid points per dimension (default: 32)")
    p.add_argument("--Re", type=float, default=125,
                   help="Reynolds number (default: 125)")
    p.add_argument("--T", type=float, default=5.0,
                   help="Final time (default: 5.0)")
    p.add_argument("--dt", type=float, default=0.004,
                   help="Time step (default: 0.004)")
    p.add_argument("--outdir", type=str, default="results",
                   help="Output directory (default: results/)")
    p.add_argument("--print-every", type=int, default=50,
                   help="Print diagnostics every N steps (default: 50)")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip plot generation (data files only)")
    return p.parse_args()


# ============================================================
# GRID SETUP
# ============================================================

def setup_grid(N):
    """Initialise wavenumbers, dealiasing mask, and coordinate arrays."""
    L = 2 * np.pi
    dx = L / N
    x = np.linspace(0, L, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

    k1d = fftfreq(N, d=1.0 / N)
    KX, KY, KZ = np.meshgrid(k1d, k1d, k1d, indexing='ij')
    K2 = KX**2 + KY**2 + KZ**2
    K2_safe = K2.copy()
    K2_safe[0, 0, 0] = 1.0

    # 2/3 dealiasing
    kmax = N // 3
    dealias = ((np.abs(KX) <= kmax) &
               (np.abs(KY) <= kmax) &
               (np.abs(KZ) <= kmax)).astype(np.float64)

    K_mag = np.sqrt(K2)

    return {
        'X': X, 'Y': Y, 'Z': Z,
        'KX': KX, 'KY': KY, 'KZ': KZ,
        'K2': K2, 'K2_safe': K2_safe,
        'K_mag': K_mag, 'dealias': dealias,
        'N': N, 'L': L,
    }


# ============================================================
# INITIAL CONDITIONS
# ============================================================

def taylor_green_ic(grid):
    """Taylor-Green vortex: analytic, divergence-free, well-studied benchmark."""
    X, Y, Z = grid['X'], grid['Y'], grid['Z']
    u0 = np.sin(X) * np.cos(Y) * np.cos(Z)
    v0 = -np.cos(X) * np.sin(Y) * np.cos(Z)
    w0 = np.zeros_like(X)
    return fftn(u0), fftn(v0), fftn(w0)


# ============================================================
# NAVIER-STOKES RHS (Rotational / Lamb form)
# ============================================================

def compute_rhs(uh, vh, wh, grid, nu):
    """
    RHS of incompressible NS using rotational (Lamb) form with Leray projection.

    Rotational form avoids computing the pressure gradient explicitly:
        du/dt = -(ω × u) - ∇p + ν∇²u
    Leray projection removes the pressure term by enforcing ∇·u = 0.
    """
    KX, KY, KZ = grid['KX'], grid['KY'], grid['KZ']
    K2, K2_safe = grid['K2'], grid['K2_safe']
    dealias = grid['dealias']

    uhd = uh * dealias
    vhd = vh * dealias
    whd = wh * dealias

    # Vorticity in Fourier space
    ox_h = 1j * (KY * whd - KZ * vhd)
    oy_h = 1j * (KZ * uhd - KX * whd)
    oz_h = 1j * (KX * vhd - KY * uhd)

    # To physical space
    u = np.real(ifftn(uhd))
    v = np.real(ifftn(vhd))
    w = np.real(ifftn(whd))
    ox = np.real(ifftn(ox_h))
    oy = np.real(ifftn(oy_h))
    oz = np.real(ifftn(oz_h))

    # Nonlinear: -(ω × u)
    nlx_h = fftn(-(oy * w - oz * v)) * dealias
    nly_h = fftn(-(oz * u - ox * w)) * dealias
    nlz_h = fftn(-(ox * v - oy * u)) * dealias

    # Leray projection: remove divergent part
    kdot = (KX * nlx_h + KY * nly_h + KZ * nlz_h) / K2_safe
    kdot[0, 0, 0] = 0.0
    nlx_h -= KX * kdot
    nly_h -= KY * kdot
    nlz_h -= KZ * kdot

    # Add viscous diffusion
    return (nlx_h - nu * K2 * uh,
            nly_h - nu * K2 * vh,
            nlz_h - nu * K2 * wh)


# ============================================================
# DIAGNOSTICS
# ============================================================

def kinetic_energy(uh, vh, wh, N):
    """Volume-averaged kinetic energy."""
    return 0.5 * np.sum(np.abs(uh)**2 + np.abs(vh)**2 + np.abs(wh)**2).real / N**6


def enstrophy(uh, vh, wh, grid):
    """Volume-averaged enstrophy = mean |ω|²."""
    KX, KY, KZ = grid['KX'], grid['KY'], grid['KZ']
    ox = np.real(ifftn(1j * (KY * wh - KZ * vh)))
    oy = np.real(ifftn(1j * (KZ * uh - KX * wh)))
    oz = np.real(ifftn(1j * (KX * vh - KY * uh)))
    return np.mean(ox**2 + oy**2 + oz**2)


def energy_spectrum(uh, vh, wh, grid):
    """Spherically averaged energy spectrum E(k)."""
    N = grid['N']
    K_mag = grid['K_mag']
    Ek_3d = 0.5 * (np.abs(uh)**2 + np.abs(vh)**2 + np.abs(wh)**2) / N**6
    kmax_spec = N // 2
    k_vals = np.arange(1, kmax_spec + 1, dtype=float)
    E_k = np.zeros(kmax_spec)
    for i in range(kmax_spec):
        shell = (K_mag >= i + 0.5) & (K_mag < i + 1.5)
        E_k[i] = np.sum(Ek_3d[shell])
    return k_vals, E_k


# ============================================================
# STRAIN DIAGNOSTICS — Papers 165-166
# ============================================================

def velocity_gradient_tensor(uh, vh, wh, grid):
    """
    Full velocity gradient tensor A_ij = ∂u_i/∂x_j at every grid point.
    Returns shape (3, 3, N, N, N).
    """
    N = grid['N']
    KX, KY, KZ = grid['KX'], grid['KY'], grid['KZ']
    dealias = grid['dealias']

    A = np.zeros((3, 3, N, N, N))
    vel_hats = [uh * dealias, vh * dealias, wh * dealias]
    K_comps = [KX, KY, KZ]
    for i in range(3):
        for j in range(3):
            A[i, j] = np.real(ifftn(1j * K_comps[j] * vel_hats[i]))
    return A


def compute_QR_invariants(A):
    """
    Q-R invariants of the velocity gradient tensor and strain-rate tensor.

    Full gradient: Q = -½ tr(A²),  R = -⅓ tr(A³)
    Strain only:   Q_S = -½ tr(S²), R_S = -⅓ tr(S³)  where S = (A+Aᵀ)/2

    The (Q_S, R_S) joint PDF is the key diagnostic for Paper 165's
    Steinbach line test.
    """
    A2 = np.einsum('ijxyz,jkxyz->ikxyz', A, A)
    A3 = np.einsum('ijxyz,jkxyz->ikxyz', A, A2)
    Q = -0.5 * (A2[0, 0] + A2[1, 1] + A2[2, 2])
    R = -(1.0 / 3.0) * (A3[0, 0] + A3[1, 1] + A3[2, 2])

    # Strain-rate tensor
    S = np.zeros_like(A)
    for i in range(3):
        for j in range(3):
            S[i, j] = 0.5 * (A[i, j] + A[j, i])

    S2 = np.einsum('ijxyz,jkxyz->ikxyz', S, S)
    S3 = np.einsum('ijxyz,jkxyz->ikxyz', S, S2)
    Q_S = -0.5 * (S2[0, 0] + S2[1, 1] + S2[2, 2])
    R_S = -(1.0 / 3.0) * (S3[0, 0] + S3[1, 1] + S3[2, 2])

    return Q, R, Q_S, R_S


def strain_eigenvalues(Q_S, R_S):
    """
    Analytic eigenvalues of the traceless symmetric strain-rate tensor.

    The characteristic polynomial is λ³ + Qλ - R = 0 (since tr(S) = 0).
    Solved via Cardano / trigonometric method.

    Returns (λ₁, λ₂, λ₃) sorted λ₁ ≥ λ₂ ≥ λ₃.
    """
    Q_flat = Q_S.flatten()
    R_flat = R_S.flatten()
    neg_Q = np.maximum(-Q_flat, 1e-30)
    m = 2.0 * np.sqrt(neg_Q / 3.0)

    arg = np.zeros_like(Q_flat)
    valid = neg_Q > 1e-30
    arg[valid] = (3.0 * R_flat[valid]) / (Q_flat[valid] * m[valid])
    arg = np.clip(arg, -1.0, 1.0)
    theta = (1.0 / 3.0) * np.arccos(arg)

    lam1 = m * np.cos(theta)
    lam2 = m * np.cos(theta - 2.0 * np.pi / 3.0)
    lam3 = m * np.cos(theta - 4.0 * np.pi / 3.0)

    eigs = np.stack([lam1, lam2, lam3], axis=0)
    eigs = np.sort(eigs, axis=0)[::-1]
    return eigs[0], eigs[1], eigs[2]


def steinbach_roots():
    """
    Roots of the Steinbach polynomial: t³ - 7t/3 + 7/27 = 0

    This cubic arises from the icosahedral geometry of S³/2I.
    Its coefficients contain only primes {3, 7}.
    """
    return np.sort(np.roots([1, 0, -7.0 / 3.0, 7.0 / 27.0]))[::-1]


# ============================================================
# TIME INTEGRATION: Classical RK4
# ============================================================

def rk4_step(uh, vh, wh, grid, nu, dt):
    """Single RK4 time step."""
    k1u, k1v, k1w = compute_rhs(uh, vh, wh, grid, nu)
    k2u, k2v, k2w = compute_rhs(
        uh + 0.5 * dt * k1u, vh + 0.5 * dt * k1v, wh + 0.5 * dt * k1w,
        grid, nu)
    k3u, k3v, k3w = compute_rhs(
        uh + 0.5 * dt * k2u, vh + 0.5 * dt * k2v, wh + 0.5 * dt * k2w,
        grid, nu)
    k4u, k4v, k4w = compute_rhs(
        uh + dt * k3u, vh + dt * k3v, wh + dt * k3w,
        grid, nu)

    uh += (dt / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
    vh += (dt / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
    wh += (dt / 6.0) * (k1w + 2 * k2w + 2 * k3w + k4w)
    return uh, vh, wh


# ============================================================
# MAIN
# ============================================================

def main():
    args = parse_args()

    N = args.N
    Re = args.Re
    nu = 1.0 / Re
    dt = args.dt
    T_final = args.T
    print_every = args.print_every
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    # Banner
    print("Bootstrap DNS Solver")
    print("=" * 60)
    print(f"Grid:        {N}³ = {N**3:,} points")
    print(f"Domain:      [0, 2π]³ periodic")
    print(f"Viscosity:   ν = {nu:.6f}  (Re = {Re:.0f})")
    print(f"Time step:   dt = {dt}")
    print(f"Final time:  T = {T_final}")
    print(f"IC:          Taylor-Green vortex")
    print(f"Output:      {outdir}/")
    print("=" * 60)

    grid = setup_grid(N)
    uh, vh, wh = taylor_green_ic(grid)

    n_steps = int(round(T_final / dt))

    # Storage
    times_hist = []
    energy_hist = []
    enstrophy_hist = []
    spectra_saved = {}
    field_at_peak = None
    peak_enstrophy = 0.0
    peak_time = 0.0

    print("\nStarting time integration...\n")
    t_wall_start = time()

    for step in range(1, n_steps + 1):
        uh, vh, wh = rk4_step(uh, vh, wh, grid, nu, dt)
        t = step * dt

        if step % print_every == 0 or step == n_steps:
            E = kinetic_energy(uh, vh, wh, N)
            W = enstrophy(uh, vh, wh, grid)
            times_hist.append(t)
            energy_hist.append(E)
            enstrophy_hist.append(W)

            if W > peak_enstrophy:
                peak_enstrophy = W
                peak_time = t
                field_at_peak = (uh.copy(), vh.copy(), wh.copy())

            elapsed = time() - t_wall_start
            rate = step / elapsed if elapsed > 0 else 0
            print(f"  step {step:5d}  t = {t:6.3f}  "
                  f"E = {E:.6f}  Ω = {W:.4f}  [{rate:.0f} steps/s]")

        # Save spectrum at integer times
        if abs(t - round(t)) < dt / 2 and round(t) > 0:
            k_s, E_s = energy_spectrum(uh, vh, wh, grid)
            spectra_saved[f"t={round(t):.0f}"] = (k_s, E_s)

    elapsed_total = time() - t_wall_start
    print(f"\n{'=' * 60}")
    print(f"Simulation complete: {elapsed_total:.1f}s  ({n_steps} steps)")
    print(f"Peak enstrophy: Ω = {peak_enstrophy:.4f} at t = {peak_time:.2f}")
    print(f"Final energy:   E = {energy_hist[-1]:.6f}")

    # Final spectrum
    k_s, E_s = energy_spectrum(uh, vh, wh, grid)
    spectra_saved["t=final"] = (k_s, E_s)

    # ========================================================
    # STRAIN DIAGNOSTICS
    # ========================================================
    print(f"\n{'=' * 60}")
    print(f"Computing strain diagnostics at peak dissipation (t={peak_time:.2f})...")

    uh_p, vh_p, wh_p = field_at_peak if field_at_peak else (uh, vh, wh)

    A = velocity_gradient_tensor(uh_p, vh_p, wh_p, grid)
    Q, R, Q_S, R_S = compute_QR_invariants(A)
    lam1, lam2, lam3 = strain_eigenvalues(Q_S, R_S)

    sroots = steinbach_roots()
    steinbach_r21 = sroots[1] / sroots[0]
    steinbach_r31 = sroots[2] / sroots[0]

    valid = np.abs(lam1) > 1e-10
    ratio_21 = lam2[valid] / lam1[valid]
    ratio_31 = lam3[valid] / lam1[valid]

    print(f"  λ₁ (extensional) mean: {np.mean(lam1):.6f}")
    print(f"  λ₂ (intermediate) mean: {np.mean(lam2):.6f}")
    print(f"  λ₃ (compressive)  mean: {np.mean(lam3):.6f}")
    print(f"  Incompressibility:      {np.mean(lam1 + lam2 + lam3):.2e}")
    print(f"")
    print(f"  Steinbach roots: τ = [{sroots[0]:.4f}, {sroots[1]:.4f}, {sroots[2]:.4f}]")
    print(f"  Steinbach λ₂/λ₁ = {steinbach_r21:.4f}   DNS median = {np.median(ratio_21):.4f}")
    print(f"  Steinbach λ₃/λ₁ = {steinbach_r31:.4f}   DNS median = {np.median(ratio_31):.4f}")

    # ========================================================
    # SAVE DATA
    # ========================================================
    print(f"\n{'=' * 60}")
    print("Saving data...")

    # Diagnostics time series
    with open(os.path.join(outdir, "diagnostics.json"), 'w') as f:
        json.dump({
            'parameters': {
                'N': N, 'Re': Re, 'nu': nu,
                'dt': dt, 'T_final': T_final,
            },
            'time': times_hist,
            'energy': energy_hist,
            'enstrophy': enstrophy_hist,
            'peak_time': peak_time,
            'peak_enstrophy': peak_enstrophy,
            'wall_time_seconds': elapsed_total,
        }, f, indent=2)

    # Spectra
    spec_data = {}
    for label, (k_s, E_s) in spectra_saved.items():
        spec_data[f'k_{label}'] = k_s
        spec_data[f'E_{label}'] = E_s
    np.savez(os.path.join(outdir, "spectrum.npz"), **spec_data)

    # Strain data
    np.savez(os.path.join(outdir, "strain_data.npz"),
             Q=Q.astype(np.float32),
             R=R.astype(np.float32),
             Q_S=Q_S.astype(np.float32),
             R_S=R_S.astype(np.float32),
             lam1=lam1.astype(np.float32),
             lam2=lam2.astype(np.float32),
             lam3=lam3.astype(np.float32))

    # Final state (for restart or further analysis)
    np.savez(os.path.join(outdir, "dns_state.npz"),
             uh=uh, vh=vh, wh=wh)

    print(f"  Saved to {outdir}/")

    # ========================================================
    # PLOTS
    # ========================================================
    if not args.no_plots:
        try:
            from plot_diagnostics import generate_all_plots
            generate_all_plots(outdir)
        except ImportError:
            # Plot module not on path — try relative import
            try:
                sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
                from plot_diagnostics import generate_all_plots
                generate_all_plots(outdir)
            except ImportError:
                print("  (Run plot_diagnostics.py separately to generate figures)")

    # Summary
    print(f"\n{'=' * 60}")
    print("RESULTS SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Resolution:       {N}³")
    print(f"  Reynolds number:  Re = {Re:.0f}")
    print(f"  Peak dissipation: t = {peak_time:.2f}")
    print(f"  Energy decay:     {energy_hist[0]:.4f} → {energy_hist[-1]:.4f}")
    print(f"  Peak enstrophy:   {peak_enstrophy:.4f}")
    print(f"  Wall time:        {elapsed_total:.1f}s")
    print(f"")
    print(f"  STEINBACH TEST (Paper 165):")
    print(f"    λ₂/λ₁  predicted: {steinbach_r21:.4f}  measured: {np.median(ratio_21):.4f}")
    print(f"    λ₃/λ₁  predicted: {steinbach_r31:.4f}  measured: {np.median(ratio_31):.4f}")
    print(f"{'=' * 60}")
    print(f"\n🐕☕⬡ DNS complete.\n")


if __name__ == "__main__":
    main()
