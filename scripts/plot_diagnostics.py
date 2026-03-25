#!/usr/bin/env python3
"""
Bootstrap DNS — Diagnostic Plots
=================================
Generates publication-quality figures from DNS output data.

Usage:
    python plot_diagnostics.py results/
    python plot_diagnostics.py results/ --dpi 300

Dr. Clifford Keeble & Claude
Bootstrap Universe Programme, March 2026
"""

import argparse
import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# ============================================================
# STEINBACH REFERENCE
# ============================================================

def steinbach_roots():
    """Roots of t³ - 7t/3 + 7/27 = 0"""
    return np.sort(np.roots([1, 0, -7.0 / 3.0, 7.0 / 27.0]))[::-1]


# ============================================================
# FIGURE 1: Energy and Enstrophy
# ============================================================

def plot_energy_enstrophy(outdir, dpi=150):
    with open(os.path.join(outdir, "diagnostics.json")) as f:
        diag = json.load(f)

    params = diag['parameters']
    Re = params['Re']
    N = params['N']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(diag['time'], diag['energy'], 'b-', linewidth=1.5)
    ax1.set_xlabel('Time', fontsize=12)
    ax1.set_ylabel('Kinetic Energy E(t)', fontsize=12)
    ax1.set_title('Energy Decay', fontsize=14)
    ax1.axhline(y=0.125, color='gray', linestyle='--', alpha=0.5, label='E₀ = 0.125')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    ax2.plot(diag['time'], diag['enstrophy'], 'r-', linewidth=1.5)
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Enstrophy Ω(t)', fontsize=12)
    ax2.set_title('Enstrophy Evolution', fontsize=14)
    ax2.axvline(x=diag['peak_time'], color='gray', linestyle='--', alpha=0.5,
                label=f'Peak t = {diag["peak_time"]:.1f}')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(f'Bootstrap DNS — Taylor-Green {N}³, Re = {Re:.0f}',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig1_energy.png')
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  [1/4] {path}")


# ============================================================
# FIGURE 2: Energy Spectrum
# ============================================================

def plot_spectrum(outdir, dpi=150):
    data = np.load(os.path.join(outdir, "spectrum.npz"))
    keys = sorted([k for k in data.keys() if k.startswith('k_')])

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(keys)))

    for key, c in zip(keys, colors):
        label = key.replace('k_', '')
        k_s = data[key]
        E_s = data[key.replace('k_', 'E_')]
        mask = E_s > 0
        ax.loglog(k_s[mask], E_s[mask], '-', color=c, linewidth=1.5,
                  label=label, alpha=0.8)

    # Reference slope
    k_ref = np.array([2, 30])
    E_ref = 0.005 * k_ref**(-5.0 / 3.0)
    ax.loglog(k_ref, E_ref, 'k--', linewidth=2, alpha=0.6, label='k⁻⁵/³')

    ax.set_xlabel('Wavenumber k', fontsize=12)
    ax.set_ylabel('E(k)', fontsize=12)
    ax.set_title('Energy Spectrum — Paper 166 Reference', fontsize=14)
    ax.legend(fontsize=10, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    path = os.path.join(outdir, 'fig2_spectrum.png')
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  [2/4] {path}")


# ============================================================
# FIGURE 3: Q-R Invariants — THE KEY DIAGNOSTIC
# ============================================================

def plot_QR_invariants(outdir, dpi=150):
    data = np.load(os.path.join(outdir, "strain_data.npz"))
    Q = data['Q'].flatten()
    R = data['R'].flatten()
    Q_S = data['Q_S'].flatten()
    R_S = data['R_S'].flatten()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # --- Left: Full gradient Q-R ---
    Q_rms = max(np.sqrt(np.mean(Q**2)), 1e-20)
    R_rms = max(np.sqrt(np.mean(R**2)), 1e-20)

    h1, xe, ye = np.histogram2d(
        R / R_rms, Q / Q_rms, bins=200, range=[[-5, 5], [-5, 5]])
    h1 = np.maximum(h1, 0.5)
    ax1.pcolormesh(xe, ye, h1.T, norm=LogNorm(vmin=1, vmax=h1.max()),
                   cmap='inferno')

    # Vieillefosse discriminant
    R_line = np.linspace(-5, 5, 500)
    disc_Q = -(np.abs(R_line * R_rms)**2 * 27.0 / 4.0)**(1.0 / 3.0) / Q_rms
    ax1.plot(R_line, disc_Q, 'c-', linewidth=1.5, alpha=0.8,
             label='Vieillefosse D = 0')

    ax1.set_xlabel('R / R_rms', fontsize=12)
    ax1.set_ylabel('Q / Q_rms', fontsize=12)
    ax1.set_title('Q-R Joint PDF (Full Gradient)', fontsize=14)
    ax1.set_xlim(-5, 5)
    ax1.set_ylim(-5, 5)
    ax1.set_aspect('equal')
    ax1.legend(fontsize=10)

    # --- Right: Strain Q_S-R_S with Steinbach line ---
    Q_S_rms = max(np.sqrt(np.mean(Q_S**2)), 1e-20)
    R_S_rms = max(np.sqrt(np.mean(R_S**2)), 1e-20)

    h2, xe2, ye2 = np.histogram2d(
        R_S / R_S_rms, Q_S / Q_S_rms, bins=200, range=[[-5, 5], [-5, 5]])
    h2 = np.maximum(h2, 0.5)
    ax2.pcolormesh(xe2, ye2, h2.T, norm=LogNorm(vmin=1, vmax=h2.max()),
                   cmap='inferno')

    # Discriminant
    R_S_line = np.linspace(-5, 5, 500)
    disc_QS = -(np.abs(R_S_line * R_S_rms)**2 * 27.0 / 4.0)**(1.0 / 3.0) / Q_S_rms
    ax2.plot(R_S_line, disc_QS, 'c-', linewidth=1.5, alpha=0.8, label='Δ = 0')

    # Steinbach constraint line
    # Steinbach: Q_S = -7/3, R_S = 7/27 → R/(-Q)^{3/2} = const
    stein_q = -7.0 / 3.0
    stein_r = 7.0 / 27.0
    stein_ratio = stein_r / (-stein_q)**1.5
    Q_sl = np.linspace(-4.5, -0.1, 300)
    R_sl = stein_ratio * (-Q_sl * Q_S_rms)**1.5 / R_S_rms
    ax2.plot(R_sl, Q_sl, 'g--', linewidth=2, alpha=0.8, label='Steinbach line')

    ax2.set_xlabel('R_S / R_S,rms', fontsize=12)
    ax2.set_ylabel('Q_S / Q_S,rms', fontsize=12)
    ax2.set_title('Q_S–R_S (Strain) — Paper 165 Test', fontsize=14)
    ax2.set_xlim(-5, 5)
    ax2.set_ylim(-5, 5)
    ax2.set_aspect('equal')
    ax2.legend(fontsize=10)

    fig.suptitle('Velocity Gradient Invariants at Peak Dissipation',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig3_QR.png')
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  [3/4] {path}")


# ============================================================
# FIGURE 4: Eigenvalue Ratio Distributions
# ============================================================

def plot_eigenvalue_ratios(outdir, dpi=150):
    data = np.load(os.path.join(outdir, "strain_data.npz"))
    lam1 = data['lam1'].flatten()
    lam2 = data['lam2'].flatten()
    lam3 = data['lam3'].flatten()

    sroots = steinbach_roots()
    steinbach_r21 = sroots[1] / sroots[0]
    steinbach_r31 = sroots[2] / sroots[0]

    valid = np.abs(lam1) > 1e-10
    ratio_21 = lam2[valid] / lam1[valid]
    ratio_31 = lam3[valid] / lam1[valid]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.hist(ratio_21, bins=200, range=(-1.5, 1.5), density=True,
             color='steelblue', alpha=0.7, label='DNS')
    ax1.axvline(x=steinbach_r21, color='red', linewidth=2, linestyle='--',
                label=f'Steinbach: {steinbach_r21:.3f}')
    ax1.set_xlabel('λ₂ / λ₁', fontsize=12)
    ax1.set_ylabel('PDF', fontsize=12)
    ax1.set_title('Intermediate / Extensional', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    ax2.hist(ratio_31, bins=200, range=(-3, 1), density=True,
             color='indianred', alpha=0.7, label='DNS')
    ax2.axvline(x=steinbach_r31, color='red', linewidth=2, linestyle='--',
                label=f'Steinbach: {steinbach_r31:.3f}')
    ax2.set_xlabel('λ₃ / λ₁', fontsize=12)
    ax2.set_ylabel('PDF', fontsize=12)
    ax2.set_title('Compressive / Extensional', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.suptitle('Strain Eigenvalue Ratios — Steinbach Test (Paper 165)',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig4_eigenvalues.png')
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  [4/4] {path}")


# ============================================================
# ENTRY POINT
# ============================================================

def generate_all_plots(outdir, dpi=150):
    """Generate all diagnostic figures from saved data."""
    print(f"\nGenerating diagnostic plots from {outdir}/...")
    plot_energy_enstrophy(outdir, dpi)
    plot_spectrum(outdir, dpi)
    plot_QR_invariants(outdir, dpi)
    plot_eigenvalue_ratios(outdir, dpi)
    print("  Done.\n")


def main():
    p = argparse.ArgumentParser(description="Generate DNS diagnostic plots")
    p.add_argument("outdir", help="Directory containing DNS output data")
    p.add_argument("--dpi", type=int, default=150, help="Plot DPI (default: 150)")
    args = p.parse_args()
    generate_all_plots(args.outdir, args.dpi)


if __name__ == "__main__":
    main()
