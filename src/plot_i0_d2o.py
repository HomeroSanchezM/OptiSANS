#!/usr/bin/env python3
"""
Extract I(0) from Pepsi-SANS .dat files and plot vs D2O fraction.

Usage:
    python plot_i0_d2o.py <folder>
    python plot_i0_d2o.py <folder> <output_prefix>

The folder should contain .out files named like:
    d2o0.dat , etc.
The D2O percentage is parsed from the filename (first occurrence of d2oXX).
"""

import sys
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def parse_d2o_percent(filename):
    """Extract D2O percentage from filename like d2o_10_deut0_d2o10.out"""
    match = re.search(r'd2o(\d+)', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return None


def extract_I0(filepath):
    """Extract I(0): first value of second column (first non-comment data line)."""
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) >= 2:
                return float(cols[1])
    raise ValueError(f"No data found in {filepath}")


def main(folder, output_prefix=None):
    # Collect all .dat files
    files = [
        os.path.join(folder, f)
        for f in os.listdir(folder)
        if f.endswith('.dat')
    ]

    if not files:
        print(f"No .dat files found in {folder}")
        sys.exit(1)

    data = []
    for filepath in files:
        d2o = parse_d2o_percent(filepath)
        if d2o is None:
            print(f"  Warning: could not parse D2O% from {os.path.basename(filepath)}, skipping.")
            continue
        try:
            I0 = extract_I0(filepath)
            data.append((d2o, I0))
            print(f"  {os.path.basename(filepath):40s}  D2O={d2o:3d}%  I(0)={I0:.6e}")
        except Exception as e:
            print(f"  Error reading {filepath}: {e}")

    if not data:
        print("No valid data found.")
        sys.exit(1)

    data.sort(key=lambda x: x[0])
    d2o_vals = np.array([d[0] for d in data])
    I0_vals  = np.array([d[1] for d in data])

    if output_prefix is None:
        output_prefix = os.path.join(folder, "I0_vs_d2o")

    # ── Figure 1 : linear scale ──────────────────────────────────────────────
    I0_abs = np.abs(I0_vals)
    zero_mask = I0_abs == 0
    
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(d2o_vals, I0_vals, 'o-', color='steelblue', linewidth=1.8,
            markersize=6, markerfacecolor='white', markeredgewidth=1.8)
    ax.set_xlabel("D₂O fraction (%)", fontsize=13)
    ax.set_ylabel("I(0)  (cm⁻¹)", fontsize=13)
    ax.set_title("I(0) vs D₂O linear scale half deuterated --d2o variation", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.axhline(0, color='black', linewidth=0.8, linestyle=':')
    #ax.axhline(0.015, color='darkorange', linewidth=1.5, linestyle='--', label="I(0) limit = 0.015")
    
     # Also annotate the minimum (match point) even if not exactly zero
    idx_min = np.argmin(I0_abs)
    ax.axvline(d2o_vals[idx_min], color='gray', linestyle='--', linewidth=1,
               label=f"Min |I(0)| at {d2o_vals[idx_min]}% D₂O = {I0_vals[idx_min]}")
    
    ax.legend(fontsize=10)
    fig.tight_layout()
    path_linear = output_prefix + "_linear.png"
    fig.savefig(path_linear, dpi=150)
    print(f"\nSaved: {path_linear}")

    # ── Figure 2 : log scale of |I(0)| ──────────────────────────────────────
    I0_abs = np.abs(I0_vals)
    zero_mask = I0_abs == 0

    fig, ax = plt.subplots(figsize=(7, 5))

    # Plot nonzero points normally
    if np.any(~zero_mask):
        ax.semilogy(d2o_vals[~zero_mask], I0_abs[~zero_mask],
                    'o-', color='firebrick', linewidth=1.8,
                    markersize=6, markerfacecolor='white', markeredgewidth=1.8,
                    label="|I(0)|")

    # Mark match point (I0 ≈ 0 or minimum) if any zeros
    if np.any(zero_mask):
        y_min = I0_abs[~zero_mask].min() if np.any(~zero_mask) else 1e-10
        ax.semilogy(d2o_vals[zero_mask], np.full(zero_mask.sum(), y_min * 0.5),
                    'v', color='orange', markersize=10, label="Match point (I=0)")

    # Also annotate the minimum (match point) even if not exactly zero
    idx_min = np.argmin(I0_abs)
    ax.axvline(d2o_vals[idx_min], color='gray', linestyle='--', linewidth=1,
               label=f"Min |I(0)| at {d2o_vals[idx_min]}% D₂O")

    #ax.axhline(0.015, color='darkorange', linewidth=1.5, linestyle='--', label="I(0) limit = 0.015")
    ax.set_xlabel("D₂O fraction (%)", fontsize=13)
    ax.set_ylabel("|I(0)|  (cm⁻¹)", fontsize=13)
    ax.set_title("I(0) vs D₂O — log|I(0)|", fontsize=14)
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=10)
    fig.tight_layout()
    path_log = output_prefix + "_log.png"
    fig.savefig(path_log, dpi=150)
    print(f"Saved: {path_log}")

    # ── Figure 3 : both panels side by side ─────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(d2o_vals, I0_vals, 'o-', color='steelblue', linewidth=1.8,
             markersize=6, markerfacecolor='white', markeredgewidth=1.8)
    ax1.set_xlabel("D₂O (%)", fontsize=12)
    ax1.set_ylabel("I(0)  (cm⁻¹)", fontsize=12)
    ax1.set_title("Linear scale", fontsize=13)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.axhline(0, color='black', linewidth=0.8, linestyle=':')
    #ax1.axhline(0.015, color='darkorange', linewidth=1.5, linestyle='--', label="I(0) limit = 0.015")
    
     # Also annotate the minimum (match point) even if not exactly zero
    idx_min = np.argmin(I0_abs)
    ax1.axvline(d2o_vals[idx_min], color='gray', linestyle='--', linewidth=1,
               label=f"Min |I(0)| at {d2o_vals[idx_min]}% D₂O = {I0_vals[idx_min]}")
    
    
    ax1.legend(fontsize=10)

    if np.any(~zero_mask):
        ax2.semilogy(d2o_vals[~zero_mask], I0_abs[~zero_mask],
                     'o-', color='firebrick', linewidth=1.8,
                     markersize=6, markerfacecolor='white', markeredgewidth=1.8)
    if np.any(zero_mask):
        ax2.semilogy(d2o_vals[zero_mask],
                     np.full(zero_mask.sum(), I0_abs[~zero_mask].min() * 0.5),
                     'v', color='orange', markersize=10)
    ax2.axvline(d2o_vals[idx_min], color='gray', linestyle='--', linewidth=1,
                label=f"Min at {d2o_vals[idx_min]}% D₂O")
    ax2.axhline(0.015, color='darkorange', linewidth=1.5, linestyle='--', label="I(0) limit = 0.015")
    ax2.set_xlabel("D₂O (%)", fontsize=12)
    ax2.set_ylabel("|I(0)|  (cm⁻¹)", fontsize=12)
    ax2.set_title("Log |I(0)|", fontsize=13)
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)
    ax2.legend(fontsize=10)

    fig.suptitle("Contrast variation — I(0) vs D₂O fraction", fontsize=14, y=1.01)
    fig.tight_layout()
    path_both = output_prefix + "_combined.png"
    fig.savefig(path_both, dpi=150, bbox_inches='tight')
    print(f"Saved: {path_both}")
    plt.close('all')


def plot_i0_vs_d2o(folder: str, output_prefix: str) -> None:
    """
    Programmatic interface: generate I(0) vs D2O% plots.

    Args:
        folder:        Directory containing the .dat SANS files.
        output_prefix: Full path prefix for output files (no extension).
                       E.g. "/path/to/output/I0_vs_d2o" will produce
                       "/path/to/output/I0_vs_d2o_linear.png", etc.
    """
    main(folder, output_prefix)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_i0_d2o.py <folder> [output_prefix]")
        sys.exit(1)
    _out_prefix = sys.argv[2] if len(sys.argv) >= 3 else None
    main(sys.argv[1], _out_prefix)
