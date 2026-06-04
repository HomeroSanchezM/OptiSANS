#!/usr/bin/env python3
"""
Plot D2O % vs Fitness Score from a fitness result CSV.

Usage:
    python plot_fitness.py results.csv
    python plot_fitness.py results.csv --output plot.png
    python plot_fitness.py results.csv --no-labels
    python plot_fitness.py results.csv --label-step 5
"""

import argparse
import re
import sys
import csv
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


# ============================================================================
#                           ARGUMENT PARSING
# ============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot D2O%% vs fitness score from a fitness result CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
CSV format expected:
  filename,fitness_score,ratio
  _d2o42.dat,0.00126857,0.13437720
  ...

The D2O percentage is extracted from the filename using the pattern d2oXX.

Examples:
  python plot_fitness.py results.csv
  python plot_fitness.py results.csv --output fitness_plot.png
  python plot_fitness.py results.csv --no-labels
  python plot_fitness.py results.csv --label-step 10
        """
    )

    parser.add_argument(
        'csv_file',
        type=str,
        help='Path to the fitness result CSV file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        metavar='FILE',
        help='Save the plot to a file instead of displaying it (e.g. plot.png, plot.pdf)'
    )
    parser.add_argument(
        '--no-labels',
        action='store_true',
        help='Disable per-point annotations (fitness + ratio labels)'
    )
    parser.add_argument(
        '--label-step',
        type=int,
        default=1,
        metavar='N',
        help='Show a label every N points to reduce clutter (default: 1 = all points)'
    )
    parser.add_argument(
        '--title',
        type=str,
        default=None,
        help='Custom plot title (default: derived from CSV filename)'
    )

    return parser.parse_args()


# ============================================================================
#                           DATA LOADING
# ============================================================================

def extract_d2o(filename):
    """Extract the D2O percentage integer from a filename like _d2o42.dat."""
    match = re.search(r'd2o(\d+)', filename, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def load_csv(path):
    """
    Load the CSV and return sorted lists of (d2o_percent, fitness, ratio).
    Rows where D2O cannot be parsed are skipped with a warning.
    """
    records = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            d2o = extract_d2o(row['filename'])
            if d2o is None:
                print(f"Warning: cannot extract D2O%% from '{row['filename']}', skipping.",
                      file=sys.stderr)
                continue
            try:
                fitness = float(row['fitness_score'])
                ratio   = float(row['ratio'])
            except (KeyError, ValueError) as e:
                print(f"Warning: bad numeric data in row {row}: {e}", file=sys.stderr)
                continue
            records.append((d2o, fitness, ratio))

    if not records:
        print("Error: no valid data rows found in CSV.", file=sys.stderr)
        sys.exit(1)

    records.sort(key=lambda r: r[0])
    d2o_vals   = [r[0] for r in records]
    fitness_vals = [r[1] for r in records]
    ratio_vals  = [r[2] for r in records]
    return d2o_vals, fitness_vals, ratio_vals


# ============================================================================
#                           PLOTTING
# ============================================================================

def plot(d2o_vals, fitness_vals, ratio_vals, args):
    fig, ax = plt.subplots(figsize=(14, 6))

    # ── colour-code: zero-fitness points (failed ratio check) in red ──
    colors = ['#d62728' if f == 0.0 else '#1f77b4' for f in fitness_vals]

    ax.scatter(d2o_vals, fitness_vals, c=colors, s=60, zorder=3, linewidths=0.5,
               edgecolors='white')
    ax.plot(d2o_vals, fitness_vals, color='#1f77b4', linewidth=1.2,
            alpha=0.5, zorder=2)

    # ── per-point annotations ──
    if not args.no_labels:
        for i, (x, y, r) in enumerate(zip(d2o_vals, fitness_vals, ratio_vals)):
            if i % args.label_step != 0:
                continue
            label = f"f={y:.4f}\nr={r:.4f}"
            txt = ax.annotate(
                label,
                xy=(x, y),
                xytext=(0, 10),
                textcoords='offset points',
                ha='center',
                va='bottom',
                fontsize=6.5,
                color='#222222',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.6,
                          ec='#aaaaaa', lw=0.5),
                arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.6)
            )
            txt.set_path_effects([pe.withStroke(linewidth=1, foreground='white')])

    # ── legend for failed points ──
    if any(f == 0.0 for f in fitness_vals):
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4',
                   markersize=8, label='Valid curve'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
                   markersize=8, label='Failed ratio check (fitness = 0)'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

    # ── axes & title ──
    ax.set_xlabel('D₂O (%)', fontsize=12)
    ax.set_ylabel('Fitness score', fontsize=12)

    title = args.title if args.title else f'Fitness vs D₂O% — {args.csv_file}'
    ax.set_title(title, fontsize=13, pad=12)

    ax.set_xlim(min(d2o_vals) - 1, max(d2o_vals) + 1)
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.5f'))
    ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f"Plot saved to '{args.output}'")
    else:
        plt.show()


def plot_fitness_from_csv(
    csv_file: str,
    output: str,
    no_labels: bool = False,
    label_step: int = 5,
    title: Optional[str] = None,
) -> None:
    """
    Programmatic interface: generate the fitness-vs-D2O% plot and save to `output`.

    Args:
        csv_file:   Path to the result CSV (columns: filename, fitness_score, ratio).
        output:     Output image path (e.g. "/path/to/fitness_vs_d2o.png").
        no_labels:  If True, suppress per-point annotations.
        label_step: Show a label every N points (default 5 to reduce clutter).
        title:      Custom plot title; if None, a default title is used.
    """
    import types
    args = types.SimpleNamespace(
        csv_file=csv_file,
        output=output,
        no_labels=no_labels,
        label_step=label_step,
        title=title if title else f"Fitness vs D₂O% — {Path(csv_file).name}",
    )
    d2o_vals, fitness_vals, ratio_vals = load_csv(csv_file)
    plot(d2o_vals, fitness_vals, ratio_vals, args)


# ============================================================================
#                           MAIN
# ============================================================================

def main():
    args = parse_arguments()
    d2o_vals, fitness_vals, ratio_vals = load_csv(args.csv_file)
    print(f"Loaded {len(d2o_vals)} data points (D2O range: {min(d2o_vals)}%–{max(d2o_vals)}%)")
    plot(d2o_vals, fitness_vals, ratio_vals, args)


if __name__ == '__main__':
    main()
