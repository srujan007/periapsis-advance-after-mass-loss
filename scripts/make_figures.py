#!/usr/bin/env python3
"""Write the two grayscale figures used in the companion note.

Usage
-----
    python scripts/make_figures.py
    python scripts/make_figures.py --outdir figures
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from massloss import AU, MSUN, bound_fraction
from massloss.precession import no_kick_advance_ratio, no_kick_century_rate_ratio

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.linewidth": 0.8,
        "legend.fontsize": 8.5,
        "mathtext.fontset": "stix",
    }
)


def fig_ratio(outdir: Path) -> None:
    mus = np.linspace(0.505, 0.995, 120)
    per_orbit = np.array([no_kick_advance_ratio(1.0, mu) for mu in mus])
    per_century = np.array([no_kick_century_rate_ratio(1.0, mu) for mu in mus])

    fig, ax = plt.subplots(figsize=(5.4, 3.6), dpi=180)
    ax.plot(mus, per_orbit, "k-", lw=1.8, label=r"Per orbit, $\Delta\phi_f/\Delta\phi_i=\mu^2$")
    ax.plot(mus, per_century, "k--", lw=1.8, label="Per Julian century")
    ax.axvline(0.5, color="0.5", ls=":", lw=1)
    ax.set_xlabel(r"Remnant mass fraction $\mu=M_f/M_i$")
    ax.set_ylabel("Residual / pre-collapse advance")
    ax.set_xlim(0.50, 1.00)
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    path = outdir / "fig1_precession_ratio.png"
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print("wrote", path)


def fig_bound_fraction(outdir: Path) -> None:
    mi = 15 * MSUN
    vks = np.linspace(0, 280, 29)
    styles = [
        (1, 8.0, "-", 1.8, "k", r"1 AU, 8 $M_\odot$"),
        (1, 6.0, "--", 1.6, "k", r"1 AU, 6 $M_\odot$"),
        (5, 8.0, "-", 1.8, "0.35", r"5 AU, 8 $M_\odot$"),
        (5, 6.0, "--", 1.6, "0.35", r"5 AU, 6 $M_\odot$"),
        (1, 1.4, ":", 1.8, "k", r"1 AU, 1.4 $M_\odot$"),
        (10, 8.0, "-.", 1.6, "0.15", r"10 AU, 8 $M_\odot$"),
    ]

    fig, ax = plt.subplots(figsize=(5.4, 3.6), dpi=180)
    for r_au, mf_msun, ls, lw, color, lab in styles:
        ys = [
            bound_fraction(mi, mf_msun * MSUN, r_au * AU, vk, n=12_000, seed=7)
            for vk in vks
        ]
        ax.plot(vks, ys, ls, color=color, lw=lw, label=lab)
    ax.set_xlabel(r"Kick speed $v_k$ (km s$^{-1}$)")
    ax.set_ylabel("Bound fraction")
    ax.set_xlim(0, 280)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, loc="upper right", handlelength=2.6)
    fig.tight_layout()
    path = outdir / "fig2_bound_fraction.png"
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print("wrote", path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig_ratio(outdir)
    fig_bound_fraction(outdir)


if __name__ == "__main__":
    main()
