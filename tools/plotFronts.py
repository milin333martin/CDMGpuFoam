#!/usr/bin/env python3
"""
plotFronts.py <case> [<case> ...]

Plots the shock and flame tracking written by CDMGpuFoam to

    <case>/postProcessing/fronts/<startTime>/fronts.dat

Three panels: the leading shock and flame positions against time, their
speeds, and the shock-flame separation (the decoupling signal).  With
--channels, a fourth panel shows how far the shock got into each channel of
an arrestor.

    ./plotFronts.py br09cht
    ./plotFronts.py br05 br09cht --channels          # compare two cases
    ./plotFronts.py br09cht --save fronts.png

Columns of fronts.dat: time, then five per group (shock x, shock speed,
flame x, flame speed, separation).  Group 0 is the whole domain, then one
group per channel.
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read(case):
    pat = os.path.join(case, "postProcessing", "fronts", "*", "fronts.dat")
    files = sorted(glob.glob(pat))
    if not files:
        # a run restarted from a later time writes a second directory
        files = sorted(glob.glob(os.path.join(case, "**", "fronts.dat"),
                                 recursive=True))
    if not files:
        raise SystemExit("no fronts.dat under " + case)

    names, data = [], []
    for f in files:
        head = [l for l in open(f) if l.startswith("#")]
        for l in head:
            if "columns" in l:
                cols = l.split(":", 1)[1].split()
                names = [c.split("_", 1)[1] for c in cols
                         if c.startswith("xShock_")]
        # read by hand: a run that was interrupted can leave a partial line,
        # and a case rebuilt with a different channel count can leave rows of
        # another width in the same file
        rows, widths = [], {}
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            widths[len(v)] = widths.get(len(v), 0) + 1
            rows.append(v)
        if not rows:
            continue
        width = max(widths, key=widths.get)     # the most common width wins
        bad = sum(n for w, n in widths.items() if w != width)
        if bad:
            print("%s: skipped %d malformed row(s) of %d" % (f, bad, len(rows)))
        good = [r for r in rows if len(r) == width]
        data.append(np.array(good, dtype=float))

    if not data:
        raise SystemExit("no usable rows in " + case)
    width = min(a.shape[1] for a in data)
    d = np.vstack([a[:, :width] for a in data])
    d = d[np.argsort(d[:, 0])]
    n = (d.shape[1] - 1)//5
    if not names:
        names = ["all"] + ["ch%d" % i for i in range(1, n)]
    return d, names


def group(d, i):
    """shock x, shock speed, flame x, flame speed, separation of group i"""
    c = 1 + 5*i
    return d[:, c], d[:, c+1], d[:, c+2], d[:, c+3], d[:, c+4]


def main():
    a = argparse.ArgumentParser()
    a.add_argument("cases", nargs="+")
    a.add_argument("--channels", action="store_true",
                   help="add a panel of the penetration into each channel")
    a.add_argument("--arrestor", type=float, nargs=2, default=None,
                   metavar=("X0", "X1"),
                   help="arrestor position in m, marked on the plots")
    a.add_argument("--save", default=None)
    args = a.parse_args()

    nPanels = 4 if args.channels else 3
    fig, ax = plt.subplots(nPanels, 1, figsize=(9, 3*nPanels), sharex=True)

    for case in args.cases:
        d, names = read(case)
        t = d[:, 0]*1e6                       # us
        xs, Ds, xf, Df, sep = group(d, 0)
        label = os.path.basename(os.path.normpath(case))

        l, = ax[0].plot(t, xs*1000, label=label + " shock")
        ax[0].plot(t, xf*1000, "--", color=l.get_color(),
                   label=label + " flame")
        ax[1].plot(t, Ds, color=l.get_color(), label=label)
        ax[2].plot(t, sep*1000, color=l.get_color(), label=label)

        if args.channels and len(names) > 1:
            x0 = args.arrestor[0] if args.arrestor else None
            for i in range(1, len(names)):
                xsi = group(d, i)[0]
                ax[3].plot(t, xsi*1000, lw=0.8,
                           label=names[i] if len(args.cases) == 1 else None)

        # the arrestor, if given
        if args.arrestor:
            for p in (0, 3) if args.channels else (0,):
                ax[p].axhspan(args.arrestor[0]*1000, args.arrestor[1]*1000,
                              color="0.85", zorder=0)

    ax[0].set_ylabel("position [mm]")
    ax[0].legend(fontsize=8, ncol=2)
    ax[1].set_ylabel("shock speed [m/s]")
    ax[1].axhline(1970.06, color="k", lw=0.6, ls=":")
    ax[1].text(ax[1].get_xlim()[0], 1970.06, " D$_{CJ}$", va="bottom",
               fontsize=8)
    ax[2].set_ylabel("shock-flame\nseparation [mm]")
    ax[2].set_yscale("symlog", linthresh=0.1)
    if args.channels:
        ax[3].set_ylabel("per-channel\nshock position [mm]")
        if len(args.cases) == 1:
            ax[3].legend(fontsize=6, ncol=4)
    ax[-1].set_xlabel("time [$\\mu$s]")

    for x in ax:
        x.grid(alpha=0.3)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150)
        print("written", args.save)
    else:
        out = "fronts.png"
        fig.savefig(out, dpi=150)
        print("written", out)
        if os.environ.get("DISPLAY"):
            plt.show()


if __name__ == "__main__":
    main()
