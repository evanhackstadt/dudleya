"""
pcangsd_visualize.py -- plot PCAngsd covariance matrices as PC scatterplots.

Evan Hackstadt / Whittall Lab Dudleya Genomics (refactored to use the sample sheet)

Reads sample metadata from config/samples.tsv and styling from config/taxa.yaml;
it never parses sample names. The .cov rows are aligned to the sheet through
population.info (the sample order that ANGSD saw), so subset runs just work.

Outputs (in outdir)
  pca_species.png, pca_taxon.png, pca_population.png   full plots, three colorings
  pca_species_labeled.png                              species plot with per-sample labels
  zoom/pca_<variant>_zoom<k>.png                       auto-detected dense clusters (optional)
  points.csv, metrics.txt

Usable from Snakemake (script: directive) or standalone for post-run exploration:
  python pcangsd_visualize.py --cov X.cov --info population.info \\
      --sheet config/samples.tsv --taxa config/taxa.yaml --outdir out/
"""

import argparse
import itertools
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless HPC nodes
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.cluster.hierarchy import fcluster, linkage  # noqa: E402

# make `import sampleinfo` work however the script is launched
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sampleinfo import Taxonomy  # noqa: E402

VARIANTS = ("species", "taxon", "population")
DPI_FULL, DPI_ZOOM = 600, 300
# Okabe-Ito colors for subspecies lacking their own color (ordered to avoid near-twins of
# the species colors in taxa.yaml; orange is last because it resembles vermillion)
SUPPLEMENT = ["#CC79A7", "#F0E442", "#56B4E9", "#999999", "#000000", "#E69F00"]
POP_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


# --- styling ----------------------------------------------------------------

def make_styles(points, tax, variant):
    """Ordered {category: {color, marker, label}} for the categories of one coloring."""
    present = set(points[variant])
    if variant == "population":
        # >10 populations: reuse tab10 colors with a different marker every 10 populations.
        # Ordered by taxon then name, so related populations sit together in the legend.
        rank = {code: i for i, code in enumerate(tax.taxa)}
        owner = points.groupby("population")["taxon"].first()
        pops = sorted(present, key=lambda p: (rank[owner[p]], p))
        cmap = plt.get_cmap("tab10")
        return {
            p: dict(color=cmap(i % 10), marker=POP_MARKERS[(i // 10) % len(POP_MARKERS)], label=p)
            for i, p in enumerate(pops)
        }

    # species: color+marker from taxa.yaml.  taxon: shape from taxa.yaml (subspecies
    # inherit their species' unless overridden); color is the taxon's own, else a distinct one.
    spare = itertools.cycle(SUPPLEMENT)
    styles = {}
    for code in (c for c in tax.taxa if c in present):
        st = tax.style(code)
        own = variant == "species" or tax.own_color(code)
        styles[code] = dict(color=st["color"] if own else next(spare),
                            marker=st["marker"], label=tax.label(code))
    return styles


# --- zoom detection ---------------------------------------------------------

def dense_boxes(x, y, max_zooms=4, link_frac=0.05, min_points=5, max_span=0.3, pad=0.2):
    """Bounding boxes of tight clusters, as [((xmin, xmax), (ymin, ymax)), ...].

    Single-linkage clusters at a distance of `link_frac` of the plot extent; a cluster
    earns a zoom if it has >= min_points and spans < `max_span` of the extent (so zooming
    reveals something). Largest first. Works on the geometry only, not on any coloring.
    """
    xy = np.column_stack([x, y])
    if max_zooms < 1 or len(xy) < max(min_points, 2):
        return []
    span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1]))
    labels = fcluster(linkage(xy, "single"), t=link_frac * span, criterion="distance")
    found = []
    for lab in np.unique(labels):
        pts = xy[labels == lab]
        w, h = np.ptp(pts, axis=0)
        if len(pts) < min_points or max(w, h) > max_span * span:
            continue
        px, py = max(w * pad, 0.01 * span), max(h * pad, 0.01 * span)
        found.append((len(pts), (pts[:, 0].min() - px, pts[:, 0].max() + px),
                      (pts[:, 1].min() - py, pts[:, 1].max() + py)))
    found.sort(key=lambda f: -f[0])
    return [(xl, yl) for _, xl, yl in found[:max_zooms]]


# --- plotting ---------------------------------------------------------------

def draw(points, variant, styles, x_pc, y_pc, ve, path, dpi, xlim=None, ylim=None, annotate=False):
    xcol, ycol = f"PC{x_pc}", f"PC{y_pc}"
    view = points
    if xlim:
        view = points[points[xcol].between(*xlim) & points[ycol].between(*ylim)]

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for cat, st in styles.items():
        sub = view[view[variant] == cat]
        if sub.empty:
            continue
        ax.scatter(sub[xcol], sub[ycol], color=st["color"], marker=st["marker"], s=30,
                   edgecolors="black", linewidths=0.3, alpha=0.9, label=f"{st['label']} (n={len(sub)})")
    if annotate:
        for x, y, lab in zip(view[xcol], view[ycol], view["label"]):
            ax.annotate(lab, (x, y), fontsize=5)
    if xlim:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    ax.set_xlabel(f"{xcol} ({ve[x_pc - 1]:.2f}% of variance)")
    ax.set_ylabel(f"{ycol} ({ve[y_pc - 1]:.2f}% of variance)")
    ax.set_title(f"{len(view)}-sample PCA by {variant}" + (" (zoom)" if xlim else ""), fontsize=13)
    n_leg = len(ax.get_legend_handles_labels()[0])
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
              fontsize=9 if n_leg <= 15 else 7, ncol=1 if n_leg <= 25 else 2)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --- driver -----------------------------------------------------------------

def run(cov_path, info_path, sheet_path, taxa_path, outdir, max_zooms=4, x_pc=1, y_pc=2):
    tax = Taxonomy.load(taxa_path)
    samples = [line.strip() for line in open(info_path) if line.strip()]
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str, keep_default_na=False).set_index("sample", drop=False)
    missing = [s for s in samples if s not in sheet.index]
    if missing:
        raise ValueError(f"{len(missing)} sample(s) in {info_path} are not in {sheet_path}, e.g. {missing[:3]}")
    meta = sheet.loc[samples]

    cov = np.loadtxt(cov_path)
    if cov.shape != (len(samples), len(samples)):
        raise ValueError(f"{cov_path} is {cov.shape}, but {info_path} lists {len(samples)} samples")

    # PCs = eigenvectors of the covariance matrix, descending variance
    evals, evecs = np.linalg.eigh(cov)
    evals, evecs = evals[::-1], evecs[:, ::-1]
    ve = 100 * evals / evals.sum()
    k = min(10, len(samples))
    if max(x_pc, y_pc) > k:
        raise ValueError(f"requested PCs {x_pc},{y_pc} but only the first {k} are kept")

    pcs = pd.DataFrame(evecs[:, :k], index=meta.index, columns=[f"PC{i + 1}" for i in range(k)])
    points = pd.concat([meta[["taxon", "species", "population", "lp_id", "du_id", "batch"]], pcs], axis=1)
    points["label"] = points["du_id"].where(points["du_id"] != "", points.index.to_series())

    os.makedirs(os.path.join(outdir, "zoom"), exist_ok=True)
    boxes = dense_boxes(points[f"PC{x_pc}"], points[f"PC{y_pc}"], max_zooms)
    for i, (xl, yl) in enumerate(boxes, 1):
        print(f"zoom {i}: PC{x_pc} {xl[0]:.4f}..{xl[1]:.4f}, PC{y_pc} {yl[0]:.4f}..{yl[1]:.4f}")

    for variant in VARIANTS:
        styles = make_styles(points, tax, variant)
        args = (points, variant, styles, x_pc, y_pc, ve)
        draw(*args, os.path.join(outdir, f"pca_{variant}.png"), DPI_FULL)
        if variant == "species":
            draw(*args, os.path.join(outdir, "pca_species_labeled.png"), DPI_FULL, annotate=True)
        for i, (xl, yl) in enumerate(boxes, 1):
            draw(*args, os.path.join(outdir, "zoom", f"pca_{variant}_zoom{i}.png"), DPI_ZOOM,
                 xlim=xl, ylim=yl, annotate=True)

    points.drop(columns="label").to_csv(os.path.join(outdir, "points.csv"), index_label="sample")
    with open(os.path.join(outdir, "metrics.txt"), "w") as fh:
        fh.write(f"samples: {len(samples)}\n\nvariance explained (top 20 PCs)\n")
        fh.writelines(f"PC{i + 1}\t{v:.3f}%\n" for i, v in enumerate(ve[:20]))
        for col in ("species", "taxon", "batch"):
            fh.write(f"\nsamples by {col}\n")
            fh.writelines(f"{k_}\t{n}\n" for k_, n in points[col].value_counts().items())
    print(f"Plots, points.csv and metrics.txt saved to {outdir}")


if "snakemake" in globals():
    sm = snakemake  # noqa: F821  (injected by Snakemake's script: directive)
    run(sm.input.cov, sm.input.info, sm.input.sheet, sm.input.taxa,
        os.path.dirname(sm.output.species), getattr(sm.params, "max_zooms", 4))
elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cov", required=True)
    ap.add_argument("--info", required=True, help="population.info: sample order of the .cov")
    ap.add_argument("--sheet", default="config/samples.tsv")
    ap.add_argument("--taxa", default="config/taxa.yaml")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-zooms", type=int, default=4, help="0 disables auto-zoom")
    ap.add_argument("--pcs", type=int, nargs=2, default=[1, 2], metavar=("X", "Y"))
    a = ap.parse_args()
    run(a.cov, a.info, a.sheet, a.taxa, a.outdir, a.max_zooms, *a.pcs)