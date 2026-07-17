#!/usr/bin/env python3
"""Build distance-based phylogenetic trees from a CSV distance matrix.

Features:
- Accepts full or triangular (upper/lower) distance matrices with first column as row labels
- Builds NJ and UPGMA trees using scikit-bio
- Optionally builds minimum-evolution style tree if `gme` or `bme` is available in scikit-bio
- Roots trees by outgroup (if provided)
- Exports Newick
- Renders PNGs in rectangular / radial / unrooted / diagonal layouts
"""

from __future__ import annotations

import argparse
import io
import math
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np
import pandas as pd

try:  # pragma: no cover
    from skbio import DistanceMatrix
    from skbio.tree import nj, upgma
except Exception:  # pragma: no cover
    DistanceMatrix = None  # type: ignore
    nj = None  # type: ignore
    upgma = None  # type: ignore

try:  # pragma: no cover
    from skbio.tree import gme  # type: ignore
except Exception:  # pragma: no cover
    gme = None

try:  # pragma: no cover
    from skbio.tree import bme  # type: ignore
except Exception:  # pragma: no cover
    bme = None

try:  # pragma: no cover
    from ete3 import NodeStyle, TextFace, Tree, TreeStyle
except Exception:  # pragma: no cover
    NodeStyle = None  # type: ignore
    TextFace = None  # type: ignore
    Tree = Any  # type: ignore
    TreeStyle = None  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build NJ/UPGMA/ME trees from a CSV distance matrix and render multiple layouts."
    )
    parser.add_argument("--input", required=True, help="Path to CSV distance matrix")
    parser.add_argument(
        "--outdir",
        default="tree_outputs",
        help="Output directory (default: tree_outputs)",
    )
    parser.add_argument(
        "--prefix", default="dist_tree", help="Output filename prefix (default: dist_tree)"
    )
    parser.add_argument(
        "--methods",
        default="nj,upgma",
        help="Comma-separated: nj,upgma,me (default: nj,upgma)",
    )
    parser.add_argument(
        "--me-algorithm",
        default="auto",
        choices=["auto", "gme", "bme"],
        help="Minimum-evolution algorithm if --methods includes me (default: auto)",
    )
    parser.add_argument(
        "--layouts",
        default="rectangular,radial,unrooted,diagonal",
        help=(
            "Comma-separated render layouts: rectangular,radial,unrooted,diagonal "
            "(default: rectangular,radial,unrooted,diagonal)"
        ),
    )
    parser.add_argument("--outgroup", default=None, help="Outgroup taxon name for rooting")
    parser.add_argument(
        "--font-size", type=int, default=10, help="Leaf label font size (default: 10)"
    )
    parser.add_argument(
        "--label-color", default="black", help="Default leaf label color (default: black)"
    )
    parser.add_argument(
        "--label-colors-csv",
        default=None,
        help="Optional CSV mapping labels to colors with columns: label,color",
    )
    parser.add_argument(
        "--branch-color", default="#444444", help="Branch line color (default: #444444)"
    )
    parser.add_argument(
        "--show-branch-lengths",
        action="store_true",
        help="Show branch lengths in ETE-rendered layouts",
    )
    parser.add_argument(
        "--show-scale", action="store_true", help="Show scale bar in ETE-rendered layouts"
    )
    parser.add_argument("--width", type=int, default=1800, help="Output image width in px")
    parser.add_argument("--height", type=int, default=1200, help="Output image height in px")
    parser.add_argument("--dpi", type=int, default=200, help="DPI for diagonal layout image")
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-12,
        help="Tolerance for symmetry checks (default: 1e-12)",
    )
    parser.add_argument(
        "--triangle-policy",
        choices=["mean", "upper", "lower"],
        default="mean",
        help="How to reconcile upper/lower entries when both present and differ (default: mean)",
    )
    return parser.parse_args()


def _require_skbio() -> None:
    if DistanceMatrix is None or nj is None or upgma is None:
        raise RuntimeError(
            "scikit-bio is required for tree building. Install with: pip install scikit-bio"
        )


def _require_ete3() -> None:
    if NodeStyle is None or TextFace is None or TreeStyle is None:
        raise RuntimeError("ete3 is required for rendering. Install with: pip install ete3")


def _coerce_numeric_matrix(df: pd.DataFrame) -> np.ndarray:
    return df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)


def _fill_triangular_matrix(mat: np.ndarray, policy: str, tol: float) -> np.ndarray:
    n = mat.shape[0]
    out = mat.copy()

    for i in range(n):
        out[i, i] = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            a = out[i, j]
            b = out[j, i]

            a_nan = np.isnan(a)
            b_nan = np.isnan(b)

            if a_nan and b_nan:
                raise ValueError(
                    f"Missing both mirrored entries at ({i},{j}) and ({j},{i}); "
                    "cannot infer full distance matrix."
                )

            if a_nan:
                out[i, j] = b
                continue
            if b_nan:
                out[j, i] = a
                continue

            if math.isclose(a, b, abs_tol=tol):
                continue

            if policy == "upper":
                out[j, i] = a
            elif policy == "lower":
                out[i, j] = b
            else:  # mean
                m = (a + b) / 2.0
                out[i, j] = m
                out[j, i] = m

    return out


def load_distance_matrix(csv_path: Path, triangle_policy: str, tol: float) -> DistanceMatrix:
    raw = pd.read_csv(csv_path)

    if raw.shape[1] < 2:
        raise ValueError("CSV must include one label column and >=1 distance column")

    row_labels = raw.iloc[:, 0].astype(str).str.strip().tolist()
    col_labels = pd.Index(raw.columns[1:]).astype(str).str.strip().tolist()

    if len(row_labels) != len(col_labels):
        raise ValueError(
            f"Matrix dimensions mismatch: {len(row_labels)} rows vs {len(col_labels)} columns"
        )

    if row_labels != col_labels:
        raise ValueError(
            "Row labels (column 1) and column headers must match in the same order"
        )

    numeric = _coerce_numeric_matrix(raw.iloc[:, 1:])

    if numeric.shape[0] != numeric.shape[1]:
        raise ValueError("Distance matrix must be square")

    numeric = _fill_triangular_matrix(numeric, policy=triangle_policy, tol=tol)

    if not np.allclose(np.diag(numeric), 0.0, atol=tol):
        raise ValueError("Matrix diagonal must be zeros")

    if not np.allclose(numeric, numeric.T, atol=tol):
        raise ValueError("Matrix is not symmetric after reconciliation")

    return DistanceMatrix(numeric, row_labels)


def _tree_to_newick(tree_obj) -> str:
    sio = io.StringIO()
    tree_obj.write(sio)
    return sio.getvalue().strip()


def build_trees(dm: DistanceMatrix, methods: Sequence[str], me_algorithm: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    chosen = {m.strip().lower() for m in methods if m.strip()}

    if "nj" in chosen:
        out["NJ"] = _tree_to_newick(nj(dm))

    if "upgma" in chosen:
        out["UPGMA"] = _tree_to_newick(upgma(dm))

    if "me" in chosen:
        me_builder = None
        me_name = None

        if me_algorithm == "gme":
            me_builder = gme
            me_name = "GME"
        elif me_algorithm == "bme":
            me_builder = bme
            me_name = "BME"
        else:
            if bme is not None:
                me_builder = bme
                me_name = "BME"
            elif gme is not None:
                me_builder = gme
                me_name = "GME"

        if me_builder is None:
            raise RuntimeError(
                "Requested ME tree, but no ME builder found in scikit-bio. "
                "Install a scikit-bio version exposing bme/gme or skip me."
            )

        out[me_name] = _tree_to_newick(me_builder(dm))

    if not out:
        raise ValueError("No valid methods requested. Use any of: nj, upgma, me")

    return out


def load_label_color_map(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    df = pd.read_csv(path)
    if not {"label", "color"}.issubset(df.columns):
        raise ValueError("Label colors CSV must contain columns: label,color")
    return {str(r["label"]): str(r["color"]) for _, r in df.iterrows()}


def root_and_style_tree(
    newick: str,
    outgroup: Optional[str],
    font_size: int,
    default_label_color: str,
    label_colors: Dict[str, str],
    branch_color: str,
) -> Tree:
    tree = Tree(newick, format=1)

    leaves = {leaf.name for leaf in tree.iter_leaves()}
    if outgroup:
        if outgroup in leaves:
            tree.set_outgroup(tree & outgroup)
        else:
            print(f"WARNING: outgroup '{outgroup}' not found; keeping original rooting")

    for node in tree.traverse():
        style = NodeStyle()
        style["size"] = 0
        style["hz_line_width"] = 2
        style["vt_line_width"] = 2
        style["hz_line_color"] = branch_color
        style["vt_line_color"] = branch_color
        node.set_style(style)

    for leaf in tree.iter_leaves():
        color = label_colors.get(leaf.name, default_label_color)
        face = TextFace(leaf.name, fsize=font_size, fgcolor=color)
        leaf.add_face(face, column=0, position="aligned")

    return tree


def render_ete_layout(
    tree: Tree,
    layout: str,
    output_png: Path,
    width: int,
    height: int,
    show_branch_lengths: bool,
    show_scale: bool,
) -> None:
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = show_branch_lengths
    ts.show_scale = show_scale

    if layout == "rectangular":
        ts.mode = "r"
    elif layout == "radial":
        ts.mode = "c"
        ts.arc_start = -180
        ts.arc_span = 360
    elif layout == "unrooted":
        ts.mode = "c"
        ts.arc_start = 0
        ts.arc_span = 360
        ts.root_opening_factor = 0.0
    else:
        raise ValueError(f"Unsupported ETE layout: {layout}")

    tree.render(str(output_png), w=width, h=height, tree_style=ts)


def _distance_from_root(tree: Tree) -> Dict[Tree, float]:
    coords: Dict[Tree, float] = {tree: 0.0}
    for node in tree.traverse("preorder"):
        parent_x = coords[node]
        for child in node.children:
            dist = child.dist if child.dist is not None else 0.0
            coords[child] = parent_x + max(0.0, float(dist))
    if max(coords.values(), default=0.0) == 0.0:
        coords = {node: float(i) for i, node in enumerate(tree.traverse("preorder"))}
    return coords


def _leaf_y_positions(tree: Tree) -> Dict[Tree, float]:
    leaves = list(tree.iter_leaves())
    y: Dict[Tree, float] = {leaf: float(i) for i, leaf in enumerate(leaves)}

    def assign(node: Tree) -> float:
        if node.is_leaf():
            return y[node]
        vals = [assign(child) for child in node.children]
        y[node] = sum(vals) / len(vals)
        return y[node]

    assign(tree)
    return y


def render_diagonal_layout(
    tree: Tree,
    output_png: Path,
    width: int,
    height: int,
    dpi: int,
    font_size: int,
    default_label_color: str,
    label_colors: Dict[str, str],
    branch_color: str,
    show_branch_lengths: bool,
) -> None:
    import matplotlib.pyplot as plt

    x = _distance_from_root(tree)
    y = _leaf_y_positions(tree)

    fig_w = max(4.0, width / dpi)
    fig_h = max(3.0, height / dpi)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    for node in tree.traverse("preorder"):
        for child in node.children:
            ax.plot(
                [x[node], x[child]],
                [y[node], y[child]],
                color=branch_color,
                linewidth=1.5,
                zorder=1,
            )

            if show_branch_lengths:
                mx = (x[node] + x[child]) / 2.0
                my = (y[node] + y[child]) / 2.0
                bl = child.dist if child.dist is not None else 0.0
                ax.text(mx, my, f"{bl:.4g}", fontsize=max(6, font_size - 2), color="#555555")

    x_span = max(x.values()) - min(x.values()) if x else 1.0
    label_offset = max(0.02 * (x_span if x_span > 0 else 1.0), 0.02)

    for leaf in tree.iter_leaves():
        color = label_colors.get(leaf.name, default_label_color)
        ax.text(
            x[leaf] + label_offset,
            y[leaf],
            leaf.name,
            va="center",
            ha="left",
            fontsize=font_size,
            color=color,
            zorder=2,
        )

    ax.set_axis_off()
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(output_png, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    _require_skbio()
    _require_ete3()

    csv_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    methods = [m.strip().lower() for m in args.methods.split(",") if m.strip()]
    layouts = [l.strip().lower() for l in args.layouts.split(",") if l.strip()]

    dm = load_distance_matrix(
        csv_path=csv_path, triangle_policy=args.triangle_policy, tol=args.tol
    )
    trees = build_trees(dm, methods=methods, me_algorithm=args.me_algorithm)

    label_colors = load_label_color_map(Path(args.label_colors_csv) if args.label_colors_csv else None)

    for method_name, newick in trees.items():
        newick_path = outdir / f"{args.prefix}.{method_name}.newick"
        newick_path.write_text(newick + "\n", encoding="utf-8")

        styled = root_and_style_tree(
            newick=newick,
            outgroup=args.outgroup,
            font_size=args.font_size,
            default_label_color=args.label_color,
            label_colors=label_colors,
            branch_color=args.branch_color,
        )

        for layout in layouts:
            png_path = outdir / f"{args.prefix}.{method_name}.{layout}.png"

            if layout in {"rectangular", "radial", "unrooted"}:
                render_ete_layout(
                    tree=styled.copy(method="deepcopy"),
                    layout=layout,
                    output_png=png_path,
                    width=args.width,
                    height=args.height,
                    show_branch_lengths=args.show_branch_lengths,
                    show_scale=args.show_scale,
                )
            elif layout == "diagonal":
                render_diagonal_layout(
                    tree=styled.copy(method="deepcopy"),
                    output_png=png_path,
                    width=args.width,
                    height=args.height,
                    dpi=args.dpi,
                    font_size=args.font_size,
                    default_label_color=args.label_color,
                    label_colors=label_colors,
                    branch_color=args.branch_color,
                    show_branch_lengths=args.show_branch_lengths,
                )
            else:
                raise ValueError(
                    f"Unknown layout '{layout}'. Use rectangular, radial, unrooted, diagonal"
                )

    print(f"Done. Wrote outputs to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
