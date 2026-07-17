#!/usr/bin/env python3
"""Core utilities for distance-matrix tree inference and plotting."""

from __future__ import annotations

import io
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import TreeNode, nj, upgma


@dataclass
class ValidationReport:
    n_taxa: int
    labels: List[str]
    triangle_violations: int
    triangle_examples: List[str]


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
                    f"Missing both mirrored entries at ({i},{j}) and ({j},{i})."
                )
            if a_nan:
                out[i, j] = b
                continue
            if b_nan:
                out[j, i] = a
                continue

            if np.isclose(a, b, atol=tol):
                continue

            if policy == "upper":
                out[j, i] = a
            elif policy == "lower":
                out[i, j] = b
            else:
                m = (a + b) / 2.0
                out[i, j] = m
                out[j, i] = m

    return out


def triangle_inequality_violations(
    mat: np.ndarray, labels: Sequence[str], tol: float, max_examples: int = 8
) -> Tuple[int, List[str]]:
    n = mat.shape[0]
    violations = 0
    examples: List[str] = []

    for i in range(n):
        for j in range(i + 1, n):
            dij = mat[i, j]
            for k in range(n):
                if k == i or k == j:
                    continue
                rhs = mat[i, k] + mat[k, j]
                if dij > rhs + tol:
                    violations += 1
                    if len(examples) < max_examples:
                        examples.append(
                            f"d({labels[i]}, {labels[j]})={dij:.6g} > "
                            f"d({labels[i]}, {labels[k]}) + d({labels[k]}, {labels[j]})={rhs:.6g}"
                        )
                    break

    return violations, examples


def parse_distance_csv(
    csv_data,
    triangle_policy: str = "mean",
    tol: float = 1e-12,
    check_triangle_inequality: bool = False,
) -> Tuple[DistanceMatrix, np.ndarray, ValidationReport]:
    raw = pd.read_csv(csv_data)

    if raw.shape[1] < 2:
        raise ValueError("CSV must contain one label column plus distance columns.")

    labels = raw.iloc[:, 0].astype(str).str.strip().tolist()
    col_labels = pd.Index(raw.columns[1:]).astype(str).str.strip().tolist()

    if len(labels) != len(col_labels):
        raise ValueError(
            f"Matrix dimensions mismatch: {len(labels)} rows vs {len(col_labels)} columns."
        )

    if labels != col_labels:
        raise ValueError(
            "Row labels (first column) and column headers must match in the same order."
        )

    if len(set(labels)) != len(labels):
        raise ValueError("Sample labels must be unique.")

    numeric = _coerce_numeric_matrix(raw.iloc[:, 1:])

    if np.isnan(numeric).all():
        raise ValueError("No numeric distances found in matrix body.")

    if numeric.shape[0] != numeric.shape[1]:
        raise ValueError("Distance matrix must be square.")

    numeric = _fill_triangular_matrix(numeric, policy=triangle_policy, tol=tol)

    if np.isnan(numeric).any():
        raise ValueError("Distance matrix still contains NaN after triangular fill.")

    if not np.allclose(np.diag(numeric), 0.0, atol=tol):
        raise ValueError("Matrix diagonal must be zero.")

    if (numeric < -tol).any():
        raise ValueError("Distance matrix contains negative distances.")

    if not np.allclose(numeric, numeric.T, atol=tol):
        raise ValueError("Distance matrix is not symmetric.")

    if check_triangle_inequality:
        violations, examples = triangle_inequality_violations(numeric, labels, tol=tol)
    else:
        violations, examples = 0, []

    dm = DistanceMatrix(numeric, labels)
    report = ValidationReport(
        n_taxa=len(labels),
        labels=labels,
        triangle_violations=violations,
        triangle_examples=examples,
    )
    return dm, numeric, report


def infer_distance_tree(dm: DistanceMatrix, algorithm: str) -> TreeNode:
    algo = algorithm.lower().strip()
    if algo == "nj":
        return nj(dm)
    if algo == "upgma":
        return upgma(dm)
    raise ValueError("Unsupported algorithm. Choose 'nj' or 'upgma'.")


def root_tree(tree: TreeNode, outgroup: Optional[str]) -> TreeNode:
    if not outgroup:
        return tree

    tip_names = {tip.name for tip in tree.tips()}
    if outgroup not in tip_names:
        raise ValueError(f"Outgroup '{outgroup}' not found among tips.")

    rooted = None
    try:
        rooted = tree.root_by_outgroup(outgroup)
    except Exception:
        rooted = tree.root_by_outgroup([outgroup])
    if rooted is None:
        return tree
    return rooted


def _iter_nonroot_nodes(tree: TreeNode):
    for node in tree.preorder():
        if node.parent is None:
            continue
        yield node


def handle_negative_branch_lengths(tree: TreeNode, mode: str) -> TreeNode:
    mode = mode.lower().strip()
    if mode == "keep":
        return tree

    for node in _iter_nonroot_nodes(tree):
        length = 0.0 if node.length is None else float(node.length)
        if length >= 0:
            continue

        deficit = -length
        node.length = 0.0

        if mode == "redistribute" and len(node.children) > 0:
            share = deficit / len(node.children)
            for child in node.children:
                current = 0.0 if child.length is None else float(child.length)
                child.length = max(0.0, current + share)

    if mode not in {"clamp", "redistribute"}:
        raise ValueError("Invalid negative branch mode. Use keep, clamp, or redistribute.")

    return tree


def enforce_ultrametric(tree: TreeNode, tol: float = 1e-12) -> TreeNode:
    root_to_tip: Dict[str, float] = {}

    def walk(node: TreeNode, dist: float) -> None:
        if node.is_tip():
            root_to_tip[node.name] = dist
            return
        for child in node.children:
            bl = 0.0 if child.length is None else max(0.0, float(child.length))
            walk(child, dist + bl)

    walk(tree, 0.0)

    if not root_to_tip:
        return tree

    target = max(root_to_tip.values())

    for tip in tree.tips():
        current = root_to_tip.get(tip.name, 0.0)
        delta = target - current
        if delta > tol:
            old = 0.0 if tip.length is None else float(tip.length)
            tip.length = max(0.0, old + delta)

    return tree


def tree_to_newick(tree: TreeNode) -> str:
    sio = io.StringIO()
    tree.write(sio)
    return sio.getvalue().strip()


def _node_positions(
    tree: TreeNode,
    use_branch_lengths: bool,
) -> Tuple[Dict[TreeNode, float], Dict[TreeNode, float], List[TreeNode]]:
    x: Dict[TreeNode, float] = {tree: 0.0}

    for node in tree.preorder():
        parent_x = x[node]
        for child in node.children:
            if use_branch_lengths:
                step = 0.0 if child.length is None else max(0.0, float(child.length))
            else:
                step = 1.0
            x[child] = parent_x + step

    tips = list(tree.tips())
    y: Dict[TreeNode, float] = {tip: float(i) for i, tip in enumerate(tips)}

    def assign(node: TreeNode) -> float:
        if node.is_tip():
            return y[node]
        vals = [assign(ch) for ch in node.children]
        y[node] = sum(vals) / len(vals)
        return y[node]

    assign(tree)

    if max(x.values(), default=0.0) == 0.0:
        # Fallback so labels do not collapse if all branch lengths are zero
        for node in tree.preorder():
            if node.parent is None:
                continue
            x[node] = x[node.parent] + 1.0

    return x, y, tips


def plot_rectangular_tree(
    tree: TreeNode,
    use_branch_lengths: bool,
    show_branch_lengths: bool,
    font_size: int,
    font_style: str,
    label_color: str,
    branch_color: str,
    line_width: float,
    figsize: Tuple[float, float] = (11.0, 7.0),
):
    x, y, tips = _node_positions(tree, use_branch_lengths=use_branch_lengths)

    fig, ax = plt.subplots(figsize=figsize)

    for node in tree.preorder():
        if len(node.children) == 0:
            continue

        child_ys = [y[ch] for ch in node.children]
        ax.plot(
            [x[node], x[node]],
            [min(child_ys), max(child_ys)],
            color=branch_color,
            linewidth=line_width,
            zorder=1,
        )

        for child in node.children:
            ax.plot(
                [x[node], x[child]],
                [y[child], y[child]],
                color=branch_color,
                linewidth=line_width,
                zorder=1,
            )

            if use_branch_lengths and show_branch_lengths:
                bl = 0.0 if child.length is None else float(child.length)
                mx = (x[node] + x[child]) / 2.0
                ax.text(
                    mx,
                    y[child] + 0.18,
                    f"{bl:.4g}",
                    fontsize=max(6, font_size - 2),
                    color="#555555",
                    ha="center",
                    va="bottom",
                )

    x_span = max(x.values(), default=1.0) - min(x.values(), default=0.0)
    x_offset = max(0.01, x_span * 0.03)

    for tip in tips:
        ax.text(
            x[tip] + x_offset,
            y[tip],
            str(tip.name),
            fontsize=font_size,
            color=label_color,
            fontstyle=font_style,
            va="center",
            ha="left",
            zorder=2,
        )

    ax.set_axis_off()
    ax.invert_yaxis()
    fig.tight_layout()
    return fig
