#!/usr/bin/env python3
"""Streamlit UI for distance-matrix tree inference (NJ/UPGMA)."""

from __future__ import annotations

import io

import streamlit as st

from tree_core import (
    enforce_ultrametric,
    handle_negative_branch_lengths,
    infer_distance_tree,
    parse_distance_csv,
    plot_rectangular_tree,
    root_tree,
    tree_to_newick,
)


st.set_page_config(page_title="Distance Matrix to Tree", layout="wide")
st.markdown(
    """
    <div style="
        background: linear-gradient(90deg, #0f766e 0%, #115e59 100%);
        padding: 18px 22px;
        border-radius: 12px;
        color: white;
        margin-bottom: 14px;
    ">
      <h1 style="margin:0; font-size: 1.9rem;">Distance Matrix to Tree</h1>
      <p style="margin:8px 0 0 0; font-size:1rem;">
        Upload a pairwise distance matrix and infer NJ or UPGMA trees with validation, outgroup rooting, and styling controls.
      </p>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.expander("Input format instructions", expanded=True):
    st.markdown(
        """
Upload a CSV distance matrix with:
- First column = sample names
- First row headers = the same sample names, same order
- Diagonal values = `0`
- You may provide a full matrix or only upper/lower triangle (leave opposite triangle blank)
        """
    )

    example_csv = (
        "Taxon,Sp1,Sp2,Sp3,Sp4\n"
        "Sp1,0,,,\n"
        "Sp2,0.12,0,,\n"
        "Sp3,0.18,0.09,0,\n"
        "Sp4,0.25,0.14,0.07,0\n"
    )
    st.download_button(
        label="Download example CSV",
        data=example_csv,
        file_name="example_distance_matrix.csv",
        mime="text/csv",
    )

uploaded = st.file_uploader("Upload distance matrix CSV", type=["csv"])

if uploaded is None:
    st.info("Upload a CSV to continue.")
    st.stop()

st.subheader("Tree Inference Settings")

left, right = st.columns(2)
with left:
    algorithm = st.radio("Algorithm", options=["nj", "upgma"], horizontal=True)
    branch_mode = st.radio(
        "Tree style",
        options=["phylogram (use branch lengths)", "cladogram (equal branch steps)"],
    )
    ultrametric = st.toggle(
        "Ultrametric enforcement",
        value=False,
        help="Forces equal root-to-tip distances by extending terminal branches.",
    )

with right:
    negative_mode = st.selectbox(
        "Negative branch lengths",
        options=["keep", "clamp", "redistribute"],
        help="Some distance methods can produce negatives. Clamp sets negatives to 0.",
    )
    triangle_policy = st.selectbox(
        "Triangular conflict policy",
        options=["mean", "upper", "lower"],
        help="If both mirrored cells are present but differ, choose how to reconcile them.",
    )
    check_triangle = st.toggle("Check triangle inequality", value=False)

st.subheader("Display Settings")
col1, col2, col3, col4 = st.columns(4)
with col1:
    font_size = st.slider("Font size", min_value=7, max_value=24, value=11)
with col2:
    font_style = st.selectbox("Font style", options=["normal", "italic", "oblique"])
with col3:
    label_color = st.color_picker("Label color", value="#111111")
with col4:
    branch_color = st.color_picker("Branch color", value="#3a3a3a")

show_branch_lengths = st.toggle("Show branch lengths", value=True)
line_width = st.slider("Branch line width", min_value=0.5, max_value=4.0, value=1.6, step=0.1)

st.subheader("Validation")

try:
    dm, mat, report = parse_distance_csv(
        uploaded,
        triangle_policy=triangle_policy,
        check_triangle_inequality=check_triangle,
    )
except Exception as exc:
    st.error(f"Matrix parsing/validation failed: {exc}")
    st.stop()

st.success(f"Parsed {report.n_taxa} samples.")
st.caption(f"Samples: {', '.join(report.labels)}")

if check_triangle:
    if report.triangle_violations == 0:
        st.success("No triangle inequality violations detected.")
    else:
        st.warning(f"Triangle inequality violations: {report.triangle_violations}")
        with st.expander("Show first violations"):
            for row in report.triangle_examples:
                st.write(f"- {row}")

outgroup_choices = ["(none)"] + report.labels
outgroup = st.selectbox("Outgroup", options=outgroup_choices)
outgroup_value = None if outgroup == "(none)" else outgroup

if st.button("Build Tree", type="primary"):
    try:
        tree = infer_distance_tree(dm, algorithm=algorithm)
        tree = root_tree(tree, outgroup_value)
        tree = handle_negative_branch_lengths(tree, mode=negative_mode)

        if ultrametric:
            tree = enforce_ultrametric(tree)

        use_branch_lengths = branch_mode.startswith("phylogram")
        fig = plot_rectangular_tree(
            tree=tree,
            use_branch_lengths=use_branch_lengths,
            show_branch_lengths=show_branch_lengths,
            font_size=font_size,
            font_style=font_style,
            label_color=label_color,
            branch_color=branch_color,
            line_width=line_width,
        )
        st.pyplot(fig, clear_figure=True)

        newick = tree_to_newick(tree)
        st.download_button(
            label="Download Newick",
            data=(newick + "\n"),
            file_name=f"distance_tree_{algorithm}.newick",
            mime="text/plain",
        )

        png_bytes = io.BytesIO()
        fig.savefig(png_bytes, format="png", dpi=220, bbox_inches="tight")
        st.download_button(
            label="Download PNG",
            data=png_bytes.getvalue(),
            file_name=f"distance_tree_{algorithm}.png",
            mime="image/png",
        )
    except Exception as exc:
        st.error(f"Tree inference failed: {exc}")
