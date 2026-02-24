#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px


def generate_data(seed: int, n: int):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "genotype": [f"G{i + 1:03d}" for i in range(n)],
            "yield_t_ha": rng.normal(6.4, 0.75, n),
            "protein_pct": rng.normal(11.9, 1.1, n),
            "disease_score": rng.normal(4.2, 0.9, n),
        }
    )


def main():
    st.set_page_config(page_title="Breeding Dashboard Demo", layout="wide")
    st.title("Breeding Dashboard Demo (Streamlit)")

    n = st.sidebar.slider("Number of genotypes", 40, 300, 120)
    seed = st.sidebar.number_input("Random seed", value=42)
    metric = st.sidebar.selectbox(
        "Color metric", ["yield_t_ha", "protein_pct", "disease_score"]
    )

    df = generate_data(int(seed), int(n))
    st.dataframe(df.head(20), use_container_width=True)

    fig = px.scatter(
        df,
        x="yield_t_ha",
        y="protein_pct",
        color=metric,
        hover_data=["genotype"],
        title="Yield vs Protein with Metric Overlay",
    )
    st.plotly_chart(fig, use_container_width=True)

    st.markdown(
        "This minimal app demonstrates an interactive pattern that can be wired to qtl-analysis or breeding-trial-management outputs."
    )


if __name__ == "__main__":
    main()
