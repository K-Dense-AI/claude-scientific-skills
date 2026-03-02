#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from shiny import App, render, ui
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


app_ui = ui.page_fluid(
    ui.h2("Shiny Trait Dashboard Demo"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_slider("n", "Number of genotypes", 40, 260, 120),
            ui.input_numeric("seed", "Seed", 42),
        ),
        ui.output_plot("trait_plot"),
        ui.output_table("preview"),
    ),
)


def server(input, output, session):
    def make_data():
        rng = np.random.default_rng(int(input.seed()))
        n = int(input.n())
        return pd.DataFrame(
            {
                "yield_t_ha": rng.normal(6.3, 0.8, n),
                "protein_pct": rng.normal(12.0, 1.0, n),
                "disease_score": rng.normal(4.1, 0.9, n),
            }
        )

    @output
    @render.plot
    def trait_plot():
        df = make_data()
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(df["yield_t_ha"], df["protein_pct"], alpha=0.7)
        ax.set_xlabel("Yield (t/ha)")
        ax.set_ylabel("Protein (%)")
        ax.set_title("Yield vs Protein")
        fig.tight_layout()
        return fig

    @output
    @render.table
    def preview():
        return make_data().head(10)


app = App(app_ui, server)
