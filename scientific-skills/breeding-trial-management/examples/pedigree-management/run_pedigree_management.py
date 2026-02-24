#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    ped = pd.DataFrame(
        {
            "id": ["P1", "P2", "F1", "BC1", "BC2", "BC3"],
            "sire": ["0", "0", "P1", "F1", "BC1", "BC2"],
            "dam": ["0", "0", "P2", "P1", "P1", "P1"],
        }
    )
    ped.to_csv(out / "pedigree_records.csv", index=False)

    inbreed = pd.DataFrame(
        {
            "id": ped["id"],
            "inbreeding_coef": [0.0, 0.0, 0.0, 0.125, 0.1875, 0.2344],
        }
    )
    inbreed.to_csv(out / "inbreeding_coefficients.csv", index=False)

    g = nx.DiGraph()
    for _, r in ped.iterrows():
        g.add_node(r["id"])
        if r["sire"] != "0":
            g.add_edge(r["sire"], r["id"])
        if r["dam"] != "0":
            g.add_edge(r["dam"], r["id"])
    pos = nx.spring_layout(g, seed=7)
    plt.figure(figsize=(6, 4))
    nx.draw(g, pos, with_labels=True, node_color="#9ecae1", arrows=True)
    plt.title("Pedigree lineage graph")
    plt.tight_layout()
    plt.savefig(out / "pedigree_lineage.png", dpi=150)
    plt.close()
    print("Saved pedigree records, inbreeding coefficients, and lineage graph")


if __name__ == "__main__":
    main()
