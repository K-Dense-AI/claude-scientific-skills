#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    accessions = pd.DataFrame(
        {
            "accession_id": ["BB001", "BB002", "BB003"],
            "name": ["Elite-A", "Elite-B", "Donor-X"],
            "status": ["active", "active", "active"],
        }
    )
    accessions.to_csv(out / "breedbase_accessions_read.csv", index=False)

    payload = pd.DataFrame(
        {
            "accession_id": ["BB004", "BB005"],
            "name": ["Line-4", "Line-5"],
            "write_result": ["mock_created", "mock_created"],
        }
    )
    payload.to_csv(out / "breedbase_accessions_write_result.csv", index=False)
    print("Saved Breedbase read/write mock outputs")


if __name__ == "__main__":
    main()
