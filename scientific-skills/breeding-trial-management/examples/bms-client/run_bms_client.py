#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    trials = pd.DataFrame(
        {
            "trial_id": ["BMS-T01", "BMS-T02"],
            "location": ["Site-A", "Site-B"],
            "season": ["2026A", "2026A"],
            "status": ["retrieved_mock", "retrieved_mock"],
        }
    )
    trials.to_csv(out / "bms_trials_read.csv", index=False)

    update = pd.DataFrame(
        {
            "trial_id": ["BMS-T03"],
            "operation": ["create_mock"],
            "result": ["success"],
        }
    )
    update.to_csv(out / "bms_trials_write_result.csv", index=False)
    print("Saved BMS client mock read/write outputs")


if __name__ == "__main__":
    main()
