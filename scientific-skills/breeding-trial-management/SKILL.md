---
name: breeding-trial-management
description: >
  Breeding trial management toolkit for design, fieldbook generation, germplasm workflows,
  selection, and crossing decisions. Complements qtl-analysis for end-to-end breeding genomics.
license: Apache-2.0
metadata:
  author: Clayton Young (borealBytes / Superior Byte Works, LLC)
  version: "0.1.0"
  category: breeding-genomics
---

# Breeding Trial Management

## Overview

This skill handles breeding operations around trial execution and decision support:
- Trial design planning
- Fieldbook preparation
- Germplasm tracking
- Selection and crossing workflows

## Installation

```bash
python scripts/check_system.py
```

## Unified CLI

```bash
python scripts/breeding_cli.py design --help
python scripts/breeding_cli.py fieldbook --help
python scripts/breeding_cli.py germplasm --help
python scripts/breeding_cli.py select --help
python scripts/breeding_cli.py cross --help
```

## Tool Selection Guide

- Use `design` for RCBD/alpha-lattice/augmented planning tasks.
- Use `fieldbook` to generate plot sheets and labels.
- Use `germplasm` for accession list and inventory operations.
- Use `select` for ranking and shortlist generation.
- Use `cross` for mate pairing and crossing plan scaffolds.
