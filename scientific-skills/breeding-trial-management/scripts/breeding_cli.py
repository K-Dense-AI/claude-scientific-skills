#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

import argparse


def cmd_design(args):
    print(f"design stub: mode={args.mode}")


def cmd_fieldbook(args):
    print(f"fieldbook stub: input={args.input}")


def cmd_germplasm(args):
    print(f"germplasm stub: action={args.action}")


def cmd_select(args):
    print(f"select stub: trait={args.trait}")


def cmd_cross(args):
    print(f"cross stub: strategy={args.strategy}")


def main():
    p = argparse.ArgumentParser(description="Breeding trial management CLI")
    sp = p.add_subparsers(dest="command", required=True)

    d = sp.add_parser("design", help="Trial design operations")
    d.add_argument(
        "--mode", choices=["rcbd", "alpha-lattice", "augmented"], default="rcbd"
    )
    d.set_defaults(func=cmd_design)

    f = sp.add_parser("fieldbook", help="Fieldbook generation")
    f.add_argument("--input", default="plots.csv")
    f.set_defaults(func=cmd_fieldbook)

    g = sp.add_parser("germplasm", help="Germplasm operations")
    g.add_argument("--action", choices=["list", "import", "export"], default="list")
    g.set_defaults(func=cmd_germplasm)

    s = sp.add_parser("select", help="Selection operations")
    s.add_argument("--trait", default="yield")
    s.set_defaults(func=cmd_select)

    c = sp.add_parser("cross", help="Cross planning operations")
    c.add_argument(
        "--strategy",
        choices=["top-cross", "factorial", "optimal-contribution"],
        default="top-cross",
    )
    c.set_defaults(func=cmd_cross)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
