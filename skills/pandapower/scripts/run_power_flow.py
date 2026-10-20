#!/usr/bin/env python3
"""Run an AC power flow and optional N-1 line/transformer screening on a JSON case.

The case file describes buses, lines, two-winding transformers, loads, static
generators and PV generators in engineering units (kV, MW, Mvar, ohm, nF).
The script builds the network with pandapower, solves it, and reports bus
voltages, branch loadings and limit violations. Exit codes:

    0  converged, no violations
    2  converged, at least one limit violation
    3  case rejected or power flow did not converge

See references/case-format.md for the field-by-field contract.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_VIOLATION = 2
EXIT_FAILED = 3


def _require(mapping: dict, key: str, where: str):
    if key not in mapping:
        raise ValueError(f"{where}: missing required field '{key}'")
    return mapping[key]


def build_network(case: dict):
    """Translate a validated case description into a pandapower network."""
    import pandapower as pp

    net = pp.create_empty_network(sn_mva=float(case.get("base_mva", 100.0)))
    index: dict[str, int] = {}

    for i, bus in enumerate(case.get("buses", [])):
        name = _require(bus, "id", f"buses[{i}]")
        if name in index:
            raise ValueError(f"buses[{i}]: duplicate bus id '{name}'")
        index[name] = pp.create_bus(
            net, vn_kv=float(_require(bus, "vn_kv", f"buses[{i}]")), name=name
        )
    if not index:
        raise ValueError("case has no buses")

    def bus_of(ref, where):
        if ref not in index:
            raise ValueError(f"{where}: unknown bus '{ref}'")
        return index[ref]

    for i, grid in enumerate(case.get("ext_grids", [])):
        where = f"ext_grids[{i}]"
        pp.create_ext_grid(
            net,
            bus=bus_of(_require(grid, "bus", where), where),
            vm_pu=float(grid.get("vm_pu", 1.0)),
            va_degree=float(grid.get("va_degree", 0.0)),
            name=grid.get("id", f"ext{i}"),
        )

    for i, gen in enumerate(case.get("gens", [])):
        where = f"gens[{i}]"
        pp.create_gen(
            net,
            bus=bus_of(_require(gen, "bus", where), where),
            p_mw=float(_require(gen, "p_mw", where)),
            vm_pu=float(gen.get("vm_pu", 1.0)),
            min_q_mvar=gen.get("min_q_mvar"),
            max_q_mvar=gen.get("max_q_mvar"),
            name=gen.get("id", f"gen{i}"),
        )

    for i, sgen in enumerate(case.get("sgens", [])):
        where = f"sgens[{i}]"
        pp.create_sgen(
            net,
            bus=bus_of(_require(sgen, "bus", where), where),
            p_mw=float(_require(sgen, "p_mw", where)),
            q_mvar=float(sgen.get("q_mvar", 0.0)),
            name=sgen.get("id", f"sgen{i}"),
        )

    for i, load in enumerate(case.get("loads", [])):
        where = f"loads[{i}]"
        pp.create_load(
            net,
            bus=bus_of(_require(load, "bus", where), where),
            p_mw=float(_require(load, "p_mw", where)),
            q_mvar=float(load.get("q_mvar", 0.0)),
            name=load.get("id", f"load{i}"),
        )

    for i, line in enumerate(case.get("lines", [])):
        where = f"lines[{i}]"
        pp.create_line_from_parameters(
            net,
            from_bus=bus_of(_require(line, "from", where), where),
            to_bus=bus_of(_require(line, "to", where), where),
            length_km=float(_require(line, "length_km", where)),
            r_ohm_per_km=float(_require(line, "r_ohm_per_km", where)),
            x_ohm_per_km=float(_require(line, "x_ohm_per_km", where)),
            c_nf_per_km=float(line.get("c_nf_per_km", 0.0)),
            max_i_ka=float(_require(line, "max_i_ka", where)),
            parallel=int(line.get("parallel", 1)),
            name=line.get("id", f"line{i}"),
        )

    for i, trafo in enumerate(case.get("trafos", [])):
        where = f"trafos[{i}]"
        pp.create_transformer_from_parameters(
            net,
            hv_bus=bus_of(_require(trafo, "hv_bus", where), where),
            lv_bus=bus_of(_require(trafo, "lv_bus", where), where),
            sn_mva=float(_require(trafo, "sn_mva", where)),
            vn_hv_kv=float(_require(trafo, "vn_hv_kv", where)),
            vn_lv_kv=float(_require(trafo, "vn_lv_kv", where)),
            vk_percent=float(_require(trafo, "vk_percent", where)),
            vkr_percent=float(_require(trafo, "vkr_percent", where)),
            pfe_kw=float(trafo.get("pfe_kw", 0.0)),
            i0_percent=float(trafo.get("i0_percent", 0.0)),
            name=trafo.get("id", f"trafo{i}"),
        )

    return net


def solve(net) -> None:
    import pandapower as pp

    # numba is optional; turning it off keeps the run deterministic on machines
    # without llvmlite and avoids the import-time warning.
    pp.runpp(net, numba=False, calculate_voltage_angles=True, init="auto")


def collect(net, vm_min: float, vm_max: float, loading_limit: float) -> dict:
    buses = []
    violations = []
    for idx, row in net.res_bus.iterrows():
        name = net.bus.at[idx, "name"]
        vm = float(row["vm_pu"])
        entry = {
            "bus": name,
            "vn_kv": round(float(net.bus.at[idx, "vn_kv"]), 3),
            "vm_pu": None if vm != vm else round(vm, 5),  # NaN marks an isolated bus
            "va_degree": round(float(row["va_degree"]), 4),
        }
        buses.append(entry)
        if vm != vm:
            violations.append(
                {
                    "element": f"bus {name}",
                    "quantity": "vm_pu",
                    "value": None,
                    "limit": vm_min,
                    "kind": "de_energized",
                }
            )
        elif not vm_min <= vm <= vm_max:
            violations.append(
                {
                    "element": f"bus {name}",
                    "quantity": "vm_pu",
                    "value": entry["vm_pu"],
                    "limit": vm_min if row["vm_pu"] < vm_min else vm_max,
                    "kind": "undervoltage" if row["vm_pu"] < vm_min else "overvoltage",
                }
            )

    branches = []
    branch_columns = {
        "line": (net.res_line, net.line, "p_from_mw", "q_from_mvar"),
        "trafo": (net.res_trafo, net.trafo, "p_hv_mw", "q_hv_mvar"),
    }
    for kind, (res, table, p_col, q_col) in branch_columns.items():
        if len(table) == 0:
            continue
        for idx, row in res.iterrows():
            loading = float(row["loading_percent"])
            entry = {
                "element": f"{kind} {table.at[idx, 'name']}",
                "loading_percent": round(loading, 3),
                "p_from_mw": round(float(row[p_col]), 4),
                "q_from_mvar": round(float(row[q_col]), 4),
            }
            branches.append(entry)
            if loading > loading_limit:
                violations.append(
                    {
                        "element": entry["element"],
                        "quantity": "loading_percent",
                        "value": entry["loading_percent"],
                        "limit": loading_limit,
                        "kind": "overload",
                    }
                )

    energized = [b["vm_pu"] for b in buses if b["vm_pu"] is not None]
    return {
        "converged": bool(net.converged),
        "iterations": int(getattr(net, "_ppc", {}).get("iterations", 0) or 0),
        "buses": buses,
        "branches": branches,
        "violations": violations,
        "summary": {
            "vm_min_pu": min(energized) if energized else None,
            "vm_max_pu": max(energized) if energized else None,
            "de_energized_buses": sum(1 for b in buses if b["vm_pu"] is None),
            "max_loading_percent": round(
                max((b["loading_percent"] for b in branches), default=0.0), 3
            ),
            "p_loss_mw": round(
                float(net.res_line.p_from_mw.sum() + net.res_line.p_to_mw.sum())
                + (
                    float(net.res_trafo.p_hv_mw.sum() + net.res_trafo.p_lv_mw.sum())
                    if len(net.res_trafo)
                    else 0.0
                ),
                4,
            ),
        },
    }


def screen_n1(case: dict, net, vm_min: float, vm_max: float, loading_limit: float,
              max_outages: int) -> dict:
    """Trip each line and transformer in turn and record the post-outage state."""
    results = []
    targets = [("line", i) for i in net.line.index] + [
        ("trafo", i) for i in net.trafo.index
    ]
    if max_outages and len(targets) > max_outages:
        targets = targets[:max_outages]
    for kind, idx in targets:
        table = getattr(net, kind)
        original = bool(table.at[idx, "in_service"])
        name = table.at[idx, "name"]
        table.at[idx, "in_service"] = False
        try:
            solve(net)
            state = collect(net, vm_min, vm_max, loading_limit)
            results.append(
                {
                    "outage": f"{kind} {name}",
                    "converged": state["converged"],
                    "violations": len(state["violations"]),
                    "vm_min_pu": state["summary"]["vm_min_pu"],
                    "de_energized_buses": state["summary"]["de_energized_buses"],
                    "max_loading_percent": state["summary"]["max_loading_percent"],
                    "worst": [
                        f"{v['element']} de-energized"
                        if v["kind"] == "de_energized"
                        else f"{v['element']} {v['quantity']}={v['value']}"
                        for v in state["violations"][:5]
                    ],
                }
            )
        except Exception as exc:  # noqa: BLE001 - a failed contingency is a result
            results.append(
                {"outage": f"{kind} {name}", "converged": False, "error": str(exc)}
            )
        finally:
            table.at[idx, "in_service"] = original
    return {
        "screenings": results,
        "non_converged": sum(1 for r in results if not r.get("converged", False)),
        "with_violations": sum(1 for r in results if r.get("violations", 0)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="AC power flow with limit checks and N-1 screening (pandapower)."
    )
    parser.add_argument("case", help="path to the JSON case file")
    parser.add_argument("--json-out", help="write the full report to this path")
    parser.add_argument("--n1", action="store_true", help="run N-1 line/transformer screening")
    parser.add_argument(
        "--max-n1", type=int, default=40, help="cap the number of N-1 screenings (default 40)"
    )
    parser.add_argument("--vm-min", type=float, default=0.95, help="lower voltage limit in pu")
    parser.add_argument("--vm-max", type=float, default=1.05, help="upper voltage limit in pu")
    parser.add_argument(
        "--loading-limit", type=float, default=100.0, help="branch loading limit in percent"
    )
    args = parser.parse_args(argv)

    try:
        case = json.loads(Path(args.case).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[FAIL] cannot read case {args.case}: {exc}")
        return EXIT_FAILED

    try:
        import pandapower
        import numpy
        import pandas

        net = build_network(case)
        solve(net)
    except ValueError as exc:
        print(f"[FAIL] invalid case: {exc}")
        return EXIT_FAILED
    except Exception as exc:  # noqa: BLE001 - solver failures are reported, not raised
        print(f"[FAIL] power flow failed: {exc}")
        return EXIT_FAILED

    report = {
        "case": case.get("name", Path(args.case).stem),
        "case_file": str(Path(args.case)),
        "versions": {
            "pandapower": pandapower.__version__,
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
        },
        "limits": {
            "vm_min_pu": args.vm_min,
            "vm_max_pu": args.vm_max,
            "loading_percent": args.loading_limit,
        },
        "power_flow": collect(net, args.vm_min, args.vm_max, args.loading_limit),
    }
    if args.n1:
        report["n1"] = screen_n1(
            case, net, args.vm_min, args.vm_max, args.loading_limit, args.max_n1
        )

    pf = report["power_flow"]
    print(
        f"[{'OK' if pf['converged'] else 'FAIL'}] {report['case']}: "
        f"pandapower {report['versions']['pandapower']}, "
        f"{pf['iterations']} iterations"
    )
    print(f"{'bus':<14}{'vn_kv':>8}{'vm_pu':>9}{'va_deg':>10}")
    for bus in pf["buses"]:
        vm_text = "n/a" if bus["vm_pu"] is None else f"{bus['vm_pu']:.4f}"
        print(f"{bus['bus']:<14}{bus['vn_kv']:>8.1f}{vm_text:>9}{bus['va_degree']:>10.3f}")
    for branch in pf["branches"]:
        print(f"  {branch['element']:<18} loading {branch['loading_percent']:>7.2f}%")
    summary = pf["summary"]
    vm_range = (
        "no energized bus"
        if summary["vm_min_pu"] is None
        else f"{summary['vm_min_pu']:.4f}..{summary['vm_max_pu']:.4f} pu"
    )
    print(
        f"-> vm {vm_range}, max loading {summary['max_loading_percent']:.2f}%, "
        f"losses {summary['p_loss_mw']:.4f} MW, "
        f"de-energized buses {summary['de_energized_buses']}"
    )
    for violation in pf["violations"]:
        if violation["kind"] == "de_energized":
            print(
                f"[FAIL] de-energized: {violation['element']} has no voltage source "
                "after the outage"
            )
            continue
        print(
            f"[FAIL] {violation['kind']}: {violation['element']} "
            f"{violation['quantity']}={violation['value']} (limit {violation['limit']})"
        )
    if args.n1:
        n1 = report["n1"]
        print(
            f"-> N-1: {len(n1['screenings'])} outages, "
            f"{n1['with_violations']} with violations, "
            f"{n1['non_converged']} non-converged"
        )
        for item in n1["screenings"]:
            if not item.get("converged", False):
                print(f"[FAIL] N-1 {item['outage']}: {item.get('error', 'diverged')}")
            elif item["violations"]:
                print(
                    f"[FAIL] N-1 {item['outage']}: {item['violations']} violation(s) "
                    f"({'; '.join(item['worst'])})"
                )
        if n1["with_violations"] or n1["non_converged"]:
            report["power_flow"]["violations"].append(
                {"element": "N-1 screening", "quantity": "contingencies",
                 "value": n1["with_violations"] + n1["non_converged"], "limit": 0,
                 "kind": "n1_failure"}
            )

    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"-> report written to {out}")

    if not pf["converged"]:
        return EXIT_FAILED
    return EXIT_VIOLATION if report["power_flow"]["violations"] else EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
