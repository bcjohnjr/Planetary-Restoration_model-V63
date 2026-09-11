#!/usr/bin/env python3
"""Synthesize V63 publication experiments without converting disagreement into a failure."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ORDER = [
    "canonical_peak15p2",
    "scaled_peak5p0",
    "scaled_peak10p0",
    "scaled_peak20p0",
    "same_total_uniform",
    "same_total_peak20_early",
    "same_total_peak20_late",
]


def load_jsons(root: Path, name: str):
    return [json.loads(p.read_text(encoding="utf-8")) for p in root.rglob(name)]


def milestone(summary: dict, year: int) -> dict:
    for row in summary["milestones"]:
        if int(row["year"]) == year:
            return row
    raise KeyError((summary.get("case"), year))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="experiment_artifacts")
    ap.add_argument("--outdir", default="nature_experiment_results")
    args = ap.parse_args()

    root = Path(args.root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cases = load_jsons(root, "schedule_case_summary.json")
    if not cases:
        raise RuntimeError("No schedule_case_summary.json files found")
    by_case = {x["case"]: x for x in cases}
    missing = [x for x in ORDER if x not in by_case]
    if missing:
        raise RuntimeError(f"Missing schedule cases: {missing}")

    oscar_list = load_jsons(root, "oscar_matched_summary.json")
    if len(oscar_list) != 1:
        raise RuntimeError(f"Expected one OSCAR matched summary, found {len(oscar_list)}")
    oscar = oscar_list[0]

    rows = []
    for case in ORDER:
        s = by_case[case]
        m = milestone(s, 2183)
        rows.append({
            "case": case,
            "cumulative_cdr_gtco2": s["cumulative_cdr_gtco2"],
            "peak_cdr_gtco2_per_year": s["peak_cdr_gtco2_per_year"],
            "first_cdr_year": s["first_positive_cdr_year"],
            "last_cdr_year": s["last_positive_cdr_year"],
            "first_persistent_divergence_year": s["first_5yr_persistent_hector_outside_fair_p05_p95"],
            "fair_response_2183": m["fair_response_fraction_p50"],
            "hector_response_2183": m["hector_response_fraction"],
            "hector_minus_fair_2183": m["hector_response_fraction"] - m["fair_response_fraction_p50"],
            "fair_delta_co2_2183_ppm": m["fair_delta_co2_p50_ppm"],
            "hector_delta_co2_2183_ppm": m["hector_delta_co2_ppm"],
        })
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "cdr_schedule_generality_summary.csv", index=False)

    canonical = by_case["canonical_peak15p2"]
    cm = milestone(canonical, 2183)
    three_model = {
        "year": 2183,
        "FaIR_response_fraction_p50": cm["fair_response_fraction_p50"],
        "Hector_response_fraction": cm["hector_response_fraction"],
        "OSCAR_response_fraction_p50": oscar["final_2183_response_fraction_p50"],
        "FaIR_delta_co2_p50_ppm": cm["fair_delta_co2_p50_ppm"],
        "Hector_delta_co2_ppm": cm["hector_delta_co2_ppm"],
        "OSCAR_delta_co2_p50_ppm": oscar["final_2183_delta_co2_p50_ppm"],
    }

    fixed_total_names = [
        "canonical_peak15p2",
        "same_total_uniform",
        "same_total_peak20_early",
        "same_total_peak20_late",
    ]
    fixed_total = table[table.case.isin(fixed_total_names)].copy()
    total_span = float(fixed_total.cumulative_cdr_gtco2.max() - fixed_total.cumulative_cdr_gtco2.min())

    result = {
        "experiment": "V63 Nature-candidate generality and mechanism suite",
        "schedule_cases": rows,
        "canonical_three_model_2183": three_model,
        "oscar_status": oscar["status"],
        "oscar_forcing_protocol": oscar["forcing_protocol"],
        "fixed_total_design_check": {
            "cases": fixed_total_names,
            "max_minus_min_cumulative_cdr_gtco2": total_span,
            "passes_exact_total_design": total_span < 1e-6,
        },
        "predeclared_questions": [
            "Does FaIR/Hector long-horizon separation persist when CDR amplitude changes?",
            "At fixed cumulative CDR, does changing removal timing change atmospheric response efficacy?",
            "Where does matched-forcing OSCAR fall relative to FaIR and Hector in the canonical case?",
            "How much of OSCAR's CDR compensation is assigned to land versus ocean sinks?",
        ],
        "claim_rule": (
            "A universal threshold is not claimed. A Nature-level mechanistic claim requires reproducible behavior across schedules "
            "and interpretable reservoir compensation, with model structural uncertainty reported explicitly."
        ),
    }
    (outdir / "nature_experiment_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# V63 Nature-candidate experiment results",
        "",
        "## CDR schedule generality",
        "",
        "| Case | Cumulative CDR (GtCO2) | Peak (GtCO2/yr) | First persistent FaIR/Hector range separation | FaIR response 2183 | Hector response 2183 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        div = "none" if r["first_persistent_divergence_year"] is None else str(r["first_persistent_divergence_year"])
        lines.append(
            f"| {r['case']} | {r['cumulative_cdr_gtco2']:.3f} | {r['peak_cdr_gtco2_per_year']:.3f} | {div} | "
            f"{r['fair_response_2183']:.4f} | {r['hector_response_2183']:.4f} |"
        )
    lines += [
        "",
        "## Canonical three-model comparison at 2183",
        "",
        f"- FaIR response fraction (median): **{three_model['FaIR_response_fraction_p50']:.4f}**",
        f"- Hector response fraction: **{three_model['Hector_response_fraction']:.4f}**",
        f"- OSCAR response fraction (median): **{three_model['OSCAR_response_fraction_p50']:.4f}**",
        "",
        "## Mechanism test",
        "",
        "See `oscar_reservoir_decomposition.csv` for annual atmospheric benefit, land compensation, ocean compensation, and budget residual.",
        "",
        "No universal threshold is inferred automatically from these experiments; model disagreement is retained as a result.",
    ]
    (outdir / "NATURE_EXPERIMENTS_V63_RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
