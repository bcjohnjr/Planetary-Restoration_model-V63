#!/usr/bin/env python3
"""Summarize one 2026-2183 FaIR/Hector CDR-schedule generality case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def nearest(df: pd.DataFrame, year: int) -> pd.Series:
    return df.iloc[int(np.argmin(np.abs(df.year.to_numpy(float) - year)))]


def first_persistent(mask: np.ndarray, years: np.ndarray, n: int = 5):
    run = 0
    for i, flag in enumerate(mask):
        run = run + 1 if flag else 0
        if run >= n:
            return int(years[i - n + 1])
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--trajectory", default="external_validation_net_co2_trajectory.csv")
    ap.add_argument("--fair", default="fair_matched_attribution.csv")
    ap.add_argument("--hector", default="hector_matched_attribution.csv")
    ap.add_argument("--fair-summary", default="fair_matched_summary.json")
    ap.add_argument("--hector-summary", default="hector_matched_summary.json")
    args = ap.parse_args()

    tr = pd.read_csv(args.trajectory)
    f = pd.read_csv(args.fair)
    h = pd.read_csv(args.hector)
    fs = json.loads(Path(args.fair_summary).read_text(encoding="utf-8"))
    hs = json.loads(Path(args.hector_summary).read_text(encoding="utf-8"))

    # Integer-year alignment for 2027-2183 only.
    ff = f.copy()
    ff["iyear"] = np.rint(ff.year.to_numpy(float)).astype(int)
    hh = h.copy()
    hh["iyear"] = np.rint(hh.year.to_numpy(float)).astype(int)
    ff = ff[(ff.iyear >= 2027) & (ff.iyear <= 2183)].drop_duplicates("iyear")
    hh = hh[(hh.iyear >= 2027) & (hh.iyear <= 2183)].drop_duplicates("iyear")
    m = ff[["iyear", "fraction_p05", "fraction_p50", "fraction_p95", "delta_co2_p50_ppm"]].merge(
        hh[["iyear", "apparent_response_fraction", "delta_co2_ppm"]], on="iyear", how="inner"
    )
    m = m[np.isfinite(m.fraction_p50) & np.isfinite(m.apparent_response_fraction)].copy()
    outside = (
        (m.apparent_response_fraction < m.fraction_p05)
        | (m.apparent_response_fraction > m.fraction_p95)
    ).to_numpy(bool)
    persistent_year = first_persistent(outside, m.iyear.to_numpy(int), 5)

    milestones = []
    for y in (2040, 2100, 2156, 2183):
        fr = nearest(f, y)
        hr = nearest(h, y)
        milestones.append({
            "year": y,
            "fair_delta_co2_p50_ppm": float(fr.delta_co2_p50_ppm),
            "hector_delta_co2_ppm": float(hr.delta_co2_ppm),
            "fair_response_fraction_p05": float(fr.fraction_p05),
            "fair_response_fraction_p50": float(fr.fraction_p50),
            "fair_response_fraction_p95": float(fr.fraction_p95),
            "hector_response_fraction": float(hr.apparent_response_fraction),
            "hector_inside_fair_p05_p95": bool(
                fr.fraction_p05 <= hr.apparent_response_fraction <= fr.fraction_p95
            ),
        })

    cdr = tr.cdr_gtco2.to_numpy(float)
    summary = {
        "experiment": "V63 CDR schedule generality — matched FaIR/Hector through 2183",
        "case": args.case,
        "trajectory_years": [int(tr.year.min()), int(tr.year.max())],
        "cumulative_cdr_gtco2": float(cdr.sum()),
        "peak_cdr_gtco2_per_year": float(cdr.max()),
        "first_positive_cdr_year": int(tr.loc[tr.cdr_gtco2 > 0, "year"].min()),
        "last_positive_cdr_year": int(tr.loc[tr.cdr_gtco2 > 0, "year"].max()),
        "first_5yr_persistent_hector_outside_fair_p05_p95": persistent_year,
        "milestones": milestones,
        "forcing_gates": {
            "fair_on_max_abs_wm2": fs["forcing_match_max_abs_on_wm2"],
            "fair_off_max_abs_wm2": fs["forcing_match_max_abs_off_wm2"],
            "hector_on_max_abs_wm2": hs["forcing_match_max_abs_on_wm2"],
            "hector_off_max_abs_wm2": hs["forcing_match_max_abs_off_wm2"],
        },
        "claim_boundary": (
            "This case tests sensitivity to CDR timing/rate while holding the V62 non-CDR source terms fixed. "
            "Cross-model disagreement is an experimental result, not a failure gate."
        ),
    }
    Path("schedule_case_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame(milestones).to_csv("schedule_case_milestones.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
