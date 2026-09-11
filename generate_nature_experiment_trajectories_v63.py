#!/usr/bin/env python3
"""Generate predeclared CDR schedule experiments for the V62.2 Nature-candidate tests.

The non-CDR source terms are held exactly as in the canonical V62.2 trajectory.
Only cdr_gtco2 is changed. Every generated CSV spans 2026-2183 and obeys
net = gross + permafrost + stored-carbon reversal - CDR.

Experiment families
-------------------
1. Dose/rate: canonical CDR shape scaled to 5, 10, 15.2 and 20 GtCO2/yr peaks.
2. Timing at fixed cumulative removal: uniform, early 20-Gt pulse, late 20-Gt pulse.

The fixed-cumulative cases all remove exactly the same total CO2 as the canonical
2026-2183 schedule. This separates timing/rate effects from total-removal effects.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED = [
    "year",
    "gross_co2_gtco2",
    "permafrost_co2_gtco2",
    "stored_carbon_reversal_gtco2",
    "cdr_gtco2",
    "net_co2_gtco2",
]


def validate(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing columns: {missing}")
    x = df.sort_values("year").reset_index(drop=True).copy()
    if int(x.year.iloc[0]) != 2026 or int(x.year.iloc[-1]) != 2183:
        raise RuntimeError("Expected canonical years 2026-2183")
    if x.year.duplicated().any():
        raise RuntimeError("Duplicate years")
    ident = (
        x.gross_co2_gtco2
        + x.permafrost_co2_gtco2
        + x.stored_carbon_reversal_gtco2
        - x.cdr_gtco2
    )
    err = float(np.max(np.abs(ident - x.net_co2_gtco2)))
    if err > 1e-9:
        raise RuntimeError(f"Canonical mass-balance error: {err}")
    if (x.cdr_gtco2 < -1e-12).any():
        raise RuntimeError("Negative canonical CDR")
    return x


def pulse_schedule(years: np.ndarray, total: float, peak: float, late: bool) -> np.ndarray:
    """Exact-total rectangular pulse with one fractional year, excluding 2026."""
    out = np.zeros(len(years), dtype=float)
    eligible = np.where(years >= 2027)[0]
    n_full = int(total // peak)
    rem = float(total - n_full * peak)
    n_needed = n_full + (1 if rem > 1e-12 else 0)
    if n_needed > len(eligible):
        raise RuntimeError("Pulse cannot fit inside 2027-2183")
    chosen = eligible[-n_needed:] if late else eligible[:n_needed]
    if late:
        # Fractional year first, then full-rate years, so the late pulse ends at 2183.
        if rem > 1e-12:
            out[chosen[0]] = rem
            out[chosen[1:]] = peak
        else:
            out[chosen] = peak
    else:
        out[chosen[:n_full]] = peak
        if rem > 1e-12:
            out[chosen[n_full]] = rem
    return out


def finalize(base: pd.DataFrame, cdr: np.ndarray) -> pd.DataFrame:
    out = base[REQUIRED].copy()
    out["cdr_gtco2"] = np.asarray(cdr, dtype=float)
    out["net_co2_gtco2"] = (
        out.gross_co2_gtco2
        + out.permafrost_co2_gtco2
        + out.stored_carbon_reversal_gtco2
        - out.cdr_gtco2
    )
    if not np.isfinite(out[REQUIRED].to_numpy(float)).all():
        raise RuntimeError("Non-finite generated trajectory")
    if (out.cdr_gtco2 < -1e-12).any():
        raise RuntimeError("Generated negative CDR")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="external_validation_net_co2_trajectory.csv")
    ap.add_argument("--outdir", default="nature_experiments/trajectories")
    args = ap.parse_args()

    source = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base = validate(pd.read_csv(source))
    years = base.year.to_numpy(int)
    canonical = base.cdr_gtco2.to_numpy(float)
    canonical_total = float(canonical.sum())
    canonical_peak = float(canonical.max())
    if abs(canonical_peak - 15.2) > 1e-9:
        raise RuntimeError(f"Expected canonical peak 15.2, found {canonical_peak}")

    variants: dict[str, tuple[np.ndarray, dict]] = {}

    def add(name: str, schedule: np.ndarray, family: str, purpose: str) -> None:
        variants[name] = (
            np.asarray(schedule, dtype=float),
            {"family": family, "purpose": purpose},
        )

    add(
        "canonical_peak15p2",
        canonical,
        "canonical",
        "Exact V62.2 canonical pathway; reference case.",
    )

    for peak in (5.0, 10.0, 20.0):
        add(
            f"scaled_peak{str(peak).replace('.', 'p')}",
            canonical * (peak / canonical_peak),
            "dose_rate",
            f"Preserve canonical time-shape while changing peak CDR to {peak:g} GtCO2/yr.",
        )

    active = years >= 2027
    uniform = np.zeros(len(years), dtype=float)
    uniform[active] = canonical_total / int(active.sum())
    add(
        "same_total_uniform",
        uniform,
        "timing_fixed_total",
        "Same canonical cumulative CDR distributed uniformly from 2027 through 2183.",
    )

    add(
        "same_total_peak20_early",
        pulse_schedule(years, canonical_total, 20.0, late=False),
        "timing_fixed_total",
        "Same canonical cumulative CDR delivered as early as possible at <=20 GtCO2/yr.",
    )
    add(
        "same_total_peak20_late",
        pulse_schedule(years, canonical_total, 20.0, late=True),
        "timing_fixed_total",
        "Same canonical cumulative CDR delivered as late as possible at <=20 GtCO2/yr.",
    )

    manifest = {
        "experiment": "V63 Nature-candidate CDR schedule generality",
        "source": source.name,
        "years": [2026, 2183],
        "canonical_total_cdr_gtco2": canonical_total,
        "canonical_peak_cdr_gtco2_per_year": canonical_peak,
        "design": {
            "dose_rate_family": "same canonical shape, different amplitudes",
            "timing_fixed_total_family": "same cumulative removal, different timing/rate",
            "held_constant": [
                "gross_co2_gtco2",
                "permafrost_co2_gtco2",
                "stored_carbon_reversal_gtco2",
            ],
        },
        "cases": [],
    }

    for name, (cdr, meta) in variants.items():
        out = finalize(base, cdr)
        path = outdir / f"{name}.csv"
        out.to_csv(path, index=False)
        row = {
            "case": name,
            **meta,
            "file": str(path),
            "cumulative_cdr_gtco2": float(out.cdr_gtco2.sum()),
            "peak_cdr_gtco2_per_year": float(out.cdr_gtco2.max()),
            "first_positive_cdr_year": int(out.loc[out.cdr_gtco2 > 0, "year"].min()),
            "last_positive_cdr_year": int(out.loc[out.cdr_gtco2 > 0, "year"].max()),
        }
        manifest["cases"].append(row)

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    pd.DataFrame(manifest["cases"]).to_csv(outdir / "manifest.csv", index=False)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
