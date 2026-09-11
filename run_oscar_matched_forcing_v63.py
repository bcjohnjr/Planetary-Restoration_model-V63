#!/usr/bin/env python3
"""OSCAR v3.3 cross-model matched-forcing + reservoir-decomposition experiment.

This is a publication experiment, not a replacement for the already-green V62.2
OSCAR validation.  It uses the same FaIR-derived non-CO2 forcing target used by
the green FaIR/Hector matched experiment.

OSCAR's RF_contr driver is deliberately repurposed as a *synthetic additive
bookkeeping forcing carrier*.  It must not be interpreted as physical contrail
forcing in this experiment.  The original scenario RF_contr is zeroed from 2027.
The synthetic carrier is iterated until the OSCAR ensemble-median
(RF_warm - RF_CO2) matches the common target.

The script also decomposes the paired CDR response into atmospheric benefit,
land-sink compensation and ocean-sink compensation.  OSCAR's endogenous
permafrost stock is disabled because the V62 trajectory already contains an
explicit exogenous permafrost-CO2 source; that source is identical in ON/OFF.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import run_oscar_validation_v62 as base

MATCH_REVISION = "V63-oscar-common-nonco2-median-match-2026-09-11"
DEFAULT_TOL = 1e-4
DEFAULT_MAX_ITERS = 6


def year_matrix(da) -> np.ndarray:
    """Return xarray DataArray as year x flattened-ensemble matrix."""
    if "year" not in da.dims:
        raise RuntimeError(f"Expected year dimension, found {da.dims}")
    order = ["year"] + [d for d in da.dims if d != "year"]
    arr = np.asarray(da.transpose(*order).values, dtype=float)
    return arr.reshape(arr.shape[0], -1)


def qrow(values) -> dict[str, float]:
    a = np.asarray(values, dtype=float).reshape(-1)
    return {
        "p05": float(np.quantile(a, 0.05)),
        "p50": float(np.quantile(a, 0.50)),
        "p95": float(np.quantile(a, 0.95)),
    }


def fit_forcing_case(
    OSCAR,
    ini_2026,
    par,
    forcing_base,
    target: np.ndarray,
    years: np.ndarray,
    label: str,
    tol: float,
    max_iters: int,
):
    """Match OSCAR ensemble-median non-CO2 warming forcing to common target."""
    import xarray as xr

    if "RF_contr" not in forcing_base:
        raise RuntimeError("OSCAR forcing lacks RF_contr synthetic-carrier slot")

    carrier = np.zeros(len(years), dtype=float)
    alpha = 1.0
    alpha_min = 0.0625
    audit: list[dict] = []
    current = None

    keep = [
        "RF_warm",
        "RF_CO2",
        "D_Focean",
        "D_Fland",
        "D_Epf",
        "D_Tg",
    ]

    def evaluate(vec: np.ndarray, iteration: int, tag: str):
        f = forcing_base.copy(deep=True)
        # Explicitly remove scenario contrail forcing and use this variable only
        # as the documented synthetic residual carrier.
        f["RF_contr"] = xr.DataArray(
            np.asarray(vec, dtype=float),
            coords={"year": f.year.values},
            dims=("year",),
        )
        out, fin = OSCAR(
            Ini=ini_2026,
            Par=par,
            For=f,
            var_keep=keep,
            get_final=True,
            nt=4,
            nt_max=48,
            adapt_nt=True,
        )
        non = year_matrix(out["RF_warm"] - out["RF_CO2"])
        med = np.median(non, axis=1)
        residual = target - med
        err = float(np.max(np.abs(residual)))
        audit.append(
            {
                "case": label,
                "iteration": int(iteration),
                "tag": tag,
                "alpha": float(alpha),
                "max_abs_median_nonco2_error_wm2": err,
            }
        )
        print(f"{label}: iteration {iteration} {tag}: max median forcing error={err:.9g} W/m2", flush=True)
        return {
            "carrier": np.asarray(vec, dtype=float),
            "out": out,
            "fin": fin,
            "non": non,
            "median": med,
            "residual": residual,
            "error": err,
        }

    current = evaluate(carrier, 0, "seed")
    if current["error"] > tol:
        # RF_contr enters RF_warm additively, so the undamped residual is the
        # natural first correction. Climate-mediated secondary feedbacks are
        # handled by subsequent accepted-state iterations.
        carrier = carrier + current["residual"]
        current = evaluate(carrier, 1, "direct-residual")

    iteration = 1
    while current["error"] > tol and iteration < max_iters:
        iteration += 1
        proposal = current["carrier"] + alpha * current["residual"]
        trial = evaluate(proposal, iteration, "proposal")
        if trial["error"] < current["error"]:
            current = trial
            alpha = min(1.0, alpha * 1.25)
        else:
            alpha *= 0.5
            if alpha < alpha_min:
                raise RuntimeError(
                    f"{label}: forcing match stalled at {current['error']:.6g} W/m2"
                )

    if current["error"] > tol:
        raise RuntimeError(
            f"{label}: forcing match failed: {current['error']:.6g} > {tol:.6g} W/m2"
        )
    return current, audit


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--oscar-root", required=True)
    ap.add_argument("--trajectory", default="external_validation_net_co2_trajectory.csv")
    ap.add_argument("--forcing-target", default="fair_matched_nonco2_forcing.csv")
    ap.add_argument("--outdir", default="oscar_matched_results")
    ap.add_argument("--nmc", type=int, default=200)
    ap.add_argument("--end-year", type=int, default=2183)
    ap.add_argument("--forcing-tol", type=float, default=DEFAULT_TOL)
    ap.add_argument("--max-iters", type=int, default=DEFAULT_MAX_ITERS)
    args = ap.parse_args()

    launch_dir = Path.cwd().resolve()
    trajectory = Path(args.trajectory)
    if not trajectory.is_absolute():
        trajectory = (launch_dir / trajectory).resolve()
    target_path = Path(args.forcing_target)
    if not target_path.is_absolute():
        target_path = (launch_dir / target_path).resolve()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = (launch_dir / outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if args.end_year != 2183:
        raise RuntimeError("Publication matched-OSCAR experiment is predeclared to end at 2183")

    pathway = base.load_pathway(trajectory, args.end_year)
    target_df = pd.read_csv(target_path).sort_values("year")
    if not {"year", "matched_nonco2_forcing_wm2"}.issubset(target_df.columns):
        raise RuntimeError("Common forcing target has wrong columns")

    (
        OSCAR,
        load_all_param,
        generate_config,
        for_hist0,
        for_scen0,
    ) = base.import_oscar(Path(args.oscar_root).resolve())

    import xarray as xr

    np.random.seed(base.SEED)
    par0 = load_all_param(mod_region="RCP_5reg")
    par = generate_config(par0, nMC=args.nmc)

    # Avoid double counting the V62 exogenous permafrost source.
    if "Cfroz_0" in par:
        par["Cfroz_0"] = 0 * par["Cfroz_0"]

    for_hist = for_hist0.copy(deep=True)
    par, for_hist = base.move_static_to_par(xr, par, for_hist)
    _, ini_2014 = OSCAR(
        Ini=None,
        Par=par,
        For=for_hist,
        get_final=True,
        nt=4,
    )

    scen = for_scen0.sel(scen=base.SCENARIO, drop=True).copy(deep=True)
    bridge = scen.sel(year=slice(2015, 2026)).copy(deep=True)
    if "D_CO2" not in bridge:
        raise RuntimeError("OSCAR bridge lacks D_CO2")
    co2_0 = float(np.asarray(par0["CO2_0"]).squeeze())
    bridge["D_CO2"].loc[dict(year=2026)] = base.COMMON_CO2_PPM - co2_0
    _, ini_2026 = OSCAR(
        Ini=ini_2014,
        Par=par,
        For=bridge,
        get_final=True,
        nt=4,
    )

    future = scen.sel(year=slice(base.START_YEAR, args.end_year)).copy(deep=True)
    if "D_CO2" in future:
        future = future.drop_vars("D_CO2")
    future = base.zero_separate_luc(future)
    years = np.asarray(future.year.values, dtype=int)

    target_series = target_df.set_index("year")["matched_nonco2_forcing_wm2"]
    missing_target = [int(y) for y in years if int(y) not in target_series.index]
    if missing_target:
        raise RuntimeError(f"Forcing target missing years: {missing_target[:5]}")
    target = np.asarray([float(target_series.loc[int(y)]) for y in years], dtype=float)

    p = pathway[pathway.year >= base.START_YEAR].set_index("year")
    on = np.asarray([float(p.loc[y, "net_co2_gtco2"]) for y in years], dtype=float)
    off = np.asarray([
        float(
            p.loc[y, "gross_co2_gtco2"]
            + p.loc[y, "permafrost_co2_gtco2"]
            + p.loc[y, "stored_carbon_reversal_gtco2"]
        )
        for y in years
    ], dtype=float)

    ref_eff = scen["Eff"].sel(year=2025, drop=True)
    for_on = base.set_global_eff(xr, future.copy(deep=True), on, ref_eff)
    for_off = base.set_global_eff(xr, future.copy(deep=True), off, ref_eff)

    # RF_contr differs by construction after fitting, so identity-test all other drivers.
    xr.testing.assert_identical(
        for_on.drop_vars([v for v in ("Eff", "RF_contr") if v in for_on]),
        for_off.drop_vars([v for v in ("Eff", "RF_contr") if v in for_off]),
    )

    print("Matching OSCAR removal-ON to common non-CO2 forcing target", flush=True)
    fit_on, audit_on = fit_forcing_case(
        OSCAR, ini_2026, par, for_on, target, years, "removal_ON",
        args.forcing_tol, args.max_iters,
    )
    print("Matching OSCAR removal-OFF to common non-CO2 forcing target", flush=True)
    fit_off, audit_off = fit_forcing_case(
        OSCAR, ini_2026, par, for_off, target, years, "removal_OFF",
        args.forcing_tol, args.max_iters,
    )

    out_on = fit_on["out"]
    out_off = fit_off["out"]
    co2_on = year_matrix(out_on["D_CO2"] + par["CO2_0"])
    co2_off = year_matrix(out_off["D_CO2"] + par["CO2_0"])
    if co2_on.shape != co2_off.shape:
        raise RuntimeError("OSCAR paired CO2 shapes differ")

    cdr_map = pathway.set_index("year")["cdr_gtco2"].to_dict()
    cdr_annual = np.asarray([float(cdr_map[int(y)]) for y in years], dtype=float)
    cumulative_cdr = np.cumsum(cdr_annual)

    delta_ppm = co2_off - co2_on
    atmospheric_benefit_gtco2 = delta_ppm * base.GTCO2_PER_PPM

    land_on = year_matrix(out_on["D_Fland"])
    land_off = year_matrix(out_off["D_Fland"])
    ocean_on = year_matrix(out_on["D_Focean"])
    ocean_off = year_matrix(out_off["D_Focean"])
    epf_on = year_matrix(out_on["D_Epf"])
    epf_off = year_matrix(out_off["D_Epf"])

    # Positive compensation means the OFF world keeps more carbon in a sink than
    # the ON world; equivalently CDR induces reduced uptake / backflux.
    land_comp_flux_pgcyr = land_off - land_on
    ocean_comp_flux_pgcyr = ocean_off - ocean_on
    land_comp_gtco2 = np.cumsum(land_comp_flux_pgcyr, axis=0) * base.GTCO2_PER_PGC
    ocean_comp_gtco2 = np.cumsum(ocean_comp_flux_pgcyr, axis=0) * base.GTCO2_PER_PGC
    closure_gtco2 = (
        cumulative_cdr[:, None]
        - atmospheric_benefit_gtco2
        - land_comp_gtco2
        - ocean_comp_gtco2
    )

    rows_attr = []
    rows_res = []
    rows_ts = []
    for i, y in enumerate(years):
        dq = qrow(delta_ppm[i])
        resp = delta_ppm[i] * base.GTCO2_PER_PPM / cumulative_cdr[i] if cumulative_cdr[i] > 0 else np.full(delta_ppm.shape[1], np.nan)
        rq = qrow(resp) if cumulative_cdr[i] > 0 else {"p05": None, "p50": None, "p95": None}
        onq = qrow(co2_on[i]); offq = qrow(co2_off[i])
        rows_attr.append({
            "year": int(y),
            "cumulative_cdr_gtco2": float(cumulative_cdr[i]),
            "delta_co2_p05_ppm": dq["p05"],
            "delta_co2_p50_ppm": dq["p50"],
            "delta_co2_p95_ppm": dq["p95"],
            "fraction_p05": rq["p05"],
            "fraction_p50": rq["p50"],
            "fraction_p95": rq["p95"],
        })
        lq = qrow(land_comp_gtco2[i]); oq = qrow(ocean_comp_gtco2[i]); aq = qrow(atmospheric_benefit_gtco2[i]); cq = qrow(closure_gtco2[i])
        rows_res.append({
            "year": int(y),
            "cumulative_cdr_gtco2": float(cumulative_cdr[i]),
            "atmospheric_benefit_p05_gtco2": aq["p05"],
            "atmospheric_benefit_p50_gtco2": aq["p50"],
            "atmospheric_benefit_p95_gtco2": aq["p95"],
            "land_compensation_p05_gtco2": lq["p05"],
            "land_compensation_p50_gtco2": lq["p50"],
            "land_compensation_p95_gtco2": lq["p95"],
            "ocean_compensation_p05_gtco2": oq["p05"],
            "ocean_compensation_p50_gtco2": oq["p50"],
            "ocean_compensation_p95_gtco2": oq["p95"],
            "budget_residual_p05_gtco2": cq["p05"],
            "budget_residual_p50_gtco2": cq["p50"],
            "budget_residual_p95_gtco2": cq["p95"],
        })
        non_on_q = qrow(fit_on["non"][i]); non_off_q = qrow(fit_off["non"][i])
        rows_ts.append({
            "year": int(y),
            "target_nonco2_wm2": float(target[i]),
            "nonco2_on_p05_wm2": non_on_q["p05"],
            "nonco2_on_p50_wm2": non_on_q["p50"],
            "nonco2_on_p95_wm2": non_on_q["p95"],
            "nonco2_off_p05_wm2": non_off_q["p05"],
            "nonco2_off_p50_wm2": non_off_q["p50"],
            "nonco2_off_p95_wm2": non_off_q["p95"],
            "co2_on_p05_ppm": onq["p05"],
            "co2_on_p50_ppm": onq["p50"],
            "co2_on_p95_ppm": onq["p95"],
            "co2_off_p05_ppm": offq["p05"],
            "co2_off_p50_ppm": offq["p50"],
            "co2_off_p95_ppm": offq["p95"],
        })

    pd.DataFrame(rows_attr).to_csv(outdir / "oscar_matched_attribution.csv", index=False)
    pd.DataFrame(rows_res).to_csv(outdir / "oscar_reservoir_decomposition.csv", index=False)
    pd.DataFrame(rows_ts).to_csv(outdir / "oscar_matched_timeseries.csv", index=False)
    pd.DataFrame(audit_on + audit_off).to_csv(outdir / "oscar_forcing_match_audit.csv", index=False)

    epf_gate = float(max(np.max(np.abs(epf_on)), np.max(np.abs(epf_off))))
    milestones = []
    attr_df = pd.DataFrame(rows_attr)
    res_df = pd.DataFrame(rows_res)
    for y in (2040, 2100, 2156, 2183):
        ar = attr_df.loc[attr_df.year == y].iloc[0].to_dict()
        rr = res_df.loc[res_df.year == y].iloc[0].to_dict()
        milestones.append({**ar, **{k: v for k, v in rr.items() if k not in {"year", "cumulative_cdr_gtco2"}}})

    gates = {
        "forcing_target_on_median_max_abs_wm2": float(fit_on["error"]) <= args.forcing_tol,
        "forcing_target_off_median_max_abs_wm2": float(fit_off["error"]) <= args.forcing_tol,
        "finite_co2": bool(np.isfinite(co2_on).all() and np.isfinite(co2_off).all()),
        "removal_lowers_2183_median_co2": bool(np.median(delta_ppm[-1]) > 0),
        "endogenous_permafrost_disabled": epf_gate < 1e-10,
    }
    passed = all(gates.values())

    final_resp = np.median(delta_ppm[-1]) * base.GTCO2_PER_PPM / cumulative_cdr[-1]
    summary = {
        "status": "PASS" if passed else "FAIL",
        "experiment": "V63_OSCAR_COMMON_NONCO2_MATCH_AND_RESERVOIR_DECOMPOSITION",
        "revision": MATCH_REVISION,
        "oscar_version": "v3.3",
        "oscar_pinned_commit": base.OSCAR_PIN,
        "ensemble_members": args.nmc,
        "seed": base.SEED,
        "years": [int(years[0]), int(years[-1])],
        "forcing_protocol": {
            "target_file": target_path.name,
            "target_definition": "FaIR calibrated-constrained medium-extension median non-CO2 ERF",
            "oscar_match_quantity": "ensemble median of RF_warm - RF_CO2",
            "synthetic_carrier": "RF_contr repurposed as additive bookkeeping forcing; not physical contrail forcing",
            "tolerance_wm2": args.forcing_tol,
            "on_max_abs_median_error_wm2": float(fit_on["error"]),
            "off_max_abs_median_error_wm2": float(fit_off["error"]),
        },
        "canonical_cdr_2027_2183_gtco2": float(cumulative_cdr[-1]),
        "final_2183_delta_co2_p50_ppm": float(np.median(delta_ppm[-1])),
        "final_2183_response_fraction_p50": float(final_resp),
        "reservoir_sign_convention": "positive compensation = OFF sink minus ON sink; CDR caused reduced uptake/backflux",
        "permafrost_note": "OSCAR endogenous permafrost is disabled because V62 supplies an exogenous permafrost source identically to ON/OFF.",
        "milestones": milestones,
        "gates": gates,
        "claim_boundary": (
            "This experiment matches the OSCAR ensemble-median non-CO2 warming forcing to the same target used by FaIR/Hector. "
            "It does not force every OSCAR ensemble member individually. Reservoir decomposition is diagnostic of OSCAR structure, not universal Earth-system attribution."
        ),
    }
    (outdir / "oscar_matched_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
