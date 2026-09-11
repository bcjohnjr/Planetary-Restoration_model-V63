# V63 predeclared Nature-candidate experiment protocol

## Purpose

Test whether the long-horizon separation in atmospheric CDR efficacy seen in the green V62.2 matched-forcing FaIR/Hector experiment is robust to CDR amplitude and timing, and identify the carbon-reservoir compensation responsible for the response in an independently structured model (OSCAR v3.3).

This protocol does **not** predeclare that the result is a universal threshold or law. Cross-model disagreement is retained and reported.

## Experiment A — CDR schedule generality

All non-CDR source terms are held exactly at the V62.2 canonical values. Only `cdr_gtco2` changes. All cases end at the canonical pathway boundary, 2183; no terminal-rate extension is permitted in this experiment.

Predeclared cases:

1. `canonical_peak15p2` — exact V62.2 canonical schedule.
2. `scaled_peak5p0` — canonical temporal shape scaled to a 5 GtCO2/yr peak.
3. `scaled_peak10p0` — canonical temporal shape scaled to a 10 GtCO2/yr peak.
4. `scaled_peak20p0` — canonical temporal shape scaled to a 20 GtCO2/yr peak.
5. `same_total_uniform` — canonical cumulative removal spread uniformly over 2027-2183.
6. `same_total_peak20_early` — canonical cumulative removal delivered as early as possible at no more than 20 GtCO2/yr.
7. `same_total_peak20_late` — canonical cumulative removal delivered as late as possible at no more than 20 GtCO2/yr.

The amplitude family tests dose/rate dependence. The fixed-total family tests timing dependence independently of total removal.

### Models

- FaIR 2.2.4, 841 calibrated/constrained configurations.
- Hector 3.5.0.

Future non-CO2 forcing is matched using the already validated V62 matched-forcing protocol. Structural model disagreement is not a failure gate.

### Primary schedule metric

For each case, report the first year in which Hector remains outside the FaIR p05-p95 atmospheric-response-fraction interval for at least five consecutive annual points. Also report 2040, 2100, 2156 and 2183 response fractions.

No special status is assigned to 2184 in advance; the schedule experiment ends at 2183.

## Experiment B — OSCAR common-forcing third-model test

Run OSCAR v3.3 at pinned commit `3ce008400e06363564e5981a35cbc32377d41d86`, 200 Monte Carlo configurations, seed 6201, under the canonical V62.2 carbon pathway through 2183.

The future non-CO2 target is the same FaIR-derived target used by the green FaIR/Hector matched experiment. OSCAR's `RF_contr` driver is repurposed only as a synthetic additive bookkeeping carrier. It is not interpreted as physical contrail forcing. The experiment iterates until the OSCAR ensemble-median `RF_warm - RF_CO2` matches the common target to <= 1e-4 W/m2 in both removal-ON and removal-OFF branches.

This is an ensemble-median forcing match for OSCAR; unlike FaIR, it does not claim every OSCAR member is individually forced to the identical target.

## Experiment C — carbon-reservoir decomposition

Within the matched OSCAR experiment, decompose the paired CDR response into:

- atmospheric CO2 benefit (`OFF - ON` atmospheric burden),
- land-sink compensation (`OFF D_Fland - ON D_Fland`, cumulatively integrated),
- ocean-sink compensation (`OFF D_Focean - ON D_Focean`, cumulatively integrated),
- carbon-budget residual.

Positive sink compensation means CDR causes reduced natural uptake or backflux relative to removal-OFF.

OSCAR endogenous permafrost is disabled because V62.2 already supplies an explicit exogenous permafrost CO2 source, identical in ON/OFF. Therefore this experiment does not claim an endogenous OSCAR permafrost attribution.

## Claim rules

1. A universal CDR-efficiency threshold is not claimed from one pathway or one model pair.
2. A mechanism claim requires schedule robustness plus interpretable reservoir compensation.
3. The climate experiments do not validate V62.2 ecological, engineering, finance, governance, health, food, legal or robotics modules.
4. Failed or contrary results are retained; experimental definitions are not changed after observing the outcomes without a new version and explicit reason.
