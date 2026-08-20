# Reproducibility closure report

## Decision

The accepted literature-constrained baseline is reproducible from the frozen checkpoint: f_tau=8.3574e-6, q_tau=0.9274, z_tau=1.4568, P(q>0.90)=0.710 and P(q>0.95)=0.3167. Reference physical collisional values remain 0.187 (ring) and 0.231 (broad); the controlled diagnostic matrix remains 0.062/0.077, 0.106/0.116, 0.432/0.472, and 0.664/0.726 for the specified scenarios.

## Gate status

| Gate | Status | Evidence |
|---|---|---|
| Accepted baseline | PASS | `final_report/accepted_baseline_verification.md`; manifest status `ok` |
| Collisional definitions/configuration | PASS | YAML configs, generated LaTeX prior table, and 200,000-draw production rerun archived in `diagnostics/collisional_full_rerun.log` and `tables/collisional_geometry_summary.csv` |
| Grain scenario traceability | PASS for analytic grid; HDF5/CSV archive | `results/grains`; parquet unavailable |
| Population four-chain MCMC diagnostics | PASS as validation sensitivity | Strict merged and DUNES-only ArviZ archives; R-hat=1.000, ESS>2,300, 0 divergences, BFMI>0.96 for strict merged |
| Tau Ceti holdout | PASS | `data/canonical/tau_ceti_holdout.json` |
| Geometry posterior/injection coverage | NOT SATISFIED | retained products are diagnostic map/visibility/Laplace outputs |
| Clean-environment reproduction | PARTIAL | existing project environment only; comparison table passes locked values |
| Public DOI | PASS | Zenodo record [10.5281/zenodo.22021388](https://doi.org/10.5281/zenodo.22021388) for the v1.0.2 release |
| Manuscript finalisation | BLOCKED | named revised manuscript file was not supplied |

## Scientific interpretation

The manuscript should remain pre-submission. The new PyMC population run is a converged validation sensitivity, not a replacement for the accepted MLE/bootstrap headline. No other diagnostic result is promoted. The broad outer scale is robust; inner edge and position angle remain model-sensitive. D_bl, D_min and D_char remain distinct, and the fixed 20 um value must remain a sensitivity assumption. The archival luminosity is not a replacement for the reference inference.

## Reproduction command

```bash
cd '/Users/romonedunlop/Documents/Tau Ceti Research-final-reproducibility-closure'
'/Users/romonedunlop/Documents/Tau Ceti Research/.venv/bin/python' run_baseline_verification.py
'/Users/romonedunlop/Documents/Tau Ceti Research/.venv/bin/python' finalisation/close_release.py
```

This is a closure package, not a claim that all requested specialist-submission gates passed.
