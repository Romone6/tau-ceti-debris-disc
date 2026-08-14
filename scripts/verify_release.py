"""Fast verification of the frozen public release."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

EXPECTED = {
    "f_tau": 8.357427655638217e-06,
    "q_tau": 0.9274083133073248,
    "z_tau": 1.456758216415282,
    "P_q_gt_090": 0.71,
    "P_q_gt_095": 0.31666666666666665,
    "ring_ratio": 0.18608025936957256,
    "broad_ratio": 0.2291458892158028,
}


def main() -> None:
    comparison = {row["quantity"]: float(row["regenerated"]) for row in csv.DictReader((ROOT / "tables/reproduction_comparison.csv").open())}
    mapping = {"f_tau": "f_tau_reference_median", "q_tau": "q_tau_reference_median", "z_tau": "z_tau_reference_median", "P_q_gt_090": "P_q_gt_090", "P_q_gt_095": "P_q_gt_095"}
    for key, quantity in mapping.items():
        assert abs(comparison[quantity] - EXPECTED[key]) < (5e-8 if key == "f_tau" else 5e-3), (key, comparison[quantity])
    collisional = {r["formulation"]: float(r["ratio_p50"]) for r in csv.DictReader((ROOT / "tables/collisional_geometry_summary.csv").open())}
    assert abs(collisional["A_characteristic_radius_ring"] - EXPECTED["ring_ratio"]) < 1e-9
    assert abs(collisional["B_continuous_broad_belt"] - EXPECTED["broad_ratio"]) < 1e-9
    grain = next(r for r in csv.DictReader((ROOT / "results/grains/blowout_summary.csv").open()) if r["scenario"] == "compact silicate" and r["Qpr"] == "1.0")
    assert float(grain["Dbl_median_um"]) > 0
    for path in [ROOT / "posterior_samples/population_strict_merged_idata.nc", ROOT / "posterior_samples/population_dunes_only_idata.nc"]:
        assert path.stat().st_size > 100_000
    manifest = json.loads((ROOT / "data/archive_manifest.json").read_text())
    assert manifest["verification"]["status"] == "ok"
    print("release verification passed")
    print(json.dumps(EXPECTED, indent=2))


if __name__ == "__main__":
    main()
