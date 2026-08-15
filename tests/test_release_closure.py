import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_baseline_locked_values():
    payload = {r["quantity"]: float(r["regenerated"]) for r in csv.DictReader((ROOT / "tables/reproduction_comparison.csv").open())}
    assert math.isclose(payload["q_tau_reference_median"], 0.9274083133073248, abs_tol=5e-4)
    assert math.isclose(payload["z_tau_reference_median"], 1.456758216415282, abs_tol=2e-2)


def test_tau_ceti_holdout():
    holdout = json.loads((ROOT / "data/processed/tau_ceti_holdout.json").read_text())
    assert holdout["tau_ceti_in_fit"] is False


def test_collisional_reference():
    rows = list(csv.DictReader((ROOT / "tables/collisional_geometry_summary.csv").open()))
    by = {r["formulation"]: r for r in rows}
    assert math.isclose(float(by["A_characteristic_radius_ring"]["ratio_p50"]), 0.18608025936957256, rel_tol=1e-9)
    assert math.isclose(float(by["B_continuous_broad_belt"]["ratio_p50"]), 0.2291458892158028, rel_tol=1e-9)


def test_grain_metadata_preserves_distinction():
    text = (ROOT / "config/grain_scenarios.yaml").read_text()
    assert "D_bl:" in text and "D_min:" in text and "D_char:" in text


def test_posterior_archives_exist():
    for name in ("population_strict_merged_idata.nc", "population_dunes_only_idata.nc"):
        assert (ROOT / "posterior_samples" / name).stat().st_size > 100_000


if __name__ == "__main__":
    for name, test in globals().copy().items():
        if name.startswith("test_"):
            test()
    print("release tests passed")
