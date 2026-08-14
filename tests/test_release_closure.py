import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
F = ROOT / "finalisation"


def test_baseline_locked_values():
    payload = json.loads((ROOT / "results/final/accepted_merged_model_summary.json").read_text())
    assert math.isclose(payload["percentile"]["median"], 0.9274083133073248, abs_tol=5e-4)
    assert math.isclose(payload["z"]["median"], 1.456758216415282, abs_tol=2e-2)


def test_tau_ceti_holdout():
    holdout = json.loads((F / "data/canonical/tau_ceti_holdout.json").read_text())
    assert holdout["tau_ceti_in_fit"] is False


def test_collisional_matrix_reference_and_diagnostic():
    rows = list(csv.DictReader((F / "tables/controlled_model_matrix.csv").open()))
    by = {r["model"]: r for r in rows}
    assert math.isclose(float(by["A_reference_physical_Dbl"]["ring_ratio_median"]), 0.18740704196554983, rel_tol=1e-6)
    assert math.isclose(float(by["F_full_diagnostic_Dmin20"]["ring_ratio_median"]), 0.6639454235362406, rel_tol=1e-6)


def test_grain_metadata_preserves_distinction():
    text = (F / "config/grain_scenarios.yaml").read_text()
    assert "D_bl:" in text and "D_min:" in text and "D_char:" in text
    assert "not radiation-pressure D_bl" in (F / "config/collisional_dmin_20um.yaml").read_text()


def test_no_stale_fmax_six_claim_in_finalisation_text():
    paths = [F / "README.md", F / "final_report/reproducibility_closure_report.md", F / "paper/tau_ceti_raw_data_revision.tex"]
    for path in paths:
        if path.exists():
            assert "f/fmax approximately 6" not in path.read_text().lower()


def run_smoke_tests():
    tests = [test_baseline_locked_values, test_tau_ceti_holdout, test_collisional_matrix_reference_and_diagnostic, test_grain_metadata_preserves_distinction, test_no_stale_fmax_six_claim_in_finalisation_text]
    for test in tests: test()
    print(f"{len(tests)} closure smoke tests passed")


if __name__ == "__main__":
    run_smoke_tests()
