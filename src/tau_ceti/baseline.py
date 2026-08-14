"""Accepted-baseline verification helpers for the raw-data reduction phase."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ACCEPTED_BASELINE = {
    "fractional_luminosity_median": 8.36e-06,
    "fractional_luminosity_p16": 7.06e-06,
    "fractional_luminosity_p84": 1.03e-05,
    "q_tau_median": 0.927,
    "q_tau_p16": 0.886,
    "q_tau_p84": 0.966,
    "q_tau_p95_low": 0.854,
    "q_tau_p95_high": 0.983,
    "z_tau_median": 1.46,
    "z_tau_p16": 1.21,
    "z_tau_p84": 1.82,
    "p_q_gt_090": 0.710,
    "p_q_gt_095": 0.317,
    "p_q_gt_099": 0.000,
    "ring_ratio_median": 0.19,
    "broad_ratio_median": 0.23,
}

RANDOM_SEEDS = {
    "sed": 20260723,
    "population_catalogue_and_bootstrap": 20260723,
    "population_tau_predictive": 20260724,
}

CHECKSUM_PATHS = [
    "results/final/accepted_merged_model_summary.json",
    "results/tables/sed_model_comparison.csv",
    "results/tables/tau_ceti_merged_population_position.csv",
    "results/tables/collisional_geometry_summary.csv",
    "report/sed_analysis.md",
    "report/collisional_geometry_audit.md",
    "report/strict_merged_population_analysis.md",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def collect_baseline_metrics(root: Path) -> dict[str, Any]:
    """Collect accepted baseline metrics from the current repository outputs."""
    root = Path(root)

    sed_rows = _read_csv_rows(root / "results/tables/sed_model_comparison.csv")
    sed_fractional = next(
        row for row in sed_rows if row["model"] == "single_mbb" and row["quantity"] == "fractional_luminosity"
    )

    collisional_rows = _read_csv_rows(root / "results/tables/collisional_geometry_summary.csv")
    collisional = {
        row["formulation"]: {
            "ratio_p16": float(row["ratio_p16"]),
            "ratio_p50": float(row["ratio_p50"]),
            "ratio_p84": float(row["ratio_p84"]),
            "p_ratio_gt_1": float(row["p_ratio_gt_1"]),
            "p_ratio_gt_10": float(row["p_ratio_gt_10"]),
            "classification": row["classification"],
        }
        for row in collisional_rows
    }

    population_summary = json.loads((root / "results/final/accepted_merged_model_summary.json").read_text(encoding="utf-8"))

    return {
        "fractional_luminosity": {
            "p16": float(sed_fractional["p16"]),
            "median": float(sed_fractional["p50"]),
            "p84": float(sed_fractional["p84"]),
        },
        "population": {
            "q_tau": {
                "p16": float(population_summary["percentile"]["ci68"][0]),
                "median": float(population_summary["percentile"]["median"]),
                "p84": float(population_summary["percentile"]["ci68"][1]),
                "p95_low": float(population_summary["percentile"]["ci95"][0]),
                "p95_high": float(population_summary["percentile"]["ci95"][1]),
                "p_q_gt_090": float(population_summary["percentile"]["P_q_gt_090"]),
                "p_q_gt_095": float(population_summary["percentile"]["P_q_gt_095"]),
                "p_q_gt_099": float(population_summary["percentile"]["P_q_gt_099"]),
            },
            "z_tau": {
                "p16": float(population_summary["z"]["ci68"][0]),
                "median": float(population_summary["z"]["median"]),
                "p84": float(population_summary["z"]["ci68"][1]),
            },
            "sample": population_summary["sample"],
        },
        "collisional": collisional,
    }


def verify_accepted_baseline(metrics: dict[str, Any]) -> dict[str, Any]:
    """Verify collected metrics against the accepted baseline values."""
    mismatches: list[str] = []

    def check(name: str, observed: float, expected: float, tolerance: float) -> None:
        if abs(observed - expected) > tolerance:
            mismatches.append(f"{name}: observed={observed} expected≈{expected} tol={tolerance}")

    check("fractional_luminosity_median", metrics["fractional_luminosity"]["median"], ACCEPTED_BASELINE["fractional_luminosity_median"], 5e-08)
    check("fractional_luminosity_p16", metrics["fractional_luminosity"]["p16"], ACCEPTED_BASELINE["fractional_luminosity_p16"], 5e-08)
    check("fractional_luminosity_p84", metrics["fractional_luminosity"]["p84"], ACCEPTED_BASELINE["fractional_luminosity_p84"], 5e-07)

    q_tau = metrics["population"]["q_tau"]
    z_tau = metrics["population"]["z_tau"]
    check("q_tau_median", q_tau["median"], ACCEPTED_BASELINE["q_tau_median"], 5e-04)
    check("q_tau_p16", q_tau["p16"], ACCEPTED_BASELINE["q_tau_p16"], 5e-04)
    check("q_tau_p84", q_tau["p84"], ACCEPTED_BASELINE["q_tau_p84"], 5e-04)
    check("q_tau_p95_low", q_tau["p95_low"], ACCEPTED_BASELINE["q_tau_p95_low"], 5e-04)
    check("q_tau_p95_high", q_tau["p95_high"], ACCEPTED_BASELINE["q_tau_p95_high"], 5e-04)
    check("p_q_gt_090", q_tau["p_q_gt_090"], ACCEPTED_BASELINE["p_q_gt_090"], 5e-03)
    check("p_q_gt_095", q_tau["p_q_gt_095"], ACCEPTED_BASELINE["p_q_gt_095"], 5e-03)
    check("p_q_gt_099", q_tau["p_q_gt_099"], ACCEPTED_BASELINE["p_q_gt_099"], 1e-12)
    check("z_tau_median", z_tau["median"], ACCEPTED_BASELINE["z_tau_median"], 2e-02)
    check("z_tau_p16", z_tau["p16"], ACCEPTED_BASELINE["z_tau_p16"], 2e-02)
    check("z_tau_p84", z_tau["p84"], ACCEPTED_BASELINE["z_tau_p84"], 2e-02)

    check(
        "ring_ratio_median",
        metrics["collisional"]["A_characteristic_radius_ring"]["ratio_p50"],
        ACCEPTED_BASELINE["ring_ratio_median"],
        2e-02,
    )
    check(
        "broad_ratio_median",
        metrics["collisional"]["B_continuous_broad_belt"]["ratio_p50"],
        ACCEPTED_BASELINE["broad_ratio_median"],
        2e-02,
    )

    return {
        "status": "ok" if not mismatches else "mismatch",
        "mismatches": mismatches,
        "accepted_baseline": ACCEPTED_BASELINE.copy(),
    }


def build_baseline_manifest(
    root: Path,
    *,
    metrics: dict[str, Any],
    verification: dict[str, Any],
    git_commit: str,
    python_version: str,
    package_versions: dict[str, str],
) -> dict[str, Any]:
    """Build a read-only manifest for the accepted baseline state."""
    root = Path(root)
    checksums = {path: _sha256(root / path) for path in CHECKSUM_PATHS}

    return {
        "git_commit": git_commit,
        "python_version": python_version,
        "package_versions": dict(package_versions),
        "random_seeds": RANDOM_SEEDS.copy(),
        "verification": verification,
        "metrics": metrics,
        "checksums": checksums,
    }


def write_baseline_artifacts(manifest_path: Path, report_path: Path, manifest: dict[str, Any]) -> None:
    """Write the accepted baseline manifest and a short verification report."""
    manifest_path = Path(manifest_path)
    report_path = Path(report_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    verification = manifest["verification"]
    q_tau = manifest["metrics"]["population"]["q_tau"]
    z_tau = manifest["metrics"]["population"]["z_tau"]
    f_tau = manifest["metrics"]["fractional_luminosity"]
    lines = [
        "# Accepted baseline verification",
        "",
        f"- Status: `{verification['status']}`",
        f"- Git commit: `{manifest['git_commit']}`",
        f"- Python: `{manifest['python_version']}`",
        "",
        "## Verified baseline values",
        "",
        f"- f_tau median: `{f_tau['median']:.12g}`",
        f"- q_tau median: `{q_tau['median']:.12g}`",
        f"- z_tau median: `{z_tau['median']:.12g}`",
        f"- P(q_tau > 0.90): `{q_tau['p_q_gt_090']:.12g}`",
        f"- P(q_tau > 0.95): `{q_tau['p_q_gt_095']:.12g}`",
        f"- P(q_tau > 0.99): `{q_tau['p_q_gt_099']:.12g}`",
        "",
    ]
    if verification["mismatches"]:
        lines.extend(["## Mismatches", ""])
        lines.extend(f"- {item}" for item in verification["mismatches"])
    else:
        lines.extend(["No mismatches detected.", ""])

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
