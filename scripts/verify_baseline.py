#!/usr/bin/env python3
"""Freeze the accepted baseline into a read-only manifest and verification note."""

from __future__ import annotations

import platform
import subprocess
from pathlib import Path

import matplotlib
import numpy
import scipy

from tau_ceti.baseline import (
    build_baseline_manifest,
    collect_baseline_metrics,
    verify_accepted_baseline,
    write_baseline_artifacts,
)


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "results/baseline/accepted_baseline_manifest.json"
REPORT = ROOT / "report/accepted_baseline_verification.md"


def git_commit(root: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        cwd=root,
    ).stdout.strip()


def main() -> None:
    metrics = collect_baseline_metrics(ROOT)
    verification = verify_accepted_baseline(metrics)
    manifest = build_baseline_manifest(
        ROOT,
        metrics=metrics,
        verification=verification,
        git_commit=git_commit(ROOT),
        python_version=platform.python_version(),
        package_versions={
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    )
    write_baseline_artifacts(MANIFEST, REPORT, manifest)


if __name__ == "__main__":
    main()
