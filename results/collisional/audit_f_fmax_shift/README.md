# Tau Ceti (f/f_{\max}) reconciliation audit

This isolated workspace decomposes the shift between the accepted reference
collisional calculation and the archival Herschel–ALMA diagnostic calculation.
It does not modify production outputs.

Run from the repository root:

```bash
PYTHONPATH=. .venv/bin/python audit_f_fmax_shift/reconcile.py
PYTHONPATH=. .venv/bin/python -m unittest discover -s audit_f_fmax_shift/tests -v
```

The script reconstructs the reference pathway from
`results/sed/single_mbb_posterior.npz`, reconstructs the omitted diagnostic
collisional random variables from the diagnostic seed and source draw order,
and verifies the stored diagnostic ratios. It writes tables, figures, a PDF
report, manuscript-ready text and a run manifest under this directory.

The accepted reference result remains authoritative. The diagnostic result is
retained as a sensitivity analysis because its Herschel/ALMA posterior is
map-level/compatibility-runtime based and boundary-sensitive.
