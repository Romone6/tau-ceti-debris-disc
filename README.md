# Tau Ceti debris-disc reproducibility release

Frozen scientific product for the Tau Ceti debris-disc analysis. This repository contains only the release configuration, processed inputs, derived tables, posterior archives, diagnostics, figures and reproduction code—not the full working directory or third-party telescope archives.

## Accepted headline values

| Quantity | Frozen value |
|---|---:|
| (f_\tau) median | (8.3574\times10^{-6}) |
| (q_\tau) median | 0.9274 |
| (z_\tau) median | 1.4568 |
| (P(q>0.90)) | 0.7100 |
| (P(q>0.95)) | 0.3167 |
| (f/f_{\max}), ring | 0.1861 |
| (f/f_{\max}), broad belt | 0.2291 |
| (D_{\rm bl}) scenario median | 0.95 μm |
| (D_{\rm bl}), 16–84% scenario range | 0.46–2.06 μm |

The PyMC files are validation sensitivities, not replacements for the accepted MLE/bootstrap headline. Geometry credible intervals are not claimed: the outer scale is robust near 50–60 AU, while the inner edge and position angle remain model-sensitive.

## Reproduce

```bash
conda env create -f environment.yml
conda activate tauceti
make verify
make reproduce
```

`make verify` is the fast frozen-value and archive-integrity check. `make reproduce` reruns the 200,000-draw collisional geometry analysis and refreshes its tables and figures. The full Herschel/ALMA archive reductions are intentionally not bundled; see [`data/README.md`](data/README.md).

## Release and citation

This is a pre-deposit release. No public DOI is claimed until the repository is connected to Zenodo and the record is published. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) and [`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md).
