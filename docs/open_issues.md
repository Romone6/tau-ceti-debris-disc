# Open closure issues

1. The named author manuscript `Tau_Ceti_MNRAS_PreCodex_Revised.tex` was not present in the supplied attachments or project tree; only `paper/tau_ceti_raw_data_revision.tex` was available.
2. The accepted headline population result remains MLE plus bootstrap; separate four-chain PyMC strict-merged and DUNES-only validation sensitivities now pass R-hat/ESS/divergence/BFMI checks and are archived, but they are not silently promoted because explicit priors shift q_tau slightly.
3. The 200,000-draw collisional geometry audit has now been rerun from the archived production script and its exact prior/configuration table is frozen; the lightweight module CLI remains an archived-output validator rather than a second implementation.
4. No full Herschel/ALMA sampled joint posterior or injection-recovery coverage archive is present; geometry remains diagnostic.
5. PyArrow/pandas are unavailable, so canonical inputs and grain samples are archived as CSV/HDF5 with explicit `.parquet.MISSING` notices.
6. The v1.0.2 release is archived at Zenodo under DOI `10.5281/zenodo.22021388`; future changes must use a new versioned release rather than modifying the frozen record.
7. A clean isolated dependency environment was not built; the reproducibility command was run using the existing project environment.
