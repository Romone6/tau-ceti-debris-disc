# Change log from accepted baseline

## Preserved

- Accepted fractional luminosity posterior.
- Accepted strict merged DEBRIS--DUNES population result.
- Corrected collisional-geometry result.
- Treatment of the ALMA central source as stellar/chromospheric.

## Added in `raw-data-reduction-phase`

- Official archive observation inventories and checksums.
- Herschel and ALMA raw/product downloads.
- Extracted and staged archive product trees.
- CASA/HIPE runtime probes.
- Herschel map-level image diagnostics.
- Deduplicated PACS map-level forward fit with injection tests.
- Modular CASA compatibility execution of the delivered 12-m and ACA calibration
  scripts, calibrated measurement-set generation, and target visibility extraction.
- Diagnostic 12-m+ACA visibility-domain broad-belt fit and joint Herschel
  map/ALMA visibility posterior with propagated physical and population metrics.
- Formal baseline-to-reanalysis comparison table.
- New-observation proposal package.

## Not claimed

- No full HIPE telemetry re-reduction.
- No original CASA 4.2.2 binary execution; modular CASA 6.7.5 compatibility
  execution was used instead.
- No full-HIPE/raw-telemetry plus ALMA joint posterior.
- No new telescope observations.
- No replacement of the accepted \(f_\tau\), \(f/f_{\max}\), or \(q_\tau\)
  values based on the simplified map-level fit.

## Current interpretation

The new archive work materially improves provenance and provides a diagnostic
map/visibility revision. Its inferred fractional luminosity is lower and its
collisional-ratio tail is broader than the accepted baseline, while the shared
geometry remains broad and several-tens of AU in scale. Because the Herschel
telemetry was not reprocessed with HIPE and the ALMA scripts ran under a CASA
6.7.5 compatibility layer, the accepted baseline remains the authoritative
headline result.
