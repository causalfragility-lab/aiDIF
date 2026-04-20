## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on the check environment

* The `mirt` package is listed in Suggests but not available on the check
  system; it is an optional dependency used only for reading mirt output
  via `read_ai_scored()`.

* The Quarto TMPDIR warning is a Windows-only cosmetic issue in devtools
  and does not affect the package.

* "DIF" (Differential Item Functioning) and "DASB" (Differential AI Scoring
  Bias) are established and novel technical acronyms in psychometrics,
  not misspellings.

## Downstream dependencies

None (new submission).
