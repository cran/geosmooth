# DGP-library frozen registry -- manifest

Frozen: 2026-06-11T07:37:11Z
Contract: LPS Tiers 1-4, Amendment 1 (consolidate the DGP library).
Plan (DGP definitions matched): lps_experimental_plan_2026-06-09.tex (sec:dgp).

## Provenance (source-tree state at freeze time)

- geosmooth git head: `10e6bafe35d418c238691e66c331fa79606acd94`
- branch: `codex/geosmooth-dgp-library`
- source tree clean at freeze: **true**
- R: R version 4.5.2 (2025-10-31)
- platform: aarch64-apple-darwin20
- BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
- LAPACK library: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib
- LAPACK version: 3.12.1
- digest: 0.6.39
- generator source: `R/dgp_library.R`  sha256 `2cff1220e5910b9f337264c8f01dca6f8f50bf691fc662ba31df9c0b4034156b`
- registry CSV: `inst/dgp_registry/dgp_registry.csv`  sha256 `da8ffa5226fcd5c9a7fffc4b52a1c1a6cac5d5058735f01037b27697e17236eb`

## Registry

- rows: 24
- G-tags covered: G1, G2, G3a, G3b, G3c, G3d, G4, G5, G6, G7

## Checksum definition (so the auditor can recompute each row)

Each row's `sha256` is `geosmooth::dgp.content.sha256(ds)` =
`digest::digest(.dgp.content(ds), algo = "sha256", serialize = TRUE)`, where
`.dgp.content(ds)` is the fixed-order list
`(dataset.id, gtag, n, p, d, seed, sigma, U, X, truth, y, region)` of the
materialized object -- environment-derived provenance (package version) is
excluded, so the checksum depends only on the data. Reproduce a row with:

```r
library(geosmooth)
ds <- dgp.materialize(<gtag>, <args from params column>)
stopifnot(dgp.content.sha256(ds) == <sha256 from CSV>)
```

Determinism is environment-relative (spec sec:rng): same R/BLAS/digest as
above -> bitwise-identical objects and checksums.

## Generating command

```sh
Rscript scripts/freeze_dgp_registry.R
```

Full environment in `inst/dgp_registry/sessionInfo.txt`.
