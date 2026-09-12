# geosmooth synthetic registry

This directory is the package-owned registry for reusable synthetic statistical
datasets. It separates reusable recipes from frozen dataset instances.

## Tables

- `recipes.csv` binds each recipe ID to one geometry, sampling, truth, and
  response component and to the complete canonical specification SHA-256.
- `geometries.csv`, `samplings.csv`, `truths.csv`, and `responses.csv` are
  normalized component ledgers. Each row contains typed, inspectable JSON
  parameters, its stable ID, kind, version, and SHA-256. Matrix parameters
  use declared dimensions and row-major nested arrays. The public resolver
  reconstructs recipes from these rows rather than dispatching on recipe-name
  fragments.
- `instances.csv` records the sample size, seed, RNG policy, and expected
  content checksum for each frozen instance. Its `checksum.id` and component
  recipe foreign keys must each resolve exactly once.
- `checksums.csv` states the checksum payload, XDR serialization, scientific
  tolerances, and an exact foreign key to `environment_fingerprints.csv`.
- `environment_fingerprints.csv` records complete exact-environment scopes,
  including compiler, operating-system machine identity, explicitly
  established RNG version and kinds, absolute BLAS/LAPACK paths and file
  SHA-256 values, math-runtime identity, `geosmooth` version, and
  registry/evaluator and dependency versions.
- `environment_fingerprint_schema.csv` is the normative field contract used by
  both runtime scope matching and registry tests.
- `fixtures/` contains one committed canonical dataset object for every
  frozen instance. Exact content checksums are enforced when the current
  environment matches the recorded scope; otherwise the object is compared
  scientifically with its fixture under the recorded tolerances.
- `gaussian_mixture_1d_shapes.csv` contains the maintained S01--S16 scientific
  shape definitions.
- `gaussian_mixture_1d_variants.csv` contains the maintained V1--V3
  sampling/noise mechanisms.
- `legacy_family_disposition.csv` records implementation, follow-up, and
  retirement decisions for legacy families considered during the `gflow`
  split. A pending row is not evidence that its legacy code is safe to delete.

## Identity and regeneration

Ordinary materializations are identified by the full specification hash,
sample size, seed, and RNG policy in `dataset.id`. Frozen instance IDs are the
only shorter exception, and `materialize.synthetic.instance()` verifies their
complete registry rows and content. Ordinary API behavior has no
process-global verification bypass; the source-owned constructor dispatcher
is private to the regeneration script.

Regenerate normalized ledgers from the package source and catalog assets with:

```sh
Rscript scripts/update_synthetic_registry.R
```

Tests require every foreign key to resolve exactly once, every JSON component
and recipe hash to reproduce, every environment field in the normative schema
to be present, both exact and controlled scientific-fallback verification to
run, every fixture to validate, and every active frozen checksum to replay. Do
not edit generated normalized ledgers without also changing their owning
component or catalog source.

The legacy G1--G7 source was frozen from geosmooth revision
`5ecf721eae2765d8cde01d7fb82b17ff8bde8599`. The shared SSRHE source is
`trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R`;
its source SHA-256 is stored on the committed parity fixture.
