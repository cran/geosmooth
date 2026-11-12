# geosmooth 0.2.1

* Reduce package-check runtime by sharing test setup and grouping reference
  comparisons, retaining all reference cases and numerical tolerances.
* Declare the `stats::predict()` and `utils::tail()` imports used by
  trend-filtering helpers.
* Add `fit.projection.trend.filtering()` and
  `cv.projection.trend.filtering()` for sums of order-2 trend-filtered curves
  on fixed or learned projections, with Gaussian and binomial responses.
* Refresh function documentation and vignettes. Include the stored PCA/PLS
  comparison results and their reproduction script.
* Reject exactly constant covariate columns before weighted centering in
  Gaussian and logistic trend filtering, independently of BLAS rounding.

# geosmooth 0.2.0

## Documentation and checks

* The PCA/PLS vignette now loads stored simulation results instead of
  repeating the simulations each time it is built. The complete results and
  the code that reproduces them are included with the package.
* The tests run faster. The synthetic-data tests share their setup, and the
  tests of graph cross-validation, binomial fits, and constant responses use
  smaller examples. Every test case and every numerical tolerance is kept.

## Trend filtering and graph methods

* Gaussian and logistic trend filtering now reject a covariate column that
  is exactly constant, and they check this before the column is centered.
  The check no longer depends on rounding in the linear algebra library, so
  a constant column cannot be mistaken for a varying one. Covariates with
  small but nonzero variation are still accepted.

* GCV for graph heat smoothing with `numerical.policy = "verified.fallback"`
  is now computed accurately on weakly connected graphs. The Schur complement
  over the observed vertices is formed as the crossproduct of a residual
  factor, which avoids subtracting two nearly equal matrices. The model, the
  independent reference calculation, and the acceptance tolerances are
  unchanged, and a smoothing time is still not selected when the GCV scores
  cannot distinguish the candidates.

* Fixed cases in which logistic trend filtering stopped making progress
  before converging, in both the banded and the dense solver. The banded
  solver now projects the dual variables of inactive knots in the same
  centered coordinates that its convergence check uses. The dense solver
  evaluates ridge penalties from their residual factors and computes the
  change of the objective in the line search directly. The objective and the
  convergence tolerances are unchanged, also for the Gaussian solver, which
  shares this code. Tests cover the 17 recorded failures and the penalty
  paths around them, with independent checks of optimality.

* Corrected how the penalty is scaled in cross-validation for trend
  filtering. Local polynomial lifting, graph, and one-dimensional Gaussian
  and logistic fits now use `cv.penalty.scale = "training_mass"` by default:
  in each fold, lambda (and, for the one-dimensional fits, any explicit ridge
  penalty) is multiplied by the ratio of the training weight to the full
  weight. `"none"` keeps the earlier behavior of using the same numerical
  penalty in every fold. The penalties used in each fold are recorded. Fits
  with a fixed lambda are not affected. The assessment functions use the
  same convention in their outer and inner folds. Fitted objects saved by
  earlier versions keep the earlier convention. LPL-TF now selects a penalty
  only among candidates that succeeded in every fold.

* Add `fit.gaussian.trend.filtering()` and
  `assess.gaussian.trend.filtering()` for continuous outcomes with order-2
  trend filtering. Native dense and banded solvers use squared-error loss;
  the assessment retunes on every permutation, retains failed replicates,
  reports nested-CV MSE/MAE improvement, and supports optional absolute-deviation
  and total-variation tests with Holm adjustment. Randomization requires
  exchangeable outcomes under a common-distribution null, not just equal means
  with arbitrary heteroskedasticity. Existing logistic behavior is preserved.

* Extend `assess.logistic.trend.filtering()` with participant-averaged absolute
  probability deviation and total probability variation, including interior
  turning points. Keep likelihood as the default test; optional shape tests
  reuse the same retuned permutation fits and report individual and Holm-adjusted
  p-values, complete diagnostics and failure bounds.

* Stabilize logistic trend-filtering fits with clustered active knots using a
  scaled banded QR projection instead of normal equations. Use compensated
  objective sums in the dense solver and shared quadratic initialization to
  prevent roundoff-induced line-search stalls. Preserve the objective and
  convergence tolerances.

* Add `assess.logistic.trend.filtering()` for a global permutation test of
  constant binary probability and separate nested-CV log-loss/Brier improvement.
  It repeats tuning, uses outcome-independent folds, retains failed replicates
  and restricts shuffling to independent unit-weight binary observations without
  covariates or offsets.

* Add a reproducible logistic trend-filtering vignette with synthetic truth,
  training-only penalty selection, binomial trial weights and held-out Pima
  classification/calibration examples.

* Add streamed intrinsic surface distances to all three neighbor probability
  estimators with `metric = "linf.simplex"`. Reuse the dgraphs C++ solver through
  LinkingTo, cache sparse training profiles per prediction/CV split, and retain
  only one query neighborhood at a time. Preserve participant exclusions and
  tie rules; reject ANN and positive epsilon for intrinsic distances.

* Add a compact, banded C++ backend for order-2 logistic trend filtering.
  Automatic selection avoids full dense knot matrices above 256 locations or
  when the dense workspace limit is exceeded. Full-knot KKT checks preserve the
  objective without thinning; prediction also uses the compact representation.

* Add exact Euclidean ANN search and opt-in approximation (`eps > 0`) to all
  three neighbor probability estimators. `search = "auto"` selects exact ANN
  for Euclidean coordinates and brute force for other inputs. Preserve grouped
  exclusions and nearest answering visits; reuse trees within prediction and CV.

* Add `fit.logistic.trend.filtering()` for order-2 univariate logistic trend
  filtering with a C++ backend, exact ties, binomial weights, optional covariates
  and offsets, explicit ridge options, fold-based penalty selection, and
  stationarity/complementarity diagnostics. Failed fits remain explicit.
* Add `fit.neighbor.probability()` with C++ ordinary kNN, raw nearest-case
  and stopping-rule corrected estimators. Shared prediction paths use
  participant exclusions, closest answering visits, deterministic ID ties,
  missing-label handling and grouped Brier selection. Validate against the
  original sPTB Python prediction and minimum-error tuning code.

* Add opt-in intrinsic simplex and Euclidean Bernoulli charts to
  `fit.local.likelihood()`, with weighted PCA/PLS, training-fold-only fitting,
  stable ridge optimization and explicit diagnostics/fallbacks. Legacy defaults
  and density fits are unchanged.
* Add `unfold.simplex.chart()` for conservative one/two-face certificates.
  Corner, tie and multi-face supports are not forced into a chart.
* Graph low-pass operators, bases and fits accept `weight.type = "conductance"`
  for direct symmetric nonnegative affinities. Zero-weight edges are removed;
  positive conductances are not transformed or clipped. Length remains the default.
* `spectrum = "full"` explicitly obtains and verifies all eigenpairs. Auto uses
  dense decomposition for requests exceeding the sparse backend capability;
  incompatible sparse-only requests now fail instead of silently returning n-2.
* Low-pass GCV now uses the conventional observed-count normalization. Exact
  zero heat time is an identity diagnostic with undefined GCV, excluded from
  selection. Positive-time heat scores use stable spectral complements.
* Full-spectrum heat fits, paths and refits accept `NA` responses, fitting
  observed squared error plus the heat penalty `exp(t L) - I`. GCV uses the
  observed fitting operator and observed count. A wholly unobserved component
  is unidentified and fails explicitly. This extension does not model assay
  selection, repeated-measurement covariance or observation-frequency weights.

* Add opt-in verified observed-heat fallback with independently checked score
  formulations and explicit numerical diagnostics. Automatic tuning refuses
  structurally flat paths and unresolved competing minima; prescribed times
  remain available. Numerical resolution bands are not statistical intervals.

## Solver reliability

* Hessian L1 ADMM uses an operator-scale initial penalty by default
  (`admm.rho = NULL`) and adapts it during its initial iterations. An explicit
  positive `admm.rho` still sets the initial value; use
  `admm.adaptive.rho = FALSE` to hold it fixed. Numerical stabilization now
  preserves the requested objective instead of adding an unintended ridge.
* ADMM convergence also checks stationarity of the original objective.
  Incomplete fixed fits warn and print their status. CV excludes candidates
  with any unsuccessful fold and errors if none qualify or the selected
  full-data ADMM refit is incomplete. Fold diagnostics are retained, including
  on the no-eligible-candidate error. Generalized-lasso warnings are retained.
* The shared PTTF L1 solver receives the same safeguards. Hessian L1 objective
  reporting now accounts for the penalty actually used with row scaling;
  the unscaled Hessian norm remains available separately.

* Generalized-lasso Hessian paths now propose penalties only: all returned
  coefficients and CV candidates are solved by ADMM on the original objective.
  Finite path coefficients can be suboptimal for nearly rank-deficient operators;
  missing-label path coefficients can include an extra ridge. The path remains
  diagnostic, with its original objectives in `solver$path.objective` and the
  refined route labeled `genlasso_admm_refined`. Default Hessian ADMM stopping
  tolerances are now `1e-6` absolute and `1e-5` relative. Fitted values and selected
  penalties can change; difficult candidates may require a larger iteration cap.

## Uncertainty and result interpretation

* `lps.pointwise.band()` now estimates noise variance using residual noise
  degrees of freedom `sum((I - S)^2) = n - 2*tr(S) + sum(S^2)`. The former
  `n - tr(S)` denominator is available as `variance.method = "legacy"`.
  `df` still reports `tr(S)`; `residual.df` is returned separately. Bands remain
  pointwise, conditional on a fixed smoother, and uncorrected for smoothing bias.
* Regression families support `fitted()` and `residuals()` with their existing
  vector/matrix and missing-label contracts. Residuals at other evaluation
  coordinates are rejected. Tracked harmonic fitted values use `fitted()`; density mass
  stays `$rho`. No new out-of-sample prediction or reusable-refit contract is implied.
* Lifting, synchronized lifting, density, prediction-synchronized, chart-kernel,
  and local-likelihood fits print bounded diagnostics. Structured `summary()`
  methods cover the main regression families. Unavailable convergence information
  is labeled explicitly, and density accounting and local fallbacks remain visible.
* The existing guides now include interpreted first-fit and graph workflows,
  embedded figures with text alternatives, and clear latent/ambient synthetic
  views. `plot.synthetic_dataset()` adds named `view`, `legend`, and `col`
  controls while preserving automatic view selection. Point-color keys now
  describe the displayed values; region coloring requires region labels.
* Graph adaptation uses dgraphs' public distance API and accepts both transitional
  graph lists and the forthcoming `dgraph` representation, without restoring
  retired geometry re-exports.

## Input validation and harmonic diagnostics

* Density methods now reject misspelled, duplicated, or unsupported controls
  and unused terminal arguments. Logical controls require `TRUE` or `FALSE`.
  Recognized controls forwarded to subject-occupation smoothers remain supported.
* Graph low-pass, density, and harmonic methods validate finite integer vertex
  indices before conversion. Fractional indices and iteration counts now error
  instead of silently selecting a different vertex or truncating a control.
* Harmonic smoothing documents the actual boundary rule: a region vertex must
  have a neighbor outside the region; graph leaves are not automatically fixed.
  The basic function still leaves a boundary-free region unchanged, while the
  tracked function still relaxes it. Both now report explicit numerical status.
* `harmonic.smoother()` reports extrema-stability detection separately from
  numerical convergence. `stable_iteration` is now `NA_integer_` when no
  stability window was detected; check `stability_detected` before using it.
  It returns update/residual histories and actual recorded iteration numbers,
  always records the final state, and uses those numbers in plots. Stability
  of recorded extrema does not establish convergence or optimal smoothing.

## API changes

This release has 71 explicit public exports, compared with 93 in version 0.1.0.
The 35 retired names are replaced by dependency-qualified calls, a unified
fitter, or two S3 generics.

* Replace the six exported `refit.*()` functions with the S3 generic
  `refit(object, y, ...)`, covering MALPS, both lifting trend filters,
  metric graph low-pass, quadratic Hessian, and Hessian L1 fits. In old
  calls, rename `fitted.model` to `object` and `y.new` to `y`.
  Method-specific controls, response-shape restrictions, selection behavior,
  and return classes are preserved. Former help names remain available.
* Replace `lps.smoother.matrix()` and `malps.smoother.matrix()` with
  `smoother.matrix(object, ...)`. LPS retains `check.tol` and its
  fixed-configuration requirements; MALPS retains `max.n`, `allow.robust`,
  and its conditional interpretation. The eight former entry points are
  no longer exported. The function-guide vignette documents all supported methods.

* API consolidation (breaking change): retire the 25 transitional synthetic
  geometry and sampling re-exports. Call the same names and arguments with
  `dgraphs::`, for example `dgraphs::synthetic.circle()` and
  `dgraphs::synthetic.sampling.uniform.interval(0, 1)`. Statistical recipe,
  truth, response, materialization, and validation functions remain in
  geosmooth, including its local `synthetic.sampling.stratified()`.
* Consolidate quadratic Hessian regression into
  `fit.ssrhe.hessian.regression()`. Existing fixed calls with `lambda1` and
  `lambda2` are unchanged. Replace the former `.cv()` and `.gcv()` fitters
  with `lambda.selection = "cv"` or `"gcv"` and penalty grids. The old
  function names are no longer exported; their help aliases explain migration.
* CV settings now use `cv.control = list(foldid = ..., cv.folds = ...,
  loss = "mse", selection = "min")`; `foldid` replaces `fold.id` and
  `cv.folds` replaces `nfolds`. GCV trace settings use
  `gcv.control = list(trace.method = "hutchinson", trace.n.probes = 50L,
  trace.seed = ...)`. Both control lists support `support.selection`,
  `support.grid`, and `support.max.candidates` (formerly the mode-specific
  `support.cv.max.candidates` / `support.gcv.max.candidates`). Nonempty controls
  for the wrong mode, unknown control entries, and simultaneous scalar/grid
  penalties are errors. CV still supports missing labels; GCV still requires
  a single fully observed response and strictly positive weights. Existing
  fit classes, CV/GCV tables, and refit support are preserved.

## Dependencies

* Require R >= 4.1.0, matching the minimum required by dgraphs.
* Require dgraphs >= 0.4.0 in both Imports and LinkingTo for synthetic geometry,
  point generation, intrinsic simplex distances and their C++ header interface.
  Statistical recipes, responses, dataset identities and legacy G4 remain
  here. Geometry Lab and maintained geometry reference tools now live in
  the dgraphs development tree.
* Require grip >= 0.3.0 for the optional GRIP and edge-KK layout methods,
  which use the dot-delimited argument interface. Grip remains in Suggests.

## Documentation

* Clarify that metric graph low-pass refits process response columns
  sequentially. The retained `n.cores` argument does not enable parallel
  execution, and per-column GCV results report `n.cores.used = 1L`.
* Add "Finding your way around geosmooth", a task-oriented vignette covering
  every exported function, supported prediction/refit workflows, and migration
  to the consolidated API. Link the guide from the package overview and README.
* Update the synthetic-datasets vignette to use geometry and sampling from
  dgraphs. Both vignettes are distributed as HTML with the built package.

# geosmooth 0.1.0

* Initial CRAN release.
* Provides local polynomial and model-averaged local polynomial smoothing,
  local polynomial lifting trend filtering, synchronized local polynomial
  lifting trend filtering, and occupation-density adapters.
* Provides metric-graph low-pass filtering and Hessian-energy regression with
  compiled backends.
* Provides reproducible synthetic geometry and response generators for method
  evaluation.
* Makes synthetic-registry identity verification stable across supported R
  versions by excluding the serializer's producer-version header from the
  checksum contract.
