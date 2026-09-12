# Maintained synthetic recipe and frozen-instance registries.

.synthetic.registry.asset <- function(file) {
  dev <- file.path("inst", "synthetic_registry", file)
  if (file.exists(dev)) return(dev)
  system.file("synthetic_registry", file, package = "geosmooth")
}

.synthetic.registry.row <- function(file, key, value) {
  path <- .synthetic.registry.asset(file)
  if (!nzchar(path) || !file.exists(path)) {
    stop("Synthetic registry table is unavailable: ", file, call. = FALSE)
  }
  table <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!key %in% names(table)) {
    stop("Synthetic registry table ", file, " lacks key ", key, ".",
         call. = FALSE)
  }
  row <- table[table[[key]] == value, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(
      "Synthetic registry key must resolve exactly once: ",
      file, "$", key, " = ", value, call. = FALSE)
  }
  row
}

.synthetic.registry.json.node <- function(value) {
  if (is.null(value)) return(list(type = "null"))
  if (is.factor(value)) value <- as.character(value)
  if (is.matrix(value)) {
    type <- typeof(value)
    if (!type %in% c("integer", "double", "logical", "character")) {
      stop("Registry JSON does not support this matrix type.", call. = FALSE)
    }
    rows <- lapply(seq_len(nrow(value)), function(i) {
      structure(unname(value[i, , drop = TRUE]), class = "AsIs")
    })
    dim.labels <- dimnames(value)
    return(list(
      type = paste0(type, ".matrix"),
      nrow = as.integer(nrow(value)),
      ncol = as.integer(ncol(value)),
      values = rows,
      row.names = if (is.null(dim.labels[[1L]])) NULL else
        structure(enc2utf8(dim.labels[[1L]]), class = "AsIs"),
      column.names = if (is.null(dim.labels[[2L]])) NULL else
        structure(enc2utf8(dim.labels[[2L]]), class = "AsIs")
    ))
  }
  if (is.array(value)) {
    stop("Registry JSON supports matrices but not higher-order arrays.",
         call. = FALSE)
  }
  if (is.list(value)) {
    return(list(
      type = "list",
      names = if (is.null(names(value))) NULL else
        structure(enc2utf8(names(value)), class = "AsIs"),
      values = lapply(value, .synthetic.registry.json.node)
    ))
  }
  type <- typeof(value)
  if (!type %in% c("integer", "double", "logical", "character")) {
    stop("Registry JSON does not support this atomic type.", call. = FALSE)
  }
  if (type == "character") value <- enc2utf8(value)
  list(
    type = type,
    names = if (is.null(names(value))) NULL else
      structure(enc2utf8(names(value)), class = "AsIs"),
    values = structure(unname(value), class = "AsIs")
  )
}

.synthetic.registry.json.text <- function(value) {
  as.character(jsonlite::toJSON(
    .synthetic.registry.json.node(value),
    auto_unbox = TRUE, null = "null", na = "null", digits = 17L,
    pretty = FALSE))
}

.synthetic.registry.json.atomic <- function(values, type) {
  convert <- switch(
    type,
    integer = function(value) {
      if (is.null(value)) NA_integer_ else as.integer(value)
    },
    double = function(value) {
      if (is.null(value)) NA_real_ else as.double(value)
    },
    logical = function(value) {
      if (is.null(value)) NA else as.logical(value)
    },
    character = function(value) {
      if (is.null(value)) NA_character_ else enc2utf8(as.character(value))
    },
    stop("Unknown registry JSON atomic type: ", type, call. = FALSE))
  prototype <- switch(
    type, integer = integer(1L), double = double(1L),
    logical = logical(1L), character = character(1L))
  vapply(values, convert, prototype)
}

.synthetic.registry.json.decode.node <- function(node) {
  if (!is.list(node) || !is.character(node$type) ||
      length(node$type) != 1L) {
    stop("Registry JSON node is malformed.", call. = FALSE)
  }
  type <- node$type
  if (identical(type, "null")) return(NULL)
  if (identical(type, "list")) {
    value <- lapply(node$values, .synthetic.registry.json.decode.node)
    if (!is.null(node$names)) {
      names(value) <- enc2utf8(unlist(node$names, use.names = FALSE))
    }
    return(value)
  }
  if (grepl("\\.matrix$", type)) {
    atomic.type <- sub("\\.matrix$", "", type)
    rows <- lapply(node$values, .synthetic.registry.json.atomic,
                   type = atomic.type)
    flat <- if (length(rows)) unlist(rows, use.names = FALSE) else
      switch(
        atomic.type, integer = integer(), double = double(),
        logical = logical(), character = character())
    value <- matrix(
      flat, nrow = as.integer(node$nrow), ncol = as.integer(node$ncol),
      byrow = TRUE)
    row.names <- if (is.null(node$row.names)) NULL else
      enc2utf8(unlist(node$row.names, use.names = FALSE))
    column.names <- if (is.null(node$column.names)) NULL else
      enc2utf8(unlist(node$column.names, use.names = FALSE))
    if (!is.null(row.names) || !is.null(column.names)) {
      dimnames(value) <- list(row.names, column.names)
    }
    return(value)
  }
  value <- .synthetic.registry.json.atomic(node$values, type)
  if (!is.null(node$names)) {
    names(value) <- enc2utf8(unlist(node$names, use.names = FALSE))
  }
  value
}

.synthetic.registry.json.value <- function(value) {
  value <- .synthetic.scalar.character(value, "registry JSON")
  node <- tryCatch(
    jsonlite::fromJSON(value, simplifyVector = FALSE),
    error = function(error) {
      stop("Registry JSON is malformed: ", conditionMessage(error),
           call. = FALSE)
    })
  .synthetic.registry.json.decode.node(node)
}

.synthetic.registry.component.subclass <- function(kind, family) {
  if (identical(kind, "geometry") &&
      identical(family, "g4.segment.rectangle")) {
    return("synthetic_stratified_geometry")
  }
  paste0(
    "synthetic_", gsub(".", "_", family, fixed = TRUE), "_", kind)
}

.synthetic.registry.component <- function(component.id, expected.kind) {
  file <- paste0(
    switch(
      expected.kind,
      geometry = "geometries", sampling = "samplings",
      truth = "truths", response = "responses"),
    ".csv")
  row <- .synthetic.registry.row(file, "component.id", component.id)
  if (!identical(row$kind, expected.kind)) {
    stop("Registry component kind is inconsistent for ",
         component.id, ".", call. = FALSE)
  }
  component <- .new.synthetic.component(
    kind = expected.kind,
    family = row$family,
    version = as.integer(row$version),
    parameters = .synthetic.registry.json.value(row$parameters),
    subclass = .synthetic.registry.component.subclass(
      expected.kind, row$family))
  if (!.is.synthetic.component(component, expected.kind) ||
      !identical(component$family, row$family) ||
      !identical(component$version, as.integer(row$version)) ||
      !identical(.synthetic.sha256(component), row$component.sha256)) {
    stop("Registry component payload failed identity validation for ",
         component.id, ".", call. = FALSE)
  }
  component
}

.synthetic.registry.resolve <- function(recipe.id) {
  row <- .synthetic.registry.row("recipes.csv", "recipe.id", recipe.id)
  if (!identical(row$status, "active") ||
      !identical(as.integer(row$version), 1L)) {
    stop("Registry recipe status or version is unsupported.",
         call. = FALSE)
  }
  spec <- synthetic.spec(
    geometry = .synthetic.registry.component(row$geometry.id, "geometry"),
    sampling = .synthetic.registry.component(row$sampling.id, "sampling"),
    truth = .synthetic.registry.component(row$truth.id, "truth"),
    response = .synthetic.registry.component(row$response.id, "response"),
    recipe.id = if (is.na(row$spec.recipe.id) ||
                       !nzchar(row$spec.recipe.id)) NULL else row$spec.recipe.id,
    registry.tag = if (is.na(row$registry.tag) ||
                         !nzchar(row$registry.tag)) NULL else row$registry.tag,
    compatibility = .synthetic.registry.json.value(
      row$compatibility.recipe),
    metadata = .synthetic.registry.json.value(row$metadata))
  if (!identical(spec$specification.sha256, row$specification.sha256)) {
    stop("Synthetic recipe registry hash mismatch for ", recipe.id, ".",
         call. = FALSE)
  }
  spec
}

.synthetic.environment.schema <- function() {
  path <- .synthetic.registry.asset("environment_fingerprint_schema.csv")
  if (!nzchar(path) || !file.exists(path)) {
    stop("Environment fingerprint schema is unavailable.", call. = FALSE)
  }
  schema <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("field", "type", "required", "description")
  if (!identical(names(schema), required) ||
      anyDuplicated(schema$field) ||
      !all(schema$type %in% c("character", "integer", "double")) ||
      !all(schema$required)) {
    stop("Environment fingerprint schema is malformed.", call. = FALSE)
  }
  schema
}

.synthetic.file.identity <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(list(path = "unavailable", sha256 = "unavailable"))
  }
  path <- normalizePath(path, mustWork = TRUE)
  list(
    path = path,
    sha256 = digest::digest(
      file = path, algo = "sha256", serialize = FALSE))
}

.synthetic.math.runtime.identity <- function() {
  candidates <- c(
    file.path(R.home("lib"), "libR.dylib"),
    file.path(R.home("lib"), "libR.so"),
    file.path(R.home("bin"), "x64", "R.dll"),
    file.path(R.home("bin"), "R.dll"))
  lib.r <- candidates[file.exists(candidates)][1L]
  if (is.na(lib.r)) return("unavailable")
  lib.r <- normalizePath(lib.r, mustWork = TRUE)
  linked <- character()
  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      nzchar(Sys.which("otool"))) {
    output <- tryCatch(
      system2("otool", c("-L", lib.r), stdout = TRUE, stderr = FALSE),
      error = function(error) character())
    linked <- trimws(grep("libSystem|libm", output, value = TRUE))
  } else if (nzchar(Sys.which("ldd"))) {
    output <- tryCatch(
      system2("ldd", lib.r, stdout = TRUE, stderr = FALSE),
      error = function(error) character())
    linked <- trimws(grep("libm", output, value = TRUE))
  }
  paste0(
    "libR.path=", lib.r,
    ";libR.sha256=",
    digest::digest(file = lib.r, algo = "sha256", serialize = FALSE),
    ";linked.math=", if (length(linked)) {
      paste(linked, collapse = " | ")
    } else {
      "none-reported"
    })
}

.synthetic.rng.version <- function() {
  as.character(getRversion())
}

.synthetic.environment.fingerprint <- function(rng.policy = "legacy") {
  rng.policy <- match.arg(rng.policy, c("legacy", "named.stream.v1"))
  soft <- extSoftVersion()
  compiler.command <- tryCatch(
    trimws(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "config", "CXX17"),
      stdout = TRUE, stderr = FALSE))[1L],
    error = function(e) "")
  compiler.version <- if (nzchar(compiler.command)) {
    tryCatch(
      trimws(system2(
        compiler.command, "--version",
        stdout = TRUE, stderr = FALSE))[1L],
      error = function(e) "")
  } else {
    ""
  }
  compiler <- paste(
    c(compiler.command, compiler.version)[
      nzchar(c(compiler.command, compiler.version))],
    collapse = " | ")
  if (!nzchar(compiler)) {
    compiler <- "unavailable"
  }
  blas <- .synthetic.file.identity(unname(soft["BLAS"]))
  lapack <- .synthetic.file.identity(La_library())
  dependencies <- c("dgraphs", "digest", "jsonlite", "MASS", "Matrix")
  dependency.versions <- paste(
    paste0(
      dependencies, "=",
      vapply(dependencies, function(package) {
        tryCatch(
          as.character(utils::packageVersion(package)),
          error = function(e) "unavailable")
      }, character(1))),
    collapse = ";")
  info <- Sys.info()
  fingerprint <- list(
    r.version = R.version.string,
    platform = R.version$platform,
    architecture = R.version$arch,
    compiler = compiler,
    os.name = unname(info["sysname"]),
    os.release = unname(info["release"]),
    os.version = unname(info["version"]),
    os.machine = unname(info["machine"]),
    endianness = .Platform$endian,
    rng.policy = rng.policy,
    rng.kind = if (rng.policy == "legacy") {
      "Mersenne-Twister"
    } else "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection",
    rng.version = .synthetic.rng.version(),
    blas.path = blas$path,
    blas.sha256 = blas$sha256,
    lapack.path = lapack$path,
    lapack.sha256 = lapack$sha256,
    lapack.version = as.character(La_version()),
    math.runtime = .synthetic.math.runtime.identity(),
    geosmooth.version = as.character(utils::packageVersion("geosmooth")),
    registry.schema.version = 1L,
    evaluator.registry.version = 1L,
    dependency.versions = dependency.versions
  )
  schema <- .synthetic.environment.schema()
  if (!identical(names(fingerprint), schema$field)) {
    stop("Environment fingerprint implementation and schema disagree.",
         call. = FALSE)
  }
  fingerprint
}

.synthetic.environment.fingerprint.id <- function(fingerprint) {
  paste0("synthetic-env-v1.", .synthetic.sha256(fingerprint))
}

.synthetic.environment.fingerprint.from.row <- function(row) {
  schema <- .synthetic.environment.schema()
  stats::setNames(lapply(seq_len(nrow(schema)), function(i) {
    value <- row[[schema$field[i]]]
    switch(
      schema$type[i],
      character = as.character(value),
      integer = as.integer(value),
      double = as.double(value))
  }), schema$field)
}

.synthetic.environment.fingerprint.row <- function(fingerprint.id) {
  row <- .synthetic.registry.row(
    "environment_fingerprints.csv",
    "environment.fingerprint.id", fingerprint.id)
  fingerprint <- .synthetic.environment.fingerprint.from.row(row)
  expected.id <- .synthetic.environment.fingerprint.id(fingerprint)
  if (!identical(row$fingerprint.sha256, .synthetic.sha256(fingerprint)) ||
      !identical(fingerprint.id, expected.id)) {
    stop("Environment fingerprint identity is inconsistent.", call. = FALSE)
  }
  row
}

.synthetic.checksum.scope.matches <- function(row, current = NULL) {
  environment.row <- .synthetic.environment.fingerprint.row(
    row$environment.fingerprint.id)
  expected <- .synthetic.environment.fingerprint.from.row(environment.row)
  if (is.null(current)) {
    current <- .synthetic.environment.fingerprint(expected$rng.policy)
  }
  schema <- .synthetic.environment.schema()
  if (!identical(names(current), schema$field)) return(FALSE)
  all(vapply(schema$field, function(field) {
    expected.value <- expected[[field]]
    actual <- current[[field]]
    if (is.numeric(actual)) {
      identical(as.numeric(expected.value), as.numeric(actual))
    } else {
      identical(as.character(expected.value), as.character(actual))
    }
  }, logical(1)))
}

.synthetic.registry.fixture.path <- function(relative.path) {
  relative.path <- .synthetic.scalar.character(
    relative.path, "fixture.path")
  dev <- file.path(relative.path)
  if (file.exists(dev)) return(dev)
  installed.relative <- sub("^inst/", "", relative.path)
  system.file(installed.relative, package = "geosmooth")
}

.synthetic.numeric.fields <- function(x) {
  as.double(strsplit(x, ";", fixed = TRUE)[[1L]])
}

.synthetic.gaussian.mixture.catalog <- function() {
  path <- .synthetic.registry.asset("gaussian_mixture_1d_shapes.csv")
  if (!nzchar(path) || !file.exists(path)) {
    stop("The one-dimensional shape registry is unavailable.", call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

.synthetic.gaussian.variant.catalog <- function() {
  path <- .synthetic.registry.asset("gaussian_mixture_1d_variants.csv")
  if (!nzchar(path) || !file.exists(path)) {
    stop("The one-dimensional variant registry is unavailable.",
         call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

.synthetic.ssrhe.recipe.ids <- function() {
  one.d <- as.vector(outer(
    sprintf("S%02d", 1:16), paste0("V", 1:3), paste, sep = "."))
  flat <- paste0("ssrhe.flat.d", 2:4, ".v1")
  curvature <- c(35L, 85L, 125L, 150L, 200L)
  quad <- unlist(lapply(2:4, function(d) {
    unlist(lapply(seq_len(d), function(k) {
      sprintf(
        "ssrhe.quadform.d%d.k%d.c%03d.v1", d, k, curvature)
    }), use.names = FALSE)
  }), use.names = FALSE)
  c(one.d, flat, quad)
}

.synthetic.one.d.recipe.spec <- function(recipe.id) {
  pieces <- strsplit(recipe.id, ".", fixed = TRUE)[[1L]]
  shapes <- .synthetic.gaussian.mixture.catalog()
  variants <- .synthetic.gaussian.variant.catalog()
  shape <- shapes[shapes$shape.id == pieces[1L], , drop = FALSE]
  variant <- variants[variants$variant.id == pieces[2L], , drop = FALSE]
  if (nrow(shape) != 1L || nrow(variant) != 1L) {
    stop("Unknown one-dimensional SSRHE recipe: ", recipe.id,
         call. = FALSE)
  }
  sampling <- switch(
    variant$variant.id,
    V1 = synthetic.sampling.uniform.interval(0, 1, order = "ascending"),
    V2 = synthetic.sampling.truncated.normal(
      mean = 0.5, sd = 0.2, lower = 0, upper = 1,
      order = "ascending", algorithm = "inverse.cdf.v1"),
    V3 = synthetic.sampling.gapped.uniform(
      intervals = rbind(c(0.02, 0.42), c(0.58, 0.98)),
      allocation = "fixed.proportion",
      probabilities = c(0.42, 0.58),
      rounding = "floor.first.remainder.last",
      order = "ascending",
      algorithm = "sequential.interval.runif.v1"))
  response <- switch(
    variant$variant.id,
    V1 = synthetic.response.gaussian(0.08),
    V2 = synthetic.response.heteroskedastic.gaussian(0.045, 0.14),
    V3 = synthetic.response.laplace.outlier(
      laplace.scale = 0.055, outlier.fraction = 0.04,
      outlier.sd = 0.45, minimum.outliers = 1L,
      count.rounding = "round",
      laplace.algorithm = "uniform.inverse.v1"))
  synthetic.spec(
    synthetic.quadform(1L, 1L, frame = "canonical"),
    sampling,
    synthetic.truth.gaussian.mixture(
      centers = matrix(.synthetic.numeric.fields(shape$means), ncol = 1L),
      scales = .synthetic.numeric.fields(shape$sds),
      weights = .synthetic.numeric.fields(shape$weights),
      normalize = "sample.max"),
    response,
    recipe.id = recipe.id,
    registry.tag = paste0(pieces[1L], "_", pieces[2L]),
    compatibility = list(
      seed.policy = "ssrhe.1d.v1",
      source = "ssrhe_order3_l1_validation_helpers.R"),
    metadata = list(
      shape.id = shape$shape.id,
      shape.name = shape$shape.name,
      variant.id = variant$variant.id,
      variant.name = variant$variant.name))
}

.synthetic.ssrhe.surface.recipe.spec <- function(recipe.id) {
  is.flat <- grepl("^ssrhe\\.flat\\.", recipe.id)
  if (is.flat) {
    match <- regexec("^ssrhe\\.flat\\.d([234])\\.v1$", recipe.id)
    pieces <- regmatches(recipe.id, match)[[1L]]
    if (!length(pieces)) return(NULL)
    d <- as.integer(pieces[2L])
    forms <- list()
    ambient.dim <- d
    metadata <- list(family = "flat", intrinsic.dim = d)
    seed.policy <- "ssrhe.flat.v1"
  } else {
    match <- regexec(
      "^ssrhe\\.quadform\\.d([234])\\.k([1-4])\\.c([0-9]{3})\\.v1$",
      recipe.id)
    pieces <- regmatches(recipe.id, match)[[1L]]
    if (!length(pieces)) return(NULL)
    d <- as.integer(pieces[2L])
    index.k <- as.integer(pieces[3L])
    curvature <- as.integer(pieces[4L]) / 100
    if (index.k > d ||
        !as.integer(pieces[4L]) %in% c(35L, 85L, 125L, 150L, 200L)) {
      return(NULL)
    }
    coefficient.base <- switch(
      as.character(d),
      "2" = c(1, 1),
      "3" = c(1, 1, 1.5),
      "4" = c(1, 1, 2, 4))
    signs <- c(rep(1, index.k), rep(-1, d - index.k))
    coefficients <- curvature * coefficient.base
    forms <- list(diag(signs * coefficients, d))
    ambient.dim <- d + 1L
    metadata <- list(
      family = "quadform", intrinsic.dim = d,
      index.k = index.k, curvature = curvature,
      coefficients = coefficients)
    seed.policy <- "ssrhe.quadform.v1"
  }
  synthetic.spec(
    synthetic.quadform(
      d, ambient.dim, forms = forms, frame = "canonical"),
    synthetic.sampling.uniform.box(-1, 1, order = "draw"),
    synthetic.truth.named("ssrhe.shared.smooth.v1"),
    synthetic.response.gaussian(0.08),
    recipe.id = recipe.id,
    registry.tag = recipe.id,
    compatibility = list(
      seed.policy = seed.policy,
      source = "ssrhe_order3_l1_validation_helpers.R"),
    metadata = metadata)
}

.synthetic.poly2.plan <- c(
  b0 = 0.5, b1 = 1.0, b2 = -0.7,
  b11 = 0.4, b12 = 0.3, b22 = -0.6)

.synthetic.recipe.defaults <- function() {
  list(
    G1 = list(D = 2L, degree = 2L, sigma = 0, poly.coef = NULL),
    G2 = list(d = 2L, D = 3L, sigma = 0,
              frame.seed = NULL, poly.coef = NULL),
    G3a = list(R = 1, truth = "smooth", sigma = 0.1, rho0 = 1),
    G3b = list(R = 2, rho0 = 1, truth = "smooth", sigma = 0.1),
    G3c = list(pitch = 0.2, sigma = 0.1),
    G3d = list(
      R1 = 1, R2 = 0.35,
      u.range = c(-pi / 2, pi / 2),
      v.range = c(-pi / 2, pi / 2),
      truth = "smooth", sigma = 0.1),
    G4 = list(fracA = 0.5, eta = 0.02, sigma = 0.1),
    G5 = list(K = 40L, m = 20L, sigma.x = 0.05,
              rho = 0.6, sigma = 0.1),
    G6 = list(prevalence = 0.5, amplitude = 1.5,
              clip = c(0.05, 0.95)),
    G7 = list(D = 5L, alpha = 1, zero.fraction = 0.5,
              zero.parts = NULL, sigma = 0.1))
}

.synthetic.recipe.spec <- function(recipe.id, args) {
  compatibility <- list(
    recipe = recipe.id,
    seed.policy = "geosmooth.g1.g7.v1")
  if (recipe.id == "G1") {
    D <- .synthetic.scalar.integer(args$D, "D", 1L)
    coefficients <- args$poly.coef
    if (is.null(coefficients)) {
      if (D != 2L) {
        stop("G1 D > 2 requires poly.coef.", call. = FALSE)
      }
      coefficients <- .synthetic.poly2.plan
    }
    return(synthetic.spec(
      synthetic.quadform(D, D, forms = list(), frame = "canonical"),
      synthetic.sampling.uniform.box(-1, 1, order = "draw"),
      synthetic.truth.polynomial(coefficients, "latent"),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g1.ambient.polynomial.v1", registry.tag = "G1",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G2") {
    d <- .synthetic.scalar.integer(args$d, "d", 1L)
    D <- .synthetic.scalar.integer(args$D, "D", d + 1L)
    coefficients <- args$poly.coef
    if (is.null(coefficients)) {
      if (d != 2L) stop("G2 d != 2 requires poly.coef.", call. = FALSE)
      coefficients <- .synthetic.poly2.plan
    }
    compatibility$frame.seed <- args$frame.seed
    return(synthetic.spec(
      synthetic.quadform(
        d, D, forms = list(), frame = "random.orthonormal",
        frame.algorithm = "legacy.base.qr.gaussian.v1"),
      synthetic.sampling.uniform.box(-1, 1, order = "draw"),
      synthetic.truth.polynomial(coefficients, "latent"),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g2.flat.embedded.v1", registry.tag = "G2",
      compatibility = compatibility, metadata = args))
  }
  truth.2d <- function(kind) {
    if (identical(kind, "linear")) {
      synthetic.truth.named("intrinsic.linear.first.v1")
    } else {
      synthetic.truth.named("intrinsic.smooth2.v1")
    }
  }
  if (recipe.id == "G3a") {
    R <- .synthetic.scalar.double(args$R, "R", 0, strict = TRUE)
    return(synthetic.spec(
      synthetic.quadform(
        2L, 3L, forms = list(diag(1 / (2 * R), 2L)),
        frame = "canonical"),
      synthetic.sampling.uniform.disk(args$rho0),
      truth.2d(args$truth),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g3a.paraboloid.v1", registry.tag = "G3a",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G3b") {
    return(synthetic.spec(
      synthetic.sphere.cap(
        radius = args$R, footprint.radius = args$rho0),
      synthetic.sampling.uniform.disk(args$rho0),
      truth.2d(args$truth),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g3b.sphere.cap.v1", registry.tag = "G3b",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G3c") {
    return(synthetic.spec(
      synthetic.helix(pitch = args$pitch),
      synthetic.sampling.uniform.interval(0, 2 * pi),
      synthetic.truth.named("helix.sin.v1"),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g3c.helix.v1", registry.tag = "G3c",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G3d") {
    return(synthetic.spec(
      synthetic.torus.patch(
        major.radius = args$R1, minor.radius = args$R2,
        u.range = args$u.range, v.range = args$v.range),
      synthetic.sampling.uniform.rectangle(
        c(args$u.range[1], args$v.range[1]),
        c(args$u.range[2], args$v.range[2])),
      truth.2d(args$truth),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g3d.torus.patch.v1", registry.tag = "G3d",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G4") {
    return(synthetic.spec(
      .synthetic.g4.segment.rectangle(args$eta),
      synthetic.sampling.stratified(
        args$fracA, "legacy.g4.stratum.sequential.v1"),
      synthetic.truth.named("g4.ridge.v1"),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g4.stratified.v1", registry.tag = "G4",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G5") {
    return(synthetic.spec(
      synthetic.quadform(2L, 2L, forms = list(), frame = "canonical"),
      synthetic.sampling.clustered(
        args$K, args$m, within.sd = args$sigma.x),
      synthetic.truth.named("g5.smooth.v1"),
      synthetic.response.clustered.gaussian(args$sigma, args$rho),
      recipe.id = "g5.clustered.v1", registry.tag = "G5",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G6") {
    return(synthetic.spec(
      synthetic.quadform(2L, 2L, forms = list(), frame = "canonical"),
      synthetic.sampling.uniform.box(-1, 1, order = "draw"),
      synthetic.truth.logit(
        amplitude = args$amplitude,
        target.prevalence = args$prevalence, clip = args$clip),
      synthetic.response.bernoulli(),
      recipe.id = "g6.binary.surface.v1", registry.tag = "G6",
      compatibility = compatibility, metadata = args))
  }
  if (recipe.id == "G7") {
    D <- .synthetic.scalar.integer(args$D, "D", 2L)
    zero.parts <- if (is.null(args$zero.parts)) D else args$zero.parts
    return(synthetic.spec(
      synthetic.simplex(D),
      synthetic.sampling.dirichlet.zeros(
        args$alpha, args$zero.fraction, zero.parts,
        zero.row.policy = "first"),
      synthetic.truth.named("g7.smooth.v1"),
      synthetic.response.gaussian(args$sigma),
      recipe.id = "g7.simplex.zeros.v1", registry.tag = "G7",
      compatibility = compatibility, metadata = args))
  }
  stop("Unknown synthetic recipe ID: ", recipe.id, call. = FALSE)
}

#' List maintained synthetic recipe IDs
#' @return Character recipe IDs.
#' @examples
#' head(synthetic.registry.ids())
#' @export
synthetic.registry.ids <- function() {
  c(names(.synthetic.recipe.defaults()), .synthetic.ssrhe.recipe.ids())
}

# Build registry rows without consulting generated ledgers. This function is
# private to registry regeneration and parameterized G-family compatibility.
.synthetic.registry.spec.from.code <- function(
    recipe.id, parameters = list()) {
  recipe.id <- .synthetic.scalar.character(recipe.id, "recipe.id")
  if (grepl("^S[0-9]{2}\\.V[1-3]$", recipe.id)) {
    if (length(parameters)) {
      stop("One-dimensional registry recipes do not accept overrides.",
           call. = FALSE)
    }
    return(.synthetic.one.d.recipe.spec(recipe.id))
  }
  if (grepl("^ssrhe\\.(flat|quadform)\\.", recipe.id)) {
    if (length(parameters)) {
      stop("SSRHE surface registry recipes do not accept overrides.",
           call. = FALSE)
    }
    spec <- .synthetic.ssrhe.surface.recipe.spec(recipe.id)
    if (is.null(spec)) {
      stop("Unknown synthetic recipe ID: ", recipe.id, call. = FALSE)
    }
    return(spec)
  }
  defaults <- .synthetic.recipe.defaults()[[recipe.id]]
  if (is.null(defaults)) {
    stop("Unknown synthetic recipe ID: ", recipe.id, call. = FALSE)
  }
  if (!is.list(parameters) || is.null(names(parameters)) && length(parameters)) {
    stop("parameters must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(parameters), names(defaults))
  if (length(unknown)) {
    stop("Unknown parameters for ", recipe.id, ": ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  args <- utils::modifyList(defaults, parameters, keep.null = TRUE)
  spec <- .synthetic.recipe.spec(recipe.id, args)
  spec
}

#' Resolve a maintained synthetic recipe
#'
#' Default recipes are reconstructed from normalized component foreign keys.
#' Parameter overrides are accepted only for the legacy G-family compatibility
#' surface and remain content-bound, non-frozen specifications.
#'
#' @param recipe.id One of `synthetic.registry.ids()`.
#' @param parameters Optional named parameter overrides.
#' @return A validated `synthetic_spec`.
#' @examples
#' synthetic.registry.spec("G1")
#' @export
synthetic.registry.spec <- function(recipe.id, parameters = list()) {
  recipe.id <- .synthetic.scalar.character(recipe.id, "recipe.id")
  if (!is.list(parameters) || is.null(names(parameters)) && length(parameters)) {
    stop("parameters must be a named list.", call. = FALSE)
  }
  if (!length(parameters)) {
    return(.synthetic.registry.resolve(recipe.id))
  }
  if (!recipe.id %in% names(.synthetic.recipe.defaults())) {
    stop("Registry parameter overrides are supported only for G1--G7.",
         call. = FALSE)
  }
  .synthetic.registry.spec.from.code(recipe.id, parameters)
}

#' Resolve the legacy seed for an SSRHE registry recipe
#'
#' @param recipe.id A one-dimensional, flat, or quadform SSRHE recipe ID.
#' @param replicate Positive replicate number.
#' @param n Sample size. Required for flat and quadform recipes.
#' @param base.seed One-dimensional suite base seed.
#' @return A validated integer seed.
#' @examples
#' synthetic.registry.seed("S16.V1", replicate = 2L)
#' @export
synthetic.registry.seed <- function(
    recipe.id, replicate = 1L, n = NULL, base.seed = 273001L) {
  recipe.id <- .synthetic.scalar.character(recipe.id, "recipe.id")
  replicate <- .synthetic.scalar.integer(replicate, "replicate", 1L)
  if (grepl("^S[0-9]{2}\\.V[1-3]$", recipe.id)) {
    pieces <- strsplit(recipe.id, ".", fixed = TRUE)[[1L]]
    shape <- as.integer(sub("^S", "", pieces[1L]))
    variant <- as.integer(sub("^V", "", pieces[2L]))
    return(.validate.synthetic.seed(
      as.double(base.seed) + shape * 1000 + variant * 100 + replicate))
  }
  n <- .synthetic.scalar.integer(n, "n", 1L)
  spec <- .synthetic.ssrhe.surface.recipe.spec(recipe.id)
  if (is.null(spec)) {
    stop("recipe.id is not an SSRHE registry recipe.", call. = FALSE)
  }
  metadata <- spec$metadata
  if (metadata$family == "flat") {
    seed <- 283000 + 1000 * metadata$intrinsic.dim + 10 * n + replicate
  } else {
    seed <- 291000 + 1000 * metadata$intrinsic.dim +
      100 * metadata$index.k + round(100 * metadata$curvature) +
      10 * n + replicate
  }
  .validate.synthetic.seed(seed)
}

.synthetic.instance.path <- function() {
  dev <- file.path("inst", "synthetic_registry", "instances.csv")
  if (file.exists(dev)) return(dev)
  system.file("synthetic_registry", "instances.csv", package = "geosmooth")
}

#' Materialize a frozen synthetic instance
#' @param instance.id Frozen instance ID.
#' @param validate Whether to validate and checksum the result.
#' @return A `synthetic_dataset` whose `dataset.id` is `instance.id`.
#' @examples
#' materialize.synthetic.instance("geosmooth.g1.default.v1")
#' @export
materialize.synthetic.instance <- function(instance.id, validate = TRUE) {
  instance.id <- .synthetic.scalar.character(instance.id, "instance.id")
  path <- .synthetic.instance.path()
  if (!nzchar(path) || !file.exists(path)) {
    stop("Synthetic instance registry is not installed.", call. = FALSE)
  }
  row <- .synthetic.registry.row("instances.csv", "instance.id", instance.id)
  checksum.row <- .synthetic.registry.row(
    "checksums.csv", "checksum.id", row$checksum.id)
  if (!identical(checksum.row$instance.id, instance.id) ||
      !identical(checksum.row$content.sha256, row$content.sha256)) {
    stop("Frozen instance checksum foreign key is inconsistent.",
         call. = FALSE)
  }
  fixture.path <- .synthetic.registry.fixture.path(row$fixture.path)
  if (!nzchar(fixture.path) || !file.exists(fixture.path)) {
    stop("Frozen instance fixture is unavailable: ", row$fixture.path,
         call. = FALSE)
  }
  spec <- synthetic.registry.spec(row$recipe.id)
  object <- materialize.synthetic(
    spec, n = row$n, seed = row$seed, rng.policy = row$rng.policy,
    validate = FALSE)
  object$dataset.id <- instance.id
  attr(object, "frozen.instance.id") <- instance.id
  if (validate) {
    validate.synthetic.dataset(object)
    object <- .synthetic.verify.frozen.instance(
      object = object,
      fixture = readRDS(fixture.path),
      instance.row = row,
      checksum.row = checksum.row)
  }
  object
}

.synthetic.verify.frozen.instance <- function(
    object, fixture, instance.row, checksum.row,
    current.fingerprint = NULL) {
  checksum <- synthetic.dataset.checksum(object)
  if (.synthetic.checksum.scope.matches(
      checksum.row, current = current.fingerprint)) {
    if (!identical(checksum, instance.row$content.sha256)) {
      stop(
        "Frozen instance checksum mismatch for ", instance.row$instance.id,
        ": expected ", instance.row$content.sha256,
        ", obtained ", checksum, call. = FALSE)
    }
    attr(object, "verification.scope") <- "exact-environment"
    return(object)
  }
  comparison <- compare.synthetic.dataset(
    object, fixture,
    tolerance = c(
      checksum.row$scientific.abs.tolerance,
      checksum.row$scientific.rel.tolerance))
  if (!comparison$equal) {
    stop(
      "Frozen instance failed cross-environment scientific parity: ",
      paste(comparison$mismatches, collapse = ", "), call. = FALSE)
  }
  attr(object, "verification.scope") <- "scientific-parity"
  object
}
