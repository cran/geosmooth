# Opt-in Bernoulli chart path; legacy density/local-likelihood code is unchanged.
.fit.local.binary.chart <- function(X, y, X.eval, support.size, degree, kernel,
    bandwidth.multiplier, coordinate.method, chart.dim, quadrature.weights,
    lambda.ridge, fallback, max.iter, tol, return.details, chart.geometry,
    unfold.faces, intercept.mass, support.grid, degree.grid, kernel.grid,
    bandwidth.multiplier.grid, lambda.ridge.grid, chart.dim.grid, foldid,
    cv.folds, cv.seed) {

    prepared <- .local.chart.prepare.X.eval(X, X.eval)
    X <- prepared$X
    X.eval <- prepared$X.eval
    n <- nrow(X)
    y <- .local.chart.validate.response(y, n)
    if (any(y < 0 | y > 1)) stop("y must be in [0, 1].", call. = FALSE)
    if (!coordinate.method %in% c("local.pca", "local.pls"))
        stop("Nonlegacy charts require local.pca or local.pls.", call. = FALSE)
    if (!fallback %in% c("degree0", "na"))
        stop("Nonlegacy charts support fallback = 'degree0' or 'na'.", call. = FALSE)
    if (!is.logical(unfold.faces) || length(unfold.faces) != 1L || is.na(unfold.faces))
        stop("unfold.faces must be TRUE or FALSE.", call. = FALSE)
    if (unfold.faces && chart.geometry != "linf.simplex")
        stop("Unfolding requires linf.simplex geometry.", call. = FALSE)
    if (chart.geometry == "linf.simplex") {
        X <- .simplex.chart.normalize(X)
        X.eval <- .simplex.chart.normalize(X.eval)
    }
    quadrature.weights <- .local.chart.validate.quadrature(quadrature.weights, n)
    support.size <- .local.chart.validate.support.size(support.size, n)
    chart.dim <- if (is.null(chart.dim)) max(1L, min(ncol(X), support.size - 1L)) else
        .local.chart.validate.positive.integer(chart.dim, "chart.dim")
    max.iter <- .local.chart.validate.positive.integer(max.iter, "max.iter")
    tol <- .local.chart.validate.positive.scalar(tol, "tol")
    intercept.mass <- .local.chart.validate.positive.scalar(intercept.mass, "intercept.mass")
    validate.controls <- function(support.size, degree, kernel, bandwidth.multiplier,
                                  lambda.ridge, chart.dim) {
        if (!is.numeric(support.size) || length(support.size) != 1L ||
            !is.finite(support.size) || support.size != as.integer(support.size) || support.size < 1 || support.size > n)
            stop("Support candidates must be integer sizes between 1 and nrow(X).", call. = FALSE)
        if (!is.numeric(degree) || length(degree) != 1L || !degree %in% c(0, 1))
            stop("Nonlegacy charts support degree 0 or 1 only.", call. = FALSE)
        if (!kernel %in% c("gaussian", "tricube", "epanechnikov", "triangular"))
            stop("Unknown kernel.", call. = FALSE)
        .local.chart.validate.positive.scalar(bandwidth.multiplier, "bandwidth.multiplier")
        .local.chart.validate.nonnegative.scalar(lambda.ridge, "lambda.ridge")
        .local.chart.validate.positive.integer(chart.dim, "chart.dim")
        invisible(NULL)
    }
    controls <- list(support.size = support.size, degree = degree, kernel = kernel,
        bandwidth.multiplier = bandwidth.multiplier, lambda.ridge = lambda.ridge,
        chart.dim = chart.dim)
    do.call(validate.controls, controls)
    fixed <- list(coordinate.method = coordinate.method, chart.geometry = chart.geometry,
                  unfold.faces = unfold.faces, fallback = fallback,
                  max.iter = max.iter, tol = tol, intercept.mass = intercept.mass)
    run <- function(train, query, controls) do.call(.local.binary.chart.evaluate,
        c(list(X = X[train, , drop = FALSE], y = y[train], X.eval = query,
               quadrature.weights = quadrature.weights[train]), controls, fixed))
    cv.table <- cv.predictions <- NULL
    requested <- !is.null(foldid) || any(!vapply(list(support.grid, degree.grid,
        kernel.grid, bandwidth.multiplier.grid, lambda.ridge.grid, chart.dim.grid), is.null, logical(1)))
    if (requested) {
        if (is.null(foldid)) {
            cv.folds <- .local.chart.validate.positive.integer(cv.folds, "cv.folds")
            if (cv.folds < 2L || cv.folds > n) stop("cv.folds must be between 2 and nrow(X).")
            # Preserve the caller's RNG state, including its absence.
            present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            if (present) old.seed <- get(".Random.seed", envir = .GlobalEnv)
            on.exit(if (present) assign(".Random.seed", old.seed, envir = .GlobalEnv) else
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                    rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
            set.seed(cv.seed)
            foldid <- sample(rep(seq_len(cv.folds), length.out = n))
        }
        if (!is.numeric(foldid) || length(foldid) != n || any(!is.finite(foldid)) ||
            any(foldid < 1 | foldid != trunc(foldid)) || length(unique(foldid)) < 2L)
            stop("foldid must assign each training row to one of at least two positive integer folds.")
        grids <- list(support.size = support.grid, degree = degree.grid,
            kernel = kernel.grid, bandwidth.multiplier = bandwidth.multiplier.grid,
            lambda.ridge = lambda.ridge.grid, chart.dim = chart.dim.grid)
        for (name in names(grids)) if (is.null(grids[[name]])) grids[[name]] <- controls[[name]]
        cv.table <- do.call(expand.grid, c(grids, list(stringsAsFactors = FALSE)))
        if (!nrow(cv.table)) stop("CV grids must not be empty.")
        predictions <- matrix(NA_real_, n, nrow(cv.table))
        cv.table$brier <- cv.table$logloss <- Inf
        cv.table$fallback.count <- NA_integer_
        cv.table$status <- "support.exceeds.training"
        minimum <- min(vapply(unique(foldid), function(f) sum(foldid != f), integer(1)))
        for (j in seq_len(nrow(cv.table))) {
            candidate <- as.list(cv.table[j, names(controls), drop = FALSE])
            do.call(validate.controls, candidate)
            if (candidate$support.size > minimum) next
            failures <- 0L
            for (f in unique(foldid)) {
                held <- which(foldid == f)
                fit <- run(which(foldid != f), X[held, , drop = FALSE], candidate)
                predictions[held, j] <- fit$values
                failures <- failures + sum(fit$details$fallback)
            }
            cv.table$fallback.count[j] <- failures
            if (all(is.finite(predictions[, j]))) {
                pr <- pmin(1 - 1e-15, pmax(1e-15, predictions[, j]))
                cv.table$brier[j] <- mean((y - pr)^2)
                cv.table$logloss[j] <- -mean(y * log(pr) + (1 - y) * log1p(-pr))
                cv.table$status[j] <- "evaluated"
            } else cv.table$status[j] <- "nonfinite.predictions"
        }
        if (!any(is.finite(cv.table$brier))) stop("No CV candidate produced complete finite predictions.")
        best <- which.min(cv.table$brier)
        controls <- as.list(cv.table[best, names(controls), drop = FALSE])
        cv.predictions <- predictions[, best]
    }
    fit <- run(seq_len(n), X.eval, controls)
    details <- fit$details
    structure(list(method.id = "local_likelihood", likelihood.family = "bernoulli",
        X = X, y = y, X.eval = X.eval, fitted.values = fit$values,
        selected = c(controls, fixed, list(likelihood.family = "bernoulli",
            cv.brier.observed = if (requested) min(cv.table$brier) else NA_real_)),
        cv.table = cv.table, foldid = foldid, cv.predictions = cv.predictions,
        quadrature.weights = quadrature.weights,
        diagnostics = list(status.counts = table(details$status),
            fallback.count = sum(details$fallback), fallback.fraction = mean(details$fallback),
            per.eval = if (isTRUE(return.details)) details else NULL)), class = "local_likelihood")
}

.local.binary.chart.evaluate <- function(X, y, X.eval, quadrature.weights,
    support.size, degree, kernel, bandwidth.multiplier, lambda.ridge, chart.dim,
    coordinate.method, chart.geometry, unfold.faces, fallback, max.iter, tol,
    intercept.mass) {
    values <- numeric(nrow(X.eval))
    details <- vector("list", nrow(X.eval))
    for (i in seq_len(nrow(X.eval))) {
        query <- X.eval[i, ]
        distances <- if (chart.geometry == "linf.simplex")
            drop(dgraphs::linf.simplex.distances(rbind(query, X), sources = 1L,
                                                targets = 1L + seq_len(nrow(X)))) else
            sqrt(rowSums(sweep(X, 2L, query)^2))
        idx <- order(distances, seq_along(distances))[seq_len(support.size)]
        local <- X[idx, , drop = FALSE]
        w <- .local.chart.kernel(distances[idx], kernel, bandwidth.multiplier)$weights * quadrature.weights[idx]
        if (sum(w) <= 0) stop("The selected kernel has zero total local weight.")
        w <- w * support.size / sum(w)
        yy <- y[idx]
        constant <- (sum(w * yy) + intercept.mass) / (support.size + 2 * intercept.mass)
        cert <- list(certified = TRUE, status = "not.requested", distance.error = NA_real_,
                     radius = max(distances[idx]), clearance = NA_real_)
        if (unfold.faces) {
            cert <- unfold.simplex.chart(local, query)
            if (cert$certified) { local <- cert$coordinates; query <- cert$query }
        }
        status <- "degree0"
        failed <- FALSE
        gradient <- NA_real_
        dimension <- 0L
        explained <- NA_real_
        values[i] <- constant
        if (degree > 0L) {
            if (!cert$certified) {
                failed <- TRUE
                status <- paste0("chart.", cert$status)
            } else {
                rotation <- tryCatch(.local.binary.rotation(local, yy, w,
                    coordinate.method, chart.dim), error = function(e) NULL)
                if (is.null(rotation)) {
                    failed <- TRUE
                    status <- "rotation.failed"
                } else {
                    centered <- sweep(local, 2L, rotation$center)
                    scores <- centered %*% rotation$rotation
                    score.query <- drop((query - rotation$center) %*% rotation$rotation)
                    dimension <- ncol(scores)
                    if (coordinate.method == "local.pca") {
                        total <- sum(w * rowSums(centered^2))
                        explained <- if (total > 0) sum(w * rowSums(scores^2)) / total else NA_real_
                    }
                    scale <- pmax(1e-12, sqrt(colSums(scores^2 * w) / sum(w)))
                    A <- cbind(1, sweep(scores, 2L, scale, "/"))
                    aq <- c(1, score.query / scale)
                    fit <- .local.binary.solve(A, yy, w, aq, lambda.ridge,
                                               intercept.mass, max.iter, tol)
                    gradient <- fit$gradient
                    status <- fit$status
                    failed <- status != "converged"
                    if (!failed) values[i] <- fit$value
                }
            }
        }
        if (failed && fallback == "na") values[i] <- NA_real_
        details[[i]] <- data.frame(status = status, fallback = failed,
            chart.status = cert$status, chart.certified = if (unfold.faces) cert$certified else NA,
            radius = cert$radius, clearance = cert$clearance, distance.error = cert$distance.error,
            dimension = dimension, pca.variance.fraction = explained,
            effective.support = sum(w)^2 / sum(w^2), gradient = gradient)
    }
    list(values = values, details = do.call(rbind, details))
}
