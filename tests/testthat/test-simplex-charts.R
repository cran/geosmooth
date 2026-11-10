test_that("unfolding flattens a seam and declines unsupported corners", {
    X <- rbind(c(1, .8, .2), c(.8, 1, .2), c(1, .9, .2))
    query <- c(1, .95, .2)
    a <- unfold.simplex.chart(X, query)
    expect_true(a$certified)
    expect_equal(a$coordinates[2, ], c(1, 1.2, .2))
    expect_lt(a$distance.error, 1e-12)
    expect_equal(unfold.simplex.chart(4 * X, 8 * query), a)
    expect_identical(unfold.simplex.chart(rbind(X, c(1, 1, 0)), query)$status, "tied.maxima")
    expect_identical(unfold.simplex.chart(rbind(X, c(.2, .3, 1)), query)$status, "multiple.faces")
    expect_false(unfold.simplex.chart(rbind(c(1, 0, .99), c(0, 1, .99)), query)$certified)
    expect_error(unfold.simplex.chart(-X, query), "nonnegative")
    expect_error(unfold.simplex.chart(X, c(0, 0, 0)), "positive maximum")
    expect_equal(unfold.simplex.chart(matrix(1, 3, 1), 1)$distance.error, 0)
})

chart.fit.test <- function(X, y, ..., X.eval = X, method = "local.pca", unfold = TRUE) {
    fit.local.likelihood(X, y, X.eval, likelihood.family = "bernoulli",
        chart.geometry = "linf.simplex", unfold.faces = unfold,
        coordinate.method = method, ...)
}

test_that("R package reproduces the independent synthetic Python predictions", {
    path <- test_path("fixtures", "simplex-charts")
    reference <- read.csv(file.path(path, "predictions.csv"))
    for (scenario in 0:4) {
        dat <- read.csv(file.path(path, paste0("scenario-", scenario, ".csv")))
        X <- as.matrix(dat[1:320, 1:8])
        query <- as.matrix(dat[321:326, 1:8])
        for (method in unique(reference$method)) for (d in 1:3) for (ridge in c(.1, 1, 10)) {
            fit <- chart.fit.test(X, dat$y[1:320], X.eval = query,
                method = if (method == "unfolded_pls") "local.pls" else "local.pca",
                unfold = method != "ambient_pca", support.size = 60, chart.dim = d,
                lambda.ridge = ridge, max.iter = 70, tol = 1e-7)
            expected <- reference[reference$scenario == scenario & reference$method == method &
                reference$dimension == d & reference$ridge == ridge, ]
            expect_equal(fit$fitted.values, expected$prediction, tolerance = 2e-7)
            expect_equal(fit$diagnostics$fallback.count, 0)
            expect_true(all(fit$diagnostics$per.eval$gradient < 1e-7))
        }
    }
})

test_that("unfolded charts preserve explicit fallback and solver failures", {
    X <- rbind(c(1, .8, .2), c(.8, 1, .2), c(1, 1, .3))
    y <- c(0, 1, 0)
    fit <- chart.fit.test(X, y, support.size = 3, chart.dim = 1)
    expect_equal(fit$diagnostics$fallback.count, 3)
    expect_true(all(is.finite(fit$fitted.values)))
    missing <- chart.fit.test(X, y, support.size = 3, chart.dim = 1, fallback = "na")
    expect_true(all(is.na(missing$fitted.values)))
    for (outcome in c(0, 1)) {
        constant <- chart.fit.test(X, rep(outcome, 3), unfold = FALSE, method = "local.pls",
                                   support.size = 3, chart.dim = 2)
        expect_equal(constant$fitted.values, rep((3 * outcome + .5) / 4, 3))
        expect_equal(constant$diagnostics$fallback.count, 0)
    }
    expect_error(chart.fit.test(X, y, degree = 2), "degree 0 or 1")
    expect_error(chart.fit.test(X, y, chart.dim = "auto"), "numeric chart dimension")
    expect_error(chart.fit.test(X, y, fallback = "zero"), "fallback")
    expect_error(fit.local.likelihood(X, y, chart.geometry = "linf.simplex"), "Bernoulli")
    expect_error(fit.local.likelihood(X, y, coordinate.method = "local.pls"), "nonlegacy")
    dat <- read.csv(test_path("fixtures", "simplex-charts", "scenario-1.csv"))
    failed <- chart.fit.test(as.matrix(dat[1:50, 1:8]), dat$y[1:50],
        unfold = FALSE, support.size = 30, chart.dim = 2, max.iter = 1, tol = 1e-14)
    expect_gt(failed$diagnostics$fallback.count, 0)
    expect_true(all(failed$diagnostics$per.eval$gradient[failed$diagnostics$per.eval$fallback] >= 1e-14))
})

test_that("CV reconstructs PCA/PLS from each training fold only", {
    dat <- read.csv(test_path("fixtures", "simplex-charts", "scenario-1.csv"))[1:36, ]
    X <- as.matrix(dat[, 1:8]); y <- dat$y; folds <- rep(1:3, each = 12)
    for (method in c("local.pca", "local.pls")) {
        fit <- chart.fit.test(X, y, method = method, unfold = FALSE,
            support.size = 15, chart.dim = 2, lambda.ridge = 1, foldid = folds)
        for (f in 1:3) {
            train <- which(folds != f); held <- which(folds == f)
            manual <- chart.fit.test(X[train, ], y[train], X.eval = X[held, ], method = method,
                unfold = FALSE, support.size = 15, chart.dim = 2, lambda.ridge = 1)
            expect_equal(fit$cv.predictions[held], manual$fitted.values)
        }
        flipped <- y; flipped[folds == 1] <- 1 - flipped[folds == 1]
        check <- chart.fit.test(X, flipped, method = method, unfold = FALSE,
            support.size = 15, chart.dim = 2, lambda.ridge = 1, foldid = folds)
        expect_equal(check$cv.predictions[folds == 1], fit$cv.predictions[folds == 1])
    }
    set.seed(123); seed <- .Random.seed
    chart.fit.test(X, y, unfold = FALSE, support.size = 15, chart.dim.grid = c(1, 2), cv.folds = 3)
    expect_identical(.Random.seed, seed)
})

test_that("legacy behavior is unchanged and new summaries print", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1); y <- rep(c(0, 1), 10)
    old <- fit.local.likelihood(X, y, likelihood.family = "bernoulli", coordinate.method = "local.pca")
    explicit <- fit.local.likelihood(X, y, likelihood.family = "bernoulli",
                                    coordinate.method = "local.pca", chart.geometry = "legacy")
    old$call <- explicit$call <- NULL
    expect_equal(old, explicit)
    new <- fit.local.likelihood(X, y, likelihood.family = "bernoulli",
                               coordinate.method = "local.pls", chart.geometry = "euclidean")
    expect_output(print(new), "Local likelihood")
})
