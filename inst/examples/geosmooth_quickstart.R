# Runnable quick-start examples for geosmooth.
#
# Run from an installed package, or from the source tree after:
# pkgload::load_all("/Users/pgajer/current_projects/geosmooth")

set.seed(1)
x <- seq(0, 1, length.out = 60)
X <- cbind(x = x)
y <- sin(2 * pi * x) + stats::rnorm(length(x), sd = 0.08)
foldid <- rep(1:5, length.out = length(y))

lps.fit <- geosmooth::fit.lps(
    X = X,
    y = y,
    foldid = foldid,
    support.grid = c(8L, 12L, 16L),
    degree.grid = 0:1,
    kernel.grid = c("gaussian", "tricube")
)
lps.pred <- stats::predict(lps.fit, X)

malps.fit <- geosmooth::fit.malps(
    X = X,
    y = y,
    degree = 1L,
    support.type = "knn",
    support.size = 12L,
    kernel = "tricube",
    support.selection = "fixed",
    coordinate.method = "coordinates"
)

lpl.op <- geosmooth::lpl.tf.operator(
    X = X,
    degree = 1L,
    support.type = "knn",
    support.size = 12L,
    kernel = "gaussian",
    coordinate.method = "coordinates"
)

slpl.op <- geosmooth::slpl.tf.operator(
    X = X,
    degree = 1L,
    support.type = "knn",
    support.size = 12L,
    kernel = "gaussian",
    coordinate.method = "coordinates"
)

if (requireNamespace("genlasso", quietly = TRUE)) {
    lpl.fit <- geosmooth::fit.lpl.tf(
        y = y,
        operator = lpl.op,
        lambda = 0.1,
        lambda.selection = "fixed"
    )

    slpl.fit <- geosmooth::fit.slpl.tf(
        y = y,
        operator = slpl.op,
        lambda1 = 0.1,
        lambda2 = 0.01,
        lambda.selection = "fixed"
    )
}

grid <- expand.grid(x = seq(0, 1, length.out = 5),
                    y = seq(0, 1, length.out = 5))
X2 <- as.matrix(grid)
y2 <- sin(2 * pi * X2[, 1]) + 0.25 * X2[, 2]

ssrhe.fit <- geosmooth::fit.ssrhe.hessian.regression(
    X = X2,
    y = y2,
    k = 12L,
    tangent.dim = 2L,
    lambda1 = 0.05,
    return.local.diagnostics = FALSE
)

list(
    lps.selected = lps.fit$selected,
    lps.rmse = sqrt(mean((lps.pred - y)^2)),
    malps.rmse = sqrt(mean((malps.fit$fitted.values - y)^2)),
    lpl.operator.rows = nrow(lpl.op$A),
    slpl.sync.rows = nrow(slpl.op$C_sync),
    ssrhe.rmse = sqrt(mean((ssrhe.fit$fitted.values - y2)^2))
)
