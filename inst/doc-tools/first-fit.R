## ---- first-fit-model
library(geosmooth)
set.seed(1)
X <- matrix(seq(0, 1, length.out = 60), ncol = 1)
truth <- sin(2 * pi * X[, 1])
y <- truth + rnorm(nrow(X), sd = 0.08)
fit <- fit.lps(
  X, y, foldid = rep(1:5, length.out = nrow(X)),
  support.grid = c(8L, 12L, 16L), degree.grid = 0:1,
  kernel.grid = "gaussian"
)
fit$selected

## ---- first-fit-plot
plot(X[, 1], y, pch = 16, col = "#666666", xlab = "Coordinate",
     ylab = "Response", main = "A local polynomial fit to noisy observations")
lines(X[, 1], truth, col = "#D55E00", lty = 2, lwd = 2)
lines(X[, 1], fitted(fit), col = "#0072B2", lty = 1, lwd = 2)
legend("topright", c("Observations", "Fitted response", "Known synthetic truth"),
       pch = c(16, NA, NA), lty = c(NA, 1, 2),
       col = c("#666666", "#0072B2", "#D55E00"), bty = "n", cex = .8)
