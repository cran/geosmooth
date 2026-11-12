## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", fig.width=7, fig.height=3.5)
library(geosmooth)

## ----fixed--------------------------------------------------------------------
set.seed(41)
X <- matrix(runif(72,-1,1),ncol=2)
colnames(X) <- c("x1","x2")
mean.response <- sin(X[,1]+0.6*X[,2]) + 0.3*X[,2]^2
y <- mean.response + 0.05*cos(seq_len(nrow(X)))
directions <- rbind(c(1,0.6), c(0,1))
fixed <- fit.projection.trend.filtering(
    X, y, lambda=0.001, directions=directions, direction.search="fixed",
    preprocessing="none", curve.ridge=0.05)
fixed
fixed$checks

## ----learned------------------------------------------------------------------
learned <- fit.projection.trend.filtering(
    X, y, lambda=0.001, nterms=1, curve.ridge=0.05, seed=17,
    control=list(starts=1, sparse.directions=2, dense.directions=2,
                 refinement.passes=1, refinement.proposals=8))
learned
learned$directions
learned$search.starts

## ----binary-------------------------------------------------------------------
set.seed(72)
probability <- plogis(-0.4 + 1.5*sin(X[,1]+0.6*X[,2]))
binary <- rbinom(nrow(X), 1, probability)
binary.fit <- fit.projection.trend.filtering(
    X, binary, lambda=0.001, family="binomial", curve.ridge=0.05,
    directions=matrix(c(1,0.6),1), direction.search="fixed",
    preprocessing="none")
head(predict(binary.fit))

## ----cv, eval=FALSE-----------------------------------------------------------
# cv <- cv.projection.trend.filtering(
#     X, y, foldid=folds, groups=participant,
#     nterms.grid=c(0,1,2), lambda.grid=c(Inf,0.01,0.001),
#     curve.ridge.grid=c(0.001,0.01), seed=17)
# cv$cv.table
# predict(cv, newdata=new.X)

## ----prediction---------------------------------------------------------------
new.X <- rbind(c(-1.5,0), c(0,0), c(1.5,0))
colnames(new.X) <- colnames(X)
prediction <- predict(binary.fit, new.X)
prediction
attr(prediction, "outside.range")

## ----diagnostic-plot, fig.alt="Estimated binary probabilities plotted against the simulated probabilities for this illustrative fixed-direction fit."----
plot(probability, as.numeric(predict(binary.fit)), pch=19,
     xlab="Simulated probability", ylab="Estimated probability",
     xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2,col="gray50")

