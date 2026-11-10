## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>")
library(geosmooth)

## ----seam---------------------------------------------------------------------
X <- rbind(c(1, .8, .2), c(.8, 1, .2), c(1, .9, .2))
query <- c(1, .95, .2)
chart <- unfold.simplex.chart(X, query)
chart[c("certified", "status", "distance.error")]
chart$coordinates

## ----pca-geometry-------------------------------------------------------------
t.seam <- c(seq(-.2, -.01, length.out=20), seq(.01, .2, length.out=20))
seam <- cbind(ifelse(t.seam < 0, 1, 1-t.seam),
              ifelse(t.seam < 0, 1+t.seam, 1), .2)
flat <- unfold.simplex.chart(seam, c(1, .95, .2))
variance.one <- function(x) {
  s <- svd(scale(x, center=TRUE, scale=FALSE), nu=0, nv=0)$d^2
  s[1]/sum(s)
}
data.frame(chart=c("Ambient PCA", "Unfolded PCA"),
           first.PC.variance=c(variance.one(seam), variance.one(flat$coordinates)))

## ----fit----------------------------------------------------------------------
set.seed(2710)
n <- 100L
t <- runif(n, -.2, .2)
U <- cbind(ifelse(t < 0, 1, 1-t), ifelse(t < 0, 1+t, 1), runif(n, .1, .3))
y <- rbinom(n, 1, plogis(-1 + 8*t))
fit <- fit.local.likelihood(U, y, X.eval=U[81:100, ],
  likelihood.family="bernoulli", chart.geometry="linf.simplex",
  coordinate.method="local.pca", unfold.faces=TRUE,
  support.size=25, chart.dim=2, lambda.ridge=1)
fit$diagnostics$status.counts
head(fit$diagnostics$per.eval)

## ----heldout------------------------------------------------------------------
pls <- fit.local.likelihood(U[1:80, ], y[1:80], X.eval=U[81:100, ],
  likelihood.family="bernoulli", chart.geometry="linf.simplex",
  coordinate.method="local.pls", unfold.faces=TRUE,
  support.size=25, chart.dim=2, lambda.ridge=1)
range(pls$fitted.values)

## ----cv-----------------------------------------------------------------------
fold <- rep(1:4, length.out=80)
selected <- fit.local.likelihood(U[1:80, ], y[1:80], X.eval=U[81:100, ],
  likelihood.family="bernoulli", chart.geometry="linf.simplex",
  coordinate.method="local.pls", unfold.faces=TRUE, support.size=20,
  chart.dim.grid=c(1,2), lambda.ridge.grid=c(.1,1), foldid=fold)
selected$cv.table

