## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", fig.width=7.2,
                    fig.height=4.5, dpi=120)
library(geosmooth)
recompute <- isTRUE(getOption("geosmooth.recompute.vignettes", FALSE))

## ----seam-geometry, fig.cap="A straight latent coordinate appears bent in ambient coordinates. Unfolding rotates the adjacent face into the query face; both axes use the same physical scale.", fig.alt="Two panels show a right-angle seam before unfolding and a straight line afterward."----
t.line <- c(seq(-.2, -.005, length.out=30), seq(.005, .2, length.out=30))
seam <- cbind(ifelse(t.line < 0, 1, 1-t.line),
              ifelse(t.line < 0, 1+t.line, 1), .2)
flat <- unfold.simplex.chart(seam, c(1, .95, .2))
stopifnot(flat$certified)
old.par <- par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(seam[,1:2], asp=1, pch=16, col="#0072B2", xlim=c(.76,1.24),
     ylim=c(.76,1.24), xlab="First coordinate", ylab="Second coordinate",
     main="Ambient coordinates")
plot(flat$coordinates[,1:2], asp=1, pch=16, col="#009E73",
     xlim=c(.76,1.24), ylim=c(.76,1.24), xlab="First coordinate",
     ylab="Unfolded second coordinate", main="Unfolded coordinates")
par(old.par)
variance.one <- function(x) {
    values <- svd(scale(x, center=TRUE, scale=FALSE), nu=0, nv=0)$d^2
    values[1]/sum(values)
}
knitr::kable(data.frame(Chart=c("Ambient", "Unfolded"),
    `Variance retained by one PC (%)`=100*c(variance.one(seam),
        variance.one(flat$coordinates)), check.names=FALSE), digits=1)

## ----generator----------------------------------------------------------------
scenarios <- c("Seam", "High variance", "Low variance", "Constant", "Nonlinear")
make.example <- function(seed, scenario, n=360L) {
    set.seed(seed)
    t <- runif(n, -.2, .2)
    X <- cbind(ifelse(t < 0, 1, 1-t), ifelse(t < 0, 1+t, 1),
               runif(n, .1, .5), runif(n, .10, .14),
               runif(n, .02, .06), runif(n, .02, .06))
    eta <- switch(scenario,
        "Seam"=-1+2*t/.2,
        "High variance"=-1+2*(X[,3]-.3)/.2,
        "Low variance"=-1+2*(X[,4]-.12)/.02,
        "Constant"=rep(-1,n),
        "Nonlinear"=-1+2.5*((t/.2)^2-1/3))
    probability <- plogis(eta)
    list(X=X, y=rbinom(n,1,probability), probability=probability, t=t)
}

## ----fit-helper---------------------------------------------------------------
methods <- c("Ambient PCA", "Unfolded PCA", "Unfolded PLS")
fit.chart <- function(X, y, query, method, dimension, degree=1L) {
    fit.local.likelihood(X, y, X.eval=query, likelihood.family="bernoulli",
        chart.geometry="linf.simplex", unfold.faces=method != "Ambient PCA",
        coordinate.method=if (method == "Unfolded PLS") "local.pls" else "local.pca",
        support.size=60L, chart.dim=dimension, degree=degree,
        kernel="gaussian", lambda.ridge=1, intercept.mass=.5,
        max.iter=70L, tol=1e-7, fallback="degree0")
}
train <- 1:180
tune <- 181:240
test <- 241:360

## ----single-fit, eval=FALSE---------------------------------------------------
# dat <- make.example(271301, "Low variance")
# fit <- fit.chart(dat$X[train,], dat$y[train], dat$X[test,],
#                  "Unfolded PLS", dimension=1L)
# head(fit$fitted.values)
# fit$diagnostics$status.counts

## ----comparison, eval=recompute-----------------------------------------------
# rows <- list()
# illustrations <- list()
# append.result <- function(scenario, repetition, method, mode, dimension,
#                           prediction, truth, response, fallbacks, tuning.fallbacks=NA_integer_) {
#     stopifnot(all(is.finite(prediction)), all(prediction > 0 & prediction < 1))
#     data.frame(scenario=scenario, repetition=repetition, method=method,
#         mode=mode, dimension=dimension,
#         probability.mse=mean((prediction-truth)^2),
#         brier=mean((prediction-response)^2),
#         log.loss=-mean(response*log(prediction)+(1-response)*log1p(-prediction)),
#         fallbacks=fallbacks, tuning.fallbacks=tuning.fallbacks)
# }
# for (s in seq_along(scenarios)) for (repetition in 1:5) {
#     scenario <- scenarios[s]
#     dat <- make.example(271000+100*s+repetition, scenario)
#     example.predictions <- list()
#     for (method in methods) {
#         fits <- lapply(1:3, function(d) fit.chart(dat$X[train,], dat$y[train],
#             dat$X[c(tune,test),], method, d))
#         tuning.loss <- vapply(fits, function(f)
#             mean((f$fitted.values[1:60]-dat$y[tune])^2), numeric(1))
#         best <- which.min(tuning.loss)
#         for (mode in c("One dimension", "Tuned dimension")) {
#             d <- if (mode == "One dimension") 1L else best
#             chosen <- fits[[d]]
#             prediction <- chosen$fitted.values[61:180]
#             rows[[length(rows)+1L]] <- append.result(scenario, repetition, method,
#                 mode, d, prediction, dat$probability[test], dat$y[test],
#                 sum(chosen$diagnostics$per.eval$fallback[61:180]),
#                 sum(chosen$diagnostics$per.eval$fallback[1:60]))
#             if (repetition == 1L) example.predictions[[paste(method,mode)]] <- prediction
#         }
#     }
#     local <- fit.chart(dat$X[train,], dat$y[train], dat$X[test,],
#                        "Ambient PCA", dimension=1L, degree=0L)
#     constants <- list("Local constant"=local$fitted.values,
#         "Global constant"=rep((sum(dat$y[train])+.5)/(length(train)+1), length(test)))
#     for (method in names(constants)) {
#         rows[[length(rows)+1L]] <- append.result(scenario, repetition, method,
#             "Constant", 0L, constants[[method]], dat$probability[test], dat$y[test], 0L)
#         if (repetition == 1L) example.predictions[[method]] <- constants[[method]]
#     }
#     if (repetition == 1L) illustrations[[scenario]] <- list(data=dat,
#         predictions=example.predictions)
# }
# results <- do.call(rbind, rows)
# stopifnot(nrow(results) == 200L)

## ----stored-comparison, include=FALSE-----------------------------------------
if (!recompute) {
    stored <- readRDS(system.file("extdata", "vignette-results",
                                 "pca-pls-comparison.rds", package="geosmooth"))
    stopifnot(identical(stored$provenance$vignette.sha256,
        digest::digest(file=knitr::current_input(dir=TRUE), algo="sha256")))
    results <- stored$results
    illustrations <- stored$illustrations
}
stopifnot(nrow(results) == 200L, length(illustrations) == 5L,
          all(table(results$scenario, results$repetition) == 8L))

## ----performance-table--------------------------------------------------------
summary.table <- aggregate(cbind(probability.mse,brier) ~ scenario+method+mode,
                           results, mean)
summary.table <- summary.table[order(match(summary.table$scenario, scenarios),
    match(summary.table$mode,c("One dimension","Tuned dimension","Constant")),
    match(summary.table$method,c(methods,"Local constant","Global constant"))),]
knitr::kable(summary.table, digits=4, row.names=FALSE,
    col.names=c("Example","Method","Chart","Probability MSE","Brier loss"))

## ----performance-figure, fig.height=8, fig.cap="Held-out probability error for all five repetitions. Small points are independent repetitions; diamonds are their means. Each panel uses its own labelled vertical scale. Lower is better. The constant baselines are identical in the two columns.", fig.alt="Ten panels compare probability mean squared error across methods, with all repetitions and means shown for fixed and tuned chart dimensions."----
colors <- c("Ambient PCA"="#0072B2", "Unfolded PCA"="#009E73",
    "Unfolded PLS"="#D55E00", "Local constant"="#666666", "Global constant"="#AA4499")
labels <- c("PCA", "Unfolded\nPCA", "Unfolded\nPLS", "Local\nconstant", "Global\nconstant")
old.par <- par(mfrow=c(5,2), mar=c(3.2,4,2,1), oma=c(0,0,1,0), cex=.8)
for (scenario in scenarios) for (mode in c("One dimension","Tuned dimension")) {
    tab <- results[results$scenario==scenario & results$mode %in% c(mode,"Constant"),]
    order.methods <- names(colors)
    plot(NA, xlim=c(.7,5.3), ylim=c(0, max(tab$probability.mse)*1.08+1e-6),
         xaxt="n", xlab="", ylab="Probability MSE", main=paste(scenario,mode,sep=" — "))
    axis(1, at=1:5, labels=labels, cex.axis=.8)
    for (j in 1:5) {
        values <- tab$probability.mse[tab$method==order.methods[j]]
        points(j+seq(-.12,.12,length.out=5), values, pch=16, cex=.7,
               col=adjustcolor(colors[j],alpha.f=.6))
        points(j, mean(values), pch=18, cex=1.3, col=colors[j])
    }
}
par(old.par)

## ----predictions-figure, fig.height=8, fig.cap="Predicted versus known true probabilities for first-repetition test observations. Columns compare unfolded PCA and PLS, each with one dimension. Points near the diagonal recover the true probabilities; all panels share the same axes.", fig.alt="Ten scatterplots show true versus predicted test probabilities for unfolded PCA and unfolded PLS, including constant and nonlinear cases."----
old.par <- par(mfrow=c(5,2), mar=c(3,3.3,2,1), cex=.8)
for (scenario in scenarios) for (method in c("Unfolded PCA","Unfolded PLS")) {
    item <- illustrations[[scenario]]
    prediction <- item$predictions[[paste(method,"One dimension")]]
    plot(item$data$probability[test], prediction, xlim=c(0,1), ylim=c(0,1),
         xlab="True probability", ylab="Predicted probability", pch=16, cex=.65,
         col=adjustcolor(colors[method],alpha.f=.55), main=paste(scenario,method,sep=" — "))
    abline(0,1,lty=2,col="grey35")
}
par(old.par)

## ----slice-figure, fig.height=7.5, fig.cap="Conditional probability slices through the same first-repetition training fits. Nuisance coordinates are fixed at their generating centers. Black is the known truth; colored curves use one-dimensional charts. All probability axes run from zero to one.", fig.alt="Five probability curves compare the known truth with ambient PCA, unfolded PCA and unfolded PLS as the signal coordinate changes."----
slices <- list()
old.par <- par(mfrow=c(3,2), mar=c(4,4,2,1), cex=.85)
for (scenario in scenarios) {
    grid <- seq(-1,1,length.out=60)
    t <- if (scenario %in% c("Seam","Nonlinear","Constant")) .2*grid else rep(.08,60)
    X.slice <- cbind(ifelse(t < 0,1,1-t), ifelse(t < 0,1+t,1), .3,.12,.04,.04)
    if (scenario=="High variance") X.slice[,3] <- .3+.2*grid
    if (scenario=="Low variance") X.slice[,4] <- .12+.02*grid
    eta <- if (scenario=="Constant") rep(-1,60) else if (scenario=="Nonlinear")
        -1+2.5*(grid^2-1/3) else -1+2*grid
    truth <- plogis(eta)
    dat <- illustrations[[scenario]]$data
    prediction <- vapply(methods, function(method)
        fit.chart(dat$X[train,],dat$y[train],X.slice,method,1L)$fitted.values,
        numeric(60))
    axis.value <- if (scenario=="High variance") X.slice[,3] else if (scenario=="Low variance")
        X.slice[,4] else t
    axis.label <- if (scenario=="High variance") "Broad coordinate u3" else if (scenario=="Low variance")
        "Narrow coordinate u4" else "Seam coordinate t"
    plot(axis.value,truth,type="l",lwd=2,col="black",ylim=c(0,1),
         xlab=axis.label,ylab="Probability",main=scenario)
    for (j in seq_along(methods)) lines(axis.value,prediction[,j],col=colors[j],lwd=1.5,lty=j)
    slices[[scenario]] <- list(X=X.slice,truth=truth,prediction=prediction)
}
plot.new()
legend("center",c("Truth",methods),col=c("black",colors[methods]),
       lty=c(1,1,2,3),lwd=c(2,1.5,1.5,1.5),bty="n",cex=.9)
par(old.par)

## ----diagnostic-tables--------------------------------------------------------
tuned <- results[results$mode=="Tuned dimension",]
selection <- as.data.frame(xtabs(~scenario+method+dimension, tuned))
selection <- selection[selection$Freq>0,]
knitr::kable(selection, row.names=FALSE,
    col.names=c("Example","Method","Selected dimension","Repetitions"))
knitr::kable(aggregate(cbind(fallbacks,tuning.fallbacks) ~ scenario+method+mode,
                      results, sum), row.names=FALSE,
    col.names=c("Example","Method","Chart","Test fallbacks","Tuning fallbacks"))

## ----reproduce-comparison, eval=FALSE-----------------------------------------
# options(geosmooth.recompute.vignettes=TRUE)
# rmarkdown::render("vignettes/pca-pls-chart-comparison.Rmd")

## ----save-validation, include=FALSE-------------------------------------------
# A developer may request a machine-readable audit copy; ordinary vignette
# builds do not write outside knitr's output directory.
audit.path <- Sys.getenv("GEOSMOOTH_VIGNETTE_AUDIT", unset="")
if (nzchar(audit.path)) saveRDS(list(results=results, illustrations=illustrations,
    slices=slices, summary=summary.table, session=sessionInfo()), audit.path)

## ----session, echo=FALSE------------------------------------------------------
if (!recompute) {
    stored$provenance$package.versions
} else {
    c(geosmooth=as.character(packageVersion("geosmooth")),
      dgraphs=as.character(packageVersion("dgraphs")))
}

