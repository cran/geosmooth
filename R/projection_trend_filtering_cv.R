#' Cross-Validate Projection Trend Filtering
#'
#' Select explicit mean-loss penalties and component counts using supplied
#' folds. Every candidate in every training split repeats preprocessing,
#' direction selection, refinement and joint fitting.
#'
#' @inheritParams fit.projection.trend.filtering
#' @param foldid Integer vector, or matrix with one column per CV repeat. Each
#'   column must contain at least two folds; labels need not be consecutive.
#' @param nterms.grid Nonnegative integer component counts. In fixed mode the
#'   nonzero count must match the supplied direction matrix.
#' @param lambda.grid Nonnegative finite penalties or `Inf`.
#' @param curve.ridge.grid Nonnegative finite component ridge coefficients.
#' @param groups Optional participant/group labels. All rows of a group must
#'   belong to the same fold in each repeat. Groups validate supplied splits;
#'   they do not automatically change observation weights.
#' @param keep.fits Retain training-split fits for inspection, default `FALSE`.
#'
#' @details Scores pool weighted held-out losses within each repeat and then
#'   average repeats. Gaussian selection uses MSE; binomial selection uses
#'   stable log loss from logits. Every required fold must converge and
#'   predict every held-out row. Fixed and pattern-search modes select the
#'   first grid row on equal scores.
#'   Gauss--Newton modes also require explicit search eligibility and resolve
#'   equal CV scores by fewer requested components, stronger penalty, then
#'   grid order. Their random streams are keyed by repeat and fold.
#'   Penalties already use mean loss, so no additional fold-mass multiplier
#'   is applied by this wrapper. A final full-data refit is performed only for
#'   the selected candidate. Its failure is reported without substituting
#'   another candidate. CV scores used for selection are not an independent
#'   assessment of the selected procedure; use outer folds or a test set.
#'
#' @return A `projection.trend.filtering.cv` containing `cv.table`, the
#'   candidate-by-split `fold.results`, held-out `predictions` (logits for
#'   binomial data), `foldid`, optional `fold.fits`, `selected.index`, `fit`,
#'   and `status`. `selected.index` is NULL if validation or the full refit
#'   failed; `cv.selected.index` retains the validation choice when applicable.
#'   For Gauss--Newton modes, `fold.results` separately records
#'   `numerical.certified`, `search.eligible`, `cv.eligible`,
#'   `search.status`, requested/effective component counts and native-fit
#'   counts, including when `keep.fits=FALSE`.
#' @examples
#' X <- matrix(seq(-1, 1, length.out = 24), ncol = 1)
#' cv <- cv.projection.trend.filtering(X, sin(2*X[,1]),
#'     foldid = rep(1:3, 8), nterms.grid = 1,
#'     lambda.grid = c(Inf, 0.001), direction.search = "fixed")
#' cv$cv.table
#' @export
cv.projection.trend.filtering <- function(
    X, y, foldid, lambda.grid, nterms.grid = c(0L, 1L, 2L),
    curve.ridge.grid = 0.001, family = c("gaussian", "binomial"),
    directions = NULL, direction.search = c("learned", "fixed", "gn", "pls"),
    weights = NULL, offset = NULL, groups = NULL,
    preprocessing = c("scale", "isotropic", "none"), seed = 1L,
    control = list(), backend = c("auto", "dense", "banded"),
    keep.fits = FALSE, warn = TRUE) {
    call <- match.call(); family <- match.arg(family)
    direction.search <- match.arg(direction.search)
    if (missing(preprocessing) && direction.search %in% c("gn", "pls")) preprocessing <- "isotropic"
    preprocessing <- match.arg(preprocessing); backend <- match.arg(backend)
    .ptf.logical(keep.fits,"keep.fits"); .ptf.logical(warn,"warn")
    control <- .ptf.control(control); .ptf.scalar(seed,"seed",integer=TRUE)
    for (nm in c("lambda.grid","nterms.grid","curve.ridge.grid")) {
        grid <- get(nm)
        if (!is.numeric(grid) || !length(grid) || !is.null(dim(grid)))
            stop(nm," must be a nonempty numeric vector.",call.=FALSE)
        for (value in grid) .ptf.scalar(value,nm,integer=nm=="nterms.grid",infinite=nm=="lambda.grid")
    }
    checked <- .ptf.problem(X,y,weights,offset,family,preprocessing,0,0,control,backend)
    X <- checked$X; weights <- checked$weights; offsets <- checked$offset; n <- nrow(X)
    if (!is.numeric(foldid) || any(!is.finite(foldid)) || any(foldid!=floor(foldid)) ||
        (is.null(dim(foldid)) && length(foldid)!=n) ||
        (!is.null(dim(foldid)) && (!is.matrix(foldid) || nrow(foldid)!=n)))
        stop("foldid must be an integer vector or matrix with one row per observation.",call.=FALSE)
    folds <- matrix(foldid,nrow=n)
    if (!ncol(folds) || any(vapply(seq_len(ncol(folds)),function(r) length(unique(folds[,r]))<2L,logical(1))))
        stop("Each repeat needs at least two folds.",call.=FALSE)
    if (!is.null(groups)) {
        if (!is.atomic(groups) || !is.null(dim(groups)) || length(groups)!=n || anyNA(groups))
            stop("groups must be nonmissing labels, one per observation.",call.=FALSE)
        rows <- split(seq_len(n),groups)
        for (r in seq_len(ncol(folds))) if (any(vapply(rows,function(ix)
            length(unique(folds[ix,r]))!=1L,logical(1))))
            stop("Every group must stay within one fold in each repeat.",call.=FALSE)
    }
    candidates <- unique(expand.grid(nterms=nterms.grid,lambda=lambda.grid,
        curve.ridge=curve.ridge.grid,KEEP.OUT.ATTRS=FALSE))
    candidates$candidate <- seq_len(nrow(candidates))
    predictions <- lapply(seq_len(ncol(folds)),function(r) matrix(NA_real_,n,nrow(candidates)))
    retained <- if (keep.fits) vector("list",nrow(candidates)) else NULL
    ledger <- list()
    for (candidate in seq_len(nrow(candidates))) {
        spec <- candidates[candidate,]
        args <- list(family=family,directions=if (spec$nterms==0) NULL else directions,
            direction.search=direction.search,nterms=spec$nterms,lambda=spec$lambda,
            curve.ridge=spec$curve.ridge,preprocessing=preprocessing,seed=seed,
            control=control,backend=backend,warn=FALSE)
        if (keep.fits) retained[[candidate]] <- list()
        for (r in seq_len(ncol(folds))) for (fold in sort(unique(folds[,r]))) {
            test <- folds[,r]==fold; train <- !test
            if (direction.search %in% c("gn", "pls"))
                args$control$stream.key <- c(control$stream.key, "repeat", as.character(r),
                                            "fold", as.character(fold))
            fit <- tryCatch(do.call(fit.projection.trend.filtering,c(list(
                X=X[train,,drop=FALSE],y=y[train],weights=weights[train],
                offset=if (is.null(offset)) NULL else offsets[train]),args)),error=.ptf.numerical.error)
            if (keep.fits) retained[[candidate]][[paste(r,fold,sep=":")]] <- fit
            status <- if (inherits(fit,"error")) conditionMessage(fit) else fit$status
            valid <- !inherits(fit,"error") && .ptf.eligible(fit)
            loss.sum <- NA_real_
            if (valid) {
                pred <- tryCatch(as.numeric(stats::predict(fit,X[test,,drop=FALSE],
                    offset=if (is.null(offset)) NULL else offsets[test],type="link")),error=.ptf.numerical.error)
                valid <- !inherits(pred,"error") && length(pred)==sum(test) && all(is.finite(pred))
                if (valid) {
                    predictions[[r]][test,candidate] <- pred
                    loss.sum <- sum(weights[test]*.ptf.loss(y[test],pred,family)) *
                        if (family=="gaussian") 2 else 1
                } else status <- if (inherits(pred,"error")) conditionMessage(pred) else "nonfinite_predictions"
            }
            row <- data.frame(candidate=candidate,repeat.id=r,fold=fold,
                training.mass=sum(weights[train]),validation.mass=sum(weights[test]),
                loss.sum=loss.sum,converged=valid,status=status)
            if (direction.search %in% c("gn", "pls")) {
                failed <- inherits(fit,"error")
                row$numerical.certified <- !failed && isTRUE(fit$checks$converged)
                row$search.eligible <- !failed && isTRUE(fit$cv.eligible)
                row$cv.eligible <- valid
                row$search.status <- if (failed) status else fit$search.status
                row$requested.nterms <- spec$nterms
                row$effective.nterms <- if (failed) NA_integer_ else fit$nterms
                row$native.fits <- if (failed) NA_real_ else fit$counters$native.fits
            }
            ledger[[length(ledger)+1L]] <- row
        }
    }
    ledger <- do.call(rbind,ledger)
    candidates$cv.loss <- vapply(candidates$candidate,function(j) {
        take <- ledger$candidate==j
        if (!all(ledger$converged[take])) return(NA_real_)
        sum(ledger$loss.sum[take])/(sum(weights)*ncol(folds))
    },numeric(1))
    candidates$eligible <- is.finite(candidates$cv.loss)
    choice <- if (any(candidates$eligible)) which.min(ifelse(candidates$eligible,candidates$cv.loss,Inf)) else NULL
    if (!is.null(choice) && direction.search %in% c("gn", "pls")) {
        take <- which(candidates$eligible)
        choice <- take[order(candidates$cv.loss[take], candidates$nterms[take],
                             -candidates$lambda[take], take)[1L]]
    }
    fit <- NULL; selected <- choice; status <- "no_eligible_candidate"
    if (!is.null(choice)) {
        spec <- candidates[choice,]
        fit <- tryCatch(fit.projection.trend.filtering(X,y,spec$lambda,family=family,
            directions=if (spec$nterms==0) NULL else directions,direction.search=direction.search,
            nterms=spec$nterms,weights=weights,offset=offset,curve.ridge=spec$curve.ridge,
            preprocessing=preprocessing,seed=seed,control=control,backend=backend,warn=FALSE),error=.ptf.numerical.error)
        if (inherits(fit,"error") || !.ptf.eligible(fit)) {
            selected <- NULL; status <- "full_refit_failed"
        } else status <- "complete"
    }
    if (warn && status!="complete") warning("Projection CV: ",status,"; inspect fold.results and fit.",call.=FALSE)
    structure(list(call=call,family=family,cv.table=candidates,fold.results=ledger,
        predictions=predictions,foldid=folds,groups=groups,fold.fits=retained,
        selected.index=selected,cv.selected.index=choice,fit=fit,status=status),
        class="projection.trend.filtering.cv")
}

#' @rdname predict.projection.trend.filtering.fit
#' @method predict projection.trend.filtering.cv
#' @export
predict.projection.trend.filtering.cv <- function(object, newdata=NULL, ...) {
    if (is.null(object$selected.index)) stop("No validated full-data CV fit is available.",call.=FALSE)
    stats::predict(object$fit,newdata=newdata,...)
}
