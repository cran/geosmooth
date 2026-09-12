make.odcv4.path.graph <- function(lengths) {
    n <- length(lengths) + 1L
    adj <- vector("list", n)
    wt <- vector("list", n)
    for (i in seq_along(lengths)) {
        j <- i + 1L
        ell <- as.double(lengths[[i]])
        adj[[i]] <- c(adj[[i]], j)
        wt[[i]] <- c(wt[[i]], ell)
        adj[[j]] <- c(adj[[j]], i)
        wt[[j]] <- c(wt[[j]], ell)
    }
    list(
        adj.list = lapply(adj, as.integer),
        weight.list = lapply(wt, as.double)
    )
}

expect.odcv4.visit.fit <- function(fit, n.visits, n.candidates) {
    expect_s3_class(fit, "density_fit")
    expect_identical(fit$method.id, "graph_random_walk")
    expect_identical(fit$theta$od.cv, "visit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), n.candidates)
    expect_equal(dim(fit$visit.cv.predicted.mass),
                 c(n.visits, n.candidates))
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_equal(
        fit$diagnostics$od.visit.cv.selection$visit.cv.neg.log.rho,
        min(fit$visit.cv.table$visit.cv.neg.log.rho),
        tolerance = 1e-12
    )
}

test_that("OD-CV4 graph random-walk visit CV searches graph-control candidates", {
    n <- 8L
    X <- cbind(seq(0, 1, length.out = n), 0)
    graph <- make.odcv4.path.graph(c(1, 2, 1, 1, 3, 1, 2))
    subject.index <- c(2L, 2L, 3L, 5L, 6L, 6L, 7L, 8L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "graph_random_walk",
        graph = graph,
        graph.control = list(
            walk.step.grid = c(0L, 1L, 2L),
            affinity.method.grid = c("exp_neg_length_over_median",
                                     "inverse_length"),
            affinity.scale.grid = c(NA_real_, 0.75),
            affinity.epsilon.grid = 1e-8,
            normalize.grid = TRUE
        ),
        od.cv = "visit",
        visit.foldid = rep(1:4, length.out = length(subject.index))
    )

    expect.odcv4.visit.fit(fit, length(subject.index), 9L)
    expect_true(all(c("walk.step", "affinity.method", "affinity.scale",
                      "affinity.epsilon", "normalize") %in%
                    names(fit$visit.cv.table)))
    expect_setequal(fit$visit.cv.table$walk.step, c(0L, 1L, 2L))
    expect_setequal(fit$visit.cv.table$affinity.method,
                    c("exp_neg_length_over_median", "inverse_length"))
    inverse.rows <- fit$visit.cv.table$affinity.method == "inverse_length"
    expect_true(all(is.na(fit$visit.cv.table$affinity.scale[inverse.rows])))
    expect_true(fit$theta$selected.walk.step %in% c(0L, 1L, 2L))
    expect_true(fit$theta$affinity.method %in%
                c("exp_neg_length_over_median", "inverse_length"))
    expect_equal(
        fit$theta$selected.walk.step,
        fit$diagnostics$od.visit.cv.selection$walk.step
    )
    expect_identical(
        fit$theta$affinity.method,
        fit$diagnostics$od.visit.cv.selection$affinity.method
    )
})

test_that("OD-CV4 graph random-walk visit CV keeps compact selected metadata", {
    n <- 7L
    X <- matrix(seq(0, 1, length.out = n), ncol = 1L)
    graph <- make.odcv4.path.graph(rep(1, n - 1L))
    subject.index <- c(2L, 3L, 3L, 5L, 6L, 7L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "graph_random_walk",
        graph = graph,
        graph.control = list(walk.step.grid = c(0L, 1L)),
        od.cv = "visit",
        visit.foldid = rep(1:3, length.out = length(subject.index)),
        return.details = FALSE
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.list(fit$diagnostics$od.visit.cv))
    expect_true(is.list(fit$diagnostics$od.visit.cv.selection))
    expect_null(fit$visit.cv.table)
    expect_null(fit$visit.foldid)
    expect_null(fit$visit.cv.predicted.mass)
})

test_that("OD-CV4 graph random-walk visit CV validates graph-control grids", {
    X <- matrix(seq(0, 1, length.out = 6), ncol = 1L)
    graph <- make.odcv4.path.graph(rep(1, 5L))
    subject.index <- c(1L, 2L, 4L, 6L)

    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "graph_random_walk",
            graph = graph,
            od.cv = "visit",
            visit.cv.folds = 2L,
            walk.step.grid = c(0L, 1L)
        ),
        "graph.control"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "graph_random_walk",
            graph = graph,
            graph.control = list(walk.step.grid = -1L),
            od.cv = "visit",
            visit.cv.folds = 2L
        ),
        "walk.step.grid"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "graph_random_walk",
            graph = graph,
            graph.control = list(affinity.method.grid = "bad"),
            od.cv = "visit",
            visit.cv.folds = 2L
        ),
        "unsupported"
    )
})
