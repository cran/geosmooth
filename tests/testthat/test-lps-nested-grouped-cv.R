# E1.10 Part A GATEs -- nested and grouped cross-validation machinery.
#
# Contract: dev/methods/lps/audit_contracts/lps_tiers1to4_contract_2026-06-11.md sB / E1.10
# (sub-item typing per sA1); plan sE1.10 + sec:paired. Spec memo:
# dev/methods/lps/audit_contracts/tiers1to4/e1_10_spec_questions_implementer_2026-06-11.md.
#
# Deterministic GATEs on inline fixtures (no DGP library):
#   E1.10A1  no selection leakage -- the held-out outer fold never enters
#            inner selection (behavioral invariance, per-fold).
#   E1.10A2  cluster integrity -- whole-cluster folds, fold/train cluster
#            disjointness, leave-cluster-out, determinism.
#   E1.10A3  paired discipline -- both arms of a comparison run on the same
#            recorded foldid; fold plumbing cannot be bypassed.
# The E1.10(a)/(b) STUDY decision rules are NOT here; they live in
# dev/methods/lps/ci/e1_10_nested_grouped_cv.R and are gated on the audited DGP
# library.

e110.fixture <- function() {
    n <- 40L
    set.seed(5201)
    X <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
    set.seed(5202)
    y <- sin(2 * pi * X[, 1L]) + 0.5 * X[, 2L]^2 +
        0.1 * stats::rnorm(n)
    list(X = X, y = y, outer.foldid = rep(1:4, length.out = n))
}

e110.fit.args <- function() {
    list(
        support.grid = c(6L, 10L),
        degree.grid = 0:1,
        kernel.grid = c("gaussian", "tricube"),
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
}

e110.cluster.fixture <- function() {
    sizes <- c(1L, 2L, 3L, 4L, 5L, 7L, 8L, 2L, 1L, 3L, 6L, 2L)
    cluster.id <- rep(paste0("c", seq_along(sizes)), times = sizes)
    n <- length(cluster.id)
    set.seed(5301)
    centers <- matrix(stats::runif(length(sizes) * 2L, -1, 1), ncol = 2L)
    X <- centers[match(cluster.id, paste0("c", seq_along(sizes))), ,
                 drop = FALSE]
    set.seed(5302)
    X <- X + 0.1 * matrix(stats::rnorm(n * 2L), ncol = 2L)
    set.seed(5303)
    y <- X[, 1L] + 0.2 * stats::rnorm(n)
    list(X = X, y = y, cluster.id = cluster.id, sizes = sizes)
}

cluster.whole <- function(foldid, cluster.id) {
    all(vapply(split(foldid, cluster.id),
               function(f) length(unique(f)) == 1L,
               logical(1L)))
}

test_that("E1.10A1 no selection leakage: held-out fold y cannot influence inner selection", {
    fx <- e110.fixture()
    base <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L
    )
    for (pos in seq_along(sort(unique(fx$outer.foldid)))) {
        label <- sort(unique(fx$outer.foldid))[[pos]]
        test.idx <- which(fx$outer.foldid == label)
        y.shift <- fx$y
        y.shift[test.idx] <- y.shift[test.idx] + 7

        shifted <- lps.nested.cv(
            X = fx$X, y = y.shift, outer.foldid = fx$outer.foldid,
            fit.args = e110.fit.args(), inner.folds = 3L
        )
        # Fold `label`'s inner selection saw inner-training rows only, so its
        # ENTIRE inner candidate table, its selected configuration, and its
        # held-out predictions are bit-identical under the shift.
        expect_identical(shifted$inner.cv.table[[pos]],
                         base$inner.cv.table[[pos]])
        expect_identical(shifted$folds[pos, c(
            "selected.support.size", "selected.degree",
            "selected.kernel", "selected.bandwidth.multiplier",
            "inner.cv.rmse"
        )], base$folds[pos, c(
            "selected.support.size", "selected.degree",
            "selected.kernel", "selected.bandwidth.multiplier",
            "inner.cv.rmse"
        )])
        expect_identical(shifted$predictions[test.idx],
                         base$predictions[test.idx])
        expect_identical(shifted$inner.foldid[[pos]],
                         base$inner.foldid[[pos]])
    }
})

test_that("E1.10A1 no selection leakage: held-out fold X cannot influence inner selection", {
    fx <- e110.fixture()
    base <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L
    )
    test.idx <- which(fx$outer.foldid == 1L)
    X.moved <- fx$X
    X.moved[test.idx, ] <- X.moved[test.idx, ] + 5

    moved <- lps.nested.cv(
        X = X.moved, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L
    )
    # Moving the held-out rows changes WHERE fold 1 predicts (X.eval), but
    # must not change WHAT fold 1 selected: inner table and selection are
    # functions of inner-training rows only.
    expect_identical(moved$inner.cv.table[[1L]], base$inner.cv.table[[1L]])
    expect_identical(moved$folds[1L, c(
        "selected.support.size", "selected.degree",
        "selected.kernel", "selected.bandwidth.multiplier",
        "inner.cv.rmse"
    )], base$folds[1L, c(
        "selected.support.size", "selected.degree",
        "selected.kernel", "selected.bandwidth.multiplier",
        "inner.cv.rmse"
    )])
})

test_that("E1.10A1 structural: inner folds partition exactly the inner-training rows", {
    fx <- e110.fixture()
    out <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L
    )
    for (pos in seq_along(out$test.index)) {
        test.idx <- out$test.index[[pos]]
        train.idx <- out$train.index[[pos]]
        expect_length(intersect(test.idx, train.idx), 0L)
        expect_identical(sort(c(test.idx, train.idx)),
                         seq_len(nrow(fx$X)))
        # one inner fold id per inner-training row; no entry for held-out rows
        expect_length(out$inner.foldid[[pos]], length(train.idx))
        expect_setequal(unique(out$inner.foldid[[pos]]), 1:3)
    }
})

test_that("E1.10A2 cluster integrity: whole-cluster folds, disjointness, leave-cluster-out", {
    fx <- e110.cluster.fixture()
    v <- 4L
    foldid <- lps.grouped.foldid(fx$cluster.id, v = v)

    expect_true(is.integer(foldid))
    expect_length(foldid, length(fx$cluster.id))
    expect_setequal(unique(foldid), seq_len(v))
    # every cluster wholly inside one fold
    expect_true(cluster.whole(foldid, fx$cluster.id))
    # no cluster appears both in a fold and in its training complement
    for (f in seq_len(v)) {
        held.out.clusters <- unique(fx$cluster.id[foldid == f])
        training.clusters <- unique(fx$cluster.id[foldid != f])
        expect_length(intersect(held.out.clusters, training.clusters), 0L)
    }
    # deterministic without a seed; reproducible with one
    expect_identical(foldid, lps.grouped.foldid(fx$cluster.id, v = v))
    expect_identical(
        lps.grouped.foldid(fx$cluster.id, v = v, shuffle.seed = 11L),
        lps.grouped.foldid(fx$cluster.id, v = v, shuffle.seed = 11L)
    )
    expect_true(cluster.whole(
        lps.grouped.foldid(fx$cluster.id, v = v, shuffle.seed = 11L),
        fx$cluster.id
    ))
    # leave-cluster-out: v = #clusters puts each cluster alone in its fold
    lco <- lps.grouped.foldid(fx$cluster.id, v = length(fx$sizes))
    expect_identical(length(unique(lco)), length(fx$sizes))
    expect_true(cluster.whole(lco, fx$cluster.id))
    counts <- table(lco, fx$cluster.id)
    expect_true(all(colSums(counts > 0) == 1L))
    expect_true(all(rowSums(counts > 0) == 1L))
    # more folds than clusters is an error, not an empty fold
    expect_error(
        lps.grouped.foldid(fx$cluster.id, v = length(fx$sizes) + 1L),
        "between 2 and the number of distinct clusters"
    )
    expect_error(lps.grouped.foldid(c("a", NA, "b"), v = 2L), "missing")
    # fold counts must be whole scalars, never silently truncated (audit P3)
    expect_error(lps.grouped.foldid(fx$cluster.id, v = 2.9), "whole number")
    expect_error(lps.grouped.foldid(fx$cluster.id, v = c(2L, 3L)),
                 "whole number")
    expect_identical(lps.grouped.foldid(fx$cluster.id, v = 4),
                     lps.grouped.foldid(fx$cluster.id, v = 4L))
})

test_that("E1.10A2 cluster integrity: grouped nested CV keeps inner folds cluster-whole", {
    fx <- e110.cluster.fixture()
    outer.foldid <- lps.grouped.foldid(fx$cluster.id, v = 3L)
    out <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L,
        cluster.id = fx$cluster.id, inner.foldid.method = "grouped"
    )
    expect_true(isTRUE(out$outer.cluster.whole))
    for (pos in seq_along(out$inner.foldid)) {
        train.idx <- out$train.index[[pos]]
        expect_true(cluster.whole(out$inner.foldid[[pos]],
                                  fx$cluster.id[train.idx]))
    }
    # a random outer foldid that splits clusters is reported as such
    split.foldid <- rep(1:3, length.out = length(fx$cluster.id))
    out.split <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = split.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L,
        cluster.id = fx$cluster.id, inner.foldid.method = "grouped"
    )
    expect_false(isTRUE(out.split$outer.cluster.whole))
})

test_that("E1.10A3 paired discipline: both arms run on the same recorded foldid", {
    fx <- e110.fixture()
    out <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L
    )
    # The selected-min arm used exactly the outer foldid (same object value).
    expect_identical(out$selected.min$foldid, out$outer.foldid)
    expect_identical(out$selected.min$fit$foldid, out$outer.foldid)
    expect_identical(out$outer.foldid, as.integer(fx$outer.foldid))
    # Each inner fit.lps recorded exactly the constructed inner foldid.
    for (pos in seq_along(out$inner.foldid)) {
        expect_identical(out$inner.foldid.used[[pos]],
                         out$inner.foldid[[pos]])
    }
    # Fold plumbing cannot be bypassed through fit.args.
    for (bad in c("foldid", "X", "y", "X.eval")) {
        bad.args <- e110.fit.args()
        bad.args[[bad]] <- if (identical(bad, "foldid")) {
            rep(1:2, length.out = nrow(fx$X))
        } else {
            fx$X
        }
        expect_error(
            lps.nested.cv(
                X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
                fit.args = bad.args, inner.folds = 3L
            ),
            "must not contain"
        )
    }
    # Unnamed fit.args are rejected (no positional smuggling).
    expect_error(
        lps.nested.cv(
            X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
            fit.args = list(10L), inner.folds = 3L
        ),
        "named list"
    )
    # Non-gaussian outcome families are out of scope for this machinery.
    expect_error(
        lps.nested.cv(
            X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
            fit.args = c(e110.fit.args(), list(outcome.family = "binomial")),
            inner.folds = 3L
        ),
        "gaussian"
    )
    # inner.folds must be a whole scalar, never silently truncated (audit P3)
    expect_error(
        lps.nested.cv(
            X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
            fit.args = e110.fit.args(), inner.folds = 2.9
        ),
        "whole number"
    )
    expect_error(
        lps.nested.cv(
            X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
            fit.args = e110.fit.args(), inner.folds = c(2L, 3L)
        ),
        "whole number"
    )
})

test_that("E1.10A3 paired discipline: randomized inner folds are seeded, recorded, reproducible", {
    fx <- e110.fixture()
    one <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L,
        inner.shuffle.seed = 71L
    )
    two <- lps.nested.cv(
        X = fx$X, y = fx$y, outer.foldid = fx$outer.foldid,
        fit.args = e110.fit.args(), inner.folds = 3L,
        inner.shuffle.seed = 71L
    )
    expect_identical(one$inner.foldid, two$inner.foldid)
    expect_identical(one$predictions, two$predictions)
    expect_identical(one$inner.shuffle.seed, 71L)
    # balanced multiset preserved by the shuffle
    for (pos in seq_along(one$inner.foldid)) {
        expect_identical(sort(table(one$inner.foldid[[pos]])),
                         sort(table(two$inner.foldid[[pos]])))
        expect_setequal(unique(one$inner.foldid[[pos]]), 1:3)
    }
})
