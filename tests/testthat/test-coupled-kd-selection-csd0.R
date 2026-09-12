test_that("CSD0 candidate table clips auto seeds and records telemetry", {
    auto.seed <- data.frame(
        support.size = c(15L, 25L),
        degree = c(2L, 2L),
        chart.dim = c(25L, 3L)
    )

    tab <- .coupled.kd.candidate.table(
        support.grid = c(15L, 25L),
        degree.grid = 2L,
        chart.dim.grid = c(1L, 2L, "auto", 8L),
        ambient.dim = 20L,
        chart.dim.max = 8L,
        design.margin = 2L,
        auto.chart.dim = auto.seed
    )

    expect_true(all(c("candidate.id", "stage", "support.size",
                      "chart.dim", "chart.dim.source", "chart.dim.raw",
                      "chart.dim.clipped", "chart.dim.seed.clipped",
                      "chart.dim.max", "design.ncol", "design.margin",
                      "feasible", "skip.reason", "reuse.key",
                      "reuse.chart.dim.max", "score", "elapsed.sec") %in%
                    names(tab)))
    auto.15 <- tab[tab$support.size == 15L &
                       tab$chart.dim.source == "auto_seed", , drop = FALSE]
    expect_equal(nrow(auto.15), 1L)
    expect_equal(auto.15$chart.dim.raw, 25L)
    expect_equal(auto.15$chart.dim.clipped, 3L)
    expect_true(auto.15$chart.dim.seed.clipped)
    expect_true(auto.15$feasible)
    expect_equal(auto.15$design.ncol, choose(3L + 2L, 2L))

    auto.25 <- tab[tab$support.size == 25L &
                       tab$chart.dim.source == "auto_seed", , drop = FALSE]
    expect_equal(auto.25$chart.dim.raw, 3L)
    expect_equal(auto.25$chart.dim.clipped, 3L)
    expect_false(auto.25$chart.dim.seed.clipped)
})

test_that("CSD0 candidate table records infeasible pairs instead of dropping them", {
    tab <- .coupled.kd.candidate.table(
        support.grid = 15L,
        degree.grid = 2L,
        chart.dim.grid = c(1L, 8L, "local.auto"),
        ambient.dim = 20L,
        chart.dim.max = 8L,
        design.margin = 2L
    )

    feasible <- tab[tab$chart.dim == "1", , drop = FALSE]
    expect_true(feasible$feasible)
    expect_true(is.na(feasible$skip.reason))

    high.dim <- tab[tab$chart.dim == "8", , drop = FALSE]
    expect_false(high.dim$feasible)
    expect_equal(high.dim$skip.reason, "design_underdetermined")
    expect_equal(high.dim$design.ncol, choose(8L + 2L, 2L))

    local.auto <- tab[tab$chart.dim.source == "local_auto_policy", ,
                      drop = FALSE]
    expect_false(local.auto$feasible)
    expect_equal(local.auto$skip.reason, "local_auto_separate_policy")
    expect_true(is.na(local.auto$reuse.key))
})

test_that("CSD0 duplicate removal is deterministic and favors numeric candidates", {
    tab <- .coupled.kd.candidate.table(
        support.grid = 15L,
        degree.grid = 2L,
        chart.dim.grid = c(1L, 2L, "auto"),
        ambient.dim = 10L,
        chart.dim.max = 8L,
        design.margin = 2L,
        auto.chart.dim = 2L
    )

    dim.2 <- tab[tab$chart.dim.clipped == 2L, , drop = FALSE]
    expect_equal(nrow(dim.2), 1L)
    expect_equal(dim.2$chart.dim.source, "numeric")
    expect_equal(tab$candidate.id, seq_len(nrow(tab)))
})

test_that("CSD0 reuse grouping uses max feasible dimension across degrees", {
    tab <- .coupled.kd.candidate.table(
        support.grid = 25L,
        degree.grid = 1:2,
        chart.dim.grid = c(1L, 2L, 4L),
        ambient.dim = 10L,
        chart.dim.max = 8L,
        design.margin = 2L,
        reuse.type = "weighted"
    )

    feasible <- tab[tab$feasible, , drop = FALSE]
    expect_true(all(feasible$reuse.chart.dim.max == 4L))
    expect_equal(length(unique(feasible$reuse.key)), 1L)
    expect_match(unique(feasible$reuse.key), "25\rgaussian\r4", fixed = TRUE)
})

test_that("CSD0 chart-method reuse keys ignore kernel but keep support and dimension", {
    tab <- .coupled.kd.candidate.table(
        support.grid = 25L,
        degree.grid = 1L,
        chart.dim.grid = c(1L, 3L),
        ambient.dim = 8L,
        kernel.grid = c("gaussian", "tricube"),
        chart.dim.max = 8L,
        design.margin = 2L,
        reuse.type = "chart"
    )

    feasible <- tab[tab$feasible, , drop = FALSE]
    expect_equal(length(unique(feasible$reuse.key)), 1L)
    expect_match(unique(feasible$reuse.key), "25\r3", fixed = TRUE)
    expect_setequal(feasible$kernel, c("gaussian", "tricube"))
})
