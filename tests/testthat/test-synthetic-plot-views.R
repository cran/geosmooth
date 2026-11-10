test_that("synthetic plot views are explicit and preserve RNG and datasets", {
    data <- materialize.synthetic(synthetic.spec(
        geometry = dgraphs::synthetic.circle(),
        sampling = dgraphs::synthetic.sampling.grid.interval(0, 2*pi),
        truth = synthetic.truth.polynomial(c(b0 = 0)),
        response = synthetic.response.gaussian(sd = .05)), n = 20, seed = 42)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file)
    on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
    withr::local_seed(17)
    state <- .Random.seed
    before <- serialize(data, NULL)
    for (view in c("auto", "geometry", "response")) {
        for (color in c("truth", "response")) {
            expect_invisible(plot(data, view = view, color = color))
        }
    }
    expect_identical(.Random.seed, state)
    expect_identical(serialize(data, NULL), before)
    expect_error(plot(data, view = "other"), "arg")
    expect_error(plot(data, legend = NA), "TRUE or FALSE")
    if (is.null(data$region)) expect_error(plot(data, color = "region"), "stored region")
    expect_invisible(plot(data, view = "geometry", col = "grey30", pch = 1,
                          xlab = "Horizontal coordinate"))
})
