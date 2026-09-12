test_that("GE6 local polynomial design helper preserves MALPS shim behavior", {
    Z <- matrix(
        c(
            0.0, 0.0,
            0.5, 0.2,
            1.0, 0.8,
            1.5, 1.1
        ),
        ncol = 2L,
        byrow = TRUE
    )

    for (degree in 0:2) {
        generic <- .local.polynomial.design.matrix(Z, degree)
        legacy <- .malps.design.matrix(Z, degree)

        expect_equal(generic, legacy)
        expect_equal(nrow(generic), nrow(Z))
        expect_equal(
            ncol(generic),
            length(.local.polynomial.design.column.names(ncol(Z), degree))
        )
        expect_equal(
            .local.polynomial.design.column.names(ncol(Z), degree),
            .malps.design.column.names(ncol(Z), degree)
        )
    }
})

test_that("GE6 local polynomial design helper matches expected low-degree columns", {
    Z <- matrix(c(2, 3, 5, 7), ncol = 2L, byrow = TRUE)
    design <- .local.polynomial.design.matrix(Z, degree = 2L)

    expect_null(colnames(design))
    expect_equal(
        .local.polynomial.design.column.names(2L, 2L),
        c("1", "z1", "z2", "z1.z1", "z1.z2", "z2.z2")
    )
    expect_equal(
        design,
        cbind(
            1,
            Z,
            c(4, 25),
            c(6, 35),
            c(9, 49)
        )
    )
})
