context("diffcomp")

test_that("diffcomp works on example", {
  expect_equal({
    ## Load example compartment data
    data("compartments")

    ## Define eigenvector data
    eigen <- compartments[,grep("EI", colnames(compartments))]

    ## Define colData describing samples and bioreps
    cData <- data.frame(
      time = c("0000", "0000", "1440", "1440"),
      br = c(1, 2, 1, 2),
      row.names = colnames(eigen)
    )

    ## Define featureData describing rows of data
    fData <- compartments[,1:3]

    ## Calculate differential compartments
    res <- diffcomp(data = eigen,
                    colData = cData,
                    featureData = fData,
                    bioreps = "br",
                    samples = "time",
                    cutoff = 0.01)

    ## Filter for significant compartments
    sigComp <- res$distdf[res$distdf$significant == T,]

    ## Number of significant
    nrow(sigComp)
  }, 13)
})
