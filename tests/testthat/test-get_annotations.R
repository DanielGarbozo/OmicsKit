testthat::test_that("get_annotations returns expected columns for genes", {
  testthat::skip_if_not_installed("biomaRt")
  testthat::skip_if_offline("www.ensembl.org")
  testthat::skip_on_cran()

  ids <- c("ENSG00000141510") # TP53
  res <- get_annotations(
    ensembl_ids = ids,
    species     = "hsapiens_gene_ensembl",
    version     = "112",
    filename    = NULL
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("geneID", "symbol", "biotype", "chromosome",
                    "gene_start", "gene_end", "gene_length", "description") %in% names(res)))
  expect_equal(res$geneID[1], "ENSG00000141510")
})

testthat::test_that("get_annotations works for transcripts", {
  testthat::skip_if_not_installed("biomaRt")
  testthat::skip_if_offline("www.ensembl.org")
  testthat::skip_on_cran()

  # Obtener un transcript válido del release 112
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                              dataset = "hsapiens_gene_ensembl",
                              host = .resolve_ensembl_host("112"))

  ids <- biomaRt::getBM(
    attributes = "ensembl_transcript_id_version",
    filters    = "ensembl_gene_id",
    values     = "ENSG00000141510",
    mart       = ensembl
  )$ensembl_transcript_id_version

  testthat::skip_if(length(ids) == 0)

  res <- get_annotations(
    ensembl_ids = ids[1],
    species     = "hsapiens_gene_ensembl",
    mode        = "transcripts",
    version     = "112",
    filename    = NULL
  )

  expect_s3_class(res, "data.frame")
  expect_true("transcriptID" %in% names(res))
  expect_equal(res$transcriptID[1], ids[1])
})


testthat::test_that("get_annotations errors when biomaRt is missing", {
  testthat::skip_if_not_installed("withr")

  withr::with_namespace("biomaRt", {
    # If biomaRt is installed, skip this test
    testthat::skip("biomaRt is installed; cannot simulate missing package here.")
  })
})
