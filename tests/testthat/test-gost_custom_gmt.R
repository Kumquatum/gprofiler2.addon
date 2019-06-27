context("test-gost_custom_gmt")

path_to_gmt <- system.file("extdata", "h.all.v6.2.entrez.gmt",
                           package = "gprofiler2.addon", mustWork = TRUE)
query <- c("26118", "5837", "6781", "23036", "694", "123", "1466", "7436",
           "23210", "2131", "2152", "5165", "55139", "7360", "229", "8614",
           "54206", "2027", "10957", "3162", "5228", "26330", "9435", "55076")
gost_custom_gmt_res <- gost_custom_gmt(query, custom_gmt = path_to_gmt)


test_that("GMT file checked correctly", {
  expect_error(gost_custom_gmt(query = letters[1:5], custom_gmt = NULL), "custom_gmt should be a string representing GMT file's path")
  expect_error(gost_custom_gmt(query = letters[1:5], custom_gmt = "/fakepath/forunittest/fakeGMT.gmt"), "GMT file path must exist")
  expect_error(gost_custom_gmt(query = letters[1:5], custom_gmt = system.file("DESCRIPTION", package = "gprofiler2.addon",
                                                                              mustWork = TRUE)), "File extension isn't GMT")
  expect_null(gost_custom_gmt(query = letters[1:5], custom_gmt = system.file("extdata", "h.all.v6.2.entrez.gmt", package = "gprofiler2.addon",
                                                                              mustWork = TRUE)))
})


test_that("Output conform", {
  # Case if enrichment has been found
  if (!is.null(gost_custom_gmt_res)) {
    expect_output(str(names(gost_custom_gmt_res)), "result|meta")
    expect_output(str(unlist(gost_custom_gmt_res, recursive = FALSE)), "result.query|result.significant|result.p_value|result.term_size|
                result.query_size|result.intersection_size|result.precision|result.recall|result.term_id|result.source|result.term_name|
                result.effective_domain_size, result.source_order|result.parents|meta.query_metadata|meta.result_metadata|meta.genes_metadata")
  }
})


test_that("Exemple works as before", {
  expect_true(identical(gost_custom_gmt_res$result, gost_custom_gmt_res_ref$result))
  # organism (which is in fact the token returned by the API) and timestamp should be the only ones changing in meta
  expect_true(identical(gost_custom_gmt_res$meta$query_metadata[names(gost_custom_gmt_res$meta$query_metadata) != "organism"],
                        gost_custom_gmt_res_ref$meta$query_metadata[names(gost_custom_gmt_res$meta$query_metadata) != "organism"]))
  expect_true(identical(gost_custom_gmt_res$meta$result_metadata[names(gost_custom_gmt_res$meta$result_metadata) != "timestamp"],
                        gost_custom_gmt_res_ref$meta$result_metadata[names(gost_custom_gmt_res$meta$result_metadata) != "timestamp"]))
})
