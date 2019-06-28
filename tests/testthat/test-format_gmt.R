context("test-format_gmt")

test_that("GMT file checked correctly", {
  expect_error(format_gmt(path_to_gmt = NULL), "custom_gmt should be a string representing GMT file's path")
  expect_error(format_gmt(path_to_gmt = "/fakepath/forunittest/fakeGMT.gmt"), "GMT file path must exist")
  expect_error(format_gmt(path_to_gmt = system.file("DESCRIPTION", package = "gprofiler2.addon",
                                                                              mustWork = TRUE)), "File extension isn't GMT")
  expect_error(format_gmt(path_to_gmt = system.file("extdata", "h.all.v6.2.entrez.gmt", package = "gprofiler2.addon",
                                                    mustWork = TRUE)), NA)
})
