context("test-gost_aggreg_res")

query <- c("26118", "5837", "6781", "23036", "694", "123", "1466", "7436",
           "23210", "2131", "2152", "5165", "55139", "7360", "229", "8614",
           "54206", "2027", "10957", "3162", "5228", "26330", "9435", "55076")
gost_classic_res <- gprofiler2::gost(query, numeric_ns = "ENTREZGENE_ACC")
aggreg_res = gost_aggreg_res(gost_classic_res, gost_custom_gmt_res_ref)


test_that("Queries content", {
  expect_error(gost_aggreg_res(gprofiler2::gost(list(q1 = query, q2 = query), numeric_ns = "ENTREZGENE_ACC"), gost_custom_gmt_res_ref), "Number of queries different.")
  expect_error(gost_aggreg_res(gprofiler2::gost(c(query, "18325", "2452"), numeric_ns = "ENTREZGENE_ACC"), gost_custom_gmt_res_ref), "Length of queries different.")
  fake_gost_classic_res <- gost_classic_res
  fake_gost_classic_res$meta$query_metadata$ordered  <- TRUE
  expect_error(gost_aggreg_res(fake_gost_classic_res, gost_custom_gmt_res_ref), "Item ordered of query_metadata isn't identical")
})


test_that("Output conform", {
  # Case if enrichment has been found
  if (!is.null(aggreg_res)) {
    expect_output(str(names(aggreg_res)), "result|meta")
    expect_output(str(unlist(aggreg_res, recursive = FALSE)), "result.query|result.significant|result.p_value|result.term_size|
                result.query_size|result.intersection_size|result.precision|result.recall|result.term_id|result.source|result.term_name|
                result.effective_domain_size, result.source_order|result.parents|meta.query_metadata|meta.result_metadata|meta.genes_metadata")
  }
})

test_that("Able to plot gostplot", {
  expect_error(gprofiler2::gostplot(aggreg_res), NA)
})
