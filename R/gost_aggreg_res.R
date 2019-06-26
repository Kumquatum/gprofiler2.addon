#' Aggregation of results from gost_standard and gost_custom_gmt
#'
#' Merging of the output from both \code{\link[gprofiler2]{gost}} and \code{\link{gost_custom_gmt}} on the same query or list of queries
#' May also be used on an output from this function with another output from \code{\link{gost_custom_gmt}}
#'
#' @param standard_output output from \code{\link[gprofiler2]{gost}}
#' @param custom_output output from \code{\link{gost_custom_gmt}}
#' @param check_content boolean specifying if content must be checked to be from same queries
#' @return A named list where 'result' contains data.frame with the enrichment analysis results and 'meta' contains metadata needed for Manhattan plot. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query'.
#'  When requesting a 'multi_query', either TRUE or FALSE, the columns of the resulting data frame differ.
#'  If 'evcodes' is set, the return value includes columns 'evidence_codes' and 'intersection'.
#'  The latter conveys info about the intersecting genes between the corresponding query and term.
#' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' @examples
#' path_to_gmt <- system.file("extdata", "h.all.v6.2.entrez.gmt",
#'                            package = "gprofiler2.addon", mustWork = TRUE)
#' query <- c("26118", "5837", "6781", "23036", "694", "123", "1466", "7436",
#'            "23210", "2131", "2152", "5165", "55139", "7360", "229", "8614",
#'            "54206", "2027", "10957", "3162", "5228", "26330", "9435", "55076")
#' gost_custom_gmt_res <- gost_custom_gmt(query, custom_gmt = path_to_gmt)
#' gost_classic_res <- gprofiler2::gost(query, numeric_ns = "ENTREZGENE_ACC")
#' gost_aggregated <- gost_aggreg_res(gost_classic_res, gost_custom_gmt_res)
#'
#' @import gprofiler2
#'
#' @export

gost_aggreg_res <- function(standard_output = NULL, custom_output = NULL, check_content = TRUE){
  # Checking format
  lapply(list(standard_output, custom_output), function(output){
    # Checking first level of output structure
    if (!("result" %in% names(output))) stop("Name 'result' not found from the input")
    if (!("meta" %in% names(output))) stop("Name 'meta' not found from the input")

    # Checking second level of output structure
    sub_names <- c("result.query", "result.significant", "result.p_value", "result.term_size", "result.query_size",
                   "result.intersection_size", "result.precision", "result.recall", "result.term_id", "result.source",
                   "result.term_name", "result.effective_domain_size", "result.source_order", "meta.query_metadata",
                   "meta.result_metadata", "meta.genes_metadata")
    lapply(sub_names, function(name){
      if (!(name %in% names(unlist(output, recursive = FALSE)))) stop(paste0("Name ", name, " not found from the input"))
    })
  })

  # Checking content is from same queries
  if (check_content) {
    meta_standard <- standard_output$meta
    meta_custom <- custom_output$meta

    lapply(names(meta_standard$query_metadata), function(name) {
      if (!(name %in% c("organism", "sources"))) {
        if (name == "queries") {
          if (length(meta_standard$query_metadata$queries) != length(meta_custom$query_metadata$queries)) stop("Number of queries different.")
          if (length(setdiff(meta_standard$query_metadata$queries, meta_custom$query_metadata$queries)) > 0) stop("Length of queries different.")
          if (mapply(function(x,y){all(x == y)}, meta_standard$query_metadata$queries, meta_custom$query_metadata$queries)) {
            warning("Different queries IDs. It may be due to different type of ID (Ensembl, Entrez, etc.). IDs from standard output will be kept.")
          }
        } else if (name == "numeric_ns") {
          if (meta_standard$query_metadata$numeric_ns != meta_custom$query_metadata$numeric_ns)
            warning("Different type of IDs. ", meta_standard$query_metadata$numeric_ns, " (standard) and ", meta_custom$query_metadata$numeric_ns, " (custom)")
        } else {
          # if (meta_standard$query_metadata[[name]] != meta_custom$query_metadata[[name]]) stop(paste0("Item ", name, " of query_metadata isn't identical"))
          if (!all.equal(meta_standard$query_metadata[[name]], meta_custom$query_metadata[[name]])) stop(paste0("Item ", name, " of query_metadata isn't identical"))
        }
      }
    })
  }

  # Merged meta based on standard output
  meta_merged <- meta_standard
  meta_merged$query_metadata$sources <- c(meta_merged$query_metadata$sources, meta_custom$query_metadata$sources)
  meta_merged$result_metadata[[names(meta_custom$result_metadata[2])]] <- meta_custom$result_metadata[[2]]
  meta_merged$result_metadata$timestamp <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", "GMT") # Note : not exactly the same format as the one from API

  # Aggregation
  merged_output <- list(
    result = rbind(standard_output$result, custom_output$result),
    meta = meta_merged
  )

  return(merged_output)
}
