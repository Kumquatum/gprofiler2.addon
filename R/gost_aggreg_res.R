#' Aggregation of results from gost_standard and gost_custom_gmt
#'
#' Merging of the output from both \code{\link{gost_standard}} and \code{\link{gost_custom_gmt}} on the same query or list of queries
#'
#' @param standard_output output from \code{\link{gost_standard}}
#' @param custom_output output from \code{\link{gost_custom_gmt}}
#' @return A named list where 'result' contains data.frame with the enrichment analysis results and 'meta' contains metadata needed for Manhattan plot. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query'.
#'  When requesting a 'multi_query', either TRUE or FALSE, the columns of the resulting data frame differ.
#'  If 'evcodes' is set, the return value includes columns 'evidence_codes' and 'intersection'.
#'  The latter conveys info about the intersecting genes between the corresponding query and term.
#' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' @examples
#'
#'
#' @export

gost_aggreg_res <- function(standard_output = NULL, custom_output = NULL){
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

  # Aggregation
  merged_output <- list(
    result = rbind(standard_output$result, custom_output$result),
    meta = "TODO"
    # meta = rbind(standard_output$rmeta, custom_output$meta)
  )

}
