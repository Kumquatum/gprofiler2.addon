#' Gene list functional enrichment with custom gmt file
#'
#' Interface to the g:Profiler tool g:GOSt for functional enrichment of a list of gene IDs based on a given GMT file.
#' In case the input 'query' is a list, results for multiple queries will be returned in the same data frame with column 'query' indicating the corresponding query name.
#' If 'multi_query' is selected, the result is a data frame for comparing multiple input lists,
#' just as in the web tool.
#'
#' @param query vector that can consist of IDs of the same type as those inside gmt file; or a (named) list of such vectors. If not, use \code{\link{gprofiler2::gconvert}} to get required IDs.
#' family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' @param ordered_query in case input gene lists are ranked this option may be
#'  used to get GSEA style p-values.
#' @param multi_query in case of multiple gene lists, returns comparison table of these lists.
#' If enabled, the result data frame has columns named 'p_values', 'query_sizes', 'intersection_sizes' with vectors showing values in the order of input queries.
#' @param significant whether all or only statistically significant results should
#'  be returned.
#' @param exclude_iea exclude GO electronic annotations (IEA).
#' @param measure_underrepresentation measure underrepresentation.
#' @param evcodes include evidence codes to the results. Note
#'  that this can decrease performance and make the query slower.
#'  In addition, a column 'intersection' is created that contains the gene id-s that intersect between the query and term.
#' @param user_threshold custom p-value threshold, results with a larger p-value are
#'  excluded.
#' @param correction_method the algorithm used for multiple testing correction, one of "gSCS" (synonyms: "analytical", "g_SCS"), "fdr" (synonyms: "false_discovery_rate"), "bonferroni".
#' @param domain_scope how to define statistical domain, one of "annotated", "known".
#' @param custom_bg vector of gene names to use as a statistical background. If given, the domain_scope is set to 'custom'.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param custom_gmt A path to the gmt file.
#' @return A named list where 'result' contains data.frame with the enrichment analysis results and 'meta' contains metadata needed for Manhattan plot. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query'.
#'  When requesting a 'multi_query', either TRUE or FALSE, the columns of the resulting data frame differ.
#'  If 'evcodes' is set, the return value includes columns 'evidence_codes' and 'intersection'.
#'  The latter conveys info about the intersecting genes between the corresponding query and term.
#' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' @examples
#' # Use GMT file with entrez gene IDs or download gmt file at
#' software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/h.all.v6.2.symbols.gmt
#' path_to_gmt <- "path to the gmt file or the downloaded one"
#' gost_custom_gmt_res <- gost_custom_gmt(c("26118", "5837", "6781", "23036", "694", "123", "1466", "7436",
#'                                          "23210", "2131", "2152", "5165", "55139", "7360", "229", "8614",
#'                                          "54206", "2027", "10957", "3162", "5228", "26330", "9435", "55076"),
#'                                          custom_gmt = path_to_gmt)
#'
#' @export

gost_custom_gmt = function(query, ordered_query = FALSE,
                           multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                           measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05,
                           correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate",
                                                 "gSCS", "analytical"), domain_scope = c("annotated", "known", "custom"),
                           custom_bg = NULL, numeric_ns = "", custom_gmt = NULL)
{
  # Variables that need to be defined here because they are globals into the package
  gp_globals = new.env()
  gp_globals$version =
    tryCatch(
      utils::packageVersion("gprofiler2"),
      error = function(e) { return("unknown_version") }
    );
  # Set SSL version to TLSv1_2 with fallback to TLSv1_1
  # CURL_SSLVERSION_SSLv3 is not used due to the SSLv3 vulnerability <https://access.redhat.com/articles/1232123>
  # CURL_SSLVERSION_TLSv1_3 is not widespread enough to have a built-in LibreSSL support yet.
  # (curl's authors may decide to change it at some point, so links to the source are provided.)
  gp_globals$CURL_SSLVERSION_TLSv1_1 <- 5L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1925>
  gp_globals$CURL_SSLVERSION_TLSv1_2 <- 6L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1926>
  gp_globals$rcurl_opts =
    RCurl::curlOptions(useragent = paste("gprofiler2/", gp_globals$version, sep=""), sslversion = gp_globals$CURL_SSLVERSION_TLSv1_2)
  gp_globals$base_url = "http://biit.cs.ut.ee/gprofiler"

  # Check gmt file
  if (!is.character(custom_gmt)) {stop("custom_gmt should be a string representing GMT file's path")}
  if (!file.exists(custom_gmt)) {stop("GMT file path must exist")}
  if (!grepl("\\.gmt$", custom_gmt, ignore.case = TRUE)) {stop("File extension isn't GMT")}
  gmt = readChar(custom_gmt, file.info(custom_gmt)$size)

  # Usual checks in the original function gost
  if (is.null(query)) {
    stop("Missing query")
  }
  else if (is.list(query)) {
    if (is.data.frame(query)) {
      stop("Query can't be a data.frame. Please use a vector or list of identifiers.")
    }
    qnames = names(query)
    if (is.null(qnames)) {
      qnames = paste("query", seq(1, length(query)), sep = "_")
      names(query) = qnames
    }
    query = lapply(query, function(x) x[!is.na(x)])
  }
  else {
    query = query[!is.na(query)]
  }
  correction_method <- match.arg(correction_method)
  if (!is.null(custom_bg)) {
    if (!is.vector(custom_bg)) {
      stop("custom_bg must be a vector")
    }
    message("Detected custom background input, domain scope is set to 'custom'")
    domain_scope <- "custom"
    t <- ifelse(length(custom_bg) == 1, custom_bg <- jsonlite::unbox(custom_bg),
                custom_bg <- custom_bg)
  }

  # Upload gmt and get token needed for second query
  url_custom <- paste0(file.path(gp_globals$base_url, "api", "gost", "custom"), "/")
  body_custom <- jsonlite::toJSON(list(gmt = jsonlite::unbox(gmt), name = jsonlite::unbox(basename(custom_gmt)),
                                       output = jsonlite::unbox("json")), auto_unbox = FALSE, null = "null")
  headers_custom <- list(Accept = "application/json", `Content-Type` = "application/json",
                         charset = "UTF-8")
  h1_custom = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2_custom = RCurl::getCurlHandle()

  r_custom = RCurl::curlPerform(url = url_custom, postfields = body_custom, httpheader = headers_custom,
                                customrequest = "POST", verbose = FALSE, ssl.verifypeer = FALSE,
                                writefunction = h1_custom$update, curl = h2_custom, .opts = gp_globals$rcurl_opts)
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2_custom)[["response.code"]]
  txt <- h1_custom$value()
  if (rescode != 200) {
    stop("GMT upload - Bad request, response code ", rescode)
  }
  res <- jsonlite::fromJSON(txt)
  token <- res$organism

  # Usual query launched with token
  url = paste0(file.path(gp_globals$base_url, "api", "gost", "profile"), "/")
  domain_scope <- match.arg(domain_scope)
  body <- jsonlite::toJSON((list(organism = jsonlite::unbox(token),
                                 query = query, user_threshold = jsonlite::unbox(user_threshold),
                                 all_results = jsonlite::unbox(!significant), ordered = jsonlite::unbox(ordered_query),
                                 no_evidences = jsonlite::unbox(!evcodes), combined = jsonlite::unbox(multi_query),
                                 measure_underrepresentation = jsonlite::unbox(measure_underrepresentation),
                                 no_iea = jsonlite::unbox(!exclude_iea), domain_scope = jsonlite::unbox(domain_scope),
                                 numeric_ns = jsonlite::unbox(numeric_ns), significance_threshold_method = jsonlite::unbox(correction_method),
                                 background = custom_bg, output = jsonlite::unbox("json"))),
                           auto_unbox = FALSE, null = "null")
  headers <- list(Accept = "application/json", `Content-Type` = "application/json",
                  charset = "UTF-8")
  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2 = RCurl::getCurlHandle()
  r = RCurl::curlPerform(url = url, postfields = body, httpheader = headers,
                         customrequest = "POST", verbose = FALSE, ssl.verifypeer = FALSE,
                         writefunction = h1$update, curl = h2, .opts = gp_globals$rcurl_opts)
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()
  if (rescode != 200) {
    stop("Results - Bad request, response code ", rescode)
  }
  res <- jsonlite::fromJSON(txt)
  df = res$result
  meta = res$meta
  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the organism is correct or set significant = FALSE")
    return(NULL)
  }
  if (multi_query) {
    df$significant <- lapply(df$p_values, function(x) x <=
                               user_threshold)
    col_names <- c("term_id", "p_values", "significant",
                   "term_size", "query_sizes", "intersection_sizes",
                   "source", "term_name", "effective_domain_size", "source_order")
  }
  else {
    col_names <- c("query", "significant", "p_value", "term_size",
                   "query_size", "intersection_size", "precision", "recall",
                   "term_id", "source", "term_name", "effective_domain_size",
                   "source_order")
    if (evcodes) {
      col_names <- append(col_names, c("evidence_codes",
                                       "intersection"))
      df$intersection <- mapply(function(evcodes, query) paste0(meta$genes_metadata$query[[query]]$ensgs[which(lengths(evcodes) >
                                                                                                                 0)], collapse = ","), df$intersections, df$query,
                                SIMPLIFY = TRUE)
      df$evidence_codes <- sapply(df$intersections, function(x) paste0(sapply(x[which(lengths(x) >
                                                                                        0)], paste0, collapse = " "), collapse = ","),
                                  USE.NAMES = FALSE)
    }
    df <- df[with(df, order(query, p_value)), ]
  }
  colnames(df)[colnames(df) == "native"] <- "term_id"
  colnames(df)[colnames(df) == "name"] <- "term_name"
  row.names(df) <- NULL
  df <- df[, col_names]
  return(list(result = df, meta = meta))
}
