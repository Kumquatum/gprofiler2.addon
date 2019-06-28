#' Format a gmt file into a list
#'
#' Take a path to a gmt file and format it into a easily readble/printing format
#'
#' @param path_to_gmt A path to the gmt file.
#' @return A list where name is the term of the gmt and content is the ids
#' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' @examples
#' path_to_gmt <- system.file("extdata", "h.all.v6.2.entrez.gmt",
#'                            package = "gprofiler2.addon", mustWork = TRUE)
#' formated_gmt <- format_gmt(path_to_gmt)
#'
#' @export

format_gmt <- function(path_to_gmt){
  if (!is.character(path_to_gmt)) {stop("custom_gmt should be a string representing GMT file's path")}
  if (!file.exists(path_to_gmt)) {stop("GMT file path must exist")}
  if (!grepl("\\.gmt$", path_to_gmt, ignore.case = TRUE)) {stop("File extension isn't GMT")}

  gmt_raw <- readChar(path_to_gmt, file.info(path_to_gmt)$size)
  gmt <- unlist(strsplit(gmt_raw, "\\n"))
  gmt <- strsplit(gmt, "\\t")
  names(gmt) <- sapply(gmt, `[[`, 1)
  gmt <- lapply(gmt, function(item) item[-c(1,2)])
}
