#' Combine query files from PubMedScrapeR into a single list object
#'
#' @param files a vector containing query file names that are to be combined into a single object
#' @param normalize boolean to indicate if query term counts should be normalized by number of abstracts
#' @param blacklist a vector containing terms to exclude, typically very common words
#'
#' @return a list object where each element corresponds to the term count of a given query
#' @importFrom pbapply pblapply
#' @importFrom dplyr data_frame
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table

#' @export
#' @examples
#' library(JubRi)
#' query_dir <- file.path("example", "test_queries")
#' files <- list.files(query_dir, pattern = "*.txt", full.names = TRUE, recursive = TRUE)
#' queries <- aggregrate_queries(files = files, normalize = TRUE)
#' # lapply(queries, head)
aggregrate_queries <- function(files, normalize = TRUE, blacklist = c()) {
  exclude <- c(
    "nlmcategory", "expression",
    "cells", "cell",
    "genes", "gene"
  )
  blacklist <- c(blacklist, exclude)

  obj <- pbapply::pblapply(files, function(x) {
    r <- utils::read.table(x, header = TRUE, row.names = 1, sep = "\t")
    if (nrow(r) == 0) {
      return(NA)
    }
    id <- sapply(strsplit(names(r), "_", fixed = TRUE), function(x) (x[1]))
    num_abstracts <- as.numeric(sapply(strsplit(names(r), "_", fixed = TRUE), function(x) (x[2])))
    r <- r[nzchar(gsub("[0-9]+", "", row.names(r))), , drop = FALSE]
    r <- r[!row.names(r) %in% blacklist, , drop = FALSE]

    if (normalize == TRUE) {
      r <- dplyr::data_frame(ID = id, TERM = row.names(r), COUNT = r[[1]] / num_abstracts)
    } else {
      r <- dplyr::data_frame(ID = id, TERM = row.names(r), COUNT = r[[1]])
    }

    return(r)
  })

  names(obj) <- tools::file_path_sans_ext(basename(files))
  obj <- obj[!is.na(obj)]

  return(obj)
}



#' Reshape a list of query entries to a dataframe
#'
#' @param x a list of query entries
#'
#' @return an aggregrated data.frame with IDs as columns and TERMS as rows
#' @importFrom reshape2 dcast
reshape_to_df <- function(x) {
  x <- do.call(rbind.data.frame, x)
  x <- unique(x)
  x <- reshape2::dcast(TERM ~ ID, value.var = "COUNT", data = x)
  x[is.na(x)] <- 0
  row.names(x) <- x$TERM
  x$TERM <- NULL
  return(x)
}




#' Computer P-value enrichment of query terms
#'
#' @param fg foreground matrix
#' @param bg background matrix
#' @param terms a vector of character terms that are to be tested
#'
#' @return p-values
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
wilcox <- function(fg, bg, terms) {
  # if (parallel == TRUE) {
  #   cl <- parallel::makeCluster(parallel::detectCores() - 1)
  #   parallel::clusterExport(cl, varlist = c("fg", "bg"), envir = environment())
  #   p <- pbapply::pbsapply(terms, function(y) {
  #     input_x <- as.numeric(fg[y, , drop = FALSE])
  #     input_y <- as.numeric(bg[y, , drop = FALSE])
  #     return(stats::wilcox.test(input_x, input_y, alternative = "greater")$p.value)
  #   }, USE.NAMES = TRUE, cl = cl)
  #   parallel::stopCluster(cl)
  # } else {
    p <- pbapply::pbsapply(terms, function(y) {
      input_x <- as.numeric(fg[y, , drop = FALSE])
      input_y <- as.numeric(bg[y, , drop = FALSE])
      return(stats::wilcox.test(input_x, input_y, alternative = "greater")$p.value)
    }, USE.NAMES = TRUE)
  # }

  return(p)
}



#' Verbose logging function
#'
#' @param msg character log statement
#'
#' @return a character log statement tagged with current system time
log <- function(msg) {
  message(sprintf("%s: %s", Sys.time(), msg))
}
