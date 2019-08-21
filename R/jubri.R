#' JubRi: a method for literature-based keyword ontology
#'
#' @param x a character vector of genes IDs thought to be found in database
#' @param db a database of literature queries
#' @param to_test_critera quantile to cut the number of terms to test. 0 if wanting to test all.
#' @param background_size number of elements to randomly pull from the database to use as a background model
#' @param p_adj_method p-value correction method
#' @param pval_threshold p-value significance threshold
#' @param log2fc_threshold log2 fold change significance threshold
#' @param seed set the initial seed for background aggregration and p-value computation
#' @param verbose speak alot
#'
#' @return a list of data.frames, each element corresponding to the foreground, background and metadata data.frames
#' @importFrom stats p.adjust quantile
#'
#' @export
#' @examples
#' genes <- c("NANOG", "DPPA5A", "POU5F1", "SOX2")
#' data(queries)
#' results <- JubRi(x = genes,
#'                  db = queries,
#'                  to_test_critera = 0.2,
#'                  background_size = 1000,
#'                  p_adj_method = "fdr",
#'                  pval_threshold = 0.05,
#'                  log2fc_threshold = 1,
#'                  seed = 42,
#'                  verbose = TRUE)
JubRi <- function(x, db, to_test_critera = 0.25, background_size = 1000,
                  p_adj_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                  pval_threshold = 0.01, log2fc_threshold = 1,
                  seed = 42, verbose = TRUE) {

  fg_terms <- toupper(intersect(tolower(x), tolower(names(db))))
  if (verbose) {
    log(sprintf("Foreground length of %s", length(x)))
    log(sprintf("Database length of %s", length(db)))
    log(sprintf("Database contains %s terms overlapping with foreground", length(fg_terms)))
  }
  stopifnot(length(fg_terms) > 0)
  fg <- db[fg_terms]

  # RESHAPE
  if (verbose) {
    log("Reshaping data")
  }
  fg <- reshape_to_df(fg)

  # DETERMINE TEST TERMS
  if (verbose) {
    log("Determining number of terms to test")
  }
  meta <- data.frame(
    TERM = row.names(fg),
    FG_MU = rowMeans(fg),
    BG_MU = NA,
    TO_TEST = NA,
    ADJ_PVAL = NA,
    LOG2FC = NA
  )
  meta$TO_TEST <- ifelse(meta$FG_MU >= stats::quantile(meta$FG_MU, to_test_critera), TRUE, FALSE)
  n_test <- nrow(meta[meta$TO_TEST == TRUE, ])
  if (verbose) {
    log(sprintf("Testing %s/%s terms", n_test, nrow(meta)))
  }

  # BACKGROUND
  if (verbose) {
    log(sprintf("Making background of length %s", background_size))
  }
  stopifnot(background_size < length(db))
  set.seed(seed)
  bg <- sample(names(db), 1000)
  bg <- db[bg]
  bg <- reshape_to_df(bg)

  # TESTING
  shared_terms <- intersect(
    as.vector(meta$TERM[meta$TO_TEST == TRUE]),
    row.names(bg)
  )
  log(sprintf("WARNING - %s terms not found in background. Consider increasing background set size", n_test - length(shared_terms)))
  bg_mu <- rowMeans(bg[row.names(bg) %in% shared_terms, ])

  log(sprintf("Computing adjusted P-values and Log2(FC) for %s terms", length(shared_terms)))
  p <- wilcox(fg = fg, bg = bg, terms = shared_terms) #  parallel = parallel
  meta$BG_MU[match(names(bg_mu), row.names(fg))] <- bg_mu
  meta$ADJ_PVAL[match(names(p), row.names(fg))] <- stats::p.adjust(p, method = p_adj_method, n = length(p))
  meta$LOG2FC <- ifelse(meta$TO_TEST == TRUE, log2((meta$FG_MU + 0.01) / (meta$BG_MU + 0.01)), NA)
  meta$SIGNIFICANT <- ifelse(meta$ADJ_PVAL <= pval_threshold & abs(meta$LOG2FC) >= log2fc_threshold, "PVAL_LOG2FC",
    ifelse(abs(meta$LOG2FC) >= log2fc_threshold, "LOG2FC",
      ifelse(meta$ADJ_PVAL <= pval_threshold, "PVAL", "NOT_SIGNIFICANT")
    )
  )
  meta$SIGNIFICANT <- factor(meta$SIGNIFICANT, levels = c("PVAL_LOG2FC", "LOG2FC", "PVAL", "NOT_SIGNIFICANT"))

  log("DONE")

  return(list(
    FOREGROUND = fg,
    BACKGROUND = bg,
    META = meta[with(meta, order(SIGNIFICANT, -LOG2FC, ADJ_PVAL)), ])
    )
}
