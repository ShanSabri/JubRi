#' Create metaplot of results
#'
#' @param x output metatable from \code{JubRi()}
#' @param top_n number of top enriched terms to show on on metaplot panel
#' @param to_label significance levels to label terms
#' @param file output file
#'
#' @return metaplot
#' @importFrom cowplot plot_grid
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils head
#' @importFrom stats complete.cases
#'
#' @export
metaplot <- function(x, top_n = 30,
                     to_label = c("PVAL_LOG2FC", "LOG2FC"),
                     file = file.path(getwd(), "metaplot.png")) {
  cols <- c(
    "PVAL_LOG2FC" = "firebrick",
    "LOG2FC" = "pink3",
    "PVAL" = "pink1",
    "NOT_SIGNIFICANT" = "grey"
  )
  sign_tbl <- table(x$META$SIGNIFICANT)

  # VOLCANO PLOT
  if (sign_tbl[1] + sign_tbl[2] >= 200) {
    repel_data <- subset(x$META, SIGNIFICANT %in% to_label)
    repel_data <- utils::head(repel_data[with(repel_data, order(SIGNIFICANT, PVAL)), ], n = top_n)
    repel <- ggrepel::geom_text_repel(data = repel_data, aes(label = TERM), size = 2)
  } else {
    repel <- ggrepel::geom_text_repel(data = subset(x$META, SIGNIFICANT %in% to_label), aes(label = TERM), size = 2)
  }
  VOLCANO <- ggplot(x$META[stats::complete.cases(x$META), ], aes(x = LOG2FC, y = -log10(PVAL))) +
    geom_point(aes(colour = SIGNIFICANT)) +
    repel +
    scale_colour_manual(values = cols) +
    labs(
      colour = "", x = "Log2(FC + 0.01)", y = "-Log10(P-value)",
      title = sprintf(
        "%s Terms satisfy Fold Change & P-value critera\n%s Terms satisfy Fold Change critera",
        sign_tbl[1], sign_tbl[2]
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )

  # ENRICHMENT
  enrichment <- subset(x$META, SIGNIFICANT %in% to_label)
  enrichment <- utils::head(enrichment[with(enrichment, order(SIGNIFICANT, -LOG2FC)), ], n = top_n)
  enrichment <- enrichment[with(enrichment, order(SIGNIFICANT, -FG_MU)), ]
  enrichment$TERM <- factor(enrichment$TERM, levels = enrichment$TERM)
  ENRICHMENT <- ggplot(enrichment, aes(x = TERM, y = FG_MU)) +
    geom_bar(aes(fill = SIGNIFICANT), stat = "identity", colour = "black", alpha = 0.8) +
    geom_bar(data = enrichment, aes(x = TERM, y = BG_MU), stat = "identity", fill = "grey", colour = "grey55", alpha = 0.5) +
    scale_fill_manual(values = cols) +
    labs(
      x = "", y = "Average\nNormalized Frequency",
      title = sprintf("Top %s most enriched words associated with query", top_n),
      subtitle = "Sorted by significance level & frequency"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  # WORDCLOUD
  wc <- subset(x$META, SIGNIFICANT %in% to_label)
  if (sum(sign_tbl[to_label]) >= 100) {
    wc <- utils::head(wc[with(wc, order(SIGNIFICANT, PVAL)), ], n = 100)
  }
  wc$FC_BIN <- as.factor(cut(wc$LOG2FC, labels = FALSE, breaks = 10))
  WORDCLOUD <- ggplot(wc) +
    aes(x = 1, y = 1, size = LOG2FC, label = TERM, colour = LOG2FC) +
    ggrepel::geom_text_repel(segment.size = 0, force = 10, max.iter = 5000) +
    scale_size(range = c(1, 5), guide = FALSE) +
    scale_color_gradientn(colours = c("grey", "red", "firebrick")) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    labs(x = "", y = "", title = "Wordcloud of significant terms by enrichment") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  # SCATTER
  scatter <- x$META[stats::complete.cases(x$META), ]
  scatter$BIN <- as.factor(cut(scatter$FG_MU, labels = FALSE, breaks = 4))
  SCATTER <- ggplot(scatter, aes(FG_MU, BG_MU, colour = BIN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed", alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw(base_size = 14) +
    labs(
      x = "Foreground", y = "Background",
      title = "Top de-enriched terms associated with query",
      subtitle = "Stratified by foreground frequency bin"
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    )

  # DE-ENRICHMENT
  de_enrichment <- subset(x$META, LOG2FC <= 0)
  de_enrichment$BIN <- as.factor(cut(de_enrichment$FG_MU, labels = FALSE, breaks = 4))
  de_enrichment %>%
    group_by(BIN) %>%
    top_n(10, desc(LOG2FC)) %>%
    arrange(desc(FG_MU)) -> top_de_enrichment
  top_de_enrichment$TERM <- factor(top_de_enrichment$TERM, levels = top_de_enrichment$TERM)
  DE_ENRICHMENT <- ggplot(top_de_enrichment, aes(x = TERM, y = FG_MU, fill = BIN)) +
    facet_wrap(~BIN, nrow = 1, scales = "free") +
    geom_bar(data = top_de_enrichment, aes(x = TERM, y = BG_MU), stat = "identity", fill = "grey", colour = "grey55", alpha = 0.5) +
    geom_bar(aes(fill = BIN), stat = "identity", colour = "black") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      x = "", y = "Average Normalized Frequency"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  METAPLOT <- cowplot::plot_grid(VOLCANO, ncol = 1, labels = "A") +
    cowplot::plot_grid(ENRICHMENT, WORDCLOUD, ncol = 1, labels = c("B", "C")) +
    cowplot::plot_grid(SCATTER, DE_ENRICHMENT, ncol = 1, labels = c("D", ""), rel_heights = c(1, 2)) +
    ggsave(file, heigh = 7, width = 25)

  return(METAPLOT)
}
