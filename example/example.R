### INIT
rm(list = ls(all = TRUE))
options(warn = -1)
library(JubRi)



### BUILD QUERY-TERM LIST
query_dir <- file.path("example", "test_queries")
files <- list.files(query_dir, pattern = "*.txt", full.names = TRUE, recursive = TRUE)
queries_norm <- aggregrate_queries(files = files, normalize = TRUE)
queries_raw  <- aggregrate_queries(files = files, normalize = FALSE)
lapply(queries_norm, head)
lapply(queries_raw, head)



### LOAD MOUSE QUERY LIST OBJ
data(queries)
length(queries)



### DEFINE SEARCH TERMS
prog <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/2019.02.08.Drop10.All.Program.List.Cleanedx4.rds")
genes <- sort(as.vector(subset(prog, Set == "CDH1_P_8")$Genes))



### JUBRI
results <- JubRi(x = genes,
                 db = queries,
                 to_test_critera = 0.25,
                 background_size = 1000,
                 p_adj_method = "fdr",
                 pval_threshold = 0.05,
                 log2fc_threshold = 1,
                 seed = 42,
                 verbose = TRUE)

# results$FOREGROUND[1:5,1:5]
# results$BACKGROUND[1:5,1:5]
# head(results$META)

