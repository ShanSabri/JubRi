### INIT
rm(list = ls(all = TRUE))
options(warn = -1)
library(JubRi)



### BUILD QUERY-TERM LIST, OR 'DATABASE'
query_dir <- file.path("example", "test_queries")
files <- list.files(query_dir, pattern = "*.txt", full.names = TRUE, recursive = TRUE)
queries_norm <- aggregrate_queries(files = files, normalize = TRUE)
queries_raw  <- aggregrate_queries(files = files, normalize = FALSE)
# lapply(queries_norm, head)
# lapply(queries_raw, head)



### LOAD MOUSE QUERY LIST OBJ
data(queries)
length(queries) # 18039 entires



### DEFINE SEARCH TERMS
prog <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/2019.02.08.Drop10.All.Program.List.Cleanedx4.rds")
genes <- sort(as.vector(subset(prog, Set == "CDH1_P_7")$Genes)) # pluripotency-related genes



### JUBRI
results <- JubRi(x = genes,
                 db = queries,
                 to_test_critera = 0.50,
                 background_size = 1000,
                 p_adj_method = "fdr",
                 pval_threshold = 0.05,
                 log2fc_threshold = 1,
                 seed = 42,
                 verbose = TRUE)
# results$FOREGROUND[1:5,1:5]
# results$BACKGROUND[1:5,1:5]
# head(results$META)



### VISUALIZE OUTPUT
mp <- metaplot(results, top_n = 30, file = file.path("man", "figures", "metaplot.png"))



### ONE STOP SHOP
JubRi(x = genes,
      db = queries,
      to_test_critera = 0.50,
      background_size = 1000,
      p_adj_method = "fdr",
      pval_threshold = 0.05,
      log2fc_threshold = 1,
      seed = 42,
      verbose = TRUE) %>%
  metaplot()



# ### JUBRI ON PROGRAMS
# prog_ls <- split(prog, f = prog$Set)
# prog_ls <- prog_ls[grep("TC1|CDH1", names(prog_ls))]
# bg_size <- sapply(sort(sapply(prog_ls, nrow)), function(s){ # bg of 10x query set, max of 1000
#   min(1000, 10*s)
# })
# STDOUT <- lapply(seq_along(prog_ls), function(x){
#   prog_id <- names(prog_ls)[x]
#   p <- as.vector(prog_ls[[x]]$Genes)
#   bg <- bg_size[[prog_id]]
#   out <- JubRi(x = p,
#                db = queries,
#                to_test_critera = 0.25,
#                background_size = bg,
#                p_adj_method = "fdr",
#                pval_threshold = 0.05,
#                log2fc_threshold = 1,
#                seed = 42,
#                verbose = TRUE)
#   saveRDS(out, compress = TRUE, file = sprintf("~/Dropbox/PlathLab/Code/Jubri/v3/program_objs/%s.rds", prog_id))
# })
