#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied.", call.=FALSE)
} 

# read in annotation file
load("../annotation/gencode27ns_flattenedwindows_width50_steps20.Rda")

suppressPackageStartupMessages({
      require(tidyverse)
      require(DESeq2)
      require(DEWSeq)
      require(cowplot)
      require(IHW)
      require(tidyr)
      #require(data.table)
})

#protein <- "SSB-HepG2"
protein <- args[1]

project_dir <- file.path("results", protein)
tmp_dir <- file.path('/tmp',protein)
fig_dir <- file.path('/tmp',protein,'figures')
dir.create(tmp_dir)
dir.create(fig_dir)
cat(paste0("processing ", protein, "\n"))
# sanity check for files
proj_files <- list.files(file.path("results",protein))
if(!"windows_countmatrix.txt.gz" %in% proj_files){
  stop("Cannot find count data file windows_countmatrix.txt.gz in project folder")
}else if(!"samples.txt" %in% proj_files){
  stop("Cannot find sample information file samples.txt in project folder")
}
WINDOWCOUNTS <- read_tsv(file.path(project_dir, "windows_countmatrix.txt.gz")) %>% column_to_rownames("ID")

SAMPLEINFO <- read_tsv(file.path(project_dir, "samples.txt")) %>% as.data.frame %>% column_to_rownames("fileID")

SAMPLEINFO <- SAMPLEINFO[ colnames(WINDOWCOUNTS), ]

SAMPLEINFO <- SAMPLEINFO %>% mutate(type = factor(type, levels = c("IP", "SMI")))

dds <- DESeqDataSetFromMatrix(countData = WINDOWCOUNTS,colData = SAMPLEINFO,design = ~ type)                                  

selection <- rowSums( counts(dds) > 0 ) >= 2

dds <- dds[ selection, ]

message("Running DESeq2")
dds <- dds %>% DESeq(betaPrior = T)
message("Running DEWSeq")
res <- results_DEWSeq(dds, contrast = c("type", "IP", "SMI"), tidy = T)
resAnn <- merge(res,ANNOTATION_WINDOWS,by='unique_id')
message("Merging windows")
resWindows <- mergeWindows(annRes=resAnn,minDist=0,padjWindow='bonferroni',ncores=5)
message("Running IHW")
resIHW <- ihw(pBonferroni ~ baseMean, data =  resWindows, alpha=0.05, nfolds=10)
# save IHW figures
# cross validation
pdf(file.path(fig_dir,'IHW_cross_validation.pdf'),width=10,height=9,title='IHW cross validation')
plot(resIHW)
dev.off()
# decision boundary
pdf(file.path(fig_dir,'IHW_decision_boundary.pdf'),width=10,height=9,title='IHW decision boundary')
plot(resIHW,what='decisionboundary')
dev.off()
resIHW <- as.data.frame(resIHW)
colnames(resIHW) <- paste('IHW',colnames(resIHW),sep='.')
res <- cbind(resWindows,resIHW[,c('IHW.adj_pvalue','IHW.weight','IHW.weighted_pvalue')]) %>%  
  mutate(significant = !is.na(IHW.adj_pvalue) & IHW.adj_pvalue < 0.05 & log2FoldChange > 1)

SIG <- res %>% dplyr::filter(significant) %>% arrange(desc(log2FoldChange)) %>%
    dplyr::left_join( y = (counts(dds) %>% as.data.frame %>% rownames_to_column(var = "unique_id")), by = "unique_id")

write_tsv(SIG %>% dplyr::select(gene_name, gene_type) %>% distinct() %>% .["gene_type"] %>% table %>% as.data.frame,
                path = file.path(tmp_dir, "_tmp_hits_type_table.txt")) 
                                  
write_tsv(SIG, path = file.path(tmp_dir, "_tmp_hits.txt"), col_names = F)

write_tsv(SIG %>% dplyr::select(gene_id, gene_name) %>% distinct(),
          path = file.path(tmp_dir, "_tmp_hits_id_name.txt"), col_names = F) 

write_tsv(SIG %>% .["gene_id"] %>% table %>% as.data.frame,
          path = file.path(tmp_dir, "_tmp_hits_id_table.txt")) 

save(list = c("dds", "res", "project_dir", "protein", "SIG"), file = file.path(project_dir, "_tmp_Workspace.Rdata"))
