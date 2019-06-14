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


#for(protein in list.files(path = "results")) {
    project_dir <- file.path("results", protein)
    
 #   if(!file.exists(file.path(project_dir, "hits.txt"))) {
        cat(paste0("processing ", protein, "\n"))

        WINDOWCOUNTS <- read_tsv(file.path(project_dir, "windows_countmatrix.txt.gz")) %>% column_to_rownames("ID")
        
        SAMPLEINFO <- read_tsv(file.path(project_dir, "samples.txt")) %>% as.data.frame %>% column_to_rownames("fileID")
        
        SAMPLEINFO <- SAMPLEINFO[ colnames(WINDOWCOUNTS), ]
        #ANNOTATION <- read_tsv(file.path(project_dir, "windows_annotation.txt.gz"))
        
        #ANNOTATION <- fread(input = "")
        
        SAMPLEINFO <- SAMPLEINFO %>% mutate(type = factor(type, levels = c("IP", "SMI")))
        #SAMPLEINFO <- SAMPLEINFO %>% rename_at(vars(oldnames), ~ newnames) 
        
        dds <- DESeqDataSetFromMatrix(countData = WINDOWCOUNTS,
                                     colData = SAMPLEINFO,
                                     design = ~ type)                                  
        
        
        selection <- rowSums( counts(dds) > 0 ) >= 2
        
        
        dds <- dds[ selection, ]
        
        #ANNOTATION <- ANNOTATION[ selection, ]
        
        dds <- dds %>% DESeq(betaPrior = T)
        
        res <- results_DEWSeq(dds, contrast = c("type", "IP", "SMI"), tidy = T)
        resAnn <- merge(res,ANNOTATION_WINDOWS,by='unique_id')
        resWindows <- mergeWindows(annRes=resAnn,minDist=0,padjWindow='bonferroni',ncores=5)
        resIHW <- ihw(pBonferroni ~ baseMean, data =  resWindows, alpha=0.05, nfolds=10)
        resIHW <- as.data.frame(resIHW)
        colnames(resIHW) <- paste('IHW',colnames(resIHW),sep='.')
        res <- cbind(resWindows,resIHW[,c('IHW.adj_pvalue','IHW.weight','IHW.weighted_pvalue')]) %>%  
          mutate(significant = !is.na(IHW.adj_pvalue) & IHW.adj_pvalue < 0.05 & log2FoldChange > 1)
        
        SIG <- res %>% dplyr::filter(significant) %>% arrange(desc(log2FoldChange)) %>%
            dplyr::left_join( y = (counts(dds) %>% as.data.frame %>% rownames_to_column(var = "unique_id")), by = "unique_id")
        
        write_tsv(SIG %>% dplyr::select(gene_name, gene_type) %>% distinct() %>% .["gene_type"] %>% table %>% as.data.frame,
                        path = file.path(project_dir, "_tmp_hits_type_table.txt")) 
            
        #res <- res %>% left_join(y = ANNOTATION, by = "featureID")
                                     
        write_tsv(SIG, path = file.path(project_dir, "_tmp_hits.txt"), col_names = F)
        
        write_tsv(SIG %>% dplyr::select(gene_id, gene_name) %>% distinct(),
                  path = file.path(project_dir, "_tmp_hits_id_name.txt"), col_names = F) 
        
        write_tsv(SIG %>% .["gene_id"] %>% table %>% as.data.frame,
                  path = file.path(project_dir, "_tmp_hits_id_table.txt")) 
        
        
        save(list = c("dds", "res", "project_dir", "protein", "SIG"), file = file.path(project_dir, "_tmp_Workspace.Rdata"))
#    }
#}
#sessionInfo()

#file.path(R.home("bin"), "R")
