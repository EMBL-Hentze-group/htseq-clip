library(pheatmap)

args = commandArgs(trailingOnly = T)
  
herv = read.table(args[1],stringsAsFactors = FALSE,header = T)

data = data.matrix(herv[,2:ncol(herv)])
#data = scale(data)
colnames(data) = colnames(herv)[2:4]
rnames = herv[,1]
rownames(data) = rnames

pdf(NULL)
pheatmap(data,fontsize_row = 7,width = strtoi(args[3]), height = strtoi(args[4]), show_rownames=T,show_colnames=T, scale= 'row', filename = args[2],cluster_rows=F,cluster_cols=F)

dev.off()