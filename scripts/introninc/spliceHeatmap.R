#!/usr/bin/Rscript

#source("https://bioconductor.org/biocLite.R")
library("gplots")

data(mtcars)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
  
heatmap.2(x, dendrogram="row")


maindir <- "/home/escott/workspace/toe1"
resfile <- file.path(maindir,"splicefigs/minusbatch1/refseq_npcnew_final.tsv")

m <- data.frame(read.table(resfile, header=T, sep="\t"))
row.names(m) <- m$test_id
msub <- m[,10:19]
m_matrix <- data.matrix(msub)
#cor_t <- 1 - cor(m_matrix)
cor_t <- 1 - abs(cor(t(m_matrix))) # edited
#cor_t <- 1 - cor(t(m_matrix))^2

distancem <- dist(t(m_matrix))
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(cor_t, Colv=NA, scale="column")

#Rowv=dendcompletem, 

distancem <- dist(m_matrix)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(m_matrix, Rowv=dendcompletem, Colv=NA, scale="column")

#heatmap.2(x)

maindir <- "/home/escott/workspace/toe1"
resfile <- file.path(maindir,"splicefigs/minusbatch1/refseq_npcnew_final.tsv")

m <- data.frame(read.table(resfile, header=T, sep="\t"))
msub <- m[,10:19]
m_matrix <- data.matrix(msub)
cor_t <- 1 - cor(t(m_matrix))^2

distancem <- dist(t(m_matrix))
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(cor_t, Colv=NA, scale="column")


