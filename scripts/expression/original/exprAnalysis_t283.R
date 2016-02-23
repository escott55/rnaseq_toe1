#!/usr/bin/Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicFeatures")
#biocLite("BiocStyle")
#biocLite("knitr")
#biocLite("rmarkdown")
#biocLite("Rsamtools")
#biocLite("GenomicFeatures")
#biocLite("GenomicAlignments")
#biocLite("BiocParallel")
#biocLite("DESeq2")
#biocLite("pheatmap")
#biocLite("RColorBrewer")
#biocLite("PoiClaClu")
#biocLite("ggplot2")

#samtools view -H file.bam > header.sam  # extract header only
#samtools reheader header.sam npc1603U3.sort.bam > file.new.bam
#samtools view -H file.bam > header.sam
#samtools view npc1603U3.sort.bam | cat header.sam - | sed s/chr// | \
    #samtools view -Sb - > npc1603U3.b37.bam
#samtools view npc1603A1.sort.bam | cat header.sam - | sed s/chr// | \
    #samtools view -Sb - > npc1603A1.b37.bam
#samtools view npc1603A3.sort.bam | cat header.sam - | sed s/chr// | \
    #samtools view -Sb - > npc1603A3.b37.bam
#samtools view /home/escott/workspace/toe1/clean/npc1603U4.sort.bam | \
    #cat header.sam - | sed s/chr// | samtools view -Sb - > npc1603U4.b37.bam


## ----style, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----
library("BiocStyle")
library("knitr")
library("rmarkdown")
library("Rsamtools")
library("GenomicFeatures")
library("ggplot2")
## ------------------------------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")
## ------------------------------------------------------------------------
library("DESeq2")

print("Running t283")

## ------------------------------------------------------------------------
#options(width=100)
#opts_chunk$set(message = FALSE, error = FALSE, warning = FALSE, fig.width=5, fig.height=5)

#dir <- system.file("extdata", package="airway", mustWork=TRUE)
maindir <- "/home/escott/workspace/toe1"
resdir <- file.path(maindir,"results/t283")

stopifnot( file.exists(resdir) )

## ------------------------------------------------------------------------
csvfile <- file.path(maindir,"raw/t283/SampleTable.txt")

sampleTable <- read.csv(csvfile,row.names=1,sep="\t",comment.char="#")

sampleTable$id <- row.names(sampleTable)
##sampleTable <- sampleTable[ sampleTable$Batch == 2, ]

filenames <- as.vector(sampleTable$File)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=1000000)

## ------------------------------------------------------------------------
seqinfo(bamfiles[1])

## ------------------------------------------------------------------------
#gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
#gtffile <- file.path(maindir,"resources/gencode.vM7.annotation_encode_ccds.gtf")
#gtffile <- file.path(maindir,"resources/Homo_sapiens.GRCh37.b37.gtf")
gtffile <- file.path(maindir,"resources/Homo_sapiens.GRCh37.b37.protcoding.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))

## ------------------------------------------------------------------------
(ebg <- exonsBy(txdb, by="gene"))


## ------------------------------------------------------------------------
register(SerialParam())

## ------------------------------------------------------------------------
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

## ------------------------------------------------------------------------
se
dim(se)
assayNames(se)
head(assay(se), 3)
colSums(assay(se))

## ------------------------------------------------------------------------
rowRanges(se)

## ------------------------------------------------------------------------
str(metadata(rowRanges(se)))

## ------------------------------------------------------------------------
colData(se)

## ------------------------------------------------------------------------
(colData(se) <- DataFrame(sampleTable))

## ------------------------------------------------------------------------
#se$cell
# Reorder levels to put affecteds first
se$Status <- relevel(se$Status, "A")

## ------------------------------------------------------------------------
round( colSums(assay(se)) / 1e6, 1 )

## ------------------------------------------------------------------------
colData(se)

## ------------------------------------------------------------------------
# Should modifie this equation to include family
#dds <- DESeqDataSet(se, design = ~ type)
#dds <- DESeqDataSet(se, design = ~ Status + Family + Family:Status)
dds <- DESeqDataSet(se, design = ~ Status)

## ------------------------------------------------------------------------
countdata <- assay(se)
head(countdata, 3)

## ------------------------------------------------------------------------
coldata <- colData(se)

## ------------------------------------------------------------------------
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ Status))
                                  #design = ~ Family + Status + Family:Status))

## ------------------------------------------------------------------------
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

## ----rlog----------------------------------------------------------------
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

objsavefile <- file.path(resdir,"deseq2_toe1.RData")
#save(dxd, dxr1, file = obj.save.file)
save.image(objsavefile) # Saves image
# load(objsavefile)

## ----rldplot, fig.width=7, fig.height=3.75-------------------------------

dds <- estimateSizeFactors(dds)
pdf(file.path(resdir,"qc2.pdf"))
par( mfrow = c( 1, 2 ) )
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
dev.off()

normcnts <- head(log2(counts(dds, normalized=TRUE)))

######################## VERIFY SAMPLE DISTANCES ##########################
source("sampleDistances.R")
cntDistHeatmap( rld, dds, resdir, typecol="Status", shapecol="Family" )
poisDistHeatmap( rld, dds, resdir, typecol="Status", shapecol="Family" )
samplePCAPlots( rld, resdir, typecol="Status", shapecol="Family" )

## ------------------------------------------------------------------------
dds <- DESeq(dds)

# Fix Ensembl IDs
row.names(dds) <- substr(row.names(dds),0,18)

## ------------------------------------------------------------------------
(res <- results(dds))

## ------------------------------------------------------------------------
mcols(res, use.names=TRUE)

## ------------------------------------------------------------------------
summary(res)

## ------------------------------------------------------------------------
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

## ------------------------------------------------------------------------
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

## ------------------------------------------------------------------------
results(dds, contrast=c("Status", "U", "A"))

## ------------------------------------------------------------------------
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

## ------------------------------------------------------------------------
sum(res$padj < 0.1, na.rm=TRUE)

## ------------------------------------------------------------------------
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

## ------------------------------------------------------------------------
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])


########################## START TOP GENE #################################

source("genePlots.R")
topGene <- rownames(res)[which.min(res$padj)]
#topGene <- "ENSMUSG00000037071.2"
pltCntsPdf( dds, topGene, resdir, typecol="Status", shapecol="Family" )

toe1Gene <- "ENSG00000132773"
pltCntsPdf( dds, toe1Gene, resdir, typecol="Status", shapecol="Family" )

#mboat7acc <- "ENSMUSG00000035596"
###trgtGene <- rownames(res)[which.isMatchingStartingAt(mboat7acc,rownames(res))]

########################## END TOP GENE ###################################

########################## Plot MA  ###################################
source("maPlots.R")
genMAPlots( res, topGene, resdir )

###################### END Plot MA  ###################################

## ----histpvalue2---------------------------------------------------------
#p values for genes with mean normalized count larger
pdf(file.path(resdir,"plothist.pdf"))
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
dev.off()

library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

## ----genescluster--------------------------------------------------------
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Status","Family")])
#df <- as.data.frame(colData(rld)[,c("Status","id")])
tfile <- file.path(resdir,"geneClustHeatmap.pdf")
#anyZero <- function(x) any(abs(x) == 0 )
#anyNull <- function(x) any(is.null(x))
#apply( mat, 1, anyZero )
pheatmap(mat, annotation_col=df, filename=tfile)
# Getting error: 'gpar' element 'fill' must not be length 0

## ----sensitivityovermean, fig.width=6------------------------------------
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
pdf(file.path(resdir,"meanNormCounts.pdf"))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
dev.off()


library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
#library("org.Mm.eg.db")
## ------------------------------------------------------------------------
#columns(org.Mm.eg.db)

## ------------------------------------------------------------------------
# remove the versioning number from the ensembl ids
#res$normids <- substr(row.names(res),0,18)
#row.names(res) <- normids

## ------------------------------------------------------------------------
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

## ------------------------------------------------------------------------
resOrdered <- res[order(res$padj),]
head(resOrdered)

## ----eval=FALSE----------------------------------------------------------
#ressig <- subset(res, padj < 0.1)
resOrderedDF <- as.data.frame(subset(resOrdered,padj<0.05)) #[1:100,]
write.csv(resOrderedDF, file=file.path(resdir,"sigset.csv"))

write.csv(as.data.frame(resOrdered), file=file.path(resdir,"fullset.csv"))

resOrderedDFClean <- as.data.frame(subset(resOrdered,complete.cases(symbol)))
toe1 <- resOrderedDFClean[ resOrderedDFClean$symbol == "TOE1", ]

#mboat7acc <- "ENSMUSG00000035596"
#topGene <- rownames(dds)[substr(rownames(dds),0,18) == mboat7acc]
#pltCntsPdf( dds, topGene, resdir )

## ----eval=FALSE----------------------------------------------------------
## library("ReportingTools")
## htmlRep <- HTMLReport(shortName="report", title="My report",
##                       reportDirectory="./report")
## publish(resOrderedDF, htmlRep)
## url <- finish(htmlRep)
## browseURL(url)

## ------------------------------------------------------------------------
(resGR <- results(dds, lfcThreshold=1, format="GRanges"))

## ------------------------------------------------------------------------
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

## ------------------------------------------------------------------------
library("Gviz")

## ------------------------------------------------------------------------
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

## ------------------------------------------------------------------------
sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))

## ----gvizplot------------------------------------------------------------
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
               type="h", name="log2 fold change", strand="+")
pdf(file.path(resdir,"groupTracks.pdf"))
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
dev.off()

q()
#################### REMOVING HIDDEN BATCH EFFECTS #######################
## ------------------------------------------------------------------------
library("sva")

## ------------------------------------------------------------------------
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
svseq$sv

## ----svaplot-------------------------------------------------------------
pdf(file.path(resdir,"svaPlot.pdf"))
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)
dev.off()

## ------------------------------------------------------------------------
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex

## ----svaDE, eval=FALSE---------------------------------------------------
## ddssva <- DESeq(ddssva)

#################### TIME COURSE EXPERIMENTS ##############################
## ------------------------------------------------------------------------
library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

## ----fissionDE-----------------------------------------------------------
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

## ----fissioncounts, fig.width=6, fig.height=4.5--------------------------
pltcnts <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("minute","strain"), returnData=TRUE)
pdf(file.path(resdir,"fissionCnts.pdf"))
ggplot2(pltcnts, aes(x=minute, y=count, color=strain, group=strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
dev.off()

## ------------------------------------------------------------------------
resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]

## ------------------------------------------------------------------------
betas <- coef(ddsTC)
colnames(betas)

## ----fissionheatmap------------------------------------------------------
#library("pheatmap")
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pdf(file.path(resdir,"fissionPHeatmap.pdf"))
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
dev.off()

## ------------------------------------------------------------------------
sessionInfo()

