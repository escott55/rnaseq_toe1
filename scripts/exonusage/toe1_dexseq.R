#!/usr/bin/Rscript

# DEPRICATED
## ----knitr, echo=FALSE, results='hide'-----------------------------------
#library("knitr")
#opts_chunk$set(tidy=FALSE,dev="pdf",fig.show="hide",
               #fig.width=4,fig.height=4.5,
               #message=FALSE)
## ----style, eval=TRUE, echo=FALSE, results='asis'------------------------
#BiocStyle::latex()
## ----options,results='hide',echo=FALSE-----------------------------------
#options(digits=2, width=80, prompt=" ", continue=" ")
## ----systemFile----------------------------------------------------------
#pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
#list.files(pythonScriptsDir)
## ----systemFileCheck,echo=FALSE,results='hide'---------------------------
#system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

# System calls:
# python /usr/local/lib/R/site-library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
#        resources/gencode.vM7.annotation.gtf resources/gencode.vM7.annotation.gff
#
# python /usr/local/lib/R/site-library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
#       resources/gencode.vM8.annotation.gtf resources/gencode.vM8.annotation.gff
#
# python /usr/local/lib/R/site-library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
#        resources/mouse_grcm38_refseq_sub.gtf resources/mouse_grcm38_refseq.gff
#
# Generate counts


library(ggplot2)
library(reshape)


ng1 = theme(panel.background = element_rect(fill = "white",colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    axis.line = element_line(size = 1.2, colour="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black",size=15),
    axis.text.x = element_text(colour="black",angle=45,hjust=1,vjust=1),
    axis.title = element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    strip.text.y = element_text(color="black",face="bold",size=15),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white")
    #legend.title = element_blank()
    )

basedir <- "/media/data/workspace/toe1"
cntDir <- file.path(basedir,"raw/newnpc/exoncnt")
flattenedFile <- file.path(basedir,"resources/Homo_sapiens.GRCh37.b37.exons.gff")
figDir <- file.path(basedir, "exonusage/newnpc")

csvfile <- file.path(basedir,"raw/newnpc/SampleTable.txt")

print(csvfile)
sampleTable <- read.csv(csvfile,row.names=1,sep="\t",comment.char="#")
sampleTable$Family <- paste("fam",sampleTable$Family,sep="")

sampleTable$id <- row.names(sampleTable)
#sampleTable <- sampleTable[ sampleTable$Family != "fam1101", ]
sampleTable <- sampleTable[ sampleTable$Batch != 1, ]

countFiles <- file.path(cntDir,paste(sampleTable$samp,"cnts.txt",sep="_"))

## ------------------------------------------------------------------------
#list.files(dir)
## ------------------------------------------------------------------------
#csvfile <- file.path(dir,"sample_table.csv")
#(sampleTable <- read.csv(csvfile,row.names=1,sep="\t"))
#filenames <- as.vector(sampleTable$File)
#file.exists(filenames)
## ----loadDEXSeq----------------------------------------------------------
#inDir = system.file("extdata", package="pasilla")
#countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
#countFiles
#flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
#flattenedFile
#countFiles = list.files(inDir, pattern="_cnts.txt$", full.names=TRUE)

## ----sampleTable---------------------------------------------------------
#sampleTable = data.frame(
   #row.names = c( "KO1", "KO17", "KO25", "KO29", "WT9", "WT10", "WT15", "WT16" ),
   #condition = c("knockout", "knockout","knockout", "knockout", "control",
                 #"control",  "control", "control" ),
   #libType = c( "single-end", "single-end", "single-end", "single-end",
               #"single-end", "single-end", "single-end", "single-end" ) )

## ----displaySampleTable--------------------------------------------------
print(sampleTable)

## ----makeecs, eval=TRUE--------------------------------------------------
suppressPackageStartupMessages( library( "DEXSeq" ) )

dxd = DEXSeqDataSetFromHTSeq(
   countFiles,
   sampleData=sampleTable,
   design= ~ sample + exon + Status:exon,
   flattenedfile=flattenedFile )

## ----start---------------------------------------------------------------
# Trim dataset to expedite this process
#genesForSubset = read.table( 
  #file.path(inDir, "geneIDsinsubset.txt"), 
  #stringsAsFactors=FALSE)[[1]]
#dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

## ----seeColData----------------------------------------------------------
colData(dxd)

## ----seeCounts-----------------------------------------------------------
head( counts(dxd), 5 )

## ----seeSplitted---------------------------------------------------------
split( seq_len(ncol(dxd)), colData(dxd)$exon )

## ----seeCounts2----------------------------------------------------------
head( featureCounts(dxd), 5 )

## ----fData---------------------------------------------------------------
head( rowRanges(dxd), 3 )
## ----pData---------------------------------------------------------------
sampleAnnotation( dxd )

## ----sizeFactors1--------------------------------------------------------
dxd = estimateSizeFactors( dxd )

## ----estDisp1------------------------------------------------------------
dxd = estimateDispersions( dxd )

## ----fitdiagnostics, dev='png', resolution=220---------------------------
plotDispEsts( dxd )

## ----testForDEU1,cache=TRUE----------------------------------------------
dxd = testForDEU( dxd )

## ----estFC,cache=TRUE----------------------------------------------------
dxd = estimateExonFoldChanges( dxd, fitExpToVar="Status")

## ----results1,cache=TRUE-------------------------------------------------
dxr1 = DEXSeqResults( dxd )

# Need a new accession here
figname <- file.path(figDir,"mboat7_expr.pdf")
pdf(figname)
targetacc <- "ENSG00000132773"
plotDEXSeq( dxr1, targetacc, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar="Status" )
dev.off()

objsavefile <- file.path(figDir,"dexseq_toe1.RData")
print(paste("Saving intermediate:",objsavefile))
#save(dxd, dxr1, file = objsavefile)
save.image() # Saves image
#unlink(obj.save.file) # Deletes file
#unlink(".RData") # Deletes Image
#load(objsavefile)

################################################################################

## ----results2,cache=TRUE-------------------------------------------------
elementMetadata(dxr1)$description


## ----tallyExons----------------------------------------------------------
table ( dxr1$padj < 0.1 )

dxr1.sub <- dxr1[ !is.na(dxr1$padj), ]

sigset <- dxr1.sub[ dxr1.sub$padj < 0.1,  ]
sigset <- sigset[ sort(sigset$padj), ]
print( head(sigset) )
#c("stat","padj","log2fold_control_knockout")
## ----tallyGenes----------------------------------------------------------
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )


## ----MvsA, dev='png', resolution=200-------------------------------------
figname <- file.path(figDir,"toe1_ma.pdf",sep="/")
pdf()
plotMA( dxr1, cex=0.8 )
dev.off()


## ----design--------------------------------------------------------------
sampleAnnotation(dxd)

## ----table2--------------------------------------------------------------
dxr1filt <- dxr1[ !is.na(dxr1$padj), ]
dxr1filt <- dxr1filt[ order(dxr1filt$padj), ]
sigset <- dxr1filt[ (dxr1filt$padj < 0.1) ,]

write.csv(as.data.frame(sigset), file=file.path(figDir,"dexseq_dxr1_sigset.csv"))

targetacc <- "ENSG00000140988"
figname <- file.path(figDir,paste("topGene",targetacc,"expr.pdf",sep="_"))
pdf(figname)
plotDEXSeq( dxr1, targetacc, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar="Status" )
dev.off()


targetacc <- "ENSG00000136682"
figname <- file.path(figDir,paste("topGene",targetacc,"expr.pdf",sep="_"))
pdf(figname)
plotDEXSeq( dxr1, targetacc, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar="Status" )
dev.off()

targetacc <- "ENSG00000128891"
figname <- file.path(figDir,paste("topGene",targetacc,"expr.pdf",sep="_"))
pdf(figname)
plotDEXSeq( dxr1, targetacc, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar="Status" )
dev.off()


write.csv(as.data.frame(dxr1), file=file.path(figDir,"dexseq_dxr1.csv"))


#q()
# More complicated Model
## ----formulas2-----------------------------------------------------------
formulaFullModel    =  ~ sample + exon + Family:exon + Status:exon
formulaReducedModel =  ~ sample + exon + Family:exon 


## ----estDisps_again, cache=TRUE, results='hide'--------------------------
dxd = estimateDispersions( dxd, formula = formulaFullModel )


## ----test_again, cache=TRUE----------------------------------------------
dxd = testForDEU( dxd, 
	reducedModel = formulaReducedModel, 
        fullModel = formulaFullModel )


## ----res_again-----------------------------------------------------------
dxr2 = DEXSeqResults( dxd )


## ----table2--------------------------------------------------------------
table( dxr2$padj < 0.1 )


## ----table3--------------------------------------------------------------
table( before = dxr1$padj < 0.1, now = dxr2$padj < 0.1 )


q()
## ----plot1, fig.height=8, fig.width=12-----------------------------------
#mboat7acc <- "ENSMUSG00000035596"
plotDEXSeq( dxr2, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


## ----checkClaim,echo=FALSE-----------------------------------------------
wh = (dxr2$groupID=="FBgn0010909")
stopifnot(sum(dxr2$padj[wh] < formals(plotDEXSeq)$FDR)==1)


## ----plot2, fig.height=8, fig.width=12-----------------------------------
plotDEXSeq( dxr2, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE,
   cex.axis=1.2, cex=1.3, lwd=2 )


## ----plot3, fig.height=8, fig.width=12-----------------------------------
plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, norCounts=TRUE,
   legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


## ----plot4, fig.height=8, fig.width=12-----------------------------------
plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, splicing=TRUE,
   legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


## ----DEXSeqHTML,cache=TRUE, eval=FALSE-----------------------------------
## DEXSeqHTML( dxr2, FDR=0.1, color=c("#FF000080", "#0000FF80") )


## ----para1,cache=TRUE,results='hide', eval=FALSE-------------------------
## BPPARAM = MultiCoreParam(workers=4)
## dxd = estimateSizeFactors( dxd )
## dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
## dxd = testForDEU( dxd, BPPARAM=BPPARAM)
## dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)


## ----alldeu, cache=TRUE--------------------------------------------------
dxr = DEXSeq(dxd)
class(dxr)


## ----buildExonCountSetLoadPacks,cache=TRUE, eval=FALSE-------------------
## library(GenomicRanges)
## library(GenomicFeatures)
## library(GenomicAlignments)


## ----buildExonCountSetDownloadAnno,cache=TRUE, eval=FALSE----------------
## hse = makeTranscriptDbFromBiomart( biomart="ensembl",
##    dataset="hsapiens_gene_ensembl" )


## ----buildExonCountSetDisjoin,cache=TRUE, eval=FALSE---------------------
## exonicParts = disjointExons( hse, aggregateGenes=FALSE )


## ----buildExonCountSet2FindBAMs,cache=TRUE, eval=FALSE-------------------
## bamDir = system.file( "extdata", package="parathyroidSE", mustWork=TRUE )
## fls = list.files( bamDir, pattern="bam$", full=TRUE )


## ----buildExonCountSet2FindBAMs2,cache=TRUE, eval=FALSE------------------
## bamlst = BamFileList( fls, index=character(), yieldSize=100000, obeyQname=TRUE )
## SE = summarizeOverlaps( exonicParts, bamlst, mode="Union", singleEnd=FALSE,
##    ignore.strand=TRUE, inter.feature=FALSE, fragments=TRUE )


## ----buildExonCountSet3,cache=TRUE, eval=FALSE---------------------------
## colData(SE)$condition = c("A", "A", "B")
## DEXSeqDataSetFromSE( SE, design= ~ sample + exon + condition:exon )


## ----acc-----------------------------------------------------------------
head( geneIDs(dxd) )
head( exonIDs(dxd) )


## ----grmethods-----------------------------------------------------------
interestingRegion = GRanges("chr2L", IRanges(start=3872658, end=3875302))
subsetByOverlaps( query=dxr, subject=interestingRegion )
findOverlaps( query=dxr, subject=interestingRegion )


## ----sessionInfo---------------------------------------------------------
sessionInfo()


