#!/usr/bin/Rscript

#library("BiocStyle")
#library("knitr")
#library("rmarkdown")
#library("Rsamtools")
#library("GenomicFeatures")
## ------------------------------------------------------------------------
#library("GenomicAlignments")
#library("BiocParallel")
## ------------------------------------------------------------------------

library("DESeq2")
library("ggplot2")

func = new.env()

################################################################################
# Data Manipulation functions
################################################################################

func.identifySignificantSet <- function(dds)
{   # Summarizes DEseq object to results
    # Returns res object
    print("Running func.identifySignificantSet")
    res <- results(dds)

    # Summary functions
    #mcols(res, use.names=TRUE)
    #summary(res)
    #res.05 <- results(dds, alpha=.05)
    #table(res.05$padj < .05)
    # Log fold change threshold
    #table(resLFC1$padj < 0.1)
    #results(dds, contrast=c("Status", "U", "A"))
    #sum(res$pvalue < 0.05, na.rm=TRUE)
    #sum(!is.na(res$pvalue))

    # Significant set by cutoff
    #resSig <- subset(res, padj < 0.1)
    #head(resSig[ order(resSig$log2FoldChange), ])
    return(res)
}

func.annotateResults <- function(res, species="Human")
{   # Returns annotated res object
    require("AnnotationDbi")
    print("Running func.annotateResults")

    if (species == "Human") {
        library("org.Hs.eg.db")
        annotdb <- org.Hs.eg.db
    } else if (species == "Mouse") {
        library("org.Mm.eg.db")
        annotdb <- org.Mm.eg.db
    } else if (species == "Fly") {
        library("org.Dr.eg.db")
        annotdb <- org.Dr.eg.db
    } else {
        stop("Error: unknown species provided")
    }

    columns(annotdb)
    res$symbol <- mapIds(annotdb,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
    res$entrez <- mapIds(annotdb,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

    resOrdered <- res[order(res$padj),]
    head(resOrdered)
    return(resOrdered)
}

func.writeResultsFile <- function(res, resdir)
{   # Write significant subset and full dataset to file
    ## ----eval=FALSE----------------------------------------------------------
    #ressig <- subset(res, padj < 0.1)
    resDF <- as.data.frame(subset(res,padj<0.05)) #[1:100,]
    write.csv(res, file=file.path(resdir,"sigset.csv"))

    write.csv(as.data.frame(res), file=file.path(resdir,"fullset.csv"))
}

################################################################################
# Plotting functions
################################################################################

func.plt.QC <- function(dds, resdir)
{
    print("Running func.plt.QC")
    rld <- rlog(dds, blind=FALSE)
    head(assay(rld), 3)

    pdf(file.path(resdir,"qc2.pdf"))
    par( mfrow = c( 1, 2 ) )
    plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
         pch=16, cex=0.3)
    plot(assay(rld)[,1:2],
         pch=16, cex=0.3)
    dev.off()

}

func.plt.cntDistHeatmap <- function(dds, tdir="results", typecol="type",
                                    shapecol="type")
{   # Generate count based heatmap of data
    require("pheatmap")
    require("RColorBrewer")

    print("Running func.plt.cntDistHeatmap")
    rld <- rlog(dds, blind=FALSE)

    sampleDists <- dist(t(assay(rld)))
    #print(sampleDists)

    ## ----distheatmap, fig.width=6, fig.height=4.5----------------------------
    sampleDistMatrix <- as.matrix(sampleDists)
    #rownames(sampleDistMatrix) <- paste(rld$type, rld$cell, sep="-")
    rownames(sampleDistMatrix) <- paste(rld[[typecol]], rld$id, sep="-")
    colnames(sampleDistMatrix) <- NULL
    currcolors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    tfile <- paste(tdir, "sampleDistHeatmap.pdf", sep="/")
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=currcolors, filename=tfile)

    ## ----mdsrlog, fig.width=6, fig.height=4.5--------------------------------
    mdsData <- data.frame(cmdscale(sampleDistMatrix))
    mds <- cbind(mdsData, as.data.frame(colData(rld)))
    mds[,c(shapecol)] <- factor(mds[,c(shapecol)])

    p <- (ggplot(mds, aes_string("X1", "X2", color=typecol, shape=shapecol)) +
          geom_point(size=3))
    pdf(paste(tdir,"ggSampleMDS.pdf",sep="/"))
    print(p)
    dev.off()
}

func.plt.poisDistHeatmap <- function(dds, tdir="results", typecol="type",
                                     shapecol="type")
{   # Plot Distance Heaptmap using Poisson Distance
    require("PoiClaClu")

    print("Running func.plt.poisDistHeatmap")
    rld <- rlog(dds, blind=FALSE)
    poisd <- PoissonDistance(t(counts(dds)))

    ## ----poisdistheatmap, fig.width=6, fig.height=4.5------------------------
    tfile <- paste(tdir,"poisonheatmap.pdf",sep="/")
    samplePoisDistMatrix <- as.matrix( poisd$dd )
    rownames(samplePoisDistMatrix) <- paste( rld[[typecol]], rld$id, sep="-" )
    colnames(samplePoisDistMatrix) <- NULL
    currcolors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(samplePoisDistMatrix,
             clustering_distance_rows=poisd$dd,
             clustering_distance_cols=poisd$dd,
             col=currcolors, filename=tfile)

    ## ----mdspois, fig.width=6, fig.height=4.5--------------------------------
    mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
    mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
    mdsPois[,c(shapecol)] <- factor(mdsPois[,c(shapecol)])
    p <- (ggplot(mdsPois, aes_string("X1","X2",color=typecol,shape=typecol)) +
          geom_point(size=3))
    pdf(paste(tdir,"ggSampleMDSPois.pdf",sep="/"))
    print(p)
    dev.off()
} # END poisDistHeatmap

func.plt.samplePCA <- function(dds, tdir="results", typecol="type",
                                    shapecol="type")
{   # Make PCA plots for all samples
    print("Running func.plt.samplePCA")

    rld <- rlog(dds, blind=FALSE)

    # Make DEseq version of PCA plot
    if( typecol == shapecol ){
        pcadf <- plotPCA(rld, intgroup = c( typecol, "id"), returnData=TRUE)
    }else{
        pcadf <- plotPCA(rld, intgroup = c( typecol, shapecol, "id"),
                         returnData=TRUE)
        pcadf[,c(shapecol)] <- factor(pcadf[,c(shapecol)])
    }

    percentVar <- round(100 * attr(pcadf, "percentVar"))

    # Make ggplot2 PCA plot
    tfile <- paste(tdir,"ggSamplePCA.pdf",sep="/")
    print(paste("Making file:",tfile,sep=" "))
    p <- (ggplot(pcadf, aes_string("PC1", "PC2", color=typecol,
                                   shape=shapecol)) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")))
    pdf(tfile)
    print(p)
    dev.off()
}

func.plt.countsByGrp <- function(dds, topGene, tdir="results",
                                 typecol="type", shapecol="type")
{   # Plot of normalized counts for an individual Gene
    print("Running func.plt.countsByGrp")

    # Should add require statements here
    fname <- paste("topGeneCnts_",topGene,".pdf",sep="")
    #if( typecol == shapecol ){
    #}else{
        #p <- plotCounts(dds, gene=topGene, intgroup=c(typecol, shapecol))
    #}
    pdf(paste(tdir,fname,sep="/"))
    plotCounts(dds, gene=topGene, intgroup=c(typecol))
    dev.off()

    ## ----ggplotcountsjitter, fig.width=6, fig.height=4.5---------------------

    if( typecol == shapecol ){
        data <- plotCounts(dds, gene=topGene, intgroup=c(typecol),
                           returnData=TRUE)
    } else {
        data <- plotCounts(dds, gene=topGene, intgroup=c(typecol,shapecol),
                           returnData=TRUE)
        data[,c(shapecol)] <- factor(data[,c(shapecol)])
    }
    #pdf(paste(tdir,"ggTopGeneJitter.pdf",sep="/"))
    #ggplot(data, aes(x=type, y=count, color=type)) +
      #scale_y_log10() +
      #geom_point(position=position_jitter(width=.1,height=0), size=3) +
      #ggtitle(topGene)
    #dev.off()

    ## ----ggplotcountsdot, fig.width=6, fig.height=4.5------------------------
    #pdf("results/ggdotplot.pdf")
    fname <- paste("ggTopGeneDot_",topGene,".pdf",sep="")
    p <- (ggplot(data, aes_string(x=typecol, y="count", fill=typecol,
                                  shape=shapecol)) +
          scale_y_log10() +
          geom_dotplot(binaxis="y", stackdir="center") +
          ggtitle(topGene))
    pdf(paste(tdir,fname,sep="/"))
    print(p)
    dev.off()

    ## ----ggplotcountsgroup, fig.width=6, fig.height=4.5----------------------
    #pdf("results/cntgrp.pdf")
    #pdf(paste(tdir,"ggTopGeneLine.pdf",sep="/"))
    #ggplot(data, aes(x=type, y=count, color=type, group=type)) +
      #scale_y_log10() + geom_point(size=3) + geom_line() + ggtitle(topGene)
    #dev.off()
}

func.plt.PvalueDistribution <- function(resdir, res)
{   # Make a distribution plot of P-values
    print("Running func.plt.PvalueDistribution")

    Pvalues <- res$pvalue[res$baseMean > 1]
    pdf(file.path(resdir,"plothist.pdf"))
    hist(Pvalues, breaks=0:20/20, col="grey50",
         border="white")
    dev.off()
}

func.plt.MA <- function(res, dds, contrast, topGene, tdir="results")
{   # Generate MA plot from results object
    # An MA-plot (R. Dudoit et al. 2002) provides a useful overview for an
    # experiment with a two-group comparison
    print("Running func.plt.MA")

    pdf(paste(tdir,"plotma.pdf",sep="/"))
    plotMA(res, ylim=c(-5,5))
    dev.off()

    resLFC1 <- results(dds, contrast, lfcThreshold=1) # Log fold threshold of 1

    print(head(resLFC1))
    #pdf("results/plotma_lab.pdf")
    pdf(paste(tdir, "plotMA_LFC.pdf", sep="/"))
    plotMA(resLFC1, ylim=c(-5,5))
    topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
    with(resLFC1[topGene, ], {
         points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
         text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
    })
    dev.off()
}

func.plt.clustHeatmap <- function(dds, resdir)
{   # Calculate Sample distance using Gene Expression
    require("genefilter")
    print("Running func.plt.clustHeatmap")

    rld <- rlog(dds, blind=FALSE)
    topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

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
}

func.plt.meanNormalizedCounts <- function(res, resdir)
{   # Distribution of normalized genes
    print("Running func.plt.meanNormalizedCounts")

    resLFC1 <- results(dds, lfcThreshold=1) # Log fold threshold of 1
    qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
    bins <- cut(resLFC1$baseMean, qs)
    levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))

    ratios <- tapply(resLFC1$pvalue, bins,
                     function(p) mean(p < .05, na.rm=TRUE))
    pdf(file.path(resdir,"meanNormCounts.pdf"))
    barplot(ratios,
            xlab="mean normalized count",
            ylab="ratio of small p values")
    dev.off()
}

func.plt.geneTracks <- function(dds, topGene, resdir, species="Human")
{   # Plot a Gene Tracks heatmap
    require("Gviz")
    print("Running func.plt.geneTracks")

    resGR <- results(dds, lfcThreshold=1, format="GRanges")

    if (species == "Human") {
        library("org.Hs.eg.db")
        annotdb <- org.Hs.eg.db
    }else if (species == "Mouse") {
        library("org.Mm.eg.db")
        annotdb <- org.Mm.eg.db
    }else if (species == "Fly") {
        library("org.Dr.eg.db")
        annotdb <- org.Dr.eg.db
    }else{
        stop("Error: provided species is not known")
    }

    columns(annotdb)
    resGR$symbol <- mapIds(annotdb, names(resGR), "SYMBOL", "ENSEMBL")

    ## ------------------------------------------------------------------------
    library("Gviz")

    ## ------------------------------------------------------------------------
    window <- resGR[topGene] + 1e6
    strand(window) <- "*"
    resGRsub <- resGR[resGR %over% window]
    naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
    resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

    ## ------------------------------------------------------------------------
    sig <- factor(ifelse(resGRsub$padj < .1 &
                         !is.na(resGRsub$padj),"sig","notsig"))

    ## ----gvizplot------------------------------------------------------------
    options(ucscChromosomeNames=FALSE)
    g <- GenomeAxisTrack()
    a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
    d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
                   type="h", name="log2 fold change", strand="+")

    pdf(file.path(resdir,"groupTracks.pdf"))
    plotTracks(list(g,d,a),
               groupAnnotation="group",
               notsig="grey",
               sig="hotpink")
    dev.off()
}

# Reattach functions on source
while("func" %in% search()){
    detach("func")
}
attach(func)
########################## END FUNC ###################################


