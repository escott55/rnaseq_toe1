#!/usr/bin/Rscript

library(ggplot2)

cntDistHeatmap <- function(rld, dds, tdir="results", typecol="type",
                           shapecol="type") {

    ## ------------------------------------------------------------------------
    sampleDists <- dist( t( assay(rld) ) )
    #print(sampleDists)

    ## ------------------------------------------------------------------------
    library("pheatmap")
    library("RColorBrewer")

    ## ----distheatmap, fig.width=6, fig.height=4.5----------------------------
    sampleDistMatrix <- as.matrix( sampleDists )
    #rownames(sampleDistMatrix) <- paste( rld$type, rld$cell, sep="-" )
    rownames(sampleDistMatrix) <- paste( rld[[typecol]], rld$id, sep="-" )
    colnames(sampleDistMatrix) <- NULL
    currcolors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    tfile <- paste(tdir,"sampleDistHeatmap.pdf",sep="/")
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=currcolors, filename=tfile)

    ## ----mdsrlog, fig.width=6, fig.height=4.5--------------------------------
    mdsData <- data.frame(cmdscale(sampleDistMatrix))
    mds <- cbind(mdsData, as.data.frame(colData(rld)))

    p <- (ggplot(mds, aes_string("X1","X2",color=typecol,shape=shapecol)) +
          geom_point(size=3))
    pdf(paste(tdir,"ggSampleMDS.pdf",sep="/"))
    print(p)
    dev.off()
}

poisDistHeatmap <- function(rld, dds, tdir="results", typecol="type",
                            shapecol="type") {
    ## ------------------------------------------------------------------------
    library("PoiClaClu")
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
    p <- (ggplot(mdsPois, aes_string("X1","X2",color=typecol,shape=typecol)) +
          geom_point(size=3))
    pdf(paste(tdir,"ggSampleMDSPois.pdf",sep="/"))
    print(p)
    dev.off()
} # END poisDistHeatmap  

samplePCAPlots <- function(rld, tdir="results", typecol="type", shapecol="type") {
    ## ------------------------------------------------------------------------
    if( typecol == shapecol ){
        pcadf <- plotPCA(rld, intgroup = c( typecol, "id"), returnData=TRUE)
    }else{
        pcadf <- plotPCA(rld, intgroup = c( typecol, shapecol, "id"),
                         returnData=TRUE)
    }
    percentVar <- round(100 * attr(pcadf, "percentVar"))

    ## ----plotpca, fig.width=6, fig.height=4.5--------------------------------
    #pdf(paste(tdir,"samplePCA.pdf",sep="/"))
    #plotPCA(rld, intgroup = c("type", "id"))
    #dev.off()

    ## ----ggplotpca, fig.width=6, fig.height=4.5------------------------------
    tfile <- paste(tdir,"ggSamplePCA.pdf",sep="/")
    print(paste("Making file:",tfile,sep=" "))
    p <- (ggplot(pcadf, aes_string("PC1", "PC2", color=typecol, shape=shapecol)) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")))
    pdf(tfile)
    print(p)
    dev.off()
}

