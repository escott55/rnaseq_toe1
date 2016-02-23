#!/usr/bin/Rscript

#An MA-plot (R. Dudoit et al. 2002) provides a useful overview for an experiment with a two-group comparison
genMAPlots <- function( res, topGene, tdir="results" ){

    ## ----plotma--------------------------------------------------------------
    pdf(paste(tdir,"plotma.pdf",sep="/"))
    plotMA(res, ylim=c(-5,5))
    dev.off()

    ## ----plotmalabel---------------------------------------------------------
    #pdf("results/plotma_lab.pdf")
    pdf(paste(tdir,"plotMA_LFC.pdf",sep="/"))
    plotMA(resLFC1, ylim=c(-5,5))
    topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
    with(resLFC1[topGene, ], {
      points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
      text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
    })
    dev.off()
}

