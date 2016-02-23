#!/usr/bin/Rscript

pltCntsPdf <- function(dds, topGene, tdir="results", typecol="type",
                       shapecol="type"){
    
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
        data <- plotCounts(dds, gene=topGene, intgroup=c(typecol), returnData=TRUE)
    } else {
        data <- plotCounts(dds, gene=topGene, intgroup=c(typecol,shapecol),
                           returnData=TRUE)
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
    p <- (ggplot(data, aes_string(x=typecol, y="count", fill=typecol, shape=shapecol)) +
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

########################## END TOP GENE ###################################

