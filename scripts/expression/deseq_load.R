#!/usr/bin/Rscript

library("BiocStyle")
library("knitr")
library("rmarkdown")
library("Rsamtools")
library("GenomicFeatures")
## ------------------------------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")
## ------------------------------------------------------------------------
library("DESeq2")

exprload = new.env()

################################################################################
# Loading functions
# probably should split into load and clean
# regenerate flag is used to force an updated copy of the object to be made
################################################################################

exprload.generateDESeqObject <- function(gtffile, bamfiles, resdir, sampleTable,
                                         design=~Status)
{
    print("Loading Data")
    se <- exprload.summarizeBams(gtffile, bamfiles, resdir, sampleTable)

    dds <- exprload.DESeqNorm(se, resdir, design=~Status)
}

exprload.summarizeBams <- function(gtffile, bamfiles, resdir, sampleTable,
                                   mode="Union", statuscol="Status",
                                   singleEnd=FALSE, debug=FALSE,
                                   regenerate=FALSE)
{   # Summarize bamfiles into se object
    print("Running exprload.summarizeBams")

    ######################
    # Verify correct input
    ######################
    if (!(file.exists(gtffile))) {
        stop("Error: gtffile does not exist")
    }

    if (length(bamfiles) == 0) {
        stop("Error: no bamfiles provided")
    }

    if (!(statuscol %in% colnames(sampleTable))) {
        stop("Error: statuscol doesnt exist in sampleTable")
    }

    #####################################################
    # if the se object exists already, load and return it
    #####################################################
    se <- exprload.loadObjectIfPresent(resdir, "se")
    if (!is.null(se) & !regenerate) {
        print("Warning: loading precomputed object")
        return(se)
    }

    #####################################################
    # Make se object if doesnt exist
    #####################################################
    bamlistobj <- BamFileList(bamfiles, yieldSize=1000000)

    txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())

    ebg <- exonsBy(txdb, by="gene") # Exons by gene

    #register(SerialParam()) # Load serial params?

    se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE)

    if (debug) {
        print("Summarized overlaps debug")
        print(se) # Summary of data object
        print(dim(se)) # Dimensions of our matrix
        print(assayNames(se))
        print(head(assay(se), 3))
        print(colSums(assay(se)))
        #rowRanges(se)
        #str(metadata(rowRanges(se)))
    }

    # add sample information
    colData(se) <- DataFrame(sampleTable)

    # Reorder levels to put affecteds first
    se$Status <- relevel(se$Status, "A")

    ## More optional debugging
    #round( colSums(assay(se)) / 1e6, 1 )
    #colData(se)
    exprload.saveObjects(se, "se", resdir) # Save SE object
    return(se)
} # END summarizeBams

exprload.DESeqNorm <- function(se, resdir, design, regenerate=FALSE)
{   # Normalize se object into DESeq counts
    print("Running exprload.DESeqNorm")

    require("DESeq2")
    # Assert input was passed correctly
    if (!file.exists(resdir)) {
        stop("Error: data directory doesnt exist")
    }

    ######################################################
    # if the dds object exists already, load and return it
    ######################################################
    dds <- exprload.loadObjectIfPresent(resdir, "dds")
    if (!is.null(dds) & !regenerate) {
        print("Warning: loading precomputed object")
        return(dds)
    }

    dds <- DESeqDataSet(se, design = ~ Status)

    # If a matrix is used as input
    #design = ~ Family + Status + Family:Status
    #ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                      #colData = coldata,
                                      #design = ~ Status)
    #nrow(dds)
    dds <- dds[ rowSums(counts(dds)) > 1, ]
    #nrow(dds)
    dds <- estimateSizeFactors(dds)

    dds <- DESeq(dds)

    row.names(dds) <- substr(row.names(dds), 0, 18)

    normcnts <- head(log2(counts(dds, normalized=TRUE)))

    exprload.saveObjects(dds, "dds", resdir) # Save SE object
    return(dds)
} # END DESeqNorm

exprload.saveObjects <- function( tobj, objectname, targetdir )
{   # Save an object provided the object and name
    print(paste("Saving object:",objectname))
    #objsavefile <- file.path(targetdir, "deseq2_toe1.RData")
    #save(dxd, dxr1, file = obj.save.file)
    #save.image(objsavefile) # Saves image

    if (is.null(tobj)){
        stop("Error: passed object is null")
    }

    objsavefile <- file.path(targetdir, paste("deseq2_toe1_", objectname,
                                              ".RData", sep=""))
    print(paste("Saving object to file:",objsavefile))
    saveRDS(tobj, objsavefile)
}

exprload.loadObjectIfPresent <- function(targetdir, objectname)
{   # Return the saved object if it exists
    objsavefile <- file.path(targetdir, paste("deseq2_toe1_", objectname,
                                              ".RData", sep=""))

    if (file.exists(objsavefile)) {
        print(paste("Loading object:", objectname))
        tobj <- readRDS(objsavefile)
        return(tobj)
    }

    return(NULL)
}

# Reattach functions on source
while("exprload" %in% search()){
    detach("exprload")
}
attach(exprload)
########################## END LOAD ###################################
