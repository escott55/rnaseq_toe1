#!/usr/bin/Rscript

# This script runs an standard case/control RNAseq expression analysis

############################################################
# TO DO
# make parameters to specify the different status types
# Pass the correct columns for family and status
# Improve user input methods
# Add a dexseq exon usage analysis
############################################################

####################################################################################

###################################
# Installation of required packages
###################################
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

#########################
# Load required libraries
#########################
library("BiocStyle")
library("knitr")
library("rmarkdown")
library("Rsamtools")
library("GenomicFeatures")
library("ggplot2")
## ------------------------------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")
library("genefilter")
## ------------------------------------------------------------------------
library("DESeq2")

#######################
# Load local namespaces
#######################
source("deseq_load.R")
source("deseq_func.R")

parseSampleInfo <- function(datasetname)
{   # Return sampletable and bamfile associated for current dataset
    ###################
    # Define User input
    ###################
    maindir <- "/home/escott/workspace/rnaseq_analysis/toe1"

    if (!file.exists(maindir)) {
        stop("Error: base directory doesnt exist")
    }
    if (datasetname == "T283") {
        # T283 data
        csvfile <- file.path(maindir,"raw/t283/SampleTable.txt")
        species <- "Human"
        contrast <- c("Status","A","F","U")
    } else if (datasetname == "Zebrafish") {
        # Zebrafish data
        csvfile <- file.path(maindir,"raw/zebrafish/SampleTable.txt")
        species <- "Zebrafish"
        contrast <- c("Status","A","U","C")
    } else if (datasetname == "NewNPC") {
        # New NPC data
        csvfile <- file.path(maindir,"raw/newnpc/SampleTable.txt")
        species <- "Human"
        contrast <- c("Status","A","U")
        # Need to filter samples by batch == 2
    } else if (datasetname == "OldNPC") {
        # New NPC data
        csvfile <- file.path(maindir,"raw/oldnpc/SampleTable.txt")
        species <- "Human"
        contrast <- c("Status", "A", "U")
    } else if (datasetname == "fibroblasts") {
        csvfile <- file.path(maindir,"raw/fibroblast/SampleTable.txt")
        species <- "Human"
        contrast <- c("Status","A","U")
    } else {
        stop("Error: Unknown datasetname provided")
    }

    resultsdir <- file.path(maindir,"results/expression",datasetname)

    # Use species to identify GTF file
    if (species == "Human") {
        gtffile <- file.path(maindir, "resources",
                             "Homo_sapiens.GRCh37.b37.protcoding.gtf")
    }else if(species == "Mouse") {
        gtffile <- file.path(maindir, "resources",
                             "gencode.vM7.annotation_encode_ccds.gtf")
    }else if(species == "Zebrafish") {
        gtffile <- file.path(maindir, "resources",
                             "Danio_rerio.GRCz10.83.protcoding.gtf")
    }else{
        stop("Error: no known gtffile for identified species")
    }

    # Parse Sample Table information
    sampleTable <- read.csv(csvfile, row.names=1, sep="\t", comment.char="#")
    sampleTable$id <- row.names(sampleTable)

    sampleTable$Family <- as.factor(sampleTable$Family)

    # Todo - limit results to target contrasts
    # sampleTable <- sampleTable[ sampleTable$Status %in% c("U","C"), ]
    # Todo - filter by batch if necessary
    # Retrieve bamfiles from SampleTable
    rawdatadir <- dirname(csvfile)
    bamfiles <- as.vector(file.path(rawdatadir, sampleTable$File))
    #datadir <- dirname(bamfiles[0]) # Figure out original data directory

    #########################
    # Check Input consistency
    #########################
    if (!(file.exists(resultsdir))) {
        # Create the results directory if it doesnt exist
        dir.create(resultsdir)
    }

    if (!all(file.exists(bamfiles))) {
        print(cbind(bamfiles, file.exists(bamfiles)))
        stop("Error: not all bamfiles exist!")
    }
    return(list(dataset = datasetname, species = species, gtffile = gtffile,
                bamfiles = bamfiles, sampleTable = sampleTable,
                resultsdir = resultsdir, contrast = contrast))
}   # parseSampleInfo


runmain <- function()
{   # Full expression analysis pipeline
    # Ex. datasets - Zebrafish, NewNPC, Fibroblasts, T283, OldNPC
    datasetname <- "Zebrafish"
    dsmeta <- parseSampleInfo(datasetname)
    # Todo - should probably just use the object instead of the following
    gtffile <- dsmeta$gtffile
    bamfiles <- dsmeta$bamfiles
    resdir <- dsmeta$resultsdir
    sampleTable <- dsmeta$sampleTable
    contrast <- dsmeta$contrast
    curspecies <- dsmeta$species

    ######################
    print("Loading Data")
    ######################
    se <- exprload.summarizeBams(gtffile, bamfiles, resdir, sampleTable)

    # Todo - figure out how to parameterize a formula
    dds <- exprload.DESeqNorm(se, resdir, design = ~Status)
    rm(se) # We dont need this object any longer

    ######################
    print("Generating QC plots")
    ######################
    func.plt.cntDistHeatmap(dds, resdir, typecol = "Status",
                            shapecol = "Family")

    # Redundant unless we need multiple distance metrics
    #func.plt.poisDistHeatmap(dds, resdir, typecol="Status", shapecol="Family")

    func.plt.samplePCA(dds, resdir, typecol = "Status", shapecol = "Family")

    res <- func.identifySignificantSet(dds)

    func.plt.PvalueDistribution(resdir, res)

    ######################
    print("Identify top genes")
    ######################
    topGene <- rownames(res)[which.min(res$padj)]
    func.plt.countsByGrp(dds, topGene, resdir, typecol = "Status",
                         shapecol = "Family")

    func.plt.MA(res, dds, contrast, topGene, resdir)
    print("done with MA")

    func.plt.clustHeatmap(dds, resdir)

    func.plt.meanNormalizedCounts(res, resdir)

    resAnnot <- func.annotateResults(res, curspecies)

    ######################
    print("Write results to file")
    ######################
    func.writeResultsFile(res, resdir)

    # Subset data for only complete symbols if required
    resSub <- subset(resAnnot, complete.cases(symbol))
    resDF <- as.data.frame(resSub)

    func.plt.geneTracks(dds, topGene, resdir, curspecies)
}

runmain() # Full workflow

# END analysis
