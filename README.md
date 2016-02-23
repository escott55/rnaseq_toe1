TOE1 Analysis
========
This repository includes three separate RNA-seq analyses, produced to
investigate different possible outcomes of a splicing defect.

TODO
-----------------
 1. Make a full processing pipeline for the samples
 2. Adjust colors in splicing plot by log-fold-change
 3. Perform a tRNA analysis
 4. Do an analysis based on intron location

Expression
-----------------
This code base uses DESeq2 to calculate normalized expression values for
all genes and identifies significantly differentially expressed genes by
a Case vs Control comparison.

Intron Inclusion
-----------------
Emulates an analysis from the U2 paper (should add link here). Briefly we
calculate read counts for intronic and extronic regions using HT-seq. Calculate
FPKM values and convert these to the proportion of transcript inclusion for
each intron, as a function of the max value from adjacent exons.

Exon Usage
----------
Run an analysis to identify exons with significantly differential expression.

Prerequisite python packages
----------------------------
 1. Rpy2 with ggplot2 installed as part of the R installation
 2. Pandas for dataframes

Running The Analysis
--------------------

Expression Analysis:
``` R
cd scripts/expression/
./run_expression_analysis.R
```

Intron Inclusion Analysis:
``` bash
python scripts/introninc/run_region_read_counts.py
python scripts/introninc/run_intron_inclusion_analysis.py
```
