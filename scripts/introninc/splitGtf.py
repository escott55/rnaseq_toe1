#!/usr/bin/python

import os
import sys
import csv
import subprocess

import pandas as pd
sys.path.append(os.path.expanduser("~/workspace/custom_utils/general"))
import housekeeping as hk

# Globals
GTFCOLS = ["chrom", "type", "class", "start", "end", "score", "strand",
           "frame", "attribute"]


def parseAttributes(attstr):
    newdict = {}
    for attpair in attstr.split(";"):
        attpair = attpair.strip()
        if attpair.count(" ") > 0:
            key, val = [x.strip().replace('"', '')
                        for x in attpair.split(" ", 1)]
            newdict[key] = val

    # print newdict
    return pd.Series(newdict)
# END parseAttributes


def calcIntrons(geneinfo, starts, ends):
    if len(starts) < 2:
        # print "Error: No introns?", starts, ends
        return None

    print "Starts:", starts
    print "Ends:", ends
    introns = []
    for i in range(len(ends)-1):
        att = ('gene_id "%s"; gene_name "%s"; transcript_id %s; '
               'intron_id "%d";' %
               (geneinfo[1], geneinfo[3], geneinfo[4], i))
        newintron = (geneinfo[0], "protein_coding", "exon",
                     str(int(ends[i])+1), str(int(starts[i+1])-1), ".",
                     geneinfo[2], ".", att)
        # print newintron
        # diff = int(ends[i])+1 + int(starts[i+1])-1
        # print "Diff:",diff
        introns.append(newintron)

    return pd.DataFrame(introns, columns=GTFCOLS)
    # return DataFrame(introns, columns=["start","end"])
# END calcIntrons


def testGtf():
    resourcedir = "/home/escott/workspace/toe1/resources"
    tfile = os.path.join(resourcedir, "Homo_sapiens.GRCh37.b37.introns.gtf")
    testdf = pd.read_csv(tfile, sep="\t", header=None, names=GTFCOLS, nrows=50)

    regionlens = testdf["end"] - testdf["start"]
    print regionlens.head()

    exonsgtf = os.path.join(resourcedir, "Homo_sapiens.GRCh37.b37.exons.gtf")
    testdf = pd.read_csv(exonsgtf, sep="\t", header=None, names=GTFCOLS,
                         nrows=50)

    regionlens = testdf["end"] - testdf["start"]
    print regionlens.head()
# END testGtf


def makeGffFromGtf(targetgtf):
    dexseqscript = os.path.join("/usr/local/lib/R/site-library", "DEXSeq",
                                "python_scripts/dexseq_prepare_annotation.py")
    newgff = targetgtf[:-3] + "gff"

    cmd = ["python", dexseqscript, targetgtf, newgff]

    print " ".join(cmd)
    out = subprocess.check_output(cmd)
    print out
    assert os.path.exists(newgff)
    return newgff
# END makeGffFromGtf


def oneGeneTest(regions_sub, geneid="ENSG00000108852"):
    generegions = regions_sub[regions_sub.gene_id == geneid]
    print "Total regions:", len(generegions)
    print generegions.groupby("gene_name").size()
    print generegions.groupby("gene_source").size()
    print generegions.groupby("gene_biotype").size()
    print generegions.groupby("transcript_id").size()
    print generegions.groupby("transcript_name").size()
    print generegions.groupby("gene_source").size()
    print generegions.groupby("tag").size()
    print generegions.groupby("exon_number").size()
    print generegions.groupby("ccds_id").size()
    print "Missing ccds id:", sum(generegions.ccds_id.isnull())
    print "Missing transcript ids:", sum(generegions.transcript_id.isnull())
# END oneGeneTest


def parseRegionsGtf(regionclass, gtffile):
    # For testing use - nrows=50
    regions = pd.read_csv(gtffile, sep="\t", header=None, names=GTFCOLS,
                          comment="#", skiprows=5, compression="gzip")
    print regions.head()
    # Subset only certain regions
    print "Full set of regions:", len(regions)
    print "Set of classes available:"
    print regions.groupby("class").size()

    # We have to use exons here otherwise the DEXSEQ script doesnt
    # identify any regions. We could change the name to "exon" after
    # selection, but this hasnt been implemented yet.
    if regionclass == "coding":
        # Use only one class of regions - Ex. CDS, exon
        regions = regions[(regions["class"] == "exon")]
        print "After filtering for exons:", len(regions)
    else:
        print "Error: unknown region class"

    # Filter for Canonical regions
    regions = regions.sort_values(by=["chrom", "start", "end",
                                      "attribute"])
    regions_sub = regions.drop_duplicates(subset=["chrom", "start",
                                                  "end"])
    # Creat annotation set and merge with raw regions
    attdf = pd.DataFrame(regions_sub.attribute.apply(parseAttributes))

    # del regions["attribute"]  # I keep this column in case I make exonsgtf

    regionsannot = pd.concat([regions_sub, attdf], axis=1)

    oneGeneTest(regionsannot)

    # Why do I remove empty transcript_ids?
    # CCDS ID seems null when the region isn't verified as coding
    # (regionsannot["transcript_id"].notnull()) and
    # r_annot_sub = regionsannot[(regionsannot["ccds_id"].notnull())]
    # print "Number of Genes", len(regionsannot.gene_name.unique())
    # print "After filtering Genes", len(r_annot_sub.gene_name.unique())
    return regionsannot
# END parseRegionsGtf


def retrieveGTFs(resourcedir, regionclass, species):
    if regionclass == "coding" and species == "Human":
        gtffile = os.path.join(resourcedir, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.protcoding.gtf.gz")
        exonsgtf = os.path.join(resourcedir, "human_b37_regions",
                                "Homo_sapiens.GRCh37.b37.exons.gtf")
        intronsgtf = os.path.join(resourcedir, "human_b37_regions",
                                  "Homo_sapiens.GRCh37.b37.introns.gtf")

    elif regionclass == "coding" and species == "Zebrafish":
        gtffile = os.path.join(resourcedir, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.protcoding.gtf.gz")
        exonsgtf = os.path.join(resourcedir, "zebra_z10_regions",
                                "Danio_rerio.GRCz10.83.exons.gtf")
        intronsgtf = os.path.join(resourcedir, "zebra_z10_regions",
                                  "Danio_rerio.GRCz10.83.introns.gtf")

    elif regionclass == "noncoding" and species == "Human":
        # Make non-coding regions set
        gtffile = os.path.join(resourcedir, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.noncode.gtf.gz")
        exonsgtf = os.path.join(resourcedir, "human_b37_regions",
                                "Homo_sapiens.GRCh37.b37.ncexons.gtf")
        intronsgtf = os.path.join(resourcedir, "human_b37_regions",
                                  "Homo_sapiens.GRCh37.b37.ncintrons.gtf")

    elif regionclass == "noncoding" and species == "Zebrafish":
        # Make non-coding regions set
        gtffile = os.path.join(resourcedir, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.noncode.gtf.gz")
        exonsgtf = os.path.join(resourcedir, "zebra_z10_regions",
                                "Danio_rerio.GRCz10.83.ncexons.gtf")
        intronsgtf = os.path.join(resourcedir, "zebra_z10_regions",
                                  "Danio_rerio.GRCz10.83.ncintrons.gtf")
    else:
        print "Error: unrecognized parameters"
        sys.exit(1)

    return gtffile, exonsgtf, intronsgtf
# END retrieveGTFs


if __name__ == "__main__":
    resourcedir = os.path.abspath("../../resources")

    regionclass = "coding"
    species = "Human"
    gtffile, exonsgtf, intronsgtf = retrieveGTFs(resourcedir, regionclass,
                                                 species)
    # intronsgff = intronsgtf[:-3] + "gff"
    # exonsgff = exonsgtf[:-3] + "gff"

    # tdir = os.path.abspath("../../resources")
    regionsannot = parseRegionsGtf(regionclass, gtffile)

    # Write out Exons GTF
    if not os.path.exists(exonsgtf):
        # Make a new exonsgtf if one doesnt exists
        regionsannot[GTFCOLS].to_csv(exonsgtf, sep="\t", index=False,
                                     header=False, quoting=csv.QUOTE_NONE)

    exongff = makeGffFromGtf(exonsgtf)

    sys.exit(1)

    # Identify Intronic Regions
    allintrons = []
    gidcols = ["chrom", "gene_id", "strand", "gene_name", "transcript_id"]
    for geneinfo, gdata in regionsannot.groupby(gidcols):
        exonstarts = sorted([start for start in
                             gdata.start[gdata["class"] == "exon"]])
        exonends = sorted([end for end in
                           gdata.end[gdata["class"] == "exon"]])
        introns = calcIntrons(geneinfo, exonstarts, exonends)
        allintrons.append(introns)

    allintrons = pd.concat(allintrons)

    # Make Introns GTF
    print allintrons.head()
    allintrons.to_csv(intronsgtf, sep="\t", index=False, header=False,
                      quoting=csv.QUOTE_NONE)

    introngff = makeGffFromGtf(intronsgtf)

    # intronsgtf = os.path.join(resourcedir, "test.GRCh37.b37.introns.gtf")
    # allintrons.to_csv( intronsgtf, sep="\t", index=False, header=False,
    # quoting=csv.QUOTE_NONE )
    # escapechar='\\', encoding='utf-8' )
# END MAIN
