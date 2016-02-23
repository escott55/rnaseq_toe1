#!/usr/bin/python

import os
import sys

import pandas as pd

sys.path.append(os.path.expanduser("~/workspace/custom_utils/general"))
import housekeeping as hk

############################
# Globals
############################
BASEDIR = os.path.abspath("../..")
RESOURCEDIR = os.path.join(BASEDIR, "resources")
FIGUREDIR = os.path.join(BASEDIR, "results")


def getGff(rtype, species, regionclass="coding"):
    if species == "Human":
        if regionclass == "coding" and rtype == "exon":
            gff = os.path.join(RESOURCEDIR, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.exons.gff")
        elif regionclass == "coding" and rtype == "intron":
            gff = os.path.join(RESOURCEDIR, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.introns.gff")
        elif regionclass == "noncoding" and rtype == "ncexon":
            gff = os.path.join(RESOURCEDIR, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.ncexons.gff")
        elif regionclass == "noncoding" and rtype == "ncintron":
            gff = os.path.join(RESOURCEDIR, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.ncintrons.gff")
        else:
            print "Error:", hk.whoami(), "Unknown rtype -", rtype
            sys.exit(1)

    elif species == "Zebrafish":
        if regionclass == "coding" and rtype == "exon":
            gff = os.path.join(RESOURCEDIR, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.exons.gff")
        elif regionclass == "coding" and rtype == "intron":
            gff = os.path.join(RESOURCEDIR, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.introns.gff")
        elif regionclass == "noncoding" and rtype == "ncexon":
            gff = os.path.join(RESOURCEDIR, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.ncexons.gff")
        elif regionclass == "noncoding" and rtype == "ncintron":
            gff = os.path.join(RESOURCEDIR, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.ncintrons.gff")
        else:
            print "Error:", hk.whoami(), "Unknown regionclass -", rtype
            sys.exit(1)

    else:
        print "Error:", hk.whoami(), ":unknown species specified"
        sys.exit(1)

    print species, regionclass, rtype
    print "GFF:", gff
    print "File exists:", os.path.exists(gff)
    assert os.path.exists(gff), "Couldn't find correct gff"+gff
    return(gff)
# END getGff


def getGtf(species, regionclass="coding"):
    if species == "Human":
        fullgtf = os.path.join(RESOURCEDIR, "human_b37_regions",
                               "Homo_sapiens.GRCh37.b37.genes.gtf.gz")
    elif species == "Zebrafish":
        fullgtf = os.path.join(RESOURCEDIR, "zebra_z10_regions",
                               "Danio_rerio.GRCz10.83.genes.gtf.gz")
    else:
        print "Error:", hk.whoami(), "Unknown species indicated -", species
        sys.exit(1)
    return fullgtf
# END getGtf


def oldGetGff(sampleset, resourcedir):
    # Identify gff files to use
    if sampleset in ["npc", "fibro", "newnpc"]:
        exongff = os.path.join(resourcedir,
                               "Homo_sapiens.GRCh37.b37.exons.gff")
        introngff = os.path.join(resourcedir,
                                 "Homo_sapiens.GRCh37.b37.introns.gff")

    if sampleset in ["npc_noncode", "fibro_noncode", "newnpc_noncode"]:
        exongff = os.path.join(resourcedir,
                               "Homo_sapiens.GRCh37.b37.ncexons.gff")
        introngff = os.path.join(resourcedir,
                                 "Homo_sapiens.GRCh37.b37.ncintrons.gff")

    if sampleset == "zebra_noncode":
        exongff = os.path.join(resourcedir,
                               "Danio_rerio.GRCz10.83.ncexons.gff")
        introngff = os.path.join(resourcedir,
                                 "Danio_rerio.GRCz10.83.ncintrons.gff")
    elif sampleset.startswith("zebra"):
        exongff = os.path.join(resourcedir,
                               "Danio_rerio.GRCz10.83.exons.gff")
        introngff = os.path.join(resourcedir,
                                 "Danio_rerio.GRCz10.83.introns.gff")

    # if sampleset.startswith("npc"):
        # sampletable = sampletable[ sampletable["samp"] != "npc1603A3" ]
        # Focus on the new data initially
        # sampletable = (sampletable[sampletable["Batch"] == 1]
        # .reset_index(drop=True))
    # if sampleset.startswith("newnpc"):
        # sampletable = (sampletable[sampletable["Batch"] == 2]
        # .reset_index(drop=True))

    return(exongff, introngff)
# END oldGetGff


def getSpecies(sampleset):
    if sampleset == "NewNPC":
        return "Human"
    elif sampleset == "OldNPC":
        return "Human"
    elif sampleset == "Fibroblast":
        return "Human"
    elif sampleset == "Zebrafish":
        return "Zebrafish"
    elif sampleset == "T283":
        return "Human"
    else:
        print "Error:", hk.whoami(), "Unknown sampleset", sampleset
# END getSpecies


def getDataDirectory(sampleset):
    # identify raw and target directories
    if sampleset == "NewNPC":
        tdir = os.path.join(BASEDIR, "raw/newnpc")
    elif sampleset == "Zebrafish":
        tdir = os.path.join(BASEDIR, "raw/zebrafish")
    elif sampleset == "Fibroblast":
        tdir = os.path.join(BASEDIR, "raw/fibroblast")
    elif sampleset == "OldNPC":
        tdir = os.path.join(BASEDIR, "raw/oldnpc")
    elif sampleset == "T283":
        # bamdir = os.path.join(BASEDIR, "raw/t283/rnaseq")
        tdir = os.path.join(BASEDIR, "raw/t283")
    else:
        print "Error:", hk.whoami(), "unknown sampleset provided"
        sys.exit(1)

    return tdir
# END getDataDirectory


def parseSampleTable(sampleset):
    rawdatadir = getDataDirectory(sampleset)
    sampletable = pd.read_csv(os.path.join(rawdatadir, "SampleTable.txt"),
                              sep="\t")
    sampletable.index = sampletable.samp
    return sampletable
# END parseSampleTable


def getRegionCountFiles(sampleset, rawdatadir, sampletable, gfftype):
    # Read in counts files
    if sampleset.count("noncode") > 0:
        exonfileslist = hk.filesInDir(os.path.join(rawdatadir, "ncexoncnt"),
                                      "txt")
        intronfileslist = hk.filesInDir(os.path.join(rawdatadir,
                                                     "ncintroncnt"), "txt")
    else:
        exonfileslist = hk.filesInDir(os.path.join(rawdatadir,
                                                   "exoncnt"), "txt")
        intronfileslist = hk.filesInDir(os.path.join(rawdatadir,
                                                     "introncnt"), "txt")

    # print exonfileslist
    # print intronfileslist

    prefix = [os.path.split(x)[1][:-9] for x in exonfileslist]
    # exonfilelist modify the file set to exclude erroneous samples
    exonfiles = pd.DataFrame({"ftype": sampleset, "cntfile": exonfileslist,
                              "prefix": prefix, "Regiontype": "exon",
                              "gfftype": gfftype})

    exonfiles = exonfiles[exonfiles.prefix.isin(sampletable.samp)]

    prefix = [os.path.split(x)[1][:-9] for x in intronfileslist]
    intronfiles = pd.DataFrame({"ftype": sampleset, "cntfile": intronfileslist,
                                "prefix": prefix, "Regiontype": "intron",
                                "gfftype": gfftype})
    intronfiles = intronfiles[intronfiles.prefix.isin(sampletable.samp)]

    print "Intron files:", len(intronfiles)
    print "Exon files:", len(exonfiles)

    assert len(exonfiles) > 0 and len(intronfiles) > 0
    return exonfiles, intronfiles
# END getRegionCountFiles


def makeFigureDir(sampleset, analysis=None):
    # Make Figure directory if necessary
    if analysis is not None:
        figuredir = os.path.join(FIGUREDIR, sampleset, analysis)
    else:
        figuredir = os.path.join(FIGUREDIR, sampleset)

    hk.makeDir(figuredir)
    return figuredir
# END makeFigureDir


if __name__ == "__main__":
    print "Exon-Human", getGff("exon", "Human")
    print "Intron-Human", getGff("intron", "Human")

# END
