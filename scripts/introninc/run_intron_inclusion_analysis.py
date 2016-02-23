#!/usr/bin/python

import os
import sys
import getopt
import math
import numpy as np
from scipy.stats import ttest_ind
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

import pandas as pd
import load_sampleinfo as sampinfo
import riiDistributions as riidist

sys.path.append(os.path.expanduser(
    "~/workspace/custom_utils/general"))
import housekeeping as hk
sys.path.append(os.path.expanduser(
    "~/workspace/custom_figures/rnaseq_intronret"))
import qcPlots

# R package loads
rstats = importr('stats', robject_translations={"format.perc": "format_perc2"})


def parseOtherIds(otherinfo):
    allinfo = {x[:x.find(" ")]: x[x.find(" ")+1:].replace('"', "")
               for x in [y.strip() for y in otherinfo.split(";")]}

    if "exonic_part_number" in allinfo:
        return allinfo["gene_id"]+":"+allinfo["exonic_part_number"]

    return allinfo["gene_id"]
# END parseOtherIds


def parseAttributes(attstr):
    newdict = {}

    for attpair in attstr.split(";"):
        attpair = attpair.strip()
        if attpair.count(" ") > 0:
            key, val = [x.strip().replace('"', '')
                        for x in attpair.split(" ", 1)]
            newdict[key] = val

    return pd.Series(newdict)
# END parseAttributes


def ensembl2gene(gff, fullgtf):
    # if gff.count("Homo") > 0:
        # fullgtf = os.path.join(sampinfo.RESOURCEDIR,
        # "Homo_sapiens.GRCh37.b37.genes.gtf.gz")

    # elif gff.count("Danio") > 0:
        # fullgtf = os.path.join(sampinfo.RESOURCEDIR,
        # "Danio_rerio.GRCz10.83.genes.gtf.gz")

    cols = ["chrom", "script", "rtype", "start", "end", "dbid",
            "strand", "unk", "attribute"]
    allgenes = pd.read_csv(fullgtf, sep="\t", header=None, names=cols,
                           compression="gzip")
    attdf = pd.DataFrame(allgenes.attribute.apply(parseAttributes))

    reflookup = {x["gene_id"]: (x["gene_name"]
                 if type(x["gene_name"]) is str else x["gene_id"])
                 for idx, x in attdf.iterrows()}
    return reflookup
# END ensembl2gene


def geneLookup(rid, reflookup):
    regionids, repid = rid.split(":")
    #  infodf = []

    if regionids.count("+") > 0:
        genes = []
        for geneid in regionids.split("+"):
            genename = reflookup[geneid] if geneid in reflookup else None
            if genename is not None:
                genes.append(genename)

        if len(genes) == 0:
            genes = [""]
        return {"geneid": geneid, "repid": repid, "uniqgene":
                ",".join(hk.uniqify(genes)), "genename": ",".join(genes)}

    else:
        geneid = regionids
        genename = reflookup[geneid] if geneid in reflookup else None
        return {"geneid": geneid, "repid": repid, "uniqgene": genename,
                "genename": genename}
# END geneLookup


def parseRegionFiles(filedf, gff, gtf):
    print "Running parseRegionFiles:", len(filedf), "files"

    alldat = None
    for idx, csamp in filedf.iterrows():
        cnts = pd.read_csv(csamp.cntfile, sep="\t", header=None,
                           names=[csamp.prefix], index_col=0)
        if alldat is None:
            alldat = cnts
        else:
            alldat[csamp.prefix] = cnts[csamp.prefix]

    alldat["ftype"] = csamp["ftype"]

    # Add GFF information
    # gffinfo = hk.parseGff(gff)

    # Replace parseGff code
    ############################################################
    assert os.path.exists(gff)
    cols = ["chrom", "script", "rtype", "start", "end", "dbid",
            "strand", "unk", "otherids"]
    gffdf = pd.read_csv(gff, sep="\t", header=None, names=cols)

    gffdf["regionid"] = gffdf.otherids.apply(parseOtherIds)
    tcols = ["chrom", "start", "end", "rtype", "strand", "regionid"]
    gffdfsub = (gffdf.loc[gffdf.rtype == "exonic_part", tcols]
                .reset_index(drop=True))
    # gffdfsub["length"] = gffdf["end"] - gffdf["start"]

    reflookup = ensembl2gene(gff, gtf)

    introninfo = pd.DataFrame([geneLookup(x, reflookup)
                              for x in gffdfsub.regionid])
    tcols = ["geneid", "uniqgene"]
    annotdf = pd.concat([gffdfsub, introninfo[tcols]], axis=1)
    ############################################################

    alldat = pd.merge(alldat, annotdf, left_index=True, right_on="regionid")
    print alldat.head()
    alldat = alldat.reset_index(drop=True)
    return alldat
    # introninfo = DataFrame([ parseRID(x, reflookup) for x in alldat.index])
    # alldat["geneid"] = [x[:9] for x in alldat.index]
    # annotdf = concat( [alldat, introninfo], axis=1 )
    # return annotdf
    # print alldat.head()
    # largedf = hk.splitByGene( annotdf )
    # annotdf["geneiduniq"] = [hk.uniqify(x.split(","))
    # if x.count(",") > 0 else x for x in annotdf.geneid]
    # largedf = hk.splitSingleColumn( targetdf, colname )
    # return largedf
# END parseRegionFiles


# FPKM - [# of mapped reads]/([length of transcript]/1000)/([total reads]/10^6)
def calcFPKM(targetdf, tcols, replace=False):
    totalreads = targetdf.loc[:, tcols].apply(sum)

    for sid, treads in totalreads.iteritems():
        fpkmvals = (targetdf[sid] / (targetdf["length"] / 1000) /
                    (treads / math.pow(10, 6)))
        if replace:
            targetdf.loc[:, sid] = fpkmvals
        else:
            targetdf[sid+"_fpkm"] = fpkmvals

    return targetdf
# END calcFPKM


def normIntronExpr(before, intron, after, tsamples):
    dist1 = intron.start - before.end
    dist2 = after.start - intron.end

    normvals = []
    if dist1 - dist2 > 0:
        # print "Before:",before
        # print "Intron:",intron
        # print "After:",after
        # print "Distances: Before:",dist1, "After:",dist2
        normvals = [None for x in tsamples]

    elif before.rclass != "exon" or after.rclass != "exon":
        # print "B:",before.rclass,"V:",intron.rclass,"A:",after.rclass
        normvals = [None for x in tsamples]

    else:
        for samp in tsamples:
            # print samp,"B:",before.loc[samp],
            # "V:",intron.loc[samp],"A:",after.loc[samp]
            maxexpr = max(before.loc[samp], after.loc[samp])
            if maxexpr > 0:
                normvals.append(min(intron.loc[samp] / float(maxexpr), 1.))
            else:
                normvals.append(0)

    # print normvals
    # print tsamples
    normvals = pd.Series({tsamples[i]: normvals[i]
                         for i in range(len(normvals))})
    tcols = ["chrom", "start", "end", "ftype", "regionid", "length",
             "geneid", "rclass", "uniqgene"]
    finalser = pd.concat([intron[tcols], normvals])
    finalser["dist"] = "%d:%d" % (dist1, dist2)
    finalser["classes"] = "%s:%s:%s" % (before.rclass, intron.rclass,
                                        after.rclass)
    # Should calculate the median expression values for each genotype here
    # finalser["maxFPKM_adj"] = maxexpr
    # print normvals
    # print intron[tcols]
    return finalser
# END normIntronExpr


def calcRII(row, aff, unaff):
    unaffmean = np.mean(row[unaff])
    affmean = np.mean(row[aff])
    if unaffmean == 0. or affmean == 0.:
        return 0.
    return math.log(np.mean(row[aff]) / unaffmean, 2)  # log10?
# END calcRII


def calcIntronRetention(fpkmdf, affset, unaffset):
    print "Running calcIntronRetention"

    tsamples = affset + unaffset

    # fpkmdf["median"] = fpkmdf[sampletable.samp].apply(np.median,axis=1)
    # fpkmdffilt = fpkmdf[ fpkmdf["median"] >= 1 ].reset_index(drop=True)
    # introndf = fpkmdf[fpkmdf["rclass"] == "intron"]
    fpkmdf["median_aff"] = fpkmdf[affset].apply(np.median, axis=1)
    fpkmdf["median_unaff"] = fpkmdf[unaffset].apply(np.median, axis=1)
    fpkmdffilt = fpkmdf[(fpkmdf["rclass"] == "exon") |
                        (fpkmdf["median_aff"] >= 2) |
                        (fpkmdf["median_unaff"] >= 2)].reset_index(drop=True)

    print fpkmdffilt.groupby("rclass").size()

    maxsamps = len(fpkmdffilt)
    limit = -1
    intronres = []
    for idx, intron in fpkmdffilt[fpkmdffilt.rclass == "intron"].iterrows():
        before = fpkmdffilt.ix[idx-1, ]
        if maxsamps <= idx+1:
            break
        after = fpkmdffilt.ix[idx+1, ]
        norm = normIntronExpr(before, intron, after, tsamples)
        intronres.append(norm)
        limit -= 1
        if limit == 0:
            break

    intronres = pd.DataFrame(intronres).reset_index(drop=True)

    assert len(intronres) > 0

    intronfilt = intronres[~(intronres.apply(lambda x:
                           any([np.isnan(x[samp])
                                for samp in tsamples]), axis=1))]

    intronfilt["median"] = intronfilt[tsamples].apply(np.median, axis=1)
    intronfilt = intronfilt[(intronfilt.length > 0) &
                            (intronfilt["median"] < 1.) &
                            (intronfilt["median"] > 0.) &
                            (intronfilt.dist == "1:1")]

    # aff = ["npc1603A1", "npc1603A3"]
    # unaff = ["npc1603U3", "npc1603U4", "npcControl"]
    intronfilt["logRII"] = intronfilt.apply(lambda row:
                                            calcRII(row, affset, unaffset),
                                            axis=1)
    intronfilt = intronfilt[intronfilt.logRII != 0]
    # print intronfilt.columns
    # print intronfilt[ intronfilt.length == 0][["uniqgene","length"]+
    # tsamples].head()
    # print intronfilt[ intronfilt.length.isnull() ].head()

    intronfilt["loglen"] = intronfilt.length.apply(lambda x: math.log10(x))

    intronfilt["pval"] = (intronfilt
                          .apply(lambda row: ttest_ind(row[affset],
                                                       row[unaffset])[1],
                                 axis=1))

    r_pvals = robjects.FloatVector(intronfilt.pval)
    intronfilt["bonf_p"] = (rstats.p_adjust(r_pvals, method='BH'))
    return intronfilt
# END calcIntronRetention


def exonAnalysis(sampletable, figuredir, alldf, aff, unaff, dtype):
    # Exon inclusion
    exondf = alldf[alldf.rclass == "exon"]
    exonfpkm = calcFPKM(exondf, sampletable.samp, replace=True)
    exonfpkm.loc["median_aff"] = exonfpkm[aff].apply(np.median, axis=1)
    exonfpkm.loc["median_unaff"] = exonfpkm[unaff].apply(np.median, axis=1)
    exonfilt = exonfpkm[(exonfpkm["median_aff"] >= 2) |
                        (exonfpkm["median_unaff"] > 2)].reset_index(drop=True)

    aff = sampletable.samp[sampletable.Status == "A"].tolist()
    unaff = sampletable.samp[sampletable.Status != "A"].tolist()

    # intronfilt["loglen"] = intronfilt.length.apply( lambda x: math.log10(x) )

    exonfilt["pval"] = exonfilt.apply(lambda row: ttest_ind(row[aff],
                                                            row[unaff])[1],
                                      axis=1)

    exonfilt["bonf_p"] = (rstats.p_adjust(robjects.FloatVector(exonfilt.pval),
                                          method='BY'))

    resfile = os.path.join(figuredir, gfftype+"_"+dtype+"_exonfinal.tsv")
    exonfilt.to_csv(resfile, sep="\t", index=False)

    # limit = 10
    # for geneid, genedf in exondf.groupby("uniqgene") :
    # print genedf.head()
    # exonusage = calcExonUsage( genedf, tsamples )
    # limit -= 1
    # if limit == 0 : break
# END exonAnalysis


def generateIntronRetentionValues(affset, unaffset, sampletable, speciesgtf,
                                  exonfiles, exongff, intronfiles,
                                  introngff, figuredir,
                                  gfftype, dtype, rerun=False):
    print "Running generateIntronRetentionValues"

    samplelist = affset + unaffset
    # sampletable.index = sampletable["samp"]

    print figuredir, gfftype, dtype
    intronretfile = os.path.join(figuredir,
                                 gfftype+"_"+dtype+"_final.tsv")

    if os.path.exists(intronretfile) and not rerun:
        print "Intronretention data already exists"
        return pd.read_csv(intronretfile, sep="\t")

    # tsamples = concat([npctsamples, fibtsamples]).reset_index()

    exondf = parseRegionFiles(exonfiles, exongff, speciesgtf)
    exondf["rclass"] = "exon"
    assert len(exondf) > 0

    introndf = parseRegionFiles(intronfiles, introngff, speciesgtf)
    introndf["rclass"] = "intron"
    assert len(introndf) > 0

    # Filter for at least 5 reads covering the intronic regions
    introndf = introndf[(introndf[samplelist].sum(axis=1) > 5)]
    # Zebrafish data - went from 265059 to 198451

    # Merge all counts into a single dataframe
    alldf = (pd.concat([exondf, introndf])
             .sort_values(by=["chrom", "uniqgene", "start"])
             .reset_index(drop=True))
    alldf["length"] = alldf.end-alldf.start

    # Filter if the intron if it has a 0 length
    alldf = alldf[alldf.length > 0]
    # Zebrafish data this was 2526 Exons, 502 introns

    print alldf.head()

    # Convert reads into fpkm values
    qcPlots.plotReadCountQc(alldf, gfftype, dtype, sampletable, figuredir)

    fpkmdf = calcFPKM(alldf, samplelist, replace=True)

    fpkmfile = os.path.join(figuredir, gfftype+"_"+dtype+"_fpkm.tsv")
    fpkmdf.to_csv(fpkmfile, sep="\t", index=False)
    qcPlots.plotFpkmQc(fpkmdf, gfftype, dtype, figuredir)

    intronret = calcIntronRetention(fpkmdf, affset, unaffset)

    intronret.sort(["bonf_p", "pval"], inplace=True)
    intronret.to_csv(intronretfile, sep="\t", index=False)

    qcPlots.plotIntronRetQc(intronret, sampletable, gfftype, dtype, figuredir)

    return intronret
# END generateIntronRetentionValues


def parseUserInput():
    # Currently not being used
    optlist, args = getopt.getopt(sys.argv[1:], "g:d:")
    optlist = dict(optlist)
    assert "-g" in optlist and "-d" in optlist, optlist
    gfftype = optlist.get("-g", None)
    dtype = optlist.get("-d", None)
    return gfftype, dtype
# END parseUserInput


def main(sampleset):
    print "Processing Sampleset:", sampleset

    regionclass = "coding"  # Alternatively - "noncoding"
    gfftype = "refseq"

    # Set variables
    figuredir = sampinfo.makeFigureDir(sampleset, analysis="intronret")

    # Example datasets - NewNPC, Zebrafish, T283
    species = sampinfo.getSpecies(sampleset)
    rawdatadir = sampinfo.getDataDirectory(sampleset)
    sampletable = sampinfo.parseSampleTable(sampleset)

    exongff = sampinfo.getGff("exon", species, regionclass)
    introngff = sampinfo.getGff("intron", species, regionclass)
    speciesgtf = sampinfo.getGtf(species, regionclass)

    exonfiles, intronfiles = (sampinfo.getRegionCountFiles(
        sampleset, rawdatadir, sampletable, gfftype))

    affset = sampletable.samp[sampletable.Status == "A"].tolist()
    unaffset = sampletable.samp[sampletable.Status != "A"].tolist()
    # unaffset = sampletable.samp[sampletable.Status == "C"].tolist()

    intronret = generateIntronRetentionValues(
        affset, unaffset, sampletable, speciesgtf, exonfiles, exongff,
        intronfiles, introngff, figuredir, gfftype, sampleset)

    # Write out introns that significantly deviate between case controls
    sigset = intronret[intronret.bonf_p < 0.05]
    sigsetfile = os.path.join(figuredir, gfftype+"_"+sampleset+"_sigset.tsv")
    sigset.to_csv(sigsetfile, sep="\t", index=False)

    qcPlots.plotIntronLenHist(intronret, gfftype, sampleset, figuredir)

    ####################################################################
    # Analyze RII distributions

    riidist.plotDensity(intronret, figuredir, gfftype, sampleset)

    if False:
        riidist.plotRowPermutations(intronret, affset, unaffset,
                                    figuredir, gfftype, sampleset)
    else:
        print "Skipping row permutations, it just takes too long"

    riidist.plotAffVersusUnaff(intronret, affset, unaffset,
                               figuredir, gfftype, sampleset)

    riidist.plotAffVersusAff(intronret, affset, unaffset,
                             figuredir, gfftype, sampleset)

    riidist.plotRiiByLen(intronret, gfftype, sampleset, figuredir)

    ####################################################################
    # irfilt = intronret[(intronret.pval < 0.05) ] #& bonf_p
    # (intronret.logRII.apply(np.abs) > 0) ]
# END MAIN


if __name__ == "__main__":
    # dtype = "npc"
    # dtype = "fibro"
    # dtype = "newnpc"
    # dtype = "zebra_AvU"
    # dtype = "zebra_UvC"
    # dtype = "newnpc_noncode"
    # dtype = "zebra_noncode"

    # NewNPC, Fibroblast
    # sampleset = "NewNPC"
    sampleset = "Fibroblast"
    sampleset = "Zebrafish"

    # main(sampleset)
# END MAIN
