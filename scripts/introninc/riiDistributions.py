#!/usr/bin/python

import os
import sys
# import getopt
import numpy as np
import math
from scipy.stats import gaussian_kde
import pandas as pd

sys.path.append(os.path.expanduser(
    "~/workspace/custom_utils/general"))
import housekeeping as hk
sys.path.append(os.path.expanduser(
    "~/workspace/custom_figures/ggplot_imports"))
from ggplotInclude import *


def plotDensity(intronret, figuredir, gfftype, dtype):
    r_dataframe = com.convert_to_r_dataframe(intronret)
    p = (ggplot2.ggplot(r_dataframe) +
         ggplot2.aes_string(x="logRII") +
         ggplot2.geom_density() +
         ggplot2.scale_y_continuous("Density") +
         ggplot2.scale_x_continuous("Log RII") +
         ggplot2.theme(**sitefreqtheme) +
         ggplot2.theme(**{'legend.position': "right"}))
    # ggplot2.facet_grid(robjects.Formula('status ~ .')) +
    # ggplot2.ggtitle("Loss of Function counts"+subtitle) +
    # ggplot2.stat_smooth(method="lm", se=False)
    # ggplot2.scale_x_continuous("External RVIS")+ \

    loffile = os.path.join(figuredir, gfftype+"_"+dtype+"_logirrdensity.pdf")
    print "Writing file %s" % loffile
    grdevices.pdf(loffile, width=7, height=6)
    p.plot()
    grdevices.dev_off()
# ENd plotDensity


# TODO
# This plot hasnt been made yet.
# Also, do we have to rerun the introncounts?
def plotRiiByIntronNumber(intronret, figuredir, gfftype, dtype):
    r_dataframe = com.convert_to_r_dataframe(intronret)
    p = (ggplot2.ggplot(r_dataframe) +
         ggplot2.aes_string(x="logRII") +
         ggplot2.geom_density() +
         ggplot2.scale_y_continuous("Density") +
         ggplot2.scale_x_continuous("Log RII") +
         ggplot2.theme(**sitefreqtheme) +
         ggplot2.theme(**{'legend.position': "right"}))
    # ggplot2.facet_grid(robjects.Formula('status ~ .')) +
    # ggplot2.ggtitle("Loss of Function counts"+subtitle) +
    # ggplot2.stat_smooth(method="lm", se=False)
    # ggplot2.scale_x_continuous("External RVIS")+ \

    loffile = os.path.join(figuredir, gfftype+"_"+dtype+"_logirrdensity.pdf")
    print "Writing file %s" % loffile
    grdevices.pdf(loffile, width=7, height=6)
    p.plot()
    grdevices.dev_off()
# END plotRiiByIntronNumber


def shuffle(df):
    cols = df.columns
    _ = [np.random.shuffle(row) for idx, row in df.iterrows()]
    return pd.DataFrame(df, columns=cols)
# END shuffle


def calcRII(row, aff, unaff):
    unaffmean = np.mean(row[unaff])
    affmean = np.mean(row[aff])
    if unaffmean == 0. or affmean == 0.:
        return 0.
    return math.log(np.mean(row[aff]) / unaffmean, 2)  # log10?
# END calcRII


def plotRowPermutations(intronret, affset, unaffset, figuredir, gfftype,
                        dtype, nboots=100):
    print "Running plotRowPermutations"
    inretmat = intronret.loc[:, affset+unaffset]

    inretmatcopy = inretmat.copy()

    # For N iterations,
    # Shuffle all rows of the intron retention Matrix
    # Calculate the RII values for each intron
    # Fit a density for the new values
    maxx = intronret["logRII"].max()
    samplings = []  # List of density distributions for each iteration
    ind = np.linspace(intronret["logRII"].min(), maxx, 512)
    for boot in range(0, nboots):
        print "Running boot:", boot
        shuffled = shuffle(inretmatcopy)
        # print shuffled.head()
        newriivals = shuffled.apply(lambda row: calcRII(row, affset, unaffset),
                                    axis=1)
        kdesub = gaussian_kde(newriivals)
        kdedf = pd.DataFrame({"subsetname": "Random%d" % boot,
                             "Density": kdesub.evaluate(ind), "logRII": ind})
        samplings.append(kdedf)

    samplings = pd.concat(samplings).reset_index(drop=True)

    # Calculate a 95% confidence interval using the
    # shuffled density distributions
    quants = (samplings.groupby(["logRII"])["Density"]
              .quantile([.025, .5, .975]).reset_index())
    quants.rename(columns={'level_1': "Quantile", 0: "Density"}, inplace=True)
    quants["linetype"] = ["Mean" if x == .5 else "95% threshold"
                          for x in quants.Quantile]

    # Calculate an Observed density using the original data
    newriivals = inretmat.apply(lambda row:
                                calcRII(row, affset, unaffset), axis=1)
    kde = gaussian_kde(newriivals)
    obsdf = pd.DataFrame({"vclass": "Observed", "Density": kde.evaluate(ind),
                          "logRII": ind})

    # Make the plot using the 3 forms of input data.
    # Sampling densities, 95% confidence interval regions,
    # and the observed values
    rsamplings = com.convert_to_r_dataframe(samplings)
    robsdf = com.convert_to_r_dataframe(obsdf)
    rquants = com.convert_to_r_dataframe(quants)
    rquants = fixRLevels(rquants, "linetype", ["Mean", "95% threshold"])
    p = (ggplot2.ggplot(robsdf) +
         ggplot2.aes_string(x="logRII", y="Density") +  # group="vclass"
         ggplot2.geom_line(ggplot2.aes_string(x="logRII", y="Density",
                                              group="factor(subsetname)"),
                           color="grey", data=rsamplings) +
         ggplot2.geom_line(ggplot2.aes_string(x="logRII", y="Density",
                                              linetype="factor(linetype)",
                                              group="factor(Quantile)"),
                           color="black", data=rquants) +
         ggplot2.geom_line(ggplot2.aes_string(color="factor(vclass)")) +
         ggplot2.scale_y_continuous("Density") +
         ggplot2.scale_x_continuous("Log RII") +
         ggplot2.scale_linetype("Confidence Interval") +
         ggplot2.scale_colour_brewer("Variant Type", palette="Set1") +
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle=45)}) +
         ggplot2.theme(**mytheme))

    figname = os.path.join(figuredir, gfftype+"_"+dtype+"_logirr_" +
                           str(nboots)+"boot.pdf")
    print "Writing file:", figname
    grdevices.pdf(figname, width=7, height=4.5)
    p.plot()
    grdevices.dev_off()
# END plotRowPermutations


def plotAffVersusUnaff(intronret, affset, unaffset, figuredir, gfftype, dtype):
    print "Running plotAffVersusUnaff"
    inretmat = intronret.loc[:, affset+unaffset]

    maxx = intronret["logRII"].max()
    AvUsamplings = []  # List of density distributions for each comparison
    ind = np.linspace(intronret["logRII"].min(), maxx, 512)

    for aff in affset:
        for unaff in unaffset:
            print "A:", aff, "versus", "U:", unaff
            newriivals = inretmat.apply(lambda row:
                                        calcRII(row, [aff], [unaff]), axis=1)
            kdesub = gaussian_kde(newriivals)
            kdedf = pd.DataFrame({"subsetname": aff+" vs "+unaff,
                                  "Affected": aff, "Unaffected": unaff,
                                  "Density": kdesub.evaluate(ind),
                                  "logRII": ind})
            AvUsamplings.append(kdedf)

    AvUsamplings = pd.concat(AvUsamplings).reset_index(drop=True)

    # Calculate an Observed density using the original data
    newriivals = inretmat.apply(lambda row:
                                calcRII(row, affset, unaffset), axis=1)
    kde = gaussian_kde(newriivals)
    obsdf = pd.DataFrame({"vclass": "Observed", "Density": kde.evaluate(ind),
                          "logRII": ind})

    rsamplings = com.convert_to_r_dataframe(AvUsamplings)
    # robsdf = com.convert_to_r_dataframe(obsdf)
    p = (ggplot2.ggplot(rsamplings) +
         ggplot2.aes_string(x="logRII", y="Density",
                            group="factor(subsetname)") +
         ggplot2.geom_vline(xintercept=0, linetype="dashed") +
         ggplot2.geom_hline(yintercept=0, linetype="solid") +
         ggplot2.geom_line(ggplot2.aes_string(color="factor(Unaffected)")) +
         ggplot2.scale_y_continuous("Density") +
         ggplot2.scale_x_continuous("Log RII") +
         ggplot2.scale_colour_brewer("Unaffected", palette="Set1") +
         ggplot2.facet_wrap(robjects.Formula('~ Affected'), ncol=3) +
         ggplot2.theme(**sitefreqtheme) +
         ggplot2.theme(**{'legend.position': "right"}))
    # ggplot2.geom_line( ggplot2.aes_string(x="logRII",y="Density",
    # group="factor(vclass)") +
    # linetype="dashed", color="black", data=robsdf ) +
    # ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +

    figname = os.path.join(figuredir, gfftype+"_"+dtype+"_logirr_AvU.pdf")
    print "Writing file:", figname
    grdevices.pdf(figname, width=10, height=8)
    p.plot()
    grdevices.dev_off()
# END plotAffVersusUnaff


def plotAffVersusAff(intronret, affset, unaffset, figuredir, gfftype, dtype):
    print "Running plotAffVersusAff"
    inretmat = intronret.loc[:, affset]

    maxx = intronret["logRII"].max()
    AvAsamplings = []  # List of density distributions for each comparison
    ind = np.linspace(intronret["logRII"].min(), maxx, 512)

    for i in range(len(affset)):
        for j in range(len(affset)):
            if i == j:
                continue
            print "A:", affset[i], "versus", "A:", affset[j]
            newriivals = inretmat.apply(lambda row:
                                        calcRII(row, [affset[i]], [affset[j]]),
                                        axis=1)
            kdesub = gaussian_kde(newriivals)
            kdedf = pd.DataFrame({"subsetname": affset[i]+" vs "+affset[j],
                                  "A1": affset[i], "A2": affset[j],
                                  "Density": kdesub.evaluate(ind),
                                  "logRII": ind})
            AvAsamplings.append(kdedf)

    AvAsamplings = pd.concat(AvAsamplings).reset_index(drop=True)

    # Calculate an Observed density using the original data
    # newriivals = inretmat.apply( lambda row: calcRII(row,affset,unaffset),
    # axis=1 )
    # kde = gaussian_kde( newriivals )
    # obsdf = pd.DataFrame({"vclass":"Observed", "Density":kde.evaluate(ind),
    # "logRII":ind})

    rsamplings = com.convert_to_r_dataframe(AvAsamplings)
    p = (ggplot2.ggplot(rsamplings) +
         ggplot2.aes_string(x="logRII", y="Density",
                            group="factor(subsetname)") +
         ggplot2.geom_vline(xintercept=0, linetype="dashed") +
         ggplot2.geom_hline(yintercept=0) +
         ggplot2.geom_line(ggplot2.aes_string(color="factor(A2)")) +
         ggplot2.scale_y_continuous("Density") +
         ggplot2.scale_x_continuous("Log RII") +
         ggplot2.scale_colour_brewer("Comparison", palette="Set1") +
         ggplot2.facet_wrap(robjects.Formula('~ A1'), ncol=3) +
         ggplot2.theme(**sitefreqtheme) +
         ggplot2.theme(**{'legend.position': "right"}))
    # robsdf = com.convert_to_r_dataframe(obsdf)
    # ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
    # ggplot2.geom_line( ggplot2.aes_string(x="logRII",y="Density") +
    # linetype="dashed", color="black", data=robsdf ) +

    figname = os.path.join(figuredir, gfftype+"_"+dtype+"_logirr_AvA.pdf")
    print "Writing file:", figname
    grdevices.pdf(figname, width=10, height=8)
    p.plot()
    grdevices.dev_off()
# END plotAffVersusAff


def plotRiiByLen(intronset, gfftype, dtype, figuredir):
    r_dataframe = com.convert_to_r_dataframe(intronset)
    p = (ggplot2.ggplot(r_dataframe) +
         ggplot2.aes_string(x="loglen", y="logRII") +
         ggplot2.geom_point(colour="#377EB8") +
         ggplot2.scale_y_continuous("Log RII ratio") +
         ggplot2.scale_x_continuous("Log intron length") +
         ggplot2.theme(**pointtheme))
    # ggplot2.ggtitle("Loss of Function counts"+subtitle) +
    # ggplot2.stat_smooth(method="lm", se=False)
    # ggplot2.scale_x_continuous("External RVIS")+ \

    loffile = os.path.join(figuredir, gfftype+"_"+dtype+"_logrii.pdf")
    print "Writing file %s" % loffile
    grdevices.pdf(loffile, width=7, height=6)
    p.plot()
    grdevices.dev_off()
# END plotRiiByLen


def runTest():
    gfftype = "refseq"
    dtype = "npcnew"
    rawdatadir = os.path.abspath("../../raw/newnpc")
    figuredir = os.path.abspath("../../results/intronret/test/")

    assert os.path.exists(rawdatadir)
    assert os.path.exists(figuredir)

    resfile = os.path.join(figuredir, gfftype+"_"+dtype+"_final.tsv")
    intronret = pd.read_csv(resfile, sep="\t")

    plotDensity(intronret, figuredir, gfftype, dtype)

    sampletable = pd.read_csv(os.path.join(rawdatadir, "SampleTable.txt"),
                              sep="\t")
    sampletable = sampletable[sampletable["Batch"] == 2]

    affset = sampletable.loc[sampletable.Status == "A", "samp"].tolist()
    unaffset = sampletable.loc[sampletable.Status == "U", "samp"].tolist()

    # plotRowPermutations( intronret, affset, unaffset,
    # figuredir, gfftype, dtype, nboots=10 )

    plotAffVersusUnaff(intronret, affset, unaffset, figuredir, gfftype, dtype)

    plotAffVersusAff(intronret, affset, unaffset, figuredir, gfftype, dtype)

    plotRiiByLen(intronret, gfftype, dtype, figuredir)
# END runTest


if __name__ == "__main__":
    runTest()
# END MAIN
