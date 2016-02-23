#!/usr/bin/python

# Given a bamfile and a region GFF, make a count file for the number of
# fragments overlapping all identified regions
import os
import sys
import subprocess
from joblib import Parallel, delayed

import pandas as pd
import load_sampleinfo as sampinfo

sys.path.append(os.path.expanduser("~/workspace/custom_utils/general"))
import housekeeping as hk


def testPipedOutput(currP):
    headcmd = ["head", "-30"]

    test = subprocess.Popen(headcmd, stdin=currP.stdout,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    currP.stdout.close()

    stdout, stderr = test.communicate()
    print "STDOUT:"
    print stdout
    print "STDERR:"
    print stderr
    sys.exit(1)
# END testPipedOutput


def CleanupBamFile(infile, prefix, targetdir, force=False):
    # Clean up Shawshanks bam files so they work with htseq
    # remove the 18th and 19 column and missing pairs
    #
    # Example command:
    # samtools view -h \
    # rnaseq/WTCaf_710_ACAGTG_L002_R1.polyATrim.adapterTrim.rmRep.sorted.bam |\
    # cut -f18,19 --complement | \
    # awk '$7 != "*"' - | \
    # samtools view -Shu - -o WTCaf_710_ACAGTG_L002_R1.clean.bam

    cleanbam = os.path.join(prefix, targetdir+".clean.bam")
    sortbam = os.path.join(prefix, targetdir+".sort.bam")

    # samtools view -h GleesonKO_25_S7_L006_R1_001.bam |
    # cut -f18,19 --complement | awk '$7 != "*"' - |
    # samtools view -Shu - -o GleesonKO_25_S7_L006_R1_001.clean.bam
    # if not os.path.exists(cleanbam) or force :
    if not os.path.exists(sortbam) or force:
        print "Making Clean Bam:", cleanbam
        first = ["samtools", "view", "-h", infile]
        second = ['cut', '-f18,19', '--complement']
        third = ['awk', '$7 != "*"', "-"]
        fourth = ['samtools', 'view', '-Shu', '-', '-o', cleanbam]

        p1 = subprocess.Popen(first, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        p2 = subprocess.Popen(second, stdin=p1.stdout, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(third, stdin=p2.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2.stdout.close()
        p4 = subprocess.Popen(fourth, stdin=p3.stdout, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        p3.stdout.close()

        # testPipedOutput(p4)
        stdout, stderr = p4.communicate()
        if stdout is not None or stderr is not None:
            print "STDOUT:\n", stdout
            print "STDERR:\n", stderr

        assert os.path.exists(cleanbam)

        # samtools sort -n -m8G \
        # GleesonKO_25_S7_L006_R1_001.bam GleesonKO_25_S7_L006_R1_001.sort
        print "Making Sort Bam:", sortbam[:-4]
        sortcmd = ["samtools", "sort", "-n", "-m8G", cleanbam, sortbam[:-4]]
        p1 = subprocess.Popen(sortcmd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        if stdout is not None or stderr is not None:
            print "STDOUT:\n", stdout
            print "STDERR:\n", stderr
            # sys.exit(1)

    assert os.path.exists(sortbam)
    print "Clean bam:", sortbam
    return sortbam
# END CleanupBamFile


def sortBamFile(infile, prefix, targetdir, force=False):
    sortbam = os.path.join(prefix, targetdir+".sort.bam")

    if not os.path.exists(sortbam) or force:
        print "Making Clean Bam:", sortbam

        print "Making Sort Bam:", sortbam[:-4]
        sortcmd = ["samtools", "sort", "-n", "-m8G", infile, sortbam[:-4]]
        print " ".join(sortcmd)
        p1 = subprocess.Popen(sortcmd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        stdout, stderr = p1.communicate()
        if stdout is not None or stderr is not None:
            print "STDOUT:\n", stdout
            print "STDERR:\n", stderr
            # sys.exit(1)

    assert os.path.exists(sortbam)
    print "Clean bam:", sortbam
    return sortbam
# END sortBamFile


def runHtseqTest(cleanbam, gff, ftype="bam"):
    if ftype == "bam":
        command = ["htseq-count", "-f", "bam", "-s", "no", cleanbam, gff]
        p1 = subprocess.Popen(command, stdout=subprocess.PIPE)
        # stderr=subprocess.PIPE )
        stdout, stderr = p1.communicate()

    elif ftype == "sam":
        p1 = subprocess.Popen(["samtools", "view", cleanbam],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        command = ["htseq-count", "-f", "bam", "-s", "no", cleanbam, gff]
        p2 = subprocess.Popen(command, stdin=p1.stdout)
        p1.stdout.close()
        stdout, stderr = p2.communicate()

    print stdout
# END runHtseqTest


def runDexSeqCount(cleanbam, gff, countfiledir, prefix, force=False):
    print "runDexSeqCount:", cleanbam
    dexoutfile = os.path.join(countfiledir, prefix+".txt")

    if not os.path.exists(dexoutfile) or force:
        # python /usr/local/lib/R/site-library/DEXSeq/python_scripts/
        # dexseq_count.py -f bam -p no -s no -a 10
        # /home/escott/workspace/mboat7/resources/gencode.vM7.annotation_clip.gff
        # GleesonKO_25_S7_L006_R1_001.bam
        # GleesonKO_25_S7_L006_R1_001.exoncnts2.txt
        dexcount = os.path.join("/usr/local/lib/R/site-library/DEXSeq",
                                "python_scripts/dexseq_count.py")

        dexcmd = ["python", dexcount, "-f", "bam", "-p", "no", "-s", "no",
                  "-a", "10", gff, cleanbam, dexoutfile]

        print "Dexcommand :", dexcmd

        p1 = subprocess.Popen(dexcmd, stdout=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        if stdout is not None or stderr is not None:
            print "STDOUT:\n", stdout
            print "STDERR:\n", stderr

    return dexoutfile
# END runDexSeqCount


def runSingleTest():
    # originalbam =
    # "raw/npc/ATCC_ATCACG_L001_R1.polyATrim.adapterTrim.rmRep.sorted.bam"
    originalbam = "raw/test.bam"
    gff = "resources/Homo_sapiens.GRCh37.hg19.gff"

    cleanbam = CleanupBamFile(originalbam, "clean", "test", False)

    # htout = runHtseqTest(cleanbam, gff, ftype="bam")

    dexoutfile = runDexSeqCount(cleanbam, gff, "clean", "test")
    return(dexoutfile)
# END runSingleTest


def processBam(cfile, gff, countfiledir):
    cpath, cbase = os.path.split(cfile)
    prefix = cbase[:cbase.find(".")]
    if not os.path.exists(os.path.join(countfiledir, prefix+"_cnts.txt")):
        print "Processing:", cfile, "Prefix:", prefix
        # cleanbam = sortBamFile(cfile, countfiledir, prefix, False )
        dexoutfile = runDexSeqCount(cfile, gff, countfiledir, prefix+"_cnts",
                                    False)
        print("Finished:", dexoutfile)
    else:
        print "File found:", cfile, "Prefix:", prefix
# END processBam


def cleanAllBamFiles(bamfiles, targetdir):
    for cfile in bamfiles:
        cpath, cbase = os.path.split(cfile)
        prefix = cbase[:cbase.find(".")]
        if not os.path.exists(os.path.join(targetdir, prefix+".sort.bam")):
            print "Processing:", cfile, "Prefix:", prefix
            # cleanbam = CleanupBamFile(cfile, targetdir, prefix, True)
            # cleanbam = sortBamFile(cfile, targetdir, prefix, False)
        else:
            print "File found:", cfile, "Prefix:", prefix
# END cleanAllBamFiles


def runTest():
    # To do: Fix this test
    cfile = os.path.join("raw/npc",
                         "1603_A3_TTAGGC_L001_R1.polyATrim."
                         "adapterTrim.rmRep.sorted.bam")
    cpath, cbase = os.path.split(cfile)
    # subname = cbase[:cbase.find(".")]
    # prefix = npcnames[subname]
    # print "Processing:",cfile,"Prefix:",prefix
    # cleanbam = CleanupBamFile( cfile, targetdir, prefix, True )
    # dexoutfile = runDexSeqCount( cleanbam, gff, targetdir,
    # prefix+"_"+gfftype, True )
    # gff = os.path.join(basedir,"resources",
    # "gencode.vM7.annotation_exonsonly.gff")
    # ("gencode.vM7.annotation_clip.gff")
    # python /usr/local/lib/R/site-library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
    # gencode.vM7.annotation_exonsonly.gtf gencode.vM7.annotation_exonsonly.gff
# END runTest


def main(sampleset):
    print "Sampleset:", sampleset
    species = sampinfo.getSpecies(sampleset)
    tdir = sampinfo.getDataDirectory(sampleset)

    regiontypes = ["exon", "intron"]  # ["ncintron", "ncexon"]
    sampletable = pd.read_csv(os.path.join(tdir, "SampleTable.txt"), sep="\t")
    print sampletable.head()
    # bamfiles = hk.filesInDir(bamdir, "bam") # alternative method
    bamfiles = [os.path.join(tdir, x) for x in sampletable.File.tolist()]
    assert len(bamfiles) > 0, "Error: No bamfiles provided?"

    print "Processing sampleset:", sampleset
    print "Using target directory:", tdir
    # print "Target bamfiles:", bamfiles
    missingbams = [x for x in bamfiles if not os.path.exists(x)]
    assert (len(missingbams) == 0), (
            "Error: not all bamfiles exist "+"; ".join(missingbams))

    # cleanAllBamFiles(bamfiles, tdir)

    # For each regionclass, generate region counts
    for rtype in regiontypes:
        countfiledir = os.path.join(tdir, rtype+"cnt")
        hk.makeDir(countfiledir)
        gff = sampinfo.getGff(rtype, species)
        print "Using GFF:", gff

        # Single File test
        # processBam(bamfiles[0], gff, countfiledir)
        # Run dexseq across all bam files
        Parallel(n_jobs=3, backend="threading")(delayed(processBam)
                                                (bfile, gff, countfiledir)
                                                for bfile in bamfiles)
# END main


if __name__ == "__main__":
    # Example datasets - NewNPC, Zebrafish, T283
    # sampleset = "NewNPC"
    # sampleset = "Fibroblast"
    # sampleset = "Zebrafish"
    # sampleset = "OldNPC"
    sampleset = "T283"

    main(sampleset)
# END MAIN
