#!/usr/bin/env python

import pysam
import argparse
import tempfile
import subprocess
import codecs
import sys
import numpy as np
import random
import time
import os

ap=argparse.ArgumentParser(description="Find support from all combinations of nonref kmers")


ap.add_argument("--vcf", help="Nucleotide frequency file.", required=True)
ap.add_argument("--k", help="Kmer size", default=22,type=int)
ap.add_argument("--jf", help="Jellyfish file.", required=True)
ap.add_argument("--gcbias", help="Optionally correct for gc bias.",default=None)
ap.add_argument("--jfstats", help="Jellyfish stats - mean, var, sd", required=True)
ap.add_argument("--sample", help="Sample being genotyped", required=True)
ap.add_argument("--ref", help="Reference sequence", required=True)
ap.add_argument("--lockfile", help="Make IO wait on the existence of this file", default=None)

args=ap.parse_args()
ref=pysam.FastaFile(args.ref)

gcRatio=None
gcPct=None
if args.gcbias is not None:

    gcBiasFile = open(args.gcbias)
    gcPct=[]
    cov=[]
    for line in gcBiasFile:
        vals=line.split()
        gcPct.append(float(vals[0]))
        cov.append(float(vals[1]))

    avgCov=np.average(cov)
    gcRatio = [c/avgCov for c in cov]

jfStatsFile = open(args.jfstats)
vals=jfStatsFile.readline().split()
jfMean = float(vals[0])
jfSD   = float(vals[2])

def GCCorrect(kmer, count, pct, ratios):
    if (len(kmer) == 0):
        return count
    kmer=kmer.upper()
    gc=kmer.count("G")+kmer.count("C")/float(len(kmer))
    # linear search is fine for now
    for i in range(0,len(pct)-1):
        if gc >= pct[i] and gc < pct[i+1]:
            return count / ratios[i]
    return count
    

vcf=open(args.vcf)
nQuery=0
k=args.k
query=""
sites=0
curChromName=""
curChrom=""
positions=[]
test=[]
import pdb

avgSupFile=open("avgSup.txt","w")


sys.stdout.write("""##fileformat=VCFv4.2
##fileDate=20190102
##source=genotype-vcf-by-kmers
##reference=/home/cmb-16/mjc/shared/references/hg38/hg38.fa
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##INFO=<ID=SA,Number=1,Type=String,Description="Names of original samples plus genotype">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=KS,Number=2,Type=Float,Description="Ref/Alt support by kmers">
##FORMAT=<ID=KN,Number=1,Type=Integer,Description="Number of kmers">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number of variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}
""".format(args.sample))

chroms         = []
positions      = []
refs           = []
alts           = []
altKmerIndices = []
smsSupport     = []
supportingSamples = []
for line in vcf:
    if len(line) > 0 and line[0] == "#":
        if line.split()[0] == "#CHROM":
            vals    = line.split()
            samples = vals[9:]
        continue
    vals    = line.split()
    pos     = int(vals[1])
    chrom   = vals[0]
    refNuc  = vals[3]
    altNuc  = vals[4]
    if len(altNuc) != 1:
        continue
    info    = vals[8].split(":")
    siIndex = -1
    spIndex = -1
    for vi in range(0,len(info)):
        if info[vi] == "SI":
            siIndex = vi
        if info[vi] == "SP":
            spIndex = vi
    # No values for genotyping
    if siIndex == -1 or spIndex==-1:
        continue
    allSI={}
    allSPList=[]
    #
    # Store supporting information for each sample 
    supSamples=[]
    sampleIndex=0
    for sample in vals[9:]:
        if sample[0] != ".":
            sampleData=sample.split(":")
            sampleSI = [ int(i) for i in sampleData[siIndex].split(",") ]
            
            for i in sampleSI:
                if i not in allSI:
                    allSI[i] = 0
                allSI[i] +=1
            allSPList.append(sampleData[spIndex])
            supSamples.append(samples[sampleIndex] + ","+sampleData[0])
        sampleIndex+=1
    supportingSamples.append("/".join(supSamples))
        
    chroms.append(chrom)
    positions.append(pos)
    alts.append(altNuc)
    refs.append(refNuc)
    altKmerIndices.append(sorted(allSI.keys()))
    if curChromName != chrom:
        sys.stderr.write("caching " + chrom + "\n")
        curChromName=chrom
        curChrom=ref.fetch(chrom)
    refSeq=curChrom[pos-args.k+1:pos+args.k].upper()

    altSeq=refSeq[0:k-1] + altNuc + refSeq[k:]

    if len(refSeq) != len(altSeq):
        print("ERROR!")
        sys.exit(1)
    query  += ">" + str(pos) + "_R" +"\n" + refSeq + "\n" + ">" + str(pos) + "_" + altNuc + "\n" + altSeq + "\n"

if args.lockfile is not None:
    while args.lockfile is not None and os.path.exists(args.lockfile):
        sleeptime=random.randint(10,20)
        sys.stderr.write("Sleeping for " + str(sleeptime) + "\n")
        time.sleep(sleeptime)

    of=open(args.lockfile, "w")
    of.close()

command = "jellyfish query --load --sequence=/dev/stdin {}".format(args.jf)
sys.stderr.write("querying jf \n")
proc        = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
proc_stdout = proc.communicate(input=bytes(query, 'utf-8'))
allLines    = codecs.latin_1_decode(proc_stdout[0])[0]
jfRes       = allLines.split("\n")
if (args.lockfile is not None):
    try:
        os.remove(args.lockfile)
    except:
        sys.stderr.write("oops probably had two processes in the same lock\n")


i=0
idx=0

while i < len(jfRes)-1:
    allSup={}
    for test in ["R", "A"]:
        sup=[]
        kmers=[]
        for j in range(0,k):
            if (i >= len(jfRes)):
                sys.stderr.write("ERROR at " + str(i) + " of " + str(len(jfRes))+"\n")
                sys.exit(1)
            vals=jfRes[i].split(" ")
            if len(vals) > 1:
                sup.append(int(vals[1]))
                kmers.append(vals[0])
            i+=1
        if gcPct is not None:            
            for s in range(0,len(kmers)):
                sup[s] = GCCorrect(kmers[s], sup[s], gcPct, gcRatio)
        allSup[test] = (sup,kmers)

    copyNumber=2
    #
    # Determine genotype only based on alt coverage.
    #
    if len(altKmerIndices[idx]) > 0:
        avgAltSupport = np.average(allSup["A"][0])
        genotype="None"
        avgRefSupport = np.average(allSup["R"][0])
        
        if avgAltSupport > jfMean + 5*jfSD:
            genotype="1/1"
            copyNumber=3
        elif avgAltSupport < max(jfMean/2 - 2*jfSD, 5):
            genotype="0/0"
            copyNumber=0
        elif abs(avgAltSupport - jfMean/2) < abs(avgAltSupport - jfMean):
            genotype="0/1"
            copyNumber=1
        else:
            genotype="1/1"
            copyNumber=2
        sys.stdout.write(chroms[idx] + "\t" + str(positions[idx]) + "\t.\t"  + refs[idx] + "\t" + alts[idx] + "\t30\tPASS\t" + "SA="+supportingSamples[idx] + "\tGT:KS:KN:CN\t" + "{}:{:2.2f},{:2.2f}:{}:{}".format(genotype, avgRefSupport,avgAltSupport,len(allSup["A"][0]), copyNumber) + "\n")

    idx+=1
sys.stderr.write(str(idx) + "\t" + str(len(positions)) + "\n")

    
    
