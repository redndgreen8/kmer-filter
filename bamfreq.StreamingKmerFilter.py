#!/usr/bin/env python
import pysam
import argparse
import tempfile
import subprocess
import codecs
import sys
import os
import random
import time

ap=argparse.ArgumentParser(description="Find support from all combinations of nonref kmers")

ap.add_argument("--ref", help="Reference sequence", required=True)
ap.add_argument("--nucfreq", help="Nucleotide frequency file.", required=True)
ap.add_argument("--k", help="Kmer size", default=32,type=int)
ap.add_argument("--f", help="alt count threshold", default=4,type=int)


args=ap.parse_args()
ref=pysam.FastaFile(args.ref)
if args.nucfreq == "stdin" or args.nucfreq == "-" or args.nucfreq == "/dev/stdin":
    nucfreq=sys.stdin
else:
    nucfreq=open(args.nucfreq)
nQuery=0
k=args.k
chroms=[]
positions=[]
test=[]
query=""
sites=0
curChromName=""
curChrom=""

queries=[]
nLines=0
faiFile=open(args.ref+ ".fai")
inFai = { line.strip().split()[0] : True for line in faiFile }


sys.stdout.write("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.9+htslib-1.9
##bcftoolsCommand=mpileup -O v -I -f /home/cmb-16/mjc/shared/references/hg38/hg38.fa -q 30 -R DUPcalls.masked_CN.AK1.clr.lra.POSTSCALE.viterout.bed -o DUPcalls.masked_CN.AK1.clr.lra.vcf /staging/mjc/jingwenr/AK1/AK1.lra.bam
##reference=file:///home/cmb-16/mjc/shared/references/hg38/hg38.fa
##contig=<ID=chr1,length=248956422>
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
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT\n""")


for line in nucfreq:
    vals=line.split()
    chrom=vals[0]
    
    if "_" in chrom or chrom not in inFai:
        continue
    pos=int(vals[1])
    counts = [(int(vals[2]), 'A'), (int(vals[3]),'C'), (int(vals[4]),'G'), (int(vals[5]),'T')]
    


    nDel = int(vals[6])


    nLines+=1
    if nLines % 1000000 == 0:
        sys.stderr.write("streamed " + str(nLines) + " lines for filter and stored " + str(len(queries)) + "\n")
    counts.sort()
    if curChromName != chrom:
        sys.stderr.write("caching " + chrom + "\n")
        curChromName=chrom
        try:
            curChrom=ref.fetch(chrom)
        except KeyError:
            continue
            

    total=counts[0][0]+counts[1][0]+counts[2][0]+counts[3][0]#+nDel
    altMult=0
    refMult=0
    refSeq=curChrom[pos-args.k+1:pos+args.k].upper()
    refNuc=refSeq[k-1].upper()
    #print( str(len(refSeq)) +" "+refNuc +" "+refSeq+" refFrominput:"+ ref.fetch(chrom,pos,pos+1 ) )
    
    #
    # Check to see if the nucleotide counts on a line support a variant.
    #

    if counts[-2][1] != refNuc:
        altNuc  = counts[-2][1]
        altMult = counts[-2][0]
        refMult = counts[-1][0]        
    else:
        altNuc = counts[-1][1]
        altMult = counts[-1][0]
        refMult = counts[-2][0]
    

    if altMult < args.f:
        continue
    else:
        sys.stdout.write(chrom+"\t"+str(pos)+"\t"+"."+"\t"+refNuc+"\t"+altNuc+"\t"+str(refMult)+","+str(altMult) +"\n")


