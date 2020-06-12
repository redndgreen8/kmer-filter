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
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Find support from all combinations of nonref kmers")

ap.add_argument("--bpseqs", help="Breakpoint sequence file, two sequences (ref,alt) per variant.", required=True, default=sys.stdin)
ap.add_argument("--k", help="Kmer size", type=int,default=32)
ap.add_argument("--jf", help="Jellyfish file of matched Illumina data (e.g./home/cmb-16/mjc/hgsvg/datasets/HG00514/Illumina/mer_counts.q0.L4.jf) ", required=True)

#ap.add_argument("--hgjf", help="Human genome jellyfish file. (e.g. /home/cmb-16/mjc/shared/references/hg38_noalts/indexes/jf/hg38.m22.jf). ", required=True)
ap.add_argument("--lockfile", help="Make IO wait on the existence of this file", default=None)
ap.add_argument("--output", help="File to write to", required=True)

args = ap.parse_args()

bpSeqFile = open(args.bpseqs)
k = args.k
jf = args.jf
#hgjf = args.hgjf
lockfile = args.lockfile
outfile = args.output
queries=[]

##
## Modify this loop to add reference/query sequences (alternating every other entry)
nEntries = 0
queryLengths = []
queryName = []

for idx,record in enumerate(SeqIO.parse(bpSeqFile, "fasta")):
   # print(str(record.name))
    queries.append(">{}\n".format(record.name) + str(record.seq))               #the "name" is already formatted with information
    nEntries +=1
    queryLengths.append(len(record.seq))
    queryName.append(">{}".format(record.name))

#    sys.stdout.write(str(idx) +" "+ str(record.seq) +" "+str(len(record.seq)) +"\n")
#time.sleep(30)

# added one too many new lines, this causes jellyfish to freak out.

if lockfile is not None:
    while lockfile is not None and os.path.exists(lockfile):
        sleeptime=random.randint(10,30)
        sys.stderr.write("Sleeping for " + str(sleeptime) + "\n")
        time.sleep(sleeptime)
    of=open(lockfile, "w")
    of.close()

sys.stderr.write("making query\n")

query = '\n'.join(queries).rstrip() + "\n"

command = "jellyfish query --load --sequence=/dev/stdin {}".format(jf)

sys.stderr.write("Submitting query " + str(len(queries)) + "\n")
proc        = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
proc_stdout = proc.communicate(input=bytes(query, 'utf-8'))
allLines    = codecs.latin_1_decode(proc_stdout[0])[0]
jfRes       = allLines.split("\n")

#command = "jellyfish query --load --sequence=/dev/stdin {}".format(hgjf)

#sys.stderr.write("querying hg jf \n")
#proc        = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
#proc_stdout = proc.communicate(input=bytes(query, 'utf-8'))
#allLines    = codecs.latin_1_decode(proc_stdout[0])[0]
#jfRefRes    = allLines.split("\n")
print("il has " + str(len(jfRes)))
if jfRes[-1] == '':
    del jfRes[-1] #delete empty value
#if jfRefRes[-1] == '':
 #   del jfRefRes[-1] #delete empty value

exp_result_len = sum([l-(k-1) for l in queryLengths])

if (lockfile is not None):
    try:
        os.remove(lockfile)
    except:
        sys.stderr.write("oops probably had two processes in the same lock\n")

#
# Parse jellyfish output
# The total k-mers sent to jellyfish can be computed from sum(queryLengths) - nEntries * (k-1)
#
#
if any([len(jfRes) != exp_result_len,               #make sure that there are the correct number of JF results
        len(queries) % 2 != 0,                      #make sure there are 2 queries for each variant
        ]):
    print("one of these is not like the other " + str(len(jfRes)) + " "  + str(exp_result_len) + "  "  + str(len(queries))        )
    sys.exit(1)

#
# Iterate over the queries 2 at a time (ref and alt for each variant)
results_idx = 0
with open(outfile,'w') as f:
    for idx in range(int(len(queries)/2)):
        ref_idx = 2*idx
        alt_idx = 2*idx + 1
        
        #both sequences include comment line so we split them
        ref_id,ref_seq = queries[ref_idx].split('\n')
        alt_id,alt_seq = queries[alt_idx].split('\n')
        
        assert str(queryName[ref_idx])==str(ref_id), "ref id not matching queryName{}{}".format(queryName[ref_idx],ref_id)
        
        assert str(queryName[alt_idx])==str(alt_id),"alt id not matching queryName{}{}".format(queryName[alt_idx],alt_id)

        #grabbing ref results
        num_ref_results = len(ref_seq) - (k-1)

        #ref_hg_counts = [result.split()[1] for result in jfRefRes[results_idx:results_idx+num_ref_results]]
        ref_il_counts = [result.split()[1] for result in jfRes[results_idx:results_idx+num_ref_results]]
        results_idx += num_ref_results
        
        #grabbing alt results
        num_alt_results = len(alt_seq) - (k-1)
       # alt_hg_counts = [result.split()[1] for result in jfRefRes[results_idx:results_idx+num_alt_results]]
        alt_il_counts = [result.split()[1] for result in jfRes[results_idx:results_idx+num_alt_results]]
        results_idx += num_alt_results
        
        #write to file
       # ref_hg_str = ','.join(ref_hg_counts)
        ref_il_str = ','.join(ref_il_counts)
     #   alt_hg_str = ','.join(alt_hg_counts)
        alt_il_str = ','.join(alt_il_counts)
        
        ##>chr1/16576202/alt/T/0/DP
        assert str(queryName[ref_idx].split("/")[3] )==ref_seq[k-1], "ref allele is not matching ref_query allele"
        assert str(queryName[alt_idx].split("/")[3] )==alt_seq[k-1], "alt allele is not matching alt_query allele"
        q=queryName[ref_idx].replace(">","")
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    #chr1    16576202        .       A       T       275,9
        contig,pos=q.split("/")[0:2]
        DP=q.split("/")[-1]
        #Output 7 tab-separated columns
        #0- id
        #1- ref sequence
        #2- alt sequence
        #3- ref seq count in human genome
        #4- ref seq count in Illumina
        #5- alt seq count in human genome
        #6- alt seq count in Illumina
        out_string = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig,pos,ref_seq,alt_seq,ref_il_str,alt_il_str,DP)
        f.write(out_string)