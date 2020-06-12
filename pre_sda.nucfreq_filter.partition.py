#!/usr/bin/env python
import pysam
import argparse
import sys



ap = argparse.ArgumentParser(description="Find support from all combinations of nonref kmers")

ap.add_argument("--fraction", "-f", default=0.25,type=float)

ap.add_argument("-sc","--second_count", type=int)

ap.add_argument("--first_sub", "-fb",  required=True)


ap.add_argument("--second_sub", "-sb",  required=True)

ap.add_argument("--median", "-med", type=int, required=True)

ap.add_argument("--mean_cov", "-m", type=int, required=True)

ap.add_argument("--region", "-r",help="list of regions, chr:s-e",required=True,default=sys.stdin)

ap.add_argument("--prefout", help="prefix root folder")

ap.add_argument("--nucfreq", "-n",help="bamFreq intersect stream",required=True,default=sys.stdin)


#ap.add_argument("--lockfile", help="Make IO wait on the existence of this file", default=None)
ap.add_argument("--output", help="File to write to", default="filtered.orig_counts.nucfreq")

args = ap.parse_args()


reg=[]
with open(args.region,'r') as rg:
    for r in rg.readlines():
        reg.append(r)

def return_nucpos(nuc):
    if nuc=="A":
        return 0
    if nuc=="C":
        return 1
    if nuc=="G":
        return 2
    if nuc=="T":
        return 3


contigs = {}
ctglist = []
with open(args.nucfreq,'r') as n:
    for line in n.readlines():
        if line[0] == '#':
            continue
        vals=line.split()

        ctg=vals[0]
        
        if ctg not in contigs:
            contigs[ctg] = []
        
        if line in contigs[ctg]:
            raise ValueError('Found Duplicate: {}'.format(line))
        contigs[ctg].append(line)
        ctglist.append(vals[0])

chrset=set([m for m in ctglist])

for line in reg:
    line=line.rstrip()
    outfile="/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/"+line+"/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq"
    filtered="/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/"+line+"/kmer_filter/"+args.output
    
    line=line.split("_")[0]
    cont=line.split(":")
    contig= cont[0]

    if contig not in chrset:
        continue
    start=int(cont[1].split("-")[0])
    stop=int(cont[1].split("-")[1])
    #print(contig,start,stop)

    with open(filtered,'w') as filt, open(outfile,'w') as out:
        for liner in contigs[contig]:
            vals=liner.split()
            
            if vals[0]!=contig:
                raise ValueError('Dictionary error: {}{}'.format(vals[0],contig))
            pos=int(vals[1])

            if pos >= start and pos <= stop :
                info=vals[5].split(";")
                #print("entered")

                #print(liner,vals,info)

                mean=int(info[0].split("=")[1])
                median=int(info[1].split("=")[1])
                depth=info[2].split("=")[1]
                first_count=int(depth.split(",")[0])
                second_count=int(depth.split(",")[1])
        
               # RKN=info[4].split("=")[1]
                #AKN=info[6].split("=")[1]


                if mean < ((args.fraction*args.mean_cov)) or mean > ((args.fraction+1)*args.mean_cov):
                    continue
        
                if median < args.median:
                    continue

                if second_count < args.second_count:
                    continue
                
                orig_nuc_count=[None]*6
                nuc_count=[None]*6#vals[3:7]
                for i in range(4):
                    if int(return_nucpos(vals[3]))==i:
                        nuc_count[i]=args.first_sub
                        orig_nuc_count[i]=first_count
                    elif int(return_nucpos(vals[4]))==i:
                        nuc_count[i]=args.second_sub
                        orig_nuc_count[i]=second_count
                    else:
                        nuc_count[i]=0
                        orig_nuc_count[i]=0
                nuc_count[4]=0
                nuc_count[5]=0
                orig_nuc_count[4]=0
                orig_nuc_count[5]=0
                
                outs=vals[0:3]+nuc_count[:]
                orig=vals[0:3]+orig_nuc_count[:]


                filt.write("\t".join(orig) +"\n")
                out.write("\t".join(outs) +"\n" )
            
