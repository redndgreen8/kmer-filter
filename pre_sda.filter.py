#!/usr/bin/env python
import pysam
import argparse
import sys



ap = argparse.ArgumentParser(description="Find support from all combinations of nonref kmers")

ap.add_argument("--first_sub", "-fb",  required=True)

ap.add_argument("--second_sub", "-sb",  required=True)

#ap.add_argument("--mean_cov", "-m", type=int, required=True)

ap.add_argument("--prefout", help="prefix root folder")

ap.add_argument("--nucfreq", "-n",help="bamFreq intersect stream",required=True,default=sys.stdin)

ap.add_argument("--region","-r",help="loci", required=True)

ap.add_argument("--output", help="File to write to", default="filtered.orig_counts.nucfreq")

args = ap.parse_args()




def return_nucpos(nuc):
    if nuc=="A":
        return 0
    if nuc=="C":
        return 1
    if nuc=="G":
        return 2
    if nuc=="T":
        return 3



filtered="/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/"+args.region+"/kmer_filter/"+args.output


with open(filtered,'w') as filt:

    for line in sys.stdin:
        line=line.rstrip()


        vals=line.split()

        ctg=vals[0]      
        start=vals[1]
        end=str(int(start)+1)
        loci = ctg + "\t" + start + "\t" + end


        info=vals[5].split(";")
        #print("entered")

        #print(liner,vals,info)

        mean=int(info[0].split("=")[1])
        median=int(info[1].split("=")[1])
        depth=info[2].split("=")[1]
        first_count=depth.split(",")[0]
        second_count=depth.split(",")[1]

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
                nuc_count[i]=str(0)
                orig_nuc_count[i]=str(0)
        nuc_count[4]=str(0)
        nuc_count[5]=str(0)
        orig_nuc_count[4]=str(0)
        orig_nuc_count[5]=str(0)
        
        outs="\t".join(nuc_count)
        orig="\t".join(orig_nuc_count)

        inf=vals[5]
        filt.write(loci +"\t"+orig +"\t"+inf+"\n")
        sys.stdout.write(loci +"\t"+outs +"\n" )
                
