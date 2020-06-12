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

ap.add_argument("--region", "-r",description="list of regions chr:s-e",required=True)

ap.add_argument("--prefout", help="prefix root folder")


#ap.add_argument("--lockfile", help="Make IO wait on the existence of this file", default=None)
ap.add_argument("--output", help="File to write to", default="filtered.orig_counts.nucfreq")

args = ap.parse_args()



reg=open(args.region)
for r in reg:
    regi=reg.split("_")[0]
    
cont=regi.split(":")
contig= cont[0]
start=int(cont[1].split("-")[0])
stop=int(cont[1].split("-")[1])
print(contig,start,stop)

outfile="/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/" + args.region + "/kmer_filter/DUP_call."+args.region+".queries.counts"





filtered=args.output

with open(filtered,'w') as filt:
    for line in sys.stdin:
        if line[0] == '#':
            continue


        vals=line.split()
        info=vals[16].split(";")
        #print(info)
        mean=int(info[0].split("=")[1])
        median=int(info[1].split("=")[1])
        depth=info[2].split("=")[1]

        first_count=int(depth.split(",")[0])
        second_count=int(depth.split(",")[1])
        
        RKN=info[4].split("=")[1]
        AKN=info[6].split("=")[1]


        if mean < ((args.fraction*args.mean_cov)) or mean > ((args.fraction+1)*args.mean_cov):
            continue
        
        if median < args.median:
            continue

        if second_count < args.second_count:
            continue
        #print("passed")
        nuc_count=vals[3:7]
    #    sorted(vals[3:7])
        for i in range(4):
            if int(nuc_count[i])==first_count:
                vals[i+3]=args.first_sub
            elif int(nuc_count[i])==second_count:
                vals[i+3]=args.second_sub
            else:
                continue
        
        filt.write(line)
        sys.stdout.write("\t".join(vals[:9]) +"\n" )
    
