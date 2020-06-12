import argparse
import math
import statistics
import sys


ap = argparse.ArgumentParser(description="Parse the queries list/counts for support.")
ap.add_argument("--counts", help="File with counts of ref/alt sequence in HG38 and Illumina.", required=True,default=sys.stdin)
# ap.add_argument("--report", help="Output report file to write to", required=True)
ap.add_argument("--mincount", help="Minimum number of times a k-mer appears to be kept (class1/2).", default=5)
ap.add_argument("--maxcount", help="Maximum number of times a k-mer appears to be kept (class2).", default=34)
ap.add_argument("--minill", help="Minimum number of times a k-mer appears to be kept for Illumina (class1/2).", default=5)
ap.add_argument("--maxill", help="Maximum number of times a k-mer appears to be kept for Illumina (class1/2).", default=120)
ap.add_argument("--mean_cov","-m", help="Mean sequencing coverage", required=True)
ap.add_argument("--f", help="filter 0/1",default=0)

ap.add_argument("--fraction", help="fraction", type=float,default=0.25)
ap.add_argument("--vcf", help="vcf header. To be generated for each genome", required=True)
ap.add_argument("--k", help="Kmer size", type=int,default=32)


ap.add_argument("--median", "-med", type=int, default=10)

ap.add_argument("-sc","--second_count", type=int,default=6)

args=ap.parse_args()

counts_filepath = args.counts
min_counts = int(args.mincount)
max_counts = int(args.maxcount)
min_ill_given = int(args.minill)
max_ill_given = int(args.maxill)

mean = float(args.mean_cov)


def _calculate_support(alt_list,min_ill,max_ill):
    return min_ill <= sum([1 for count in alt_list if count > 0]) <= max_ill
    #change max_il

c0=[]

with open(counts_filepath) as f:
    for line in f:
        contig,pos,ref_seq, alt_seq, ref_hg_list, alt_hg_list, ref_ill_list,  alt_ill_list, DP = line.split()
        
        alt_ill_list = [int(v) for v in alt_ill_list.split(',')]
        len_alt_list=len(alt_ill_list)      
   
#//edits
        altsum=0
        for i in range(len(alt_ill_list )):
            altsum=altsum+alt_ill_list[i]
        mean_alt_support=altsum/len_alt_list

        med2=statistics.median(alt_ill_list)
        
#edits//
        alt_hg_list = [int(v) for v in alt_hg_list.split(',')]
        alt_hg_sup = [1 for count in alt_hg_list if count == 0]


        line=line.rstrip()+"\t"+str(mean_alt_support)+"\t"+str(med2)+"\n"
        
        #class0/all calls logic
        alt_hg_support = sum(alt_hg_sup) >= min_counts
       
        if _calculate_support(alt_ill_list, min_ill_given, max_ill_given) and alt_hg_support:
            ref_hg_list = [int(v) for v in ref_hg_list.split(',')]
            ref_hg_sup1 = [1 for count in ref_hg_list if count == 1]
            ref_hg_sup2 = [1 for count in ref_hg_list if 1 < count ]# max_counts]
            
            #class 1
            cl=0
            if sum(ref_hg_sup1) >= min_counts:
                cl=str(1)

            #class 2
            elif sum(ref_hg_sup2) >= min_counts:
                cl=str(2)    


            if mean_alt_support < ((0.1*float(args.mean_cov) )) or mean_alt_support > ((args.fraction+1)*float(args.mean_cov) ):
                continue
            
            if med2 < args.median:
                continue

            if int(DP.split(",")[1]) < args.second_count:
                continue

            c0.append(line.rstrip()+"\t"+cl+"\n")


print(str(len(c0)))

VCFheader = open(args.vcf)

outlines = []
for line in VCFheader:
    if line[0] == '#':
        outlines.append(line.rstrip())
VCFheader.close()


svIndex=-1
#0      1     2       3       4             5             6             7               8   9       10  11
#contig,pos,ref_seq, alt_seq, ref_hg_list, alt_hg_list,ref_ill_list,  alt_ill_list , mean , median , DP,class
b=0
for line in c0:
    line=line.rstrip()
    contig,pos,ref_seq, alt_seq, ref_hg_list, alt_hg_list, ref_ill_list,  alt_ill_list,DP,meann,mediann, clas = line.split()
    if b<2:
        print(line)
        b=b+1

    outp=[]
    outp.append(contig)
    outp.append(pos)
    outp.append(".")
    outp.append(ref_seq[args.k-1])
    outp.append(alt_seq[args.k-1])
    outp.append(".")
    outp.append(".")



    svIndex+=1
    ref_k_vals = []
    alt_k_vals = []
    
    mean_alt=0
    med_alt=0

    if int(clas)==1:
        #counts returned
        ref_kmers = ref_hg_list
        alt_kmers = alt_hg_list

        mean_alt=meann
        med_alt=mediann

        #get the idx of the kmers for c1 (ref, alt) that contributed to classification
        #for c1, ref count should be exactly 1
        #        alt count should be exactly 0
        for idx,count in enumerate(ref_kmers):
            if count == 1:
                ref_k_vals.append(idx)
        for idx,count in enumerate(alt_kmers):
            if count == 0:
                alt_k_vals.append(idx)
                
        #get the idx of the kmers for c2 (alt only) that contributed to classification
    elif int(clas)==2:
        mean_alt=meann
        med_alt=mediann
    

    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    #chr1    16576388        .       G       G       .       .       mean=38;med=0;DP=205,65.;REFKMERIDX=;RKN=22;ALTKMERIDX=;AKN=22

    med_alt=round(float(med_alt))
    mean_alt=round(float(mean_alt))


    info="mean="+str(mean_alt)+";"+"med="+str(med_alt)+";"+"DP="+DP+";"+"class="+clas+";ref_alt_kmer="+ref_seq+","+alt_seq
    
  #  if len(ref_k_vals) > 0 or len(alt_k_vals) > 0:
   #     if len(ref_k_vals) > 0:
    info = info + ';{}{}'.format('REFKMERIDX=',','.join([str(val) for val in ref_k_vals]))
    info = info + ';{}{}'.format('RKN=',len(ref_k_vals))
    if len(alt_k_vals) > 0:
        info = info + ';{}{}'.format('ALTKMERIDX=',','.join([str(val) for val in alt_k_vals]))
        info = info + ';{}{}'.format('AKN=',len(alt_k_vals))

    
    outp.append(info)
    outlines.append('\t'.join(outp))


print(args.counts + '.filtered.vcf')

with open(args.counts + '.filtered.vcf','w') as f:
    f.write( '\n'.join(outlines))
#print(outlines)


