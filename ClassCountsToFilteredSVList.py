#!/usr/bin/env python


import argparse
import pysam
from collections import Counter

ap = argparse.ArgumentParser(description="Create VCF of filtered SVs.")
ap.add_argument("--sv", help="SV vcf. This should have the sequence of each SV.", required=True)
ap.add_argument("--c1count", help="Class 1 counts file.", required=True)
ap.add_argument("--c2count", help="Class 2 counts file.", required=True)
ap.add_argument("--f", help="filter 0/1")
ap.add_argument("--median", help="median filter", type=int,default=5)

args=ap.parse_args()

#read both count files into dictionary
c1path = args.c1count
c2path = args.c2count

c1dict = {}
c1dict_a_il = {}
c1dict_med = {}

c2dict = {}
c2dict_a_il = {}
c2dict_med = {}

with open(c1path) as f:
    for line in f.read().splitlines():
        if line == '':
            continue
        s_line = line.split()
        r_hg = [int(v) for v in s_line[2].split(',')]
        #r_il = [int(v) for v in s_line[3].split(',')]
        a_hg = [int(v) for v in s_line[4].split(',')]

        a_il = [int(v) for v in s_line[5].split(',')]
        altsum=0
        for i in range(len(a_il )):
            altsum=altsum+a_il[i]
        mean_alt_support=altsum/len(a_il)
        c1dict_a_il[s_line[0]+','+s_line[1]] = round(mean_alt_support)

        c1dict[s_line[0]+','+s_line[1]] = [r_hg,a_hg]
        c1dict_med[s_line[0]+','+s_line[1]]= s_line[-1]
with open(c2path) as f:
    for line in f.read().splitlines():
        if line == '':
            continue
        s_line = line.split()
        #r_hg = [int(v) for v in s_line[2].split(',')]
        #r_il = [int(v) for v in s_line[3].split(',')]
        a_hg = [int(v) for v in s_line[4].split(',')]

        a_il = [int(v) for v in s_line[5].split(',')]
        altsum=0
        for i in range(len(a_il )):
            altsum=altsum+a_il[i]
        mean_alt_support=altsum/len(a_il)
        c2dict_a_il[s_line[0]+','+s_line[1]] = round(mean_alt_support)

        c2dict[s_line[0]+','+s_line[1]] = a_hg
        c2dict_med[s_line[0]+','+s_line[1]]= s_line[-1]
svVCF     = open(args.sv)

#issue: ref/alt fields are not necessarily the ones being used throughout the pipeline
#sometimes the alt sequence does not match (as seen earlier)
#Fix:
#should write new vcf with the 'updated' ref and alt sequences the first time SVs are parsed
#this will save time later so we can just parse the counts and compare sequences to the updated sequences
#and write to a final vcf file

def GetKVPairs(s):
    kvps={}
    v=s.split(";")
    for i in v:
        kvp=i.split("=")
        kvps[kvp[0]] = kvp[1]
    return kvps

svIndex=-1
found = []
not_found = []

outlines = []
for line in svVCF:
    if line[0] == '#':
        outlines.append(line.rstrip('\n'))
        continue
    svIndex+=1
    ref_k_vals = []
    alt_k_vals = []
    
    vals=line.split()
    ref_alt_seq_id = vals[3]+','+vals[4]
    mean_alt=0
    med_alt=0
    if ref_alt_seq_id in c1dict.keys():
        found.append(ref_alt_seq_id)
        ref_kmers = c1dict[ref_alt_seq_id][0]
        alt_kmers = c1dict[ref_alt_seq_id][1]
        mean_alt=c1dict_a_il[ref_alt_seq_id]
        med_alt=c1dict_med[ref_alt_seq_id]
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
    elif ref_alt_seq_id in c2dict.keys():
        found.append(ref_alt_seq_id)
        mean_alt=c2dict_a_il[ref_alt_seq_id]
        med_alt=c2dict_med[ref_alt_seq_id]
    else:
        not_found.append(ref_alt_seq_id)
   # print(str(len(vals)))
    vals.append(".")
    vals.append(".")
    med_alt=round(float(med_alt))
    if args.f:
        if med_alt<args.median:
            #print("med filtered")
            continue


    vals[7]="mean="+str(mean_alt)+";"+"med="+str(med_alt)+";"+"DP="+vals[5]+vals[7]
    
    if len(ref_k_vals) > 0 or len(alt_k_vals) > 0:
        if len(ref_k_vals) > 0:
            vals[7] = vals[7] + ';{}{}'.format('REFKMERIDX=',','.join([str(val) for val in ref_k_vals]))
            vals[7] = vals[7] + ';{}{}'.format('RKN=',len(ref_k_vals))
        if len(alt_k_vals) > 0:
            vals[7] = vals[7] + ';{}{}'.format('ALTKMERIDX=',','.join([str(val) for val in alt_k_vals]))
            vals[7] = vals[7] + ';{}{}'.format('AKN=',len(alt_k_vals))
        




        vals[5]="."
        ref=vals[3]
        alt=vals[4]
        reff=ref[21]
        altt=alt[21]
        vals[3]=reff
        vals[4]=altt
        #print(ref +" "+ref[19:23]+" "+reff +" "+str(len(ref))+" "+alt+" "+alt[19:23]+" "+altt+" "+str(len(alt)))
        outlines.append('\t'.join(vals))#print(reff + " " + altt)


print('Num found: {}'.format(len(found)))
print('Num set found: {}'.format(len(set(found))))
print('Num not found: {}'.format(len(not_found)))
print('Num set not found: {}'.format(len(set(not_found))))

print(args.sv + '.filtered')

with open(args.sv + '.filtered','w') as f:
    f.write( '\n'.join(outlines))
    # info=vals[7]
    # info=vals[7]
    # kv=GetKVPairs(info)
    
    # seq=None
    # qName=None
    # qPos=None
    # if "QNAME" in kv:
        # qName=kv["QNAME"]
    # if "CONTIG" in kv:
        # qName=kv["CONTIG"]

    # if "QSTART" in kv:
        # qPos = int(kv["QSTART"])-1
    # if "CONTIG_START" in kv:
        # qPos = int(kv["CONTIG_START"])

    # if "SEQ" in kv:
        # seq = kv["SEQ"]
    # else:
        # seq=vals[4]
    # svType=None
    # svLen=None
    # chrom=vals[0]
    # if "SVTYPE" in kv:
        # svType = kv["SVTYPE"]
    # if "SVLEN" in kv:
        # svLen = int(kv["SVLEN"])

    # refLen=len(vals[3])
    # altLen=len(vals[4])

    # if svType is None:
        # if refLen > altLen and (refLen >= 50 or altLen >= 50):
            # svType="DEL"
        # elif altLen > refLen and (refLen >= 50 or altLen >= 50):
            # svType ="INS"
            
    # if svLen is None:
        # svLen=abs(refLen-altLen)

    # if seq is None or svType is None:
        # continue
        
    # #svName=str(svIndex)+"/"+vals[0]+"/"+vals[1]+"/"+svType+"/"+str(svLen) #not used, kept for reference
    # pos=int(vals[1])-1
    # alt=None
    
    # if svType == "DEL":
        # if pos < args.kmer:
            # continue
            
        # if qPos is not None and qName is not None:
            # left=max(0,qPos-(args.kmer-1))
            # alt=asm.fetch(qName, left, qPos+args.kmer-1)
        # else:
            # pre=ref.fetch(chrom, pos-(args.kmer-2), pos)
            # suf=ref.fetch(chrom, pos+refLen,pos+refLen+args.kmer-1)
            # alt=pre+vals[4]+suf
        # gen=ref.fetch(chrom,pos-(args.kmer-2),pos+svLen+args.kmer).upper()
    # elif svType == "INS":
        # if qPos is not None and qName is not None:
            # left=max(0,qPos-(args.kmer-1))
            # alt=asm.fetch(qName, left, qPos+(args.kmer-1)+svLen)
        # else:
            # pre=ref.fetch(chrom, pos-(args.kmer-2), pos)
            # suf=ref.fetch(chrom, pos+refLen,pos+refLen+args.kmer-1)
            # alt=pre+vals[4]+suf
        # gen=ref.fetch(chrom,pos-(args.kmer-2),pos+args.kmer).upper()
        
    # alt=alt.upper()
    # alt=alt.replace("N", "A")
    # gen=gen.replace("N", "A")
    
    # queries.write(">" + chrom + "/" + str(pos) + "/ref/" + str(svIndex) + "\n")
    # queries.write(gen + "\n")
    # queries.write(">" + chrom + "/" + str(pos) + "/alt/" + str(svIndex) + "\n")
    # queries.write(alt + "\n")