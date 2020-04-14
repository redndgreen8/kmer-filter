#!/usr/bin/env python


import argparse
import pysam

ap = argparse.ArgumentParser(description="Create kmer table from SVs.")
ap.add_argument("--sv", help="SV vcf. This should have the sequence of each SV.", required=True)
ap.add_argument("--asm", help="Assembly calling sv", required=True)
ap.add_argument("--ref", help="Reference", required=True)
ap.add_argument("--kmer", help="Kmer size (22)",type=int, default=22)
ap.add_argument("--queries", help="File to write queries to.", required=True)

args=ap.parse_args()

svVCF     = open(args.sv)
ref       = pysam.FastaFile(args.ref)
asm       = pysam.FastaFile(args.asm)

def GetKVPairs(s):
    kvps={}
    v=s.split(";")
    for i in v:
        kvp=i.split("=")
        kvps[kvp[0]] = kvp[1]
    return kvps

svIndex=-1

queries = None
if args.queries is not None:
    queries = open(args.queries, 'w')

    
correctedSeq_lines = []
for line in svVCF:
    if line[0] == '#':
        correctedSeq_lines.append(line.rstrip('\n')) #remove newlines from header lines
        continue
    svIndex+=1
    
    vals=line.split()
    info=vals[7]
    kv=GetKVPairs(info)
    
    seq=None
    qName=None
    qPos=None
    if "QNAME" in kv:
        qName=kv["QNAME"]
    if "CONTIG" in kv:
        qName=kv["CONTIG"]

    if "QSTART" in kv:
        qPos = int(kv["QSTART"])-1
    if "CONTIG_START" in kv:
        qPos = int(kv["CONTIG_START"])

    if "SEQ" in kv:
        seq = kv["SEQ"]
    else:
        seq=vals[4]
    svType=None
    svLen=None
    chrom=vals[0]
    if "SVTYPE" in kv:
        svType = kv["SVTYPE"]
    if "SVLEN" in kv:
        svLen = int(kv["SVLEN"])

    refLen=len(vals[3])
    altLen=len(vals[4])

    if svType is None:
        if refLen > altLen and (refLen >= 50 or altLen >= 50):
            svType="DEL"
        elif altLen > refLen and (refLen >= 50 or altLen >= 50):
            svType ="INS"
            
    if svLen is None:
        svLen=abs(refLen-altLen)

    if seq is None or svType is None:
        continue
        
    #svName=str(svIndex)+"/"+vals[0]+"/"+vals[1]+"/"+svType+"/"+str(svLen) #not used, kept for reference
    pos=int(vals[1])-1
    alt=None
    
    if svType == "DEL":
        if pos < args.kmer:
            continue
            
        if qPos is not None and qName is not None:
            left=max(0,qPos-(args.kmer-1))
            alt=asm.fetch(qName, left, qPos+args.kmer-1)
        else:
            pre=ref.fetch(chrom, pos-(args.kmer-2), pos)
            suf=ref.fetch(chrom, pos+refLen,pos+refLen+args.kmer-1)
            alt=pre+vals[4]+suf
        gen=ref.fetch(chrom,pos-(args.kmer-2),pos+svLen+args.kmer).upper()
    elif svType == "INS":
        if qPos is not None and qName is not None:
            left=max(0,qPos-(args.kmer-1))
            alt=asm.fetch(qName, left, qPos+(args.kmer-1)+svLen)
        else:
            pre=ref.fetch(chrom, pos-(args.kmer-2), pos)
            suf=ref.fetch(chrom, pos+refLen,pos+refLen+args.kmer-1)
            alt=pre+vals[4]+suf
        gen=ref.fetch(chrom,pos-(args.kmer-2),pos+args.kmer).upper()
        
    alt=alt.upper()
    alt=alt.replace("N", "A")
    gen=gen.replace("N", "A")
    
    queries.write(">" + chrom + "/" + str(pos) + "/ref/" + str(svIndex) + "\n")
    queries.write(gen + "\n")
    queries.write(">" + chrom + "/" + str(pos) + "/alt/" + str(svIndex) + "\n")
    queries.write(alt + "\n")
    
    new_vals = vals
    new_vals[3] = gen
    new_vals[4] = alt
    correctedSeq_lines.append('\t'.join(new_vals))
    
vcf_path = args.sv
# print('\n'.join(correctedSeq_lines))
with open(vcf_path+'.correctedseq','w') as f:
    f.write('\n'.join(correctedSeq_lines))