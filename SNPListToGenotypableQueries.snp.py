#!/usr/bin/env python


import argparse
import pysam

ap = argparse.ArgumentParser(description="Create kmer table from SVs.")
ap.add_argument("--snp", help="PSV vcf. This should have the sequence of each PSV.", required=True)
#ap.add_argument("--asm", help="Assembly calling sv", required=True)
ap.add_argument("--ref", help="Reference", required=True)
ap.add_argument("--kmer", help="Kmer size (32)",type=int,default=32)
ap.add_argument("--queries", help="File to write queries to.", required=True)

args=ap.parse_args()

snpVCF     = open(args.snp)
ref       = pysam.FastaFile(args.ref)
#asm       = pysam.FastaFile(args.asm)

def GetKVPairs(s):
    kvps={}
    v=s.split(";")
    for i in v:
        kvp=i.split("=")
        kvps[kvp[0]] = kvp[1]
    return kvps

Index=-1

queries = None
if args.queries is not None:
    queries = open(args.queries, 'w')

    
correctedSeq_lines = []
for line in snpVCF:
    if line[0] == '#':
        correctedSeq_lines.append(line.rstrip('\n')) #remove newlines from header lines
        continue
    Index+=1
    
    vals=line.split()

    
    chrom=vals[0]

#    refLen=len(vals[3])# multiple entries?? , separated
    alt_allele=vals[4].split(",")
    altLen=len( alt_allele )

            
  
        
    #svName=str(svIndex)+"/"+vals[0]+"/"+vals[1]+"/"+svType+"/"+str(svLen) #not used, kept for reference
    
    qpos=int(vals[1])-1
    pos=int(vals[1])
    alt=None


    if altLen>1:
      continue
    else:
        pre=ref.fetch(chrom, pos-(args.kmer-1), pos)
        suf=ref.fetch(chrom, pos+1,pos+1+args.kmer-1)

        gen=ref.fetch(chrom,pos-(args.kmer-1),pos+args.kmer).upper()
        gen=gen.replace("N", "A")
      #for i in  range(altLen):

        alt=pre+vals[4]+suf
        

        alt=alt.upper()
        alt=alt.replace("N", "A")

        #print( str(len(pre)) +" "+str(len(suf)) + " "+ pre +" " + vals[4] + " " + suf +" " +str(len(gen)) + " " +gen[0:21]+" "+gen[21]+" "+gen[22:] +" "+ vals[3]+" refFrominput:"+ ref.fetch(chrom,pos,pos+1 )  )
        assert (len(gen)-((2*args.kmer)-1)) ==0, "ref seq doesnt match kmer seq size"
        assert (len(alt)- ((2*args.kmer)-1)) ==0, "alt seq doesnt match kmer seq size"

        queries.write(">" + chrom + "/" + str(qpos+1) + "/ref/" +gen[args.kmer-1]+"/"+str(Index) +"/"+str(vals[5])+ "\n")
        queries.write(gen + "\n")
        queries.write(">" + chrom + "/" + str(qpos+1) + "/alt/" +vals[4]+"/"+str(Index) + "/"+str(vals[5])+ "\n")
        queries.write(alt + "\n")
      

        new_vals = vals
        new_vals[3] = gen
        new_vals[4] = alt
        correctedSeq_lines.append('\t'.join(new_vals))
queries.close()

vcf_path = args.snp
#print('\n'.join(correctedSeq_lines))
with open(vcf_path+'.correctedseq','w') as f:
    f.write('\n'.join(correctedSeq_lines))