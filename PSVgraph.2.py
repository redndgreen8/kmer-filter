import sys
import pysam
import argparse
import networkx as nx

ap = argparse.ArgumentParser(description="Build PSV graph from reads and PSV kmers [(k*2)-1 mer]")

ap.add_argument("--reads",required=True,default=sys.stdin)
ap.add_argument("--psv",required=True)
ap.add_argument("--kmer",default=32,type=int)
ap.add_argument("--out",default="psv.gml")
ap.add_argument("--type",type=int)
args = ap.parse_args()

k = args.kmer

PSV = dict()
MG = nx.DiGraph()
MG.add_node("s")
MG.add_node("t")

#DUP_call.queries.fasta
#>chr6/160609166/ref/T/2/19,18
#GCAATTCTTCTATTCTTCTTATACCTAATTCTTTCGGGTCACCTCTGTTGCTACATAAGCCTG
#>chr6/160609155/alt/C/1/16,17
#CAAATTTTCATGCAATTCTTCTATTCTTCTTCTACCTAATTCTTTCGGGTCACCTCTGTTGCT

def storePSVraw(PSV,PSVlist,MG,k):
    l=0
    r = ""
    a = ""
    posi=0

    for n in PSVlist:
        if l%4==0:
            r=n.split("/")[3]
            posi=n.split("/")[1]
            l=l+1
            continue
        if l%4==1:
            l=l+1
            continue
        if l%4==2:
            a=n.split("/")[3]
            l=l+1
            continue
            
        #n = n.rstrip
        MG.add_node(n.rstrip(), pos=posi, ref=r, alt=a )
        PSV[n.rstrip()] = []
        kmer = ""
#        assert k == len(line)
        for i in range(k):
            kmer = n.rstrip()[i:i+k]
            PSV[n.rstrip()].append(kmer)
        l=l+1
    print("stored raw PSV list",str(l))

#filtered.orig_counts.nucfreq    
#mean=37;med=36;DP=24,6;class=1;ref_alt_kmer=ATTCTTCTATTCTTCTTATACCTAATTCTTTCGGGTCACCTCTGTTGCTACATAAGCCTGCTG,ATTCTTCTATTCTTCTTATACCTAATTCTTTTGGGTCACCTCTGTTGCTACATAAGCCTGCTG;REFKMERIDX=;RKN=0

def storePSVfilt(PSV,PSVlist,MG,k):
    l=0
    _refnode = ""
    _node = ""
    info=""
    pre=""
    for n in PSVlist:
        if n[0]=="#":
            continue  
        n = n.rstrip()
        n = n.split()
        info = n[9]
        kmer = info.split(";")[4]
        pre = kmer.split("=")[1].split(",")
        _refnode = pre[0]
        _node = pre[1]
        MG.add_node(_node , pos=n[1], ref=_refnode[k-1], alt=_node[k-1] )
        PSV[_node] = []
        kmer = ""
#        assert k == len(line)
        for i in range(k):
            kmer = _node[i:i+k]
            PSV[_node].append(kmer)
        l=l+1
    print("stored filt PSV list:",str(l))


    
#set but need to keep track of multiplicity
PSVlist=[]
with open(args.psv) as p:
    PSVlist=p.readlines()
#print(PSVlist)
print(MG.number_of_nodes())

if args.type:
    storePSVfilt(PSV, PSVlist, MG,k)
else:
    storePSVraw(PSV, PSVlist, MG,k)

print(MG.number_of_nodes())
#print(PSV)

read_name = ""
offset = 0
starter=0
ender=0
un_read=0
read_stream = open(args.reads)


read_count = 0
for line in read_stream:
    if line[0] == ">":
        read_name = line.rstrip()[1:]
        continue
    read = line.rstrip()
    read_len = len(read)
    first_match = 1
    read_count = read_count+1
    for i in range((read_len - k)+ 1):
        kmer = read[i:i+k]
        for key in PSV.keys():
            current_off=-1
            if kmer in PSV[key]:
                offset = PSV[key].index(kmer)
                
                if first_match and i == read_len - k:
                    MG.add_edge("s", key, rn=read_name, of=offset, ii=i, rl= read_len)
                    MG.add_edge(key, "t", rn=read_name, ii=i, rl= read_len)
                    continue

                if first_match:
                    MG.add_edge("s", key, rn=read_name, of=offset, ii=i,rl= read_len)
                    current_match = key
                    current_off = offset
                    first_match = 0
                    starter = starter +1 
                    continue

                if  i == read_len - k:
                    MG.add_edge(key, "t", rn=read_name, ii=i,rl= read_len)
                    ender = ender +1 
                    continue
                                                
                MG.add_edge(current_match, key, rn=read_name, ii=i, iof = current_off,rl= read_len, of=offset)
                current_match = key
                current_off = offset
                
            elif not first_match and i == read_len - k:
                MG.add_edge(current_match, "t", rn=read_name, ii=i , iof = current_off,rl= read_len)
                ender = ender +1 

    if first_match:
        un_read = un_read+1            
    if read_count%10 == 0:
        #print(read_count)
        print(MG.number_of_edges())
    if read_count%100 == 0:
        print(read_count)

output_dir = args.out 
nx.write_gml(MG,output_dir)


print(MG.number_of_nodes())
print(MG.number_of_edges())
print("reads: ",str(read_count),"discard reads: ",str(un_read))
print("source: ",str(starter),"\nsink: ",str(ender))




