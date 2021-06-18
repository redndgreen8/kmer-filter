import sys
import math
import argparse
import networkx as nx
from collections import OrderedDict

ap = argparse.ArgumentParser(description="Build PSV graph from reads and PSV kmers [(k*2)-1 mer]")

ap.add_argument("--reads",required=True,default=sys.stdin)
ap.add_argument("--cov",required=True,type=int)

ap.add_argument("--psv",required=True,default=sys.stdin)
ap.add_argument("--kmer",default=32,type=int)
ap.add_argument("--out",default="psv.gml")
ap.add_argument("--type",type=int)#raw 0 , filt
args = ap.parse_args()

k = args.kmer
cov = args.cov
PSV = OrderedDict()
MG = nx.MultiDiGraph()
MG.add_node("s")
#MG.add_node("t")


#>chr6/160609166/ref/T/2/19,18
#GCAATTCTTCTATTCTTCTTATACCTAATTCTTTCGGGTCACCTCTGTTGCTACATAAGCCTG
#>chr6/160609155/alt/C/1/16,17
#CAAATTTTCATGCAATTCTTCTATTCTTCTTCTACCTAATTCTTTCGGGTCACCTCTGTTGCT

def storePSVraw(PSV,PSVlist,MG,k,cov,node_multi):
    print("storing PSV Dict!!!")
   # print(PSVlist)
    l=0
    r = ""
    a = ""
    posi=0
    p = ""
    ct=0
    gate=False
    s=0
    for n in PSVlist:
      #  n=n.rstrip()
        if l%4==0:
            p = n.split("/")[5]
            ct = p.split(",")[1]
           # print(ct)
            if int(ct) > math.floor(0.15 * cov):
                r=n.split("/")[3]
                posi=n.split("/")[1]
                gate=True            
                s=s+1                
            l=l+1
            continue
        if l%4==1:
            l=l+1
            continue
        if gate:
            if l%4==2:
                a=n.split("/")[3]
                l=l+1
                continue   
       #     MG.add_node(n.rstrip(), pos=posi, ref=r, alt=a )
            PSV[n.rstrip()] = []
            kmer = ""
            for i in range(k):
                kmer = n.rstrip()[i:i+k]
                PSV[n.rstrip()].append(kmer)
            l=l+1
        else:
            l=l+1
            continue
            
    print("stored raw PSV list",str(l),str(s))

    
#chr6	160609169	160609170	0	24	0	6	0	0	mean=37;med=36;DP=24,6;class=1;ref_alt_kmer=ATTCTTCTATTCTTCTTATACCTAATTCTTTCGGGTCACCTCTGTTGCTACATAAGCCTGCTG,ATTCTTCTATTCTTCTTATACCTAATTCTTTTGGGTCACCTCTGTTGCTACATAAGCCTGCTG;REFKMERIDX=;RKN=0
def sstorePSVfilt(PSV,PSVlist,MG,k):
    l=0
    _node = ""
    po = ""
    for n in PSVlist:
        if n[0]=="#":
            continue  
        po = n.split()[1]
        g = n.split()[9]
        f = g.rstrip()
        s = f.split(";")[4]
        m = s.split("=")[1]

        r = m.split(",")[0]
        _node = m.split(",")[1]
      #  print("\n".join([po,g,f,s,m,r,_node]))


        #MG.add_node(_node , pos=n[1], ref=n[3].split()[k-1], alt=_node.split()[k-1] )
        PSV[_node] = []
        kmer = ""
#        assert k == len(line)
        for i in range(k):
            kmer = _node.split()[i:i+k]
            PSV[_node].append(kmer)
  #          MG.add_node(kmer , pos=po)#, ref=r[k-1], alt=_node[k-1] )
        l=l+1
    print("stored filt PSV list",str(l),PSV)

def ssstorePSVfilt(PSV,PSVlist,MG,k):
    l=0
    _node = ""    
    r = ""
    a = ""
    for n in PSVlist:
        if n[0]=="#":
            continue  
        n = n.rstrip()
        n = n.split()
        #print(n)
        _node = n[4]
        r = n[3].split()[k-1]
        a = _node.split()[k-1]
        MG.add_node(_node , pos=n[1], ref=r, alt=a )
        PSV[_node] = []
        kmer = ""
#        assert k == len(line)
        for i in range(k):
            kmer = _node.split()[i:i+k]
            PSV[_node].append(kmer)
        l=l+1
    print("stored filt PSV list",str(l))
    
def storePSVfilt(PSV,PSVlist,MG,k,node_multi):
    l=0
    _refnode = ""
    _node = ""
    info=""
    pre=""
    for n in PSVlist:
        if n[0]=="#":
            continue  
       # print(n)
        n = n.rstrip()
        n = n.split()
      #  print(n)
        info = n[9]
        kmer = info.split(";")[4]
        pre = kmer.split("=")[1].split(",")
        _refnode = pre[0]
        _node = pre[1]
        PSV[_node] = []
        kmer = ""
#        assert k == len(line)
        if MG.has_node(_node):
            node_multi[_node] = node_multi[_node] + 1
            MG.add_node(_node, mult=node_multi[_node])#, alt=a )
        else:
            node_multi[_node] = 1
            MG.add_node(_node, mult=1 ,pos=n[1], ref=_refnode[k-1], alt=_node[k-1])
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

node_multi = {}


if args.type:
    storePSVfilt(PSV, PSVlist, MG,k,node_multi)
else:
    storePSVraw(PSV, PSVlist, MG,k,cov,node_multi)

print(MG.number_of_nodes())
#print(PSV)

read_name = ""
offset = 0
starter=0
ender=0
un_read=0
read_stream = open(args.reads)


read_count = 0
for key in PSV.keys():
    for kmer in PSV[key]:
        node_multi[kmer]=0

for line in read_stream:
    if line[0] == ">":
        line = line.lstrip()
        read_name = "/".join(line.split("/")[:2])
        continue
    read = line.rstrip()
    read_len = len(read)
    first_match = 1
    read_count = read_count+1

    sentinel = "_".join(["t", read_name ])
    
    for i in range((read_len - k)+ 1):
        kmer = read[i:i+k]
     #   a=""

        for key in PSV.keys():
            current_off=-1
            
            if kmer in PSV[key]:
                offset = PSV[key].index(kmer)
            #    print(kmer)



                if first_match and i == read_len - k:
                    MG.add_node(sentinel)
                    MG.add_edge("s", key, rn=read_name, of=offset, ii=i, rl= read_len)
                    MG.add_edge(key, sentinel, rn=read_name, ii=i, rl= read_len)
                    continue

                if first_match:
                    MG.add_edge("s", key, rn=read_name, of=offset, ii=i,rl= read_len)
                    current_match = key
                    current_off = offset
                    first_match = 0
                    starter = starter +1 
                    continue

                if  i == read_len - k:
                    MG.add_node(sentinel)
                    MG.add_edge(key, sentinel, rn=read_name, ii=i,rl= read_len)
                    ender = ender +1 
                    continue
                                                
                MG.add_edge(current_match, key, rn=read_name, ii=i, iof = current_off,rl= read_len, of=offset)
                current_match = key
                current_off = offset
                
            elif not first_match and i == read_len - k:
                MG.add_node(sentinel)
                MG.add_edge(current_match, sentinel, rn=read_name, ii=i , iof = current_off,rl= read_len)
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


print("\nnodes:",MG.number_of_nodes())
print("edges:",MG.number_of_edges())
print("reads: ",str(read_count),"discard reads: ",str(un_read))
print("source: ",str(starter),"\nsink: ",str(ender))




