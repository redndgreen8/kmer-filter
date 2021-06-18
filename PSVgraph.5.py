import sys
import math
import argparse
import networkx as nx
from collections import OrderedDict

ap = argparse.ArgumentParser(description="Build PSV graph from reads and PSV kmers [(k*2)-1 mer]")

ap.add_argument("--reads",required=True,default=sys.stdin)
#unmapped

ap.add_argument("--cov",required=True,type=int)

ap.add_argument("--psv",required=True)
ap.add_argument("--kmer",default=32,type=int)
ap.add_argument("--out",default="psv.gml")
#ap.add_argument("--type",type=int)#raw 0 , filt
args = ap.parse_args()

k = args.kmer
cov = args.cov
PSV = OrderedDict()
MG = nx.MultiDiGraph()
#MG.add_node("s")
#MG.add_node("t")





def storePSVfilt(PSV,PSVlist,MG,k):
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
        MG.add_node(_node , pos=n[1], ref=_refnode[k-1], alt=_node[k-1],t=0 )

        for i in range(k):
            kmer = _node[i:i+k]
            PSV[_node].append(kmer)
   #         MG.add_node(kmer , pos=n[1], ref=_refnode[k-1], alt=_node[k-1] )
            l=l+1
    print("stored filt PSV list:",str(l))



    
#set but need to keep track of multiplicity
PSVlist=[]
with open(args.psv) as p:
    PSVlist=p.readlines()
#print(PSVlist)
print(MG.number_of_nodes())

storePSVfilt(PSV, PSVlist, MG,k)


print(MG.number_of_nodes())
#print(PSV)

read_name = ""
offset = 0
starter=0
ender=0
un_read=0
read_stream = open(args.reads)


read_count = 0
node_multi = {}
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
    _starter = "_".join(["s", read_name ])


    for i in range((read_len - k)+ 1):
        kmer = read[i:i+k]


        for key in PSV.keys():
            current_off=-1
            if kmer in PSV[key]:
                offset = PSV[key].index(kmer)

                if MG.has_node(kmer):
                    node_multi[kmer] = node_multi[kmer] + 1
                    MG.add_node(kmer, mult=node_multi[kmer])
                else:
                    node_multi[kmer] = 1
                    MG.add_node(kmer, mult=1)
                    
                    
                if first_match and i == read_len - k:
                    MG.add_node(_starter,t=1)
                    MG.add_node(sentinel,t=2)
                    MG.add_edge(_starter, key, rn=read_name, of=offset, ii=i, rl= read_len)
                    MG.add_edge(key, sentinel, rn=read_name, ii=i, rl= read_len)
                    continue

                if first_match:
                    MG.add_node(_starter,t=1)  
                    MG.add_edge(_starter, key, rn=read_name, of=offset, ii=i,rl= read_len)
                    current_match = key
                    current_off = offset
                    first_match = 0
                    starter = starter +1 
                    continue

                if  i == read_len - k:
                    MG.add_node(sentinel,t=2)
                    MG.add_edge(key, sentinel, rn=read_name, ii=i,rl= read_len)
                    ender = ender +1 
                    continue
                                                
                MG.add_edge(current_match, key, rn=read_name, ii=i, iof = current_off,rl= read_len, of=offset)
                current_match = key
                current_off = offset
                
            elif not first_match and i == read_len - k:
                MG.add_node(sentinel,t=2)
                MG.add_edge(current_match, sentinel, rn=read_name, ii=i , iof = current_off,rl= read_len)
                ender = ender +1 

    if first_match:
        un_read = un_read+1            
    if read_count%10 == 0:
        #print(read_count)
        print(MG.number_of_edges())
    if read_count%100 == 0:
        print(read_count)


rem=[]

for node in MG.nodes():
    if MG.degree(node) == 0:
        rem.append(node)
MG.remove_nodes_from(rem)
print(len(rem))

output_dir = args.out 
nx.write_gml(MG,output_dir)


print(MG.number_of_nodes())
print(MG.number_of_edges())


print("reads: ",str(read_count),"discard reads: ",str(un_read))
print("source: ",str(starter),"\nsink: ",str(ender))




