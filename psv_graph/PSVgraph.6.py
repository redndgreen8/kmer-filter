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

        if _node not in node_multi.keys():
            node_multi[_node]=1
        else:
            node_multi[_node] = node_multi[_node] + 1
    
        MG.add_node(_node , pos=n[1], ref=_refnode[k-1], alt=_node[k-1],t=0, mult= node_multi[_node] )
    


#        assert k == len(line)
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

node_multi = {}

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
prev_edge=0
    

edge_list={}    
for line in read_stream:
    if line[0] == ">":
        line = line.lstrip()
        read_name = "/".join(line.split("/")[:2])
        continue
    read = line.rstrip()
    read_len = len(read)
    first_match = 1
    read_count = read_count+1
    prev_match = ""

    sentinel = "_".join(["t", read_name ])
    _starter = "_".join(["s", read_name ])

    edge_list[read_name]=[]

    for i in range((read_len - k)+ 1):
        kmer = read[i:i+k]


        for key in PSV.keys():
            prev_off = -1
            prev_kmer = ""
            if kmer in PSV[key]:
                offset = PSV[key].index(kmer)


                edge_list[read_name].append(key)    
                

                if first_match and i == read_len - k:
                    MG.add_node(_starter,t=1)
                    MG.add_node(sentinel,t=2)
                    first_match = 0
                    MG.add_edge(_starter, key, rn=read_name, of=offset, ri=i, rl= read_len)
                    MG.add_edge(key, sentinel, rn=read_name, of=offset, ri=i, rl= read_len)
                    continue

                if first_match:
                    MG.add_node(_starter,t=1)  
                    MG.add_edge(_starter, key, rn=read_name, of=offset, ri=i,rl= read_len)
                    prev_match = key
                    prev_kmer = kmer
                    prev_off = offset
                    first_match = 0
                    starter = starter +1 
                    continue

                if  i == read_len - k:
                    MG.add_node(sentinel,t=2)
                    MG.add_edge(key, sentinel, rn=read_name, of=offset, ri=i,rl= read_len)
                    ender = ender +1 
                    continue
                
                if prev_match==key:
                    if int(node_multi[key]) > 1:
                        MG.add_edge(prev_match, key)#, rn=read_name, ri=i, iof = prev_off,rl= read_len, of=offset)
                    else:
                        continue
                    continue


                MG.add_edge(prev_match, key, rn=read_name, ri=i, iof = prev_off,rl= read_len, of=offset)

                prev_kmer = kmer
                prev_match = key
                prev_off = offset
                
            elif not first_match and i == read_len - k:
                MG.add_node(sentinel,t=2)
                MG.add_edge(prev_match, sentinel, rn=read_name, ri=i , iof = prev_off,rl= read_len)
                ender = ender +1 

    if first_match:
        un_read = un_read+1            
    x = MG.number_of_edges()
    if int(x) - prev_edge >50:
        print(read_count , MG.number_of_edges())

    prev_edge = x

rem=[]
node_list={}
for node in MG.nodes():
    node_list[node]=[MG.in_degree(node),MG.out_degree(node)]
    if MG.degree(node) == 0:
        rem.append(node)
MG.remove_nodes_from(rem)
print("removed nodes:",len(rem))

output_dir = args.out 
nx.write_gml(MG,output_dir)

odir="/" + "/".join(output_dir.split("/")[:-2])
file_name=output_dir.split("/")[-1]


edge_o=odir+"/"+ ".".join(file_name.split(".")[:-2]) + "." +"edgelist"
node_o=odir+"/"+ ".".join(file_name.split(".")[:-2]) + "." +"nodelist"

print(edge_o,"\n",node_o)


with open(edge_o,'w') as e:
    e.write(str(edge_list))
with open(node_o,'w') as n:
    n.write(str(node_list))


print(MG.number_of_nodes())
print(MG.number_of_edges())


print("reads: ",str(read_count),"discard reads: ",str(un_read))
print("source: ",str(starter),"\nsink: ",str(ender))




