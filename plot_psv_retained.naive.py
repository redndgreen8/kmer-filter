import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--bamFreqOut")
parser.add_argument("--filtered",default=sys.stdin)
parser.add_argument("--retained",default=None)
parser.add_argument("--mean")
parser.add_argument("--prefix")

args = parser.parse_args()


w=open(args.bamFreqOut)

#f=open("/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/chr9:63809200-63827100_17/kmer_filter/nucplot_PSVcandidate.txt2")
#b=open("/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/chr9:63809200-63827100_17/kmer_filter/retained_PSV.txt")
#w=open("/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/chr9:63809200-63827100_17/kmer_filter/chr9:63809200-63827100_17.bamToFreq.DUP_CALL.out")


mean=int(args.mean)
STD = np.sqrt(mean)
COV=mean
MINCOV = int( max(COV - 4*STD, COV/2.0) )
MAXCOV = int( COV + 1*STD )
MINTOTAL = int( 2*COV - 3*STD )


firstpass  = []
secondpass=[]
truepospass=[]
thirdp=[]
fourp=[]

#alt_il = []

firstf  = []
secondf = []
trueposf= []
thirdf=[]
fourf=[]

firstr  = []
secondr = []
trueposr= []


#bamfreqout
i=0
for line in w:
    line = line.split()
    i=i+1
    trueposf.append(int(line[1]))
    bases = []
    for basepair in line[2:6]:
        bases.append(int(basepair))
    bases = sorted(bases, reverse=True)
    firstf.append(bases[0])
    secondf.append(bases[1])
    thirdf.append(bases[2])
    fourf.append(bases[3])
w.close()

print( str(trueposf[0]),str(i) )

if(args.retained is not None):
            b=open(args.retained)
            retain=b.readlines()
            b.close()   

#filtered psv candidates
if(args.filtered is not None):
    f=open(args.filtered)
    for line in f:
        line = line.split()
        truepospass.append(int(line[1]))
        bases = []
        for basepair in line[3:7]:
            bases.append(int(basepair))
        bases = sorted(bases, reverse=True)
        firstpass.append(bases[0])
        secondpass.append(bases[1])
        thirdp.append(bases[2])
        fourp.append(bases[3])

   #     alt=line[9].split(";")[0].split("=")[1]
    #    alt_il.append( alt)
        if(args.retained is not None):
            for liner in retain:
                liner = liner.split()
              # print( liner[0] , curpos,trueposf[0],int(line[1]),trueposf[0]+ int(liner[0]) )
                if trueposf[0]+ int(liner[0]) == int(line[1]):
             #   print( liner[0] , trueposf[0],int(line[1]),trueposf[0]+ int(liner[0]) )
                    trueposr.append(int(line[1]))
                    secondr.append(bases[1])
                    firstr.append(bases[0])
                else:
                    continue

if(args.retained is not None):
    print(trueposr[0],firstr[0],secondr[0])
    print(len(retain))
    print(len(trueposr),len(firstr),len(secondr))


secondnp=np.array(secondf)
trueposnp=np.array(trueposf)

if(args.filtered is not None):
    print(str(firstpass[0]))
 #   print(str(alt_il[0]))
 #   alt_ilnp=np.array(alt_il,dtype=float)
    fp=np.array(firstpass)
    sp=np.array(secondpass)
    tp=np.array(thirdp)
    fop=np.array(fourp)

firstnp=np.array(firstf)


print(secondnp.max())
print(firstnp.max())

#fig = plt.figure()

#plt.hist(secondf, bins=(secondnp.max()-1))
#plt.xlabel('alt alleles binned count')
#plt.ylabel('Frequency')
#plt.axis([0, mean*1.25, 0,10 ])
pre=args.prefix
prehist=pre+".alt_il_kmer_hist.png"
#plt.savefig(prehist)

#0.01 threshold

ff=np.array(firstf)
sf=np.array(secondf)
tf=np.array(thirdf)
fof=np.array(fourf)


fig, ax = plt.subplots( figsize=(16,9) )
#cov, = plt.plot(trueposf, ff+sf, 'o', color="brown", markeredgewidth=0.0, markersize=1.5, label = "coverage profile")

if(args.filtered is not None):
#    cov2p, = plt.plot(truepospass, fp+sp, 'o', color="purple", markeredgewidth=0.0, markersize=1.5)
  #  cov3p, = plt.plot(truepospass, fp+sp+tp, 'o', color="yellow", markeredgewidth=0.0, markersize=1.5)
    cov4p, = plt.plot(truepospass, fp+sp+tp+fop, 'o', color="black", markeredgewidth=0.0, markersize=1.5)
    tri, = plt.plot(truepospass, secondpass , 'o', color="blue",   markeredgewidth=0.0, markersize=2, label = "second most pass")
   # for i in range(len(truepospass)):
    #    round_alt=round(float(alt_il[i]))
     #   if round_alt>mean*1.25 or round_alt<min(4,0.1*mean):
      #      fnt=3
       # else:
        #    fnt=7.5
        #plt.text(truepospass[i], secondpass[i],str(round_alt), color="blue",  fontsize=fnt)



sec, = plt.plot(trueposf, secondf , 'o', color="red",    markeredgewidth=0.0, markersize=1.5, label = "second most fail")
cov4, = plt.plot(trueposf, ff+sf+tf+fof, 'o', color="black", markeredgewidth=0.0, markersize=1.5, label="total coverage")
#cov3, = plt.plot(trueposf, ff+sf+tf, 'o', color="yellow", markeredgewidth=0.0, markersize=1.5)
#cov2, = plt.plot(trueposf, ff+sf, 'o', color="purple", markeredgewidth=0.0, markersize=1.5)

#prime, = plt.plot(trueposf, firstf, 'o', color="black", markeredgewidth=0.0, markersize=1.5, label = "most frequent base pair")

if(args.retained is not None):
    psec, = plt.plot(trueposr, secondr , 's', color="green",    markeredgewidth=0.0, markersize=5,label = "second most retained")






ax.axhline(y=MINCOV,linewidth=0.5)
ax.axhline(y=MAXCOV,linewidth=0.5)
ax.axhline(y=MINTOTAL,linewidth=0.5)
ax.axhline(y=mean,linewidth=0.5,color="purple")






ax.set_xlabel('BP Position')
ax.set_ylabel('Depth')

ylabels = [format(label, ',.0f') for label in ax.get_yticks()]
xlabels = [format(label, ',.0f') for label in ax.get_xticks()]
ax.set_yticklabels(ylabels)
ax.set_xticklabels(xlabels)

# Hide the right and top spines
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.legend()
preplot=pre+"psv_plot.png"
plt.savefig(preplot)
