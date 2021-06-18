import os
import tempfile
import numpy as np

#
# Locate the tempdir for grid processing
#
#configfile: "SVGenotyping.json"

#if "TMPDIR" in os.environ:
#    TMPDIR = os.environ['TMPDIR']
#elif "TMPDIR" in config:
#    TMPDIR = config['TMPDIR']
#else:
#    TMPDIR = tempfile.gettempdir()

SD =config['SD']#/home/cmb-16/mjc/rdagnew/summerproj/kmer_filter"#os.path.dirname(workflow.snakefile)
WD=config['WD']#"/panfs/qcb-panasas/rdagnew/HG00514/ont"

#config['CD']

#SDA="/home/cmb-16/mjc/rdagnew/summerproj/SDA/SDA"
#
# 
# samples = config["samples"].keys()
#samples = ["HG00514"]

loci_list=config['ll']
#filtering criteria
min_kmer_count = 5
max_kmer_count = 34
min_ill_kmer_count = 5
max_ill_kmer_count = 120

listl=[]
#listfull=[]
i=0
with open(loci_list,'r') as l:
    for r in l.readlines():
        if i<250 and i>=0:
            r=r.rstrip()
            r=r.split()
            r=r[0] + ":" + r[1] +"-"+r[2]
            listl.append(r)
        i=i+1

        

means=config['mean']
COV=means
# standard deviation of read depth
STD = np.sqrt(int(means))
MINCOV = int( max(COV - 4*STD, COV/2.0) )
MAXCOV = int( COV + 1*STD )
MINTOTAL = int( 2*COV - 3*STD )
MINREADS = int(MINCOV/2.0)
        
sbs= (MINCOV+2)
fss=(MINTOTAL - sbs) +3





#when input gets too large 
# no. queries/psv sites * k  >50 million
# no. queries >1.5 million
# no. of loci >20

part=0 
if len(listl)>20:
    part=1
    print("partition il and hg38")

rtype=config['read_type']
#bams=config['bam']

print(rtype)
#count directory
CD=WD+"/"+rtype+"/"




#def partitionss(listl):
rule all:
    input:
        freqout=expand("{WD}/{loci}/{r_type}/kmer_filter/bamToFreq.DUP_CALL.out", WD=WD,loci=listl,r_type=rtype),
        vcff=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.vcf",loci=listl ,WD=WD,r_type=rtype),
        queryFasta=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta", loci=listl, WD=WD,r_type=rtype),
       # preNucfreqdone=expand("{WD}/pre.assembly.consensus.nucfreq.done",loci=listl, WD=WD),

        resultsFile=expand("{WD}/DUP_call.queries.counts",loci=listl ,WD=WD,r_type=rtype),
        resultsFileccs=expand("{WD}/DUP_call.queries.counts.ccs",loci=listl ,WD=WD,r_type=rtype),

       # resultsFileil=expand("{WD}/DUP_call.queries.counts.il",loci=listl ,WD=WD),
       # resultsFilehg=expand("{WD}/DUP_call.queries.counts.hg",loci=listl ,WD=WD),

        loci_bed=expand("{WD}/{loci}/{r_type}/loci.bed",loci=listl ,WD=WD,r_type=rtype),
        preNucfreq=expand("{WD}/{loci}/{r_type}/kmer_filter/pre.assembly.nucfreq",loci=listl ,WD=WD,r_type=rtype),



rule bamToFreq:
    input:
        bam=config['bam'],
    output:
        freqout="{WD}/{loci}/{r_type}/kmer_filter/bamToFreq.DUP_CALL.out",
    params:
        lo="{loci}",
    shell:"""
$sum/SegDupSNV/bamToFreq {input.bam} <(echo {params.lo}|tr "_" "\t"|cut -f 1) $hg38/hg38.fa > {output.freqout}
    
    echo bamfreqdone

    """

rule vcf:
    input:
        freqout="{WD}/{loci}/{r_type}/kmer_filter/bamToFreq.DUP_CALL.out",
    output:
        vcff="{WD}/{loci}/{r_type}/kmer_filter/DUP_call.vcf",
    params:
        sd=SD,
    shell:"""
cat {input.freqout} | $py34 {params.sd}/bamfreq.StreamingKmerFilter.py --nucfreq /dev/stdin --ref $hg38/hg38.fa > {output.vcff}

    """


rule GenerateQueriesFasta:
    input:
        vcff="{WD}/{loci}/{r_type}/kmer_filter/DUP_call.vcf",
    output:
        queryFasta="{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta",
    params:
        sd=SD,
        wd=WD,
        rt=rtype,
    shell:"""
python3 {params.sd}/SNPListToGenotypableQueries.snp.py --snp {input.vcff} --ref $hg38/hg38.fa --queries {output.queryFasta}
echo snplisttogenotypedone


mkdir -p {params.wd}/{params.rt}
"""
#when input gets too large 
# no. queries/psv sites * k  >50 million
#no. queries >1.5 million
if part:
    rule GenerateResultsFromQuerieshg:
        input:
            queryFasta=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD,r_type=rtype),
        output:
            resultsFilehg="{WD}/{r_type}/DUP_call.queries.counts.hg",
        threads:16
        params:
            hgjf='/panfs/qcb-panasas/rdagnew/hg38.no_alts.32.jf',
            sd=SD,
        shell:"""
        python3 {params.sd}/SVStreamingFilter.hg.py --bpseqs <(cat {input.queryFasta}) --hgjf {params.hgjf} --output {output.resultsFilehg}
        """
    rule GenerateResultsFromQueriesil:
        input:
            queryFasta=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD,r_type=rtype),
            resultsFilehg="{WD}/{r_type}/DUP_call.queries.counts.hg",
        output:
            resultsFileil="{WD}/{r_type}/DUP_call.queries.counts.il",
            resultsFile="{WD}/{r_type}/DUP_call.queries.counts",
        threads:16
        params:
            jf=config['iljf'],
            sd=SD,
        shell:"""
        python3 {params.sd}/SVStreamingFilter.il.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --output {output.resultsFileil}
   
        paste <(cat {input.resultsFilehg} |cut -f 1-6) <(cat {output.resultsFileil}|cut -f 5-7) > {output.resultsFile}
        
        wc -l {input.resultsFilehg} {output.resultsFileil} {output.resultsFile}
        """
else:


    rule GenerateResultsFromQueries:
        input:
            queryFasta=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD,r_type=rtype),
        output:
            resultsFile="{WD}/{r_type}/DUP_call.queries.counts",
        threads:16
        params:
            jf=config['iljf'],
            sd=SD,
            hgjf='/panfs/qcb-panasas/rdagnew/hg38.no_alts.32.jf',
        shell:"""
        python3 {params.sd}/SVStreamingFilter.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --hgjf {params.hgjf} --output {output.resultsFile}

        """

rule GenerateResultsFromQueriesccs:
    input:
        queryFasta=expand("{WD}/{loci}/{r_type}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD,r_type=rtype),
        resultsFile="{WD}/{r_type}/DUP_call.queries.counts",
    output:
        resultsFileccs="{WD}/{r_type}/DUP_call.queries.counts.ccs",
#        resultsFile="{WD}/DUP_call.queries.counts",
    threads:16
    params:
        jf=config['ccsjf'],
        sd=SD,
    shell:"""
    python3 {params.sd}/SVStreamingFilter.il.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --output {output.resultsFileccs}

    """
        
        

rule combineCounts:
    input:    
        resultsFileccs="{WD}/{r_type}/DUP_call.queries.counts.ccs",
        resultsFilehg="{WD}/{r_type}/DUP_call.queries.counts.hg",
        resultsFileil="{WD}/{r_type}/DUP_call.queries.counts.il",
    output:
        combineResult="{WD}/{r_type}/DUP_call.queries.counts.ALL",
    shell:"""
                    
paste <( awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$7,$5,$6,$8,$9}}' {input.resultsFilehg} ) <(cut -f 5,6,8,9 {input.resultsFileil}  ) <(cut -f 5,6,8,9 {input.resultsFileccs}   > {output.combineResult}

"""


rule VCFtonucfreq_filter:
    input:
        counts="{WD}/{r_type}/DUP_call.queries.counts.ALL",
      #  freqout=expand("{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",WD=WD,loci=listl),
    output:
        preNucfreq="{WD}/{loci}/{r_type}/kmer_filter/pre.assembly.nucfreq"
    params:
        wd=WD,
        lo="{loci}",
        sd=SD,
        fs=fss,
        sb=sbs,
        rt=rtype,
    shell:"""
mkdir -p {params.wd}/{params.lo}/kmer_filter
intersectBed -a <(grep -v "#" {input.counts} |awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}}'|python $sum/rstrip.py ) -b {params.wd}/{params.lo}/loci.bed > {output.preNucfreq}

    echo preassemblydone
    """












