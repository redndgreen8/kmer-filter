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

SD ="/home/cmb-16/mjc/rdagnew/summerproj/kmer_filter"#os.path.dirname(workflow.snakefile)
WD="/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm"
SDA="/home/cmb-16/mjc/rdagnew/summerproj/SDA/SDA"
#
# 
# samples = config["samples"].keys()
#samples = ["HG00514"]

#filtering criteria
min_kmer_count = 5
max_kmer_count = 34
min_ill_kmer_count = 5
max_ill_kmer_count = 120

listl=[]
#listfull=[]
i=0
with open("/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/dup_list.HG00514.txt") as l:
    for r in l.readlines():
        if i<250 and i>=0:
            r=r.rstrip()
            listl.append(r)
        i=i+1
#parts=range(0,len(listfull),10)

#when input gets too large 
# no. queries/psv sites * k  >50 million
# no. queries >1.5 million
# no. of loci >20

part=0 
if len(listl)>20:
    part=1
    print("partition il and hg38")


#def partitionss(listl):
rule all:
    input:
        freqout=expand("{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out", WD=WD,loci=listl),
        vcff=expand("{WD}/{loci}/kmer_filter/DUP_call.vcf",loci=listl ,WD=WD),
        queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta", loci=listl, WD=WD),
        preNucfreqdone=expand("{WD}/pre.assembly.consensus.nucfreq.done",loci=listl, WD=WD),

        resultsFile=expand("{WD}/DUP_call.queries.counts",loci=listl ,WD=WD),
       # resultsFileil=expand("{WD}/DUP_call.queries.counts.il",loci=listl ,WD=WD),
       # resultsFilehg=expand("{WD}/DUP_call.queries.counts.hg",loci=listl ,WD=WD),
        counts=expand("{WD}/DUP_call.queries.counts.filtered.vcf",loci=listl ,WD=WD),

        preNucfreq=expand("{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq",loci=listl ,WD=WD),
        loci_bed=expand("{WD}/{loci}/loci.bed",loci=listl ,WD=WD),
        sda=expand("{WD}/{loci}/sda_out/snvs/sda.assembly.consensus.nucfreq",loci=listl ,WD=WD),


rule bamToFreq:
    input:
        bam="/panfs/qcb-panasas/rdagnew/HG00514.clr.pbmm.bam",
    output:
        freqout="{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",
    params:
        lo="{loci}",
    shell:"""
$sum/SegDupSNV/bamToFreq {input.bam} <(echo {params.lo}|tr "_" "\t"|cut -f 1) $hg38/hg38.fa > {output.freqout}
    
    echo bamfreqdone

    """

rule vcf:
    input:
        freqout="{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",
    output:
        vcff="{WD}/{loci}/kmer_filter/DUP_call.vcf",
    params:
        sd=SD,
    shell:"""
cat {input.freqout} | $py34 {params.sd}/bamfreq.StreamingKmerFilter.py --nucfreq /dev/stdin --ref $hg38/hg38.fa > {output.vcff}

    """


rule GenerateQueriesFasta:
    input:
        vcff="{WD}/{loci}/kmer_filter/DUP_call.vcf",
    output:
        queryFasta="{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",
    params:
        sd=SD,
    shell:"""
python3 {params.sd}/SNPListToGenotypableQueries.snp.py --snp {input.vcff} --ref $hg38/hg38.fa --queries {output.queryFasta}
echo snplisttogenotypedone

"""
#when input gets too large 
# no. queries/psv sites * k  >50 million
#no. queries >1.5 million
if part:
    rule GenerateResultsFromQuerieshg:
        input:
            queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD),
        output:
            resultsFilehg="{WD}/DUP_call.queries.counts.hg",
        threads:16
        params:
            hgjf='/panfs/qcb-panasas/rdagnew/hg38.no_alts.32.jf',
            sd=SD,
        shell:"""
        python3 {params.sd}/SVStreamingFilter.hg.py --bpseqs <(cat {input.queryFasta}) --hgjf {params.hgjf} --output {output.resultsFilehg}
        """
    rule GenerateResultsFromQueriesil:
        input:
            queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD),
            resultsFilehg="{WD}/DUP_call.queries.counts.hg",
        output:
            resultsFileil="{WD}/DUP_call.queries.counts.il",
            resultsFile="{WD}/DUP_call.queries.counts",
        threads:16
        params:
            jf="/panfs/qcb-panasas/rdagnew/HG00514.IL.32.mer_counts.jf",
            sd=SD,
        shell:"""
        python3 {params.sd}/SVStreamingFilter.il.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --output {output.resultsFileil}
   
        paste <(cat {input.resultsFilehg} |cut -f 1-6) <(cat {output.resultsFileil}|cut -f 5-7) > {output.resultsFile}
        
        wc -l {input.resultsFilehg} {output.resultsFileil} {output.resultsFile}
        """
else:


    rule GenerateResultsFromQueries:
        input:
            queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD),
        output:
            resultsFile="{WD}/DUP_call.queries.counts",
        threads:16
        params:
            jf="/panfs/qcb-panasas/rdagnew/HG00514.IL.32.mer_counts.jf",
            sd=SD,
            hgjf='/panfs/qcb-panasas/rdagnew/hg38.no_alts.32.jf',
        shell:"""
        python3 {params.sd}/SVStreamingFilter.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --hgjf {params.hgjf} --output {output.resultsFile}

        """


rule FilterQueryResults:
    input:
        vcf_header="{WD}/vcf_header.txt",
        resultsFile2="{WD}/DUP_call.queries.counts",
    output:
        counts="{WD}/DUP_call.queries.counts.filtered.vcf",
    params:
        min_count = min_kmer_count,
        max_count = max_kmer_count,
        min_ill_count = min_ill_kmer_count,
        max_ill_count = max_ill_kmer_count,
        #lo="{loci}",
        sd=SD,
    shell:"""

    python3 {params.sd}/GenotypableQueriesToFilterdQueries.NEW.py --counts {input.resultsFile2} \
                                              --mincount {params.min_count} \
                                              --maxcount {params.max_count} \
                                              --minill {params.min_ill_count} \
                                              --maxill {params.max_ill_count} \
                                              --mean_cov 95 --vcf {input.vcf_header}
    touch {output.counts} 

    """



rule VCFtonucfreq_filter:
    input:
        counts="{WD}/DUP_call.queries.counts.filtered.vcf",
      #  freqout=expand("{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",WD=WD,loci=listl),
    output:
        preNucfreq="{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq"
    params:
        wd=WD,
        lo="{loci}",
        sd=SD,
    shell:"""
intersectBed -a <(grep -v "#" {input.counts} |awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$2+1,$4,$5,$8}}' ) -b {params.wd}/{params.lo}/loci.bed |     python {params.sd}/pre_sda.filter.py --nucfreq /dev/stdin/ --region {params.lo} -fb 112 -sb 105  > {output.preNucfreq}

    echo preassemblydone
    """


rule filterNucfreqToSDAformat:
    input:
        #preNucfreqdone="{WD}/pre.assembly.consensus.nucfreq.done",
        preNucfreq="{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq",
        #loci_bed="{WD}/{loci}/loci.bed",
    output:
        sda="{WD}/{loci}/sda_out/snvs/sda.assembly.consensus.nucfreq",
    params:
        wd=WD,
        lo="{loci}",
    shell:"""
paste <(cut -f 1-9 {params.wd}/{params.lo}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq) <(cat {params.wd}/{params.lo}/loci.bed )|awk 'BEGIN{{OFS="\t";b=0;c=0}} {{if(NR==1) {{b=$11;c=$12; print $1":"b"-"c,$2-b+1,$4,$5,$6,$7,$8,$9 }} else {{print $1":"b"-"c,$2-b+1,$4,$5,$6,$7,$8,$9}}  }}'  >{output.sda}
    echo sdaconsensusdone

    """

#def main():
#    for i in range(len(parts)):
#        if i ==len(parts)-1:
#            continue
#        listl=listfull[parts[i]:parts[i+1]]
#        print(listl,WD,SD,SDA)
#        partitionss(listl)


#if __name__ == '__main__':
#    main()



#rule SDA:
  #  input:
 #       bamr="{WD}/{loci}/reads.orig.bam",
#        sda="{WD}/{loci}/sda_out/snvs/sda.assembly.consensus.nucfreq",

   # output:
   #     done="{WD}/{loci}/sda_out/sda.done",
   # threads:8
   # params:
   #     wd=WD,
  #      lo="{loci}",
 #       sdap=SDA,
#shell:"""
##echo enteringsda
#cd {params.wd}/{params.lo}
#touch {input.sda}

#{params.sdap} collapse --coverage 95 --reads {input.bamr} --ref ref.fasta -t 8 --debug
#echo sdadone
#touch {output.done}

#"""
     #   done=expand("{WD}/{loci}/sda.done",loci=listl ,WD=WD),



#rule plotSDAout:
#shell:"""
#if [ -s $pan/HG00514/$t.$g/$r/sda_out/CC/sda.mi.gml.sites ];then
#                echo $r
#                echo SDApass
#                cat $pan/HG00514/$t.$g/$r/sda_out/CC/sda.mi.gml.sites|tr "\t" "\n" |sort -k1,1n>$pan/HG00514/$t.$g/$r/kmer_filter/retained_PSV.txt;
#                python $sum/kmer_filter/plot_psv_retained.py --filtered $pan/HG00514/$t.$g/$r/kmer_filter/filtered.orig_counts.nucfreq --retained $pan/HG00514/$t.$g/$r/kmer_filter/retained_PSV.txt --bamFreqOut $pan/HG00514/$t.$g/$r/kmer_filter/$r.local.bamToFreq.DUP_CALL.out --mean 95 --prefix $pan/HG00514/$t.$g/$r/sda_out/$r.final.
#            else
#                echo $r
#                echo SDAfail
#                python $sum/kmer_filter/plot_psv_retained.py --filtered $pan/HG00514/$t.$g/$r/kmer_filter/filtered.orig_counts.nucfreq --bamFreqOut $pan/HG00514/$t.$g/$r/kmer_filter/$r.local.bamToFreq.DUP_CALL.out --mean 95 --prefix $pan/HG00514/$t.$g/$r/sda_out/$r.final.
#            fi
#"""










