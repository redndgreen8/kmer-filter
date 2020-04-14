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
samples = ["HG00514"]

#filtering criteria
min_kmer_count = 10
max_kmer_count = 24
min_ill_kmer_count = 5
max_ill_kmer_count = 120

listl=[]
i=0
with open("/panfs/qcb-panasas/rdagnew/HG00514/clr.pbmm/dup_list.HG00514.txt") as l:
    for r in l.readlines():
        if i<2:
            r=r.rstrip()
            listl.append(r)
        i=i+1
#listl=listl[0]

rule all:
    input:
       # done=expand("{WD}/{loci}/sda_out/sda.done",loci=listl ,WD=WD),
        #bamr=expand("{WD}/{loci}/reads.orig.bam",loci=listl ,WD=WD),
        freqout=expand("{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",loci=listl ,WD=WD),
        sda=expand("{WD}/{loci}/sda_out/snvs/sda.assembly.consensus.nucfreq",loci=listl ,WD=WD),
        preNucfreq=expand("{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq",loci=listl ,WD=WD),
        annotatedVcf=expand("{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq.filtered",loci=listl ,WD=WD),
        counts=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.counts",loci=listl ,WD=WD),
        filteredCounts1=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class1.filtered",loci=listl ,WD=WD),
        filteredCounts2=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class2.filtered",loci=listl ,WD=WD),
        comb=expand("{WD}/DUP_call.full.queries.counts",loci=listl ,WD=WD),
        combvcf=expand("{WD}/DUP_call.queries.counts.vcf",loci=listl ,WD=WD),
        resultsFileil=expand("{WD}/DUP_call.queries.counts.il",loci=listl ,WD=WD),
        resultsFilehg=expand("{WD}/DUP_call.queries.counts.hg",loci=listl ,WD=WD),
        queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl ,WD=WD),
        corseq_vcf=expand("{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq",loci=listl ,WD=WD),
        vcff=expand("{WD}/{loci}/kmer_filter/DUP_call.vcf",loci=listl ,WD=WD),
        
                


#<(awk '{{print $1":"$2"-"$3}}' {input.region} )



rule bamToFreq:
    input:
        bam="/panfs/qcb-panasas/rdagnew/HG00514.clr.pbmm.bam",
    output:
        freqout="{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",
    params:
        wd=WD,
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
        wd=WD,
        lo="{loci}",
        sd=SD,
    shell:"""
cat {input.freqout} | $py34 {params.sd}/bamfreq.StreamingKmerFilter.py --nucfreq /dev/stdin --ref $hg38/hg38.fa > {output.vcff}

    """


rule GenerateQueriesFasta:
    input:
        vcff="{WD}/{loci}/kmer_filter/DUP_call.vcf",
    output:
        queryFasta="{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",
        corseq_vcf="{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq",
    params:
        sd=SD,
        #k=config["k"],
    shell:"""
python3 {params.sd}/SNPListToGenotypableQueries.snp.py --snp {input.vcff} --ref $hg38/hg38.fa --queries {output.queryFasta}
touch {output.corseq_vcf}
echo snplisttogenotypedone

"""

rule GenerateResultsFromQuerieshg:
    input:
        queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD),
    output:
        resultsFilehg="{WD}/DUP_call.queries.counts.hg",
    threads:16
    params:
        hgjf='/staging/mjc/rdagnew/hg38.no_alts.32.jf',
        sd=SD,
       # k=config["k"],
    shell:"""
python3 {params.sd}/SVStreamingFilter.hg.py --bpseqs <(cat {input.queryFasta}) --hgjf {params.hgjf} --output {output.resultsFilehg}
"""


rule GenerateResultsFromQueriesil:
    input:
        queryFasta=expand("{WD}/{loci}/kmer_filter/DUP_call.queries.fasta",loci=listl,WD=WD),
    output:
        resultsFileil="{WD}/DUP_call.queries.counts.il",
    threads:16
    params:
        jf="/staging/mjc/rdagnew/HG00514.IL.32.mer_counts.jf",
        sd=SD,
     #   k=config["k"],
    shell:"""
python3 {params.sd}/SVStreamingFilter.il.py --bpseqs <(cat {input.queryFasta}) --jf {params.jf} --output {output.resultsFileil}
"""


rule combineQuerycounts:
    input:
        resultsFile="{WD}/DUP_call.queries.counts.hg",
        resultsFile2="{WD}/DUP_call.queries.counts.il",
        corseq=expand("{WD}/{loci}/kmer_filter/DUP_call.vcf",loci=listl,WD=WD),
        vcf_header="{WD}/vcf_header.txt",

    output:
        comb="{WD}/DUP_call.full.queries.counts",
        combvcf="{WD}/DUP_call.queries.counts.vcf",

    shell:"""
paste <(cat {input.resultsFile}) <(cut -f 3-4 {input.resultsFile2}) > {output.comb}

#cat {input.vcf_header} <(paste <( cat {input.corseq}|grep -v "#") <(cat {output.comb} |tr "\t" ";")) > {output.combvcf}    
#tail {output.combvcf}

    """















rule FilterQueryResults:
    input:
        combvcf="{WD}/DUP_call.queries.counts.vcf",

    output:
        counts="{WD}/{loci}/kmer_filter/DUP_call.queries.counts",
        filteredCounts1="{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class1.filtered",
        filteredCounts2="{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class2.filtered",
        #report="{sample}.filtered.report"
    params:
        min_count = min_kmer_count,
        max_count = max_kmer_count,
        min_ill_count = min_ill_kmer_count,
        max_ill_count = max_ill_kmer_count,
        lo="{loci}",
        sd=SD,
        wd=WD,

    shell:"""
    python {params.wd}/partition.full_vcf.py --vcf {input.combvcf} --region {params.lo} 


    python3 {params.sd}/GenotypableQueriesToFilterdQueries.py --counts {output.counts} \
                                              --mincount {params.min_count} \
                                              --maxcount {params.max_count} \
                                              --minill {params.min_ill_count} \
                                              --maxill {params.max_ill_count} \
                                              --mean_cov 95
    touch {output.filteredCounts1} {output.filteredCounts2}

"""

        
rule WriteFilteredVCF:
    input:
        unannotatedVcf="{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq",
        class_1_counts="{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class1.filtered",
        class_2_counts="{WD}/{loci}/kmer_filter/DUP_call.queries.counts.class2.filtered",
    output:
        annotatedVcf="{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq.filtered",
    params:
        sd=SD,
       # k=config["k"],
    shell:"""
python3 {params.sd}/ClassCountsToFilteredSVList.py --sv {input.unannotatedVcf} --c1count {input.class_1_counts} --c2count {input.class_2_counts}
touch {output.annotatedVcf}
echo correctedseqdone

"""






rule VCFtonucfreq:
    input:
        annotatedVcf="{WD}/{loci}/kmer_filter/DUP_call.vcf.correctedseq.filtered",
        freqout="{WD}/{loci}/kmer_filter/bamToFreq.DUP_CALL.out",
    output:
        preNucfreq="{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq",
    params:
        wd=WD,
        lo="{loci}",
        sd=SD,
    shell:"""

    intersectBed -wb -a <(awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9,$10}}' {input.freqout}) -b <(grep -v "#" {input.annotatedVcf}|awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$2+1,$4,$5,$8}}'| tr -d "." ) | \
                                            python {params.sd}/pre_sda.nucfreq_filter.py \
                                            -sc 6 -fb 115 -sb 102 -med 10 -m 95 --output {params.wd}/{params.lo}/kmer_filter/filtered.orig_counts.nucfreq > {output.preNucfreq}
        echo preassemblydone
    """


rule nucfreqToSDAformat:
    input:
        preNucfreq="{WD}/{loci}/kmer_filter/mean_median_cov.filter.pre.assembly.consensus.nucfreq",
        loci_bed="{WD}/{loci}/loci.bed",
    output:
        sda="{WD}/{loci}/sda_out/snvs/sda.assembly.consensus.nucfreq",
    params:
        wd=WD,
        lo="{loci}",
    shell:"""
paste <(cut -f 1-9 {input.preNucfreq}) <(cat {input.loci_bed} )|awk 'BEGIN{{OFS="\t";b=0;c=0}} {{if(NR==1) {{b=$11;c=$12; print $1":"b"-"c,$2-b+1,$4,$5,$6,$7,$8,$9 }} else {{print $1":"b"-"c,$2-b+1,$4,$5,$6,$7,$8,$9}}  }}'  >{output.sda}
    echo sdaconsensusdone

    """

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










