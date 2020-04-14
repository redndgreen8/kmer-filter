import os
import tempfile
import numpy as np

#
# Locate the tempdir for grid processing
#
#configfile: "SVGenotyping.json"

if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
else:
    TMPDIR = tempfile.gettempdir()

SD ="/home/cmb-16/mjc/rdagnew/summerproj/kmer_filter/"#os.path.dirname(workflow.snakefile)

#
# 
# samples = config["samples"].keys()
samples = ["HG00514"]

#filtering criteria
min_kmer_count = 10
max_kmer_count = 24
min_ill_kmer_count = 5
max_ill_kmer_count = 120

def GetVcfPaths(sample_id):
    paths = glob_wildcards('{filename}.vcf')
    return paths

def GetAsm(sampleID):
    return config["samples"][sampleID]["ASM"]

def GetJF(sampleID):
    return config["samples"][sampleID]["JF"]
    
def _get_lines_in_count_file(filepath,output_stats):
    assert output_stats in ['all','class1','class2'], 'output_stats value not recognized'
    total = 0
    ref_hg_means = []
    alt_hg_means = []
    ref_hg_len = []
    alt_hg_len = []
    with open(filepath) as f:
        for line in f:
            s_line = line.split()
            if len(s_line): #don't count empty lines
                total += 1
                
                ref_hg_len.append(len(s_line[2]))
                alt_hg_len.append(len(s_line[4]))
                if output_stats == 'all':
                    ref_counts = [int(v) for v in s_line[2].split(',')]
                    alt_counts = [int(v) for v in s_line[4].split(',')]
                    ref_hg_means.append(np.mean(ref_counts))
                    alt_hg_means.append(np.mean(alt_counts))
                elif output_stats == 'class1':
                    ref_counts = [1 for v in s_line[2].split(',') if int(v) == 1]
                    alt_counts = [1 for v in s_line[4].split(',') if int(v) == 0]
                    ref_hg_means.append(sum(ref_counts))
                    alt_hg_means.append(sum(alt_counts))
                elif output_stats == 'class2':
                    ref_counts = [1 for v in s_line[2].split(',') if 1 < int(v) <max_kmer_count]
                    alt_counts = [1 for v in s_line[4].split(',') if int(v) == 0]
                    ref_hg_means.append(sum(ref_counts))
                    alt_hg_means.append(sum(alt_counts))
                
    return total,np.mean(ref_hg_len),np.mean(ref_hg_means),np.mean(alt_hg_len),np.mean(alt_hg_means)
    
report_base_str = ""
#ops=["DEL","INS"]
#classes=["tr", "non"]
#class_types=[1,2]

rule all:
    input:
        queryFasta=expand("svcall.{sample}.{op}.{cl}.queries.fasta", sample=samples, op=ops, cl=classes),
        queryResults=expand("svcall.{sample}.{op}.{cl}.queries.counts", sample=samples, op=ops, cl=classes),
        filteredQueryResults=expand("svcall.{sample}.{op}.{cl}.queries.counts.class1.filtered", sample=samples, op=ops, cl=classes),
        sampleFilteredReport=expand("{sample}.filtered.report", sample=samples,),
        combinedFilteredReport="all_samples.filtered.report",
        vcfFilteredByClasses=expand("svcall.{sample}.{op}.{cl}.vcf.correctedseq.filtered", sample=samples, op=ops, cl=classes),
        
rule WriteFilteredVCF:
    input:
        unannotatedVcf="svcall.{sample}.{op}.{cl}.vcf.correctedseq",
        class_1_counts="svcall.{sample}.{op}.{cl}.queries.counts.class1.filtered",
        class_2_counts="svcall.{sample}.{op}.{cl}.queries.counts.class2.filtered",
    output:
        annotatedVcf="svcall.{sample}.{op}.{cl}.vcf.correctedseq.filtered"
    params:
        ref=config["ref"],
        asm=lambda wildcards: GetAsm(wildcards.sample),
        sd=SD,
        k=config["k"],
    shell:"""
python3 {params.sd}/ClassCountsToFilteredSVList.py --sv {input.unannotatedVcf} --c1count {input.class_1_counts} --c2count {input.class_2_counts}
"""
    
rule GenerateCombinedReports:
    input:
        sample_list = ["{}.filtered.report".format(sample) for sample in samples]
    output:
        combined_report = "all_samples.filtered.report",
        # tr_report = "tr.filtered.report",
        # non_report = "non.filtered.report",
    run:
        out_lines = []
        tr_lines = []
        non_lines = []
        
        header = report_base_str.format('Sample Name','Op','Cl.',
                                        'Calls','Cl.1','#r kmers','#r=1','#a kmers','#a=0',
                                        'Cl.2','#r kmers','#r<{}'.format(max_kmer_count),'#a kmers','#a=0') + '\n'
        out_lines.append(header)
        tr_lines.append(header)
        non_lines.append(header)
        
        for file in input.sample_list:
            with open(file) as f:
                print(file)
                for line in f:
                    out_lines.append(line)
                    if line == '\n':
                        continue
                    
                    if line.split()[2] == 'tr':
                        tr_lines.append(line)
                    else:
                        non_lines.append(line)
        with open(output.combined_report,'w') as f:
            for line in out_lines:
                f.write(line)
        # with open(output.tr_report,'w') as f:
            # for line in tr_lines:
                # f.write(line)
        # with open(output.non_report,'w') as f:
            # for line in non_lines:
                # f.write(line)
                
rule GenerateSampleReport:
    input:
        queryResults = ["svcall.{sample}"+".{}.{}.queries.counts".format(op,cl) for op in ops for cl in classes] \
                      +["svcall.{sample}"+".{}.{}.queries.counts".format(op,cl)+".class1.filtered" for op in ops for cl in classes] \
                      +["svcall.{sample}"+".{}.{}.queries.counts".format(op,cl)+".class2.filtered" for op in ops for cl in classes],
    output:
        report = "{sample}.filtered.report"
    run:
        out_lines = []
        sample = wildcards.sample
        for op in ops:
            for cl in classes:
                all_counts = "svcall.{}.{}.{}.queries.counts".format(sample,op,cl)
                class1_counts = all_counts + ".class1.filtered"
                class2_counts = all_counts + ".class2.filtered"
                
                unfiltered_count, uf_mean_ref_len, uf_ref_mean, uf_mean_alt_len,uf_alt_mean = _get_lines_in_count_file(all_counts,'all')
                filtered_class1, c1_mean_ref_len, c1_ref_mean, c1_mean_alt_len, c1_alt_mean = _get_lines_in_count_file(class1_counts,'class1')
                filtered_class2, c2_mean_ref_len, c2_ref_mean, c2_mean_alt_len, c2_alt_mean = _get_lines_in_count_file(class2_counts,'class2')
                
                out_lines.append(report_base_str.format(sample,op,cl,
                                                        '{:d}'.format(unfiltered_count),'{:d}'.format(filtered_class1),'{:.2f}'.format(c1_mean_ref_len),
                                                        '{:.2f}'.format(c1_ref_mean),'{:.2f}'.format(c1_mean_alt_len),'{:.2f}'.format(c1_alt_mean),
                                                        '{:d}'.format(filtered_class2),'{:.2f}'.format(c2_mean_ref_len),'{:.2f}'.format(c2_ref_mean),
                                                        '{:.2f}'.format(c2_mean_alt_len),'{:.2f}'.format(c2_alt_mean))
                                                        )
                                            
        with open(output.report,'w') as f:
            for line in out_lines:
                f.write(line)
                
rule FilterQueryResults:
    input:
        counts="svcall.{sample}.{op}.{cl}.queries.counts"
    output:
        filteredCounts1="svcall.{sample}.{op}.{cl}.queries.counts.class1.filtered",
        filteredCounts2="svcall.{sample}.{op}.{cl}.queries.counts.class2.filtered",
        #report="{sample}.filtered.report"
    params:
        min_count = min_kmer_count,
        max_count = max_kmer_count,
        min_ill_count = min_ill_kmer_count,
        max_ill_count = max_ill_kmer_count,
    shell:"""python3 GenotypableQueriesToFilterdQueries.py --counts {input.counts} \
                                              --mincount {params.min_count} \
                                              --maxcount {params.max_count} \
                                              --minill {params.min_ill_count} \
                                              --maxill {params.max_ill_count} \
"""

rule GenerateResultsFromQueries:
    input:
        queryFasta="svcall.{sample}.{op}.{cl}.queries.fasta"
    output:
        resultsFile="svcall.{sample}.{op}.{cl}.queries.counts"
    params:
        hgjf='/home/cmb-16/mjc/shared/references/hg38_noalts/indexes/jf/hg38.22.jf',
        jf=lambda wildcards: GetJF(wildcards.sample),
        sd=SD,
        k=config["k"],
    shell:"""
python3 {params.sd}/SVStreamingFilter.py --bpseqs {input.queryFasta} --jf {params.jf} --hgjf {params.hgjf} --output {output.resultsFile}
"""

rule GenerateQueriesFasta:
    input:
        unannotatedVcf="svcall.{sample}.{op}.{cl}.vcf"
    output:
        queryFasta="svcall.{sample}.{op}.{cl}.queries.fasta",
        corseq_vcf="svcall.{sample}.{op}.{cl}.vcf.correctedseq"
    params:
        ref=config["ref"],
        asm=lambda wildcards: GetAsm(wildcards.sample),
        sd=SD,
        k=config["k"],
    shell:"""
python3 {params.sd}/SVListToGenotypableQueries.py --sv {input.unannotatedVcf} --ref {params.ref}/indexes/bwa/hg38.fa --kmer {params.k} --queries {output.queryFasta} --asm {params.asm}
"""

# FIRST part to complete:
rule AnnotateGenotypingKmers:
    input:
        unannotatedVcf="svcall.{sample}.vcf"
        
    output:
        kmerAnnotatedVcf="svcall.{sample}.annotated.vcf"
    params:
        ref=config["ref"]
    shell:"""
.... add commands (probaby some script that combines /home/cmb-16/mjc/mchaisso/projects/SVKmerGenotyping/SVListToGenotypableKmers.py and /home/cmb-16/mjc/mchaisso/projects/SegDupDiversity/code/SegDupSNV/StreamingKmerFilter.py
"""
