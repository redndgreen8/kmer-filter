import argparse
import math
import statistics

ap = argparse.ArgumentParser(description="Parse the queries list/counts for support.")
ap.add_argument("--counts", help="File with counts of ref/alt sequence in HG38 and Illumina.", required=True)
# ap.add_argument("--report", help="Output report file to write to", required=True)
ap.add_argument("--mincount", help="Minimum number of times a k-mer appears to be kept (class1/2).", default=10)
ap.add_argument("--maxcount", help="Maximum number of times a k-mer appears to be kept (class2).", default=24)
ap.add_argument("--minill", help="Minimum number of times a k-mer appears to be kept for Illumina (class1/2).", default=5)
ap.add_argument("--maxill", help="Maximum number of times a k-mer appears to be kept for Illumina (class1/2).", default=120)
ap.add_argument("--mean_cov", help="Mean sequencing coverage", required=True)
ap.add_argument("--f", help="filter 0/1",default=0)
ap.add_argument("--fraction", help="fraction", type=float,default=0.25)

args=ap.parse_args()

counts_filepath = args.counts
min_counts = int(args.mincount)
max_counts = int(args.maxcount)
min_ill_given = int(args.minill)
max_ill_given = int(args.maxill)

mean = float(args.mean_cov)


class_1_out = counts_filepath+'.class1.filtered'
class_2_out = counts_filepath+'.class2.filtered'

report_base_str = "\n{:13}{:5}{:4}{:>6}{:>7}{:>10}{:>10}{:>10}{:>9}{:>8}{:>9}{:>9}{:>10}{:>8}"

def _calculate_support(alt_list,min_ill,max_ill):
    return min_ill <= sum([1 for count in alt_list if count > 0]) <= max_ill
    #change max_il

c0 = {'lines':[],
     }
c1 = {'lines':[],
      'mean_r_len':[],
      'mean_r_support':[],
      'mean_a_len':[],
      'mean_a_support':[],
     }
c2 = {'lines':[],
      'mean_r_len':[],
      'mean_r_support':[],
      'mean_a_len':[],
      'mean_a_support':[],
     }








with open(counts_filepath) as f:
    for line in f:
        ref_seq, alt_seq, ref_hg_list, ref_ill_list, alt_hg_list, alt_ill_list = line.split()
        
        alt_ill_list = [int(v) for v in alt_ill_list.split(',')]
        len_alt_list=len(alt_ill_list)
        
   
#//edits
        altsum=0
        for i in range(len(alt_ill_list )):
            altsum=altsum+alt_ill_list[i]
        mean_alt_support=altsum/len(alt_ill_list)
       # print(str(altsum) +"\t"+str(len(alt_ill_list)) +"\t"+ str(mean_alt_support) )
        if args.f is not None:
            if mean_alt_support < max(4,(args.fraction*mean)) or mean_alt_support > ((args.fraction+1)*mean):
                continue
        #median filter     #0-based!!!
        #sorted_alt_list=sorted(alt_ill_list)
        med=0
        #if len_alt_list%2!=0:
         #   med=sorted_alt_list[floor( len_alt_list/2 )]
        #else:
         #   med=(sorted_alt_list[int(len_alt_list/2)] + sorted_alt_list[int(len_alt_list/2)-1])/2
        med2=statistics.median(alt_ill_list)
        
#edits//



        alt_hg_list = [int(v) for v in alt_hg_list.split(',')]
        alt_hg_sup = [1 for count in alt_hg_list if count == 0]
       


        alt_hg_support = sum(alt_hg_sup) >= min_counts
        



        line=line.rstrip()+"\t"+str(med2)+"\n"
        
        #class0/all calls logic
        c0['lines'].append(line)
        
        if _calculate_support(alt_ill_list, min_ill_given, max_ill_given) and alt_hg_support:
            ref_hg_list = [int(v) for v in ref_hg_list.split(',')]
            ref_hg_sup1 = [1 for count in ref_hg_list if count == 1]
            ref_hg_sup2 = [1 for count in ref_hg_list if 1 < count < max_counts]
            
            #class 1
            if sum(ref_hg_sup1) >= min_counts:
                    c1['lines'].append(line)
                    c1['mean_r_len'].append(len(ref_hg_list))
                    c1['mean_r_support'].append(len(ref_hg_sup1))
                    c1['mean_a_len'].append(len(alt_hg_list))
                    c1['mean_a_support'].append(len(alt_hg_sup))
            #class 2
            elif sum(ref_hg_sup2) >= min_counts:
                    c2['lines'].append(line)
                    c2['mean_r_len'].append(len(ref_hg_list))
                    c2['mean_r_support'].append(len(ref_hg_sup2))
                    c2['mean_a_len'].append(len(alt_hg_list))
                    c2['mean_a_support'].append(len(alt_hg_sup))
with open(class_1_out,'w') as f:
    f.write('\n'.join(c1['lines']))
with open(class_2_out,'w') as f:
    f.write('\n'.join(c2['lines']))
