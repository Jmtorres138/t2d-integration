#!/usr/bin/python -O
# Jason Matthew Torres
'''
Subset best fgwas model output to snps in significant blocks
Note: First must run "Get Significant Blocks" section of file functional_credible_sets.Rmd
Usage: python subset_to_seg.py pre sig_block_file 0.9 
'''
# libraries
import sys,os,gzip


# globals
server_dir = "/well/got2d/jason/"
work_dir = server_dir + "projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
#pre = work_dir + "fgwas_output/" + "drop-islet_state11+islet_state6+islet_state12+utr_5+islet_stretch_enhancer+promoter+distance_tss+islet_state9"
pre = sys.argv[1]
#sig_block_file = work_dir + "sig_blocks_ppa90.txt"
sig_block_file = sys.argv[2]
cum_thresh = sys.argv[3]
# functions

def build_dic():
    fin = open(sig_block_file,'r')
    header = fin.readline()
    dic = {}
    for line in fin:
        l = line.strip().split()
        seg, chrom, st, sp, ppa = l[0],l[1],l[2],l[3],l[4]
        try:
            dic[chrom].append([int(st),int(sp),seg,ppa])
        except:
            dic[chrom] = [[int(st),int(sp),seg,ppa]]
    fin.close()
    return(dic)

def subset_to_seg():
    print("Building dictionary...")
    ref_dic =  build_dic()
    fin = gzip.open(pre+".bfs.gz",'rb')
    head_list = fin.readline().strip().split()
    head_list = ["SEGNUMBER"] + head_list
    fout = gzip.open(work_dir+"sig_block_snps.bfs"+"_"+str(cum_thresh)+".txt.gz",'wb')
    print("Writing output file...")
    fout.write(" ".join(head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        snpid, chrom, pos = l[0],l[1],int(l[2])
        try:
            for li in ref_dic[chrom]:
                start, stop, seg = li[0],li[1],li[2]
                if pos >= start and pos <= stop:
                    print ("\n" + snpid + "\n")
                    write_list = [seg] + l
                    fout.write(" ".join(write_list)+"\n")
        except:
            pass
    fout.close()
    fin.close()



def main():
    subset_to_seg()

if (__name__=="__main__"): main()
