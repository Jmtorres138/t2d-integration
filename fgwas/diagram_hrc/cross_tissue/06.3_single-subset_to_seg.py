#!/usr/bin/python -O
# Jason Matthew Torres
'''
Subset best fgwas model output to snps in significant blocks
Note: First must run "Get Significant Blocks" section of file functional_credible_sets.Rmd
Usage: python subset_to_seg.py pre block_file
'''
# libraries
import sys,os,gzip


# globals
server_dir = "/well/got2d/jason/"
work_dir = server_dir + "projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/"

# functions

def build_dic(sig_block_file):
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

def subset_to_seg(sig_block_file,pre,out_dir):
    print("Building dictionary...")
    ref_dic =  build_dic(sig_block_file)
    fin = gzip.open(pre+".bfs.gz",'rb')
    head_list = fin.readline().strip().split()
    head_list = ["SEGNUMBER"] + head_list
    fout = gzip.open(out_dir+"loci_block_snps.bfs.txt.gz",'wb')
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
    tiss_list = ["adipose","islet","liver","muscle"]
    for tiss in tiss_list:
        pre = work_dir + "fgwas_output_" + tiss + "/" + "SINGLE_fgwas_run_loci-partition"
        sig_block_file = work_dir + "fgwas_output_" + tiss + "/" + "fgwas_blk-t2d_loci_151.txt"
        out_dir = work_dir + "fgwas_output_" + tiss + "/"
        subset_to_seg(sig_block_file,pre,out_dir)

if (__name__=="__main__"): main()
