#!/usr/bin/python -O
# Jason Matthew Torres
'''
Subset best fgwas model output to snps in significant blocks
Note: First must run "Get Significant Blocks" section of file functional_credible_sets.Rmd
Usage: python subset_to_seg.py
'''
# libraries
import sys,os,gzip


# globals
server_dir = "/well/got2d/jason/"
work_dir = server_dir + "projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/"
cred_dir = work_dir + "credible_sets/"
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

def subset_to_seg(pre,sig_block_file,out_file):
    print("Building dictionary...")
    ref_dic =  build_dic(sig_block_file)
    fin = gzip.open(pre+".bfs.gz",'rb')
    head_list = fin.readline().strip().split()
    head_list = ["SEGNUMBER"] + head_list
    fout = gzip.open(out_file,'wb')
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

def run_all(roundd,ppthr="0.9",tiss_list=["islet","liver","adipose","muscle"]):
    print "Round: " + str(roundd)
    for tiss in tiss_list:
        print tiss
        if roundd == 1:
            pre = work_dir + "fgwas_output_" + tiss +"/" + "best-joint-model"
        elif roundd ==2:
            pre = work_dir + "fgwas_output_" + tiss +"/round2/" + "best-joint-model"
        else:
            raise TypeError("Must enter valid round value: 1 or 2")
        sig_block_file = cred_dir + "Round"+str(roundd)+"."+tiss+".sig_blocks_ppa"+str(ppthr)+".txt"
        out_file = cred_dir + "Round"+str(roundd)+"."+tiss+".loci_block_snps.bfs.txt.gz"
        subset_to_seg(pre,sig_block_file,out_file)


def main():
    #run_all(roundd=1)
    run_all(roundd=2,tiss_list=["liver","adipose","muscle"])

if (__name__=="__main__"): main()
