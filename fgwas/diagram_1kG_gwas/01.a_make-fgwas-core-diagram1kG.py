#!/usr/bin/python -O
# Jason Matthew Torres
'''
Prepare fgwas core file
Usage: python 01_make-fgwas-core-diagram1kG.py
'''
# libraries
import sys,os,gzip
import subprocess as sp


# globals
home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/"
out_dir=home_dir+"intermediate_files/"
#out_file = out_dir + "diagram1kG-fgwas-core.txt.gz"
in_dir = "/well/got2d/jason/reference/gwas/diagram_1Kgenomes/"
input_file=in_dir+"DIAGRAM_T2D_1000G.2017.gz" # "1ktest.txt.gz" #
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)

mychrom = int(sys.argv[1])

def make_core(chromosome,threshold=0.01):
    '''
    Read input file and write to output file with proper format
    SNPID CHR POS F Z NCASE NCONTROL SE
    '''
    fin = gzip.open(input_file,'rb')
    out_file = out_dir + "chr"+str(chromosome)+".fgwas-core.txt.gz"
    fout = gzip.open(out_file,'wb')
    head_list = ["SNPID","CHR","POS","F","Z","BETA","SE","PVAL","NCASE","NCONTROL"]
    fout.write("\t".join(head_list)+"\n")
    fin.readline() # header
    count = 0
    for line in fin:
        #count+=1
        #sys.stdout.write("\r%d" % count)
        #sys.stdout.flush()
        l = line.strip().split()
        snpid = "chr"+l[0]
        f = l[3]
        chrom = snpid.split(":")[0]
        pos = snpid.split(":")[1]
        ncase, ncont = "NA","NA"
        beta, se, pval = l[7],l[8], l[9]
        z = str(float(beta)/float(se))
        maf = l[16]
        if float(f) >= threshold and float(f) <= (1-threshold):
            if chrom == ("chr"+str(chromosome)):
                write_list = [snpid,chrom,pos,f,z,beta,se,pval,ncase,ncont]
                fout.write("\t".join(write_list)+"\n")
    fout.close()
    fin.close()
    sys.stdout.write("\n")

def main():
    #for c in range(1,23):
        #print ("\nChromsome: %d\n" % c)
    make_core(mychrom,0.01)

if (__name__=="__main__"): main()
