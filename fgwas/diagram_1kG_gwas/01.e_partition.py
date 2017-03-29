#!/usr/bin/python -O
# Jason Matthew Torres
'''
Prepare fgwas core file
Usage: python 01.e_partition.py
'''
# libraries
import sys,os,gzip
import subprocess as sp


# globals
home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/"
out_dir=home_dir+"intermediate_files/"
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)

def part(mychrom, part=20):
    '''
    Partition a core file into smaller files
    '''
    fname = out_dir + "chr"+str(mychrom) + ".fgwas-core.txt.gz"
    fin = gzip.open(fname,'rb')
    head = fin.readline()
    nrows = 0
    for line in fin:
        nrows += 1
    fin.close()
    num = nrows/part
    rem = nrows%part
    fin = gzip.open(fname,'rb')
    head = fin.readline().strip()
    track=1
    outfile = out_dir + "chr"+str(mychrom) + ".p"+str(track)+".fgwas-core.txt.gz"
    fout = gzip.open(outfile,'wb')
    fout.write(head+"\n")
    count=0
    for line in fin:
        count+=1
        if track != part and count < num:
            fout.write(line)
        elif track==part:
            fout.write(line)
        elif count >= num:
            count=0
            track+=1
            fout.close()
            outfile = out_dir + "chr"+str(mychrom) + ".p"+str(track)+".fgwas-core.txt.gz"
            fout = gzip.open(outfile,'wb')
            fout.write(head+"\n")
            fout.write(line)
        else:
            print "Error: Please inspect"
    fout.close()
    fin.close()

def main():
    for c in range(1,23):
        print "Chromosome: %d" % c
        part(c)

if (__name__=="__main__"): main()
