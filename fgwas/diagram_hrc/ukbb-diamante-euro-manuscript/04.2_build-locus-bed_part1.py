#!/usr/bin/python -O
# Jason Matthew Torres
'''
Make file of chromosome position boundries (needed for constructing locus delimited bed file for subsequent fgwas analysis)
Usage: python 05.2_build-locus-bed_part1.py
'''
# libraries
import sys,os,gzip


# globals
server_dir = "/well/got2d/jason/"
work_dir = server_dir + "projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
input_dir = work_dir + "fgwas_input/"
in_file = input_dir + "ukbb_diamante-euro.fgwas.gz"

# functions

def write_file():
    fin = gzip.open(in_file,'rb')
    fout = open(work_dir+"chromosome_limits.txt",'w')
    header = fin.readline()
    chrom_list = []
    start_list = []
    end_list = []
    dic = {}
    prev = "NA"
    count = 0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        chrom, pos = l[0], l[2]
        if chrom not in chrom_list:
            print ("\n" + chrom)
            chrom_list.append(chrom)
            start_list.append(pos)
            end_list.append(prev)
        prev = pos
    fin.close()
    end_list.append(pos)
    print chrom_list
    print start_list
    print end_list
    end_list.pop(0)
    for i in range(0,len(chrom_list)):
        chrom = chrom_list[i]
        start = start_list[i]
        end = end_list[i]
        write_list = [chrom,start,str(int(end)+1)]
        fout.write("\t".join(write_list)+"\n")
    fout.close()


def main():
    write_file()

if (__name__=="__main__"): main()
