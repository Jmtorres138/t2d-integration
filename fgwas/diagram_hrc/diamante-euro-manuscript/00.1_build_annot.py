#!/usr/bin/python -O
# Jason Matthew Torres
'''
Generate input file for alternative fgwas run
Usage: python JTbuild_annot.py
'''
# libraries
import sys,os,gzip
import subprocess as sp


#globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"
sortBed = "/apps/well/bedtools/2.24.0/sortBed"
atac_bed = "/well/got2d/jason/reference/islet/atac_peaks/all_successful_ISL_peaks.sorted.merged.bed"
chromseg_bed = "/well/got2d/jason/reference/islet/chromatin_states/Pancreatic_islet_11_segments.bed"
genom_bed = "/well/got2d/jason/reference/genomic/all_genomic.bed"

cur_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/diamante-euro-manuscript/"
input_dir = cur_dir + "fgwas_input/"

# functions

def append_to_column(fname, pycolnum, text):
    fin = open(fname,'r')
    fout = open(cur_dir+"temp.bed",'w')
    for line in fin:
        l = line.strip().split()
        l[pycolnum] = l[pycolnum]+text
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    os.rename(cur_dir+"temp.bed",fname)

def bed_intersect(chromseg_bed,atac_bed):
    sys.stdout.write("Intersecting chromatin states and atac peaks...\n")
    command_list1 = [bedtools, "intersect", "-a", chromseg_bed, "-b", atac_bed,
                    ">", input_dir+"chrom_atac_intersect"]
    command1 = " ".join(command_list1)
    sp.check_call(command1,shell=True)
    append_to_column(input_dir+"chrom_atac_intersect",3,"_ATAC")
    sys.stdout.write("Getting chromatin states that DO NOT overlap atac peaks...\n")
    command_list2 = [bedtools, "intersect", "-a", chromseg_bed, "-b", atac_bed,
                    "-v",">", input_dir+"chrom_atac_no-overlap"]
    command2 = " ".join(command_list2)
    sp.check_call(command2,shell=True)
    append_to_column(input_dir+"chrom_atac_no-overlap",3,"_No_ATAC")

    sys.stdout.write("Write input file...\n")
    command_list3 = ["cat", input_dir+"chrom_atac_intersect", input_dir+"chrom_atac_no-overlap",
                     #genom_bed, atac_bed, " | ", sortBed, ">", input_dir+"anno_input.bed"]
                     genom_bed,  " | ", sortBed, ">", input_dir+"anno_input.bed"]

    command3 = " ".join(command_list3)
    sp.check_call(command3,shell=True)
    command4 = " ".join(["awk '!seen[$0]++'",input_dir+"anno_input.bed",">",
                         input_dir+"temp"])
    sp.check_call(command4,shell=True)
    os.rename(input_dir+"temp",input_dir+"anno_input.bed")



def main():
    bed_intersect(chromseg_bed,atac_bed)

if (__name__=="__main__"): main()
