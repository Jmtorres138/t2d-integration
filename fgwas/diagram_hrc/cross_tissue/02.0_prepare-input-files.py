#!/usr/bin/python -O
# Jason Matthew Torres
'''
Prepare annotation input files for second round
Usage: python 02.0_prepare-input-files.py
'''
# libraries
import sys,os,gzip
import subprocess as sp

#globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"
sortBed = "/apps/well/bedtools/2.24.0/sortBed"

cur_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/"

genom_bed = "/well/got2d/jason/reference/genomic/all_genomic.bed"

param_file_prefix = "best-joint-model"

path_dic = {'liver':["/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/bed_files/Liver.chromatinStates.bed",
"/well/got2d/jason/reference/encode/liver/liver.hg19.bed"],
'islet':["/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/bed_files/Islets.chromatinStates.bed",
"/well/got2d/jason/reference/islet/atac_peaks/all_successful_ISL_peaks.sorted.merged.bed"],
'muscle':["/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/bed_files/SkeletalMuscle.chromatinStates.bed",
"/well/got2d/jason/reference/encode/muscle/muscle.hg19.bed"],
'adipose':["/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/bed_files/Adipose.chromatinStates.bed",
"/well/got2d/jason/reference/encode/adipose/adipose.hg19.bed"]}

split_numbers = ["1","2","3","9","10","11"] # based on Varshney scheme

#param_file = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/fgwas_output_liver/best-joint-model.params"

#functions

def get_param_list(param_file):
    fin = open(param_file,'r')
    fin.readline()
    param_list = []
    for line in fin:
        l = line.strip().split()
        annot = l[0].split("_ln")[0]
        param_list.append(annot)
    fin.close()
    param_list.remove('pi_region')
    param_list.remove('distance_tss_0_5000')
    return param_list


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

def bed_intersect(chromseg_bed,atac_bed,annot,input_dir_2):
    sys.stdout.write("Intersecting chromatin states and atac peaks...\n")
    command_list1 = [bedtools, "intersect", "-a", chromseg_bed, "-b", atac_bed,
                    ">", input_dir_2+annot+"_open"]
    command1 = " ".join(command_list1)
    sp.check_call(command1,shell=True)
    append_to_column(input_dir_2+annot+"_open",3,"_Open")
    sys.stdout.write("Getting chromatin states that DO NOT overlap atac peaks...\n")
    command_list2 = [bedtools, "intersect", "-a", chromseg_bed, "-b", atac_bed,
                    "-v",">", input_dir_2+annot+"_closed"]
    command2 = " ".join(command_list2)
    sp.check_call(command2,shell=True)
    append_to_column(input_dir_2+annot+"_closed",3,"_Closed")
    command_list3 = ["cat", input_dir_2+annot+"_open", input_dir_2+annot+"_closed",
                     " | ", sortBed, ">", input_dir_2+annot+".bed"]
    command3 = " ".join(command_list3)
    sp.check_call(command3,shell=True)
    command4 = " ".join(["awk '!seen[$0]++'",input_dir_2+annot+".bed",">",
                         input_dir_2+"temp"])
    sp.check_call(command4,shell=True)
    os.rename(input_dir_2+"temp",input_dir_2+annot+".bed")

def build_anno_input(tissue):
    chromseg_bed = path_dic[tissue][0]
    atac_bed = path_dic[tissue][1]
    input_dir_1 = cur_dir + "fgwas_input_" + tissue + "/"
    output_dir_1 = cur_dir + "fgwas_output_" + tissue + "/"

    input_dir_2 = input_dir_1 + "round2/"
    output_dir_2 = output_dir_1 + "round2/"
    if os.path.isdir(input_dir_2)==False:
        os.mkdir(input_dir_2)
    if os.path.isdir(output_dir_2)==False:
        os.mkdir(output_dir_2)
    param_file = output_dir_1 + param_file_prefix + ".params"
    param_list = get_param_list(param_file)

    for annot in param_list:
        print("\t"+annot)
        fout = open(input_dir_2 + annot+".bed",'w')
        for f in [genom_bed,chromseg_bed]:
            fin = open(f,'r')
            for line in fin:
                l = line.strip().split()
                a = l[3]
                a = a.replace("/","_")

                if a==annot:
                    l[3] = a
                    fout.write("\t".join(l)+"\n")
            fin.close()
        fout.close()

    for annot in param_list:
        if annot.split("_")[0] in split_numbers:
            print annot
            annot_bed = input_dir_2 + annot+".bed"
            bed_intersect(annot_bed,atac_bed,annot,input_dir_2)

    command_list1 = ["cat", input_dir_2+"*.bed", ">",input_dir_2+"anno_input"]
    command1 = " ".join(command_list1)
    sp.check_call(command1,shell=True)
    command_list2 = ["rm", input_dir_2+"*.bed", input_dir_2+"*open",input_dir_2+"*closed"]
    command2 = " ".join(command_list2)
    sp.check_call(command2,shell=True)
    command_list3 = ["mv", input_dir_2+"anno_input",input_dir_2+"anno_input.bed"]
    command3 = " ".join(command_list3)
    sp.check_call(command3,shell=True)

def main():
    tiss_list = ["islet"] # ["muscle"] #"liver","adipose"]
    for tiss in tiss_list:
        print ("\n"+tiss)
        build_anno_input(tiss)

if (__name__=="__main__"): main()
