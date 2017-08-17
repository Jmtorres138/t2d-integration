#!/usr/bin/python -O
# Jason Matthew Torres
'''
Generate input file for varshney_2016 islet fgwas run
Usage: python JTbuild_annot.py
'''
# libraries
import sys,os,gzip
import subprocess as sp


#globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"
sortBed = "/apps/well/bedtools/2.24.0/sortBed"

chromseg_dir = "/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/bed_files/"#
genom_bed = "/well/got2d/jason/reference/genomic/all_genomic.bed"
cur_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/atac/"
input_dir = cur_dir + "fgwas_input/"
annot_dir = input_dir + "annot/"
if os.path.exists(annot_dir)==False:
    os.makedirs(annot_dir)
if os.path.exists(input_dir)==False:
    os.makedirs(input_dir)

isl_atac_bed = "/well/got2d/jason/reference/islet/atac_peaks/all_successful_ISL_peaks.sorted.merged.bed" # Matthias published
adi_atac_bed = "/well/got2d/jason/reference/encode/adipose/adipose.hg19.bed"
liv_atac_bed = "/well/got2d/jason/reference/encode/liver/liver.hg19.bed"
mus_atac_bed = "/well/got2d/jason/reference/encode/muscle/muscle.hg19.bed"

# functions

def create_beds():
    sys.stdout.write("\nAnnotation: coding\n")
    fin = open(genom_bed,'r')
    fout = open(annot_dir+"coding.bed",'w')
    count = 0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d"%count)
        sys.stdout.flush()
        l = line.strip().split()
        if l[3] == "coding":
            fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()

    annot_list = ["Islets","SkeletalMuscle","Adipose","Liver"]
    for annot in annot_list:
        sys.stdout.write("\nAnnotation: %s\n" % annot)
        fin = open(chromseg_dir+annot+".chromatinStates.bed",'r')
        fout = open(annot_dir+annot+"_enhancer.bed",'w')
        count = 0
        for line in fin:
            count+=1
            sys.stdout.write("\r%d"%count)
            sys.stdout.flush()
            l = line.strip().split()
            if "Active_enhancer" in l[3] or "Weak_enhancer" in l[3]:
                l[3] = annot+"_enhancer"
                fout.write("\t".join(l)+"\n")
        fin.close()
        fout.close()
        # Disregard SNPs that are in the coding set
        sp.check_call(" ".join(["mv",annot_dir+annot+"_enhancer.bed",annot_dir+"temp.bed"]),shell=True)
        command = [bedtools, "intersect","-v", "-a", annot_dir+"temp.bed",
                   "-b", annot_dir+"coding.bed", ">", annot_dir+"temp2.bed"]
        sp.check_call(" ".join(command),shell=True)

        if annot == "Islets":
            atac_bed = isl_atac_bed
        elif annot == "SkeletalMuscle":
            atac_bed = mus_atac_bed
        elif annot == "Adipose":
            atac_bed = adi_atac_bed
        elif annot == "Liver":
            atac_bed = liv_atac_bed
        else:
            raise NameError("Not a valid tissue")
        sys.stdout.write("\nIntersecting with ATAC...\n")
        command = [bedtools, "intersect","-wa", "-a", annot_dir+"temp2.bed",
                   "-b", atac_bed, ">", annot_dir+"temp3.bed"]
        sp.check_call(" ".join(command),shell=True)
        command = ["uniq", annot_dir+"temp3.bed",">",annot_dir+annot+"_enhancer.bed"]#
        sp.check_call(" ".join(command),shell=True)

    print("\nCreated bed files")

def create_anno_input():
    sys.stdout.write("Creating annotation input bed file for fgwas...\n")
    command = ["cat",annot_dir+"Liver"+"_enhancer.bed",
               annot_dir+"SkeletalMuscle"+"_enhancer.bed",
               annot_dir+"Adipose"+"_enhancer.bed","|","cut -f 1,2,3",">",annot_dir+"temp.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "sort","-i", annot_dir+"temp.bed",">",annot_dir+"temp2.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "merge", "-i",  annot_dir+"temp2.bed",">",annot_dir+"temp3.bed"]
    sp.check_call(" ".join(command),shell=True)

    fin = open(annot_dir+"temp3.bed",'r')
    fout = open(annot_dir+"peripheral.bed",'w')
    for line in fin:
        l = line.strip().split()
        l.append("peripheral")
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    command = ["rm", annot_dir+"temp.bed",annot_dir+"temp2.bed",annot_dir+"temp3.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "intersect","-v", "-a", annot_dir+"Islets_enhancer.bed",
               "-b", annot_dir+"peripheral.bed", ">", annot_dir+"islet-specific_enhancer.bed"]
    sp.check_call(" ".join(command),shell=True)
    command = [bedtools, "intersect","-v", "-a", annot_dir+"peripheral.bed",
               "-b", annot_dir+"Islets_enhancer.bed", ">", annot_dir+"peripheral-specific_enhancer.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "intersect","-u", "-a", annot_dir+"Islets_enhancer.bed",
               "-b", annot_dir+"peripheral.bed", "|", "cut -f 1,2,3",">", annot_dir+"temp1.bed"]
    sp.check_call(" ".join(command),shell=True)
    command = [bedtools, "intersect","-u", "-a", annot_dir+"peripheral.bed",
               "-b", annot_dir+"Islets_enhancer.bed", "|", "cut -f 1,2,3",">", annot_dir+"temp2.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = ["cat",annot_dir+"temp1.bed",annot_dir+"temp2.bed","|","uniq",">",annot_dir+"temp3.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "sort","-i", annot_dir+"temp3.bed",">",annot_dir+"temp4.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = [bedtools, "merge", "-i",  annot_dir+"temp4.bed",">",annot_dir+"temp5.bed"]
    sp.check_call(" ".join(command),shell=True)

    fin = open(annot_dir+"temp5.bed",'r')
    fout = open(annot_dir+"shared_enhancer.bed",'w')
    for line in fin:
        l = line.strip().split()
        l.append("shared")
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()

    command = ["rm", annot_dir+"temp*"]
    sp.check_call(" ".join(command),shell=True)

    command = ["cat", annot_dir+"coding.bed",annot_dir+"islet-specific_enhancer.bed",
               annot_dir+"peripheral-specific_enhancer.bed",annot_dir+"shared_enhancer.bed",
               ">" ,input_dir+"anno_input.bed"]
    sp.check_call(" ".join(command),shell=True)

    command = ["rm -r", annot_dir]
    sp.check_call(" ".join(command),shell=True)


def main():
    create_beds()
    create_anno_input()

if (__name__=="__main__"): main()
