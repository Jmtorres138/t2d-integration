#!/usr/bin/python -O
# Jason Matthew Torres
'''
Profile annotations in the fgwas input file
Usage: python JTprofile-annot.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import time

bedtools = "/apps/well/bedtools/2.24.0/bedtools"
fgwas_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/diamante-euro-manuscript/"
fgwas_input_dir = fgwas_dir + "fgwas_input/"
bed_dir = fgwas_input_dir + "bed_files/"
if os.path.isdir(bed_dir)==False:
    os.makedirs(bed_dir)
annot_bed = fgwas_input_dir + "anno_input.bed"
gwas_bed = "/well/got2d/jason/reference/gwas/diagram_hrc/diagram_hrc.bed"
fgwas_file = fgwas_input_dir  + "diagram_hrc.fgwas.gz"
start_index=9



def write_beds(annot_file):
    fin = open(annot_bed,'r')
    count = 0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        annot = l[3]
        fname = bed_dir + annot +".bed"
        if os.path.exists(fname)==False:
            fout = open(fname,'w')
        else:
            fout = open(fname,'a')
        fout.write("\t".join(l)+"\n")
        fout.close()
    fin.close()
    print "\nJob Complete"

def three_intersect_slow(gwas_bed,a1,a2): # faster ~ 57 seconds
    #sys.stdout.write("Intersecting...\n")
    command = [bedtools,"intersect","-u","-a",gwas_bed,"-b", bed_dir+a1+".bed",">",bed_dir+"temp.bed"]
    sp.check_call(" ".join(command),shell=True)
    command = ["cut -f 1,2,3", bed_dir+"temp.bed", ">", bed_dir+"check.bed"]
    sp.check_call(" ".join(command),shell=True)
    command = ["uniq", bed_dir+"check.bed", "|", "wc" ,"-l"]
    p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    a1_count = output.strip()
    os.remove(bed_dir+"temp.bed")
    os.remove(bed_dir+"check.bed")

    if a1 != a2:
        command = [bedtools,"intersect","-u","-a",gwas_bed,"-b",bed_dir+a2+".bed",">",bed_dir+"temp.bed"]
        sp.check_call(" ".join(command),shell=True)
        command = ["cut -f 1,2,3", bed_dir+"temp.bed", ">", bed_dir+"check.bed"]
        sp.check_call(" ".join(command),shell=True)
        command = ["uniq", bed_dir+"check.bed", "|", "wc" ,"-l"]
        p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        a2_count = output.strip()
        os.remove(bed_dir+"temp.bed")
        os.remove(bed_dir+"check.bed")

        command = [bedtools,"intersect","-c","-a",gwas_bed,"-b", bed_dir+a1+".bed",
                   bed_dir+a2+".bed",">",bed_dir+"temp.bed"]
        sp.check_call(" ".join(command),shell=True)
        command = ["cut -f 10",bed_dir+"temp.bed","|","grep 2","|","wc -l"]
        p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        overlap = output.strip()
        os.remove(bed_dir+"temp.bed")
    else:
        a2_count, overlap = a1_count, a1_count

    return [a1_count, a2_count, overlap]


def three_intersect(gwas_bed,a1,a2): # slower ~62 seconds
    fin = gzip.open(fgwas_file,'rb')
    head_list = fin.readline().strip().split()
    fin.close()
    #print head_list
    a1_index = head_list.index(a1) + 1 # 1-based reference
    a2_index = head_list.index(a2) + 1 # 1-based reference

    command = ["zcat", fgwas_file, "|", "cut -f",str(a1_index),"|", "grep 1","|","wc -l"]
    p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    a1_count = output.strip()

    if a1 != a2:

        command = ["zcat", fgwas_file, "|", "cut -f",str(a2_index),"|", "grep 1","|","wc -l"]
        p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        a2_count = output.strip()

        command = [bedtools,"intersect","-c","-a",gwas_bed,"-b", bed_dir+a1+".bed",
                   bed_dir+a2+".bed",">",bed_dir+"temp.bed"]
        sp.check_call(" ".join(command),shell=True)
        command = ["cut -f 10",bed_dir+"temp.bed","|","grep 2","|","wc -l"]
        p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        overlap = output.strip()
        os.remove(bed_dir+"temp.bed")
    else:
        a2_count, overlap = a1_count, a1_count

    return [a1_count, a2_count, overlap]


def write_profile_file():
    sys.stdout.write("Writing output file...\n")
    annot_list = [x.split(".bed")[0] for x in os.listdir(bed_dir) if ".bed" in x and "temp" not in x]
    print annot_list 
    fout = open(fgwas_input_dir+"gwas-annotation-profile.txt",'w')
    fout.write("\t".join(["Name1","Count1","Name2","Count2","Overlap"])+"\n")
    track_dic = {}
    for annot in annot_list:
        iter_list = list(annot_list)
        for a in iter_list:
            sys.stdout.write("%s : %s\n" % (annot,a))
            annot_key = annot+":"+a
            try:
                count_list = track_dic[annot_key]
            except:
                count_list = three_intersect(gwas_bed,annot,a)
                track_dic[annot_key] = count_list
                alt_key = a+":"+annot
                track_dic[alt_key] = [count_list[1],count_list[0],count_list[2]]
            write_list = [annot,str(count_list[0]),a,str(count_list[1]),str(count_list[2])]
            fout.write("\t".join(write_list)+"\n")
    fout.close()

def main():
    #write_beds(annot_bed)
    #start = time.time()
    #print three_intersect(gwas_bed,"exon","coding")
    #end = time.time()
    #print(end - start)
    #start = time.time()
    #print three_intersect_v2(gwas_bed,"exon","coding")
    #end = time.time()
    #print(end - start)
    write_profile_file()

if (__name__=="__main__"): main()
