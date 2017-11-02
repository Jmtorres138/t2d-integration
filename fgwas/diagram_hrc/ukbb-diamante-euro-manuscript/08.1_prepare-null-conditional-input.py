#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run fGWAS across annotations
Usage: python 08.1_null-fgwas-wrapper.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import operator
import time
from select import select
from math import ceil
from math import floor
import moniter_rescomp_jobs

# globals
fgwas = "LD_LIBRARY_PATH=/apps/well/gsl/2.2.1-gcc4.9.3/lib /users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
ref_dir = home_dir + "conditional/fgwas_input_files/"
out_dir=home_dir+"null/fgwas_input_files/"
uncond_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/fgwas_input/"
uncond_name = "ukbb_diamante-euro.fgwas.gz"

def prepare_inputs():
    print uncond_name
    infile = uncond_dir + uncond_name
    outfile = out_dir + "null_" + uncond_name
    outfile = outfile.split(".gz")[0]
    command_list = ["zcat",infile," | cut -f 1,2,3,4,5,6,7,8,9",">", outfile]
    command = " ".join(command_list)
    sp.check_call(command,shell=True)
    command_list = ["gzip",outfile]
    command = " ".join(command_list)
    sp.check_call(command,shell=True)
    file_list = os.listdir(ref_dir)
    file_list = [x for x in file_list if ".fgwas.gz" in x]
    for f in file_list:
        print f
        infile = ref_dir + f
        outfile = out_dir + "null_" + f
        outfile = outfile.split(".gz")[0]
        command_list = ["zcat",infile," | cut -f 1,2,3,4,5,6,7,8,9",">", outfile]
        command = " ".join(command_list)
        sp.check_call(command,shell=True)
        command_list = ["gzip",outfile]
        command = " ".join(command_list)
        sp.check_call(command,shell=True)

def main():
    prepare_inputs()

if (__name__=="__main__"): main()
