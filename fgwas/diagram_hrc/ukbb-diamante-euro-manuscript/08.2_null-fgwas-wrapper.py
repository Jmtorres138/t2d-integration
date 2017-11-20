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
python = "/apps/well/python/2.7.11/bin/python"
rscript = "/apps/well/R/3.3.1/bin/Rscript" #"/apps/well/R/3.3.1/bin/Rscript"
fgwas = "LD_LIBRARY_PATH=/apps/well/gsl/2.2.1-gcc4.9.3/lib /users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
in_dir=home_dir+"null/fgwas_input_files/"
out_dir = home_dir + "null/fgwas_output_files/"
input_file=in_dir+"null_ukbb_diamante-euro.fgwas.gz" # Optional: Run 01.1 script and use this for abbreviated annotations: in_dir+"diagram_hrc.renamed.fgwas.gz"
part_dir = home_dir+"partition_bedfiles/"

job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)
if os.path.isdir(job_dir)==False:
    os.mkdir(job_dir)
if os.path.isdir(log_dir)==False:
    os.mkdir(log_dir)
#start_index = 9 #0-based index of column in fgwas input file where annotations start
job_prefix = "null_"



manual_list = ["210_1","210_2","86_1","87_1", "132_1", "132_2", "132_3", "132_4", "132_5",
                "133_1","133_2","133_3","133_4","133_5","133_6","133_7","133_8","133_9","133_10"]
missing_list = ["20_1","86_2","87_2","163_2","164_2"]

gwas_bed_dir = "/well/got2d/jason/reference/gwas/diamante-ukbb_hrc/conditioned/gwas_bedfiles/"
loc_ref_file = "/well/got2d/jason/reference/gwas/diamante-ukbb_hrc/conditioned/list.for.credible.sets.ALL.txt"

def fgwas_run():
    bed_file = home_dir + "loci-partition.bed"
    job_file = job_dir+"null_fgwas_run.sh"
    fout=open(job_file,'w')
    command_list = [fgwas, "-i", input_file, "-cc", "-bed", bed_file,
                    #"-w", "+".join(annot_list), "-print","-o", out_dir+"null_fgwas_run_loci-partition"]
                    "-print","-o", out_dir+"null_fgwas_run_loci-partition"]
    command = " ".join(command_list)
    script='''
    #$ -N %s%s
    #$ -pe shmem 1
    #$ -P mccarthy.prjc
    #$ -q short.qc
    #$ -e %s%s.error
    #$ -o %s%s.out
    echo "start time" `date`
    %s
    echo "end time" `date`
    ''' % (job_prefix,"loci-partition", log_dir,"loci-partition",
    log_dir,"loci-partition", command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)


def run_loc_job(loc_id):
    out_dir=home_dir+"null/conditional/fgwas_output_files/"
    if os.path.isdir(out_dir)==False:
        os.mkdir(out_dir)

    # Note: loc_id can be locus id (e.g. "188_2") or condition ref name (e.g. "cond2")
    if loc_id in manual_list:
        fgwas_input_file = in_dir + "null_ukbb_diamante-euro." + loc_id + ".fgwas.gz"
    elif loc_id in missing_list:
        cond = loc_id.split("_")[1]
        fgwas_input_file = in_dir + "null_ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    elif "_" in loc_id:
        cond = loc_id.split("_")[1]
        fgwas_input_file = in_dir + "null_ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    else:
        cond = loc_id.split("cond")[1]
        fgwas_input_file = in_dir + "null_ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    bed_part_file = part_dir + "loci-partition-"+loc_id + ".bed"

    loc_out_dir = out_dir + loc_id + "/"
    if os.path.isdir(loc_out_dir)==False:
        os.mkdir(loc_out_dir)

    job_file = job_dir+"job."+loc_id+".sh"
    fout=open(job_file,'w')

    #command_list1 = [rscript,"--vanilla",home_dir+"06.0_conditional-create-bed.R",loc_id]
    #command1 = " ".join(command_list1)


    command_list2 = [fgwas, "-i", fgwas_input_file, "-cc", "-bed", bed_part_file,
                    "-print","-o", loc_out_dir+"null_fgwas_run_loci-partition"]
    command2 = " ".join(command_list2)


    command_list3 = [rscript, "--vanilla",home_dir+"08.0_get-cond-block.R",loc_id]
    command3 = " ".join(command_list3)

    command_list4 = [python,home_dir+"08.0_subset_to_seg.py",loc_id]
    command4 = " ".join(command_list4)

    #command_list5 = [rscript, "--vanilla",home_dir+"06.0_functional_credible_sets.R",loc_id]
    #command5 = " ".join(command_list5)

    #command_list6 = ["rm", loc_out_dir+"fgwas_run_loci-partition.bfs.gz"]
    #command6 = " ".join(command_list6)

    script='''
#$ -N null_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V
echo "start time" `date`
%s
%s
%s
echo "end time" `date`
        ''' % (loc_id, log_dir,loc_id,log_dir,loc_id,
        command2, command3, command4)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)

def run_all_loci():
    loc_list = []
    for i in [1,2,3,4,5,6,7,8]:
        loc_list.append("cond"+str(i))
    loc_list = loc_list + manual_list
    loc_list = loc_list + missing_list
    for loc in loc_list:
        run_loc_job(loc)
    print loc_list

def run_cred_sets():
    loc_list = []
    #for i in [1,2,3,4,5,6,7,8]: # TEMPORARY comment out
    #    loc_list.append("cond"+str(i)) # TEMPORARY comment out
    loc_list = loc_list + manual_list
    loc_list = loc_list + missing_list
    for loc in loc_list:
        command_list5 = [rscript, "--vanilla",home_dir+"08.0_functional_credible_sets.R",loc]
        sp.check_call(command_list5)
    print loc_list


def main():
    #fgwas_run()
    #run_all_loci()
    run_cred_sets()

if (__name__=="__main__"): main()
