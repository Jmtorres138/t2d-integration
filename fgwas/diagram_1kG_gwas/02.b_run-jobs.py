#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run 01a jobs
Usage: python 01.d_run-jobs.py
'''
# libraries
import sys,os
import gzip
import subprocess as sp

# module load R/3.3.1

work_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/"

def run_job(chrom,part):
    command = ["Rscript", work_dir+"02.a_run-annotations.R", str(chrom),str(part)]
    command = "\t".join(command)
    jobfile = work_dir + "jobs/chr"+str(chrom)+".p"+str(part)+".job.sh"
    fout = open(jobfile,'w')
    script = '''
#$ -N %s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %slogs/%s.error
#$ -o %slogs/%s.out
#$ -V

%s

    ''' % ("chr"+str(chrom)+".p"+str(part),
            work_dir,"chr"+str(chrom)+".p"+str(part),
            work_dir,"chr"+str(chrom)+".p"+str(part),
            command)

    fout.write(script)
    fout.close()
    comm = ["qsub", jobfile]
    sp.check_call(comm)


def run_jobs():
    for c in range(21,23):
        for p in range(1,21):
            run_job(c,p)

def main():
    run_jobs()

if (__name__=="__main__"): main()
