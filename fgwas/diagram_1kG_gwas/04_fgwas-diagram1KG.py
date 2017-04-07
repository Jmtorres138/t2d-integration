#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run fGWAS across annotations
Usage: python 04_fgwas-diagram1KG.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import operator
import time
from select import select
from math import ceil
from math import floor


# globals
fgwas = "/users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/"
in_dir=home_dir+"fgwas_input/"
out_dir = home_dir + "fgwas_output/"
input_file=in_dir+"diagram_1kG_fgwas_24annot.txt.gz"
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)

def step1(start_index=10):
    '''
    Start index is the first index for an annotation in the file
    Here, the start index is 10 (column 11)
    '''
    print("Running each annotation separately....")
    fin = gzip.open(input_file,'rb')
    annot_list = fin.readline().strip().split()[start_index:]
    fin.close()
    for annot in annot_list:
        print annot
        job_file = job_dir+"job_"+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",  "-dists",
                            annot+":"+home_dir+"dist_model", "-o", out_dir+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            annot, "-o", out_dir+annot]
        command = " ".join(command_list)
        script='''
#$ -N job_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (annot, log_dir,"job_"+annot,log_dir,"job_"+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        sp.check_call(call)
    raw_input('Press enter to continue when jobs are complete: ')

def step2(start_index=10,first_go=True):
    print("Finding annotation with highest model likelihood..")
    fin1 = gzip.open(input_file,'rb')
    annot_list = fin1.readline().strip().split()[start_index:]
    fin1.close()
    track_dic = {}
    for annot in annot_list:
        f = out_dir+annot+".llk"
        fin=open(f,'r')
        ## Use next line for ln(lk) directly
        l = fin.readline().strip().split()
        ## Use next 3 line for AIC
        ##fin.readline()
        ##fin.readline()
        ##l = fin.readline().strip().split()
        fin.close()
        if l[0]=="ln(lk):":
        #if l[0]=="AIC:":
            track_dic[annot] = float(l[1])
    sorted_annot = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_annot.reverse() # use for ln(lk)
    top_annot = sorted_annot[0][0]
    top_val = sorted_annot[0][1]
    print "Top annotation: %s ; ln(lk): %f" % (top_annot, top_val)
    #print "Top annotation: %s ; AIC: %f" % (top_annot, top_val)
    for annot in sorted_annot[1:]:
        annot = annot[0]
        #print annot
        job_file = job_dir+"job_"+top_annot+"-"+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", annot+":"+home_dir+"dist_model",
                            "-w", annot, "-o", out_dir+top_annot+"+"+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            top_annot+"+"+annot, "-o", out_dir+top_annot+"+"+annot]
        command = " ".join(command_list)
        script='''
#$ -N job_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (top_annot+"-"+annot, log_dir,"job_"+top_annot+"-"+annot,
        log_dir,"job_"+top_annot+"-"+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        if first_go==True:
            sp.check_call(call) # Run this only on the first attempt
    if first_go==True:
        raw_input('Press enter to continue when jobs are complete: ')
    top = [top_annot,top_val]
    return(top)

def run_models(fixed_list,eval_list):
    print "Fixed annotations: " + ", ".join(fixed_list)
    fixed_name = "-".join(fixed_list)
    fixed = "+".join(fixed_list)
    for annot in eval_list:
        job_file = job_dir+"job_"+fixed_name+"-"+annot+".sh"
        fout=open(job_file,'w')
        if annot == "distance_tss":
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", annot+":"+home_dir+"dist_model",
                            "-w", fixed, "-o", out_dir+fixed+"+"+annot]
        elif "distance_tss" in fixed_list:
            #print True
            temp_list = list(fixed_list)
            temp_list.remove("distance_tss")
            fixed_sub = "+".join(temp_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", fixed_sub+"+"+annot, "-o", out_dir+fixed+"+"+annot]
        else:
            command_list = [fgwas, "-i", input_file, "-cc", "-w",
                            fixed+"+"+annot, "-o", out_dir+fixed+"+"+annot]
        #print command_list
        command = " ".join(command_list)
        script='''
#$ -N job_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (fixed+"-"+annot, log_dir,"job_"+fixed+"-"+annot,
        log_dir,"job_"+fixed+"-"+annot, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        sp.check_call(call)


def step3(top_annot, top_val,start_index=10,first_pass=True):
    fin1 = gzip.open(input_file,'rb')
    annot_list = fin1.readline().strip().split()[start_index:]
    fin1.close()

    print top_annot,top_val
    top_annot_list = top_annot.split("+")
    annot_list = [x for x in annot_list if x not in top_annot_list]
    out_list = [top_annot+"+"+ x for x in annot_list]
    if first_pass == True:
        run_models(top_annot_list,annot_list)

    raw_input('Press enter to continue when jobs are complete: ')

    track_dic = {}
    for name in out_list:
        f = out_dir+name+".llk"
        fin=open(f,'r')
        #fin.readline() #use for AIC
        #fin.readline() #use for AIC
        l = fin.readline().strip().split() # get 3rd line
        fin.close()
        if l[0]=="ln(lk):":
        #if l[0]=="AIC:":
            track_dic[name] = float(l[1])
    sorted_annot = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_annot.reverse() #used when comparing ln(lk)
    print str(top_val) #+ " : " + str(floor(top_val))
    ## Next line(s) used for comparing ln(lk)
    sig_list = [x for x in sorted_annot if float(x[1]) > float(top_val)]
    #sig_list = [x for x in sorted_annot if float(x[1]) > ceil(float(top_val))]
    ## Next line used for comparing AIC, comment out if unecessary
    #sig_list = [x for x in sorted_annot if float(x[1]) < floor(float(top_val))]
    #sig_list = [x for x in sorted_annot if float(x[1]) < float(top_val)]

    print sig_list
    print len(sig_list)
    print sig_list[0]
    return sig_list[0]


def step4(model_list,first_go=True):
    print "Finding the penalty with the best cross-validation likelihood..."
    p_list = ["0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50"]
    model_name = "-".join(model_list)
    model = "+".join(model_list)
    for p in p_list:
        job_file = job_dir+"job_"+model_name+"-"+p+".sh"
        fout=open(job_file,'w')
        if "distance_tss" in model_list:
            #print True
            temp_list = list(model_list)
            temp_list.remove("distance_tss")
            model_sub = "+".join(temp_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", model_sub, "-p", p, "-xv", "-print",
                            "-o", out_dir+model+"-p"+p]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",  "-w",
                            model, "-p", p, "-xv", "-print",
                            "-o", out_dir+model+"-p"+p]
        #print command_list
        command = " ".join(command_list)
        script='''
#$ -N job_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q long.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (model_name+"-p"+p, log_dir,"job_"+model_name+"-p"+p,
        log_dir,"job_"+model_name+"-p"+p, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        if first_go==True:
            sp.check_call(call)
    if first_go==True:
        raw_input('Press enter to continue when jobs are complete: ')
    print "Finding best parameter value..."
    track_dic = {}
    for p in p_list:
        fin = open(out_dir+model+"-p"+p+".ridgeparams",'r')
        line_list = fin.readlines()
        fin.close()
        line = line_list[-1]
        llk = line.strip().split()[-1]
        track_dic[p]=llk
        print (p + ": " + str(llk))
    sorted_p = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_p.reverse() #used when comparing ln(lk)
    best = sorted_p[0]
    print "Optimal parameter value evaluated: %s"   % best[0]
    return best


def step5_6(model_list,best_p,best_llk,first_go=True):
    print "Test dropping each annotation from the model, using cross-validation likelihood"
    print "Keep dropping annotations as long as the cross-validation likelihood keeps increasing"
    for mod in model_list:
        job_file = job_dir+"job_drop-"+mod+".sh"
        fout=open(job_file,'w')
        if "distance_tss" in model_list:
            #print True
            temp_list = list(model_list)
            temp_list.remove("distance_tss")
            model_sub = "+".join(temp_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+mod]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-w", mod, "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+mod]
        #print command_list
        command = " ".join(command_list)
        script='''
#$ -N job_drop-%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (mod, log_dir,"job_drop-"+mod,
        log_dir,"job_drop-"+mod, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        if first_go==True:
            sp.check_call(call)
    if first_go==True:
        raw_input('Press enter to continue when jobs are complete: ')
    print "The best likelihood in full model: %s" % str(best_llk)
    track_dic = {}
    for mod in model_list:
        fin = open(out_dir+"drop-"+mod+".ridgeparams",'r')
        line_list = fin.readlines()
        fin.close()
        line = line_list[-1]
        llk = line.strip().split()[-1]
        track_dic[mod]=llk
        print (mod + ": " + str(llk))
    sorted_mods = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_mods.reverse() #used when comparing ln(lk)
    check_list = [x for x in sorted_mods if float(x[1]) > float(best_llk)]
    print check_list

def step5_6_sequential(model_list,best_p,best_llk,first_go=True):
    print "Test dropping each annotation from the model, using cross-validation likelihood"
    print "Keep dropping annotations as long as the cross-validation likelihood keeps increasing"
    for mod in model_list:
        keep_list = list(model_list)
        dropped_mod = mod
        keep_list.remove(mod)
        keep_mods = "+".join(keep_list)
        job_file = job_dir+"job_drop-"+mod+".sh"
        fout=open(job_file,'w')
        if "distance_tss" in keep_list:
            #print True
            #temp_list = list(model_list)
            keep_list.remove("distance_tss")
            model_sub = "+".join(keep_list)
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-dists", "distance_tss"+":"+home_dir+"dist_model",
                            "-w", keep_mods,
                            "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+mod]
        else:
            command_list = [fgwas, "-i", input_file, "-cc",
                            "-w", keep_mods, "-p", best_p, "-xv", "-print",
                            "-o", out_dir+"drop-"+mod]
        #print command_list
        command = " ".join(command_list)
        script='''
#$ -N job_drop-%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
        ''' % (mod, log_dir,"job_drop-"+mod,
        log_dir,"job_drop-"+mod, command)
        fout.write(script)
        fout.close()
        call = ["qsub", job_file]
        if first_go==True:
            sp.check_call(call)
    if first_go==True:
        raw_input('Press enter to continue when jobs are complete: ')
    print "The best likelihood in full model: %s" % str(best_llk)
    track_dic = {}
    for mod in model_list:
        fin = open(out_dir+"drop-"+mod+".ridgeparams",'r')
        line_list = fin.readlines()
        fin.close()
        line = line_list[-1]
        llk = line.strip().split()[-1]
        track_dic[mod]=llk
        print ("dropped " + mod + ": " + str(llk))
    sorted_mods = sorted(track_dic.items(),key=operator.itemgetter(1))
    sorted_mods.reverse() #used when comparing ln(lk)
    check_list = [x for x in sorted_mods if float(x[1]) > float(best_llk)]
    try:
        best = check_list[0]
        best_dropped_mod = best[0]
        best_dropped_llk = best[1]
        report_list = list(model_list)
        report_list.remove(best_dropped_mod)
        print ("Best dropped model: %s" % best_dropped_mod)
        print ("Best dropped llk: %s" % best_dropped_llk)
        print ("Annotations to keep: %s" % ",".join(report_list))
        return best_dropped_mod,best_dropped_llk, report_list
    except:
        print ("Dropping models didn't improve cross-validated likelihood")
        print ("Keep the current model!")
        return False,False,model_list




def main():
    ### Steps 1-2
    #step1()
    #top = step2()
    ##top = step2(first_go=False)
    #top_annot, top_val = top[0], top[1]
    #print "Top annotation: " + str(top_annot) + " Value: " + str(top_val)
    #Top annotation: islet_state8 Value: 781.647


    ### Step 3
    #top_annot,top_val = "islet_state8", 781.647 # determined from steps 1,2
    #('islet_state8+islet_state9+cds+islet_state6+distance_tss+intron', 800.90999999999997)
    #top_annot,top_val = "islet_state8+islet_state9+cds+islet_state6+distance_tss+intron", 800.91
    #top_annot,top_val = "islet_state8+islet_state9+cds+islet_state6+distance_tss+intron+transcript+islet_state7", 802.614
    top_annot,top_val = "islet_state8+islet_state9+cds+islet_state6+distance_tss+intron+transcript+islet_state7+utr_5+islet_state13+islet_state15+utr_3+promoter+islet_atac+exon+islet_state11+islet_state3+islet_state10+islet_state4+islet_state14", 806.045

    #iter1 = step3(top_annot, top_val,first_pass=False)
    #iter2 = step3(iter1[0], iter1[1])
    #iter3 = step3(iter2[0], iter2[1])

    ## Note: keep running iterations until additional annotations no longer
    ## improve the model (using ln(lk) or AIC criterion above)

    ### Step 4
    ## There are 20 annotations in the resulting model (shown below)
    model_list = ["islet_state8","islet_state9","cds","islet_state6","distance_tss",
                  "intron","transcript","islet_state7","utr_5","islet_state13","islet_state15",
                  "utr_3","promoter","islet_atac","exon","islet_state11","islet_state3","islet_state10",
                  "islet_state4","islet_state14"]

    best_p, best_llk = step4(model_list)
    #best_p, best_llk = step4(model_list,first_go=False)

    # Step 5
    #mod,llk,keep = step5_6_sequential(model_list,best_p,best_llk,first_go=True)
    #
    # Need to keep iterating
    #print keep
    #model_list = ["intron","islet_state9","utr_3"]
    #best_p, best_llk = step4(model_list)
    #mod,llk,keep = step5_6_sequential(model_list,best_p,best_llk,first_go=True)

    #model_list = ["intron","utr_3"]
    #best_p, best_llk = step4(model_list)
    #mod,llk,keep = step5_6_sequential(model_list,best_p,best_llk,first_go=True)

    # The best model is ["intron","utr_3"]
if (__name__=="__main__"): main()
