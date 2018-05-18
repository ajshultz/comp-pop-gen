#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg') #Added to get working on the cluster
from matplotlib.backends.backend_pdf import PdfPages

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename):
    print("Opening %s"%config_filename)
    config_file = open(config_filename,"r")
    config_info = {}
    sample_ncbi_dict = {}
    sample_ena_dict = {}
    sample_dict = {}
    
    for line in config_file:
        if line[0] == "#":
            pass
        elif line == "\n":
            pass
        else:
            line=line.strip()
            line = line.split(" ")
            if line[0] == "--ABBV":
                config_info["abbv"] = line[1]
            elif line[0] == "--SAMPLE_NCBI":
                if line[1] not in sample_ncbi_dict:
                    sample_ncbi_dict[line[1]] = [line[2]]
                elif line[1] in sample_ncbi_dict:
                    sample_ncbi_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_ncbi_dict"] = sample_ncbi_dict
                config_info["sample_dict"] = sample_dict
            elif line[0] == "--SAMPLE_ENA":
                if line[1] not in sample_ena_dict:
                    sample_ena_dict[line[1]] = [line[2]]
                elif line[1] in sample_ena_dict:
                    sample_ena_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_ena_dict"] = sample_ena_dict
                config_info["sample_dict"] = sample_dict

            elif line[0] == "--GENOME_NCBI":
                config_info["genome_ncbi"] = line[1]
            elif line[0] == "--GENOME_LOCAL":
                config_info["genome_local"] = line[1]
            elif line[0] == "--SAMPLE_LOCAL":
                sys.exit("Local sample files are not supported at this time")
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
            elif line[0] == "--HETEROZYGOSITY":
                config_info["het"] = line[1]
            elif line[0] == "--PIPELINE":
                config_info["pipeline"] = line[1]
            elif line[0] == "--COVERAGE":
                config_info["coverage"] = line[1]
            elif line[0] == "--NINTERVALS":
                config_info["nintervals"] = line[1]
            elif line[0] == "--MEMORY_DS":
                config_info["memory_ds"] = line[1]
            elif line[0] == "--TIME_DS":
                config_info["time_ds"] = line[1]
            elif line[0] == "--MEMORY_HC":
                config_info["memory_hc"] = line[1]
            elif line[0] == "--TIME_HC":
                config_info["time_hc"] = line[1]
            elif line[0] == "--MEMORY_GG":
                config_info["memory_gg"] = line[1]
            elif line[0] == "--TIME_GG":
                config_info["time_gg"] = line[1]
    config_file.close()
    
    #Make sure all necessary inputs are present
    try:
        config_info["abbv"]
    except NameError:
        sys.exit("Oops, you forgot to specify a species abbreviation with --ABBV")
    
    try:
        config_info["out_dir"]
    except NameError:
        config_info["out_dir"] = "."
        
    if "genome_ncbi" not in config_info and "genome_local" not in config_info:
        sys.exit("Oops, you forgot to specify a reference genome!")
        
    if len(config_info["sample_dict"]) == 0:
        sys.exit("Oops, you forgot to specify samples!")
        
    if "het" not in config_info:
        config_info["het"] = "0.001"
        print("No heterozygosity specified, using default human value (0.001)")
    
    if "pipeline" not in config_info:
        config_info["pipeline"] = "lowcoverage"
        print("No pipeline specified (highcoverage or lowcoverage), using lowcoverage.")
        
    if config_info["pipeline"] != "highcoverage" and config_info["pipeline"] != "lowcoverage":
        sys.exit("Pipeline must be set to either 'highcoverage' or 'lowcoverage")
    
    if "coverage" not in config_info:
        config_info["coverage"] = "30"
        print("No desired coverage specified, using 30X")
    
    if "nintervals" not in config_info:
        config_info["nintervals"] = "10"
        print("No specification for the number of intervals to analyze, using 10")
    
    if "memory_ds" not in config_info:
        config_info["memory_ds"] = "8"
        print("No specification of how much memory to use for downsampling. Using 8GB by default")

    if "time_ds" not in config_info:
        config_info["time_ds"] = "8"
        print("No specification of how much time to use for downsampling, using 8 hours by default")
    
    if "memory_hc" not in config_info:
        config_info["memory_hc"] = "8"
        print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")
    
    if "time_hc" not in config_info:
        config_info["time_hc"] = "12"
        print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")
    
    if "memory_gg" not in config_info:
        config_info["memory_gg"] = "8"
        print("No specification of how much memory to use for GenotypeGVCF, using 8GB by default")
    
    if "time_gg" not in config_info:
        config_info["time_gg"] = "12"
        print("No specification of how much time to use for GenotypeGVCF, using 12 hours by default")
    
    #Return objects    
    return(config_info)


#Check if a directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):

        os.makedirs(dir)

#Create generic slurm script
def script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%j.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)
    
#Create generic slurm script for arrays (-o and -e for arrays)
def array_script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%A_%a.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%A_%a.err\n\n{cmd}
    '''
    return(slurm_script)

#Submit filename to slurm with sbatch with a given amount of time and memroy, returns job id number
def sbatch_submit(filename,memory,timelimit):
    proc = Popen('sbatch --mem %s --time %s:00:00 %s '%(memory,timelimit,filename),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    stdout = stdout.decode("utf-8","ignore")
    stdout = stdout.strip()
    stdout = stdout.strip('Submitted batch job ')
    return(stdout)

#Submit filename to slurm with sbatch with a given amount of time and memroy, returns job id number - includes first argument to include memory limit in java opts (use memory - 2)
def sbatch_submit_array(filename,memory,timelimit, array_nums):
    proc = Popen('sbatch --mem %s000 --time %s:00:00 --array=%s %s %d'%(memory,timelimit,array_nums,filename,(int(memory)-2)),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    stdout = stdout.decode("utf-8","ignore")
    stdout = stdout.strip()
    stdout = stdout.strip('Submitted batch job ')
    return(stdout)

#Check job status of specific jobid: returns job status
def jobid_status(jobid,date):
    proc = Popen('sacct --format state --noheader -j %d -S %s'%(int(jobid),date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split()
    return(lines[0].decode("utf-8","ignore"))
    
    
#Check the status of all jobs, returns dictionary of jobid:status. Ignores jobs with ".ba+" and ".ex+" file extensions.
def all_jobs_status(date):
    proc = Popen('sacct --format jobid,state --noheader -S %s'%(date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    stdout=stdout.decode("utf-8","ignore")
    stderr=stderr.decode("utf-8","ignore")
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split("\n")
    status_dict = {}
    for line in lines:
        line = line.split()
        if len(line) > 1:
            if not re.search('.ba+',line[0]) and not re.search('.ex+',line[0]):
                status_dict[line[0]]=line[1]
    return(status_dict)
    
#Count the number of pending jobs
def num_pend_run(job_id_list,date):
    count = 0
    status_dict = all_jobs_status(date)
    for job in job_id_list:
        if status_dict[job] == "PENDING" or status_dict[job] == "RUNNING":
            count += 1
    return(count)

#Check for missing gvcf interval files for a given sample in list of files. Will return a list of the missing intervals
def check_missing_gvcfs(arraystart,arrayend,sample_files,sample,coverage):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s.%sX.%s.g.vcf.gz"%(sample,coverage,str(i)) not in sample_files or "%s.%sX.%s.g.vcf.gz.tbi"%(sample,coverage,str(i)) not in sample_files:
            missing_ints.append(str(i))            
    return(missing_ints)
    
#Check for missing vcf interval files for a given sample in list of files. Will return a list of the missing intervals
def check_missing_vcfs(arraystart,arrayend,vcf_files,sp_abbr,coverage):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s.%sX.%s.vcf.gz"%(sp_abbr,coverage,str(i)) not in vcf_files or "%s.%sX.%s.vcf.gz.tbi"%(sp_abbr,coverage,str(i)) not in vcf_files:
            missing_ints.append(str(i))            
    return(missing_ints)


#Count number of intervals in split genome file
def count_intervals(nintervals,outputdir):
    if nintervals != "CHROMOSOME":
        nintervalfiles = nintervals
    else:
        filelist = os.listdir(outputdir)
        nintervalfiles = len(filelist)
    return(nintervalfiles)


###Create sbatch scripts
#Create a genotypegvcf sbatch file for a sample
def genotypegvcf_sbatch(sp_dir,sp_abbr,sample_list,coverage,het,nintervals,memory_gg):
    slurm_script = array_script_create()
    nintervals = str(nintervals)
    
    #Need to create string of --variant sample1.g.vcf.gz --variant sample2.g.vcf.gz for each sample in the list
    new_sample_list = []
    for i in range(0,len(sample_list)):
        new_sample_list.append('--variant %s/gvcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.g.vcf.gz'%(sp_dir,sample_list[i],coverage))

    all_sample_variant_call = " ".join(new_sample_list)

    
    #Load modules and get versions for all programs used
    ##For now, using my own installation of GATK as it is not yet installed on the cluster
    cmd_1 = 'module load java/1.8.0_45-fasrc01'
    
    cmd_2 = 'MEM=$1'    
    
    #Before running GenotypeGVCFs, need to combine GVCFs for all individuals into a single file      
    cmd_3 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" CombineGVCFs -R %s/genome/%s.fa %s -O %s/gvcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,all_sample_variant_call,sp_dir,sp_abbr,coverage,sp_dir,nintervals,sp_abbr)
    
    cmd_4 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenotypeGVCFs -R %s/genome/%s.fa -V %s/gvcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -O %s/vcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.vcf.gz --heterozygosity %s --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,sp_dir,sp_abbr,coverage,sp_dir,sp_abbr,coverage,het,sp_dir,nintervals,sp_abbr)
    
    cmd_5 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" VariantsToTable -V %s/vcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.vcf.gz -F CHROM -F POS -F TYPE -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F NCALLED -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -R %s/genome/%s.fa --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,coverage,sp_dir,sp_abbr,sp_dir,nintervals,sp_abbr)
        
    cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    genotypegvcf_script = slurm_script.format(partition="shared",cores="2",nodes="1",jobid="gg_%sX"%coverage,sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/07_genotypegvcf_%sX_array.sbatch"%(sp_dir,coverage)
    out_file = open(out_filename,"w")
    out_file.write(genotypegvcf_script)
    out_file.close

    return(out_filename)


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    now = datetime.datetime.now()
    print('Staring work on script 03: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get attributes
    config_info = extract_config(config_filename)

    #####Check if working directories exist, if not creates them
    if config_info["pipeline"] == "lowcoverage":
        sp_dir = "%s/%s_%sX_LC"%(config_info["out_dir"],config_info["abbv"],config_info["coverage"])
    elif config_info["pipeline"] == "highcoverage":
        sp_dir = "%s/%s_%sX_HC"%(config_info["out_dir"],config_info["abbv"],config_info["coverage"])
    
    print("\nOutput will be written to %s\n"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    dedup_dir = "%s/dedup"%(sp_dir)
    gvcf_dir = "%s/gvcf"%(sp_dir)
    vcf_dir = "%s/vcf"%(sp_dir)
       
    directory_create(sp_dir)    
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(genome_dir)
    directory_create(stats_dir)
    directory_create(dedup_dir)
    directory_create(gvcf_dir)
    directory_create(vcf_dir)


    #Need to count intervals from genome split in previous step. If the genome was split by a number of intervals, this is trivial (it is just the number), but if it was split by chromosome, need to count files.
        
    nintervalfiles = count_intervals(config_info["nintervals"],"%s/%s_splits_interval_lists/"%(genome_dir,config_info["nintervals"]))
    
    #####Run GenotypeGVCF

    #First, iterate through and make sure all HaplotypeCaller result files are present
    missing_dict = {}
    all_gvcfs = os.listdir(gvcf_dir)
    for sample in config_info["sample_dict"]:
        missing_ints = check_missing_gvcfs(1,int(nintervalfiles),all_gvcfs,sample,config_info["coverage"])
        if len(missing_ints) > 0:
            missing_dict[sample] = missing_ints
        else:
            pass
    
    #If there are any missing files, print statement and kill script
    if len(missing_dict) > 0:
        print('There are missing gvcf files for given samples and interval set:')
        for sample in missing_dict:
            missing = ", ".join(missing_dict[sample])
            print('%s: %s'%(sample,missing))
        sys.exit("Rerun missing intervals or change input samples")
    else:
        pass

    
    #Create GenotypeGVCF file - here we include all samples, so only a single slurm array is needed
    gg_filename = genotypegvcf_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample_list=list(config_info["sample_dict"].keys()),coverage=config_info["coverage"],het=config_info["het"],nintervals=config_info["nintervals"],memory_gg=config_info["memory_gg"])
    
    #Get number of finished files
    vcf_files = os.listdir(vcf_dir)
    finished_files = len([name for name in vcf_files if ".tbi" in name])
    
    #Submit file the first time
    if finished_files == 0:
        #Submit job
        base_jobid = sbatch_submit_array(gg_filename,memory=config_info["memory_gg"],timelimit=config_info["time_gg"], array_nums="1-%s"%str(nintervalfiles))

        #Expand job IDs
        all_jobids = []
        for i in range(1,int(nintervalfiles)+1):
            all_jobids.append("%s_%d"%(base_jobid,i))
    
    elif finished_files < nintervalfiles:
        #Check each interval, see if it has both a .vcf.gz and .tbi file
        
        missing = check_missing_vcfs(arraystart=1,arrayend=nintervalfiles,vcf_files=vcf_files,sp_abbr=config_info["abbv"],coverage=config_info["coverage"])
        missing_vec = ",".join(missing)
        
         #Submit job, get base jobid for array
        base_jobid = sbatch_submit_array(gg_filename,memory=config_info["memory_gg"],timelimit=config_info["time_gg"], array_nums=missing_vec)
        sleep(1)
    
        #Add jobids for array to dictionary with jobid as key and sample as value
        for i in missing:
            all_jobids.append("%s_%s"%(base_jobid,i))
            
    elif finished_files == nintervalfiles:
        sys.exit("All vcf files present, exiting")
        
    else:
        sys.exit("More vcf files present than expected, check")
        
    #Give sacct a chance to catch up       
    sleep(20)
    
    #Then, enter while loop that will continue until the number of completed jobs matches the number of intervals
    #Create dictionary of completed jobids and completion statuses
    completed_jobids = {}
    rerun_jobids = []
    successful_intervals = []
    failed_intervals = []
    
    while len(completed_jobids) < len(all_jobids):
        job_statuses = all_jobs_status(start_date)
        current_jobs = list(all_jobids.keys())
        for job in current_jobs:
            if job not in completed_jobids:
                if job in job_statuses:#Have to add this because array jobs may be delayed
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        array_id = job.split("_")[1]
                        
                        #If job_id is "COMPLETED", check to make sure both the .vcf.gz file and .tbi file are both present. If they are, print and add to successful_samples dictionary (sample:[intervals])
                        if job_statuses[job] == "COMPLETED":
                            if os.path.isfile("%s/vcf/%s.%sX.%s.vcf.gz"%(sp_dir,config_info["abbv"],config_info["coverage"],array_id)) and os.path.isfile("%s/vcf/%s.%sX.%s.vcf.gz.tbi"%(sp_dir,config_info["abbv"],config_info["coverage"],array_id)):
                                print("Job %s completed for interval %s"%(job, array_id))
                                successful_intervals.append(array_id)

                        #If job_id is not COMPLETED, it means there was some sort of failure in the job. Resubmit with 2x time (up to 7 days, or 168 hours) and 2x memory
                        elif job_statuses[job] != "COMPLETED" and job not in rerun_jobids:
                            new_mem = str(int(config_info["memory_hc"])*2)
                            new_time =  int(config_info["time_hc"])*2
                            if new_time > 168:
                               new_time = '168'
                            else:
                                new_time = str(new_time)
                            #Submit array with only that interval
                            resubmitted_jobid = sbatch_submit_array(gg_filename,memory=new_mem,timelimit=new_time, array_nums=array_id)
                            sleep(1)
                            
                            #Add job id (including array number) to both rerun_jobids and all_jobids
                            rerun_jobids.append('%s_%s'%(resubmitted_jobid,array_id))
                            
                            all_jobids.append('%s_%s'%(resubmitted_jobid,array_id))
                            
                            print("Job %s failed, retrying interval %s with %s memory and %s time"%(job,array_id,new_mem,new_time))
                        
                        #If just doesn't finish and already resubmitted, do not submit again, print failure to log file, and add to failed_intervals list
                        elif job_statuses[job] != "COMPLETED" and job in rerun_jobids:
                            print("GenotypeGVCF job %s failure 2x for and interval %s"%(array_id))
                            failed_intervals.append(array_id)
                                
                        else:
                            print("Error with GenotypeGVCF job checking and resubmissions")
                        
        sleep(30)
    
    #After all jobs have finished, report intervals failed twice
    failed_intervals = ",".join(failed_intervals)
    print("Failed for intervals: %s"%(failed_intervals))
     

    #Check that the final vcf is available, if so, remove intermediate files (combined gvcf)
    '''
    for i in range(1,nintervalsfiles+1):
        if os.path.isfile('%s/vcf/%s.%sX.%s.vcf.gz'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i))) and os.path.isfile('%s/vcf/%s.%sX.%s.vcf.gz.tbi'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i))):
            if os.path.isfile('%s/gvcf/%s.%sX.%s.g.vcf.gz'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i))):
                proc = Popen('rm %s/gvcf/%s.%sX.%s.g.vcf.gz'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i)),shell=True)
            if os.path.isfile('%s/gvcf/%s.%sX.%s.g.vcf.gz.tbi'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i))):
                proc = Popen('rm %s/gvcf/%s.%sX.%s.g.vcf.gz.tbi'%(sp_dir,config_info["abbv"],config_info["coverage"],str(i)),shell=True)
    '''         
    now = datetime.datetime.now()
    print('Finished script 04: %s'%now)

    
if __name__ == "__main__":
    main()