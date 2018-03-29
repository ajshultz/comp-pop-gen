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
            elif line[0] == "--MEMORY_HC":
                config_info["memory_hc"] = line[1]
            elif line[0] == "--TIME_HC":
                config_info["time_hc"] = line[1]
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
    
    if "memory_hc" not in config_info:
        config_info["memory_hc"] = "8"
        print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")
    
    if "time_hc" not in config_info:
        config_info["time_hc"] = "12"
        print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")
    
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
    proc = Popen('sbatch --mem %s --time %s %s '%(memory,timelimit,filename),shell=True,stdout=PIPE,stderr=PIPE)
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
        if "%s.%sX.%s.g.vcf.gz"%(sample,coverage,str(i)) not in sample_files and "%s.%sX.%s.g.vcf.gz.tbi"%(sample,coverage,str(i)) not in sample_files:
            missing_ints.append(str(i))            
    
    return(missing_ints)


def split_genome(sp_dir,sp_abbr,nintervals,outputdir):
	#Open input		
    fai = open("%s/genome/%s.fa.fai"%(sp_dir,sp_abbr),"r")
    outHandle = sp_abbr

    if nintervals != "CHROMOSOME":
        faiList = []
        totLen = 1
        cumStart = []
        cumEnd = []
    
        n = int(nintervals)
    
        #For each line in the .fai file, append data to a list, and create two additional lists with cumulative start and end positions.
        for line in fai:
            line = line.strip().split("\t")
            faiList.append(line)
            cumStart.append(totLen)
            totLen = totLen + int(line[1])
            cumEnd.append(totLen)

        #Create interval based on total length of scaffolds, and create list of start and end values
        interval = (totLen/n)+1
        intervalNums = range(1,n+1)
        intervalEnd = [x * interval for x in intervalNums]
        intervalStart = [(x - interval) + 1 for x in intervalEnd]
    
        #Create list of which file each scaffold should go to, depending on where it falls on the interval list. Note that scaffolds that span intervals will be moved to the preceeding interval.
        fileDes = []
        for i in range(len(faiList)):
            fileDes.append(0)
            for j in range(n):
                if cumStart[i] >= intervalStart[j] and cumEnd[i] < intervalEnd[j]:
                    fileDes[i] = (j+1)
                elif j < (n-1) and cumStart[i] >= intervalStart[j] and cumEnd[i] < intervalEnd[(j+1)] and cumStart[i] < intervalStart[(j+1)]:
                    fileDes[i] = (j+1)
    
        #Create list of output files based on species abbreviation
        outList = []
    
        for i in range(n):
            outFile = (outputdir+outHandle+"_"+(str(i+1))+".interval_list")
            outList.append(open(outFile,"w"))
        
        #Append scaffold name to appropriate output file.						
        for i in range(len(faiList)):
            outList[fileDes[(i)]-1].write(faiList[i][0]+"\n")
        
        for i in range(n):
            outList[i].close()
        
        nintervalfiles = int(nintervals)
    
    else:
        faiList = []
        for line in fai:
            line = line.strip().split("\t")
            faiList.append(line)
        #Get number of chromosomes, will produce that many interval files + 1
        nchr = len([name for name in faiList if "NC_" in name[0]])
        nintervalfiles=(nchr+1)
        filenum = 1
        #Open file to write all non-chromosomes
        randoutFile = open((outputdir+outHandle+"_"+(str(nintervalfiles))+".interval_list"),"w")
        for i in range(len(faiList)):
            if "NC_" in faiList[i][0]:
                outFile = open((outputdir+outHandle+"_"+(str(filenum))+".interval_list"),"w")
                outFile.write(faiList[i][0])
                outFile.close()
                filenum += 1
            else:
                randoutFile.write("%s\n"%faiList[i][0])
                
        randoutFile.close()

    fai.close()
	
    return(nintervalfiles)



###Create sbatch scripts
#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def downsample_sbatch(sp_dir,sp_abbr,sample_dict,coverage,coverage_dict,memory_ds):
    slurm_script = script_create()
    downsample_sbatch_filenames = []
    
    for sample in sample_dict.keys():
        #First check if dedup file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
        dedup_filename = '%s/dedup/%s.%sX.dedup.bam'%(sp_dir,sample,coverage)
        dedup_sorted_filename = '%s/dedup/%s.%sX.dedup.sorted.bam'%(sp_dir,sample,coverage)
        if os.path.isfile(dedup_filename) or os.path.isfile(dedup_sorted_filename):
            print('%s.%sX.dedup.bam already present, skipping'%(sample,coverage))
        else:
            print('Will donwsample BAM for sample %s'%(sample))

            #Original directory with full BAMs
            orig_dir = '../SPECIES_DATASETS/%s/dedup'%(sp_abbr)

            #Load modules and get versions for all programs used
            ##For now, using my own installation of GATK as it is not yet installed on the cluster
            cmd_1 = 'module load java/1.8.0_45-fasrc01'
            
            #Command to donwsample if proportion <0.95, if >0.95, just copy
            if coverage_dict[sample] < 0.95:
                cmd_2 = 'gatk --java-options "-Xmx%dg -XX:ParallelGCThreads=1" DownsampleSam -I %s/%s.dedup.sorted.bam -O %s/dedup/%s.%sX.dedup.bam -S ConstantMemory -P %f -A 0.0001 --COMPRESSION_LEVEL 5'%((int(memory_ds)-2),orig_dir,sample,sp_dir,sample,coverage,coverage_dict[sample])
                
                #Sort and index dedup bam file
                cmd_3 = 'gatk --java-options "-Xmx%dg -XX:ParallelGCThreads=1" SortSam -I %s/dedup/%s.%sX.dedup.bam -O %s/dedup/%s.%sX.dedup.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --COMPRESSION_LEVEL 5'%((int(memory_ds)-2),sp_dir,sample,coverage,sp_dir,sample,coverage)
            
            else:
                cmd_2 = "cp %s/%s.dedup.sorted.bam %s/dedup/%s.%sX.dedup.sorted.bam"%(orig_dir,sample,sp_dir,sample,coverage)
                
                #Re-create index
                cmd_3 = 'gatk --java-options "-Xmx%dg -XX:ParallelGCThreads=1" BuildBamIndex -I %s/dedup/%s.%sX.dedup.sorted.bam'%((int(memory_ds)-2),sp_dir,sample,coverage)
             
            #Compute coverage histogram of sorted bam
            cmd_4 = 'bedtools genomecov -ibam %s/dedup/%s.%sX.dedup.sorted.bam -g %s/genome/%s.fa > %s/stats/%s.%sX.coverage'%(sp_dir,sample,coverage,sp_dir,sp_abbr,sp_dir,sample,coverage)
            
            #Grab only genome output
            cmd_5 = r"""awk '$1 == "genome" {print $0}' %s/stats/%s.%sX.coverage > %s/stats/%s.%sX.genome.coverage"""%(sp_dir,sample,coverage,sp_dir,sample,coverage)
    
            cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5]

            final_cmd = "\n\n".join(cmd_list)

            #Format sbatch script
            downsample_script = slurm_script.format(partition="shared",cores="2",nodes="1",jobid="ds_%sX"%coverage,sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/05_downsample_%sX_%s.sbatch"%(sp_dir,coverage,sample)
            out_file = open(out_filename,"w")
            out_file.write(downsample_script)
            out_file.close
            downsample_sbatch_filenames.append(out_filename)
    
    return(downsample_sbatch_filenames)

#Create a haplotypecaller sbatch file for a sample
def haplotypecaller_sbatch(sp_dir,sp_abbr,sample,coverage,het,memory_hc,nintervals,pipeline):
    slurm_script = array_script_create()
    nintervals = str(nintervals)

    #Load modules and get versions for all programs used
    ##For now, using my own installation of GATK as it is not yet installed on the cluster
    cmd_1 = 'module load java/1.8.0_45-fasrc01'
    
    cmd_2 = 'MEM=$1'
    
    if pipeline == "highcoverage":
    #Command to donwsample if proportion <0.95, if >0.95, just copy
        cmd_3 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" HaplotypeCaller -I %s/dedup/%s.%sX.dedup.sorted.bam -O %s/gvcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -R %s/genome/%s.fa --heterozygosity %s --ERC GVCF --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sample,coverage,sp_dir,sample,coverage,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    
    elif pipeline == "lowcoverage":
        cmd_3 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" HaplotypeCaller -I %s/dedup/%s.%sX.dedup.sorted.bam -O %s/gvcf/%s.%sX.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -R %s/genome/%s.fa --heterozygosity %s --ERC GVCF --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list --minDanglingBranchLength 1 --minPruning 1'%(sp_dir,sample,coverage,sp_dir,sample,coverage,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    
    cmd_list = [cmd_1,cmd_2,cmd_3]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    haplocaller_script = slurm_script.format(partition="shared",cores="2",nodes="1",jobid="hc_%sX"%coverage,sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/06_haplotypecaller_%sX_%s_array.sbatch"%(sp_dir,coverage,sample)
    out_file = open(out_filename,"w")
    out_file.write(haplocaller_script)
    out_file.close

    return(out_filename)



#Given a sample name, reads in histogram of genome coverage, calculates mean and median, and returns dictionary with coverage histogram, mean and median for that sample
def collect_coverage_metrics(sp_dir,sample,coverage):
    coverage_file = '%s/stats/%s.%sX.genome.coverage'%(sp_dir,sample,coverage)
    coverage_stats = {}
    hist_bins = []
    hist_vals = []
    site_total = 0
    cov_total = 0
    try:
        full_coverage_stats = open(coverage_file,"r")
        for line in full_coverage_stats:
            split_line = line.strip().split("\t")
            med_cutoff = (int(split_line[3])+1)/2
            site_total += int(split_line[2])
            cov_total += (int(split_line[1])*int(split_line[2]))
            
            hist_bins.append(int(split_line[1]))
            hist_vals.append(float(split_line[4]))
            
            #Test whether site total is greater than midpoint if median not yet discovered - if yes, add to median in dictionary.
            if "median" not in coverage_stats and site_total > med_cutoff:
                coverage_stats["median"] = int(split_line[1])
            else:
                pass
        
        #Calculate mean from cov_total and site_total
        coverage_stats["mean"] = round(cov_total/site_total,2)
        
        #Add histogram to dictionary
        coverage_stats["hist_vals"] = hist_vals
        coverage_stats["hist_bins"] = hist_bins
    except:
        print("No genome coverage file for sample: %s"%(sample))
        
    return(coverage_stats)
    


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
    
    sp_dir = "%s/%s_%sX"%(config_info["out_dir"],config_info["abbv"],config_info["coverage"])
    
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


    #####Get genome file from previous runs, copy to this genome directory, create interval file to split up jobs; create jobs to subset original dedup files to desired coverage (based on median original coverage) **Assumes working in directory on same level as SPECIES_DATASETS and comp-pop-gen

    #Copy genome + fai and dict file from original directory to new directory
    orig_genome = "../SPECIES_DATASETS/%s/genome/%s"%(config_info["abbv"],config_info["abbv"])
    
    if not os.path.isfile('%s/genome/%s.fa'%(sp_dir,config_info["abbv"])):
        proc = Popen('cp %s.fa %s'%(orig_genome,genome_dir),shell=True)
    if not os.path.isfile('%s/genome/%s.fa.fai'%(sp_dir,config_info["abbv"])):
        proc = Popen('cp %s.fa.fai %s'%(orig_genome,genome_dir),shell=True)
        sleep(60)
    if not os.path.isfile('%s/genome/%s.dict'%(sp_dir,config_info["abbv"])):
        proc = Popen('cp %s.dict %s'%(orig_genome,genome_dir),shell=True)
        
    #Split into N specified intervals (if CHROMOSOME instead of number, will split into chromosomes, with all unplaced scaffolds in one file), create directory to store interval files. Function returns the number of interval files (nintervalfiles)
    #Added sleep to give proc time to copy the .fai file

    directory_create('%s/%s_splits_interval_lists/'%(genome_dir,config_info["nintervals"]))
    nintervalfiles = split_genome(sp_dir,config_info["abbv"],config_info["nintervals"],"%s/%s_splits_interval_lists/"%(genome_dir,config_info["nintervals"]))
    
    #Read in file with coverage info, create a dictionary with sample name as the key, and the proportion to downsample, calculated as the desired coverage/median coverage
    
    #test_sample_info_file = open("../comp-pop-gen/variant_calling_testing/test_sample_info.csv","r")
    test_sample_info_file = open("test_sample_info.csv","r")
    
    test_sample_info_dict={}
    
    for line in test_sample_info_file:
        line = line.strip().split(",")
        if line[0] == config_info["abbv"]:
            test_sample_info_dict[line[1]] = (float(config_info["coverage"])/int(line[3]))

    #####Downsampling files
    
    #Create sbatch script to grab deduped BAM files, downsample to desired proportion (just copy if proportion < 0.95)
    downsample_filenames = downsample_sbatch(sp_dir,sp_abbr = config_info["abbv"],sample_dict = config_info["sample_dict"],coverage = config_info["coverage"],coverage_dict = test_sample_info_dict,memory_ds = config_info["memory_ds"])

    #Submit sbatch files, including memory and time requirements
    downsample_jobids = []
    completed_jobids = {}
    for i in range(0,len(downsample_filenames)):
        downsample_jobids.append(sbatch_submit(downsample_filenames[i],memory="%s000"%config_info["memory_ds"],timelimit="0-08:00"))
        sleep(1)
    #Add an extra sleep to give sacct a chance to catch up
    sleep(20)
    #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
    while len(completed_jobids) < len(downsample_jobids):
        num_running = num_pend_run(downsample_jobids,start_date)
        job_statuses = all_jobs_status(start_date)
        for job in downsample_jobids:
            if job not in completed_jobids:
                if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                    completed_jobids[job] = job_statuses[job]
                    print("Job %s completed"%job)
        sleep(30)
    
    #After all jobs have finished, report which jobs failed
    for job in completed_jobids:
        if completed_jobids[job] != "COMPLETED":
            print("Downsample %s failed with code: %s"%(job,completed_jobids[job]))
 
 
    #####Collate coverage summary statistics from downsampled files
    #Get stats directory filenames
    stat_files = os.listdir("%s/stats"%(sp_dir))
    
    #Collect coverage histograms, calculate mean and median coverage
    all_coverage_stats = {}
    for sample in config_info["sample_dict"]:
        all_coverage_stats[sample] = collect_coverage_metrics(sp_dir,sample,coverage=config_info["coverage"])   

    summary_stat_file = '%s/stats/_%s_downsampled_%sX_coverage_stats.txt'%(sp_dir,config_info["abbv"],config_info["coverage"])
    sum_stat = open(summary_stat_file,"w")
    
    #Write header
    sum_stat.write("SAMPLE\tMEAN_COVERAGE\tMEDIAN_COVERAGE\n")
    
    #Iterate through samples and add results from all three dictionaries if present
    for sample in config_info["sample_dict"]:
        samp_sum = []
        samp_sum.append(sample)
        #coverage
        if "mean" in all_coverage_stats[sample]:
            samp_sum.append(str(all_coverage_stats[sample]["mean"]))
            samp_sum.append(str(all_coverage_stats[sample]["median"]))
        else:
            samp_sum.append("")
            samp_sum.append("")            
        samp_sum = "\t".join(samp_sum)
        sum_stat.write('%s\n'%samp_sum)
    sum_stat.close()
   
    #Copy this file to centralized location        
    general_dir = "_ALL_SPECIES_SUMMARIES"
    directory_create(general_dir)
    
    try:
        proc = Popen('cp %s %s/%s_all_summary_stats.txt'%(summary_stat_file,general_dir,config_info["abbv"]),shell=True)
    except:
        print("There was an error copying summary stat files")
    

    #####Run HaplotypeCaller
    
    #Submit all jobs the first time
    #hc_filenames is a dictionary with the sample as key and filename as value
    hc_filenames = {}
    #all_jobids is a dictionary with jobid (including array numbers as key and sample as value)
    all_jobids = {}

    for sample in config_info["sample_dict"]:
        sample_files = [name for name in os.listdir(gvcf_dir) if sample in name]
        
        #If no files exist, submit full array
        if len(sample_files) == 0:

            #Create sbatch file, add to filename dictionary with sample as key and filename as value
            hc_filename = haplotypecaller_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample=sample,coverage=config_info["coverage"],het=config_info["het"],memory_hc=config_info["memory_hc"],nintervals=nintervalfiles,pipeline=config_info["pipeline"])
            hc_filenames[sample] = hc_filename
        
            #Submit job, get base jobid for array
            base_jobid = sbatch_submit_array(hc_filename,memory=config_info["memory_hc"],timelimit=config_info["time_hc"], array_nums="1-%d"%nintervalfiles)
            sleep(1)
        
            #Add jobids for array to dictionary with jobid as key and sample as value
            for i in range(1,nintervalfiles+1):
                all_jobids["%s_%d"%(base_jobid,i)] = sample
        
        #If the number of sample files is less than the the number of interval files x2 (because of vcf and index), that means some intervals are missing. Only submit those intervals that don't have .tbi (index) files.
        elif len(sample_files) < 2*nintervalfiles:
            #Check each interval, see if it has both a .vcf.gz and .tbi file
            hc_filename = haplotypecaller_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample=sample,coverage=config_info["coverage"],het=config_info["het"],memory_hc=config_info["memory_hc"],nintervals=nintervalfiles,pipeline=config_info["pipeline"])
            hc_filenames[sample] = hc_filename
            
            missing = check_missing_gvcfs(arraystart=1,arrayend=nintervalfiles,sample_files=sample_files,sample=sample,coverage=config_info["coverage"])
            
            print(missing)
            
            missing_vec = ",".join(missing)
            
             #Submit job, get base jobid for array
            base_jobid = sbatch_submit_array(hc_filename,memory=config_info["memory_hc"],timelimit=config_info["time_hc"], array_nums=missing_vec)
            sleep(1)
        
            #Add jobids for array to dictionary with jobid as key and sample as value
            for i in missing:
                all_jobids["%s_%s"%(base_jobid,i)] = sample
            
        elif len(sample_files) == 2*nintervalfiles:
            print("Sample %s has all gvcf files, skipping HaplotypeCaller"%sample)
        
        else:
            print("Sample %s has more gvcfs than expected, check"%sample)
            
    
    #Give sacct a chance to catch up       
    sleep(20)
    
    #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
    #Create dictionary of completed jobids and completion statuses
    completed_jobids = {}
    rerun_jobids = {}
    successful_samples = {}
    failed_samples = {}
    
    while len(completed_jobids) < len(all_jobids):
        job_statuses = all_jobs_status(start_date)
        for job in all_jobids:
            if job not in completed_jobids:
                if job in job_statuses:#Have to add this because array jobs may be delayed
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        array_id = "_".split(job)[1]
                        
                        #If job_id is "COMPLETED", check to make sure both the .vcf.gz file and .tbi file are both present. If they are, print and add to successful_samples dictionary (sample:[intervals])
                        if job_statuses[job] == "COMPLETED":
                            if os.path.isfile("%s/gvcf/%s.%sX.%s.g.vcf.gz"%(sp_dir,all_jobids[job],config_info["coverage"])) and os.path.isfile("%s/gvcf/%s.%sX.%s.g.vcf.gz.tbi"%(sp_dir,all_jobids[job],config_info["coverage"])):
                                print("Job %s completed for sample %s"%(job, all_jobids[job]))
                                if all_jobids[job] in successful_samples:
                                    successful_sample[all_jobids[job]].append(array_id)
                                else:
                                    successful_sample[all_jobids[job]] = [array_id]
                        #If job_id is not COMPLETED, it means there was some sort of failure in the job. Resubmit with 2x time (up to 7 days, or 168 hours) and 2x memory
                        elif job_statuses[job] != "COMPLETED" and job not in rerun_jobids:
                            new_mem = str(int(config_info["memory_hc"])*2)
                            new_time =  int(config_info["time_hc"])*2
                            if new_time > 168:
                               new_time = '168'
                            else:
                                new_time = str(new_time)
                            #Submit array with only that interval
                            resubmitted_jobid = sbatch_submit_array(hc_filenames[all_jobids[job]],memory=new_mem,timelimit=new_time, array_nums=array_id)
                            sleep(1)
                            
                            #Add job id (including array number) to both rerun_jobids and all_jobids
                            rerun_jobids['%s_%s'%(resubmitted_jobid,array_id)] = [all_jobids[job]]
                            
                            all_jobids['%s_%s'%(resubmitted_jobid,array_id)] = [all_jobids[job]]
                        
                        #If just doesn't finished and already resubmitted, do not submit again, print failure to log file, and add to failed_samples dictionary
                        elif job_statuses[job] != "COMPLETED" and job in rerun_jobids:
                            print("HaplotypeCaller failure 2x for sample %s and interval %s"%(all_jobids[job],array_id))    
                            if all_jobids[job] in failed_samples:
                                failed_samples[all_jobids[job]].append(array_id)
                            else:
                                failed_samples[all_jobids[job]] = [array_id]
                                
                        else:
                            print("Error with HaplotypeCaller job checking and resubmissions")
                        
        sleep(30)
    
    #After all jobs have finished, report which samples and intervals failed twice
    for sample in failed_samples:
        failed_intervals = ",".join(failed_samples[sample])
        print("Sample %s, failed for intervals: %s"%(sample,failed_intervals))
     

    #Check that the final downsampled sorted bam and index are available, if so, remove intermediate files (unsorted dedup)
    for sample in config_info["sample_dict"]:
        if os.path.isfile('%s/dedup/%s.%sX.dedup.sorted.bam'%(sp_dir,sample,config_info["coverage"])) and os.path.isfile('%s/dedup/%s.%sX.dedup.sorted.bai'%(sp_dir,sample,config_info["coverage"])):
            if os.path.isfile('%s/dedup/%s.%sX.dedup.bam'%(sp_dir,sample,config_info["coverage"])):
                proc = Popen('%s/dedup/%s.%sX.dedup.bam'%(sp_dir,sample,config_info["coverage"]),shell=True)
               
    now = datetime.datetime.now()
    print('Finished script 03: %s'%now)


if __name__ == "__main__":
    main()