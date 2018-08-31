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
    sample_local_dict = {}
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
            elif line[0] == "--SAMPLE_LOCAL":
                if line[1] not in sample_local_dict:
                    sample_local_dict[line[1]] = [line[2]]
                elif line[1] in sample_local_dict:
                    sample_local_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_local_dict"] = sample_local_dict
                config_info["sample_dict"] = sample_dict


            elif line[0] == "--GENOME_NCBI":
                config_info["genome_ncbi"] = line[1]
            elif line[0] == "--GENOME_LOCAL":
                config_info["genome_local"] = line[1]
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
            elif line[0] == "--HETEROZYGOSITY":
                config_info["het"] = line[1]
            elif line[0] == "--PIPELINE":
                config_info["pipeline"] = line[1]
            elif line[0] == "--NINTERVALS":
                config_info["nintervals"] = line[1]
            elif line[0] == "--MEMORY_HC":
                config_info["memory_hc"] = line[1]
            elif line[0] == "--TIME_HC":
                config_info["time_hc"] = line[1]
            elif line[0] == "--MEMORY_GG":
                config_info["memory_gg"] = line[1]
            elif line[0] == "--TIME_GG":
                config_info["time_gg"] = line[1]
            elif line[0] == "--COMBINE_GVCF_PROGRAM":
                config_info["combine_gvcf_program"] = line[1]
            elif line[0] == "--BYPASS_INTERVAL":
                config_info["bypass_interval"] = line[1]
                
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
        #print("No heterozygosity specified, using default human value (0.001)")
    
    if "pipeline" not in config_info:
        config_info["pipeline"] = "lowcoverage"
        #print("No pipeline specified (highcoverage or lowcoverage), using lowcoverage.")
        
    #if config_info["pipeline"] != "highcoverage" and config_info["pipeline"] != "lowcoverage":
        #sys.exit("Pipeline must be set to either 'highcoverage' or 'lowcoverage")
    
    if "nintervals" not in config_info:
        config_info["nintervals"] = "10"
        #print("No specification for the number of intervals to analyze, using 10")
    
    if "memory_hc" not in config_info:
        config_info["memory_hc"] = "8"
        #print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")
    
    if "time_hc" not in config_info:
        config_info["time_hc"] = "12"
        #print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")
    
    if "memory_gg" not in config_info:
        config_info["memory_gg"] = "8"
        #print("No specification of how much memory to use for GenotypeGVCF, using 8 GB by default")
    
    if "time_gg" not in config_info:
        config_info["time_gg"] = "12"
        #print("No specification of how much time to use for GenotypeGVCF, using 12 hours by default")
        
    if "combine_gvcf_program" not in config_info:
        config_info["combine_gvcf_program"] = "GenomicsDBImport"
        #print("No specification of which program to use to combine gvcf files, using GenomicsDBImport by default")
        
    if "bypass_interval" not in config_info:
        config_info["bypass_interval"] = "FALSE"
        #print("BYPASS_INTERVAL set to FALSE, will require gvcfs from all samples for all intervals to proceed")
    
    
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
    proc = Popen('sbatch --mem %s000 --time %s:00:00 %s '%(memory,timelimit,filename),shell=True,stdout=PIPE,stderr=PIPE)
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
def check_missing_gvcfs(arraystart,arrayend,sample_files,sample):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s.%s.g.vcf.gz"%(sample,str(i)) not in sample_files or "%s.%s.g.vcf.gz.tbi"%(sample,str(i)) not in sample_files:
            missing_ints.append(str(i))            
    
    return(missing_ints)
    
#Check for missing vcf interval files for a given sample in list of files. Will return a list of the missing intervals
def check_missing_vcfs(arraystart,arrayend,vcf_files,sp_abbr):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s.%s.vcf.gz"%(sp_abbr,str(i)) not in vcf_files or "%s.%s.vcf.gz.tbi"%(sp_abbr,str(i)) not in vcf_files:
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



#Create an sbatch file for a given set of samples to create genome coverage graphs with bedtools - will use deduplicated BAM file and will write results to stats_coverage

def sample_coverage_sbatch(sp_dir,sp_abbr,sample):
    slurm_script = script_create()
     
    #First check if dedup file is present, if it is, continue and if not print statement.
    dedup_sorted_filename = '%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)
    genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
    if os.path.isfile(dedup_sorted_filename):
        if os.path.isfile(genome_cov_filename):
            print("Coverage file for sample %s already exists, will not recalculate"%(sample))
        else:
            #Load modules and get versions for all programs used
            cmd_1 = 'module load bedtools2/2.26.0-fasrc01'
        
            #Create bedgraph with bedtools of coverage (include regions with 0 coverage)
            cmd_2 = 'bedtools genomecov -bga -ibam %s/dedup/%s.dedup.sorted.bam -g %s/genome/%s.fa > %s/stats_coverage/%s.bg'%(sp_dir,sample,sp_dir,sp_abbr,sp_dir,sample)
    
    
            cmd_list = [cmd_1,cmd_2]

            final_cmd = "\n\n".join(cmd_list)


    #Format sbatch script
            sample_coverage_script = slurm_script.format(partition="shared",time="0-12:00",mem="8000",cores="1",nodes="1",jobid="genomecov",sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/07_coverage_bedgraph_%s.sbatch"%(sp_dir,sample)
            out_file = open(out_filename,"w")
            out_file.write(sample_coverage_script)
            out_file.close
            return(out_filename)

    else:
        print('No sorted dedup file for %s, cannot compute coverage bedgraph.'%(sample))
        
        
#Create an sbatch file for a given set of samples to create genome coverage graphs with bedtools - will use deduplicated BAM file and will write results to stats_coverage

def union_coverage_sbatch(sp_dir,sp_abbr,sample_bedgraph_file_list):
    slurm_script = script_create()
    
    sample_bedgraph_files = " ".join(sample_bedgraph_file_list)
    #First check if dedup file is present, if it is, continue and if not print statement.
            #Load modules and get versions for all programs used
    cmd_1 = 'module load bedtools2/2.26.0-fasrc01'

    #Create bedgraph with bedtools of coverage (include regions with 0 coverage)
    cmd_2 = 'bedtools unionbedg -header -empty -g %s/genome/%s.fa -i %s > %s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir,sp_abbr,sample_bedgraph_files,sp_dir,sp_abbr)


    cmd_list = [cmd_1,cmd_2]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    sample_coverage_script = slurm_script.format(partition="shared",time="0-12:00",mem="8000",cores="1",nodes="1",jobid="genomecov",sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/08_union_coverage_bedgraph_all_samples.sbatch"%(sp_dir,sample)
    out_file = open(out_filename,"w")
    out_file.write(sample_coverage_script)
    out_file.close
    return(out_filename)
    
#Takes in one line from a union coverage file, and returns a list with the [chromosome,start,end,interval_length,summed_coverage].
def compute_coverage_sum(union_cov_line):
    union_cov_line = union_cov_line.strip()
    union_cov_list = union_cov_line.split()
    
    #Grab positional info, calculate length of interval
    chrom = union_cov_list[0]
    start = union_cov_list[1]
    end = union_cov_list[2]
    interval_length = int(end) - int(start)
    
    #Turn coverage into integers, sum
    cov_list = [int(i) for i in union_cov_list[3:]]
    summed_coverage = sum(cov_list)
    
    return([chrom,start,end,interval_length,summed_coverage])


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    now = datetime.datetime.now()
    print('Staring work on script 05: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get attributes
    config_info = extract_config(config_filename)

    #####Check if working directories exist, if not creates them
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("\nOutput will be written to %s\n"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    dedup_dir = "%s/dedup"%(sp_dir)
    stats_coverage_dir = "%s/stats_coverage"%(sp_dir)
    
       
    directory_create(sp_dir)    
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(stats_coverage_dir)

    
    
    #Submit all jobs
    #Coverage_filenames is a dictionary with the sample as key and filename as value
    coverage_filenames = {}
    #all_jobids is a dictionary with jobid (including array numbers as key and sample as value)
    all_jobids = {}

    for sample in config_info["sample_dict"]:
        coverage_filenames[sample] = sample_coverage_sbatch(sp_dir,config_info["abbv"],sample)
    if len(coverage_filenames) > 0:
        #Submit dedup read sbatch files
        coverage_jobids = {}
        completed_jobids = {}
        for sample in coverage_filenames:
            all_jobids[sbatch_submit(coverage_filenames[sample],memory=8,timelimit=8)] = sample
            sleep(1)
        #Add an extra sleep to give sacct a chance to catch up
        sleep(20)
        #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
        while len(completed_jobids) < len(coverage_filenames):
            num_running = num_pend_run(all_jobids,start_date)
            job_statuses = all_jobs_status(start_date)
            for job in all_jobids:
                if job not in completed_jobids:
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        print("Job %s completed for sample %s"%(job,all_jobids[job]))
            sleep(30)
    
        #After all jobs have finished, report which jobs failed, also add samples that failed to list to exclude from whole-genome calculations
        failed_samples = []
        for job in completed_jobids:
            if completed_jobids[job] != "COMPLETED":
                print("Coverage bedgraph job %s failed with code: %s for sample %s"%(job,completed_jobids[job],all_jobids[job]))
                failed_samples.append(all_jobids[job])


    print("Now computing whole-genome coverage bedgraph")
    
    ########## Create union bedgraph with all sample coverages included
    #Make a list of all sample bedgraphs to include in genome coverage calculations. Iterates through samples names, makes sure coverage bedgraph files exist (if not, skips sample), and that sample completed above command successfully.
    sample_bedgraph_file_list = []
    for sample in config_info["sample_dict"]:
        genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
        if sample in failed_samples:
            print("Will not include sample %s in genome coverage calculation, job failed")
        elif os.path.isfile(genome_cov_filename):
            sample_bedgraph_file_list.append(genome_cov_filename)
    
    #Create and submit file for union bedgraph job        
    union_sbatch_file = union_coverage_sbatch(sp_dir,config_info["abbv"],sample_bedgraph_file_list)
    union_job_id = sbatch_submit(union_sbatch_file,memory=8,timelimit=8)
    sleep(30)

    #Only check on union job if actually submitted. 
    if union_job_id is not None: 
        dones = ['COMPLETED','CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL']
        #Check job id status of union job. If not in one of the 'done' job status categories, wait 30 seconds and check again.
        while jobid_status(union_job_id,start_date) not in dones:
            sleep(30)
    
        #Check to make sure job completed, and that all necessary files are present. If not, exit and give information.
        union_job_completion_status = jobid_status(union_job_id,start_date)
        if union_job_completion_status != 'COMPLETED':
            sys.exit("There was a problem creating the union bedgraph file. The job exited with status %s. Please diagnose and fix before moving on"%union_job_completion_status)
            
            
    #Read through resulting bedgraph file, create new bedgraph from sum across all sample coverages, and compute median
    #Dictionary to create histogram as we iterate each line
    coverage_histogram = {}
    total_sites = 0
    union_cov_filename = '%s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir, config_info["abbv"])
    summary_bedgraph = open('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,config_info["abbv"]), 'w')
    
    union_cov_file = open(union_cov_filename,'r')
    for line in union_cov_file:
        summed_cov = compute_coverage_sum(line)
        #Returns [chrom,start,end,interval_length,summed_coverage]

        #Write to new summary bedgraph
        summary_bedgraph.write('%s/t%s/t%s/t%d\n'%(summed_cov[0],summed_cov[1],summed_cov[2],summed_cov[4]))
        #Add coverage to histogram dictionary
        if summed_cov[4] not in coverage_histogram:
            coverage_histogram[summed_cov[4]] = summed_cov[3]
        else:
            coverage_histogram[summed_cov[4]] = coverage_histogram[summed_cov[4]] + summed_cov[3]
        
        #Add interval length to total sites
        total_sites += summed_cov[3]
 
    union_cov_file.close()
    summary_bedgraph.close()
       
    #Order and write histogram to file
    #Calculate median value
    ordered_hist_bins = sorted(coverage_histogram.keys())
    summed_hist_file = open('%s/stats_coverage/_%s_all_samples_summed_cov_histogram.txt','w')
    summed_hist_file.write('SUMMED_COVERAGE\tN_SITES\n')
    median_cov = None
    median_cov_cutoff = (total_sites+1)/2
    
    for bin in ordered_hist_bins:
        summed_hist_file.write('%d\t%d\n'%(bin,coverage_histogram[bin]))
        if median_cov is None and coverage_histogram[bin] > median_cov_cutoff:
            median_cov = bin
        else:
            pass
    
    #Calculate upper and lower limits of coverage and print to log file
    upper_lim_cov = median_cov * 2
    lower_lim_cov = median_cov * 0.5
    print('\nMedian coverage of all samples is %d, will exclude all sites > %d and < %d'%(median_cov,upper_lim_cov,lower_lim_cov))
    
    #Create new bedfile of sites < 2X median and > 0.5X median
    summary_bedgraph = open('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,config_info["abbv"]), 'r')
    clean_cov_bedfile = open('%s/stats_coverage/_%s_clean_coverage_sites.bed','w')
    low_cov_bedfile = open('%s/stats_coverage/_%s_too_low_coverage_sites.bed','w')
    high_cov_bedfile = open('%s/stats_coverage/_%s_too_high_coverage_sites.bed','w')
    for line in summary_bedgraph:
        line = line.strip()
        split_line = line.split()
        if split_line[3] > lower_lim_cov and split_line[3] < upper_lim_cov:
            clean_cov_bedfile.write('%s\t%s\t%s\n'%(split_line[0],split_line[1],split_line[2]))
        elif split_line[3] < lower_lim_cov:
            low_cov_bedfile.write('%s\t%s\t%s\t%d\n'%(split_line[0],split_line[1],split_line[2],split_line[3]))
        elif split_line[3] > upper_lim_cov:
            high_cov_bedfile.write('%s\t%s\t%s\t%d\n'%(split_line[0],split_line[1],split_line[2],split_line[3]))
        else:
            print("Something wrong with coverage calculation")


    summary_bedgraph.close()
    summed_hist_file.close()
    clean_cov_bedfile.close()
    low_cov_bedfile.close()
    high_cov_bedfile.close()
    
    now = datetime.datetime.now()
    print('Finished script 05: %s'%now)


if __name__ == "__main__":
    main()