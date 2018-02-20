#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import datetime

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename):
    print("Opening %s"%config_filename)
    config_file = open(config_filename,"r")
    config_info = {}
    sample_ncbi_dict = {}
    
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
                config_info["sample_ncbi_dict"] = sample_ncbi_dict
            elif line[0] == "--GENOME_NCBI":
                config_info["genome_ncbi"] = line[1]
            elif line[0] == "--GENOME_LOCAL":
                config_info["genome_local"] = line[1]
            elif line[0] == "--SAMPLE_LOCAL":
                sys.exit("Local sample files are not supported at this time")
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
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
        
    if len(config_info["sample_ncbi_dict"]) == 0:
        sys.exit("Oops, you forgot to specify samples!")
    
    #Return objects    
    return(config_info)


#Check if a directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):

        os.makedirs(dir)

#Create generic slurm script
def script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -t {time}\n#SBATCH --mem {mem}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%j.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)

    
#Submit filename to slurm with sbatch, returns job id number
def sbatch_submit(filename):
    proc = Popen('sbatch %s'%filename,shell=True,stdout=PIPE,stderr=PIPE)
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
    

def num_pend_run(job_id_list,date):
    count = 0
    status_dict = all_jobs_status(date)
    for job in job_id_list:
        if status_dict[job] == "PENDING" or status_dict[job] == "RUNNING":
            count += 1
    return(count)
    
    
##########Above functions same as script 01, may just wish to source that script in future.


#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def dedup_sbatch(sp_dir,sp_abbr,sample_ncbi_dict):
    slurm_script = script_create()
    dedup_sbatch_filenames = []
     
    for sample in sample_ncbi_dict.keys():
        #First check if dedup file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
        dedup_filename = '%s/dedup/%s.dedup.bam'%(sp_dir,sample)
        if os.path.isfile(dedup_filename):
            print('%s.dedup.bam already present, skipping'%(sample))
        else:
            print('Will dedup sras for sample %s'%(sample))
            

            align_dir = '%s/alignment/'%(sp_dir)

            #Load modules and get versions for all programs used
            ##For now, using my own installation of GATK as it is not yet installed on the cluster
            cmd_1 = 'module load java/1.8.0_45-fasrc01'
            

            #Add PL read group info (can remove once all initial libraries done with script 2)
            rg_cmd_list = []
            for sra in sample_ncbi_dict[sample]:
                
                rg_cmd = 'gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" AddOrReplaceReadGroups -I %s/alignment/%s.sorted.bam -O %s/alignment/%s.sorted.rg.bam -ID %s -SM %s -PU %s.%s -LB %s -PL illumina --COMPRESSION_LEVEL 5 --CREATE_INDEX true'%(sp_dir,sra,sp_dir,sra,sra,sample,sra,sample,sample)
                
                rg_cmd_list.append(rg_cmd)
            
            rg_cmd = "\n\n".join(rg_cmd_list)
            
            #Create GATK mark duplicates command
            gatk_cmd_1 = 'gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" MarkDuplicatesGATK '
            gatk_cmd_2 = ('-I '+align_dir+'%s.sorted.rg.bam ')*len(sample_ncbi_dict[sample])%tuple(sample_ncbi_dict[sample])
            gatk_cmd_3 = '-O %s/dedup/%s.dedup.bam '%(sp_dir,sample)
            gatk_cmd_4 = '--METRICS_FILE %s/stats/%s.dedup.metrics.txt '%(sp_dir,sample)
            gatk_cmd_5 = '--COMPRESSION_LEVEL 5'
            
            #Combine into single line
            cmd_2 = "".join([gatk_cmd_1,gatk_cmd_2,gatk_cmd_3,gatk_cmd_4,gatk_cmd_5])
            
            #Sort and index dedup bam file
            cmd_3 = 'gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" SortSam -I %s/dedup/%s.dedup.bam -O %s/dedup/%s.dedup.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --COMPRESSION_LEVEL 5'%(sp_dir,sample,sp_dir,sample)
            
            cmd_4 = 'gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" CollectAlignmentSummaryMetrics -I %s/dedup/%s.dedup.sorted.bam -R %s/genome/%s.fa --METRIC_ACCUMULATION_LEVEL=SAMPLE -O %s/stats/%s.alignment_metrics.txt'%(sp_dir,sample,sp_dir,sp_abbr,sp_dir,sample)
            
            #Validate sorted bam
            cmd_5 = 'gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" ValidateSamFile -I %s/dedup/%s.dedup.sorted.bam -O %s/stats/%s.validate.txt'%(sp_dir,sample,sp_dir,sample)
            
            #Compute coverage histogram of sorted bam
            cmd_6 = 'bedtools genomecov -ibam %s/dedup/%s.dedup.sorted.bam -g %s/genome/%s.fa > %s/stats/%s.coverage'%(sp_dir,sample,sp_dir,sp_abbr,sp_abbr,sample)
            
            #Grab only genome output
            cmd_7 = r"""awk '$1 == "genome" {print $0}' %s/stats/%s.coverage > %s/stats/%s.genome.coverage"""%(sp_dir,sample,sp_dir,sample)
        
            cmd_list = [cmd_1,rg_cmd,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7]

            final_cmd = "\n\n".join(cmd_list)


    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-0:00",mem="10000",cores="2",nodes="1",jobid="dedup",sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/04_dedup_sort_validate_%s.sbatch"%(sp_dir,sample)
            out_file = open(out_filename,"w")
            out_file.write(sra_script)
            out_file.close
            dedup_sbatch_filenames.append(out_filename)
    
    return(dedup_sbatch_filenames)

#Collect all deduplication metrics for a list of sample IDs, writes all metrics to a file and returns a dictionary of % duplication for each sample
def collect_dedup_metrics(stats_dir,sample_list):
    

#Collect all alignment summary metrics for a list of sample IDs, writes all metrics to a file and returns a dictionary of important metrics for each sample
def collect_alignment_metrics(stats_dir,sample_list):



def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    now = datetime.datetime.now()
    print('Staring work on script 02: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get Sample, SRA and Genome attributes - use same config as for pipeline script 01
#     '''
#     Example config file format (separate by spaces), only include one genome option: 
#     --ABBV <Species_Abbr>
#     --OUT_DIR <output directory>
#     --SAMPLE_NCBI <SAMPLE_1> <SRA_ID,SRA_ID,SRA_ID>
#     --SAMPLE_NCBI <SAMPLE_2> <SRA_ID,SRA_ID>
#     --GENOME_NCBI <NCBI Genome Accession>
#     --GENOME_LOCAL <Local genome fasta file>
#     '''
    
    config_info = extract_config(config_filename)

    #####Check if dedup directory exists, if not creates it
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    dedup_dir = "%s/dedup"%(sp_dir)
    directory_create(dedup_dir)


    #####Create sbatch files to dedup SRA files and combine if multiple SRAs for each sample, sort and index resulting file, validate, and compute coverage histogram
    
    #Create sbatch files
    dedup_filenames = dedup_sbatch(sp_dir,config_info["abbv"],config_info["sample_ncbi_dict"])
    '''
    #Submit dedup read sbatch files
    dedup_jobids = []
    completed_jobids = {}
    for i in range(0,len(dedup_filenames)):
        dedup_jobids.append(sbatch_submit(dedup_filenames[i]))
        sleep(1)
    #Add an extra sleep to give sacct a chance to catch up
    sleep(20)
    #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
    while len(completed_jobids) < len(dedup_filenames):
        num_running = num_pend_run(dedup_jobids,start_date)
        job_statuses = all_jobs_status(start_date)
        for job in dedup_jobids:
            if job not in completed_jobids:
                if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                    completed_jobids[job] = job_statuses[job]
                    print("Job %s completed"%job)
        sleep(30)
    
    #After all jobs have finished, report which jobs failed
    for job in completed_jobids:
        if completed_jobids[job] != "COMPLETED":
            print("Dedup, sort and validate job %s failed with code: %s"%(job,completed_jobids[job]))
    
    #####Collate summary statistics   
    #Collect alignment metrics and dedup metrics into single files for all samples
    
    
    
    #Collect coverage histograms, plot in pdf (species_date.pdf), mean and median coverage
    
    #Produce file of all most important summary stats. 
    
    #Copy this file and pdf of coverage to centralized location        
    
    #Check that the final sorted bam and index is available, if so, remove intermediate files (dedup)
    '''
    '''
    for sample in config_info["sample_ncbi_dict"]:
        for sra in config_info["sample_ncbi_dict"][sample]:
            if os.path.isfile('%s/%s.sorted.bam'%(alignment_dir,sra)) and os.path.isfile('%s/%s.sorted.bai'%(alignment_dir,sra)):
                    proc = Popen('rm %s/%s*'%(fastq_dir,sra),shell=True)
                    proc = Popen('rm %s/%s.sra'%(sra_dir,sra),shell=True)
            else:
                print("Something happened with SRA: %s for sample: %s"%(sra,sample))        
    '''
if __name__ == "__main__":
    main()