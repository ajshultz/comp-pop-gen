#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep

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
                sra_list = line[2].split(",")
                sample_ncbi_dict[line[1]] = sra_list
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
    return stdout.strip('Submitted batch job ')


#Check job status of specific jobid: returns job status
def jobid_status(jobid):
    proc = Popen('sacct --format state --noheader -j %d'%int(jobid),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split()
    return(lines[0].decode("utf-8","ignore"))
    
    
#Check the status of all jobs, returns dictionary of jobid:status. Ignores jobs with ".ba+" and ".ex+" file extensions.
def all_jobs_status():
    proc = Popen('sacct --format jobid,state --noheader',shell=True,stdout=PIPE,stderr=PIPE)
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
    

#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def sra_download_sbatch(sp_dir,sample_ncbi_dict):
    slurm_script = script_create()
    sra_dl_sbatch_filenames = []
    
    #Paths to various required software
    path_to_ascp="/n/home13/ashultz/.aspera/connect/bin/ascp"
    path_to_ascp_openssh="/n/home13/ashultz/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    path_to_sratools = "/n/home13/ashultz/sw/progs/sratoolkit.2.8.2-1-centos_linux64/bin/"
    path_to_sra_dls = "/n/holylfs/LABS/informatics/ashultz/CompPopGen/raw_sra_files"
    
    for sample in sample_ncbi_dict.keys():
        for sra in sample_ncbi_dict[sample]:
            print('Will download %s for sample %s, split, and run through fastqc'%(sra,sample))
    
            #Load modules and get versions for all programs used
            cmd_1 = 'module load fastqc/0.11.5-fasrc01'
            cmd_2 = r'%sprefetch --force all --max-size 500000000 -a "%s|%s" --ascp-options "-QT -l 10G" %s'%(path_to_sratools,path_to_ascp,path_to_ascp_openssh,sra)
            cmd_3 = r'%sfastq-dump --outdir %s/fastq --gzip --split-files %s/sra/%s.sra'%(path_to_sratools,sp_dir,path_to_sra_dls,sra)
            cmd_4 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra) 
            
            final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/01_sra_download_parse_%s.sbatch"%(sp_dir,sra)
            out_file = open(out_filename,"w")
            out_file.write(sra_script)
            out_file.close
            sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)
            


#Use FTP to download species genome fasta file. Will obtain FTP location by first downloading the current genbank assembly summary report. Returns name of sbatch script created.
def get_ncbi_genome(sp_dir,species_name,sp_abbr):

    #Recreate genome directory
    genome_dir = "%s/genome"%(sp_dir)
    
    #Species are named with spaces instead of "_" in NCBI records, so convert underscores to spaces first.
    species_name_spaces = re.sub("_"," ",species_name)
    
    #Download current genbank assembly summary report
    #wget_ncbi_summary = 'wget -O %s/assembly_summary_genbank.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'%genome_dir
    wget_ncbi_summary = 'curl ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt > %s/assembly_summary_genbank.txt'%genome_dir
    proc = Popen(wget_ncbi_summary,shell=True,stdout=PIPE,stderr=PIPE)
    proc.wait()
    
    #Use grep to pull in relevant lines for this genome
    grep_command = r'grep "%s" %s/assembly_summary_genbank.txt'%(species_name_spaces,genome_dir)
    proc = Popen(grep_command,shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    print(stdout)
    genome_opts = stdout.decode("utf-8","ignore")
    genome_opts = genome_opts.strip()
    genome_opts = genome_opts.split("\n")
    
    #If there are no entries print an error message.
    if len(genome_opts[0]) == 0:
        sys.exit("Are you sure the NCBI genome name is correct? It is not in the genbank assembly summary table (see file in genome directory).")
    
    #If there is one entry, take FTP from that. If there are more than 1, take "representative genome" for FTP. If more than one "representative genome", break and give an error
    if len(genome_opts) == 1:
        genome_opts = genome_opts[0].split("\t")
        genome_ftp_path = genome_opts[19] 
        genome_filename = '%s_genomic.fna.gz'%genome_ftp_path.split("/")[-1]
        full_genome_path = '%s/%s'%(genome_ftp_path,genome_filename)
        print("Downloading %s genome accession %s from: %s"%(species_name,genome_opts[0],full_genome_path))

    elif len(genome_opts) > 1:
        possible_genome_opts = []
        for i in range(0,len(genome_opts)):
            possible = genome_opts[i].split("\t")
            if possible[4] == "representative genome":
                possible_genome_opts.append(possible)
        
        #Double check single entry
        if len(possible_genome_opts) == 1:
            genome_ftp_path = possible_genome_opts[0][19]                
            genome_filename = '%s_genomic.fna.gz'%genome_ftp_path.split("/")[-1]
            full_genome_path = '%s/%s'%(genome_ftp_path,genome_filename)
            print("\nDownloading %s genome accession %s from: %s\n\nCopying fasta file to %s.fa and indexing with samtools faidx and bwa index, and creating a sequence dictionary with samtools dict\n"%(species_name,possible_genome_opts[0][0],full_genome_path,sp_abbr))

        elif len(possible_genome_opts) > 1:
            sys.exit("There seems to be more than one representative genome for the species you provided, which is problematic. See the genbank assebmly summary table in the genome directory to refine genome name before proceeding, or supply a local fasta.")
        else:
            sys.exit("There does not seem to be a 'representative genome' for your species. See the genbank assebmly summary table in the genome directory to refine genome name before proceeding, or supply a local fasta.") 
    
    slurm_script = script_create()
    
    #Load modules, also print samtools and bwa versions
    cmd_1 = 'module load samtools/1.5-fasrc01\nmodule load bwa/0.7.15-fasrc01'
    cmd_2 = 'wget -P %s %s'%(genome_dir,full_genome_path)
    cmd_3 = 'gunzip %s/%s'%(genome_dir,genome_filename)
    cmd_4 = 'mv %s/%s %s/%s.fa'%(genome_dir,genome_filename[:-3],genome_dir,sp_abbr)
    cmd_5 = 'samtools faidx %s/%s.fa'%(genome_dir,sp_abbr)
    cmd_6 = 'bwa index %s/%s.fa'%(genome_dir,sp_abbr)
    cmd_7 = 'samtools dict -o %s/%s.dict %s/%s.fa'%(genome_dir,sp_abbr,genome_dir,sp_abbr)
    
    final_cmd = "%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7)    
    
    #Format sbatch script
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_DL_Index",sp_dir=sp_dir,cmd=final_cmd)

    out_filename = "%s/scripts/02_genome_download_index_%s.sbatch"%(sp_dir,sp_abbr)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close
    
    return(out_filename)


#Make sure that genome and genome indexes are set up in genome directory for mapping. Returns name of sbatch script created.
def process_local_genome(sp_dir,genome_local,sp_abbr,genome_present,bwa_index_present,faidx_index_present,dict_index_present):
    genome_dir = "%s/genome"%(sp_dir)
    
    #Load modules, also print samtools and bwa versions
    cmd_1 = 'module load samtools/1.5-fasrc01\nmodule load bwa/0.7.15-fasrc01'
    
    #Copy genome if necessary
    if genome_present == False:
        print("\nCopying %s to %s/%s.fa"%(genome_local,genome_dir,sp_abbr))
        cmd_2 = 'cp %s %s/%s.fa'%(genome_local,genome_dir,sp_abbr)
    else:
        cmd_2 = ''
    
    #Index with samtools faidx if necessary
    if faidx_index_present == False:
        print("\nIndexing %s/%s.fa with Samtools faidx"%(genome_dir,sp_abbr))
        cmd_3 = 'samtools faidx %s/%s.fa'%(genome_dir,sp_abbr)
    else:
        cmd_3 = ''
        
    #Index with BWA if necessary
    if bwa_index_present == False:
        print("\nIndexing %s/%s.fa with bwa index"%(genome_dir,sp_abbr))
        cmd_4 = 'bwa index %s/%s.fa'%(genome_dir,sp_abbr)
    else:
        cmd_4 = ''
        
    if dict_index_present == False:
        print("\nCreating sequence dictionary file for %s/%s.fa"%(genome_dir,sp_abbr))
        cmd_5 = 'samtools dict -o %s/%s.dict %s/%s.fa'%(genome_dir,sp_abbr,genome_dir,sp_abbr)
    else:
        cmd_5 = ''
        
    final_cmd = "%s\n\n%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4,cmd_5)
    
    #Format sbatch script and write file
    slurm_script = script_create()
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_CP_Index",sp_dir=sp_dir,cmd=final_cmd)

    out_filename = "%s/scripts/02_genome_cp_index_%s.sbatch"%(sp_dir,sp_abbr)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close

    return(out_filename)


#Create sbatch files for trimming, mapping, sorting, indexing and alignment stat creation for each set of paired SRA fastq files.
def fastq_trim_align_stats(sp_dir,sra,sp_abbr,sample):
    print('Will trim, map, sort, index, and collect alignment stats for %s from sample %s'%(sra,sample))

    #Set read group info: ID = SRA number, SM = sample id, PU: SRA.sample, LB: sample id 
    read_group_info = '@RG\\tID:%s\\tSM:%s\\tPU:%s.%s\\tLB:%s'%(sra,sample,sra,sample,sample)
    path_to_picard = '~/sw/bin/picard.jar'

    #Load necessary modules
    cmd_1 = 'module load NGmerge/0.2-fasrc01\nmodule load bwa/0.7.15-fasrc01\nmodule load java/1.8.0_45-fasrc01'
    
    #If paired end data:
    #Trim adapters with NGmerge
    cmd_2 = 'NGmerge -1 %s/fastq/%s_1.fastq.gz -2 %s/fastq/%s_2.fastq.gz -a -v -o %s/fastq/%s_trimmed -n 8 -z -u 60'%(sp_dir,sra,sp_dir,sra,sp_dir,sra)
    
    #Map to genome with BWA mem
    cmd_3 = r"bwa mem -M -t 8 -R '%s' %s/genome/%s.fa %s/fastq/%s_trimmed_1.fastq.gz %s/fastq/%s_trimmed_2.fastq.gz > %s/alignment/%s_bwa.sam"%(read_group_info,sp_dir,sp_abbr,sp_dir,sra,sp_dir,sra,sp_dir,sra)
    
    cmd_4 = 'java -Xmx7g -XX:ParallelGCThreads=6 -jar %s SortSam I=%s/alignment/%s_bwa.sam O=%s/alignment/%s.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true'%(path_to_picard,sp_dir,sra,sp_dir,sra)
    
    #Calculate alignment summary stats
    cmd_5 = 'java -Xmx7g -XX:ParallelGCThreads=6 -jar %s CollectAlignmentSummaryMetrics I=%s/alignment/%s.sorted.bam R=%s/genome/%s.fa METRIC_ACCUMULATION_LEVEL=SAMPLE O=%s/stats/%s.%s.alignment_metrics.txt'%(path_to_picard,sp_dir,sra,sp_dir,sp_abbr,sp_dir,sample,sra)

    cmd_6 = 'if [ -f %s/stats/%s.%s.alignment_metrics.txt ]\nthen\n rm %s/alignment/%s_bwa.sam \nfi'%(sp_dir,sample,sra,sp_dir,sra)
    
    cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6]
    
    final_cmd = "\n\n".join(cmd_list)
    
    #Format sbatch script and write file
    slurm_script = script_create()
    genome_script = slurm_script.format(partition="shared",time="1-0:00",mem="16000",cores="8",nodes="1",jobid="Trim_Map",sp_dir=sp_dir,cmd=final_cmd)

    out_filename = "%s/scripts/03_trim_map_stats_%s.sbatch"%(sp_dir,sra)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close

    return(out_filename)





def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    #Open config file and get Sample, SRA and Genome attributes
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

    #####Check if species directory, logs, scripts, fastq, fastqc, genome, alignment, and stats directories, if not creates them
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("\nOutput will be written to %s\n"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    fastq_dir = "%s/fastq"%(sp_dir)
    fastqc_dir = "%s/fastqc"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    alignment_dir = "%s/alignment"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    
    directory_create(sp_dir)
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(fastq_dir)
    directory_create(fastqc_dir)
    directory_create(genome_dir)
    directory_create(alignment_dir)
    directory_create(stats_dir)
    
    
    #####Prepare genome
    #Create sbatch file to download genome if not already present (checks for abbv.fa), mv to abbv.fa.
    
    genome_path = "%s/%s.fa"%(genome_dir,config_info["abbv"])
    index_path_bwa = "%s/%s.fa.bwt"%(genome_dir,config_info["abbv"])
    index_path_faidx = "%s/%s.fa.fai"%(genome_dir,config_info["abbv"])
    index_path_dict = "%s/%s.dict"%(genome_dir,config_info["abbv"])
    genome_jobnum = []
    
    if "genome_ncbi" in config_info:
        #Create sbatch script
        genome_sbatch_name = get_ncbi_genome(sp_dir,config_info["genome_ncbi"],config_info["abbv"])
        
        #Submit sbatch script
        #genome_jobid = sbatch_submit(genome_sbatch_name)
        
    elif "genome_local" in config_info:
        #Check if genome, BWA indexes, and faidx indexes already exist in genome directory. If not, create script to copy and index.
        genome_present = os.path.isfile(genome_path)
        bwa_index_present = os.path.isfile(index_path_bwa)
        faidx_index_present = os.path.isfile(index_path_faidx)
        dict_index_present = os.path.isfile(index_path_dict)

        #Create sbatch script if any missing elements (genome, faidx or bwa index)
        if genome_present == False or bwa_index_present == False or faidx_index_present == False or dict_index_present == False:  
            genome_sbatch_name = process_local_genome(sp_dir,config_info["genome_local"],config_info["abbv"],genome_present,bwa_index_present,faidx_index_present,dict_index_present)
            
            #Submit sbatch script
            #genome_job_id = sbatch_submit(genome_sbatch_name)
        
    #Sleep 30 seconds to give job status time to get to sacct before starting to check    
    sleep(30)
    
    dones = ['COMPLETED','CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL']
    #Check job id status of genome job. If not in one of the 'done' job status categories, wait 30 seconds and check again.
    while jobid_status(genome_job_id) not in dones:
        print("Genome not done yet")
        sleep(30)
    
    #Check to make sure job completed, and that all necessary files are present. If not, exit and give information.
    genome_job_completion_status = jobid_status(genome_job_id)
    if genome_job_completion_status != 'COMPLETED':
        sys.exit("There was a problem downloading and indexing the genome. The job exited with status %s. Please diagnose and fix before moving on"%genome_job_completion_status)
        
    genome_present = os.path.isfile(genome_path)
    bwa_index_present = os.path.isfile(index_path_bwa)
    faidx_index_present = os.path.isfile(index_path_faidx)
    dict_index_present = os.path.isfile(index_path_dict)
    if genome_present == False or bwa_index_present == False or faidx_index_present == False or dict_index_present == False:
        sys.exit("The genome job finished but not all files and indexes are present. Please check, create missing indexes, and resubmit with local genome")

    print("\nGenome download and indexing successfully finished\n")

    #####Create sbatch files to download SRA files and use fastq-dump to split
    
    #Create sbatch files
    sra_dl_sbatch_filenames = sra_download_sbatch(sp_dir,config_info["sample_ncbi_dict"])
        
    #Submit SRA read sbatch files            
    sra_dl_jobids = []
    


    #Once SRA dl and genome jobs are finished, check that all fastq files are there. If so, continue below. If not, 
    
    #Trim fastq files with NGmerge
    #Set Read Group information
    #Map fastq files to genome with BWA
    #Sort (convert to BAM), Index
    #Calculate alignment stats
    #Remove SAM file if stats file is there
    mapping_jobids = []
    
    for sample in config_info["sample_ncbi_dict"]:
        for sra in config_info["sample_ncbi_dict"][sample]:
            sra_map_sbatchfile = fastq_trim_align_stats(sp_dir,sra,config_info["abbv"],sample)
            #sra_map_jobid = sbatch_submit(sra_map_sbatchfile)
            #mapping_jobids.append(sra_map_jobid)

if __name__ == "__main__":
    main()