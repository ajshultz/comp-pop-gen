#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE

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
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -t {time}\n#SBATCH --mem {mem}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o logs/{jobid}_%j.out\n#SBATCH -e logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)
    
#Submit filename to slurm with sbatch, returns job id number
def sbatch_submit(filename):
    proc = Popen('sbatch %s'%filename,shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    return stdout.strip('Submitted batch job ')

#Create an sbatch file for a given set of SRAs and split into fastq files
def sra_download_sbatch(sp_dir,sample_ncbi_dict):
    slurm_script = script_create()
    
    #Paths to various required software
    path_to_ascp="/n/home13/ashultz/.aspera/connect/bin/ascp"
    path_to_ascp_openssh="/n/home13/ashultz/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    path_to_sratools = "/n/home13/ashultz/sw/progs/sratoolkit.2.8.2-1-centos_linux64/bin/"
    path_to_sra_dls = "/n/holylfs/LABS/informatics/ashultz/CompPopGen/raw_sra_files"
    
    for sample in sample_ncbi_dict.keys():
        for sra in sample_ncbi_dict[sample]:
    
            cmd_1 = 'module load fastqc'
            cmd_2 = r'%sprefetch --force all --max-size 500000000 -a "%s|%s" --ascp-options "-QT -l 10G" %s'%(path_to_sratools,path_to_ascp,path_to_ascp_openssh,sra)
            cmd_3 = r'%sfastq-dump --outdir %s/fastq --gzip --split-files %s/sra/%s.sra'%(path_to_sratools,sp_dir,path_to_sra_dls,sra)
            cmd_4 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra) 
            
            final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",cmd=final_cmd)

            out_file = open("%s/scripts/sra_download_parse_%s.sbatch"%(sp_dir,sra),"w")
            out_file.write(sra_script)
            out_file.close
            #print(sra_script)


#Use FTP to download species genome fasta file. Will obtain FTP location by first downloading the current genbank assembly summary report.
def get_ncbi_genome(sp_dir,species_name,sp_abbr):

    #Recreate genome directory
    genome_dir = "%s/genome"%(sp_dir)
    
    #Species are named with spaces instead of "_" in NCBI records, so convert underscores to spaces first.
    species_name_spaces = re.sub("_"," ",species_name)
    
    #Download current genbank assembly summary report
    #wget_ncbi_summary = 'wget -O %s/assembly_summary_genbank.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'%genome_dir
    wget_ncbi_summary = 'curl -O %s/assembly_summary_genbank.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'%genome_dir
    Popen(wget_ncbi_summary,shell=True,stdout=PIPE,stderr=PIPE)
    
    #Use grep to pull in relevant lines for this genome
    proc = Popen(r'grep "%s" %s/assembly_summary_genbank.txt'%(species_name_spaces,genome_dir),shell=True, stdout=PIPE, stderr=PIPE)
    stdoutstr,stderrstr = proc.communicate()
    genome_opts = stdoutstr.decode("utf-8","ignore")
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
            print("\nDownloading %s genome accession %s from: %s\n\nCopying fasta file to %s.fa and indexing with samtools faidx and bwa index"%(species_name,possible_genome_opts[0][0],full_genome_path,sp_abbr))

        elif len(possible_genome_opts) > 1:
            sys.exit("There seems to be more than one representative genome for the species you provided, which is problematic. See the genbank assebmly summary table in the genome directory to refine genome name before proceeding, or supply a local fasta.")
        else:
            sys.exit("There does not seem to be a 'representative genome' for your species. See the genbank assebmly summary table in the genome directory to refine genome name before proceeding, or supply a local fasta.") 
    
    slurm_script = script_create()
    
    cmd_1 = 'module load samtools\nmodule load bwa'
    cmd_2 = 'wget -P %s %s'%(genome_dir,full_genome_path)
    cmd_3 = 'gunzip %s/%s'%(genome_dir,genome_filename)
    cmd_4 = 'mv %s/%s %s/%s.fa'%(genome_dir,genome_filename[:-3],genome_dir,sp_abbr)
    cmd_5 = 'samtools faidx %s/%s.fa'%(genome_dir,sp_abbr)
    cmd_6 = 'bwa index %s/%s.fa'%(genome_dir,sp_abbr)
    
    final_cmd = "%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6)    
    
    #Format sbatch script
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_DL_Index",cmd=final_cmd)

    out_file = open("%s/scripts/genome_download_index_%s.sbatch"%(sp_dir,sp_abbr),"w")
    out_file.write(genome_script)
    out_file.close


#Make sure that genome and genome indexes are set up in genome directory for mapping
def process_local_genome(sp_dir,genome_local,sp_abbr,genome_present,bwa_index_present,faidx_index_present):
    genome_dir = "%s/genome"%(sp_dir)
    
    cmd_1 = 'module load samtools\nmodule load bwa'
    
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
        
    final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    #Format sbatch script and write file
    slurm_script = script_create()
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_CP_Index",cmd=final_cmd)

    out_file = open("%s/scripts/genome_cp_index_%s.sbatch"%(sp_dir,sp_abbr),"w")
    out_file.write(genome_script)
    out_file.close








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

    #####Check if species directory, logs, scripts, fastq, fastqc, genome, and alignment directories, if not creates them one
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("\nOutput will be written to %s"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    fastq_dir = "%s/fastq"%(sp_dir)
    fastqc_dir = "%s/fastqc"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    alignment_dir = "%s/alignment"%(sp_dir)
    
    directory_create(sp_dir)
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(fastq_dir)
    directory_create(fastqc_dir)
    directory_create(genome_dir)
    directory_create(alignment_dir)
    
    #####Download SRA files and use fastq-dump to split
    
    #Create sbatch files
    sra_download_sbatch(sp_dir,config_info["sample_ncbi_dict"])
    
    #Submit sbatch files - only allow up to 50 jobs to be running at one time.
    #use subprocess
    

    #####Prepare genome
    #Create sbatch file to download genome if not already present (checks for abbv.fa), mv to abbv.fa.
    #e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.2_ASM73873v2/GCF_000738735.2_ASM73873v2_genomic.fna.gz
    
    if "genome_ncbi" in config_info:
        #Create sbatch script
        get_ncbi_genome(sp_dir,config_info["genome_ncbi"],config_info["abbv"])
        #Submit sbatch script
        
    elif "genome_local" in config_info:
        #Check if genome, BWA indexes, and faidx indexes already exist in genome directory. If not, create script to copy and index.
        genome_path = "%s/%s.fa"%(genome_dir,config_info["abbv"])
        genome_present = os.path.isfile(genome_path)
        index_path_bwa = "%s/%s.fa.bwt"%(genome_dir,config_info["abbv"])
        bwa_index_present = os.path.isfile(index_path_bwa)
        index_path_faidx = "%s/%s.fa.fai"%(genome_dir,config_info["abbv"])
        faidx_index_present = os.path.isfile(index_path_faidx)

        #Create sbatch script if any missing elements (genome, faidx or bwa index)
        if genome_present == False or bwa_index_present == False or faidx_index_present == False:  
                process_local_genome(sp_dir,config_info["genome_local"],config_info["abbv"],genome_present,bwa_index_present,faidx_index_present)
        
            #Submit sbatch script



    #Once SRA dl jobs are finished, check that all fastq files are there. If so, continue below. If not, 
    
    #Trim fastq files with NGmerge
    #Map fastq files to genome with BWA
    #Set Read Group information
    
    

    #Calculate alignment stats

if __name__ == "__main__":
    main()