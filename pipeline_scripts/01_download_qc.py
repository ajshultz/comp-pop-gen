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
        
    try:
        config_info["genome_ncbi"] or config_info["genome_local"]
    except UnboundLocalError:
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
            cmd_4 = 'fastqc -o %s/fastqc %s_1.fastq.gz %s_2.fastq.gz'%(sp_dir,sra,sra)
            
            
            final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",cmd=final_cmd)

            out_file = open("%s/scripts/sra_download_parse_%s.sbatch"%(sp_dir,sra),"w")
            out_file.write(sra_script)
            out_file.close
            #print(sra_script)



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

    #Check if species directory, logs, scripts, fastq, fastqc, genome, and alignment directories, if not creates them one
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("Output will be written to %s"%sp_dir)
    
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
    
    #Download SRA files and use fastq-dump to split
    sra_download_sbatch(sp_dir,config_info["sample_ncbi_dict"])
    


    #Download genome if not already present
    #e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.2_ASM73873v2/GCF_000738735.2_ASM73873v2_genomic.fna.gz

    #Run fastqc on fastq files

    #Index genome if not already present
    #Submit job once genome is downloaded

    #Trim fastq files with NGmerge

    #Map fastq files to genome with BWA
    #Set Read Group information

    #Calculate alignment stats

if __name__ == "__main__":
    main()