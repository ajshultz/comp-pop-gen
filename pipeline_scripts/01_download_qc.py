#!/usr/bin/env python

#For use with Python 3.

#Load standard modules
import re, sys, os, sets, getopt

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values.
def extract_config(config_filename):
    print("Opening %s"%config_filename)
    config_file = open(config_filename,"r")
    sample_dict = {}
    
    for line in config_file:
        if line[0] == "#"
            pass
        line = line.strip.split(" ")
        if line[0] == "--ABBV":
            abbv = line[1]
        elif line[0] == "--SAMPLE":
            sra_list = line[2].split(",")
            sample_dict[line[1]] = sra_list
        elif line[0] == "--GENOME_NCBI":
            genome_ncbi = line[1]
        elif line[0] == "--GENOME_LOCAL":
            genome_local = line[1]
    config_file.close()
    
    if genome_ncbi:
        return(abbv,sample_dict,genome_ncbi)
    elif genome_local:
        return(abbv,sample_dict,genome_local)


def main():

#Open config file and get Sample, SRA and Genome attributes

'''
Example config file format (separate by spaces), only include one genome option: 
--ABBV <Species_Abbr>
--SAMPLE <SAMPLE_1> <SRA_ID,SRA_ID,SRA_ID>
--SAMPLE <SAMPLE_2> <SRA_ID,SRA_ID>
--GENOME_NCBI <NCBI Genome Accession>
--GENOME_LOCAL <Local genome fasta file>
'''


#Download SRA files


#!/bin/bash

list=`cat accessions.lst`
njobs=24
path_to_ascp="/n/home13/sayshwarya/.aspera/connect/bin/ascp"
path_to_ascp_openssh="/n/home13/sayshwarya/.aspera/connect/etc/asperaweb_id_dsa.openssh"

for i in $list ; do
 echo -n "--> STARTING ASSESSION $i DOWNLOAD @ " ; date
 /n/holylfs/EXTERNAL_REPOS/SOFTWARE/sratoolkit.2.5.5-centos_linux64/bin/prefetch --force yes --max-size 500000000 -a "${path_to_ascp}|${path_to_ascp_openssh}" --ascp-options "-QT -l 10G" $i ./ > $i.log & 
 running=`jobs -r | grep -c Running`
 while [ "$running" -gt "$njobs" ] ; do
        sleep 1
        running=`jobs -r | grep -c Running`
        done
done



#Download genome if not already present
e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.2_ASM73873v2/GCF_000738735.2_ASM73873v2_genomic.fna.gz

#Split SRA files into fastqs

#Run fastqc on fastq files

#Index genome if not already present
#Submit job once genome is downloaded

#Trim fastq files with NGmerge

#Map fastq files to genome with BWA
#Set Read Group information

#Calculate alignment stats

if __name__ == "__main__":
    main()