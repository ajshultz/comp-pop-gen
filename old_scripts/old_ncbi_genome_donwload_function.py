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
            if possible[4] == "representative genome" or possible[4] == "reference genome":
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

    out_filename = "%s/scripts/01_genome_download_index_%s.sbatch"%(sp_dir,sp_abbr)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close
    
    return(out_filename)
