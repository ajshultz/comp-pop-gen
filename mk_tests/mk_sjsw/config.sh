#!/usr/bin/bash

# modules required
module load bcftools/1.5-fasrc02 vcftools/0.1.14-fasrc01 htslib/1.5-fasrc02 jdk/10.0.1-fasrc01 bedtools2/2.26.0-fasrc01 python/3.6.3-fasrc02 Anaconda3/5.0.1-fasrc02 
export R_LIBS_USER=$HOME/path/to/R/packages

# for example: 
# export R_LIBS_USER=$HOME/apps/R_3.6.1:$R_LIBS_USER

# create an environment that has all the other packages you'll need 
conda create -n mk python=3.6 anaconda cyvcf2 tqdm r-base r-tidyverse r-rjags r-r2jags r-lme4 r-arm
source activate mk

# set up project directory as the structure outlined in directory_tree.pdf 

export INSHORT=ingroup_spp_name (six letter code)
export OUTSHORT=outgroup_spp_name (six letter code)
export INLONG=ingroup_spp_name (longform spp name, with leading underscore)
export OUTLONG=outgroup_spp_name (longform spp name, with leading underscore)

# for example: 

# export INSHORT=corCor
# export OUTSHORT=corMon
# export INLONG=_Ccornix
# export OUTLONG=_Cmonedula

export PATHW=$HOME/path/to/working/directory

wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip 
rm -r clinEff/
cd snpEff/
mkdir -p data/$INSHORT.ncbi/
# ensure reference sequence (FASTA) and genome annotation (GFF3) are in the appropriate data directory
# rename reference sequence to sequences.fa and gzip
# rename genome annotation to genes.gff and gzip

# add the following into the snpEff.config file under the Databases & Genomes section: 

# Common name genome, NCBI version __
$INSHORT.ncbi.genome : ncbi_genome_name

# for example: 

# # Hooded crow genome, NCBI version 2
# corCor.ncbi.genome : Corvus_cornix_cornix


export PATHS=$HOME/path/to/snpEff

# build database (from working directory)
java -jar snpEff/snpEff.jar build -gff3 -v $INSHORT.ncbi 

# in working directory, will need: 
# ingroup missingness
# outgroup missingness
# ingroup coverage sites
# outgroup coverage sites
# genes.gff
# genenames.py
# gff2bed.awk
# parser_nov.py
# pi2bed.awk
# my.jag2.R
# SnIPRE_source.R
# ingroup.vcfs directory
# outgroup.vcfs directory

