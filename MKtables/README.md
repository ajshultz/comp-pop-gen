Information about running scripts to calculate MK tables.

1.) Create a directory for a particular ingroup species. All subsequent analyses
will be done within this directory.

2.) A configuration file must be created for each species that contains the location 
of directories and files that are used by the various scripts. An example for T.
gutatta may be found in this directory entitled "ExampleConfigFile.txt".

You may also find a file named "Example.sh" that is an example of the how you
may submit the next parts of the pipeline as a job, which consist of the name of a
script, followed by the name of the configuration file that it takes as input on 
the command line.

3.) Make gff database file using '00_make_gffDatabase.py'. This database  is used by
the 'gffutils' python module. This can take a little bit of time, but no more than
an hour. While this is running I tried to set up directories for additional species.

4.) Run the 'Combined_wrapper.py' script. While this may be submitted as a job, I 
ran this on bioinf01 so that there was no time restriction. However, if jobs are
not pending long and running quickly, this is not necessary.

This script will create 3 directories: "sbatchScripts" which contains the jobs
submitted by 'Combined_wrapper.py', and two directories named "*_AlleleTables" and
"*_MKtables", which contain the allele tables and MK tables per scaffold. These
directories may contain many files, so use 'ls' with caution.

5.) Run the '03_collect_results.py' script to combine MK table results across
scaffolds into a single file named "PolymorphismDivergenceStats_combined.txt".


