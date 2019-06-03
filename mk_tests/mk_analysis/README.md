Analyzing MK Tables
============

Authors:
* Allison Shultz (Assistant Curator of Ornithology, Natural History Museum of Los Angeles County; ashultz@nhm.org)
* Brian Arnold (Senior Bioinformatics Scientist, Informatics Group, Harvard University; barnold@g.harvard.edu)
* Sara Wuitchik (Postdoc, Harvard University & Boston University; sjswuit@bu.edu)
* Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)


This directory contains code to analyze MK tables produced by mk_pipeline.

**Note that these results are a work in progress, and should be considered 

SnIPRE_TestSpecies_2019_Arguments.R Runs take the MK table and generates SnIPRE output. Three arguments are required 1) Species abbreviation 2) the number of individuals sequenced in the ingroup and 3) the number of individuals sequenced in the outgroup. Note it currently assumes the MK table is in a file called SpAbbr_PolymorphismDivergenceStats_combined.txt.

SnIPRE_code_JAGS contains all source code required to run SnIPRE.

mk_res_analyses.R contains R code to analyze SnIPRE output and produce some simple summaries across species.