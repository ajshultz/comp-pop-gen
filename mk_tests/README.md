Generating MK tables from VCF files
============

Authors:
* Allison Shultz (Assistant Curator of Ornithology, Natural History Museum of Los Angeles County; ashultz@nhm.org)
* Brian Arnold (Senior Bioinformatics Scientist, Informatics Group, Harvard University; barnold@g.harvard.edu)
* Sara Wuitchik (Postdoc, Harvard University & Boston University; sjswuit@bu.edu)
* Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)


This directory contains code to create MK tables from VCF files.

mk_python is highly experimental and not-guaranteed-to-work Python code for this task.

mk_pipeline uses existing tools for most of the work but makes some simplifying assumptions as outlined in the directory.

mk_analysis contains the code to run SnIPRE once the MK table has been generated, and summarize results across spcies

mk_sjsw contains final draft of vcf2mk pipeline
