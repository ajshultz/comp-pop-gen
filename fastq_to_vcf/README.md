FASTQ to VCF pipeline for reanalysis of SRA datasets
============

Authors:
* Allison Shultz (Assistant Curator of Ornithology, Natural History Museum of Los Angeles County; ashultz@nhm.org)
* Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)


Our objective is to create a semi-automated, reproducible, robust pipeline to take raw sequencing data (user supplied or from NCBI) and generate VCF files for downstream analyses with some quality filtering. 

The pipeline will have four steps:

1. Obtain data (if necessary), perform basic quality control, map to a reference genome using BWA (including downloading and indexing a reference genome if necessary), and calculate alignment statistics with PicardTools.

2. Perform additional data preprocessing (e.g. sorting, indexing, deduplication).

3. Use HaplotypeCaller and Genotype GVCF to call variants

4. Filter variants (using hard filtering) and produce additional qc stats (e.g. depth of coverage, PCA)

At this time, our scripts are built to work with the slurm architecture for job submission to a computing cluster. However, we hope to generalize this in the future.