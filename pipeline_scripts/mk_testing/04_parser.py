# to set up conda env for cyvcf2 on cluster

module load Anaconda3/5.0.1-fasrc02
conda create -n cyvcf2 python=3.6 anaconda
source activate cyvcf2
conda install -c bioconda -n cyvcf2 cyvcf2




#!/usr/bin/python -tt

import re
import sys
import os
from cyvcf2 import VCF
# see https://brentp.github.io/cyvcf2/docstrings.html

def main():

	vcf = VCF('/scratch/swuitchik/CompPopGen/consensus_ingroups/corCor/corCor.ann.vcf', gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles

	outFile = open("corCor.effects.txt", 'w')
  
	for variant in vcf:
		# check if this site passed GATK VCF filtering, and only look at bi-allelic sites
		if len(variant.ALT) == 1:
			print(variant.CHROM, variant.end, sep="\t", end="\t", file=outFile)
			for field in variant.INFO:
				if field[0] == 'ANN':
					annotations = field[1].split(",")
					effects = []
					for i in annotations:
						i = i.split('|')
						effects.append(i[1])
					print(",".join(effects), file=outFile)

	outFile.close()
if __name__ == '__main__':
  main()


# to deactivate cyvcf2 environment, use: source deactivate

## Written with Brian J Arnold (barnold@harvard.edu)

