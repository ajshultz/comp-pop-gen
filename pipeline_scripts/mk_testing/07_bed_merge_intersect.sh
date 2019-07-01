## Merge and intersect callable regions

bedtools merge -i corCor.onlyCDS.gff.bed -c 10 -o distinct > corCor.onlyCDS.merge.gff.bed
bedtools intersect -a corCor.onlyCDS.merge.gff.bed -b vcfout.bed > corCor.consensus.bed

