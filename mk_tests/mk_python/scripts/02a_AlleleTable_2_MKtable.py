#!/usr/bin/python -tt

""" 

"""

import sys
import re
import os
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sharedFunctions as shFn

def areExonsAtLeast3bp(db):
	cdsTooSmall = []
	for gene in db.features_of_type('gene'):
		for cds in db.children(gene, featuretype='CDS', order_by='start'):
			if cds.end - cds.start + 1 < 2:
				cdsTooSmall.append(cds.id)
	"""
	if len(cdsTooSmall)>0:
		print("WARNING, FOLLOWING CDS ARE TOO SMALL FOR YOUR ALGORITHM:")
		for i in cdsTooSmall:
			print(i)
	"""


def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 5:
		# have comma separated lists for each field, split on comma
		gff3DB = sys.argv[1]
		refGenomeFile = sys.argv[2]
		alleleTableFilename = sys.argv[3]
		chromToGenotype = sys.argv[4]
	else:
		print("Didn't specify right number of input arguments, exiting...")
		sys.exit()
	
	db = gffutils.FeatureDB(gff3DB)					# this database was created in previous step
	# Assumption in analyses below is that exons are at least 1 a.a. residue long, 3 bp. confirm this
	areExonsAtLeast3bp(db)

	# create data structs to convert position -> functional effect
	Pos_2_CDS, CDSinfo = shFn.constructCdsDicts(db, chromToGenotype)  	# Pos_2_CDS[pos] = [ sorted CDS ids ]
	# Pos_2_CDS[pos] = [ sorted CDS ids ]
	# CDSinfo[firstCDS] = (0:geneName, 1:strand, 2:frame, 3:startPos, 4:endPos, 5:RefSeqScaffold, 6:UniqueGeneName) 

	refGenomeSeqDict = SeqIO.to_dict(SeqIO.parse(refGenomeFile, "fasta")) # refGenomeSeqDict[ scaffoldName ].seq = ATG...GGC, sequence of that scaffold
	translateCodonDict = shFn.constructTranslateCodonDict()

	mkDict = {} # mkDict[gene] = [SynPoly, NonsynPoly, SynFixed, NonsynFixed]
	for gene in db.features_of_type('gene'):
		if gene.seqid == chromToGenotype:
			mkDict[gene['Name'][0]] = [0,0,0,0]
	
	if not os.path.isfile(alleleTableFilename):
		print("Alleletable -> MKtable script cannot locate the allele table! exiting...")
		sys.exit()
	f = open(alleleTableFilename, 'r')
	# File structure:
	# Col 1: position
	# Col 2: Refbase
	# Col 3: Species1 allele
	# Col 4: Species1 allele frequency
	# ...
	# Col n-1: SpeciesX allele
	# Col n: SpeciesX allele frequency
	total=0
	biallelic=0
	biCoding=0
	numOutgroupPolymorphic = 0
	numOutgroupFixed = 0
	for line in f:
		if not re.match(r'^Position', line):
			mutInfo = line.split()
			posInScaff = int(mutInfo[0])
			refBase = str(mutInfo[1])
			#########
			# NOTE: only analysing ingroup and a single outgroup, for now
			#########
			ingroupBase = str(mutInfo[2])
			ingroupAF = float(mutInfo[3])
			outgroupBase = str(mutInfo[4])
			outgroupAF = float(mutInfo[5])

			allelesAllDict = {}		# used to check how many alleles at this site, across all species
			allelesAllDict[refBase] = 1	# also use ref base to see if site biallelic
			allelesAllDict[ingroupBase] = 1	
			allelesAllDict[outgroupBase] = 1	
			
			#allelesOutgroupDict = {}	# used to check how many alleles at this site, ONLY in outgroup species
			#i = 2				# ALT alleles start at index 2 (Col 3), and proceed until end, skipping a col for frequency each time
			#while i < len(mutInfo)-1:
			#	allelesAllDict[mutInfo[i]] = 1
			#	allelesOutgroupDict[mutInfo[i]] = 1
			#	i+=2
			#if len(allelesAllDict.keys()) == 2:
			#	biallelic+=1
			#if len(allelesOutgroupDict.keys()) == 1:		# IMPORTANT FILTERING STEP: OUTGROUP SPECIES MUST AGREE ON ANCESTRAL BASE, CHANGE ME
			#outgroupBase = next(iter(allelesOutgroupDict.keys()))	

			if posInScaff in Pos_2_CDS.keys():		# is position in coding region
				biCoding+=1
				inPolymorphic = None
				inFixed = None
				outPolymorphic = None
				outFixed = None
				# is ingroup polymorphic or fixed?
				if ingroupAF > 0.0 and ingroupAF < 1.0:
					inPolymorphic = 1
					inFixed = 0
				else:
					inPolymorphic = 0
					inFixed = 1
				# is outgorup polymorphic or fixed?
				if outgroupAF > 0.0 and outgroupAF < 1.0:
					outPolymorphic = 1
					outFixed = 0
					numOutgroupPolymorphic += 1
				else:
					outPolymorphic = 0
					outFixed = 1
					numOutgroupFixed += 1
				firstCDS = Pos_2_CDS[posInScaff][0]	# first of sorted mRNAs position corresponds to
				# Reminder: CDSinfo[firstCDS] = (0:geneName, 1:strand, 2:frame, 3:startPos, 4:endPos, 5:RefSeqScaffold, 6:UniqueGeneName)
				geneName = CDSinfo[firstCDS][6]
				strand = CDSinfo[firstCDS][1]
				refseqScaffold = CDSinfo[firstCDS][5]
				if strand == "+":
					p1,p2,p3 = shFn.getCodonPositionsInReference_forwardStrand(posInScaff, firstCDS, CDSinfo)
				elif strand == "-":
					p1,p2,p3 = shFn.getCodonPositionsInReference_reverseStrand(posInScaff, firstCDS, CDSinfo)

				# ignore sites in which flanking bases could not be found, i.e. flanking exon not located from partial gene model
				if p1 == None:
					continue

				posInCodonToChange = shFn.getPosInCodon(p1,p2,p3,posInScaff)

				# scaffolds in gff3 RefSeq IDs, scaffolds in reference genome Genbank IDs
				b1 = str(refGenomeSeqDict[refseqScaffold].seq[(p1-1):(p1)]).upper()
				b2 = str(refGenomeSeqDict[refseqScaffold].seq[(p2-1):(p2)]).upper()
				b3 = str(refGenomeSeqDict[refseqScaffold].seq[(p3-1):(p3)]).upper()
				codonRef = [b1, b2,  b3]
				if "N" not in codonRef:
					codonIngroup = codonRef.copy() # need to use .copy method, otherwise it creates new pointer to same object
					codonOutgroup = codonRef.copy()
					codonIngroup[posInCodonToChange] = ingroupBase.upper()
					codonOutgroup[posInCodonToChange] = outgroupBase.upper()
					
					codonRefStr = ''.join(codonRef)
					codonIngroupStr = ''.join(codonIngroup)
					codonOutgroupStr = ''.join(codonOutgroup)

					if strand == "-":
						ref = Seq(codonRefStr, generic_dna)
						ref = ref.reverse_complement()
						inG = Seq(codonIngroupStr, generic_dna)
						inG = inG.reverse_complement()
						outG = Seq(codonOutgroupStr, generic_dna)
						outG = outG.reverse_complement()
						codonRefStr = str(ref)
						codonIngroupStr = str(inG)
						codonOutgroupStr = str(outG)

					if inFixed:
						if outFixed and codonIngroupStr != codonOutgroupStr:
							#ingroup and outGroup fixed, can only be 2 bases
							if translateCodonDict[codonIngroupStr] == translateCodonDict[codonOutgroupStr]:
								mkDict[geneName][2] += 1	#synonymous
							else:
								mkDict[geneName][3] += 1	#nonsynonymous
						elif outPolymorphic:
							continue
							# WE ARE IGNORING THIS FOR NOW; 
							# this could represent a fixed diff in ingroup that's still polymorphic in outgroup
							# OR, it could be a a polymorphic in the outgroup where NOTHING happened in ingroup (unmutated, monomorphic)
						else:
							continue
							#inPolymorphic
							#outPolymorphic
							#codonIngroupStr == codonOutgroupStr, b/c the second outgroup has a different base
					
					if inPolymorphic:
						#is site biallelic?
						if len(allelesAllDict.keys()) == 2:
							if translateCodonDict[codonIngroupStr] == translateCodonDict[codonOutgroupStr]:
								mkDict[geneName][0] += 1	#synonymous
							else:
								mkDict[geneName][1] += 1	#nonsynonymous
						# if multiallelic, then Ingroup alternate allele != outgroup allele
						# these sites count as 1 polymorphism, 1 divergence
						elif len(allelesAllDict.keys()) > 2:
							# add 1 to polymorphism	
							if translateCodonDict[codonIngroupStr] == translateCodonDict[codonRefStr]:
								mkDict[geneName][0] += 1	#synonymous
							else:
								mkDict[geneName][1] += 1	#nonsynonymous
							# add 1 to divergence	
							if translateCodonDict[codonIngroupStr] == translateCodonDict[codonOutgroupStr]:
								mkDict[geneName][2] += 0.5	#synonymous
							else:
								mkDict[geneName][3] += 0.5	#nonsynonymous
							if translateCodonDict[codonRefStr] == translateCodonDict[codonOutgroupStr]:
								mkDict[geneName][2] += 0.5	#synonymous
							else:
								mkDict[geneName][3] += 0.5	#nonsynonymous

	mkFile = open('MKtable_%s.txt' % chromToGenotype, 'w')
	print("GeneName\tSynPolymorphic\tNonsynPolymorphic\tSynFixed\tNonsynFixed", file=mkFile)
	for gene in mkDict.keys():
		print(gene, "\t", end="", file=mkFile)
		for i in (0,1,2,3):	
			print(mkDict[gene][i], end = "\t", file=mkFile)
		print(file=mkFile)

	outFile = open('out_%s.txt' % chromToGenotype, 'w')
	print("Coding biallelic polymorphic sites: :", biCoding, file=outFile)
	print("Number of sites with polymorphic outgroup: ", numOutgroupPolymorphic, file=outFile)
	print("Number of sites with fixed outgroup: ", numOutgroupFixed, file=outFile)
	sys.exit()


if __name__ == '__main__':
  main()
