#!/usr/bin/python -tt

""" 
This script calculates an allele table for each chromosome only for sites that are callable across all species included.
It also calculates the number of callable sites for each gene, as this information may be useful for downstream analyses.

WARNING: loading all sites in BED file into memory takes a LOT of memory, thus we analyse by chromosome
"""

import sys
import gzip
import re
import os
import gffutils
import array as arr
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.Alphabet import generic_dna
import sharedFunctions as shFn

def getChromList(bedList):
	ChromList = []
	# Only go through first file, as I assume all BED files will have all chromosomes represented
	bedFH = open(bedList[0], 'r')
	for line in bedFH:
		l = line.split()
		chrom = l[0]
		if chrom not in ChromList:
			ChromList.append(chrom)
	return(ChromList)

def getBedIntervalsList(bedList, chromToGenotype):
	# Open bed files
	# NOTE: bedfiles have 0-based coordinates, half open, meaning it doesn't include ending coordinate
	# convert to 1-based closed interval by adding 1 to start position and taking rest of coorinates as is
	
	bedIntervalsListOfArrays = [] # using list of array to save on memory, as opposed to 2D list
	for index in range(0,len(bedList)):
		bedIntervalsListOfArrays.append(arr.array('i'))
	index = 0
	for bedFN in bedList:
		bedFH = open(bedFN, 'r')
		for line in bedFH:
			l = line.split()
			chrom = l[0]
			if chrom == chromToGenotype:
				start = int( l[1] ) + 1 
				end = int( l[2] )
				bedIntervalsListOfArrays[index].extend(range(start, end+1))		
		bedFH.close()
		index += 1
	return(bedIntervalsListOfArrays)

def getCallableSites(bedList, bedIntList):
	CallableSitesDict = {}
	TempSet = set()
	for index in range(0,len(bedList)):
		if index==0:
			TempSet = set(bedIntList[index])
		else:
			TempSet = TempSet.intersection(set(bedIntList[index]))
	for x in TempSet:
		CallableSitesDict[x]=1	

	return(CallableSitesDict)


def printCoverageStatistics(CallableSitesDict, db, CDSinfo, refGenomeSeqDict, translateCodonDict, chromToGenotype):
	CallableSitesPerGene = {}
	CallableSynSitesPerGene = {}
	CallableNsynSitesPerGene = {}
	SitesPerGene = {}
	for gene in db.features_of_type('gene'):
		if gene.seqid == chromToGenotype:
			geneName = gene['Name'][0]
			CallableSitesPerGene[geneName] = 0
			CallableSynSitesPerGene[geneName] = 0
			CallableNsynSitesPerGene[geneName] = 0
			SitesPerGene[geneName] = 0
			largestmRNAsize = 0
			for mRNA in db.children(gene, featuretype='mRNA', order_by='start'):
				mRNAsize = mRNA.end - mRNA.start
				if mRNAsize >= largestmRNAsize:
					largestmRNA = mRNA.id
					largestmRNAsize = mRNAsize
			if largestmRNAsize > 0:
				for cds in db.children(largestmRNA, featuretype='CDS', order_by='id'):
					SitesPerGene[geneName] += (cds.end - cds.start + 1)
					# remove "cds" from name to get numerical ID, need for CDSinfo dict
					cdsIdNew = cds.id.replace('cds','')     # cdsIdNew = ##_##: left of _ represent mRNA it came from (may not 
										# be same # as rnd.id, but cds with same # came from same rna),
										# right of _ represent exon number
					geneName = CDSinfo[cdsIdNew][6]
					strand = CDSinfo[cdsIdNew][1]
					refseqScaffold = CDSinfo[cdsIdNew][5]
					for pos in range(cds.start, cds.end+1):
						if pos in CallableSitesDict:
							# GET FLANKING POSITIONS FOR CODON
							if strand == "+":
								p1,p2,p3 = shFn.getCodonPositionsInReference_forwardStrand(pos, cdsIdNew, CDSinfo)
							elif strand == "-":
								p1,p2,p3 = shFn.getCodonPositionsInReference_reverseStrand(pos, cdsIdNew, CDSinfo)
							# ignore sites in which flanking bases could not be found, i.e. flanking exon not located from partial gene model
							if p1 == None:
								continue
							# get index within codon to change to determine functional effect
							posInCodonToChange = shFn.getPosInCodon(p1,p2,p3,pos) 

							# GET BASES FROM REFERENCE
							b1 = str(refGenomeSeqDict[refseqScaffold].seq[(p1-1):(p1)]).upper()
							b2 = str(refGenomeSeqDict[refseqScaffold].seq[(p2-1):(p2)]).upper()
							b3 = str(refGenomeSeqDict[refseqScaffold].seq[(p3-1):(p3)]).upper()
							codonRef = [b1, b2,  b3]
							if "N" not in codonRef:
								# ignore all sites that cannot be classified as SYN/NSYN
								CallableSitesPerGene[geneName] += 1
								codons = []
								for base in ["A", "G", "C", "T"]:
									codonTemp = codonRef.copy()
									codonTemp[posInCodonToChange] = base
									codonTempStr = ''.join(codonTemp)
									codons.append(codonTempStr)
								
								if strand == "-":
									for i in range(0,len(codons)):
										codonTemp = Seq(codons[i], generic_dna)
										codonTemp = codonTemp.reverse_complement()
										codonTempStr = str(codonTemp)
										codons[i] = codonTempStr

								# 12 possible comparisons, use 12 as denominator (4 codons, 3 possible ocmparisons each)
								for i in range(0,len(codons)):
									for j in range(0,len(codons)):
										if i != j:
											if translateCodonDict[codons[i]] == translateCodonDict[codons[j]]:
												CallableSynSitesPerGene[geneName] += (1/12)
											else:
												CallableNsynSitesPerGene[geneName] += (1/12)
	outFile = open("coverageStats_%s.txt" % chromToGenotype, 'w')
	print("GeneName\tCallableSites\tCallableSynSites\tCallableNsynSites\tTotalSites", file=outFile)
	for gene in CallableSitesPerGene.keys():
		print(gene, " ", CallableSitesPerGene[gene], " ", CallableSynSitesPerGene[gene], " ", CallableNsynSitesPerGene[gene], " " , SitesPerGene[gene], file=outFile)

def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 6:
		# have comma separated lists for each field, split on comma
		vcfTempList = sys.argv[1]
		bedTempList = sys.argv[2]
		gff3DB = sys.argv[3]
		refGenomeFile = sys.argv[4]
		chromToGenotype = sys.argv[5]

		vcfList = vcfTempList.split(",")
		bedList = bedTempList.split(",")

		#Check to make sure files exist!
		if not os.path.isfile(gff3DB):
			print("can't locate gff3 file")
			sys.exit()
		for vcf in vcfList:
			if not os.path.isfile(vcf):
				print("can't locate one of your VCF files")
				sys.exit()
		for bed in bedList:
			if not os.path.isfile(bed):
				print("can't locate one of your BED files")
				sys.exit()
	else:
		print("USAGE: Specify two things, comma-separated (no spaces) list of vcfs, comma-separated (no spaces) list of bed files")
		sys.exit()
	
	# access information in gff3 file
	db = gffutils.FeatureDB(gff3DB)		# this database was created in previous step
	# get list of chromosomes
	ChromList = getChromList(bedList) 	# list of chromosome names

	# Open bed files
	# NOTE: bedfiles have 0-based coordinates, half open, meaning it doesn't include ending coordinate
	# convert to 1-based closed interval by adding 1 to start position and taking rest of coorinates as is
	bedIntervalsListOfArrays = getBedIntervalsList(bedList, chromToGenotype) # 2D Dict of sets, bedIntervalsDict[bedFN][chrom] = set(range(start, end+1))

	# for each chrom, find intersection of all intervals across all input bed files
	# load into 2D dict CallableSitesDict, CallableSitesDict[chrom][pos]=1
	CallableSitesDict = getCallableSites(bedList, bedIntervalsListOfArrays) 
	del bedIntervalsListOfArrays[:] # reduce memory footprint

	# create data structs to convert position -> functional effect, also store information for CDSs
	Pos_2_CDS, CDSinfo = shFn.constructCdsDicts(db, chromToGenotype)	
	# Pos_2_CDS[pos] = [ sorted CDS ids ]
	# CDSinfo[firstCDS] = (0:geneName, 1:strand, 2:frame, 3:startPos, 4:endPos, 5:RefSeqScaffold, 6:UniqueGeneName) 

	refGenomeSeqDict = SeqIO.to_dict(SeqIO.parse(refGenomeFile, "fasta")) # refGenomeSeqDict[ scaffoldName ].seq = ATG...GGC, sequence of that scaffold
	translateCodonDict = shFn.constructTranslateCodonDict()
	
	# Calculate and print number of callable sites per gene, using the largest mRNA for each gene
	# Also print number of synonymous and nonsynonymous sites per gene
	printCoverageStatistics(CallableSitesDict, db, CDSinfo, refGenomeSeqDict, translateCodonDict, chromToGenotype)

	RefDict = {}	# 2D dict, RefDict[chromosome][site] = reference allele
	AltDict = {}	# 3D dict, AltDict[vcfFile][chromosome][site] = (alternate allele, frequency)
	MultiDict = {} 	# Dict recording multiallelic sites and indels
	# Initialize Dicts
	for vcfFN in vcfList:
		AltDict[vcfFN] = {}

	chr_regex = r"^" + re.escape(chromToGenotype)
	# Go through each VCF file
	for vcfFN in vcfList:	
		myVcf = gzip.open(vcfFN,'rt')
		#myVcf = vcf.Reader(filename=vcfFN)
		# Go through ea line in VCF
		#for record in myVcf:
		for record in myVcf:
			if re.match(chr_regex, record): 
				lineList = re.split(r'\t+', record)
				# col 7 is FILTER, col 4 is REF, col 5 is ALT
				# This ensures PASSing site, and ref and alt only single base
				if lineList[6] == "PASS" and re.match(r'^[AGCT]{1}$', lineList[3]) and re.match(r'^[AGCT]{1}$', lineList[4]):	
					pos=int(lineList[1]) 	# 2nd col of vcf is POSITION
					if pos in CallableSitesDict:	
						ref=lineList[3]			# 4th col of vcf is REFERENCE base
						RefDict[pos] = ref 		# position added if at least 1 sample passed, occasionally gets overwritten
						alt = lineList[4] 		# 5th col of vcf is ALTERNATE allele
						INFO = re.split(';', lineList[7]) # parse INFO field, ; delimited
						AltDict[vcfFN][pos] = (alt, INFO[1].replace('AF=','') ) # 2nd col of INFO field contains allele frequency AF=#.###
	# Print output	
	outFile = open('InfoForMKtable_%s.txt' % chromToGenotype, 'w')
	print("Position", "RefAllele", end=" ", file=outFile)
	for vcfFN in sorted(vcfList):
		vcfFNsplit = vcfFN.split("/")
		name = vcfFNsplit[len(vcfFNsplit)-1]
		nameSplit = name.split("_")
		speciesName = nameSplit[0]
		print(speciesName, "_Allele ", speciesName, "_Freq",sep="", end=" ", file=outFile)
	print(file=outFile)
	for pos in RefDict:	# sites in which at least one sample has passing SNP
		# ignore multiallelic sites and indels
		#if pos not in MultiDict:
		print(pos, end=" ", file=outFile)
		print(RefDict[pos], end=" ", file=outFile)
		for vcfFN in sorted(vcfList):
			if pos in AltDict[vcfFN]:
				# remove square brackets PyVCF puts on ALT and AF
				x = re.sub('[\[\]]', '', str(AltDict[vcfFN][pos][0]))
				y = re.sub('[\[\]]', '', str(AltDict[vcfFN][pos][1]))
				print(x, y, file=outFile, end=" ")
			else:   # otherwise was still in CallableSitesDict
				print(RefDict[pos], 1.0, file=outFile, end=" ")
		print("", file=outFile)

if __name__ == '__main__':
  main()
