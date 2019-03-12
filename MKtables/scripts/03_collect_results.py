#!/usr/bin/python -tt

"""
This wrapper takes in a parent directory and comma-separated list of species, and invokes VcfBed_2_AlleleTable_ByChrom.py script
by locating the correct vcf files to use for a particular chromosome, and locating the bed files to see which sites had sufficient
coverage
IMPORTANT: This code assumes that the set of species supplied to the wrapper function were all aligned to the same reference genome,
since the same chromosome needs to be present in the vcf files of all species analyzed.
"""

import sys
import re
import os
import sharedFunctions as shFn


def get_outdirs(speciesList, outDir):
	outdirNameMK = ""
	outdirNameAT = ""
	for index in range(len(speciesList)):
		tmp = speciesList[index].split("_")
		if index <= len(speciesList)-1:
			outdirNameMK += tmp[0] + "_"
			outdirNameAT += tmp[0] + "_"
	outdirNameMK += "MKtables"
	outdirNameAT += "AlleleTables"
	outdirNameMK = outDir + "/" + outdirNameMK
	outdirNameAT = outDir + "/" + outdirNameAT
	return(outdirNameMK, outdirNameAT)

def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 2:
		configFile = sys.argv[1]
	else:
		print("USAGE: Please specify right number of arguments")
		sys.exit()

	configInfo = shFn.parse_config(configFile)
	# get directory names of output files
	outDir_mkTable, outDir_alleleTable = get_outdirs(configInfo["speciesList"], configInfo["outDir"])

	# get list of chromosomes, as these are in the names of files
	bedList = shFn.getBedFiles(configInfo["bedVcfDir"], configInfo["speciesList"])
	ChromList = shFn.getChromList(bedList) # list of chromosome names

	covStatsDict = {}
	mkTableDict = {}
	for chrom in sorted(ChromList):
		covStats_fileName = outDir_alleleTable + "/coverageStats_" + chrom + ".txt"
		alleleTable_fileName = outDir_alleleTable + "/InfoForMKtable_" + chrom + ".txt"
		mkTable_fileName = outDir_mkTable + "/MKtable_" + chrom + ".txt"
		existsCovStats = os.path.isfile(covStats_fileName)
		existsMK = os.path.isfile(mkTable_fileName)
		if existsCovStats and existsMK:
			f = open(covStats_fileName, 'r')
			for line in f:
				line = line.split()
				if re.match(r'^GeneName', line[0]):
					covStatsDict["header"] = line[1:len(line)]
				else:
					covStatsDict[line[0]] = line[1:len(line)]
			f.close()
			f = open(mkTable_fileName, 'r')
			for line in f:
				line = line.split()
				if re.match(r'^GeneName', line[0]):
					mkTableDict["header"] = line[1:len(line)]
				else:
					mkTableDict[line[0]] = line[1:len(line)]
			f.close()

		else:
			print("Coverage stats and/or MK table didn't exist for " + chrom)
			continue

	sp_abbr = configInfo["speciesList"][0]

	outFile = open("%s_PolymorphismDivergenceStats_combined.txt"%(sp_abbr), 'w')
	print("Gene", end="\t", file=outFile)
	for i in covStatsDict["header"]:
		print(i, end="\t", file=outFile)
	for i in mkTableDict["header"]:
		print(i, end="\t", file=outFile)
	print(file=outFile)
	for gene in covStatsDict.keys():
		if gene != "header":
			print(gene, end="\t", file=outFile)
			for i in covStatsDict[gene]:
				print(i, end="\t", file=outFile)
			for i in mkTableDict[gene]:
				print(i, end="\t", file=outFile)
			print(file=outFile)

	sys.exit()


if __name__ == '__main__':
  main()
