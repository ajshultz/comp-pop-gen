#!/usr/bin/env python

""" 
The first step with working with gffutils involves importing the gff file into a local sqlite3 file-based database.
It also checks whether the gff file contains gene-level information and whether the scaffold names in the gff file 
match those in the reference sequence.
"""

import sys
import re
import os
import gzip
import gffutils
import sharedFunctions as shFn

def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 2:
		configFile = sys.argv[1]
	else:
		print("need at least 1 command line arg, config file, exiting...")
		sys.exit()

	configInfo = shFn.parse_config(configFile)
	# get name of reference from refFile, used to rename gff file
	refName = re.search(r'\/(\w+)\.fa', configInfo["refFile"]) 
	refName = refName.group(1)

	"""
	if not os.path.isdir(configInfo["outDir"] + "/annotation_info"):
		os.mkdir(configInfo["outDir"] + "/annotation_info")
	"""
	#
	# store names of scaffolds in gff File
	#
	gffScaffNames = []
	gff = gzip.open(configInfo["gffFile"],'r')
	geneMatch = 0 	# used to see if there is any gene level information
	for line in gff:
		line = line.decode("utf-8")
		if not re.match(r'^#', line):
			lineList = re.split(r'\t+', line)
			name = lineList[0]
			gffScaffNames.append(name)
			if "mRNA" in line:
				#if re.match(r'mRNA', line):
				geneMatch += 1
	if geneMatch == 0:
		print("WARNING: \"mRNA\" not detected, this is not a proper gene-level gff, fix this before moving forwards, exiting...")
		sys.exit()

	#
	# store names of scaffolds in reference sequence
	#
	refScaffNames = []
	ref = open(configInfo["refFile"], 'r')
	for line in ref:
		if re.match(r'^>', line):
			line = line.strip()
			lineList = line.split()
			#lineList = re.split(r'\w+', line)
			chName = lineList[0]
			chName = chName.replace('>','')
			refScaffNames.append(chName)
	
	gffNameInRefseq=0
	for name in gffScaffNames:
		if name in refScaffNames:
			gffNameInRefseq += 1		
	f = open("numScaffsInRef.txt",'w')
	print(len(refScaffNames), file=f)
	f.close()

	if gffNameInRefseq == len(gffScaffNames): 
		db = gffutils.create_db(configInfo["gffFile"], dbfn="gff3Database.db", force=True, keep_order=False, merge_strategy='merge', sort_attribute_values=False)
		#os.system("mv gff3Database.db " + configInfo["outDir"] + "/annotation_info")
	else:
		print("There was a discrepancy between scaffold names in the gff file and the reference sequence.")
		print("Create a new gff file with matching scaffold names and continue. Exiting...")
		sys.exit()
	
	#
	# store different scaffold names, Refseq/Genbank
	#
	"""
	refseqNames = []
	genbankNames = []
	scaff = open(scaffNamesFile, 'r')
	for line in scaff:
		if not re.match(r'^#Assembly', line):
			lineList = re.split(r'\t+', line)
			refseqNames.append(lineList[2])
			genbankNames.append(lineList[3])
	"""
	
	#
	# do the scaffold names in the gff file match those in refseq or genbank?
	#
	"""
	gffInRefseqNames=0
	gffInGenbankNames=0
	# just check first name
	if gffScaffNames[0] in refseqNames:
		gffInRefseqNames += 1
	if gffScaffNames[0] in genbankNames:
		gffInGenbankNames += 1

	refInRefseqNames=0
	refInGenbankNames=0
	if refScaffNames[0] in refseqNames:
		refInRefseqNames += 1
	if refScaffNames[0] in genbankNames:
		refInGenbankNames += 1
	# Check if gff names are refseq or genbank, otherwise you need to expand list of possible naming schemes
	if not gffInRefseqNames and not gffInGenbankNames:
		print("WARNING: gff scaffold names do not follow Refseq or Genbank format here, cannot proceed...")
		sys.exit()
	if not refInRefseqNames and not refInGenbankNames:
		print("WARNING: ref scaffold names do not follow Refseq or Genbank format here, cannot proceed...")
		sys.exit()

	if gffInGenbankNames and refInGenbankNames:
		print("both in genbank")
	elif gffInRefseqNames and refInGenbankNames:
		print("gff in refseq but ref in genbank")
	sys.exit()
	"""

	"""
	# if gff scaffold names match those in reference sequence, make a renamed copy
	if not flag==0:
		if configInfo["gffFile"].endswith(".gz"):
			copyCommand = "cp " + configInfo["gffFile"] + " " + configInfo["outDir"] + "/annotation_info/" + refName + ".gff.gz"
		else:
			copyCommand = "cp " + configInfo["gffFile"] + " " + configInfo["outDir"] + "/annotation_info/" + refName + ".gff"
		os.system(copyCommand)
	else:
		print("gff file and reference sequence had different scaffold names")

	# 	IF GENBANK NAMES DONT MATCH, MAKE NEW GFF FILE WITH MATCHING NAMES, MAKE NOTE OF THIS IN FILE
	"""

if __name__ == '__main__':
  main()
