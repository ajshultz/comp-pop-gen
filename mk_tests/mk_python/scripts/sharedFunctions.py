# functions shared across multiple scripts

import sys
import gzip
import re
import os
import gffutils
import array as arr
from subprocess import Popen,PIPE

def parse_config(configFile):
	configInfo = {}
	f = open(configFile, 'r')
	for line in f:
		if re.match(r'\-\-', line): 
			line = line.strip()
			line = line.split(" ")
			# arg values specified by "--"
			if line[0] == "--OUT_DIR":
				configInfo["outDir"] = line[1]
			if line[0] == "--GFF":
				configInfo["gffFile"] = line[1]
			if line[0] == "--REF":
				configInfo["refFile"] = line[1]
			if line[0] == "--BED_VCF_DIR":
				configInfo["bedVcfDir"] = line[1]
			if line[0] == "--SPECIES":
				speciesList = line[1].split(',')
				configInfo["speciesList"] = speciesList
			if line[0] == "--SCRIPT_DIR":
				configInfo["scriptDir"] = line[1]
			if line[0] == "--SCAFF_PER_JOB":
				configInfo["scaffPerJob"] = int(line[1])
	f.close()
	return(configInfo)

def parseScaffNameFile(scaffNamesFile):
	rvDict = {}
	fh = open(scaffNamesFile,'r')
	for line in fh:
		if not re.match('^#', line):	# ignore headers
			l = line.split()
			RefseqName = l[2]
			GenbankName = l[3]
			rvDict[RefseqName] = GenbankName 
	return(rvDict)

def constructCdsDicts(db, chromToGenotype):
	CDSinfo = {}
	Pos_2_CDS = {}
	for cds in db.features_of_type('CDS'):
		if cds.seqid == chromToGenotype:
			#print(cds.id, " ", len(cds.id.split("_")))
			cdsIdNew = cds.id.replace('cds','')	# cdsIdNew = ##_##: left of _ represent mRNA it came from (may not 
								# be same # as rnd.id, but cds with same # came from same rna), 
								# right of _ represent exon number	
			x = (list(db.parents(cds, featuretype='gene')))
			if len(x) > 1:
				print("CDS ", cds.id, " had more than one gene ascribed to it")
			CDSinfo[cdsIdNew] = (x[0].id, cds.strand, cds.frame, cds.start, cds.end, x[0].seqid, x[0].attributes['Name'][0])
			#print(CDSinfo[cdsIdNew])
			for pos in range(cds.start, cds.end+1):	# +1 bc range inherently doesn't include specified endpoint
				try:
					Pos_2_CDS[pos].append(cdsIdNew)
				except:
					Pos_2_CDS[pos] = []
					Pos_2_CDS[pos].append(cdsIdNew)
	
	for pos in Pos_2_CDS.keys():
		Pos_2_CDS[pos].sort()
	"""
	for cds in CDSinfo.keys():
		print(cds, " ", CDSinfo[cds])
	for pos in Pos_2_CDS.keys():
		print(type(pos), pos, end=" ")
		for cds in Pos_2_CDS[pos]:
			print(cds, end=" ")
		print()
	"""
	return(Pos_2_CDS, CDSinfo)

def constructTranslateCodonDict():
    translateCodonDict = {}
    translateCodonDict["ATT"] = "I" ;
    translateCodonDict["ATC"] = "I" ;
    translateCodonDict["ATA"] = "I" ;
    translateCodonDict["CTT"] = "L" ;
    translateCodonDict["CTC"] = "L" ;
    translateCodonDict["CTA"] = "L" ;
    translateCodonDict["CTG"] = "L" ;
    translateCodonDict["TTA"] = "L" ;
    translateCodonDict["TTG"] = "L" ;
    translateCodonDict["GTT"] = "V" ;
    translateCodonDict["GTC"] = "V" ;
    translateCodonDict["GTA"] = "V" ;
    translateCodonDict["GTG"] = "V" ;
    translateCodonDict["TTT"] = "F" ;
    translateCodonDict["TTC"] = "F" ;
    translateCodonDict["ATG"] = "M" ;
    translateCodonDict["TGT"] = "C" ;
    translateCodonDict["TGC"] = "C" ;
    translateCodonDict["GCT"] = "A" ;
    translateCodonDict["GCC"] = "A" ;
    translateCodonDict["GCA"] = "A" ;
    translateCodonDict["GCG"] = "A" ;
    translateCodonDict["GGT"] = "G" ;
    translateCodonDict["GGC"] = "G" ;
    translateCodonDict["GGA"] = "G" ;
    translateCodonDict["GGG"] = "G" ;
    translateCodonDict["CCT"] = "P" ;
    translateCodonDict["CCC"] = "P" ;
    translateCodonDict["CCA"] = "P" ;
    translateCodonDict["CCG"] = "P" ;
    translateCodonDict["ACT"] = "T" ;
    translateCodonDict["ACC"] = "T" ;
    translateCodonDict["ACA"] = "T" ;
    translateCodonDict["ACG"] = "T" ;
    translateCodonDict["TCT"] = "S" ;
    translateCodonDict["TCC"] = "S" ;
    translateCodonDict["TCA"] = "S" ;
    translateCodonDict["TCG"] = "S" ;
    translateCodonDict["AGT"] = "S" ;
    translateCodonDict["AGC"] = "S" ;
    translateCodonDict["TAT"] = "Y" ;
    translateCodonDict["TAC"] = "Y" ;
    translateCodonDict["TGG"] = "W" ;
    translateCodonDict["CAA"] = "Q" ;
    translateCodonDict["CAG"] = "Q" ;
    translateCodonDict["AAT"] = "N" ;
    translateCodonDict["AAC"] = "N" ;
    translateCodonDict["CAT"] = "H" ;
    translateCodonDict["CAC"] = "H" ;
    translateCodonDict["GAA"] = "E" ;
    translateCodonDict["GAG"] = "E" ;
    translateCodonDict["GAT"] = "D" ;
    translateCodonDict["GAC"] = "D" ;
    translateCodonDict["AAA"] = "K" ;
    translateCodonDict["AAG"] = "K" ;
    translateCodonDict["CGT"] = "R" ;
    translateCodonDict["CGC"] = "R" ;
    translateCodonDict["CGA"] = "R" ;
    translateCodonDict["CGG"] = "R" ;
    translateCodonDict["AGA"] = "R" ;
    translateCodonDict["AGG"] = "R" ;
    translateCodonDict["TAA"] = "Stop" ;
    translateCodonDict["TAG"] = "Stop" ;
    translateCodonDict["TGA"] = "Stop" ;
    
    return (translateCodonDict) ;

def getCodonPositionsInReference_forwardStrand(posInScaff, firstCDS, CDSinfo):
	rf = int(CDSinfo[firstCDS][2]) # rf = reading frame
	cdsStartPos = int(CDSinfo[firstCDS][3])
	cdsEndPos = int(CDSinfo[firstCDS][4])
	flankBase1 = None
	flankBase2 = None
	posInCDS = (posInScaff - cdsStartPos)+1
	posInCodon = (posInCDS - rf)%3	# = 0,1,2 if 3rd,1st,2nd codon position, respectively
	# you have codon position of mutation, now you need to get flanking bases from reference
	# if codon is at border of exon, need to account for introns
	if posInCodon == 0:	# 3rd codon position
				# are previous 2 bases in CDS? if so take prev. 2 bases from reference
		if posInScaff-2 >= cdsStartPos:
			#take posInScaff - 1, posInScaff - 2 from ref genome
			flankBase1 = posInScaff-1
			flankBase2 = posInScaff-2
			#print("bothInScaff ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
		else:
			# need to get positions from previous exon (5')
			cdsPrev = getFlankingCDS(firstCDS, "previous", CDSinfo)
			if cdsPrev == None:
				return(None, None, None)
			cdsPrevStartPos = int(CDSinfo[cdsPrev][3]) 
			cdsPrevEndPos = int(CDSinfo[cdsPrev][4])
			#print(cdsStartPos, " ", posInScaff, " ", firstCDS, " ", cdsPrev)
			if posInScaff-1 < cdsStartPos:	
				# 3rd codon position at leftmost point in exon
				# get previous position from previous exon 
				flankBase1 = cdsPrevEndPos 
				# check if next position also in previous codon; may not be if 1bp long!
				if cdsPrevEndPos-1 >= cdsPrevStartPos:
					flankBase2 = cdsPrevEndPos-1 
				else:
					cdsPrev2 = getFlankingCDS(cdsPrev, "previous", CDSinfo)
					if cdsPrev2 == None:
						return(None, None, None)
					cdsPrev2EndPos = int(CDSinfo[cdsPrev2][4])
					flankBase2 = cdsPrev2EndPos 
				#print("3rdCodonLeftmostPos: ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
			elif posInScaff-2 < cdsStartPos:
				# 2nd codon pos in curr. exon, only 1st codon pos in previous exon to (5')
				flankBase1 = posInScaff-1
				flankBase2 = int(CDSinfo[cdsPrev][4])
				#print("1stCodonPosInAnotherExon: ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
			else:
				print("WARNING: CODON MAYHEM")
				sys.exit()
		# codon = flankBase2, flankBase1, posInScaff 
		if flankBase2 < flankBase1 and flankBase1 < posInScaff:
			return(flankBase2, flankBase1, posInScaff)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_forwardStrand, exiting...")
			sys.exit()
	elif posInCodon == 1:	# 1st codon position
		if posInScaff+2 <= cdsEndPos:	
			# other bases in codon in same exon
			flankBase1 = posInScaff+1
			flankBase2 = posInScaff+2
		else:
			cdsNext = getFlankingCDS(firstCDS, "next", CDSinfo)
			if cdsNext == None:
				return(None, None, None)
			cdsNextStartPos = int(CDSinfo[cdsNext][3])
			cdsNextEndPos = int(CDSinfo[cdsNext][4])
			if posInScaff+1 > cdsEndPos:
				#1st codon position at rightmost point in exon
				flankBase1 = cdsNextStartPos
				if cdsNextStartPos+1 <= cdsNextEndPos:
					flankBase2 = cdsNextStartPos+1
				else:
					cdsNext2 = getFlankingCDS(cdsNext, "next", CDSinfo)
					if cdsNext2 == None:
						return(None, None, None)
					cdsNext2StartPos = int(CDSinfo[cdsNext2][3])
					flankBase2 = cdsNext2StartPos
			elif posInScaff+2 > cdsEndPos:
				flankBase1 = posInScaff+1
				flankBase2 = cdsNextStartPos 
			else:
				print("WARNING: CODON MAYHEM")
				sys.exit()
		# codon = posInScaff, flankBase1, flankBase2
		if posInScaff < flankBase1 and flankBase1 < flankBase2:
			return(posInScaff, flankBase1, flankBase2)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_forwardStrand, exiting...")
			sys.exit()
	elif posInCodon == 2:	# 2nd codon position
		if posInScaff-1 >= cdsStartPos:
			# entire codon in current exon
			flankBase1 = posInScaff-1
		else:
			cdsPrev = getFlankingCDS(firstCDS, "previous", CDSinfo)
			if cdsPrev == None:
				return(None, None, None)
			cdsPrevEndPos = int(CDSinfo[cdsPrev][4])
			flankBase1 = cdsPrevEndPos 
		if posInScaff+1 <= cdsEndPos:
			flankBase2 = posInScaff+1
		else:
			cdsNext = getFlankingCDS(firstCDS, "next", CDSinfo)
			if cdsNext == None:
				return(None, None, None)
			cdsNextStartPos = int(CDSinfo[cdsNext][3]) 
			flankBase2 = cdsNextStartPos

		if flankBase1 == None or flankBase2 == None:
			print("WARNING: CODON MAYHEM")
			sys.exit()
		# codon = flankBase1, posInScaff, flankBase2
		if flankBase1 < posInScaff and posInScaff < flankBase2:
			return(flankBase1, posInScaff, flankBase2)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_forwardStrand, exiting...")
			sys.exit()


def getCodonPositionsInReference_reverseStrand(posInScaff, firstCDS, CDSinfo):
	rf = int(CDSinfo[firstCDS][2]) # rf = reading frame
	cdsStartPos = int(CDSinfo[firstCDS][3])
	cdsEndPos = int(CDSinfo[firstCDS][4])
	flankBase1 = None
	flankBase2 = None
	posInCDS = (cdsEndPos - posInScaff)+1
	posInCodon = (posInCDS - rf)%3	# = 0,1,2 if 3rd,1st,2nd codon position, respectively
	if posInCodon == 0:	# 3rd codon position
				# are next 2 bases in CDS? These are other bases in codon since gene seq reversed on - strand 
		if posInScaff+2 <= cdsEndPos:	# where cdsEndPos is the start of the exon, b/c gene on - strand
			#take posInScaff - 1, posInScaff - 2 from ref genome
			flankBase1 = posInScaff+1
			flankBase2 = posInScaff+2
			#print("bothInScaff ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
		else:
			# need to get positions from exon to left (5')
			cdsPrev = getFlankingCDS(firstCDS, "previous", CDSinfo)	# left with respect to exon sequence, but actually moving right along ref genome
			if cdsPrev == None:
				return(None, None, None)
			cdsPrevStartPos = int(CDSinfo[cdsPrev][3]) 
			cdsPrevEndPos = int(CDSinfo[cdsPrev][4])
			if posInScaff+1 > cdsEndPos:	
				# 3rd codon position at very beginning of exon
				# get other 2 positions of codon from the start coord. of previous exon (corresponds to the 3' end of that exon) 
				flankBase1 = cdsPrevStartPos 
				if cdsPrevStartPos+1 <= cdsPrevEndPos:
					flankBase2 = cdsPrevStartPos+1
				else:
					cdsPrev2 = getFlankingCDS(cdsPrev, "previous", CDSinfo)
					if cdsPrev2 == None:
						return(None, None, None)
					cdsPrev2StartPos = int(CDSinfo[cdsPrev2][3])
					flankBase2 = cdsPrev2StartPos 
					
				#print("3rdCodonLeftmostPos: ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
			elif posInScaff+2 > cdsEndPos:
				# 2nd codon pos in curr. exon, 1st codon pos in previous exon 
				flankBase1 = posInScaff+1
				flankBase2 = cdsPrevStartPos 
				#print("1stCodonPosInAnotherExon: ", flankBase2, " ", flankBase1, " ", "BASE:", posInScaff)
			else:
				print("WARNING: CODON MAYHEM")
				sys.exit()
		# codon = posInScaff, flankBase1, flankBase2
		if posInScaff < flankBase1 and flankBase1 < flankBase2:
			return(posInScaff, flankBase1, flankBase2)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_reverseStrand, exiting...")

	elif posInCodon == 1:	# 1st codon position
		if posInScaff-2 >= cdsStartPos:	
			# other bases in codon in same exon
			flankBase1 = posInScaff-1
			flankBase2 = posInScaff-2
		else:
			cdsNext = getFlankingCDS(firstCDS, "next", CDSinfo)
			if cdsNext == None:
				return(None, None, None)
			cdsNextStartPos = int(CDSinfo[cdsNext][3])
			cdsNextEndPos = int(CDSinfo[cdsNext][4])
			if posInScaff-1 < cdsStartPos:
				#1st codon position at very end position (3') of exon 
				flankBase1 = cdsNextEndPos
				if cdsNextEndPos-1 >= cdsNextStartPos:
					flankBase2 = cdsNextEndPos-1
				else:
					cdsNext2 = getFlankingCDS(cdsNext, "next", CDSinfo)
					if cdsNext2 == None:
						return(None, None, None)
					cdsNext2EndPos = int(CDSinfo[cdsNext2][3])
					flankBase2 = cdsNext2EndPos
			elif posInScaff-2 < cdsStartPos:
				flankBase1 = posInScaff-1
				flankBase2 = cdsNextEndPos 
			else:
				print("WARNING: CODON MAYHEM")
				sys.exit()
		# codon = flankBase2, flankBase1, posInScaff
		if flankBase2 < flankBase1 and flankBase1 < posInScaff:
			return(flankBase2, flankBase1, posInScaff)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_reverseStrand, exiting...")
	elif posInCodon == 2:	# 2nd codon position
		if posInScaff+1 <= cdsEndPos:
			flankBase1 = posInScaff+1
		else:
			cdsPrev = getFlankingCDS(firstCDS, "previous", CDSinfo)	
			if cdsPrev == None:
				return(None, None, None)
			cdsPrevStartPos = int(CDSinfo[cdsPrev][3]) 
			flankBase1 = cdsPrevStartPos 
		if posInScaff-1 >= cdsStartPos:
			flankBase2 = posInScaff-1
		else:
			cdsNext = getFlankingCDS(firstCDS, "next", CDSinfo)	
			if cdsNext == None:
				return(None, None, None)
			cdsNextEndPos = int(CDSinfo[cdsNext][4]) 
			flankBase2 = cdsNextEndPos 

		if flankBase1 == None or flankBase2 == None:
			print("WARNING: CODON MAYHEM")
			sys.exit()
		# codon = flankBase2, posInScaff, flankBase1
		if flankBase2 < posInScaff and posInScaff < flankBase1:
			return(flankBase2, posInScaff, flankBase1)
		else:
			print("Unsorted codon positions in getCodonPositionsInReference_reverseStrand, exiting...")

def getFlankingCDS(cds, direction, CDSinfo):
	cdsName = cds.split("_")
	cdsRNA = cdsName[0]
	cdsNum = None
	if len(cdsName) == 2:
		cdsNum = int(cdsName[1])	# convert to int to increment/decrement
	else:
		cdsNum = 0			# no "_" present, signaling first exon of mRNA

	if direction == "previous":
		cdsNum -=1
	elif direction == "next":
		cdsNum +=1
	else:
		print("getFlankingCDS function not used properly, exiting...")
		sys.exit()

	newcds = cdsRNA + "_" + str(cdsNum)
	if cdsNum == 0:
		return(cdsRNA)	# 1st exon not followed by "_#"
	elif cdsNum >= 1 and newcds in CDSinfo: # if annotation errors, cdsNum could be greater than 3' most exon; opposite error of getting negative cdsNum
		return(newcds)
	else:
		print("WARNING: couldn't get flanking bases for polymorphic position, not able to locate flanking exon")
		return(None)

def getPosInCodon(p1,p2,p3,posOfInterest):
	posInCodonToChange = None
	index = -1
	for i in (p1,p2,p3):
		index+=1
		if i == posOfInterest:
			posInCodonToChange = index
	return(posInCodonToChange)

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

def getBedFiles(vcfBedDir, speciesList):
	bedList = []
	for species in speciesList:
		fn = vcfBedDir + "/" + species + "/stats_coverage/" + "_" + species + "_clean_coverage_sites_merged.bed"
		exists = os.path.isfile(fn)
		if exists:
			bedList.append(fn)	
		else:
			print("one of your bed files doesn't exist")
			sys.exit()
	return(bedList)

def create_job_script(name, outDir, queue, nCores, time, mem, command):
	outFile = open('job_%s.sh' % name , 'w')
	print("#!/bin/bash", file=outFile)
	print("#SBATCH -J "+ str(name), file=outFile)
	print("#SBATCH -o out." + str(name), file=outFile)
	print("#SBATCH -e err." + str(name), file=outFile)
	print("#SBATCH -p " + queue, file=outFile)
	print("#SBATCH -n " + str(nCores), file=outFile)
	print("#SBATCH -t " + str(time), file=outFile)
	print("#SBATCH --mem=" + str(mem), file=outFile)
	print(command, file=outFile)
	outFile.close()
	jobId = sbatch_submit(outFile.name)
	print(jobId)
	os.system("mv job_" + str(name) + ".sh " + outDir)
	return(jobId)

#Submit filename to slurm with sbatch, returns job id number
def sbatch_submit(filename):
	proc = Popen('sbatch %s'%filename,shell=True,stdout=PIPE,stderr=PIPE)
	stdout,stderr = proc.communicate()
	stdout = stdout.decode("utf-8","ignore")
	stdout = stdout.strip()
	stdout = stdout.strip('Submitted batch job ')
	return(stdout)

#Check job status of specific jobid: returns job status
def jobid_status(jobid,date):
	print(jobid)
	proc = Popen('sacct --format state --noheader -j %d'%(int(jobid)),shell=True,stdout=PIPE,stderr=PIPE)
	#proc = Popen('sacct --format state --noheader -j %d -S %s'%(int(jobid),date),shell=True,stdout=PIPE,stderr=PIPE)
	stdout,stderr=proc.communicate()
	#print(stdout)
	#print(stderr)
	if proc.returncode != 0:
		raise Exception('Error running sacct: %s'%stderr)
	if stdout.strip() == '':
		return("No Status")
	lines = stdout.split()
	print(lines)
	if not lines:
		rv = None
		print("InFunction ", jobid, " ", rv)
		return(rv)
	else:
		rv = lines[0].decode("utf-8","ignore")
		print("InFunction ", jobid, " ", rv)
		return(rv)


