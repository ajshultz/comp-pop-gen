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
from time import sleep
import datetime
import sharedFunctions as shFn

def make_outdir(speciesList, outDir):
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
	outdirNameSBATCH = outDir + "/" + "sbatchScripts"
	if not os.path.isdir(outdirNameAT):
		try:
			os.mkdir(outdirNameAT)
		except:
			print("Couldn't make output directory, exiting...")
			sys.exit()
	if not os.path.isdir(outdirNameMK):
		try:
			os.mkdir(outdirNameMK)
		except:
			print("Couldn't make output directory, exiting...")
			sys.exit()
	if not os.path.isdir(outdirNameSBATCH):
		try:
			os.mkdir(outdirNameSBATCH)
		except:
			print("Couldn't make output directory, exiting...")
			sys.exit()
	return(outdirNameMK, outdirNameAT, outdirNameSBATCH)

def GetChrom2IntervalNum(speciesList, vcfBedDir):
	Chrom2IntervalMap = {}
	for species in speciesList:
		Chrom2IntervalMap[species] = {}
		# list interval files
		intervalsDir = vcfBedDir + "/" + species + "/genome/*" + "_splits_interval_lists/*"
		lsCommand = "ls " + intervalsDir
		intervalFilesStr = os.popen(lsCommand).read()
		intervalFilesList = intervalFilesStr.split()
		# go thru each interval file, get interval number from file name, and chromosomes for that interval number
		for fn in intervalFilesList:
			fh = open(fn, 'r')
			intervalNum = re.search(r'_(\d+)\.interval_list', fn)
			for chrom in fh:
				chrom = chrom.strip()
				Chrom2IntervalMap[species][chrom] = intervalNum.group(1)
	return(Chrom2IntervalMap)

#Create job script


def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 2:
		configFile = sys.argv[1]
	else:
		print("USAGE: Please specify VCF/BED directory and comma-separated (no spaces) list of species")
		sys.exit()

	configInfo = shFn.parse_config(configFile)
	script1 = configInfo["scriptDir"] + "/01a_VcfBed_2_AlleleTable_ByChrom_LowerMem.py"
	script2 = configInfo["scriptDir"] + "/02a_AlleleTable_2_MKtable.py"
	refSpecies = configInfo["speciesList"][0]
	gffDB = configInfo["outDir"] + "/gff3Database.db"
	if not os.path.isfile(script1):
		print("can't locate ", script1)
		sys.exit()
	if not os.path.isfile(script2):
		print("can't locate ", script2)
		sys.exit()
	if not os.path.isfile(gffDB):
		print("can't locate ", gffDB)
		sys.exit()
	if not os.path.isfile(configInfo["refFile"]):
		print("can't locate ", configInfo["refFile"])
		sys.exit()

	now = datetime.datetime.now()
	print('Staring work: %s'%now)
	start_date = now.strftime("%Y-%m-%d")
	# resources to request for jobs, SLURM directives
	#queue = "serial_requeue,shared"
	queue = "serial_requeue"
	nCores = 1
	time = 4800
	mem = 15000
	chromPerJob = configInfo["scaffPerJob"]

	# Make directory for output files
	outDir_mkTable, outDir_alleleTable, outDir_sbatchScript = make_outdir(configInfo["speciesList"], configInfo["outDir"])
	# get bed files
	bedList = shFn.getBedFiles(configInfo["bedVcfDir"], configInfo["speciesList"])

	# For each chromosome, find out which VCF it's in, put it in Chrom2IntervalMap
	# get interval files for each vcf (i.e. which chromosomes a particular vcf contains)
	Chrom2IntervalMap = GetChrom2IntervalNum(configInfo["speciesList"], configInfo["bedVcfDir"])
	# get list of chromosome names
	ChromList = shFn.getChromList(bedList)


	#Submit one job per chrom
	bedCommaSep = ','.join(bedList)
	chromTotal = 0
	chromCount = 0
	command = ""
	tmpCount = 0
	jobNum = 0
	jobIdCommandDict = {}
	currentCoreMemAllocDict = {}
	currentMemAllocDict = {}
	for chrom in sorted(ChromList):
		#if tmpCount <=50:
		#tmpCount +=1
		chromTotal += 1
		chromCount += 1
		vcfList = []
		alleleTableFile = outDir_alleleTable + "/InfoForMKtable_" + chrom + ".txt" # this gets fed into MK table script, below
		# get appropriate vcfs, for first script that computes allele tables
		for species in configInfo["speciesList"]:
			#use Chrom2IntervalMap Dict to see which VCF to use for this chrom for this species
			vcf = configInfo["bedVcfDir"] + "/" + species + "/vcf/" + species + "_hardfilters." + Chrom2IntervalMap[species][chrom] + ".vcf.gz"
			vcfList.append(vcf)
			existsVcf = os.path.isfile(vcf)
			if not existsVcf:
				print("one of your vcf or bed files doesn't exist")
				sys.exit()
		vcfCommaSep = ','.join(vcfList)

		if chromCount <= chromPerJob:	# add another chromosome to current job
			# First calculate allele table, then MK table
			command += "cd " + outDir_alleleTable + "\n"
			command += "python3 " + script1 + " " + vcfCommaSep + " " + bedCommaSep + " " + gffDB + " " + configInfo["refFile"] + " " + chrom + "\n"
			command += "cd " + outDir_mkTable + "\n"
			command += "python3 " + script2 + " " + gffDB + " " + configInfo["refFile"] + " " + alleleTableFile + " " + chrom + "\n"
		if chromCount > chromPerJob or chromTotal == len(ChromList):
			# create SBATCH script and submit
			#print(command)
			#sys.exit()
			jobNum += 1
			jobId = shFn.create_job_script(jobNum, outDir_sbatchScript, queue, nCores, time, mem, command)
			jobIdCommandDict[jobId] = command
			currentCoreMemAllocDict[jobId] = nCores
			currentMemAllocDict[jobId] = mem
			# reset command to submit, chromCount
			command = "cd " + outDir_alleleTable + "\n"
			command += "python3 " + script1 + " " + vcfCommaSep + " " + bedCommaSep + " " + gffDB + " " + configInfo["refFile"] + " " + chrom + "\n"
			command += "cd " + outDir_mkTable + "\n"
			command += "python3 " + script2 + " " + gffDB + " " + configInfo["refFile"] + " " + alleleTableFile + " " + chrom + "\n"
			chromCount = 1
	ongoingJobs = ['RUNNING', 'PENDING']
	failedJobs = ['CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL','OUT_OF_ME+']
	sleep(300)	# allow time for jobs to submit, otherwise they have no status and program exits

	numSuccessJobsDict = {}
	while len(numSuccessJobsDict.keys()) < len(jobIdCommandDict.keys()):
		for jobId in jobIdCommandDict:
			status = shFn.jobid_status(jobId, start_date)
			print(status)
			if status and status not in ongoingJobs:
				if status == 'COMPLETED':
					numSuccessJobsDict[jobId] = 1
				elif status in failedJobs:
					if status == 'OUT_OF_ME+':
						# Submit job with 2X more requested resources as last time, memory and cores
						chromTotal +=1
						nCoresMore = currentCoreMemAllocDict[jobId]*1
						memMore = currentMemAllocDict[jobId]*2
						command = jobIdCommandDict[jobId]
						jobIdNew = shFn.create_job_script(jobNum, outDir_sbatchScript, queue, nCoresMore, time, memMore, command)
						# update jobId parameters, so program keeps requesting more cores and memory
						currentCoreMemAllocDict[jobIdNew] = nCoresMore
						currentMemAllocDict[jobIdNew] = memMore
						jobIdCommandDict[jobIdNew] = command
						# since this job failed, remove details
						del currentCoreMemAllocDict[jobId]
						del currentMemAllocDict[jobId]
						del jobIdCommandDict[jobId]
					else:
						# just resubmit job with same resource requests
						chromTotal +=1
						command = jobIdCommandDict[jobId]
						jobIdNew = shFn.create_job_script(jobNum, outDir_sbatchScript, queue, nCores, time, mem, command)
						jobIdCommandDict[jobIdNew] = command
						del jobIdCommandDict[jobId]
				else:
					print("UNIDENTIFIED JOB STATUS: ", status)

		sleep(120)
	print("FINISHED")
	sys.exit()


if __name__ == '__main__':
  main()
