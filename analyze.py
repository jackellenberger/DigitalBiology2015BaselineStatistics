#! /usr/bin/env python

import string
import os
import sys
import cStringIO
import math
import numpy as np
from pprint import pprint
from prody import *
from pylab import *
import random
ion()

################################################################################
# CONSTANTS / PARAMETERS #
##########################

HBOND_CA_DISTANCE = 6.5

polarity = {"ALA": ["CB"],
			"ARG": ["CB", "CG"],
			"ASN": ["CB"],
			"ASP": ["CB"],
			"CYS": ["CB"],
			"GLN": ["CB", "CG"],
			"GLU": ["CB", "CG"],
			"HIS": ["CB"],
			"ILE": ["CB", "CG1", "CG2", "CD1"],
			"LEU": ["CB", "CG", "CD1", "CD2"],
			"LYS": ["CB", "CG", "CD"],
			"MET": ["CB", "CG", "CE"],
			"PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
			"PRO": ["CB", "CG"],
			"THR": ["CG2"],
			"TRP": ["CB", "CG", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
			"TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2"],
			"VAL": ["CB", "CG1", "CG2"]}


################################################################################
# EDITED PRODY FUNCTIONS #
##########################

################################################################################
# Parsing #
###########

# "filename,ca1X,ca1Y,ca1Z,ca2X,ca2Y,ca2Z,centerX,centerY,centerZ,atomsInDD,totalAtoms"
# "filename,HBcount,HBcountLessSD,dehydronCount"
def analyze(filename,perDehydronOut,perHbondOut,perFileModel,perFileOutput,modeledAtomsInDD):
	slash = filename.rfind("/") + 1
	if slash != -1:
		truefilename=filename[slash:]
	else: truefilename = filename
	ppdb = parsePDB(filename,model=1)
	if ppdb.numAtoms() > 10000:
		print "abort! too many atoms! ain't nobody got time fo dat"
		return
	CAs = ppdb.calpha
	CAneighbors = findNeighbors(CAs,HBOND_CA_DISTANCE)
	if len(CAneighbors) > 4000:
		print "abort, too many neighbors implies too much computation time"
		return
	HBcount = 0
	HBList = []
	dehydronCount = 0
	AtomCountInDDD = 0
	AtomCountInHBDD = 0
	atomCountInDDD = 0
	uniqueAtomCountInHBDD = 0
	totalNeighbors=[]
	totalDDDneighbors=[]
	dehydronCAstring = []

	for (CA1, CA2, CAdist) in CAneighbors:

		CA1res = Atom.getResnum(CA1)
		if CA1res < 1:
			break
		CA1group = ppdb.select('resnum '+str(CA1res))
		CA1Hydrogens = CA1group.select('name HA or name 1HA or name 2HA or name H')
		CA1Nitrogens = (CA1group.select('name N'))
		CA1Oxygens = (CA1group.select('name O'))
		CA1B = None
		for atom in CA1Oxygens:
			CA1B = (CA1group.select('index '+str(Atom.getIndex(atom)-1)+\
								' or index '+str(Atom.getIndex(atom)+1)))
		
		CA2res = Atom.getResnum(CA2)
		if CA2res < 1:
			break
		CA2group = ppdb.select('resnum '+str(CA2res))
		CA2Hydrogens = (CA2group.select('name HA or name 1HA or name 2HA or name H'))
		CA2Nitrogens = (CA2group.select('name N'))
		CA2Oxygens = (CA2group.select('name O'))
		CA2B = None
		for atom in CA2Oxygens:
			CA2B = (CA2group.select('index '+str(Atom.getIndex(atom)-1)+\
								' or index '+str(Atom.getIndex(atom)+1)))
		initialHBCount = 0
		NPCcount = 0

		if (CA1Hydrogens and CA1Nitrogens and CA2Oxygens and CA2B):
			for D in CA1Nitrogens:
				for A in CA2Oxygens:
					DAdist = calcDistance(D,A)
					if DAdist < 3.5:
						for H in CA1Hydrogens:
							HAdist = calcDistance(H,A)
							if HAdist < 2.5:
								DHAangle = calcAngle(D,H,A)
								if DHAangle > 90:
									for B in CA2B:
										DABangle = calcAngle(D,A,B)
										if DABangle > 90:
											HABangle = calcAngle(H,A,B)
											if HABangle > 90:
												NPCcount,uniqueAtomCount,HBDDneighbors = countNPC(ppdb,CA1,CA2,totalNeighbors)
												AtomCountInHBDD += uniqueAtomCount
												initialHBCount += 1
												HBList.append((CA1,CA2))
												if NPCcount < 19:
													totalDDDneighbors += HBDDneighbors
													atomCountInDDD += uniqueAtomCount
													dehydronCount += 1
													printPerDehydron(perDehydronOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount)
													printPerHydrogenBond(perHbondOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount,True)
												else:
													printPerHydrogenBond(perHbondOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount,False)
												break

		if (CA2Hydrogens and CA2Nitrogens and CA1Oxygens and CA1B):
			for D in CA2Nitrogens:
				for A in CA1Oxygens:
					DAdist = calcDistance(D,A)
					if DAdist < 3.5:
						for H in CA2Hydrogens:
							HAdist = calcDistance(H,A)
							if HAdist < 2.5:
								DHAangle = calcAngle(D,H,A)
								if DHAangle > 90:
									for B in CA1B:
										DABangle = calcAngle(D,A,B)
										if DABangle > 90:
											HABangle = calcAngle(H,A,B)
											if HABangle > 90:
												NPCcount,uniqueAtomCount,HBDDneighbors = countNPC(ppdb,CA1,CA2,totalNeighbors)
												AtomCountInHBDD += uniqueAtomCount
												initialHBCount += 1
												HBList.append((CA1,CA2))
												if NPCcount < 19:
													totalDDDneighbors += HBDDneighbors
													atomCountInDDD += uniqueAtomCount
													dehydronCount += 1
													printPerDehydron(perDehydronOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount)
													printPerHydrogenBond(perHbondOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount,True)
												else:
													printPerHydrogenBond(perHbondOut,truefilename,CA1,CA2,NPCcount,uniqueAtomCount,False)
												break
		HBcount += initialHBCount
	uniqueAtomCountInHBDD = len(set(totalNeighbors))
	uniqueAtomCountInDDD = len(set(totalDDDneighbors))
	printPerFile(perFileOutput,truefilename,HBcount,dehydronCount,len(ppdb),AtomCountInHBDD,uniqueAtomCountInHBDD,atomCountInDDD,uniqueAtomCountInDDD)
	modelSampledHBList(truefilename,perFileModel,ppdb,HBcount,HBList,dehydronCount,modeledAtomsInDD,uniqueAtomCountInDDD)


def modelSampledHBList(filename,perFileModel,pdb,HBcount,HBList,dehydronCount,modeledAtomsInDD,uniqueAtomCountInDDD):
	sampledHBList = random.sample(HBList, dehydronCount)
	uniqueAtomsInDDD = []
	for (CA1,CA2) in sampledHBList:
		CA1coords = Atom.getCoords(CA1)
		CA2coords = Atom.getCoords(CA2)
		bondCenter = [(CA1coords[0]+CA2coords[0])/2.0,
				  	  (CA1coords[1]+CA2coords[1])/2.0,
				  	  (CA1coords[2]+CA2coords[2])/2.0]
		NPCcount, uniqueAtomCount, atomsInDDD = countNPC(pdb,CA1,CA2,[])

		perFileModel.write(str(filename)+","+\
						   str(Atom.getResnum(CA1))+","+\
						   str(CA1coords[0])+","+str(CA1coords[1])+","+str(CA1coords[2])+","+\
						   str(Atom.getResnum(CA2))+","+\
						   str(CA2coords[0])+","+str(CA2coords[1])+","+str(CA2coords[2])+","+\
						   str(bondCenter[0])+","+str(bondCenter[1])+","+str(bondCenter[2])+","+\
						   str(Atom.getResname(CA1))+","+\
						   str(Atom.getResname(CA2))+","+\
						   str(NPCcount)+","+\
				 		   str(uniqueAtomCount)+"\n")
		uniqueAtomsInDDD += atomsInDDD
	numUniqueAtomsInDDD = len(set(uniqueAtomsInDDD))
	modeledAtomsInDD.write(filename+","+str(numUniqueAtomsInDDD)+","+str(uniqueAtomCountInDDD)+"\n")

	



def printPerDehydron(output,filename,CA1,CA2,NPCcount,uniqueAtomCount):
	CA1coords = Atom.getCoords(CA1)
	CA2coords = Atom.getCoords(CA2)
	bondCenter = [(CA1coords[0]+CA2coords[0])/2.0,
				  (CA1coords[1]+CA2coords[1])/2.0,
				  (CA1coords[2]+CA2coords[2])/2.0]
	output.write(str(filename)+","+\
				 str(Atom.getResnum(CA1))+","+\
				 str(CA1coords[0])+","+str(CA1coords[1])+","+str(CA1coords[2])+","+\
				 str(Atom.getResnum(CA2))+","+\
				 str(CA2coords[0])+","+str(CA2coords[1])+","+str(CA2coords[2])+","+\
				 str(bondCenter[0])+","+str(bondCenter[1])+","+str(bondCenter[2])+","+\
				 str(Atom.getResname(CA1))+","+\
				 str(Atom.getResname(CA2))+","+\
				 str(NPCcount)+","+\
				 str(uniqueAtomCount)+"\n")

def printPerHydrogenBond(output,filename,CA1,CA2,NPCcount,uniqueAtomCount,isDehydron):
	CA1coords = Atom.getCoords(CA1)
	CA2coords = Atom.getCoords(CA2)
	bondCenter = [(CA1coords[0]+CA2coords[0])/2.0,
				  (CA1coords[1]+CA2coords[1])/2.0,
				  (CA1coords[2]+CA2coords[2])/2.0]
	output.write(str(filename)+","+\
				 str(Atom.getResnum(CA1))+","+\
				 str(CA1coords[0])+","+str(CA1coords[1])+","+str(CA1coords[2])+","+\
				 str(Atom.getResnum(CA2))+","+\
				 str(CA2coords[0])+","+str(CA2coords[1])+","+str(CA2coords[2])+","+\
				 str(bondCenter[0])+","+str(bondCenter[1])+","+str(bondCenter[2])+","+\
				 str(Atom.getResname(CA1))+","+\
				 str(Atom.getResname(CA2))+","+\
				 str(NPCcount)+","+\
				 str(uniqueAtomCount)+","+\
				 str(isDehydron)+"\n")

def printPerFile(output,filename,HBcount,dehydronCount,numAtoms,AtomCountInHBDD,numUniqueAtomsInHBDD,numAtomsInDDD,numUniqueAtomsInDDD):
	HBless1SD = .15865 * HBcount
	HBless11SD = .136 * HBcount #That's 1.1 SD's below the mean
	HBless15SD = .067 * HBcount #That's 1.5 SD's below the mean
	output.write(str(filename)+","+\
				 str(HBcount)+","+\
				 str(HBless1SD)+","+\
				 str(HBless11SD)+","+\
				 str(HBless15SD)+","+\
				 str(dehydronCount)+","+\
				 str(numAtoms)+","+\
				 str(AtomCountInHBDD)+","+\
				 str(numAtomsInDDD)+","+\
				 str(numUniqueAtomsInHBDD)+","+\
				 str(numUniqueAtomsInDDD)+"\n")

def countNPC(pdb,CA1,CA2,totalNeighbors):
	CA1name = Atom.getResname(CA1)
	CA2name = Atom.getResname(CA2)
	neighbors = pdb.select('within '+str(HBOND_CA_DISTANCE)+' of CA1 or within '+str(HBOND_CA_DISTANCE)+' of CA2',CA1=CA1, CA2=CA2)
	uniqueAtomCount = neighbors.numAtoms()
	DDDneighbors = []
	for i in Selection.getIndices(neighbors):
		totalNeighbors.append(i)
		DDDneighbors.append(i)

	NPCcount = 0
	for atom in neighbors:
		atomSeqName = Atom.getResname(atom)
		atomName = Atom.getName(atom)
		if atomSeqName in polarity and atomName in polarity[atomSeqName]:
			# if Atom.getResnum(CA1) == 17:
			# 	print atomSeqName, atomName
			NPCcount +=1
	
	return NPCcount, uniqueAtomCount,DDDneighbors

def writeHeaders(perDehydronOut,perHbondOut,perFileOutput,perFileModel,modeledAtomsInDD):
	perDehydronOut.write("filename,ca1ResNum,ca1X,ca1Y,ca1Z,ca2ResNum,ca2X,ca2Y,ca2Z,centerX,centerY,centerZ,ca1resName,ca2resName,NPCcount,uniqueAtomsInDDD\n")
	perHbondOut.write("filename,ca1ResNum,ca1X,ca1Y,ca1Z,ca2ResNum,ca2X,ca2Y,ca2Z,centerX,centerY,centerZ,ca1resName,ca2resName,NPCcount,uniqueAtomsInHBDD,isDehydron\n")
	perFileOutput.write("filename,HBcount,HBless1SD,HBless1.1SD,HBless1.5SD,dehydronCount,totalAtomsInPDB,totalAtomsInHBDD,totalAtomsInDDD,totalUniqueAtomsInHBDD,totalUniqueAtomsInDDD\n")
	perFileModel.write("filename,ca1ResNum,ca1X,ca1Y,ca1Z,ca2ResNum,ca2X,ca2Y,ca2Z,centerX,centerY,centerZ,ca1resName,ca2resName,NPCcount,uniqueAtomsInDDD\n")
	modeledAtomsInDD.write("filename,modeledUniqueAtomsInDDD,actualUniqueAtomsInDDD\n")

################################################################################
# MAIN #
########
inputFile = sys.argv[1]

#analyze directory full of pdb's 
if inputFile[-1:] == "/":
	perDehydronOut = open(inputFile[:-1]+"_perDehydron_output.csv","w")
	perHbondOut = open (inputFile[:-1]+"_perHBond_output.csv","w")
	perFileOutput = open(inputFile[:-1]+"_perFile_output.csv","w")	
	perFileModel = open(inputFile[:-1]+"_perFileModel_output.csv","w")
	modeledAtomsInDD = open(inputFile[:-1]+"_modeledAtomsInDD_output.csv","w")

	writeHeaders(perDehydronOut,perHbondOut,perFileOutput,perFileModel,modeledAtomsInDD)
	
	for filename in os.listdir(inputFile):
		if filename.endswith(".pdb"):
			#print "###",filename,"###"
			analyze(inputFile+"/"+filename,perDehydronOut,perHbondOut,perFileModel,perFileOutput,modeledAtomsInDD)

#analyze single pdb
else:
	perDehydronOut = open(inputFile+"_single_perDehydron_output.csv","w")
	perHbondOut = open (inputFile+"_single_perHBond_output.csv","w")
	perFileOutput = open(inputFile+"_single_perFile_output.csv","w")
	perFileModel = open(inputFile+"_single_perFileModel_output.csv","w")
	modeledAtomsInDD = open(inputFile+"_single_modeledAtomsInDD_output.csv","w")

	writeHeaders(perDehydronOut,perHbondOut,perFileOutput,perFileModel,modeledAtomsInDD)
	
	if inputFile[-4:] == ".pdb":
		''' 
		### If we want to download new pdb files, uncomment this code ###
		slash = inputFile.rfind("/") + 1
		underscore = inputFile.rfind("_")
		if underscore != -1:
			inputFile=inputFile[:(underscore)]
		else:
			inputFile=inputFile[:-4]
		if slash != -1:
			inputFile=inputFile[slash:]
		'''
	print "###",inputFile,"###"
	analyze(inputFile,perDehydronOut,perHbondOut,perFileModel,perFileOutput,modeledAtomsInDD)
perDehydronOut.close()
perHbondOut.close()
perFileModel.close()
perFileOutput.close()
