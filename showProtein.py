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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from prody.utilities import openFile, showFigure
import random

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
mainChainColors = ["blue","green","purple","blue","green","purple","blue","green","purple"]


################################################################################
# Parsing #
###########

def parse(filename,arguments):
	ppdb=parsePDB(filename,model=1)
	if ppdb.numAtoms() > 10000:
		print "abort! too many atoms! ain't nobody got time fo dat"
		return
	plt.ion()

	calpha = ppdb.calpha

	CAneighbors = findNeighbors(calpha,HBOND_CA_DISTANCE)
	if len(CAneighbors) > 4000:
		print "abort, too many neighbors implies too much computation time"
		return
	

	dehydronList = []
	hbList = []
	totalCaptured = []
	for (CA1, CA2, CAdist) in CAneighbors:
		CA1res = Atom.getResnum(CA1)
		if CA1res < 1:
			break
		CA1group = ppdb.select('resnum '+str(CA1res))
		CA1Hydrogens = CA1group.select('name HA or name 1HA or name 2HA or name H')
		CA1Nitrogens = (CA1group.select('name N'))
		CA1Oxygens = (CA1group.select('name O'))
		CA1B = None
		if CA1Oxygens != None:
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
		if CA2Oxygens != None:
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
												NPCcount,captured = countNPC(ppdb,CA1,CA2)
												if NPCcount < 19:
													dehydronList.append((CA1,CA2))
													totalCaptured += captured
												hbList.append((CA1,CA2))
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
												NPCcount,captured = countNPC(ppdb,CA1,CA2)
												if NPCcount < 19:
													dehydronList.append((CA1,CA2))
													totalCaptured += captured
												hbList.append((CA1,CA2))
												break
	sampledDehydronList = modelSampledHBList(len(hbList),hbList,len(dehydronList))
	if "significance" in arguments:
		twentySamples = []
		for i in range(0,20):
			sampleDehydrons = modelSampledHBList(len(hbList),hbList,len(dehydronList))
			twentySamples.append(sampleDehydrons)

	if "ticks" in arguments:
		print "painting ticks"
		plt.hlines([1,1,1,1,1],[1,1.2,1.4,.8,.6],max(ppdb.getResnums()))  # Draw a horizontal line
		plt.xlim(0,max(ppdb.getResnums())+5)
		plt.ylim(0.4,1.5)
		ca1s=[]
		ca2s=[]
		cas=[]
		for (ca1,ca2) in dehydronList:
			ca1Res = ca1.getResnum()
			ca2Res = ca2.getResnum()
			#if abs(ca1Res - ca2Res) < 21:
			cas.append(abs(ca1Res - ca2Res) + min([ca1Res,ca2Res]))
			ca1s.append(ca1Res)
			ca2s.append(ca2Res)
		mca1s = []
		mca2s = []
		mcas = []
		for (ca1,ca2) in sampledDehydronList:
			ca1Res = ca1.getResnum()
			ca2Res = ca2.getResnum()
			if abs(ca1Res - ca2Res) < 21:
				mcas.append(abs(ca1Res - ca2Res) + min([ca1Res,ca2Res]))
			mca1s.append(ca1Res)
			mca2s.append(ca2Res)
		hbondRes =[]
		hbondcenters=[]
		for (ca1,ca2) in hbList:
			ca1Res = ca1.getResnum()
			ca2Res = ca2.getResnum()
			if abs(ca1Res - ca2Res) < 21:
				hbondcenters.append(abs(ca1Res - ca2Res) + min([ca1Res,ca2Res]))
			hbondRes.append(ca1Res)
			hbondRes.append(ca2Res)

		y = np.ones(np.shape(ca1s))
		y1 = np.ones(np.shape(hbondRes))
		y14=y*1.4
		y12= y*1.2
		y8=y*.8
		y6=y*.6
		yHbondCenters = np.ones(np.shape(hbondcenters))
		ycas = np.ones(np.shape(cas))*1.2
		ymcas = np.ones(np.shape(mcas))*.8
		actualIODstring = "observed index of dispersion: "+str((np.std(cas)**2)/np.mean(cas))
		modelIODstring = "modeled index of dispersion: "+str((np.std(mcas)**2)/np.mean(mcas))
		if "significance" in arguments:
			modelSum = 0.0
			for model in twentySamples:
				singleModel = []
				for (ca1,ca2) in model:
					ca1Res = ca1.getResnum()
					ca2Res = ca2.getResnum()
					if abs(ca1Res - ca2Res) < 21:
						singleModel.append(abs(ca1Res - ca2Res) + min([ca1Res,ca2Res]))
				modelSum += (np.std(singleModel)**2 / np.mean(singleModel))
			sigModelIODString = "average of 20 models' IoD: "+str(modelSum / len(twentySamples))
			plt.text(1,.5,sigModelIODString)
		#plt.plot(hbondRes,y1,".",color="black",ms=10)
		#plt.plot(ca1s,y12,".",color="red",ms=10)
		#plt.plot(ca2s,y14,".",color="orange",ms=10)
		#plt.plot(mca1s,y6,".",color="blue",ms=10)
		#plt.plot(mca2s,y8,".",color="green",ms=10)
		plt.text(1,1.4,actualIODstring)
		plt.text(1,.6,modelIODstring)
		plt.plot(hbondcenters,yHbondCenters,".",color="black",ms=10)
		plt.plot(cas,ycas,".",color="red",ms=10)
		plt.plot(mcas,ymcas,".",color="blue",ms=10)
		plt.axis('off')
		plt.show()
		pause(2**31-1)

	else:
		cf = plt.figure()
		show = Axes3D(cf)
		if calpha and ("alphas" in arguments or "mainchain" in arguments):
			for ch in HierView(calpha, chain=True):
				xyz = ch._getCoords()
				chid = ch.getChid()
				show.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2],
						  label=filename,
						  color=mainChainColors.pop(0),
						  lw=1)
		if "sampled" in arguments or "model" in arguments:
			for (ca1,ca2) in sampledDehydronList:
				xyz1 = ca1._getCoords()
				xyz2 = ca2._getCoords()
				show.plot([xyz1[0],xyz2[0]], [xyz1[1],xyz2[1]], [xyz1[2],xyz2[2]],
						     color="g",
						     lw=2)
				if "domains" in arguments or "ddd" in arguments and "small" not in arguments:
					(xs,ys,zs) = drawSphere(xyz1[0],xyz1[1],xyz1[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="g",alpha=.4)
					(xs,ys,zs) = drawSphere(xyz2[0],xyz2[1],xyz2[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="g",alpha=.4)
				if "small" in arguments:
					show.scatter([xyz1[0]],[xyz1[1]],[xyz1[2]],
						c="r", s=1000, alpha=.2)
					show.scatter([xyz2[0]],[xyz2[1]],[xyz2[2]],
						c="r", s=1000, alpha=.2)
		if dehydronList and "dehydrons" in arguments:
			print "painting dehydrons..."
			for (ca1,ca2) in dehydronList:
				xyz1 = ca1._getCoords()
				xyz2 = ca2._getCoords()
				show.plot([xyz1[0],xyz2[0]], [xyz1[1],xyz2[1]], [xyz1[2],xyz2[2]],
						     color="red",
						     lw=2)
				if "domains" in arguments or "ddd" in arguments and "small" not in arguments:
					(xs,ys,zs) = drawSphere(xyz1[0],xyz1[1],xyz1[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="r",alpha=.4)
					(xs,ys,zs) = drawSphere(xyz2[0],xyz2[1],xyz2[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="r",alpha=.4)
				if "small" in arguments:
					show.scatter([xyz1[0]],[xyz1[1]],[xyz1[2]],
						c="r", s=1000, alpha=.2)
					show.scatter([xyz2[0]],[xyz2[1]],[xyz2[2]],
						c="r", s=1000, alpha=.2)
		if hbList and "hbonds" in arguments:
			print "painting hbonds..."
			for (ca1,ca2) in hbList:
				xyz1 = ca1._getCoords()
				xyz2 = ca2._getCoords()
				show.plot([xyz1[0],xyz2[0]], [xyz1[1],xyz2[1]], [xyz1[2],xyz2[2]],
						     color="blue",
						     lw=2)
				if "domains" in arguments or "hbdd" in arguments:
					(xs,ys,zs) = drawSphere(xyz1[0],xyz1[1],xyz1[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="blue",alpha=.3)
					(xs,ys,zs) = drawSphere(xyz2[0],xyz2[1],xyz2[2],HBOND_CA_DISTANCE)
					show.plot_wireframe(xs,ys,zs,color="blue",alpha=.3)
		if "atoms" in arguments:
			print "painting atoms..."
			for atom in ppdb:
				xyz = atom._getCoords()
				show.scatter(xyz[0],xyz[1],xyz[2],color="black",s=.5)
		# Don't use these they're SLOOOOOOWWW
		if "uncaptured" in arguments and "captured" not in arguments:
			print "painting uncaptured atoms..."
			totalCaptured = list(set(totalCaptured))
			for atom in ppdb:
				if atom not in totalCaptured:
					xyz = atom._getCoords()
					show.scatter(xyz[0],xyz[1],xyz[2],color="black",s=.5)
		if "captured" in arguments and "uncaptured" not in arguments:
			print "painting captured atoms..."
			totalCaptured = list(set(totalCaptured))
			for atom in ppdb:
				if atom not in totalCaptured:
					xyz = atom._getCoords()
					show.scatter(xyz[0],xyz[1],xyz[2],color="black",s=.5)
		print "plotting"
		#show.set_aspect('equal','box')
		showFigure()
		cf.show()
		pause(2**31-1)
		print "done"	

def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)

def countNPC(pdb,CA1,CA2):
	CA1name = Atom.getResname(CA1)
	CA2name = Atom.getResname(CA2)
	neighbors = pdb.select('within '+str(HBOND_CA_DISTANCE)+' of CA1 or within '+str(HBOND_CA_DISTANCE)+' of CA2',CA1=CA1, CA2=CA2)

	NPCcount = 0
	captured = []
	for atom in neighbors:
		atomSeqName = Atom.getResname(atom)
		atomName = Atom.getName(atom)
		if atomSeqName in polarity and atomName in polarity[atomSeqName]:
			NPCcount +=1
		captured.append(atom)
	
	return NPCcount,captured

def modelSampledHBList(HBcount,HBList,dehydronCount):
	sampledHBList = random.sample(HBList, dehydronCount)
	modeledDehydronList = []
	for (CA1,CA2) in sampledHBList:
		modeledDehydronList.append((CA1,CA2))
	return modeledDehydronList

################################################################################
# MAIN #
########
inputFile = sys.argv[1]
parse(inputFile, sys.argv)
#analyze single pdb



