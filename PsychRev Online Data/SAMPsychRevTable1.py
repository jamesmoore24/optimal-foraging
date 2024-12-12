#!/usr/bin/env python
# encoding: utf-8
"""
Created by Thomas Hills on 2008-10-31.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.

This file creates the data for Table 1

run in python

"""

import sys
import os
import time
import datetime
import copy
from pylab import *
from numpy import *
from scipy.optimize import fmin

startat = 0

def get_labels():
	ofile = open(labelfile, 'rU')
	for line in ofile:
		words = line.split()
		for word in words:
			wordlist.append(word)
	ofile.close()
	
def get_matrix():
	ofile = open(matfile1, 'rU')		### for optimal for.  matfile1 is the fluidcat matrix
	j = 0
	for line in ofile:
		sims = line.split()
		i = 0
		for sim in sims:
			catmat[i,j] = sim
			i += 1
		j += 1
	ofile.close()

def get_matrix2():						## for optimal for. matfile2 is the ancosim matrix
	ofile = open(matfile2, 'rU')
	j = 0
	for line in ofile:
		sims = line.split()
		i = 0
		for sim in sims:
			cosmat[i,j] = sim
			i += 1
		j += 1
	ofile.close()
		
def get_globalcue():
	ofile = open(globalcuefile, 'rU')
	for line in ofile:
		freqms = line.split()
		globallabs.append(freqms[0])
		globalcue.append(float(freqms[1]))
	ofile.close()
		
def fill_sidentrylist():
	ofile = open(entryfile, 'rU')
	for line in ofile:
		dt = line.split()
		sidlist.append(dt[0])
		entrylist.append(dt[1])
	ofile.close()

def get_groups():
	ofile = open(groupfile, 'rU')
	for line in ofile:
		dt = line.split()
		groupwords.append(dt[0])
		groupgroups.append(int(dt[1]))
	ofile.close()
		
def fill_sidswitchlist():
	ofile = open(switchfile, 'rU')
	for line in ofile:
		dt = line.split()
		sidlist.append(dt[0])
		switchlist.append(int(dt[1]))
	ofile.close()
		
def unique(inlist, keepstr=True):
  	typ = type(inlist)
  	if not typ == list:
		inlist = list(inlist)
  	i = 0
  	while i < len(inlist):
		try:
			del inlist[inlist.index(inlist[i], i + 1)]
		except:
			i += 1
	if not typ in (str, unicode):
		inlist = typ(inlist)
	else:
		if keepstr:
			inlist = ''.join(inlist)
	return inlist

def samfit(beta, globn, globd, simn, simd):
	ct = 0
	for colent in range(startat, len(globn)):	
		if colent == 0:
			numrat = pow(globn[colent],beta[0])
			denrat = sum(pow(globd[colent],beta[0]))
		else:
			numrat = pow(globn[colent],beta[0])*pow(simn[colent],beta[1])
			denrat = sum(pow(globd[colent],beta[0])*pow(simd[colent],beta[1]))		
		ct += -log(numrat/denrat)
	return	ct

def samfitglobal(beta, globn, globd):		## keeper
	ct = 0
	for colent in range(startat, len(globn)):	
		numrat = pow(globn[colent],beta[0])
		denrat = sum(pow(globd[colent],beta[0]))
		#print "start" + str(colent)
		#print numrat
		#print denrat		
		ct += -log(numrat/denrat)
	return	ct
	
def samfitlocal(beta, simn, simd):			## keeper
	ct = 0
	for colent in range(1, len(simn)):	
		numrat = pow(simn[colent],beta[0])
		denrat = sum(pow(simd[colent],beta[0]))	
		ct += -log(numrat/denrat)
	return	ct

def getfluidswitches(beta, globn, globd, productions):
	ct = 0
	subgroups = set(groupgroups)	# unique groups 
	subw = groupwords[:]			# words in groups
	groupn = array(groupgroups)		# group for each word in subw
	words = subwordlist[:]			# overall wordlist
	#basearray = zeros(len(words))	# overall in-out group array that goes with words
	switches = zeros(len(globn))
	for colent in range(0, len(globn)):	
		wordgroups = []		# this words groups

		for j in range(0,len(groupn)):
			if subw[j] == productions[colent]:
				wordgroups.append(groupn[j])	# assign groups from subw
		#print "entry " + str(colent)
		##print productions[colent]
		#print wordgroups
		subgroups = set.intersection(set(subgroups), set(wordgroups))	# the intersection of old with new groups
		#print subgroups
		if colent == 0 or len(subgroups) == 0:	# first entry or switch
			numrat = pow(globn[colent],beta[0])
			denrat = sum(pow(globd[colent],beta[0]))
			subgroups = wordgroups
			#print switches
			#print colent
			switches[colent] = 1
		else:
			numrat = pow(globn[colent],beta[0])
			basearray = zeros(len(words))
			for wg in range(0,len(groupn)):	# go through groups 
				if groupn[wg] in subgroups:	# if group for for word in subw list is among groups that are appropriate to this patch
					if subw[wg] in words:	# if the corresponding word is in the list of producible words
						basearray[words.index(subw[wg])] = 1	# assign word as in group
			denrat = sum(pow(globd[colent],beta[0])*power(basearray, beta[1]))
			subgroups = wordgroups

		ct += -log(numrat/denrat)
	#print "switches groups: " + str(switches)
	return	switches
			
						
def samswitchfit(beta, globn, globd, simn, simd, switch):
	ct = 0
	for colent in range(startat, len(globn)):	
		if colent == 0 or switch[colent] == 1:
			numrat = pow(globn[colent],beta[0])
			denrat = sum(pow(globd[colent],beta[0]))
		else:
			numrat = pow(globn[colent],beta[0])*pow(simn[colent],beta[1])
			denrat = sum(pow(globd[colent],beta[0])*pow(simd[colent],beta[1]))		
		ct += -log(numrat/denrat)
	return	ct

def samfindswitchfit(beta, globn, globd, simn, simd, sjmn):  ## similarity drop
	ct = 0
	for colent in range(startat, len(globn)):	
		if colent == 0 :
			numrat = pow(globn[colent],beta[0])
			denrat = sum(pow(globd[colent],beta[0]))
		elif colent > 0 and colent < (len(globn)-1) and sjmn[colent+1] > sjmn[colent] and sjmn[colent-1] > sjmn[colent]:
			numrat = pow(globn[colent],beta[0])
			denrat = sum(pow(globd[colent],beta[0]))			
		else:
			numrat = pow(globn[colent],beta[0])*pow(simn[colent],beta[1])
			denrat = sum(pow(globd[colent],beta[0])*pow(simd[colent],beta[1]))		
		ct += -log(numrat/denrat)
	return	ct

def writeData(datab):
	print "writing data to file"
	wfile1 = open(outfilename, 'wb')
	wfile1.write('e1\ta1\tl1\tr1\ta2\tl2\tr2\ta3\tb3\tl3\tr3\ta4\tb4\tl4\tr4\ta5\tb5\tl5\tr5\n')
	for row in datab:
		for entry in row:
			wfile1.write('%12.5f\t' % (entry))
		wfile1.write('\n')
	wfile1.close()

def writeSwitches(datab):
	print "writing data to file"
	wfile1 = open("Switchfile", 'wb')
	for row in datab:
		wfile1.write('%12.5f\n' % (row))
	wfile1.close()

				
if __name__=="__main__":
	
	labelfile = "datamatrixlabels.txt"	## labels for matrix
	matfile2 = "dataancosim.txt"	## similarity matrix
	matfile1 = "datafluidcat.txt"		## similarity matrix
	entryfile = "dataentrylist.txt"		## entry list (only has id and entry (entries should match labels--if not, error))
	globalcuefile = "datafreqlistlog.txt"	## global cue
	switchfile = "dataswitches5177.txt"	## switch cue
	
	group = 1
	if group:
		groupfile = "datagroups.txt"
		groupwords = []
		groupgroups = []
		get_groups()
		grouplist = []
		
	switchlist = []
	sidlist = []
	fill_sidswitchlist()
	
	allswitches = []
	
	t = datetime.datetime.now()
	outfilename = "output" + str(int(time.mktime(t.timetuple()))) + ".txt"	## output file

	wordlist = []
	get_labels()	## assign labels
	cosmat = zeros((len(wordlist), len(wordlist)))
	catmat = zeros((len(wordlist), len(wordlist)))
	get_matrix()	## assign sim matrix
	get_matrix2()
	
	for i in range(0,len(cosmat)):		## makes matrix > 0, because log of 0 is -inf
		for j in range(0,len(cosmat)):
			if cosmat[i][j] <= 0:
				cosmat[i][j] = 0.001
			catmat[i][j] = catmat[i][j] + 0.1


					
	globalcue = []
	globallabs = []
	get_globalcue()
	globalcue = array(globalcue)
	for i in range(0, len(globalcue)):
		if globalcue[i] <= 0:
			globalcue[i] = 0
		globalcue[i] = globalcue[i] + 0.001	## make it greater than one
		
	sidlist = []
	entrylist = []
	fill_sidentrylist()		## gets sids and entries
	unisidlist = sidlist[:]	# this is the long list of ids cut down to uniques
	unique(unisidlist)	
	

	
	betafits = zeros([2,len(unisidlist)], float)
	loglfits = []
	randomlfits = []

	databed = []

	#for vid in range(1, 2):#			TO 	run through each participant
	for vid in range(0, len(unisidlist)):
		id = unisidlist[vid]		# take first sid
		print "id is " + str(id)	
		
		usedlist = []	## this will keep track of produced wods
		sublist = []	## this is the words the subject produces
		subswitchlist = []
		subglob = []	## this is the global cue array
		
		subcosmat = cosmat.copy()	## local cue
		subcatmat = catmat.copy()	## second local cue
		#subglob = array(log(globalcue))	## global cue
		subglob = array(globalcue)	## global cue
		subwordlist = wordlist[:]	## local labels list

		print "matrix size" + str(len(subcosmat[1,:])) + ", " + str(len(subcosmat[:,1]))
		
		i = 0
		for sid in sidlist:
			if id == sid:
				sublist.append(entrylist[i])	### fills up sublist with entries for sid
				subswitchlist.append(switchlist[i])
			i += 1
		ct = 0
		
		globallisti = []		## initialize parameters
		globallistk = []	# global
		simlist1i = []		# local
		simlist1k = []		# local	- 1
		simlist2i = []		# local
		simlist2k = []		# local - 2
		simlist12i = []		# second local
		simlist12k = []		# local	- 1
		simlist22i = []		# second local
		simlist22k = []		# local	- 1
		productions = []
		rowsums = []
		
		## count through items in sublist (entries for his subject)
		## THIS INITIALIZES THE IMPORTANT DATA FOR THIS PARTICIPANT
		for ent in range(0,len(sublist)):
		##  get numerator: things you need sublist, wordlist, subcosmat, freqlist = globallist
		##  sim to animal
		##  freq
			if sublist[ent] not in usedlist:	# only use words that haven't been said before by this sid
				# global sim
				globallisti.append(float(subglob[globallabs.index(sublist[ent])]))
				globallistk.append(array(subglob))
				
				## one back sim
				if ent > 0:
					simlist1i.append(float(subcosmat[wordlist.index(sublist[ent-1]),subwordlist.index(sublist[ent]) ]))
					simlist1k.append(array(subcosmat[wordlist.index(sublist[ent-1]),:]))
					simlist12i.append(float(subcatmat[wordlist.index(sublist[ent-1]),subwordlist.index(sublist[ent]) ]))
					simlist12k.append(array(subcatmat[wordlist.index(sublist[ent-1]),:]))
				else:
					simlist1i.append(0)	
					simlist1k.append(array(subcosmat[wordlist.index(sublist[ent]),:]))
					simlist12i.append(0)	
					simlist12k.append(array(subcatmat[wordlist.index(sublist[ent]),:]))
				## two back sim
				if ent > 1:
					simlist2i.append(array(subcosmat[wordlist.index(sublist[ent-2]),subwordlist.index(sublist[ent]) ]))
					simlist2k.append(array(subcosmat[wordlist.index(sublist[ent-2]),:]))
					simlist22i.append(array(subcatmat[wordlist.index(sublist[ent-2]),subwordlist.index(sublist[ent]) ]))
					simlist22k.append(array(subcatmat[wordlist.index(sublist[ent-2]),:]))
				else:
					simlist2i.append(0)	
					simlist2k.append(array(subcosmat[wordlist.index(sublist[ent]),:]))
					simlist22i.append(0)	
					simlist22k.append(array(subcatmat[wordlist.index(sublist[ent]),:]))
			
					
					
			#  add word to used list, and set used words to zero in the global and local arrays
				usedlist.append(sublist[ent])
				subglob[globallabs.index(sublist[ent])] = 0.00000001
				subcosmat[:,subwordlist.index(sublist[ent]) ] = 0.00000001
				subcatmat[:,subwordlist.index(sublist[ent]) ] = 0
				productions.append(sublist[ent])
				rowsums.append(mean(simlist1k[-1]))
			
		## set inputs to arrays				::  this doesn't appear to be needed
		globallisti = array(globallisti)
		simlist1i = array(simlist1i)
		simlist2i = array(simlist2i)
		simlist12i = array(simlist12i)
		simlist22i = array(simlist22i)
		
		v1 = rand()
		v2 = rand()
		v3 = rand()
		v4 = rand()
		v5 = rand()
		datasheet = []
		
		datasheet.append(len(globallisti))
		doit = 1

		##  these are the fluid switches computed within this program
		## use them, they are correct
		sidfswitchlist = getfluidswitches([0,0],globallisti, globallistk, productions)  # captures list of switches
		
		
		if doit:		## this can be for Tabl1 if Frequency or Animal is the global cue
			# SAM 1 back
			v = fmin(samfitglobal, [v1], args=(globallisti, globallistk),ftol = 0.001)	
			datasheet.append(float(v[0]))
			## get error
			fitt = samfitglobal([float(v[0])],globallisti, globallistk)
			datasheet.append(float(fitt))
			## get random
			fitt = samfitglobal([0],globallisti, globallistk)
			datasheet.append(float(fitt))
			
		if doit:		#  for Table 1, simlist1i is the item-1 in BEAGLE
			# SAM 1 back
			v = fmin(samfitlocal, [v1], args=(simlist1i, simlist1k),ftol = 0.001)	
			datasheet.append(float(v[0]))
			print "got here"
			## get error
			fitt = samfitlocal([float(v[0])],simlist1i, simlist1k)
			datasheet.append(float(fitt))
			## get random
			fitt = samfitlocal([0],simlist1i, simlist1k)
			datasheet.append(float(fitt))
			
		if doit:
			# SAM local _ global				########################  this is the local + global static , table 1
			v = fmin(samfit, [v1,v2], args=(globallisti, globallistk, simlist1i, simlist1k),ftol = 0.001)	
			datasheet.append(float(v[0]))
			datasheet.append(float(v[1]))
			beta0 = float(v[0])
			beta1 = float(v[1])
			## get error
			fitt = samfit([float(v[0]),float(v[1])],globallisti, globallistk, simlist1i, simlist1k)
			datasheet.append(float(fitt))
			bestfitt = fitt
			## get random
			fitt = samfit([0,0],globallisti, globallistk, simlist1i, simlist1k)
			datasheet.append(float(fitt))
			
		if doit:		## Table 2
			# SAM cat switch glob + loc
			v = fmin(samswitchfit, [v1,v2], args=(globallisti, globallistk, simlist1i, simlist1k, sidfswitchlist),ftol = 0.001)	
			datasheet.append(float(v[0]))
			datasheet.append(float(v[1]))
			## get error
			fitt = samswitchfit([float(v[0]),float(v[1])],globallisti, globallistk, simlist1i, simlist1k, sidfswitchlist)
			datasheet.append(float(fitt))
			## get random
			fitt = samswitchfit([0,0],globallisti, globallistk, simlist1i, simlist1k, sidfswitchlist)
			datasheet.append(float(fitt))
						
		if doit:		## similarity drop, table 3
			# SAM find ass switch -- 1st local cue
			v = fmin(samfindswitchfit, [v1,v2,v3], args=(globallisti, globallistk, simlist1i, simlist1k, simlist1i),ftol = 0.001)	
			datasheet.append(float(v[0]))
			datasheet.append(float(v[1]))
			## get error
			fitt = samfindswitchfit([float(v[0]),float(v[1])],globallisti, globallistk, simlist1i, simlist1k, simlist1i)
			datasheet.append(float(fitt))
			## get random
			fitt = samfindswitchfit([0,0],globallisti, globallistk, simlist1i, simlist1k, simlist1i)
			datasheet.append(float(fitt))	
							
		databed.append(array(datasheet))

	writeData(databed)
	writeSwitches(allswitches)

	print "done"
	
