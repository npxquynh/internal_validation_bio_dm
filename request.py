#!/usr/bin/python


import os
import helper
import parse_xml as px
import parse_lgn as pl
import parse_expanded_network as pen
import write_expansion as we
from internal_validation import *
from internal_validation_2 import *
from subprocess import call
from sys import argv
from math import ceil,pow
from random import shuffle


#DEFAULT_OUTPUT_FILENAME = 'iter-results.xml'
DEFAULT_OUTPUT_FILENAME = argv[5] + '.xml' 
DEFAULT_OUTPUT_FOLDER_NAME = 'results'
DEFAULT_TILE_FILENAME = 'tile'

def getLgnProbes(lgnPath,completePath):
	'''Given the path of a file with the LGN and the path of the complete data file, 
	returns the list of probe indices and the total number of probes in the complete file.'''
	lgn = []

	with open(lgnPath) as f :
		for row in f.readlines()[1:] :
			lgn += row.replace('\n','').split(',')

	lgnNames = set(lgn)	# filters repetitions on probe names
	lgnIndices = []
	with open(completePath) as f :
		rowNumber = 1
		for row in f.readlines()[1:] :
			probeId = row.split(',')[0]
			if probeId in lgnNames :	# the row correspond to a LGN node
				lgnIndices.append(rowNumber)
			rowNumber += 1

	return (lgnIndices,rowNumber)

def filterEdges(edgeList,nodeList) :
	'''Filters the edgeList, returning those which are connecting at least a node of nodeList.'''
	return [e for e in edgeList if e[0] in nodeList or e[1] in nodeList]

def subsetToString(subList):
	s = str(subList[0])	# subList shouldn't be empty
	for x in subList[1:] :
		s += " " + str(x)
	return s

def edge2str(t) :
	'''Just converts a pair (edge) to a string in order to be used as a index.'''
	return str(edge[0]) + ',' + str(edge[1])

def pcim(execPath,lgnPath,completePath,processId,alpha=0.05,nIterations=100,subsetSize=1000) :
	
	# gets the LGN indices from the LGN network file, also returns the total number of probes from the 'complete' file
	lgnProbes,nProbes = getLgnProbes(lgnPath,completePath)
	allProbes = range(1,nProbes+1)

	# nodeCount = dict([(x,0) for x in allProbes]) # initializes the node counter
	# edgeCount = {}

	# creates a working directory accordingly to the processId
	workingDirectory = "./" + str(processId) +'/'
	call(["mkdir",workingDirectory])
	outputFolder = workingDirectory + DEFAULT_OUTPUT_FOLDER_NAME + '/'
	call(["mkdir",outputFolder])
	tilePath = workingDirectory + DEFAULT_TILE_FILENAME

	for i in range(nIterations) :

		# indices of non-LGN nodes, shuffled in order to take different subset on each iteration
		otherProbes = list(set(range(1,nProbes+1)) - set(lgnProbes))
		shuffle(otherProbes)
		
		tileFile = open(tilePath,'w')
		probeBag = list(otherProbes) # probes not in the LGN to be extracted
	
		while probeBag != [] :
			subset,probeBag = probeBag[:subsetSize],probeBag[subsetSize:] # takes a subset of non-LGN indices
			subset = subset + lgnProbes # adds the LGN to the subset

			# if the subset size is smaller, just add some random nodes outside the Local Gene Network
			if len(subset) < subsetSize : 
				shuffle(otherProbes)
				subset = subset + list(set(otherProbes).difference(set(subset)))[:(len(lgnProbes)+subsetSize-len(subset))]

			shuffle(subset)

			#for x in subset :	# increase the node count for each element of the considered subset
			#	nodeCount[x] += 1
			
			row = subsetToString(subset)
			if probeBag != [] :
				row += '\n'
			tileFile.write(row)
		
		tileFile.close()

		# execute the PC algorithm
		call([execPath,completePath,tilePath,outputFolder + DEFAULT_OUTPUT_FILENAME + '.' + str(i), str(alpha)])

def postproc(pcResultFolder, pcimFolder, lgnFilePath) :

	code = DEFAULT_OUTPUT_FILENAME
	
	genes_in_lgn, edges_in_lgn = pl.read_lgn(lgnFilePath)

	list_of_expanded_network_file = pen.merge_different_pc_run(code, pcResultFolder)
	genes_in_tiles, blocks = pen.read_expanded_network(list_of_expanded_network_file)

	# 1st internal validation
	IV = InternalValidation(genes_in_tiles, blocks, genes_in_lgn, edges_in_lgn)

	expansion_list = IV.expansion_list()

	expansion_filepath = helper.generate_filepath(pcimFolder, code)
	we.write_expansion_list(expansion_filepath, expansion_list)

	# 2nd internal validation
#	IV2 = InternalValidationRls(genes_in_tiles, blocks, genes_in_lgn, edges_in_lgn)

#	expansion_list_2 = IV2.expansion_list()
#	expansion_filepath = helper.generate_filepath(pcimFolder, code)
#	we.write_expansion_list(expansion_filepath, expansion_list_2)





if len(argv) != 9 :
	print "Usage: python request.py [PARAMS]"
	print "PARAMS:"
	print "- C++ pcim executable path"
	print "- the LGN genes file path"
	print "- the LGN rls file path"
	print "- the obs gene file path"
	print "- the name of a folder"
	print "- alpha"
	print "- number of iterations on the entire genome"
	print "- size of the subset (lgn included)"
else :
	pcim(argv[1], argv[2], argv[4], argv[5], float(argv[6]), int(argv[7]), int(argv[8]))
	#postproc(argv[5] + '/' + DEFAULT_OUTPUT_FOLDER_NAME, argv[5], argv[3])
	#postproc(argv[5] + '/' + DEFAULT_OUTPUT_FOLDER_NAME, '/home/manfredil', argv[3])
	postproc(argv[5] + '/' + DEFAULT_OUTPUT_FOLDER_NAME, '.', argv[3])
	call(["rm","-rf",argv[5] ])
