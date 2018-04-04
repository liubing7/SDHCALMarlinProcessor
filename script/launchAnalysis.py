#!/usr/bin/env python

import os
import sys 

import AnalysisProcessor as ap

import argparse

if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('runNumber' , help='number of the run' , type=str)

	args = parser.parse_args()

	runNumber = str(args.runNumber) 


	#ilcsoft
	os.environ["MARLIN_DLL"] = '/home/liu/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so'

	inputDir = '/home/liu/files/DATA/TRIVENT/SPS_Sep2017'

	inputFile = 'TDHCAL_' + runNumber + '.slcio'
	fileList = [ inputDir + '/' + inputFile ]

	a = ap.Params()
	a.collectionName = 'SDHCAL_HIT'
	a.runNumber = runNumber

	a.outputFileName = runNumber + '.root'

	ap.launch(a , fileList)


  	outputDir = inputDir.replace('TRIVENT' , 'Analysis')

  	os.system('mkdir -p ' + outputDir) 

	os.system('mv ' + a.outputFileName + ' ' + outputDir)
