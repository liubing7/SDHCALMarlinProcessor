#!/usr/bin/env python

import gfal2

import os
import sys
from os import path
import time

sys.path.insert(0 , '/gridgroup/ilc/garillot/SDHCALMarlinProcessor/script')

import AnalysisProcessor

def download(input , output) :

	# Instantiate gfal2
	ctx = gfal2.creat_context()

	# Set transfer parameters
	params = ctx.transfer_parameters()
	params.overwrite = True
	params.timeout = 300

	dlCounter = 0
	isOK = False

	print 'Try to Download ' + input
	while not isOK and dlCounter < 50 :
		source = input
		destination = output
		try :
			r = ctx.filecopy(params, source, destination)
			isOK = True
		except Exception, e :
			print "Download failed : %s" % str(e)
			isOK = False
			time.sleep(20)

		dlCounter = dlCounter + 1

	if isOK :
		print 'Download succeeded !'
	else :
		print 'Download failed !'
		sys.exit(1)


def upload(input , output) :

	# Instantiate gfal2
	ctx = gfal2.creat_context()

	# Set transfer parameters
	params = ctx.transfer_parameters()
	params.overwrite = True
	params.timeout = 300

	dlCounter = 0
	isOK = False

	print 'Try to Upload ' + input
	while not isOK and dlCounter < 50 :
		source = input
		destination = output
		try:
			r = ctx.filecopy(params, source, destination)
			isOK = True
		except Exception, e:
			print "Upload failed : %s" % str(e)
			isOK = False
			time.sleep(20)

		dlCounter = dlCounter + 1

	if isOK :
		print 'Upload succeeded !'
	else :
		print 'Upload failed !'
		sys.exit(1)



if __name__ == '__main__' :

	if len(sys.argv) < 5 :
		sys.exit('Error : too few arguments')

	os.environ["LFC_HOST"] = 'grid-lfc.desy.de'

	particle = sys.argv[1]
	energy = sys.argv[2]
	model = sys.argv[3]
	version = sys.argv[4]

	versionPath = 'Geant4.' + version

	dir = 'root://lyogrid06.in2p3.fr/dpm/in2p3.fr/home/calice/garillot/CalorimeterHit/' + versionPath + '/' + model

	print ('Searching files in ' + dir)

	ctxt = gfal2.creat_context()
	dirp = ctxt.opendir(dir)

	#try to access dir
	counter = 0
	while not dirp and counter < 50 :
		time.sleep(5)
		dirp = ctxt.opendir(dir)
		counter = counter + 1
	if not dirp :
		print dir + ' not ok or not accessible'
		sys.exit(1)


	#list files
	fileList = []

	while(True) :
		(dirent, stat) = dirp.readpp()
		if (dirent is None) :
			break ;

		if 'single_' + particle + '_' + energy + 'GeV' in dirent.d_name :
			fileList.append(dirent.d_name)

	print 'File List :'
	print fileList

	for file in fileList :
		download(dir + '/' + str(file) , 'file:' + str(file))


	os.environ["MARLIN"] = '/gridgroup/gridsoft/ipnls/ilc/v01-17-09/Marlin/v01-08'
	os.environ["PATH"] = '/gridgroup/gridsoft/ipnls/ilc/v01-17-09/Marlin/v01-08/bin' + ':' + os.environ["PATH"]
	os.environ["MARLIN_DLL"] = '/gridgroup/ilc/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so'


	a = AnalysisProcessor.Params()

	a.outputFileName = 'single_' + particle + '_' + energy + 'GeV.root'

	AnalysisProcessor.launch(a , fileList)

  	outputDir = dir.replace('CalorimeterHit' , 'Analysis')

  	os.system('gfal-mkdir -p ' + outputDir) 

	upload('file:' + a.outputFileName , outputDir + '/' + a.outputFileName)

	os.system('rm ' + a.outputFileName)

	for file in fileList :
		os.system('rm ' + file)













