#!/usr/bin/env python

import os
import sys

class Params :
	def __init__(self) :
		self.collectionName = 'HCALEndcap'
		self.outputFileName = 'analysis.root'
		self.runNumber = 0
		self.maxRecordNumber = 0

def launch(a , files) :

	fileList = ''
	for name in files :
		fileList += name + ' '

	pid = os.getpid()

	xmlFileName = str(pid) + '.xml'
	tempOutputFile = str(pid) + '.root'

	xml = '''<marlin>
 <execute>
  <processor name="AnalysisProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">''' + fileList + '''</parameter>
  <!-- limit the number of processed records (run+evt): -->
  <parameter name="MaxRecordNumber" value="''' + str(a.maxRecordNumber) + '''"/>
  <!--parameter name="SkipNEvents" value="18000" /-->
  <parameter name="SupressCheck" value="false" />
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter> 
 </global>

 <processor name="AnalysisProcessor" type="AnalysisProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit">''' + a.collectionName + '''</parameter>
  <!--Name of the root output file-->

  <parameter name="RootFileName" type="string" >''' + tempOutputFile + '''</parameter>
  <parameter name="InteractionFinder::PrintDebug" type="bool"> false </parameter>
  <parameter name="Tracking::PrintDebug" type="bool"> false </parameter>

  <parameter name="nRun" type="int">''' + str(a.runNumber) + '''</parameter>

 </processor>

</marlin>'''

	xmlFile = open(xmlFileName , 'w')
	xmlFile.write(xml)
	xmlFile.close()

	os.system('Marlin ' + xmlFileName)
	os.system('rm ' + xmlFileName)
	os.system('mv ' + tempOutputFile + ' ' + a.outputFileName)



if __name__ == '__main__' :

	if len(sys.argv) < 2 :
		sys.exit('Error : too few arguments')


	fileList = [ sys.argv[1] ]

	print 'File List :'
	print fileList

	os.environ["MARLIN"] = '/home/garillot/ilcsoft/v01-19-03/Marlin/v01-12'
	os.environ["PATH"] = '/home/garillot/ilcsoft/v01-19-03/Marlin/v01-12/bin' + ':' + os.environ["PATH"]
	os.environ["MARLIN_DLL"] = '/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so'

	a = Params()

	a.outputFileName = sys.argv[1].replace('.slcio' , '.root')

	launch(a , fileList)

