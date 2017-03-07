#!/usr/bin/env python

import os

class Params :
	def __init__(self) :
		self.collectionName = 'HCALEndcap'
		self.outputFileName = 'analysis.root'
		self.runNumber = 0

def launch(a , files) :

	fileList = ''
	for name in files :
		fileList += name + ' '

	pid = os.getpid()

	xmlFileName = str(pid) + '.xml'

	xml = '''<marlin>
 <execute>
  <processor name="AnalysisProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">''' + a.filelist + '''</parameter>
  <!-- limit the number of processed records (run+evt): -->
  <!--parameter name="MaxRecordNumber" value="200"/-->
  <!--parameter name="SkipNEvents" value="18000" /-->
  <parameter name="SupressCheck" value="false" />
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter> 
 </global>

 <processor name="AnalysisProcessor" type="AnalysisProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit">''' + a.collectionName + '''</parameter>
  <!--Name of the root output file-->
  <parameter name="RootFileName" type="string" >''' + a.outputFileName + '''</parameter>
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

