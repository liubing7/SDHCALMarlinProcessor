#!/usr/bin/env python

import os

class Params :
	def __init__(self) :
		self.collectionName = 'HCALEndcap'
		self.cosThetaLimit = 0.9
		self.thresholds = '1.0 2.0 3.0'
		self.outputFileName = 'map.root'

def launch(a , files) :

	fileList = ''
	for name in files :
		fileList += name + ' '

	pid = os.getpid()

	xmlFileName = str(pid) + '.xml'

	xml = '''<marlin>
 <execute>
  <processor name="EfficiencyProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles"> ''' + fileList + ''' </parameter>
  <!-- limit the number of processed records (run+evt): -->
  <!--parameter name="MaxRecordNumber" value="100000"/-->
  <!--parameter name="SkipNEvents" value="18000" /-->
  <parameter name="SupressCheck" value="false" />
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE  </parameter>
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="EfficiencyProcessor" type="EfficiencyProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> ''' + a.collectionName + '''</parameter>
  <!--Name of the root output file-->
  <parameter name="RootFileName" type="string" > ''' + a.outputFileName + ''' </parameter>
  <parameter name="InteractionFinder::PrintDebug" type="bool"> false </parameter>
  <parameter name="Tracking::CosThetaLimit" type="float"> ''' + str(a.cosThetaLimit) + ''' </parameter>
  <parameter name="Tracking::PrintDebug" type="bool"> false </parameter>

  <parameter name="Thresholds" type="FloatVec"> ''' + a.thresholds + ''' </parameter>

 </processor>

</marlin>'''

	xmlFile = open(xmlFileName , 'w')
	xmlFile.write(xml)
	xmlFile.close()


	os.system('Marlin ' + xmlFileName)
	os.system('rm ' + xmlFileName)

