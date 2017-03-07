#!/bin/bash

source /home/garillot/ilcsoft/v01-17-08/init_ilcsoft.sh
export MARLIN_DLL=/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so

run=$1
#energy=$1

#particle=$1
#energy=$2
#model=$3
#model=FTFP_BERT_HP

#energy=$1

echo " .... Marlin process SDHCAL analysis processor"

#FILEPATH="/home/garillot/files/grid/CalorimeterHit/Geant4.9.6/${model}"
#FILELIST=`ls ${FILEPATH} | grep -F single_${particle}_${energy}GeV`

#FILEPATH="/home/garillot/files/DATA/TRIVENT/electrons"
FILEPATH="/home/garillot/files/DATA/TRIVENT/SPS_Oct2015/163"
FILELIST=`ls ${FILEPATH} | grep -F DHCAL_${run}`
COLLECTION="SDHCAL_HIT"

#FILEPATH="/home/garillot/files/grid/CalorimeterHit/Geant4.9.6/FTF_BIC"
#FILELIST=`ls ${FILEPATH} | grep -F single_pi-_80GeV.slcio`

#FILEPATH="/home/garillot/files/local/CalorimeterHit/Geant4.10.01/QGSP_BERT_HP"
#FILELIST=`ls ${FILEPATH} | grep -F single_e-_${energy}GeV.slcio`

#FILEPATH="/home/garillot/files/muonsForBing/CalorimeterHit"
#FILELIST=`ls ${FILEPATH} | grep -F single_mu-_${energy}GeV.slcio`

#COLLECTION="HCALEndcap"



#OUTPUT="single_${particle}_${energy}GeV.root"

#OUTPUT="single_e-_${energy}GeV.root"

#OUTPUT="${energy}Test.root"
OUTPUT="${run}.root"

echo $FILELIST

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="AnalysisProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   ${FILEPATH}/${FILELIST}
  </parameter>
  <!-- limit the number of processed records (run+evt): -->
  <!--parameter name="MaxRecordNumber" value="200"/-->
  <!--parameter name="SkipNEvents" value="18000" /-->
  <parameter name="SupressCheck" value="false" />
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="AnalysisProcessor" type="AnalysisProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> ${COLLECTION} </parameter>
  <!--Name of the root output file-->
  <parameter name="RootFileName" type="string" > ${OUTPUT} </parameter>
  <parameter name="InteractionFinder::PrintDebug" type="bool"> false </parameter>
  <parameter name="Tracking::PrintDebug" type="bool"> false </parameter>
  
  <parameter name="nRun" type="int"> ${run} </parameter>
  
 </processor>

</marlin>

EOF

Marlin LCIO.xml
rm LCIO.xml

#mv ${OUTPUT} /home/garillot/files/local/RootFile/Geant4.10.01/QGSP_BERT_HP
#mv ${OUTPUT} /home/garillot/files/local/DATA/RootFile/electrons
#mv ${OUTPUT} /home/garillot/files/local/RootFile/newAnalysis/Geant4.9.6/${model}/${OUTPUT}


