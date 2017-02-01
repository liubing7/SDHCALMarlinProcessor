#!/bin/bash

source /home/garillot/ilcsoft/v01-17-08/init_ilcsoft.sh
export MARLIN_DLL=/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so

#run=$1
particle=$1
energy=$2
#model=$3
model=FTFP_BERT_HP

echo " .... Marlin process SDHCAL analysis processor"

FILEPATH="/home/garillot/files/grid/CalorimeterHit/Geant4.9.6/${model}"
FILELIST=`ls ${FILEPATH} | grep -F single_${particle}_${energy}GeV`

#FILEPATH="/home/garillot/files/DATA/TRIVENT/SPS_Oct2015/163"
#FILELIST=`ls ${FILEPATH} | grep -F TDHCAL_${run}`

#FILEPATH="/home/garillot/files/grid/CalorimeterHit/Geant4.9.6/FTF_BIC"
#FILELIST=`ls ${FILEPATH} | grep -F single_pi-_80GeV.slcio`

OUTPUT="single_${particle}_${energy}GeV.root"

#OUTPUT="${run}.root"

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
  <!--parameter name="MaxRecordNumber" value="20000"/-->  
  <!--parameter name="SkipNEvents" value="18000" /-->  
  <parameter name="SupressCheck" value="false" />  
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->  
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="AnalysisProcessor" type="AnalysisProcessor">
  <!--Name of the CalorimeterHit collection-->
  <!--parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> SDHCAL_HIT </parameter-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HCALEndcap </parameter>
  <!--Name of the root output file-->
  <parameter name="RootFileName" type="string" > ${OUTPUT} </parameter>
  <parameter name="InteractionFinder::PrintDebug" type="bool"> false </parameter>
  <parameter name="Tracking::PrintDebug" type="bool"> false </parameter>
 </processor>

</marlin>

EOF

Marlin LCIO.xml
rm LCIO.xml

#mv ${OUTPUT} /home/garillot/files/local/DATA/RootFile/${OUTPUT}
mv ${OUTPUT} /home/garillot/files/local/RootFile/newAnalysis/Geant4.9.6/${model}/${OUTPUT}


