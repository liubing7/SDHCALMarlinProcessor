#!/bin/bash

source /home/garillot/ilcsoft/v01-17-08/init_ilcsoft.sh
export MARLIN_DLL=/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so


run=$1
echo " .... Marlin process SDHCAL asic processor"

#FILEPATH="/home/garillot/files/local/CalorimeterHit/PolyaMap/Geant4.10.01/FTF_BIC"
#FILEPATH="/home/garillot/files/local/CalorimeterHit/Geant4.10.01/FTF_BIC"
#FILELIST=`ls ${FILEPATH} | grep -F single_mu-_50GeV.slcio`
#COLLECTION="HCALEndcap"

#FILEPATH="/home/garillot/files/DATA/TRIVENT/thrScan"
#FILEPATH="/home/garillot/files/DATA/TRIVENT/SPS_Oct2015/163"
FILEPATH="/home/garillot/files/DATA/TRIVENT/SPS_Oct2015/214"
FILELIST=`ls ${FILEPATH} | grep -F TDHCAL_${run}`
COLLECTION="SDHCAL_HIT"


OUTPUT="map_${run}.root"
XMLFILENAME=$$.xml

cat > ${XMLFILENAME} <<EOF
<marlin>
 <execute>
  <processor name="EfficiencyProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   ${FILEPATH}/${FILELIST}
  </parameter>
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
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> $COLLECTION </parameter>
  <!--Name of the root output file-->
  <parameter name="RootFileName" type="string" > ${OUTPUT} </parameter>
  <parameter name="InteractionFinder::PrintDebug" type="bool"> false </parameter>
  <parameter name="Tracking::CosThetaLimit" type="float"> 0.9 </parameter>
  <parameter name="Tracking::PrintDebug" type="bool"> false </parameter>
 </processor>

</marlin>

EOF

Marlin ${XMLFILENAME}
rm ${XMLFILENAME}

#mv ${OUTPUT} /home/garillot/files/PolyaScan/DATA


