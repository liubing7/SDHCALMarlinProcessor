#!/bin/bash

source /home/garillot/ilcsoft/v01-17-08/init_ilcsoft.sh
export MARLIN_DLL=/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so

polyaQ=$1
polyaD=$2
d=$3

echo " .... Marlin process SDHCAL asic processor"

# FILEPATH="/home/garillot/files/local/CalorimeterHit/Geant4.10.01/FTF_BIC"
# FILELIST=`ls ${FILEPATH} | grep -F single_mu-_50GeV.slcio`
# COLLECTION="HCALEndcapAnalog"

FILEPATH="/home/garillot/files/PolyaScan/CalorimeterHit"
FILELIST=`ls ${FILEPATH} | grep -F "${polyaQ}_${polyaD}_${d}"`
COLLECTION="HCALEndcapAnalog"

OUTPUT="${polyaQ}_${polyaD}_${d}.root"

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

  <parameter name="Thresholds" type="FloatVec"> 0.114 0.14 0.155714 0.171429 0.187143 0.202857 0.218571 0.234286 0.25 0.265714 0.281429 0.298571 0.314286 0.33 0.345714 0.361429 0.377143 0.392857 0.408571 0.424286 0.4 0.6125 0.825 1.0375 1.2375 1.45 1.6625 1.875 2.0875 2.3 2.5 2.7125 2.925 3.1375 3.35 3.5625 3.7625 3.975 4.1875 4.29448 5.33742 6.38037 7.48466 8.52761 9.57055 10.6135 11.6564 12.6994 13.8037 14.8466 15.8896 16.9325 17.9755 19.0184 20.0613 21.1656 22.2086 23.2515 </parameter>

 </processor>

</marlin>

EOF

Marlin ${XMLFILENAME}
rm ${XMLFILENAME}

mv ${OUTPUT} /home/garillot/files/PolyaScan/MulResults
