#!/usr/bin/env python


import os
import sys
from os import path
import time

sys.path.insert(0 , '/home/garillot/SDHCALMarlinProcessor/script')

import EfficiencyProcessor


#if len(sys.argv) < 4 :
#	sys.exit('Error : too few arguments')

#qbar = sys.argv[1]
#delta = sys.argv[2]
#d = sys.argv[3]

#inputFileName = qbar + '_' + delta + '_' + d + '.slcio'
#inputDir = '/home/garillot/files/PolyaScan/CalorimeterHit'
#inputFilePath = inputDir + '/' + inputFileName


#inputFileName = 'digitMapTest.slcio'
#inputFileName = 'digitMapTest728359.slcio'
#inputFileName = 'digitAnalog.slcio'
#inputFileName = 'digitElecTest.slcio'
#inputDir = '/home/garillot/testMarlinReco/trunk'

inputDir = '/home/garillot/files/muonsForBing/CalorimeterHit'
inputFileName = 'single_mu-_50GeV.slcio'


inputFilePath = inputDir + '/' + inputFileName


fileList = [ inputFilePath ]
print 'Filelist : '
print fileList


os.environ["MARLIN"] = '/home/garillot/ilcsoft/v01-17-08/Marlin/v01-07'
os.environ["PATH"] = '/home/garillot/ilcsoft/v01-17-08/Marlin/v01-07/bin' + ':' + os.environ["PATH"]
os.environ["MARLIN_DLL"] = '/home/garillot/SDHCALMarlinProcessor/lib/libsdhcalMarlin.so'


a = EfficiencyProcessor.Params()
#a.collectionName = 'HCALEndcapAnalog'
#a.thresholds='0.114 6.12 16.83'
#a.thresholds = '0.114 0.14 0.155714 0.171429 0.187143 0.202857 0.218571 0.234286 0.25 0.265714 0.281429 0.298571 0.314286 0.33 0.345714 0.361429 0.377143 0.392857 0.408571 0.424286 0.4 0.6125 0.825 1.0375 1.2375 1.45 1.6625 1.875 2.0875 2.3 2.5 2.7125 2.925 3.1375 3.35 3.5625 3.7625 3.975 4.1875 4.29448 5.33742 6.38037 7.48466 8.52761 9.57055 10.6135 11.6564 12.6994 13.8037 14.8466 15.8896 16.9325 17.9755 19.0184 20.0613 21.1656 22.2086 23.2515'

#a.outputFileName = 'map_' + qbar + '_' + delta + '_' + d + '.root'
#a.outputFileName = 'map_digitMapTest.root'
#a.outputFileName = 'map_digitMapTestAnalog.root'
#a.outputFileName = 'map_digitMapElecTest.root'
#a.outputFileName = 'map_digitMapTest728359.root'

a.outputFileName = 'muonsBing.root'

EfficiencyProcessor.launch(a , fileList)


#outputDir = '/home/garillot/files/PolyaScan/MulResults'
#os.system('mv ' + a.outputFileName + ' ' + outputDir + '/' + a.outputFileName)
