#include "AnalysisProcessor.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>

#include "EnergyOfRun.h"

using namespace lcio ;
using namespace marlin ;

AnalysisProcessor aAnalysisProcessor ;

AnalysisProcessor::AnalysisProcessor() : Processor("AnalysisProcessor")
{

	// modify processor description
	_description = "AnalysisProcessor" ;


	std::vector<std::string> hcalCollections;
	hcalCollections.push_back(std::string("HCALBarrel"));
	registerInputCollections( LCIO::CALORIMETERHIT,
							  "CollectionName" ,
							  "HCAL Collection Names"  ,
							  _hcalCollections  ,
							  hcalCollections);

	registerProcessorParameter( "RootFileName" ,
								"File name for the root output",
								outputRootName,
								std::string("toto.root") );

	registerProcessorParameter( "nRun" ,
								"Number of run",
								runNumber,
								0 ) ;

	registerProcessorParameter( "NActiveLayers" ,
								"Number of active layers",
								_nActiveLayers,
								int(48) );

	registerProcessorParameter( "N_ASIC" ,
								"Number of ASIC per layer in x direction",
								_nAsicX,
								int(12) );

	registerProcessorParameter( "N_ASIC" ,
								"Number of ASIC per layer in y direction",
								_nAsicY,
								int(12) );

	int difTab[]={ 181,94,30, 174,175,176, 158,142,141, 129,118,119, 164,152,151,  74,61,75,
				   156,111,110, 102,177,103,  133,136,134,  128,120,121,  65,64,58,  148,72,73,
				   78,79,60,  44,43,113,  243,242,241,   186,127,154,  147,70,71,   47,139,140,
				   143,77,76,   159,91,36,   179,178,183,  41,42,67,  137,46,138,  131,173,144,
				   189,184,160,  172,167,171,  146,135,145,  185,170,180,  187,188,190,  169,165,166,
				   155,57,50,  153,108,25,   51,56,109,   107,150,116,  126,124,49,  117,149,115,
				   48,45,114,   98,93,40,   92,97,100,  62,106,132,  101,35,99,  122,123,130,
				   163,161,162,  104,29,112,  59,53,54,  96,90,27,  95,8,5,  63,87,18 };
	std::vector<int> difVec(difTab, difTab + sizeof(difTab) / sizeof(int) );

	registerProcessorParameter( "DifList" ,
								"Vector of dif number",
								_difList,
								difVec );


	std::vector<float> vec ;
	vec.push_back(499.584f) ;
	vec.push_back(499.584f) ;
	vec.push_back(0) ;
	// registerProcessorParameter( "PositionShift" ,
	// 			      "3 Vector to shift to have the right (0,0,0) position",
	// 			      _posShift,
	// 			      vec );
	posShift = CLHEP::Hep3Vector( 0.0, 0.0, 0.0 ) ;

	float thr[] = { 1.0f , 2.0f , 3.0f } ; //default value for semi-digital hits
	std::vector<float> thresholdsVec(thr, thr + sizeof(thr) / sizeof(float) ) ;

	std::sort( thresholdsVec.begin() , thresholdsVec.end() ) ;

	registerProcessorParameter( "Thresholds" ,
								"Vector of thresholds",
								thresholdsFloat,
								thresholdsVec ) ;

	AlgorithmRegistrationParameters();
}


void AnalysisProcessor::AlgorithmRegistrationParameters()
{
	/*------------algorithm::Cluster------------*/
	registerProcessorParameter( "Cluster::MaxTransversalCellID" ,
								"Maximum difference between two hits cellID (0 and 1) to build a cluster",
								m_ClusterParameterSetting.maxTransversal,
								1 );

	registerProcessorParameter( "Cluster::MaxLongitudinalCellID" ,
								"Maximum difference between two hits cellID (2) to build a cluster",
								m_ClusterParameterSetting.maxLongitudinal,
								0 );

	registerProcessorParameter( "Cluster::UseDistanceInsteadCellID" ,
								"Boolean to know if clustering algorithms uses distance instead of cellID to cluster hits together",
								m_ClusterParameterSetting.useDistanceInsteadCellID,
								false );

	registerProcessorParameter( "Cluster::MaxTransversalDistance" ,
								"Maximum transversal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxTransversalDistance,
								11.0f );

	registerProcessorParameter( "Cluster::MaxLongitudinalDistance" ,
								"Maximum longitudinal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxLongitudinalDistance,
								27.0f );

	/*------------algorithm::ClusteringHelper------------*/
	registerProcessorParameter( "ClusteringHelper::LongitudinalDistanceForIsolation" ,
								"Minimum longitudinal distance (in mm) between one hits and its neighbours to decide if it is isolated",
								m_ClusteringHelperParameterSetting.longitudinalDistance,
								100.0f );

	registerProcessorParameter( "ClusteringHelper::TransversalDistanceDistanceForIsolation" ,
								"Minimum transversal distance (in mm) between one hits and its neighbours to decide if it is isolated",
								m_ClusteringHelperParameterSetting.transversalDistance,
								200.0f );

	/*------------algorithm::Tracking-----------*/
	registerProcessorParameter( "Tracking::ChiSquareLimit" ,
								"Maximum value of tracking fit chi2 to construct a track",
								m_TrackingParameterSetting.chiSquareLimit,
								100.0f );

	registerProcessorParameter( "Tracking::MaxTransverseRatio" ,
								"Maximum value of transverse ratio to construct a track",
								m_TrackingParameterSetting.maxTransverseRatio,
								0.05f );
	//								(float) 1 );

	registerProcessorParameter( "Tracking::CosThetaLimit" ,
								"Minimum value of cos(Theta) to accept the track",
								m_TrackingParameterSetting.cosThetaLimit,
								0.0f );

	registerProcessorParameter( "Tracking::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_TrackingParameterSetting.printDebug,
								false );

	/*------------algorithm::Efficiency-----------*/
	registerProcessorParameter( "Efficiency::MaxRadius" ,
								"Maximum distance parameter to find a hit to consider the layer as efficient",
								m_EfficiencyParameterSetting.maxRadius,
								25.0f );

	registerProcessorParameter( "Efficiency::SDHCALReadout" ,
								"Boolean to know if the detector used the semi digital readout",
								m_EfficiencyParameterSetting.semiDigitalReadout,
								true );

	registerProcessorParameter( "Efficiency::PrintDebug" ,
								"If true, Efficiency algorithm will print some debug information",
								m_EfficiencyParameterSetting.printDebug,
								false );

	m_EfficiencyParameterSetting.trackingParams=m_TrackingParameterSetting;

	/*------------algorithm::InteractionFinder-----------*/
	registerProcessorParameter( "InteractionFinder::MinSize" ,
								"Minimum cluster size for to define an interaction point",
								m_InteractionFinderParameterSetting.minSize,
								4 );

	registerProcessorParameter( "InteractionFinder::MaxRadius" ,
								"Maximum transversal distance to look for clusters",
								m_InteractionFinderParameterSetting.maxRadius,
								50.0f );

	registerProcessorParameter( "InteractionFinder::MaxDepth" ,
								"Maximum depth (number of layers) to look for clusters",
								m_InteractionFinderParameterSetting.maxDepth,
								4 );

	registerProcessorParameter( "InteractionFinder::MinNumberOfCluster" ,
								"Minimum number of found clusters (big enough) after the interaction point",
								m_InteractionFinderParameterSetting.minNumberOfCluster,
								3 );

	registerProcessorParameter( "InteractionFinder::UseAnalogEnergy" ,
								"Boolean to know if interaction finder algo should use cluster energy of cluster number of hits",
								m_InteractionFinderParameterSetting.useAnalogEnergy,
								false );

	registerProcessorParameter( "InteractionFinder::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_InteractionFinderParameterSetting.printDebug,
								false );

	/*------------caloobject::CaloGeom------------*/
	registerProcessorParameter( "Geometry::NLayers" ,
								"Number of layers",
								m_CaloGeomSetting.nLayers,
								48 );
	registerProcessorParameter( "Geometry::NPixelsPerLayer" ,
								"Number of pixels per layer (assume square geometry)",
								m_CaloGeomSetting.nPixelsPerLayer,
								96 );
	registerProcessorParameter( "Geometry::PixelSize" ,
								"Pixel size (assume square pixels)",
								m_CaloGeomSetting.pixelSize,
								10.408f );

	std::vector<float> vec;
	vec.push_back(10.408f) ;
	vec.push_back(97*10.408f) ;
	registerProcessorParameter( "Geometry::DetectorTransverseSize" ,
								"Define the detector transverse size used by efficiency algorithm (vector size must be 2 or 4; if 2 -> first value is min, second value is max; if 4 -> two first values define x edges , two last values define y edges) ",
								edges,
								vec );
	if( edges.size()==2 ){
		m_CaloGeomSetting.xmin=edges[0];
		m_CaloGeomSetting.ymin=edges[0];
		m_CaloGeomSetting.xmax=edges[1];
		m_CaloGeomSetting.ymax=edges[1];
	}
	else if( edges.size()==4 ){
		m_CaloGeomSetting.xmin=edges[0];
		m_CaloGeomSetting.xmax=edges[1];
		m_CaloGeomSetting.ymin=edges[2];
		m_CaloGeomSetting.ymax=edges[3];
	}
	else{
		std::cout << "WARING : Wrong number of values in paramater Geometry::DetectorTransverseSize => will use default value -500.0, +500.0" << std::endl;
	}
	/*--------------------------------------------*/

	/*------------algorithm::AsicKeyFinder-----------*/
	std::vector<int> asicKeyFactor;
	asicKeyFactor.push_back(1);
	asicKeyFactor.push_back(12);
	asicKeyFactor.push_back(1000);
	registerProcessorParameter( "AsicKeyFinder::KeyFactor" ,
								"Define the factor to build asic keys; default for sdhcal : 1000*K + 12*J + I",
								m_AsicKeyFinderParameterSetting.keyFactors,
								asicKeyFactor );

	registerProcessorParameter( "AsicKeyFinder::NPadX" ,
								"Number of pads in x direction per layer",
								m_AsicKeyFinderParameterSetting.nPadX,
								96 );

	registerProcessorParameter( "AsicKeyFinder::NPadY" ,
								"Number of pads in x direction per layer",
								m_AsicKeyFinderParameterSetting.nPadY,
								96 );

	registerProcessorParameter( "AsicKeyFinder::AsicNPad" ,
								"number of pads in x or y direction per asic (assuming a square)",
								m_AsicKeyFinderParameterSetting.asicNPad,
								8 );

	registerProcessorParameter( "AsicKeyFinder::LayerGap" ,
								"Gap size (in mm) between 2 layers",
								m_AsicKeyFinderParameterSetting.layerGap,
								26.131f );

	registerProcessorParameter( "AsicKeyFinder::PadSize" ,
								"Size of one pad in mm",
								m_AsicKeyFinderParameterSetting.padSize,
								10.408f );

	registerProcessorParameter( "AsicKeyFinder::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_AsicKeyFinderParameterSetting.printDebug,
								false );
}

void AnalysisProcessor::init()
{
	printParameters() ;
	energy = energyOfRun( runNumber ) ;

	file = new TFile(outputRootName.c_str(),"RECREATE") ;

	tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
	if( !tree )
	{
		std::cout << "tree creation" << std::endl ;
		tree = new TTree("tree","Shower variables") ;
	}

	tree->Branch("eventNumber", &eventNumber) ;
	tree->Branch("eventTime" , &evtTime) ;
	tree->Branch("spillEventTime" , &spillEvtTime) ;
	tree->Branch("cerenkovTag" , &cerenkovTag) ;
	tree->Branch("energy" , &energy) ;

	tree->Branch("nHit" , &nHit) ;
	tree->Branch("nHit1" , &nHit1) ;
	tree->Branch("nHit2" , &nHit2) ;
	tree->Branch("nHit3" , &nHit3) ;

	tree->Branch("nHough" , &nHit) ;
	tree->Branch("nHough1" , &nHough1) ;
	tree->Branch("nHough2" , &nHough2) ;
	tree->Branch("nHough3" , &nHough3) ;

	tree->Branch("nLayer" , &nLayer) ;
	tree->Branch("nInteractingLayer" , &nInteractingLayer) ;
	tree->Branch("nCluster" , &nCluster) ;
	tree->Branch("nMipCluster" , &nMipCluster) ;
	tree->Branch("nTrack" , &nTrack) ;

	tree->Branch("tracksClusterSize" , "std::vector<int>" , &tracksClusterSize) ;
	tree->Branch("tracksClusterNumber" , "std::vector<int>" , &tracksClusterNumber) ;

	tree->Branch("begin" , &begin) ;

	tree->Branch("transverseRatio" , &transverseRatio) ;
	tree->Branch("reconstructedCosTheta" , &reconstructedCosTheta) ;

	tree->Branch("thrust" , &thrust , "thrust[4]/F") ;
	tree->Branch("longiProfile" , "std::vector<double>" , &longiProfile) ;
	tree->Branch("radiProfile" , "std::vector<double>" , &radiProfile) ;

	tree->Branch("first5LayersRMS" , &first5LayersRMS) ;
	tree->Branch("neutral" , &neutral) ;

	tree->Branch("I" , "std::vector<int>" , &iVec) ;
	tree->Branch("J" , "std::vector<int>" , &jVec) ;
	tree->Branch("K" , "std::vector<int>" , &kVec) ;
	tree->Branch("thr" , "std::vector<int>" , &thrVec) ;

	_timeCut = 5e9 ; //20 sec
	_prevBCID = 0 ;
	_bcidRef = 0 ;
	firstShowerInSpill = 1 ;
	firstSpillEvtFound = true ;

	_nRun = 0 ;
	_nEvt = 0 ;

	/*--------------------Algorithms initialisation--------------------*/
	algo_Cluster=new algorithm::Cluster() ;
	algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting) ;

	algo_ClusteringHelper=new algorithm::ClusteringHelper();
	algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting) ;

	algo_Tracking=new algorithm::Tracking();
	algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting) ;

	algo_Hough = new algorithm::Hough() ;
	algo_Hough->SetHoughParameterSetting(m_HoughParameterSetting) ;

	algo_InteractionFinder=new algorithm::InteractionFinder();
	algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderParameterSetting) ;

	m_EfficiencyParameterSetting.geometry=m_CaloGeomSetting;
	algo_Efficiency=new algorithm::Efficiency();
	algo_Efficiency->SetEfficiencyParameterSetting(m_EfficiencyParameterSetting) ;

	algo_AsicKeyFinder=new algorithm::AsicKeyFinder();
	algo_AsicKeyFinder->SetAsicKeyFinderParameterSetting(m_AsicKeyFinderParameterSetting) ;

	m_ShowerAnalyserParameterSetting.energyCalibrationOption = std::string("sdhcal") ;
	m_ShowerAnalyserParameterSetting.interactionFinderParams = m_InteractionFinderParameterSetting ;
	m_ShowerAnalyserParameterSetting.geometry = m_CaloGeomSetting ;

	algo_ShowerAnalyser = new algorithm::ShowerAnalyser() ;
	algo_ShowerAnalyser->SetShowerAnalyserParameterSetting(m_ShowerAnalyserParameterSetting) ;

}

void AnalysisProcessor::processRunHeader( LCRunHeader* run)
{
	_nRun++ ;
	_nEvt = 0 ;
}

void AnalysisProcessor::findEventTime(LCEvent* evt , LCCollection* col)
{
	unsigned int hitTime = 0 ;
	EVENT::CalorimeterHit* hit = NULL ;
	if ( col->getNumberOfElements() != 0 )
	{
		try
		{
			hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt(0) ) ;
			hitTime = uint (hit->getTime() ) ;

		}
		catch (std::exception e)
		{
			std::cout << "No hits " << std::endl ;
			return ;
		}
	}

	unsigned long long _bcid ;
	unsigned long long _bcid1 ;
	unsigned long long _bcid2 ;
	std::stringstream pname1("") ;
	pname1 << "bcid1";
	_bcid1 = evt->parameters().getIntVal( pname1.str() ) ;
	std::stringstream pname2("") ;
	pname2 << "bcid2";
	_bcid2 = evt->parameters().getIntVal( pname2.str() ) ;

	unsigned long long Shift = 16777216ULL;
	_bcid=_bcid1*Shift+_bcid2;
	streamlog_out( DEBUG ) << "event : " << _nEvt+1 << " ; bcid: " << _bcid << " ; hitTime: " << hitTime <<std::endl;
	evtTime = _bcid - hitTime ;
}

void AnalysisProcessor::findSpillEventTime(LCEvent* evt,LCCollection* col)
{
	unsigned long long _bcid;
	unsigned long long _bcid1;
	unsigned long long _bcid2;
	std::stringstream pname1("");
	pname1 << "bcid1";
	_bcid1=evt->parameters().getIntVal(pname1.str());
	std::stringstream pname2("");
	pname2 << "bcid2";
	_bcid2=evt->parameters().getIntVal(pname2.str());

	unsigned long long Shift=16777216ULL;
	_bcid=_bcid1*Shift+_bcid2; //trigger time

	unsigned int hitTime = 0 ;
	EVENT::CalorimeterHit* hit=NULL;
	if (col->getNumberOfElements()!=0)
	{
		try
		{
			hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt(0) ) ;
			hitTime = uint( hit->getTime() ) ;

		}
		catch (std::exception e)
		{
			std::cout<<"No hits "<<std::endl;
			return ;
		}
	}

	//_bcidRef = absolute bcid of 1st pysical event in spill
	if(_prevBCID==0)
	{
		spillEvtTime=_bcid;
		_bcidRef=0;
		streamlog_out( DEBUG ) << "event : " << _nEvt+1
							   << " ; first event time : " << spillEvtTime
							   << " ; first reference : " << _bcidRef
							   << std::endl;
		if(numElements<400)
		{
			firstShowerInSpill=0;
			firstSpillEvtFound=false;
		}
	}
	else if( (_bcid-_prevBCID)*200 < _timeCut )
	{
		spillEvtTime=_bcid-_bcidRef;
		streamlog_out( DEBUG ) << "event : " << _nEvt+1
							   << " ; reference : " << _bcidRef
							   << " ; time to the spill start : " << spillEvtTime
							   << std::endl;
		if(firstSpillEvtFound==false&&numElements>400)
		{
			firstShowerInSpill=1;
			firstSpillEvtFound=true;
		}
		else
			firstShowerInSpill=0;
	}
	else
	{
		_bcidRef = _bcid ;
		spillEvtTime = hitTime ;
		streamlog_out( DEBUG ) << "event : " << _nEvt+1
							   << " ; distance to the previous event : " << _bcid-_prevBCID
							   << " ; New reference : " << _bcidRef
							   << " ; time to the spill start : " << spillEvtTime
							   << std::endl;
		if(numElements>400)
		{
			firstShowerInSpill=1;
			firstSpillEvtFound=true;
		}
		else{
			firstShowerInSpill=0;
			firstSpillEvtFound=false;
		}
	}
	_prevBCID = _bcid ;
}

double AnalysisProcessor::getFirst5LayersRMS()
{
	std::sort(clusterVec.begin() , clusterVec.end() , algorithm::ClusteringHelper::SortClusterByLayer) ;

	double xsum = 0 , x2sum = 0 , ysum = 0 , y2sum = 0 ;
	int n = 0 ;

	for( std::vector<caloobject::CaloCluster2D*>::const_iterator it = clusterVec.begin() ; it != clusterVec.end() ; ++it )
	{
		if ( (*it)->getLayerID() > 5 )
			break ;

		double x = (*it)->getPosition()[0] ;
		double y = (*it)->getPosition()[1] ;
		xsum += x ;
		x2sum += x*x ;
		ysum += y ;
		y2sum += y*y ;

		n++ ;
	}

	if ( n < 2 )
		return 0 ;

	return sqrt( x2sum/(n-1) - (xsum*xsum)/(n*(n-1)) + y2sum/(n-1) - (ysum*ysum)/(n*(n-1)) )  ;
}

int AnalysisProcessor::getNInteractingLayer()
{
	int toReturn = 0 ;

	for( std::map<int,std::vector<caloobject::CaloHit*> >::iterator it = hitMap.begin() ; it!=hitMap.end() ; ++it )
	{
		double xsum = 0 , x2sum = 0 , ysum = 0 , y2sum = 0 ;

		for( std::vector<caloobject::CaloHit*>::iterator jt = (it->second).begin() ; jt != (it->second).end() ; ++jt )
		{
			double x = (*jt)->getPosition()[0] ;
			double y = (*jt)->getPosition()[1] ;
			xsum += x ;
			x2sum += x*x ;
			ysum += y ;
			y2sum += y*y ;
		}
		double n = (it->second).size() ;
		if ( n <= 5 )
			continue ;

		double rmsLay = sqrt( x2sum/n - (xsum*xsum)/(n*n) + y2sum/n - (ysum*ysum)/(n*n) )  ;
		if ( rmsLay > 10 )
			toReturn++ ;
	}
	return toReturn ;
}

void AnalysisProcessor::processEvent( LCEvent * evt )
{
	//
	// * Reading HCAL Collections of CalorimeterHits*
	//
	//std::string initString;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6") ;


	for (unsigned int i(0); i < _hcalCollections.size(); ++i)
	{
		std::string colName =  _hcalCollections[i] ;
		try
		{
			col = evt->getCollection( _hcalCollections[i].c_str() ) ;
			//initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
			numElements = col->getNumberOfElements();
			//      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);

			int NHIT = 0 ;
			for (int j=0; j < numElements; ++j)
			{
				CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
				CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
				int cellID[] = { static_cast<int>( IDdecoder(hit)["I"]) , static_cast<int>( IDdecoder(hit)["J"]) , static_cast<int>( IDdecoder(hit)["K-1"]) } ;

				if ( cellID[2] > _nActiveLayers )
					continue ;

				if ( cellID[0] < 1 || cellID[0] > 96 || cellID[1] < 1 || cellID[1] > 96)
					continue ;

				caloobject::CaloHit* aHit = new caloobject::CaloHit( cellID , vec , hit->getEnergy() , hit->getTime() , posShift ) ;
				hitMap[ cellID[2] ].push_back(aHit) ;
				NHIT++ ;

				/*
				if( hit->getPosition()[0]<0 || hit->getPosition()[1]<0 || hit->getPosition()[2]<0 )
				{
					std::cout << "WARNING : hit at "
							  << hit->getPosition()[0] << ",\t"
							  << hit->getPosition()[1] << ",\t"
							  << hit->getPosition()[2] << std::endl;
					getchar() ;
				}
				*/
			}

			if ( NHIT < 5 )
			{
				std::cout << "Too few Hits : only " << NHIT << std::endl ;
				clearVec() ;
				continue ;
			}


			findEventTime(evt,col) ;
			findSpillEventTime(evt,col) ;
			std::stringstream bifTag("") ;
			bifTag << "cerenkovTag" ;
			cerenkovTag = evt->parameters().getIntVal( bifTag.str() ) ;

			clusterVec.clear() ;

			for ( int ii = 0 ; ii < _nActiveLayers ; ++ii )
				algo_Cluster->Run(hitMap[ii] , clusterVec) ;

			std::sort(clusterVec.begin() , clusterVec.end() , algorithm::ClusteringHelper::SortClusterByLayer) ;

			caloobject::SDHCALShower* shower = new caloobject::SDHCALShower(clusterVec) ;

			algo_ShowerAnalyser->Run(shower) ;

			nHit = shower->getNhit() ;
			nHit1 = shower->getSDNHits()[0] ;
			nHit2 = shower->getSDNHits()[1] ;
			nHit3 = shower->getSDNHits()[2] ;
			nLayer = shower->getNlayer() ;

			iVec.clear() ;
			jVec.clear() ;
			kVec.clear() ;
			thrVec.clear() ;

			for ( HitVec::const_iterator it = shower->getHits().begin() ; it != shower->getHits().end() ; ++it )
			{
				iVec.push_back( (*it)->getCellID()[0] ) ;
				jVec.push_back( (*it)->getCellID()[1] ) ;
				kVec.push_back( (*it)->getCellID()[2] ) ;
				thrVec.push_back( (*it)->getEnergy() ) ;
			}

			if ( !shower->findInteraction() )
				begin = -10 ;
			else
				begin = shower->getStartingPosition()[2] ;

			transverseRatio = shower->getTransverseRatio() ;
			reconstructedCosTheta = shower->getReconstructedCosTheta() ;

			for ( unsigned int jj = 0 ; jj < 4 ; jj++ )
				thrust[jj] = shower->getThrust().at(jj) ;

			longiProfile = shower->getLongitudinal() ;
			radiProfile = shower->getTransverse() ;

			std::vector<caloobject::CaloTrack*> trackVec ;
			algo_Hough->runHough(clusterVec , trackVec , algo_Tracking) ;

			std::vector<caloobject::CaloCluster2D*> clusterMipVec ;
			algo_Hough->selectNonDensePart(clusterVec , clusterMipVec) ;

			nMipCluster = static_cast<int>( clusterMipVec.size() ) ;

			nTrack = static_cast<int>( trackVec.size() ) ;

			tracksClusterSize.clear() ;
			tracksClusterNumber.clear() ;

			HitVec houghHitVec ;

			for ( std::vector<caloobject::CaloTrack*>::const_iterator it = trackVec.begin() ; it != trackVec.end() ; ++it )
			{
				std::vector<caloobject::CaloCluster2D*> cl = (*it)->getClusters() ;
				tracksClusterNumber.push_back( static_cast<int>( cl.size() ) ) ;

				for ( std::vector<caloobject::CaloCluster2D*>::const_iterator jt = cl.begin() ; jt != cl.end() ; ++jt )
				{
					houghHitVec.insert( houghHitVec.end() , (*jt)->getHits().begin() , (*jt)->getHits().end() ) ;
					tracksClusterSize.push_back( static_cast<int>( (*jt)->getHits().size() ) ) ;
				}
			}

			nHough = nHough1 = nHough2 = nHough3 = 0 ;
			for ( HitVec::const_iterator it = houghHitVec.begin() ; it != houghHitVec.end() ; ++it )
			{
				nHough++ ;
				if ( (*it)->getEnergy() >= thresholdsFloat.at(2) )
					nHough3++ ;
				else if ( (*it)->getEnergy() >= thresholdsFloat.at(1) )
					nHough2++ ;
				else
					nHough1++ ;
			}

			nInteractingLayer = getNInteractingLayer() ;

			nCluster = static_cast<int>( clusterVec.size() ) ;

			if ( clusterVec.at(0)->getPosition()[2] > 4*m_CaloGeomSetting.layerGap )
				neutral = true ;
			else
				neutral = false ;

			first5LayersRMS = getFirst5LayersRMS() ;

			eventNumber = _nEvt ;


			if ( evt->getParameters().getNFloat( std::string("ParticleEnergy") ) != 0 )
				energy = evt->getParameters().getFloatVal( std::string("ParticleEnergy") ) ;


			tree->Fill() ;

			delete shower ;

			for( std::vector<caloobject::CaloTrack*>::iterator it = trackVec.begin() ; it != trackVec.end() ; ++it )
				delete *it ;
			for( std::vector<caloobject::CaloCluster2D*>::iterator it = clusterVec.begin() ; it != clusterVec.end() ; ++it )
				delete *it ;

			clearVec() ;
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << "Exception " << std::endl ;
		}
	}
	_nEvt ++ ;
	std::cout << "Event processed : " << _nEvt << std::endl ;
}

void AnalysisProcessor::clearVec()
{
	for( std::map<int,HitVec>::iterator it = hitMap.begin() ; it!=hitMap.end() ; ++it )
		for( HitVec::iterator jt = (it->second).begin() ; jt != (it->second).end() ; ++jt )
			delete *jt ;

	hitMap.clear() ;
}


void AnalysisProcessor::check( LCEvent * evt )
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AnalysisProcessor::end()
{
	file->cd() ;

	delete algo_Cluster ;
	delete algo_ClusteringHelper ;
	delete algo_Tracking ;
	delete algo_InteractionFinder ;
	delete algo_Efficiency ;
	delete algo_AsicKeyFinder ;

	file->Write() ;
	file->Purge() ;
	file->Close() ;
}
