#include "AnalysisProcessor.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>

#include <ctime>
#include <array>
#include <numeric>

#include "Algorithm/Distance.h"
#include "EnergyOfRun.h"

using namespace lcio ;
using namespace marlin ;

AnalysisProcessor aAnalysisProcessor ;

AnalysisProcessor::AnalysisProcessor()
	: Processor("AnalysisProcessor") ,
	  _hcalCollections() ,
	  hitMap() ,
	  clusterVec() ,
	  edges() ,  //vector to recover geometry parameters
	  posShift() ,
	  algo_Cluster() ,
	  algo_ClusteringHelper() ,
	  algo_Tracking() ,
	  algo_Hough() ,
	  algo_InteractionFinder() ,
	  algo_density() ,
	  m_ClusterParameterSetting() ,
	  m_ClusteringHelperParameterSetting() ,
	  m_TrackingParameterSetting() ,
	  m_HoughParameterSetting() ,
	  m_InteractionFinderParameterSetting() ,
	  m_CaloGeomSetting() ,
	  m_ShowerAnalyserParameterSetting() ,
	  algo_ShowerAnalyser() ,
	  tracksClusterSize() ,
	  tracksClusterNumber() ,
	  longiProfile() ,
	  radiProfile()  ,
	  iVec() ,
	  jVec() ,
	  kVec() ,
	  thrVec()
{

	// modify processor description
	_description = "AnalysisProcessor" ;


	std::vector<std::string> hcalCollections;
	hcalCollections.push_back(std::string("HCALBarrel")) ;
	registerInputCollections( LCIO::CALORIMETERHIT,
							  "CollectionName" ,
							  "HCAL Collection Names"  ,
							  _hcalCollections  ,
							  hcalCollections);

	registerProcessorParameter( "RootFileName" ,
								"File name for the root output",
								outputRootName,
								std::string("toto.root") ) ;

	registerProcessorParameter( "nRun" ,
								"Number of run",
								runNumber,
								0 ) ;

	registerProcessorParameter( "recoverMissingDifsFile" ,
								"Xml file name for the missing difs to recover",
								recoverXmlFile,
								std::string("/home/garillot/SDHCALMarlinProcessor/recoverXML/H2Sept2017.xml") ) ;

	registerProcessorParameter( "NActiveLayers" ,
								"Number of active layers",
								_nActiveLayers,
								int(48) );

	posShift = CLHEP::Hep3Vector( 0.0 , 0.0 , 0.0 ) ;

	std::vector<float> thresholdsVec = { 1.0f , 2.0f , 3.0f } ;

	std::sort( thresholdsVec.begin() , thresholdsVec.end() ) ;

	registerProcessorParameter( "Thresholds" ,
								"Vector of thresholds",
								thresholds ,
								thresholdsVec ) ;

	AlgorithmRegistrationParameters() ;
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
	tree->Branch("computingTime", &computingTime) ;
	tree->Branch("eventTime" , &evtTime) ;
	tree->Branch("spillEventTime" , &spillEvtTime) ;
	tree->Branch("cerenkovTag" , &cerenkovTag) ;
	tree->Branch("energy" , &energy) ;

	tree->Branch("nHit" , &nHit) ;
	tree->Branch("nHit1" , &nHit1) ;
	tree->Branch("nHit2" , &nHit2) ;
	tree->Branch("nHit3" , &nHit3) ;

	tree->Branch("nHough" , &nHough) ;
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
	tree->Branch("end" , &_end) ;
	tree->Branch("density" , &density) ;

	tree->Branch("transverseRatio" , &transverseRatio) ;
	tree->Branch("reconstructedCosTheta" , &reconstructedCosTheta) ;

	tree->Branch("thrust" , &thrust , "thrust[4]/F") ;

	tree->Branch("meanRadius" , &meanRadius) ;
	tree->Branch("longiProfile" , "std::vector<double>" , &longiProfile) ;
	tree->Branch("radiProfile" , "std::vector<double>" , &radiProfile) ;

	tree->Branch("first5LayersRMS" , &first5LayersRMS) ;
	tree->Branch("neutral" , &neutral) ;

	tree->Branch("emFraction" , &emFraction) ;

	tree->Branch("I" , "std::vector<int>" , &iVec) ;
	tree->Branch("J" , "std::vector<int>" , &jVec) ;
	tree->Branch("K" , "std::vector<int>" , &kVec) ;
	tree->Branch("thr" , "std::vector<float>" , &thrVec) ;

	tree->Branch("nHitRecover" , &nHitRecover) ;
	tree->Branch("nHitRecover1" , &nHit1Recover) ;
	tree->Branch("nHitRecover2" , &nHit2Recover) ;
	tree->Branch("nHitRecover3" , &nHit3Recover) ;


	_timeCut = 5e9 ; //20 sec
	_prevBCID = 0 ;
	_bcidRef = 0 ;
	firstShowerInSpill = 1 ;
	firstSpillEvtFound = true ;

	_nRun = 0 ;
	_nEvt = 0 ;

	/*--------------------Algorithms initialisation--------------------*/
	algo_Cluster=new algorithm::Clustering() ;
	algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting) ;

	algo_ClusteringHelper=new algorithm::ClusteringHelper();
	algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting) ;

	algo_Tracking=new algorithm::Tracking();
	algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting) ;

	algo_Hough = new algorithm::Hough() ;
	algo_Hough->SetHoughParameterSetting(m_HoughParameterSetting) ;

	algo_InteractionFinder=new algorithm::InteractionFinder();
	algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderParameterSetting) ;

	m_ShowerAnalyserParameterSetting.energyCalibrationOption = std::string("sdhcal") ;
	m_ShowerAnalyserParameterSetting.interactionFinderParams = m_InteractionFinderParameterSetting ;
	m_ShowerAnalyserParameterSetting.geometry = m_CaloGeomSetting ;

	algo_ShowerAnalyser = new algorithm::ShowerAnalyser() ;
	algo_ShowerAnalyser->SetShowerAnalyserParameterSetting(m_ShowerAnalyserParameterSetting) ;

	algo_density = new algorithm::Density() ;

	processRecoverXmlFile() ;

}

void AnalysisProcessor::processRecoverXmlFile()
{
	TiXmlDocument doc(recoverXmlFile) ;
	if( !doc.LoadFile() )
	{
		std::cerr << "erreur lors du chargement" << std::endl ;
		std::cerr << "error #" << doc.ErrorId() << " : " << doc.ErrorDesc() << std::endl ;
		return ;
	}

	TiXmlHandle hdl(&doc);
	TiXmlElement* elem = hdl.FirstChildElement().FirstChildElement().Element() ;

	if( !elem )
	{
		std::cerr << "le noeud à atteindre n'existe pas" << std::endl ;
		return ;
	}

	while ( elem )
	{
		int id = std::atoi(elem->Attribute("id")) ;

		std::string s = elem->Attribute("dif") ;
		std::istringstream iss(s) ;
		std::vector<std::string> result{ std::istream_iterator<std::string>(iss) , {} } ;

		for ( const auto& i : result )
		{
			if ( i == std::string("left") ) recoverList.push_back( {id , left} ) ;
			else if ( i == std::string("center") ) recoverList.push_back( {id , center} ) ;
			else if ( i == std::string("right") ) recoverList.push_back( {id , right} ) ;
			else if ( i == std::string("all") ) recoverList.push_back( {id , all} ) ;
			else { std::cerr << "Error in recoverXmlFile : don't understand '" << i << "' for layer " << id << std::endl ; continue ; }
		}
		elem = elem->NextSiblingElement() ;
	}

	for ( const auto& i : recoverList )
		std::cout << "Will try to recover dif " << i.second << " of Layer " << i.first << std::endl ;
}


void AnalysisProcessor::processRunHeader(LCRunHeader*)
{
	_nRun++ ;
	_nEvt = 0 ;
}

void AnalysisProcessor::findEventTime(LCEvent* evt , LCCollection* _col)
{
	unsigned int hitTime = 0 ;
	EVENT::CalorimeterHit* hit = nullptr ;
	if ( _col->getNumberOfElements() != 0 )
	{
		try
		{
			hit = dynamic_cast<EVENT::CalorimeterHit*>( _col->getElementAt(0) ) ;
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
	_bcid = _bcid1*Shift + _bcid2 ;
	streamlog_out( DEBUG ) << "event : " << _nEvt+1 << " ; bcid: " << _bcid << " ; hitTime: " << hitTime <<std::endl;
	evtTime = _bcid - hitTime ;
}

void AnalysisProcessor::findSpillEventTime(LCEvent* evt , LCCollection* _col)
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
	if (_col->getNumberOfElements()!=0)
	{
		try
		{
			hit = dynamic_cast<EVENT::CalorimeterHit*>( _col->getElementAt(0) ) ;
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

		double rmsLay = sqrt( x2sum/n - (xsum*xsum)/(n*n) + y2sum/n - (ysum*ysum)/(n*n) ) ;
		if ( rmsLay > 10 )
			toReturn++ ;
	}
	return toReturn ;
}

std::vector<float> AnalysisProcessor::recoverHits() const
{
	std::vector<float> modifs(thresholds.size() + 1) ;

	for ( const auto& it : recoverList )
	{
		int layerBefore = std::max(0 , it.first - 1) ;
		int layerAfter = std::min(_nActiveLayers-1 , it.first + 1) ;

		const HitVec& before = hitMap.at(layerBefore) ;
		const HitVec& current = hitMap.at(it.first) ;
		const HitVec& after = hitMap.at(layerAfter) ;

		std::vector<int> nBefore(thresholds.size() + 1) ;
		std::vector<int> nCurrent(thresholds.size() + 1) ;
		std::vector<int> nAfter(thresholds.size() + 1) ;

		int min = recoverLimits.at(it.second).first ;
		int max = recoverLimits.at(it.second).second ;

		for ( caloobject::CaloHit* hit : before )
			if ( hit->getCellID()[1] >= min && hit->getCellID()[1] <= max )
				nBefore.at( static_cast<unsigned int>(hit->getEnergy()) )++ ;

		for ( caloobject::CaloHit* hit : after )
			if ( hit->getCellID()[1] >= min && hit->getCellID()[1] <= max )
				nAfter.at( static_cast<unsigned int>(hit->getEnergy()) )++ ;

		for ( caloobject::CaloHit* hit : current )
			if ( hit->getCellID()[1] >= min && hit->getCellID()[1] <= max )
				nCurrent.at( static_cast<unsigned int>(hit->getEnergy()) )-- ;

		for ( unsigned int i = 1 ; i < modifs.size() ; ++i )
		{
			modifs.at(0) += nCurrent.at(i) ;
			modifs.at(i) += nCurrent.at(i) ;

			modifs.at(0) += 0.5f*(nBefore.at(i) + nAfter.at(i)) ;
			modifs.at(i) += 0.5f*(nBefore.at(i) + nAfter.at(i)) ;
		}
	}
	return modifs ;
}


void AnalysisProcessor::processEvent( LCEvent * evt )
{
	clock_t beginClock = clock() ;

	for ( int i = 0 ; i < _nActiveLayers ; ++i )
		hitMap.insert( {i , {}} ) ;

	UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6") ;


	for (unsigned int i(0); i < _hcalCollections.size() ; ++i)
	{
		std::string colName =  _hcalCollections[i] ;
		try
		{
			col = evt->getCollection( _hcalCollections[i].c_str() ) ;
			numElements = col->getNumberOfElements();

			int NHIT = 0 ;
			for ( int j = 0 ; j < numElements ; ++j )
			{
				CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
				CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
				int cellID[] = { static_cast<int>( IDdecoder(hit)["I"]) , static_cast<int>( IDdecoder(hit)["J"]) , static_cast<int>( IDdecoder(hit)["K-1"]) } ;

				if ( cellID[2] >= _nActiveLayers )
					continue ;

				if ( cellID[0] < 1 || cellID[0] > 96 || cellID[1] < 1 || cellID[1] > 96 )
					continue ;

				caloobject::CaloHit* aHit = new caloobject::CaloHit( cellID , vec , hit->getEnergy() , hit->getTime() , posShift ) ;
				hitMap.at(cellID[2]).push_back(aHit) ;
				NHIT++ ;
			}

			if ( NHIT < 5 )
			{
				std::cout << "Too few Hits : only " << NHIT << std::endl ;
				clearVec() ;
				continue ;
			}

			findEventTime(evt,col) ;
			findSpillEventTime(evt,col) ;

			cerenkovTag = evt->parameters().getIntVal( "cerenkovTag" ) ;

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

			std::vector<float> modif = recoverHits() ;
			nHitRecover = nHit + modif.at(0) ;
			nHit1Recover = nHit1 + modif.at(1) ;
			nHit2Recover = nHit2 + modif.at(2) ;
			nHit3Recover = nHit3 + modif.at(3) ;

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
				begin = shower->getFirstIntCluster()->getLayerID() ;

			_end = shower->getLastClusterPosition() ;

			transverseRatio = shower->getTransverseRatio() ;
			reconstructedCosTheta = shower->getReconstructedCosTheta() ;

			for ( unsigned int jj = 0 ; jj < 4 ; jj++ )
				thrust[jj] = shower->getThrust().at(jj) ;

			algorithm::Distance<caloobject::CaloHit,CLHEP::Hep3Vector> dist ;

			meanRadius = 0 ;
			for ( int ii = 0 ; ii < _nActiveLayers ; ++ii )
			{
				for ( auto& cl : hitMap[ii] )
				{
					CLHEP::Hep3Vector vec(thrust[0] + thrust[1]*cl->getPosition().z() ,
							thrust[2] + thrust[3]*cl->getPosition().z() ,
							cl->getPosition().z() ) ;

					meanRadius += dist.getDistance(cl , vec) ;
				}
			}
			meanRadius /= nHit ;


			longiProfile = shower->getLongitudinal() ;
			for ( unsigned int lp = 0 ; lp < longiProfile.size() ; ++lp )
				longiProfile.at(lp) /= nHit ;

			radiProfile = shower->getTransverse() ;
			for ( unsigned int rp = 0 ; rp < radiProfile.size() ; ++rp )
				radiProfile.at(rp) /= nHit ;

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

			nHough = static_cast<int>( houghHitVec.size() ) ;
			nHough1 = 0 ;
			nHough2 = 0 ;
			nHough3 = 0 ;

			for ( HitVec::const_iterator it = houghHitVec.begin() ; it != houghHitVec.end() ; ++it )
			{
				if ( (*it)->getEnergy() >= thresholds.at(2) )
					nHough3++ ;
				else if ( (*it)->getEnergy() >= thresholds.at(1) )
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


			//density
			HitVec vec ;
			for( std::map<int,HitVec>::const_iterator it = hitMap.begin() ; it != hitMap.end() ; ++it )
				vec.insert(vec.end() , it->second.begin() , it->second.end() ) ;

			density = algo_density->compute(vec) ;

			eventNumber = _nEvt ;


			if ( evt->getParameters().getNFloat( std::string("ParticleEnergy") ) != 0 )
				energy = evt->getParameters().getFloatVal( std::string("ParticleEnergy") ) ;

			emFraction = evt->getParameters().getFloatVal( std::string("EMFraction") ) ;

			computingTime = 1.0*( clock() - beginClock )/CLOCKS_PER_SEC ;

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
	for( std::map<int,HitVec>::iterator it = hitMap.begin() ; it != hitMap.end() ; ++it )
		for( HitVec::iterator jt = (it->second).begin() ; jt != (it->second).end() ; ++jt )
			delete *jt ;

	hitMap.clear() ;
}


void AnalysisProcessor::check(LCEvent*)
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

	file->Write() ;
	file->Purge() ;
	file->Close() ;
}
