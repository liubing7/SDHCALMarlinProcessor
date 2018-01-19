#include "AnalysisProcessorILD.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>

#include <ctime>

//#include "EnergyOfRun.h"

using namespace lcio ;
using namespace marlin ;

AnalysisProcessorILD aAnalysisProcessorILD ;

AnalysisProcessorILD::AnalysisProcessorILD()
	: Processor("AnalysisProcessorILD") ,
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
	  thresholdsFloat() ,
	  thresholds() ,
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
	_description = "AnalysisProcessorILD" ;


	std::vector<std::string> hcalCollections;
	hcalCollections.push_back(std::string("HCALOther"));
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

	registerProcessorParameter( "energy" ,
								"energy",
								energy ,
								0.0f ) ;

	registerProcessorParameter( "NActiveLayers" ,
								"Number of active layers",
								_nActiveLayers,
								48 );

	posShift = CLHEP::Hep3Vector( 0.0 , 0.0 , 0.0 ) ;

	std::vector<float> thresholdsVec = { 1.0f , 2.0f , 3.0f } ;

	std::sort( thresholdsVec.begin() , thresholdsVec.end() ) ;

	registerProcessorParameter( "Thresholds" ,
								"Vector of thresholds",
								thresholdsFloat,
								thresholdsVec ) ;

	AlgorithmRegistrationParameters() ;
}


void AnalysisProcessorILD::AlgorithmRegistrationParameters()
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
								true );

	registerProcessorParameter( "Cluster::MaxTransversalDistance" ,
								"Maximum transversal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxTransversalDistance,
								11.0f );

	registerProcessorParameter( "Cluster::MaxLongitudinalDistance" ,
								"Maximum longitudinal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxLongitudinalDistance,
								0.0f );

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
								10.0f );

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

void AnalysisProcessorILD::init()
{
	printParameters() ;
//	energy = energyOfRun( runNumber ) ;

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
	tree->Branch("longiProfile" , "std::vector<double>" , &longiProfile) ;
	tree->Branch("radiProfile" , "std::vector<double>" , &radiProfile) ;

	tree->Branch("first5LayersRMS" , &first5LayersRMS) ;
	tree->Branch("neutral" , &neutral) ;

	tree->Branch("emFraction" , &emFraction) ;

	tree->Branch("I" , "std::vector<int>" , &iVec) ;
	tree->Branch("J" , "std::vector<int>" , &jVec) ;
	tree->Branch("K" , "std::vector<int>" , &kVec) ;
	tree->Branch("thr" , "std::vector<int>" , &thrVec) ;
	tree->Branch("hitDensity" , "std::vector<float>" , &densityVec) ;

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

	algo_density = new algorithm::Density() ;

}

void AnalysisProcessorILD::processRunHeader(LCRunHeader*)
{
	_nRun++ ;
	_nEvt = 0 ;
}

void AnalysisProcessorILD::findEventTime(LCEvent* evt , LCCollection* _col)
{
	unsigned int hitTime = 0 ;
	EVENT::CalorimeterHit* hit = NULL ;
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
	_bcid=_bcid1*Shift+_bcid2;
	streamlog_out( DEBUG ) << "event : " << _nEvt+1 << " ; bcid: " << _bcid << " ; hitTime: " << hitTime <<std::endl;
	evtTime = _bcid - hitTime ;
}

void AnalysisProcessorILD::findSpillEventTime(LCEvent* evt , LCCollection* _col)
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

double AnalysisProcessorILD::getFirst5LayersRMS()
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

void AnalysisProcessorILD::processEvent( LCEvent * evt )
{
	clock_t beginClock = clock() ;


	//
	// * Reading HCAL Collections of CalorimeterHits*
	//
	//std::string initString;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("system:0:5,module:5:3,stave:8:3,tower:11:5,layer:16:6,slice:22:4,x:32:-16,z:48:-16") ;


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
				int cellID[] = { static_cast<int>( IDdecoder(hit)["x"]) , static_cast<int>( IDdecoder(hit)["z"]) , static_cast<int>( IDdecoder(hit)["layer"]-1) } ;

				if ( cellID[2] > _nActiveLayers )
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

			caloobject::DigitalShower* shower = new caloobject::DigitalShower(clusterVec) ;
			shower->setInteractionSettings(m_InteractionFinderParameterSetting) ;
			shower->setGeometrySettings(m_CaloGeomSetting) ;

			shower->computePCA() ;
			shower->computeThrust() ;
			shower->computeInteraction() ;

			shower->computeProfile() ;

			nHit = shower->getNHits().at(0) ;
			nHit1 = shower->getNHits().at(1) ;
			nHit2 = shower->getNHits().at(2) ;
			nHit3 = shower->getNHits().at(3) ;
			nLayer = static_cast<int>( shower->getFiredLayers().size() ) ;
			nInteractingLayer = static_cast<int>( shower->getInteractingLayers().size() ) ;

			iVec.clear() ;
			jVec.clear() ;
			kVec.clear() ;
			thrVec.clear() ;
			densityVec.clear() ;

			//density
			HitVec hitVec ;
			for( auto it : hitMap )
				hitVec.insert(hitVec.end() , it.second.begin() , it.second.end() ) ;

			density = algo_density->compute(hitVec) ;

			auto densityPerHit = algo_density->getDensityPerHit() ;

			for ( auto hit : shower->getHits() )
			{
				iVec.push_back( hit->getCellID()[0] ) ;
				jVec.push_back( hit->getCellID()[1] ) ;
				kVec.push_back( hit->getCellID()[2] ) ;
				thrVec.push_back( hit->getEnergy() ) ;

				densityVec.push_back( densityPerHit.at(hit) ) ;
			}

			if ( !shower->getFirstIntCluster() )
				begin = -10 ;
			else
				begin = shower->getFirstIntCluster()->getLayerID() ;
			//				begin = shower->getStartingPosition()[2] ;

			_end = shower->getLastClusterLayer() ;

			transverseRatio = shower->getTransverseRatio() ;
			reconstructedCosTheta = shower->getReconstructedCosTheta() ;

			for ( unsigned int jj = 0 ; jj < 4 ; jj++ )
				thrust[jj] = shower->getThrust().at(jj) ;


			longiProfile = shower->getLongitudinalProfile() ;
			for ( unsigned int lp = 0 ; lp < longiProfile.size() ; ++lp )
				longiProfile.at(lp) /= nHit ;

			radiProfile = shower->getTransverseProfile() ;
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
				if ( (*it)->getEnergy() >= thresholdsFloat.at(2) )
					nHough3++ ;
				else if ( (*it)->getEnergy() >= thresholdsFloat.at(1) )
					nHough2++ ;
				else
					nHough1++ ;
			}

			nCluster = static_cast<int>( clusterVec.size() ) ;

			if ( clusterVec.at(0)->getPosition()[2] > ( 4*m_CaloGeomSetting.layerGap + 2670 ) )
				neutral = true ;
			else
				neutral = false ;

			first5LayersRMS = getFirst5LayersRMS() ;


			eventNumber = _nEvt ;


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

void AnalysisProcessorILD::clearVec()
{
	for( std::map<int,HitVec>::iterator it = hitMap.begin() ; it!=hitMap.end() ; ++it )
		for( HitVec::iterator jt = (it->second).begin() ; jt != (it->second).end() ; ++jt )
			delete *jt ;

	hitMap.clear() ;
}


void AnalysisProcessorILD::check(LCEvent*)
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AnalysisProcessorILD::end()
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
