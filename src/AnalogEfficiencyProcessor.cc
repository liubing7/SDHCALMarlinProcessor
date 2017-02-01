
#include "AnalogEfficiencyProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include <marlin/Global.h>
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

AnalogEfficiencyProcessor aAnalogEfficiencyProcessor ;

AnalogEfficiencyProcessor::AnalogEfficiencyProcessor() : Processor("AnalogEfficiencyProcessor")
{

	// modify processor description
	_description = "sdhcalAsicProcessor calculates a SDHCAL efficiency and pad multiplcity for each SHDCAL layer" ;


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


	std::vector<float> vec;
	std::vector<float> _posShift;
	vec.push_back(499.584);
	vec.push_back(499.584);
	vec.push_back(0);
	// registerProcessorParameter( "PositionShift" ,
	// 			      "3 Vector to shift to have the right (0,0,0) position",
	// 			      _posShift,
	// 			      vec );
	posShift=CLHEP::Hep3Vector( 0.0,
								0.0,
								0.0 );


	float thr[] = {
		0.114 , 0.14 , 0.155714 , 0.171429 , 0.187143 , 0.202857 , 0.218571 , 0.234286 , 0.25 , 0.265714 , 0.281429 , 0.298571 , 0.314286 , 0.33 , 0.345714 , 0.361429 , 0.377143 , 0.392857 , 0.408571 , 0.424286 ,
		0.4 , 0.6125 , 0.825 , 1.0375 , 1.2375 , 1.45 , 1.6625 , 1.875 , 2.0875 , 2.3 , 2.5 , 2.7125 , 2.925 , 3.1375 , 3.35 , 3.5625 , 3.7625 , 3.975 , 4.1875 ,
		4.29448 , 5.33742 , 6.38037 , 7.48466 , 8.52761 , 9.57055 , 10.6135 , 11.6564 , 12.6994 , 13.8037 , 14.8466 , 15.8896 , 16.9325 , 17.9755 , 19.0184 , 20.0613 , 21.1656 , 22.2086 , 23.2515
	} ;
	std::vector<float> thresholdsVec(thr, thr + sizeof(thr) / sizeof(float) ) ;

	//		std::vector<double> thresholdsVec ;
	//		double t = 0.1 ;
	//		while ( t <= 26 )
	//		{
	//			thresholdsVec.push_back(t) ;
	//			t += 0.1 ;
	//		}

	std::sort(thresholdsVec.begin() , thresholdsVec.end() ) ;

//	thresholds = thresholdsVec ;

	registerProcessorParameter( "Thresholds" ,
								"Vector of thresholds",
								thresholdsFloat,
								thresholdsVec );


	AlgorithmRegistrationParameters() ;
}


void AnalogEfficiencyProcessor::AlgorithmRegistrationParameters()
{
	/*------------algorithm::Cluster------------*/
	registerProcessorParameter( "Cluster::MaxTransversalCellID" ,
								"Maximum difference between two hits cellID (0 and 1) to build a cluster",
								m_ClusterParameterSetting.maxTransversal,
								(int) 1 );

	registerProcessorParameter( "Cluster::MaxLongitudinalCellID" ,
								"Maximum difference between two hits cellID (2) to build a cluster",
								m_ClusterParameterSetting.maxLongitudinal,
								(int) 0 );

	registerProcessorParameter( "Cluster::UseDistanceInsteadCellID" ,
								"Boolean to know if clustering algorithms uses distance instead of cellID to cluster hits together",
								m_ClusterParameterSetting.useDistanceInsteadCellID,
								(bool) false );

	registerProcessorParameter( "Cluster::MaxTransversalDistance" ,
								"Maximum transversal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxTransversalDistance,
								(float) 11.0 );

	registerProcessorParameter( "Cluster::MaxLongitudinalDistance" ,
								"Maximum longitudinal distance (in mm) between two hits to gathered them in one cluster",
								m_ClusterParameterSetting.maxLongitudinalDistance,
								(float) 27.0 );

	/*------------algorithm::ClusteringHelper------------*/
	registerProcessorParameter( "ClusteringHelper::LongitudinalDistanceForIsolation" ,
								"Minimum longitudinal distance (in mm) between one hits and its neighbours to decide if it is isolated",
								m_ClusteringHelperParameterSetting.longitudinalDistance,
								(float) 100.0 );

	registerProcessorParameter( "ClusteringHelper::TransversalDistanceDistanceForIsolation" ,
								"Minimum transversal distance (in mm) between one hits and its neighbours to decide if it is isolated",
								m_ClusteringHelperParameterSetting.transversalDistance,
								(float) 200.0 );

	/*------------algorithm::Tracking-----------*/
	registerProcessorParameter( "Tracking::ChiSquareLimit" ,
								"Maximum value of tracking fit chi2 to construct a track",
								m_TrackingParameterSetting.chiSquareLimit,
								(float) 100.0 );

	registerProcessorParameter( "Tracking::MaxTransverseRatio" ,
								"Maximum value of transverse ratio to construct a track",
								m_TrackingParameterSetting.maxTransverseRatio,
								(float) 0.05 );

	registerProcessorParameter( "Tracking::CosThetaLimit" ,
								"Minimum value of cos(Theta) to accept the track",
								m_TrackingParameterSetting.cosThetaLimit,
								(float) 0.0 );

	registerProcessorParameter( "Tracking::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_TrackingParameterSetting.printDebug,
								(bool) false );

	/*------------algorithm::Efficiency-----------*/
	registerProcessorParameter( "Efficiency::MaxRadius" ,
								"Maximum distance parameter to find a hit to consider the layer as efficient",
								m_EfficiencyParameterSetting.maxRadius,
								(float) 25.0 );

	registerProcessorParameter( "Efficiency::SDHCALReadout" ,
								"Boolean to know if the detector used the semi digital readout",
								m_EfficiencyParameterSetting.semiDigitalReadout,
								(bool) true );

	registerProcessorParameter( "Efficiency::PrintDebug" ,
								"If true, Efficiency algorithm will print some debug information",
								m_EfficiencyParameterSetting.printDebug,
								(bool) false );

	m_EfficiencyParameterSetting.trackingParams=m_TrackingParameterSetting;

	/*------------algorithm::InteractionFinder-----------*/
	registerProcessorParameter( "InteractionFinder::MinSize" ,
								"Minimum cluster size for to define an interaction point",
								m_InteractionFinderParameterSetting.minSize,
								(int) 4 );

	registerProcessorParameter( "InteractionFinder::MaxRadius" ,
								"Maximum transversal distance to look for clusters",
								m_InteractionFinderParameterSetting.maxRadius,
								(float) 50.0 );

	registerProcessorParameter( "InteractionFinder::MaxDepth" ,
								"Maximum depth (number of layers) to look for clusters",
								m_InteractionFinderParameterSetting.maxDepth,
								(int) 4 );

	registerProcessorParameter( "InteractionFinder::MinNumberOfCluster" ,
								"Minimum number of found clusters (big enough) after the interaction point",
								m_InteractionFinderParameterSetting.minNumberOfCluster,
								(int) 3 );

	registerProcessorParameter( "InteractionFinder::UseAnalogEnergy" ,
								"Boolean to know if interaction finder algo should use cluster energy of cluster number of hits",
								m_InteractionFinderParameterSetting.useAnalogEnergy,
								(bool) false );

	registerProcessorParameter( "InteractionFinder::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_InteractionFinderParameterSetting.printDebug,
								(bool) false );
	/*------------caloobject::CaloGeom------------*/
	registerProcessorParameter( "Geometry::NLayers" ,
								"Number of layers",
								m_CaloGeomSetting.nLayers,
								(int) 48 );
	registerProcessorParameter( "Geometry::NPixelsPerLayer" ,
								"Number of pixels per layer (assume square geometry)",
								m_CaloGeomSetting.nPixelsPerLayer,
								(int) 96 );
	registerProcessorParameter( "Geometry::PixelSize" ,
								"Pixel size (assume square pixels)",
								m_CaloGeomSetting.pixelSize,
								(float) 10.408 );

	std::vector<float> vec;
	vec.push_back(0.0);
	vec.push_back(1000.0);
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
		std::cout << "WARNING : Wrong number of values in paramater Geometry::DetectorTransverseSize => will use default value -500.0, +500.0" << std::endl;
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
								(int) 96 );

	registerProcessorParameter( "AsicKeyFinder::NPadY" ,
								"Number of pads in x direction per layer",
								m_AsicKeyFinderParameterSetting.nPadY,
								(int) 96 );

	registerProcessorParameter( "AsicKeyFinder::AsicNPad" ,
								"number of pads in x or y direction per asic (assuming a square)",
								m_AsicKeyFinderParameterSetting.asicNPad,
								(int) 8 );

	registerProcessorParameter( "AsicKeyFinder::LayerGap" ,
								"Gap size (in mm) between 2 layers",
								m_AsicKeyFinderParameterSetting.layerGap,
								(float) 26.131 );

	registerProcessorParameter( "AsicKeyFinder::PadSize" ,
								"Size of one pad in mm",
								m_AsicKeyFinderParameterSetting.padSize,
								(float) 10.408 );

	registerProcessorParameter( "AsicKeyFinder::PrintDebug" ,
								"Boolean to know if debug if printed",
								m_AsicKeyFinderParameterSetting.printDebug,
								(bool) false );
}

void AnalogEfficiencyProcessor::init()
{
	std::sort(thresholdsFloat.begin() , thresholdsFloat.end() ) ;

	thresholds.clear() ;
	for ( unsigned int i = 0 ; i < thresholdsFloat.size() ; ++i )
		thresholds.push_back( thresholdsFloat.at(i) ) ;


	printParameters() ;

	file = new TFile(outputRootName.c_str() , "RECREATE") ;
	file->cd() ;
	TDirectory* dir = file->mkdir("Graphs") ;
	dir->cd() ;


	_nRun = 0 ;
	_nEvt = 0 ;
	_goodTrackCounter = 0 ;

	/*--------------------Algorithms initialisation--------------------*/
	algo_Cluster=new algorithm::Cluster();
	algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting);

	algo_ClusteringHelper=new algorithm::ClusteringHelper();
	algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting);

	algo_Tracking=new algorithm::Tracking();
	algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting);

	algo_InteractionFinder=new algorithm::InteractionFinder();
	algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderParameterSetting);

	m_EfficiencyParameterSetting.geometry=m_CaloGeomSetting;
	algo_Efficiency=new algorithm::Efficiency();
	algo_Efficiency->SetEfficiencyParameterSetting(m_EfficiencyParameterSetting);

	algo_AsicKeyFinder=new algorithm::AsicKeyFinder();
	algo_AsicKeyFinder->SetAsicKeyFinderParameterSetting(m_AsicKeyFinderParameterSetting);

	for ( int k = 0 ; k < _nActiveLayers ; k++ )
	{
		caloobject::Layer* aLayer = new caloobject::SDHCALLayer(k , _difList.at(3*k+2) , _difList.at(3*k+1) , _difList.at(k)) ;
		aLayer->setPosition( CLHEP::Hep3Vector(10.408 , 10.408 , (k+1)*m_AsicKeyFinderParameterSetting.layerGap) ) ;
		aLayer->buildAsics() ;
		aLayer->setThresholds(thresholds) ;
		layers.push_back(aLayer) ;
	}
}

void AnalogEfficiencyProcessor::processRunHeader( LCRunHeader* run )
{
	_nRun++ ;
	_nEvt = 0 ;
}

void AnalogEfficiencyProcessor::DoTracking()
{
	std::vector<caloobject::CaloCluster2D*> clusters ;
	for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it = hitMap.begin() ; it != hitMap.end() ; ++it)
		algo_Cluster->Run(it->second , clusters) ;

	std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer) ;
	for(std::vector<caloobject::CaloCluster2D*>::iterator it = clusters.begin() ; it != clusters.end() ; ++it)
	{
		if ( algo_ClusteringHelper->IsIsolatedCluster(*it,clusters) )
		{
			delete *it ;
			clusters.erase(it) ;
			it-- ;
		}
	}
	caloobject::CaloTrack* track = NULL ;
	algo_Tracking->Run(clusters,track) ;

	if( NULL != track )
	{
		_goodTrackCounter++ ;
		algo_InteractionFinder->Run(clusters , track->getTrackParameters()) ;

		if( algo_InteractionFinder->FindInteraction() == false )
			LayerProperties(clusters) ;
	}
	file->cd() ;

	delete track ;
	for ( std::vector<caloobject::CaloCluster2D*>::iterator it = clusters.begin() ; it != clusters.end() ; ++it )
		delete (*it) ;

}

void AnalogEfficiencyProcessor::LayerProperties(std::vector<caloobject::CaloCluster2D*> &clusters)
{
	int trackBegin = (*clusters.begin())->getLayerID() ;
	int trackEnd = (*(clusters.end()-1))->getLayerID() ;
	if ( trackBegin == 1 )
		trackBegin = 0 ;

	if ( trackEnd == _nActiveLayers-2 )
		trackEnd = _nActiveLayers-1 ;


	for( int K = trackBegin ; K <= trackEnd ; K++ )
	{
		algo_Efficiency->Run(layers.at(K) , clusters) ;

		if ( algo_Efficiency->isTrack() )
		{
			layers.at(K)->update( algo_Efficiency->getExpectedPosition() , algo_Efficiency->getGoodCluster() ) ;
			//			trackPositionHist->Fill( algo_Efficiency->getExpectedPosition().x() , algo_Efficiency->getExpectedPosition().y() ) ;
		}
	}
}

void AnalogEfficiencyProcessor::processEvent( LCEvent * evt )
{
	//
	// * Reading HCAL Collections of CalorimeterHits*
	//
	//std::string initString;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");

	for (unsigned int i(0); i < _hcalCollections.size(); ++i)
	{
		std::string colName = _hcalCollections[i] ;
		try
		{
			col = evt->getCollection( _hcalCollections[i].c_str() ) ;
			//initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
			numElements = col->getNumberOfElements();
			//      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);

			for (int j=0; j < numElements; ++j)
			{
				CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
				CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
				int cellID[] = { IDdecoder(hit)["I"] , IDdecoder(hit)["J"] , IDdecoder(hit)["K-1"] } ;

				if ( cellID[2] > _nActiveLayers )
					continue ;

				if ( cellID[0] < 1 || cellID[0] > 96 || cellID[1] < 1 || cellID[1] > 96)
					continue ;

				if (hit->getEnergy() < thresholds.at(0) )
					continue ;

				caloobject::CaloHit* aHit = new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime() , posShift) ;
				hitMap[cellID[2]].push_back(aHit) ;

			}

			DoTracking() ;
			clearVec() ;
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << "Exeption " << std::endl;
		}
	}
	_nEvt ++ ;
	std::cout << "Event processed : " << _nEvt << std::endl;
}

void AnalogEfficiencyProcessor::clearVec()
{
	for( std::map<int,std::vector<caloobject::CaloHit*> >::iterator it = hitMap.begin() ; it != hitMap.end() ; ++it )
		for( std::vector<caloobject::CaloHit*>::iterator jt = (it->second).begin() ; jt != (it->second).end() ; ++jt )
			delete *jt ;

	hitMap.clear() ;
}


void AnalogEfficiencyProcessor::check( LCEvent * evt )
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AnalogEfficiencyProcessor::end()
{
	file->cd() ;
	TDirectory* graphDir = file->GetDirectory("Graphs") ;
	graphDir->cd() ;

	int nOkAsics = 0 ;
	int nTracksTotal = 0 ;
	std::vector<double> effTotal( thresholds.size() , 0 ) ;
	std::vector<double> effErrTotal( thresholds.size() , 0 ) ;

	for ( std::vector<caloobject::Layer*>::const_iterator it = layers.begin() ; it != layers.end() ; it++ )
	{
		if ( (*it)->getNTrack() == 0 )
			continue ;
		if ( (*it)->getEfficiency() < 0.1 )
			continue ;


		std::vector<double> effLayer = (*it)->getEfficiencies() ;
		std::vector<double> effErrLayer = (*it)->getEfficienciesError() ;
		TGraphErrors* graph = new TGraphErrors( static_cast<int>( thresholds.size() ) , &thresholds[0] , &effLayer[0] , NULL , &effErrLayer[0] ) ;

		int layerID = (*it)->getID() ;

		std::stringstream dirName ; dirName << "Layer" << layerID ;
		TDirectory* layerDir = graphDir->GetDirectory( dirName.str().c_str() ) ;
		if ( !layerDir )
			layerDir = graphDir->mkdir( dirName.str().c_str() ) ;

		layerDir->cd() ;
		std::stringstream graphName ;
		graphName << "Layer" << layerID ;
		graph->SetMarkerStyle(20) ;
		graph->Write(graphName.str().c_str() ) ;

		AsicMap asics = (*it)->getAsics() ;

		for ( AsicMap::const_iterator asicIt = asics.begin() ; asicIt != asics.end() ; asicIt++ )
		{
			if ( asicIt->second->getNTrack() == 0 )
				continue ;
			if ( asicIt->second->getEfficiency() < 0.1 )
				continue ;

			nOkAsics++ ;
			std::vector<double> effAsic = asicIt->second->getEfficiencies() ;
			std::vector<double> effErrAsic = asicIt->second->getEfficienciesError() ;
			TGraphErrors* graph = new TGraphErrors( static_cast<int>( thresholds.size() ) , &thresholds[0] , &effAsic[0] , NULL , &effErrAsic[0] ) ;

			for ( unsigned int i = 0 ; i < thresholds.size() ; ++i )
				effTotal.at(i) += effAsic.at(i) ;
			nTracksTotal += asicIt->second->getNTrack() ;

			int difID = asicIt->second->getDifID() ;
			int asicID = asicIt->second->getID() ;

			std::stringstream graphName ;
			graphName << difID << "," << asicID ;

			graph->SetMarkerStyle(20) ;
			graph->Write(graphName.str().c_str() ) ;
		}

	}

	for ( unsigned int i = 0 ; i < thresholds.size() ; ++i )
	{
		effTotal.at(i) /= nOkAsics ;

		double error = std::sqrt(effTotal.at(i)*(1-effTotal.at(i))/nTracksTotal) ;

		if ( error < std::numeric_limits<double>::epsilon() )
		{
			if ( effTotal.at(i) < std::numeric_limits<double>::epsilon() )
			{
				double falseEff = 1.0/nTracksTotal ;
				error = std::sqrt(falseEff*(1-falseEff)/nTracksTotal) ;
			}
			else
			{
				double falseEff = (nTracksTotal-1.0)/nTracksTotal ;
				error = std::sqrt(falseEff*(1-falseEff)/nTracksTotal) ;
			}
		}
		effErrTotal.at(i) = error ;
	}

	TGraphErrors* globalGraph = new TGraphErrors( static_cast<int>( thresholds.size() ) , &thresholds[0] , &effTotal[0] , NULL , &effErrTotal[0] ) ;

	graphDir->cd() ;
	std::stringstream graphName ;
	graphName << "Global" ;

	globalGraph->SetMarkerStyle(20) ;
	globalGraph->Write(graphName.str().c_str() ) ;


	delete algo_Cluster ;
	delete algo_ClusteringHelper ;
	delete algo_Tracking ;
	delete algo_InteractionFinder ;
	delete algo_Efficiency ;
	delete algo_AsicKeyFinder ;

	for(std::vector<caloobject::Layer*>::iterator it = layers.begin(); it!=layers.end(); ++it)
		delete (*it) ;
	layers.clear() ;


	file->Write() ;
	file->Close() ;

	std::cout << "_goodTrackCounter " << _goodTrackCounter << std::endl ;
}
