#ifndef AnalysisProcessor_h
#define AnalysisProcessor_h

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>
#include <limits>

#include "CaloObject/CaloGeom.h"
#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloCluster.h"
#include "CaloObject/Asic.h"
#include "CaloObject/Shower.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/Hough.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Efficiency.h"
#include "Algorithm/AsicKeyFinder.h"
#include "Algorithm/ShowerAnalyser.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace lcio ;
using namespace marlin ;

class AnalysisProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new AnalysisProcessor ; }


		AnalysisProcessor() ;

		/** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
		virtual void init() ;

		/** Called for every run.
   */
		virtual void processRunHeader( LCRunHeader* run ) ;

		/** Called for every event - the working horse.
   */
		virtual void processEvent( LCEvent * evt ) ;


		virtual void check( LCEvent * evt ) ;


		/** Called after data processing for clean up.
   */
		virtual void end() ;

		void AlgorithmRegistrationParameters();

		void clearVec();

		void findEventTime(LCEvent* evt,LCCollection* col) ;

		void findSpillEventTime(LCEvent* evt,LCCollection* col) ;

		double getFirst5LayersRMS() ;
		int getNInteractingLayer() ;


	protected:

		int _nRun ;
		int _nEvt ;
		/** Input collection name.
   */
		std::vector<std::string> _hcalCollections;

	private:
		std::map<int,std::vector<caloobject::CaloHit*> > hitMap ;
		std::vector<caloobject::CaloCluster2D*> clusterVec ;

		/*--------------------Global parameters--------------------*/
		int _nActiveLayers;
		int numElements;
		LCCollection * col;
		int _nAsicX;
		int _nAsicY;
		std::vector<int> _difList;
		std::vector<float> edges; //vector to recover geometry parameters
		CLHEP::Hep3Vector posShift;
		/*------------------------------------------------------------------------------*/

		/*--------------------Algorithms list to initialise--------------------*/
		algorithm::Cluster *algo_Cluster;
		algorithm::ClusteringHelper *algo_ClusteringHelper;
		algorithm::Tracking *algo_Tracking;
		algorithm::Hough* algo_Hough ;
		algorithm::InteractionFinder *algo_InteractionFinder;
		algorithm::Efficiency *algo_Efficiency;
		algorithm::AsicKeyFinder *algo_AsicKeyFinder;

		/*------------------------------------------------------------------------------*/

		/*--------------------Algorithms setting parameter structure--------------------*/
		algorithm::clusterParameterSetting m_ClusterParameterSetting;
		algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting;
		algorithm::TrackingParameterSetting m_TrackingParameterSetting;
		algorithm::HoughParameterSetting m_HoughParameterSetting ;
		algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting;
		algorithm::EfficiencyParameterSetting m_EfficiencyParameterSetting;
		algorithm::AsicKeyFinderParameterSetting m_AsicKeyFinderParameterSetting;

		/*------------------------------------------------------------------------------*/

		/*--------------------CaloObject setting parameter structure--------------------*/
		caloobject::GeomParameterSetting m_CaloGeomSetting;
		/*------------------------------------------------------------------------------*/

		algorithm::ShowerAnalyserParameterSetting m_ShowerAnalyserParameterSetting ;
		algorithm::ShowerAnalyser* algo_ShowerAnalyser ;

		/*--------------------CaloObject list to initialise--------------------*/


		/*---------------------------------------------------------------------*/

		double _timeCut ;
		unsigned long long _prevBCID ;
		unsigned long long _bcidRef ;
		int firstShowerInSpill ;
		bool firstSpillEvtFound ;

		int _timeDif_minus_bif ;

		/*--------------------Root output object--------------------*/
		std::string outputRootName ;
		TFile *file ;
		TTree* tree ;

		int eventNumber ;
		unsigned long long evtTime ;
		unsigned long long spillEvtTime ;

		bool cerenkovTag ;

		int nHit ;
		int nHit1 ;
		int nHit2 ;
		int nHit3 ;
		int nLayer ;
		int nInteractingLayer ;
		int nCluster ;
		int nMipCluster ;
		int nTrack ;

		double first5LayersRMS ;

		std::vector<int> tracksClusterSize ;
		std::vector<int> tracksClusterNumber ;

		double begin ;

		double reconstructedCosTheta ;

		float transverseRatio ;

		float thrust[4] ;

		bool neutral ;

		std::vector<double> longiProfile ;
		std::vector<double> radiProfile ;

		std::vector<int> iVec ;
		std::vector<int> jVec ;
		std::vector<int> kVec ;
		std::vector<int> thrVec ;


} ;


#endif //AnalysisProcessor_h
