#ifndef AnalysisProcessorILD_h
#define AnalysisProcessorILD_h

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
#include "Algorithm/Density.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace lcio ;
using namespace marlin ;

class AnalysisProcessorILD : public Processor
{

	public :

		virtual Processor* newProcessor() { return new AnalysisProcessorILD ; }


		AnalysisProcessorILD() ;

		/** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
		virtual void init() ;

		/** Called for every run.
   */
		virtual void processRunHeader(LCRunHeader* ) ;

		/** Called for every event - the working horse.
   */
		virtual void processEvent(LCEvent* evt) ;


		virtual void check(LCEvent* ) ;


		/** Called after data processing for clean up.
   */
		virtual void end() ;

		void AlgorithmRegistrationParameters() ;

		void clearVec() ;

		void findEventTime(LCEvent* evt , LCCollection* _col) ;

		void findSpillEventTime(LCEvent* evt , LCCollection* _col) ;

		double getFirst5LayersRMS() ;
		int getNInteractingLayer() ;

		AnalysisProcessorILD(const AnalysisProcessorILD &toCopy) = delete ;
		void operator=(const AnalysisProcessorILD &toCopy) = delete ;


	protected:

		int _nRun = 0 ;
		int _nEvt = 0 ;
		/** Input collection name.
   */
		std::vector<std::string> _hcalCollections ;

	private:
		std::map<int,HitVec> hitMap ;
		std::vector<caloobject::CaloCluster2D*> clusterVec ;

		/*--------------------Global parameters--------------------*/
		int _nActiveLayers = 0 ;
		int numElements = 0 ;
		LCCollection* col = nullptr ;
		std::vector<float> edges ; //vector to recover geometry parameters
		CLHEP::Hep3Vector posShift ;
		/*------------------------------------------------------------------------------*/

		/*--------------------Algorithms list to initialise--------------------*/
		algorithm::Clustering* algo_Cluster;
		algorithm::ClusteringHelper* algo_ClusteringHelper;
		algorithm::Tracking* algo_Tracking;
		algorithm::Hough* algo_Hough ;
		algorithm::InteractionFinder* algo_InteractionFinder;
		algorithm::Density* algo_density ;

		/*------------------------------------------------------------------------------*/

		/*--------------------Algorithms setting parameter structure--------------------*/
		algorithm::clusterParameterSetting m_ClusterParameterSetting;
		algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting;
		algorithm::TrackingParameterSetting m_TrackingParameterSetting;
		algorithm::HoughParameterSetting m_HoughParameterSetting ;
		algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting;

		/*------------------------------------------------------------------------------*/

		/*--------------------CaloObject setting parameter structure--------------------*/
		caloobject::GeomParameterSetting m_CaloGeomSetting;
		/*------------------------------------------------------------------------------*/

		algorithm::ShowerAnalyserParameterSetting m_ShowerAnalyserParameterSetting ;
		algorithm::ShowerAnalyser* algo_ShowerAnalyser ;

		/*--------------------CaloObject list to initialise--------------------*/


		/*---------------------------------------------------------------------*/

		std::vector<float> thresholdsFloat ;
		std::vector<double> thresholds ;

		double _timeCut = 0 ;
		unsigned long long _prevBCID = 0 ;
		unsigned long long _bcidRef = 0 ;
		int firstShowerInSpill = 0 ;
		bool firstSpillEvtFound = false ;

		int _timeDif_minus_bif = 0 ;

		/*--------------------Root output object--------------------*/
		std::string outputRootName = "" ;
		TFile *file = nullptr ;
		TTree* tree = nullptr ;

		double computingTime = 0 ;

		int eventNumber = 0 ;
		unsigned long long evtTime  = 0 ;
		unsigned long long spillEvtTime = 0 ;

		bool cerenkovTag = false ;

		int runNumber = 0 ;
		float energy = 0 ;

		int nHit = 0 ;
		int nHit1 = 0 ;
		int nHit2 = 0 ;
		int nHit3 = 0 ;
		int nHough = 0 ;
		int nHough1 = 0 ;
		int nHough2 = 0 ;
		int nHough3 = 0 ;
		int nLayer = 0 ;
		int nInteractingLayer = 0 ;
		int nCluster = 0 ;
		int nMipCluster = 0 ;
		int nTrack = 0 ;

		double first5LayersRMS = 0 ;

		std::vector<int> tracksClusterSize ;
		std::vector<int> tracksClusterNumber ;

		double begin = 0 ;
		double _end = 0 ;

		double reconstructedCosTheta = 0 ;

		float transverseRatio = 0 ;

		double density = 0 ;

		float thrust[4] = {0,0,0,0} ;

		bool neutral = 0 ;

		std::vector<double> longiProfile ;
		std::vector<double> radiProfile ;

		std::vector<int> iVec ;
		std::vector<int> jVec ;
		std::vector<int> kVec ;
		std::vector<int> thrVec ;

		double emFraction = 0 ;

} ;


#endif //AnalysisProcessorILD_h
