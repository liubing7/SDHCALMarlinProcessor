#ifndef EfficiencyVsAngleProcessor_h
#define EfficiencyVsAngleProcessor_h

#include <marlin/Processor.h>
#include <lcio.h>
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>
#include <limits>

#include "CaloObject/CaloGeom.h"
#include "CaloObject/CaloHit.h"
#include "CaloObject/Asic.h"

#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Efficiency.h"
#include "Algorithm/AsicKeyFinder.h"

#include "CaloObject/Layer.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace lcio ;
using namespace marlin ;


class EfficiencyVsAngleProcessor : public Processor {

	public :

		virtual Processor*  newProcessor() { return new EfficiencyVsAngleProcessor ; }


		EfficiencyVsAngleProcessor() ;

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
		void LayerProperties(std::vector<caloobject::CaloCluster2D*> &clusters) ;
		void clearVec() ;
		void DoTracking() ;
		inline int findDifID(int key) { return _difList.at( key/1000*3 + 2-key%1000/12/4 ) ; }


	protected :

		int _nRun ;
		int _nEvt ;
		int _goodTrackCounter ;
		/** Input collection name.
   */
		std::vector<std::string> _hcalCollections ;

	private:
		std::map<int,std::vector<caloobject::CaloHit*> > hitMap ;

		/*--------------------Global parameters--------------------*/
		int _nActiveLayers ;
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
		algorithm::InteractionFinder *algo_InteractionFinder;
		algorithm::Efficiency* algo_Efficiency ;
		algorithm::AsicKeyFinder *algo_AsicKeyFinder;
		/*------------------------------------------------------------------------------*/

		/*--------------------Algorithms setting parameter structure--------------------*/
		algorithm::clusterParameterSetting m_ClusterParameterSetting;
		algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting;
		algorithm::TrackingParameterSetting m_TrackingParameterSetting;
		algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting;
		algorithm::EfficiencyParameterSetting m_EfficiencyParameterSetting;
		algorithm::AsicKeyFinderParameterSetting m_AsicKeyFinderParameterSetting;
		/*------------------------------------------------------------------------------*/

		/*--------------------CaloObject setting parameter structure--------------------*/
		caloobject::GeomParameterSetting m_CaloGeomSetting;
		/*------------------------------------------------------------------------------*/

		/*--------------------CaloObject list to initialise--------------------*/
		std::vector<caloobject::Layer*> layers ;
		/*---------------------------------------------------------------------*/

		std::vector<float> thresholdsFloat ;
		std::vector<double> thresholds ;

		/*--------------------Root output object--------------------*/
		std::string outputRootName ;
		TFile* file ;
		TTree* tree ;


		TH2D* mulAngleHist ;
		TH2D* eff1AngleHist ;
		TH2D* eff2AngleHist ;
		TH2D* eff3AngleHist ;

		std::vector<int> nTracksAngleVec ;
		std::vector<double> mulAngleVec ;
		std::vector<double> eff1AngleVec ;
		std::vector<double> eff2AngleVec ;
		std::vector<double> eff3AngleVec ;

		TH2D* trackPositionHist ;

		int layerID ;
		int difID ;
		int asicID ;
		int padID ;
		double multiplicity ;
		double multiplicityError ;

		std::vector<double> efficiencies ;
		std::vector<double> efficienciesError ;

		int nTracks ;

		std::vector<double> position ;
} ;

#endif //EfficiencyVsAngleProcessor_h
