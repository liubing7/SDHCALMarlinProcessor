#ifndef LcioToRootProcessor_h
#define LcioToRootProcessor_h

#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>

#include <string>

#include <TFile.h>
#include <TTree.h>

class LcioToRootProcessor : public marlin::Processor
{

	public :

		virtual marlin::Processor* newProcessor() { return new LcioToRootProcessor ; }


		LcioToRootProcessor() ;

		virtual void init() ;

		virtual void processEvent(LCEvent* evt) ;

		virtual void end() ;

		void clear() ;

		LcioToRootProcessor(const LcioToRootProcessor &toCopy) = delete ;
		void operator=(const LcioToRootProcessor &toCopy) = delete ;


	protected :

		std::vector<std::string> inputHitCollection {} ;
		std::vector<std::string> inputRelCollection {} ;
		std::string inputMCParticleCollection = "" ;

		std::string outputRootFileName = "" ;

		TFile* file = nullptr ;
		TTree* tree = nullptr ;

		int particle1PDG = 0 ;
		int particle2PDG = 0 ;
		int particle3PDG = 0 ;

		std::vector<int> iVec {} ;
		std::vector<int> jVec {} ;
		std::vector<int> kVec {} ;
		std::vector<float> xVec {} ;
		std::vector<float> yVec {} ;
		std::vector<float> zVec {} ;
		std::vector<int> thrVec {} ;
		std::vector<float> timeVec {} ;
		std::vector<float> particle1Vec {} ;
		std::vector<float> particle2Vec {} ;
		std::vector<float> particle3Vec {} ;

		int eventNumber = 0 ;
} ;

#endif //LcioToRootProcessor_h
