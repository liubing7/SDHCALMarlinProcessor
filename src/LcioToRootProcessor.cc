#include "LcioToRootProcessor.h"

#include <cassert>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/MCParticle.h>

LcioToRootProcessor aLcioToRootProcessor ;

LcioToRootProcessor::LcioToRootProcessor()
	: marlin::Processor("LcioToRootProcessor")
{
	registerProcessorParameter( "RootFileName" ,
								"output ROOT file name"  ,
								outputRootFileName  ,
								std::string{"rootFile.root"} ) ;

	std::vector<std::string> _inputHitCollection { "HCALEndcap" } ;
	std::vector<std::string> _inputRelCollection { "RelationParticleToHit" } ;
	std::string _inputMCParticleCollection = " primaryParticles" ;

	registerInputCollections( LCIO::CALORIMETERHIT,
							  "inputHitCollection" ,
							  "input CalorimeterHit collections" ,
							  inputHitCollection  ,
							  _inputHitCollection) ;

	registerInputCollections( LCIO::LCRELATION ,
							  "inputRelCollection" ,
							  "input CalorimeterHit to MCParticle Relation collections" ,
							  inputRelCollection  ,
							  _inputRelCollection) ;

	registerInputCollection( LCIO::MCPARTICLE ,
							 "inputMCParticleCollection" ,
							 "MCParticle collection name"  ,
							 inputMCParticleCollection  ,
							 _inputMCParticleCollection) ;
}


void LcioToRootProcessor::init()
{
	assert ( inputHitCollection.size() == inputRelCollection.size() ) ;
	assert ( !inputMCParticleCollection.empty() ) ;

	file = new TFile(outputRootFileName.c_str() , "RECREATE") ;
	tree = new TTree("tree","tree") ;

	tree->Branch("eventNumber" , &eventNumber) ;

	tree->Branch("particle1PDG" , &particle1PDG) ;
	tree->Branch("particle2PDG" , &particle2PDG) ;
	tree->Branch("particle2PDG" , &particle3PDG) ;

	tree->Branch("I" , &iVec) ;
	tree->Branch("J" , &jVec) ;
	tree->Branch("K" , &kVec) ;
	tree->Branch("x" , &xVec) ;
	tree->Branch("y" , &yVec) ;
	tree->Branch("z" , &zVec) ;

	tree->Branch("thr" , &thrVec) ;
	tree->Branch("time" , &timeVec) ;

	tree->Branch("particle1Tag" , &particle1Vec) ;
	tree->Branch("particle2Tag" , &particle2Vec) ;
	tree->Branch("particle3Tag" , &particle3Vec) ;
}

void LcioToRootProcessor::clear()
{
	iVec.clear() ;
	jVec.clear() ;
	kVec.clear() ;
	xVec.clear() ;
	yVec.clear() ;
	zVec.clear() ;
	thrVec.clear() ;
	timeVec.clear() ;
	particle1Vec.clear() ;
	particle2Vec.clear() ;
	particle3Vec.clear() ;
}

void LcioToRootProcessor::processEvent(EVENT::LCEvent* evt)
{
	clear() ;



	LCCollection* mcParticleCol = evt->getCollection( inputMCParticleCollection ) ;

	std::map<MCParticle* , unsigned int> particleVec {} ;
	for ( int iCol = 0 ; iCol < mcParticleCol->getNumberOfElements() ; ++iCol )
	{
		MCParticle* particle = dynamic_cast<MCParticle*>( mcParticleCol->getElementAt(iCol) ) ;
		particleVec.insert( {particle, iCol} ) ;

		if ( iCol == 0 )
			particle1PDG = particle->getPDG() ;
		else if ( iCol == 1 )
			particle2PDG = particle->getPDG() ;
		else if ( iCol == 2 )
			particle3PDG = particle->getPDG() ;
	}

	assert ( !particleVec.empty() ) ;

	for ( unsigned int iCol = 0 ; iCol < inputHitCollection.size() ; ++iCol )
	{
		LCCollection* hitCol = evt->getCollection( inputHitCollection[iCol] ) ;
		LCCollection* relCol = evt->getCollection( inputRelCollection[iCol] ) ;

		LCRelationNavigator navi(relCol) ;

		UTIL::CellIDDecoder<EVENT::CalorimeterHit> decoder(hitCol) ;

		for ( int iHit = 0 ; iHit < hitCol->getNumberOfElements() ; ++iHit )
		{
			CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( hitCol->getElementAt(iHit) ) ;

			int i = static_cast<int>( decoder(hit)["I"] ) ;
			int j = static_cast<int>( decoder(hit)["J"] ) ;
			int k = static_cast<int>( decoder(hit)["K-1"] ) ;

			if ( i < 1 || i > 96 || j < 1 || j > 96 || k < 0 || k > 47 )
				continue ;

			iVec.push_back( i ) ;
			jVec.push_back( j ) ;
			kVec.push_back( k ) ;

			xVec.push_back( hit->getPosition()[0] ) ;
			yVec.push_back( hit->getPosition()[1] ) ;
			zVec.push_back( hit->getPosition()[2] ) ;

			thrVec.push_back( static_cast<int>(hit->getEnergy() ) ) ;
			timeVec.push_back( hit->getTime() ) ;


			std::vector<float> weightVec( particleVec.size() , 0.f ) ;

			for ( unsigned int iRelPart = 0 ; iRelPart < navi.getRelatedFromObjects( hit ).size() ; ++iRelPart )
			{
				MCParticle* particle = dynamic_cast<MCParticle*>( navi.getRelatedFromObjects(hit)[iRelPart] ) ;
				float weight = navi.getRelatedFromWeights(hit)[iRelPart] ;

				const auto it = particleVec.find(particle) ;

				assert ( it != particleVec.end() ) ;

				unsigned int index = it->second ;

				weightVec.at(index) += weight ;
			}

			particle1Vec.push_back( weightVec.at(0) ) ;

			if ( weightVec.size() > 1 )
				particle2Vec.push_back( weightVec.at(1) ) ;
			if ( weightVec.size() > 2 )
				particle3Vec.push_back( weightVec.at(2) ) ;
		}
	}

	tree->Fill() ;
	eventNumber++ ;
}

void LcioToRootProcessor::end()
{
	tree->Write("tree") ;
	file->Close() ;
}

