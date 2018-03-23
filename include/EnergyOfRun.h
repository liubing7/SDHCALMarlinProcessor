#ifndef EnergyOfRun_h
#define EnergyOfRun_h

#include <iostream>
#include <map>


float energyOfRun(int nRun)
{
	std::map<int , float> runMap ;

	//electrons SPS Aout2012
	runMap.insert( std::make_pair(715724 , 10) ) ;
	runMap.insert( std::make_pair(715725 , 10) ) ;
	runMap.insert( std::make_pair(715715 , 20) ) ;
	runMap.insert( std::make_pair(715714 , 30) ) ;
	runMap.insert( std::make_pair(715713 , 40) ) ;
	runMap.insert( std::make_pair(715716 , 50) ) ;

	//electrons SPS Avril2015
	runMap.insert( std::make_pair(728348 , 10) ) ;
	runMap.insert( std::make_pair(728350 , 15) ) ;
	runMap.insert( std::make_pair(728203 , 20) ) ;
	runMap.insert( std::make_pair(728352 , 25) ) ;
	runMap.insert( std::make_pair(728354 , 30) ) ;
	runMap.insert( std::make_pair(728357 , 40) ) ;
	runMap.insert( std::make_pair(728359 , 50) ) ;
	runMap.insert( std::make_pair(728360 , 85) ) ;

	//pions SPS Nov2012 H2
	runMap.insert( std::make_pair(716321 , 10) ) ;
	runMap.insert( std::make_pair(716315 , 20) ) ;
	runMap.insert( std::make_pair(716308 , 30) ) ;
	runMap.insert( std::make_pair(716307 , 40) ) ;
	runMap.insert( std::make_pair(716305 , 50) ) ;
	runMap.insert( std::make_pair(716298 , 60) ) ;
	runMap.insert( std::make_pair(716290 , 70) ) ;
	runMap.insert( std::make_pair(716319 , 80) ) ;

	//pions SPS Oct2015
	//163
	runMap.insert( std::make_pair(730716 , 10) ) ;
	runMap.insert( std::make_pair(730656 , 20) ) ;
	runMap.insert( std::make_pair(730634 , 30) ) ;
	runMap.insert( std::make_pair(730648 , 40) ) ;
	runMap.insert( std::make_pair(730657 , 40) ) ;
	runMap.insert( std::make_pair(730651 , 50) ) ;
	runMap.insert( std::make_pair(730655 , 60) ) ;
	runMap.insert( std::make_pair(730659 , 70) ) ;
	runMap.insert( std::make_pair(730668 , 70) ) ;
	runMap.insert( std::make_pair(730677 , 80) ) ;
	runMap.insert( std::make_pair(730661 , 90) ) ;
	runMap.insert( std::make_pair(730676 , 90) ) ;

	//208  (perlayer)
	runMap.insert( std::make_pair(730823 , 10) ) ;
	runMap.insert( std::make_pair(730678 , 20) ) ;
	runMap.insert( std::make_pair(730816 , 30) ) ;
	runMap.insert( std::make_pair(730819 , 40) ) ;
	runMap.insert( std::make_pair(730821 , 50) ) ;
	runMap.insert( std::make_pair(730824 , 60) ) ;
	runMap.insert( std::make_pair(730842 , 70) ) ;
	runMap.insert( std::make_pair(730846 , 80) ) ;

	//214  (perasic)
	runMap.insert( std::make_pair(730903 , 10) ) ;
	runMap.insert( std::make_pair(730888 , 20) ) ;
	runMap.insert( std::make_pair(730886 , 30) ) ;
	runMap.insert( std::make_pair(730882 , 40) ) ;
	runMap.insert( std::make_pair(730861 , 50) ) ;
	runMap.insert( std::make_pair(730858 , 60) ) ;
	runMap.insert( std::make_pair(730851 , 70) ) ;
	runMap.insert( std::make_pair(730847 , 80) ) ;



	//pions SPS Oct2016
	// 1 pC 2thr
	runMap.insert( std::make_pair(733660 , 20) ) ;
	runMap.insert( std::make_pair(733665 , 30) ) ;
	runMap.insert( std::make_pair(733683 , 40) ) ;
	runMap.insert( std::make_pair(733686 , 50) ) ;
	runMap.insert( std::make_pair(733689 , 60) ) ;
	runMap.insert( std::make_pair(733693 , 70) ) ;
	runMap.insert( std::make_pair(733696 , 80) ) ;


	// 5 pC 2thr
	runMap.insert( std::make_pair(733756 , 20) ) ;
	runMap.insert( std::make_pair(733750 , 30) ) ;
	runMap.insert( std::make_pair(733724 , 40) ) ;
	runMap.insert( std::make_pair(733728 , 50) ) ;
	runMap.insert( std::make_pair(733742 , 60) ) ;
	runMap.insert( std::make_pair(733743 , 70) ) ;
	runMap.insert( std::make_pair(733754 , 80) ) ;


	//pions SPS Sep2017
	runMap.insert( std::make_pair(736525 , 20) ) ;
	runMap.insert( std::make_pair(736523 , 30) ) ;
	runMap.insert( std::make_pair(736511 , 40) ) ;
	runMap.insert( std::make_pair(736517 , 50) ) ;
	runMap.insert( std::make_pair(736519 , 60) ) ;
	runMap.insert( std::make_pair(736520 , 70) ) ;
	runMap.insert( std::make_pair(736522 , 80) ) ;


	std::map<int , float>::const_iterator it = runMap.find(nRun) ;

	if ( it == runMap.end() )
	{
		std::cout << "Energy of run " << nRun << " not found : return 0" << std::endl ;
		return 0 ;
	}
	else
	{
		std::cout << "Energy of run " << nRun << " found : " << it->second << " GeV" << std::endl ;
		return it->second ;
	}
}


#endif //EnergyOfRun_h
