//============================================================================
// Name        : IonTimes.cpp
// Author      : Sz5
// Version     : 0.2
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

///GSL Headers
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_specfunc.h>

///Standard C++ Headers
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <vector>

///Program Specific Headers
#include "DefaultAction.hpp"
#include "I_Action.hpp"
#include "ResultGrid.hpp"
#include "SaddlePointSolver.hpp"
#include "SaddlePointSolver.hpp"
#include "Xtra_Funcs.hpp"
#include "inputParameters.h"
#include "ActionAnalytics.hpp"
#include "AggregateChannels.hpp"

///Additional Headers
#include <complex>
#include <complex_bessel.h>

const double Pi = 3.1415926535897932384626433832795028841971693993751;

int main(int argc, char* argv[]) {
	const std::string PATH="results/";
	const std::string PATH2="Data/Amplitudes/";
	///Make directories if they don't exist
	struct stat st;
	if(stat(PATH.c_str(),&st)==-1)
	{
		mkdir(PATH.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	if(stat("Data/",&st)==-1)
	{
		mkdir("Data/",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	if(stat(PATH2.c_str(),&st)==-1)
	{
		mkdir(PATH2.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

	std::string inputFilename;
	if(argc > 1 ){
		inputFilename = argv[1];
	}
	else{
		std::cout<<"No config filename given, using default config filename config.inp\n";
		inputFilename = "config.inp";
	}
	inputParameters Input1(inputFilename);
	Input1.ExtractValues();

	std::cout << Input1.getNtheta() << " " << Input1.getNrad() << " " << Input1.getThetaStart() << " " << Input1.getThetaEnd() << " " << Input1.getThetaStep()
			<< " " << Input1.getRStart() << " " << Input1.getREnd() << " " <<  Input1.getRStep() << std::endl;

	int NumofPulses = Input1.getNumofPulses();
	for(int iPulse = 0; iPulse < NumofPulses; ++iPulse)
	{
		std::cout<<"Starting pulse Loop: "<<iPulse+1<<" of "<<NumofPulses<<std::endl;
		std::vector<double> LaserData = Input1.getLaserValues().at(iPulse);

		double up = LaserData.at(0);
		double omega = LaserData.at(1);

		//std::vector<std::vector<int> > OrbitNums = Input1.getChannelOrbitalNums();
		std::vector<std::vector<double> > IonPots = Input1.getChannelIonisationPotentials();

		std::vector< gsl_complex > times_short = Input1.getTimesShort();
		std::vector< gsl_complex > times_long = Input1.getTimesLong();


		double rStart, rEnd, thetaStart, thetaEnd, rStep, thetaStep;
		int  Nrad, Ntheta;

		rStart = Input1.getRStart();
		rEnd = Input1.getREnd();
		thetaStart = Input1.getThetaStart();
		thetaEnd = Input1.getThetaEnd();
		Nrad = Input1.getNrad();
		Ntheta = Input1.getNtheta();
		rStep = Input1.getRStep();
		thetaStep = Input1.getThetaStep();

		std::vector<double> EnergyRanges;
		EnergyRanges.push_back(thetaStart);
		EnergyRanges.push_back(thetaEnd);
		EnergyRanges.push_back(thetaStep);
		EnergyRanges.push_back(rStart);
		EnergyRanges.push_back(rEnd);
		EnergyRanges.push_back(rStep);

		int NPref = Input1.getNPref();
		for(int iPref=0; iPref<NPref; iPref++)
		{
			std::cout<<"Prefactor Sum: "<<iPref+1<<" of "<<NPref<<std::endl;
			std::vector<bool> PrefactorChoice = Input1.getPrefactors().at(iPref);

			int prefactors_sum=0;
			for (int i=0;i<3;i++)
			{
				int bintemp=PrefactorChoice.at(i)? pow(2,i):0;
				prefactors_sum+=bintemp;
			}
			const int orig_prefactors_sum = prefactors_sum;

			///Output File setup
			std::ostringstream AmpFileNameStub;
			AmpFileNameStub << PATH2 << inputFilename << "Pre_" << orig_prefactors_sum <<"_Range_" << Ntheta<<"X"<<Nrad<<"Up_"<<iPulse<<"_channel";


			std::vector< std::vector< ResultElement > > sum_table_short;
			std::vector< std::vector< ResultElement > > sum_table_long;


			std::cout<<"Starting Channel Loop"<<std::endl;
			int Nchan = Input1.getNumofChannels();
			for(int ichan=0; ichan<Nchan; ++ichan)
			{
				std::cout<<"Doing "<<ichan+1<<" of "<<Nchan<<" Channels"<<std::endl;
				double ip1 = IonPots.at(ichan).at(0);
				double ip2 = IonPots.at(ichan).at(1);
				double ip12 = IonPots.at(ichan).at(2);

				I_Action * f = NULL;
				f = new DefaultAction(  up, omega, ip1, ip2, ip12 );

				//ToDo: Orbital nums and prefactors initialisedhere, using Usecoulomb
				//Dummy declaration needed to initialise Aggregator, not filled
				std::vector< const AbstractPrefactor* > prefactors;

				ResultGrid gridShort(f, rStart, rEnd, thetaStart, thetaEnd, rStep, thetaStep, Nrad, Ntheta, times_short);

				try {
					std::cout<<"trying to populate short table"<< std::endl;
					gridShort.populateTable();
				}
				catch (ResultElement &elem)
				{
					double plong = elem.getPLong();
					double ptrans = elem.getPTrans();
					std::complex<double> t1 = GSL2C(times_short.at(0));
					std::complex<double> t2 = GSL2C(times_short.at(1));
					std::cout<<"Error Short Orbit Finding Failed with momentum values"<<std::endl<<"p1long = "<<plong<<", p1trans = "<<ptrans<<std::endl;
					std::cout<<"and times:"<<std::endl<<"t1 = "<<t1<<", t2 = "<<t2<<std::endl;
				}
				std::cout<<"short table populated"<<std::endl;

				ResultGrid gridLong(f, rStart, rEnd, thetaStart, thetaEnd, rStep, thetaStep, Nrad, Ntheta, times_long);
				try{
					gridLong.populateTable();
				}
				catch (ResultElement &elem)
				{
					double plong = elem.getPLong();
					double ptrans = elem.getPTrans();
					std::complex<double> t1 = GSL2C(times_long.at(0));
					std::complex<double> t2 = GSL2C(times_long.at(1));
					std::cout<<"Error Long Orbit Finding Failed with momentum values"<<std::endl<<"p1long = "<<plong<<", p1trans = "<<ptrans<<std::endl;
					std::cout<<"and times:"<<std::endl<<"t1 = "<<t1<<", t2 = "<<t2<<std::endl;
	    		}
				std::cout<<"tables populated"<<std::endl;

				//ToDo: add pi
				/*
				if (Input1.getAddPi() || true) {
					std::cout<<"Shifting times by half a cycle\n";
					try {
						gridShort.PiTransformTable();
						gridLong.PiTransformTable();
					} catch (std::exception& e) {
						std::cerr << "Failed to shift times by halfa cycle\n";
						std::cerr<<e.what()<<std::endl;
						return -1;
					}

				}
				 */
				const std::vector<std::vector<ResultElement>>& tableShort = gridShort.getResultTable();

				const std::vector< std::vector< ResultElement > >& tableLong = gridLong.getResultTable();

				std::cout<<"Channel: "<<ichan+1<<" Short table size: "<<tableShort.at(0).size()<<std::endl<<"Long table size: "<<tableLong.size()<<std::endl;

				/// dump the orbits.
				std::ostringstream o_short_orbit;

				o_short_orbit << PATH <<inputFilename<<"_Range_"<< Ntheta<<"X"<<Nrad<<"_channel:"<<ichan<< ".shortorbit";

				std::cout <<"Channel: "<<ichan+1<< " Saving short orbit times" << std::endl;
				std::ofstream short_orbit_file;
				short_orbit_file.open( o_short_orbit.str().c_str() );

				std::vector< std::vector< ResultElement > >::const_iterator short_grid_it = tableShort.begin();

				for(int i=0 ; short_grid_it != tableShort.end() ; ++short_grid_it, ++i )
				{

					std::vector< ResultElement >::const_iterator it2 = short_grid_it ->begin();
					//short_orbit_file << std::endl;
					for(int j = 0; it2 != short_grid_it ->end() ; ++it2, ++j )
					{
						short_orbit_file << std::fixed << std::setprecision(8) << i*thetaStep+thetaStart << " " << j*rStep+rStart;

						const std::vector< gsl_complex >& tmp = it2->getValues();
						std::vector< gsl_complex >::const_iterator it3 = tmp.begin();

						for(int k=0 ; it3 != tmp.end() ; ++it3, ++k )
						{
							double imag = GSL_IMAG( *it3 );
							short_orbit_file<< std::setprecision(12) << " " << GSL_REAL( *it3 ) << " " << imag;
						}

						short_orbit_file << std::endl;
					}
				}

				short_orbit_file.close();

				std::cout << "Channel: "<<ichan+1<<" Saving long orbit times" <<  std::endl;
				std::ostringstream o_long_orbit;

				o_long_orbit << PATH <<inputFilename<<"_Range_"<< Ntheta<<"X"<<Nrad<<"_channel:"<<ichan<<".longorbit";

				std::ofstream long_orbit_file;
				long_orbit_file.open( o_long_orbit.str().c_str() );

				std::vector< std::vector< ResultElement > >::const_iterator long_grid_it = tableLong.begin();

				/*
				 * Not sure why this is here oof
				if(ichan==0)
				{
					sum_table_long = tableLong;
					sum_table_short = tableShort;
				}
				*/

				for(int i=0 ; long_grid_it != tableLong.end() ; ++long_grid_it, ++i )
				{
					std::vector< ResultElement >::const_iterator it2 = long_grid_it ->begin();
					//long_orbit_file << std::endl;
					for(int j=0 ; it2 != long_grid_it ->end() ; ++it2, ++j )
					{
						long_orbit_file << std::fixed << std::setprecision(8) << i*thetaStep+thetaStart << " " << j*rStep+rStart;

						const std::vector< gsl_complex >& tmp = it2->getValues();
						std::vector< gsl_complex >::const_iterator it3 = tmp.begin();
						for(int k=0 ; it3 != tmp.end() ; ++it3, ++k )
						{
							double imag = GSL_IMAG( *it3 );
							long_orbit_file <<std::setprecision(12)<< " " << GSL_REAL( *it3 ) << " " << imag;
						}
						long_orbit_file << std::endl;
					}
				}

				long_orbit_file.close();



				AggregateChannels aggro1(f, tableShort, tableLong, prefactors, EnergyRanges);

				///Amplitude File Out
				std::ostringstream AmpFileName;
				AmpFileName <<AmpFileNameStub.str().c_str()<<ichan<< ".AmpOut";


				std::cout <<std::endl<<"Channel: "<<ichan+1<< " Running Aggregator" << std::endl;
				///Create/Load/Save amplitude related data

				//Filling up aggregator, caching times and f3 values
				aggro1.LoadTimes();


				/// dump the amplitudes.
				std::ostringstream o_amps;
				o_amps << PATH2 <<inputFilename<<"_Range_"<< Ntheta<<"X"<<Nrad<<"_channel:"<<ichan<< ".outamps";
				std::cout <<"Channel: "<<ichan+1<< " Saving amplitudes" << std::endl;
				std::ofstream o_amps_file;
				o_amps_file.open( o_amps.str().c_str() );

				//Function has been augmented to print calculated amplitudes to a file,
				//in a 2d array, (ppar, pper, Re[action], Im[action])
				aggro1.SaveAmplitudes(Input1.getPrintCount(), o_amps_file, true);



				delete f;
				std::cout<<std::endl;
			}
		}
	}
}
