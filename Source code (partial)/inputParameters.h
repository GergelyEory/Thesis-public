/*
 * inputParameters.h
 *
 *  Created on: 25 Aug 2019
 *      Author: geory
 */

#ifndef INPUTPARAMETERS_H_
#define INPUTPARAMETERS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <libconfig.h++>
#include <gsl/gsl_complex.h>
#include <vector>
#include <sstream>

using namespace libconfig;

class inputParameters {
	public:
		//ctor
		inputParameters(std::string inputFilename);

		//dtor
		virtual ~inputParameters();

		void openFile(std::string inputFilename);

		std::string getFilename() const;
		std::string getPATH() const;

		void ExtractValues();
		void ExtractValues(std::string Input_filename);
		///Get specific values
		std::vector<std::vector<bool> > &getPrefactors();
		int getNPref();
		std::vector<std::vector<double> > &getLaserValues();
		int getNumofPulses();
		int getNumofChannels();
		std::vector<std::vector<double> > &getChannelIonisationPotentials();
		std::vector<std::vector<int> > &getChannelOrbitalNums();
		std::vector<gsl_complex> &getTimesShort();
		std::vector<gsl_complex> &getTimesLong();
		int getNrad();
		int getNtheta();
		double getREnd();
		double getRStart();
		double getRStep();
		double getThetaEnd();
		double getThetaStart();
		double getThetaStep();
		bool getTesting();
		bool getAddPi();
		bool getSaveAmplitudes();
		bool getLoadAmplitudes();
		bool getPrintCount();

	private:
		template<typename TYPE>
		void lookup(const std::string& varName, TYPE& var);
		void lookupString(const std::string& varName, std::string& varStr);

		std::string configFilename;
		std::string PATH="Data/";


		///Vectors of Values and their lengths
		std::vector<std::vector<bool> > Prefactors;
		int NPref;
		std::vector<std::vector<double> > LaserValues;
		int NumofPulses;
		///channels values
		int NumofChannels;
        std::vector<std::vector<int> > ChannelOrbitalNums;
		std::vector<std::vector<double> > ChannelIonisationPotentials;
		///initial start time values
		std::vector<gsl_complex> TimesShort;
		std::vector<gsl_complex> TimesLong;
		///energy ranges
		double rStart, rEnd, thetaStart, thetaEnd, rStep, thetaStep;
		int  Nrad, Ntheta;
		bool testing;
		bool addPiQ;

		std::ifstream inputFile;
		Config cfg;
};

#endif /* INPUTPARAMETERS_H_ */
