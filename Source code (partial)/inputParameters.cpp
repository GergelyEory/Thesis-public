/*
 * inputParameters.cpp
 *
 *  Created on: 25 Aug 2019
 *      Author: geory
 */

#include "inputParameters.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <gsl/gsl_complex_math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <libconfig.h++>
#include <gsl/gsl_complex.h>
#include <vector>
#include <sstream>

constexpr double pi=3.14159265358979323846264338327950288419716939937510;



inputParameters::inputParameters(std::string inputFilename) { // @suppress("Class members should be properly initialized")
	///Default values for variables that are necessary to set
	NumofChannels=1;
	thetaStart = 0.0001 * pi;
	thetaEnd = 1.9999 * pi;

	//read variables from config file

	try{
		openFile(inputFilename);
	}
	catch(const FileIOException &fioex){
		std::cerr << "I/O error while reading file." << std::endl;
		std::terminate();
	}
	catch(const ParseException &pex){
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
			            		  << " - " << pex.getError() << std::endl;
		std::terminate();
	}
	configFilename = inputFilename;

}

//dtor
inputParameters::~inputParameters() {
}

void inputParameters::openFile(std::string inputFilename){
	try{
		cfg.readFile(inputFilename.c_str());
	}
	catch(const FileIOException &fioex)
	{
		std::cerr<<"Error reading configuration file\n";
	}
}
template<typename TYPE>
void inputParameters::lookup(const std::string& varName, TYPE& var){
	try{
		var = cfg.lookup(varName);
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << varName<<" not set in configuration file\n";
	}
}
void inputParameters::lookupString(const std::string& varName, std::string & varStr)
{
	try{
		cfg.lookupValue(varName,varStr);
	}
	catch(const SettingNotFoundException &nfex){
			std::cerr << varName<<" not set in configuration file\n";
	}
}

void inputParameters::ExtractValues() {
	//Type taken int rather than size_t by libconfig parser so a type error is returned when setting these variables
	//So AutoConvert is turned on to fix this issue
	cfg.setAutoConvert(true);
	lookup("NumofPulses", NumofPulses);
	lookup("NPref", NPref);
	lookup("NumofChannels", NumofChannels);
	cfg.setAutoConvert(false);//It is turned off again to fix against incorrect types being input


	try {
		const Setting& Laser = cfg.lookup("LaserValues");
		LaserValues.clear();
		for (int i = 0; i < NumofPulses;++i) {
			const Setting& pulse = Laser[i];
			double Up, w, phi;
			pulse.lookupValue("Up", Up);
			pulse.lookupValue("w", w);
			pulse.lookupValue("phi", phi);
			std::vector<double> pulsevalues = {Up, w, phi};
			LaserValues.push_back(pulsevalues);
		}
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "LaserValues"<<" not set in configuration file\n";
	}

	try {
		const Setting& preflist = cfg.lookup("preflist");
		Prefactors.clear();
		for (int i = 0; i < NPref;++i) {
			const Setting& pref = preflist[i];
			bool a, b, c;
			pref.lookupValue("Up", a);
			pref.lookupValue("w", b);
			pref.lookupValue("phi", c);
			std::vector<bool> pulsevalues = {a,b,c};
			Prefactors.push_back(pulsevalues);
		}
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "preflist"<<" not set in configuration file\n";
	}

	try {
		const Setting& ChannelParams = cfg.lookup("ChannelParams");
		ChannelOrbitalNums.clear();
		ChannelIonisationPotentials.clear();
		for (int i = 0; i < NumofChannels;++i) {
			const Setting& channel = ChannelParams[i];
			int ng, lg, mg, ne, le, me;
			double I1g, I2g, I2e;

			cfg.setAutoConvert(true);
			channel.lookupValue("ng", ng);
			channel.lookupValue("lg", lg);
			channel.lookupValue("mg", mg);
			channel.lookupValue("ne", ne);
			channel.lookupValue("le", le);
			channel.lookupValue("me", me);
			cfg.setAutoConvert(false);

			channel.lookupValue("I1g", I1g);
			channel.lookupValue("I2g", I2g);
			channel.lookupValue("I2e", I2e);

			std::vector<int> OrbitalNums = {ng, lg, mg, ne, le, me};
			ChannelOrbitalNums.push_back(OrbitalNums);
			std::vector<double> IonPots = {I1g, I2g, I2e};
			ChannelIonisationPotentials.push_back(IonPots);
		}
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "ChannelParams"<<" not set in configuration file\n";
	}

	try {
		double t1rs, t1is, t2rs, t2is, t3rs, t3is;
		TimesShort.clear();
		cfg.lookupValue("t1Rshort",t1rs);
		cfg.lookupValue("t1Ishort",t1is);
		cfg.lookupValue("t2Rshort",t2rs);
		cfg.lookupValue("t2Ishort",t2is);
		cfg.lookupValue("t3Rshort",t3rs);
		cfg.lookupValue("t3Ishort",t3is);

		gsl_complex t1, t2, t3;
		t1 = gsl_complex_rect(t1rs, t1is);
		t2 = gsl_complex_rect(t2rs, t2is);
		t3 = gsl_complex_rect(t3rs, t3is);
		TimesShort.push_back(t1);
		TimesShort.push_back(t2);
		TimesShort.push_back(t3);
	}
	catch (const SettingNotFoundException &nfex){
		std::cerr << "TimesShort"<<" not set in configuration file\n";
	}

	try {
		cfg.setAutoConvert(true);
		cfg.lookupValue("testing", testing);
		cfg.setAutoConvert(false);
	} catch (const SettingNotFoundException &nfex) {
		testing = false;
		std::cerr << "testing" << "not set in configuration file\n";
	}

	try {
		cfg.setAutoConvert(true);
		cfg.lookupValue("addPiQ", addPiQ);
		cfg.setAutoConvert(false);
	} catch (const SettingNotFoundException &nfex) {
		addPiQ = false;
		std::cerr << "addPiQ" << "not set in configuration file\n";
	}

	try {
		double t1rl, t1il, t2rl, t2il, t3rl, t3il;
		TimesLong.clear();
		cfg.lookupValue("t1Rlong",t1rl);
		cfg.lookupValue("t1Ilong",t1il);
		cfg.lookupValue("t2Rlong",t2rl);
		cfg.lookupValue("t2Ilong",t2il);
		cfg.lookupValue("t3Rlong",t3rl);
		cfg.lookupValue("t3Ilong",t3il);

		gsl_complex t1, t2, t3;
		t1 = gsl_complex_rect(t1rl, t1il);
		t2 = gsl_complex_rect(t2rl, t2il);
		t3 = gsl_complex_rect(t3rl, t3il);
		TimesLong.push_back(t1);
		TimesLong.push_back(t2);
		TimesLong.push_back(t3);
	}
	catch (const SettingNotFoundException &nfex){
		std::cerr << "TimesLong"<<" not set in configuration file\n";
	}

	try {
		cfg.setAutoConvert(true);
		cfg.lookupValue("Nrad",Nrad);
		cfg.lookupValue("Ntheta",Ntheta);
		cfg.setAutoConvert(false);
		cfg.lookupValue("rStart",rStart);
		cfg.lookupValue("rEnd",rEnd);
		Ntheta = 4*Ntheta + 1;
		//cfg.lookupValue("thetaStart",thetaStart);
		//cfg.lookupValue("thetaEnd",thetaEnd);

		thetaStep= (thetaEnd-thetaStart)/static_cast<double>(Ntheta-1);
		rStep = (rEnd-rStart)/static_cast<double>(Nrad-1);


	}
	catch (const SettingNotFoundException &nfex){
		std::cerr << "Grid values"<<" not set in configuration file\n";
	}
}



void inputParameters::ExtractValues(std::string Input_filename) {
}


std::string inputParameters::getFilename() const {	return configFilename;}
std::string inputParameters::getPATH() const {	return PATH;}

std::vector<std::vector<bool> >& inputParameters::getPrefactors() {
	return Prefactors;
}

int inputParameters::getNPref() {
	return NPref;
}

std::vector<std::vector<double> >& inputParameters::getLaserValues() {
	return LaserValues;
}

int inputParameters::getNumofPulses() {
	return NumofPulses;
}

int inputParameters::getNumofChannels() {
	return NumofChannels;
}
std::vector<std::vector<double> >& inputParameters::getChannelIonisationPotentials() {
	return ChannelIonisationPotentials;
}

std::vector<gsl_complex>& inputParameters::getTimesShort() {
	return TimesShort;
}

std::vector<gsl_complex>& inputParameters::getTimesLong() {
	return TimesLong;
}

std::vector<std::vector<int> >& inputParameters::getChannelOrbitalNums() {
	return ChannelOrbitalNums;
}

int inputParameters::getNrad() {
	return Nrad;
}

int inputParameters::getNtheta() {
	return Ntheta;
}

double inputParameters::getREnd() {
	return rEnd;
}

double inputParameters::getRStart() {
	return rStart;
}

double inputParameters::getRStep() {
	return rStep;
}

double inputParameters::getThetaEnd() {
	return thetaEnd;
}

double inputParameters::getThetaStart() {
	return thetaStart;
}

double inputParameters::getThetaStep() {
	return thetaStep;
}

bool inputParameters::getTesting() {
	return testing;
}

bool inputParameters::getAddPi() {
	return addPiQ;
}

bool inputParameters::getSaveAmplitudes() {
	return true;
}

bool inputParameters::getLoadAmplitudes() {
	return false;
}

bool inputParameters::getPrintCount() {
	//Always prints percentage data
	return true;
}
