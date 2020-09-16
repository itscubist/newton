/********
*Author: Baran Bodur
*Date: 2019-09-11
*Description: To read, store and organize flux info
*
********/
#ifndef FLUX_INCLUDED
#define FLUX_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TString.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH1D.h"


class Flux{
public:
	// Constructor
	Flux(std::string cardName, TFile* outFile);
	//~Flux();
	
	// Functions Related To Getting Info
	double fluxAtEnergy(double energy, unsigned int flavor);
	double fluxAtDirection(double energy, double cosAngle, double azimuth, unsigned int flavor);
	double randomZenithAtEnergy(double energy, unsigned int flavor);
	double randomAzimuthAtEnergyAndZenith(double energy, double cosAngle, unsigned int flavor);

	// Variables Filled From Card File
	std::string strName; // name
	double scale; // to scale dimensions of flux
	double energyScale; // to scale dimensions of energy axis
	
	std::string storageType; // In what why flux is given
	std::string dirType; // Direction Type
	double nuFluxes[6]; // Flux value if construting flux
	double energyRange[2]; // Energy values for constructing flux
	double cosAngleRead; // For constructing flux coming from a specific direction
	double azimuthRead; // For constructing flux from a specific direction

	// Number of bins
	unsigned int angBins;
	unsigned int aziBins;
	unsigned int energyBins;

	// Where Data is Stored
	std::vector<TH1D*> fluxVsEnergy;
	std::vector<TH3D*> fluxVsEnergyDirectional;

private:
	// RNG
	TRandom3 fluxRand;
	
	void readFlux();
	void buildSingleEnergyFlux();
	void buildConstantAlongEnergyFlux();

	// File and Histogram Names (From Card File
	std::string strFluxFileName;
	std::vector<std::string> strFluxHistName;
};


#endif
