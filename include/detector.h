/********
*Author: Baran Bodur
*Date: 2019-09-12
*Description: To read, store and organize detector info
*
********/
#ifndef DET_INCLUDED
#define DET_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include "TF1.h"


class Detector{
public:
	// Constructor
	Detector(std::string cardName, TFile* outFile);
	
	// Function To Give a Random Point in the Detector Considering Detector Geometry
	TVector3 randomInteractionVertex();

	// Variables Filled From Card File
	std::string strName; // name
	double nAvogadro; // avogadros number *10^-23
	std::string material; // main detector material (will be passed to xscn)
	double molarMass; // of the material g/mol
	double density; // in g/cm3 or metricton/m^3
	double detMass; // in tons mass of the active volume
	double detYears; // in years (the time that is simulated)
	double yearsToSeconds; // *10^-7 
	double unitFixer; // fixes the units (see card file
	double overallCoeff; // so that this times*xscn*flux will give event count at given energy,angle...
	
	std::string detGeom; // Whether detector is rectangular,cylindirical or spherical
	double detX, detY, detZ, detR;

private:
	TRandom3 detRand;
	TF1 *randCyl;
	TF1 *randSph;
};


#endif
