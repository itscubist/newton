/********
*Author: Baran Bodur
*Date:  2019-09-10
*Description: To read, store and organize access to xscn data
*
********/
#ifndef XSCN_INCLUDED
#define XSCN_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

// Struct For Excited Levels
struct exStates{
	float energy; // Energy From The Ground State of the Final Product
	int spin; // Spin of Excited State
	int parity; // Parity of Excited State
};


class Xscn{
public:
	// Functions
	Xscn(std::string cardName, std::string material); // Constructor
	TH1D* readXscn(); // Read Single Diff Xscn
	std::vector<TH2D*> readXscnDouble(); // Read Double Diff Xscn Per Exc Energy
	TH2D* readExcProb(); // Read Exc Prob vs Energy
	std::vector<exStates> readExcStates(); // Read Exc State Info
	
	str::string strName; 
	double targetPerMolecule;
	unsigned int nFinalStates;
	int pdgTarget;
	int pdgFinal0;
	int pdgLepton;
	double mTarget;
	double mFinal0;
	double mLepton;

	bool writeFinal0;
	bool decayFinal0;
	bool intNu[6]; // Neutrinos participate in the interaction
	bool doubleDiff; // Whether double differential xscn exists or not
	bool excData; // Whether Excited state data exists
	
	unsigned int angBins; // Number of Angular Bins
	unsigned int energyBins; // Number of Energy Bins
	// Filenames to Read Info From	
	std::string strXscnVsEnergy;
	std::string strExcStates;
	std::string strExcProb;
	std::vector<std::string> strXscnVsEnergyAndAngle;

	// Where Info is Stored
	TH1D* xscnVsEnergy;
	std::vector<TH2D*> xscnVsEnergyAngle;
	std::vector<exStates> excLevels;
	TH2D* excProbVsEnergy;

private:

};


#endif
