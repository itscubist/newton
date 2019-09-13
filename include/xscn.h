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
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TRandom3.h"

// Project Headers
#include "global.h"
#include "event.h"
#include "flux.h"
#include "detector.h"


// Struct For Excited Levels
struct exStates{
	double energyInitial; // Energy From The Initial State
	double energyGnd; // Energy From The Ground State of the Final Product
	int spin; // Spin of Excited State
	int parity; // Parity of Excited State
};


class Xscn{
public:
	// Constructor
	Xscn(std::string cardName, std::string material, TFile* outFile);

	// Functions To Get Xscn At Asked Values
	double xscnAtEnergy(double energy);
	double xscnAtEnergyAngle(double energy, double cosAngle, unsigned int exLevel);
	unsigned int randomExLevelAtEnergy(double energy); 
	double randomAngleAtEnergy(double energy, unsigned int exLevel);
	double randomAzimuth(); //  Uniform Accross Azimuth
	// Function To Generate Event Counts with Given Flux
	int generateEventCountPerEnergy(Flux inFlux, Detector inDet);

	// Info related to cross seciton
	std::string strName; // name
	double targetPerMolecule; // how many targets per molecule
	unsigned int nFinalStates; // number of possible final states (nuclear excitation and so on)
	int pdgTarget; // pdg of target
	int pdgFinal0; // pdg of produced final nucleus, hadron...
	int pdgLepton; // pdg of procuded lepton
	double mTarget; // mass of target
	double mFinal0; // ground state mass of produced final nucleus, hadron...
	double mLepton; // mass of produced lepton

	bool writeFinal0; // Whether final hadron should be written to kin file 
	bool decayFinal0; // Whether final hadron's decay from excited state should be simulated
	bool intNu[6]; // Neutrinos participate in the interaction
	bool doubleDiff; // Whether double differential xscn exists or not
	bool excData; // Whether Excited state data exists
	
	unsigned int angBins; // Number of Angular Bins
	unsigned int energyBins; // Number of Energy Bins

	// Where Info is Stored
	TH1D* xscnVsEnergy;
	std::vector<TH2D*> xscnVsEnergyAngle;
	std::vector<exStates> excLevels;
	TH2D* excProbVsEnergy;
	
	// Generated Events
	TH1D* eventsVsEnergy;
	std::vector<double> genEnergies;
	TH3D* leptonDirEnergyDist;

private:
	// RNG
	TRandom3 xscnRand; 

	// Functions for Initialization
	TH1D* readXscn(); // Read Single Diff Xscn
	std::vector<TH2D*> readXscnDouble(); // Read Double Diff Xscn Per Exc Energy
	TH2D* readExcProb(); // Read Exc Prob vs Energy
	std::vector<exStates> readExcStates(); // Read Exc State Info
	
	// Filenames to Read Info From	
	std::string strXscnVsEnergy;
	std::string strExcStates;
	std::string strExcProb;
	std::vector<std::string> strXscnVsEnergyAndAngle;

};


#endif
