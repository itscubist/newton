/********
*Author: Baran Bodur
*Date:  2019-09-13
*Description: Store and Organize Access to Properties of a Specific Event
*
********/
#ifndef EVENT_INCLUDED
#define EVENT_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TFile.h"
#include "TVector3.h"
#include "TString.h"
#include "TRandom3.h"

// Project Headers
#include "global.h"
#include "flux.h"
#include "xscn.h"
#include "detector.h"


// Struct For Excited Levels
struct Particle{
	double energy; // Energy From The Initial State
	double mass; // Energy From The Ground State of the Final Product
	int pdgId; // Pdg ID
	TVector3 direction; // Parity of Excited State
	bool isSimulated;
};


class Event{
public:
	// Constructor
	Event(double inEnergy, Xscn &inXscn, Flux &inFlux, Detector &inDet);

	// Functions To Write Event to Kin File
	// void writeEvent(ofstream &outKinFile);
	// Lepton is the 0th indice, hadron is the 1st indice, de-excitation stuff is 2nd, 3rd...
	void fillLeptonDirAndEnergy(); // call this first
	//void fillHadronDirAndEnergy();
	// A handy kinematic function
	TVector3 combineZenAzi(double cosZen1,double azi1, double cosZen2, double azi2);

	// Info related to Event
	std::string xscnName; // Xscn Name
	std::string fluxName; // Flux Name
	double nuEnergy;
	unsigned int nuFlv;
	TVector3 vertex;
	double time;

	// Info Needed To Calculate Direction and Energy of Output Lepton 
	double fluxCosZen, fluxAzi; // zenith and azimuth from flux
	unsigned int xscnExState; // Decided Excited State No
	double initialTargetMass; // Mass of the target
	double finalHadMass; // Mass of the outgoing final (with excitation energy)
	double xscnCosAngle; // cos(ThetaNuE) from xscn
	double xscnAzi; // uniform azi

	
	// Output Particles Related To This Event
	std::vector<Particle> particles;

private:
};


#endif
