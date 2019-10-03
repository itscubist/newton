/********
*Author: Baran Bodur
*Date:  2019-10-01
*Description: Excited States Class, Initialized in Xscn Class, Controls Talys Interface 
* 						and other tools for decays
********/
#ifndef EXSTATES_INCLUDED
#define EXSTATES_INCLUDED

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
#include "xscn.h"
#include "detector.h"

const unsigned int ZMAX = 10;
const unsigned int NMAX = 10;
const unsigned int PMAX = 7;
const unsigned int EXINMAX = 150;
const unsigned int EXOUTMAX = 150;


// order: gamma, neutron,  proton, deuteron, triton, he-3, alpha
enum particle {g,n,p,d,t,h,a};
const int PDGID[7] = {22,2112,2212,1000010020,1000010030,1000020030,1000020040}; 
const int ZRED[7] = {0,0,1,1,1,2,2};
const int NRED[7] = {0,1,0,1,2,1,2};

// Will use TALYS Interface Functions So Extern Them
extern "C" {
	void machine_();
	void constants_();
	void talysinputcread_(char*,char*,int*,int*,int*,int*,float*);
	void talysinitial_(); 
	void talysreaction_(); 
	void gettalysresults_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,
			float*,float*,float*,float*,float*,float*); 
}

// Class For Excited Levels
class Exstates{
public:
	// Constructor, Destructor
	Exstates(double inEnergyGnd=0.0);
	~Exstates();

	// Functions
	void initHists(); // initialize huge arrays of histograms
	void runTalys(); // runs TALYS for this exited state
	void fillHists(); // fills histograms by questioning TALYS 
	TString getHName(int z, int n, int exc); // get a string for given values

	// Variables
	double energyInitial; // Energy From The Initial State
	double energyGnd; // Energy From The Ground State of the Final Product
	int spin; // Spin of Excited State
	int parity; // Parity of Excited State
	int nState; // The excited state no
	Xscn* xscn;
	TString fileName;
	TFile* talysFile;

	// Actual Upper Limits Based on TALYS output
	int nZ, nN, nPar, nExc, nExcout;
 		
	TH2D* decayHists[ZMAX][NMAX][EXINMAX];
	double sepEnergies[ZMAX][NMAX][PMAX];
	double excEnergies[ZMAX][NMAX][EXINMAX];
	double binW[ZMAX][NMAX][EXINMAX];
	
};

#endif
