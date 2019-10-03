/********
*Author: Baran Bodur
*Date: 2019-10-01
*Description: Functions For TALYS interface and More for decay of excited states 
*
********/

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TString.h"

// Project Headers
#include "global.h"
#include "xscn.h"
#include "exstates.h"

using namespace std;

Exstates::Exstates(double inEnergyGnd) { // Constructor
	energyGnd = inEnergyGnd;
}

Exstates::~Exstates() { // Destructor
	if(nState<1) return;
	for (int z = 0; z < nZ; z++)
		for (int n = 0; n < nN; n++)
			for (int e = 0; e < nExc; e++) {
				delete decayHists[z][n][e]; // delete histograms
			}
}
	

void Exstates::runTalys() { // runs TALYS for this exited state
	if(nState<1) return;
	//TALYS hacking stuff
	char inputfilename[64] = "talysInterface/talysinput";
	char projectile = '0';
	int tarZ = xscn->zFinalNuc, tarA = xscn->aFinalNuc;
	float exEne = static_cast<float>(energyGnd);
	machine_();
	constants_();
	talysinputcread_(inputfilename,&projectile,&tarZ,&tarA,&spin,&parity,&exEne);
	talysinitial_(); 
	talysreaction_();
	return;
}

// Functions
void Exstates::initHists() { // initialize huge arrays of histograms
	if(nState<1) return;
	fileName = "talysInterface/talysOutputs/"+(TString)xscn->strName+Form("_state%d.root",nState);
	talysFile = new TFile(fileName,"RECREATE");
	// To get size of histograms
	int tZ=0, tN=0, tPar=0, tExc=0, tExcout=0; // temporary placeholders
	float decProb,sepEnergy,excEnergy1,excBinWidth1,excEnergy2,excBinWidth2;
	gettalysresults_(&tZ,&tN,&tPar,&tExc,&tExcout,&nZ,&nN,&nPar,&nExc,&nExcout,
			&decProb,&excEnergy1,&excBinWidth1,&excEnergy2,&excBinWidth2,&sepEnergy);
	// Z and N loop means how many protons or neutron less than the initial nucleus
	for (int z = 0; z < nZ; z++) {// loop over Z non zero histograms
		for (int n = 0; n < nN; n++) {// loop over N
			for (int e = 0; e < nExc; e++) { // loop over excited states of given Z and N
				cout << "Z: " << z << " N: " << n << " Exc: " << e << endl; 
				// x axis particle, y axis decayed excitation level
				decayHists[z][n][e] = new TH2D(getHName(z,n,e),"decayHist",7,-0.5,6.5,nExcout,0,nExcout);
				cout << "Zm: " << ZMAX << " Nm: " << NMAX << " Exc m: " << EXINMAX << endl; 
				cout << "*** in Init Hists: " << getHName(z,n,e) << endl;
				for (int p = 0; p < nPar; p++) { // loop over excited states of given Z and N
					for (int eo = 0; eo < nExcout; eo++) { // loop over excited states of given Z and N
						//if(n==1) cout << "pType: " << p << " OutputEx: " << eo << endl;
						//gettalysresults_(&z,&n,&p,&e,&eo,&tZ,&tN,&tPar,&tExc,&tExcout,
						//	&decProb,&excEnergy1,&excBinWidth1,&excEnergy2,&excBinWidth2,&sepEnergy);
						// Fill histogram
						//decayHists[z][n][e]->SetBinContent(p+1,eo+1, static_cast<double>(decProb));
					} // end of output excited states
					//sepEnergies[z][n][p]  = static_cast<double>(sepEnergy);
				} // end of particle type
				//excEnergies[z][n][e] = static_cast<double>(excEnergy1);
				//binW[z][n][e] = static_cast<double>(excBinWidth1);
				decayHists[n][z][e]->Write(); // save histograms
			} // end of excited states
		} // end of neutron number
	} // end of proton number
		
	return;
}

void Exstates::fillHists() { // fills histograms by questioning TALYS 
	if(nState<1) return;
	int jZ=0, jN=0, jPar=2, jExc=4, jExcout=0; 
	float decProb,sepEnergy,excEnergy1,excBinWidth1,excEnergy2,excBinWidth2;
	gettalysresults_(&jZ,&jN,&jPar,&jExc,&jExcout,&nZ,&nN,&nPar,&nExc,&nExcout,
			&decProb,&excEnergy1,&excBinWidth1,&excEnergy2,&excBinWidth2,&sepEnergy);

	cout << " ************ Test of TALYS READ ********** " << endl;
	cout << "Z looped: " << nZ << endl;
	cout << "N looped: " << nN << endl;
	cout << "Particle looped: " << nPar << endl;
	cout << "Initial Exc State looped: " << nExc << endl;
	cout << "Final Exc State looped: " << nExcout << endl;

	cout << "Prob: " << decProb << endl;
	cout << "Separation Energy: " << sepEnergy << endl;
	cout << "Starting Excitation Energy: " << excEnergy1 << endl;
	cout << "Starting Excitation Energy Width: " << excBinWidth1 << endl;
	cout << "Final Excitation Energy: " << excEnergy2 << endl;
	cout << "Final Excitation Energy Width: " << excBinWidth2 << endl;
	cout << " ************ END OF TALYS TEST *********** " << endl;
	return;
}

TString Exstates::getHName(int z, int n, int exc) { // get a string for given values
	int zCur = xscn->zFinalNuc - z; 
	int nCur = (xscn->aFinalNuc - xscn->zFinalNuc) - n; 
	TString name = "h"+(TString)xscn->strName+Form("_state%d_z%d_n%d_exc%d",nState,zCur,nCur,exc);
	return name;
}
