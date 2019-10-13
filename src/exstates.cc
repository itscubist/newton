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

TalysData::TalysData(Xscn* inXscn) { // Constructor
	xscn = inXscn;
	decayRand.SetSeed(0);
	decayHists.assign(ZMAX, vector<vector<TH2D*> >(NMAX, vector<TH2D*>(EXINMAX) ) );
	sepEnergies.assign(ZMAX, vector<vector<double> >(NMAX, vector<double>(EXINMAX) ) );
	excEnergies.assign(ZMAX, vector<vector<double> >(NMAX, vector<double>(EXINMAX) ) );
	binW.assign(ZMAX, vector<vector<double> >(NMAX, vector<double>(EXINMAX) ) );
	
}

// Set filename
void TalysData::setFileName(TString inString) {
	talysFileName = inString;
}


TalysData::~TalysData() { // Destructor
	for (int z = 0; z < ZMAX; z++)
		for (int n = 0; n < NMAX; n++)
			for (int e = 0; e < EXINMAX; e++) {
			  //if(decayHists[z][n][e] == NULL)
			  	//cout << "Z: " << z << " N: " << n  << " Ex: " << e << " is NULL" << endl; 
			  	//delete decayHists[z][n][e]; // delete histograms
			}
}

// Run TALYS for all excited states excited with equal weight,
// We will only read decay probabilities from TALYS, hence we do not have to
// run each excited state alone, we can instead run all of them together
void TalysData::runTalys() { // runs TALYS for with given configuration
	//TALYS hacking stuff
	char inputfilename[64] = "talysInterface/talysinput";
	char projectile = '0';
	int tarZ = xscn->zFinalNuc, tarA = xscn->aFinalNuc;
	machine_();
	constants_();
	talysinputcread_(inputfilename,&projectile,&tarZ,&tarA);
	for(unsigned int sCtr = 0; sCtr < xscn->excLevels.size(); sCtr++) {
		int parity = xscn->excLevels[sCtr].parity;	
		float energy = static_cast<float>(xscn->excLevels[sCtr].energyGnd);	
		int spin = xscn->excLevels[sCtr].spin;
		float popValue = 1.0; // this is just a normalization issue...
		inputsetexenes_(&energy, &spin, &parity, &popValue);
	}
	talysinitial_(); 
	talysreaction_();
	return;
}



// Initialize all decay histograms via TALYS, for each nucleus and excited state
// find which nucleus and excited state it can decay into with what probability
void TalysData::initHists() {
	TFile* talysFile = new TFile(talysFileName,"RECREATE");
	// To get size of histograms
	int tZ=0, tN=0, tPar=0, tExc=0, tExcout=0; // temporary placeholders
	float decProb,sepEnergy,excEnergy1,excBinWidth1,excEnergy2,excBinWidth2;
	gettalysresults_(&tZ,&tN,&tPar,&tExc,&tExcout,&nZ,&nN,&nPar,&nExc,&nExcout,
			&decProb,&excEnergy1,&excBinWidth1,&excEnergy2,&excBinWidth2,&sepEnergy);
	// Z and N loop means how many protons or neutron less than the initial nucleus
	for (int z = 0; z < nZ; z++) {// loop over Z non zero histograms
		for (int n = 0; n < nN; n++) {// loop over N
			for (int e = 0; e < nExc; e++) { // loop over excited states of given Z and N
				// x axis particle, y axis decayed excitation level
				decayHists[z][n][e]
					= new TH2D(getHName(z,n,e),"decayHist",7,-0.5,6.5,nExcout,-0.5,(double)nExcout-0.5);
				for (int p = 0; p < nPar; p++) { // loop over excited states of given Z and N
					for (int eo = 0; eo < nExcout; eo++) { // loop over excited states of given Z and N
						gettalysresults_(&z,&n,&p,&e,&eo,&tZ,&tN,&tPar,&tExc,&tExcout,
							&decProb,&excEnergy1,&excBinWidth1,&excEnergy2,&excBinWidth2,&sepEnergy);
						/*if(decProb !=0) {
							cout << "*** in Init Hists: " << getHName(z,n,e) << endl;
							cout << "Dec Prob Via: " << p << " to: " << eo << " is: " << decProb << endl;	
						} */
						// Fill histogram
						decayHists[z][n][e]->SetBinContent(p+1,eo+1, static_cast<double>(decProb));
					} // end of output excited states
					sepEnergies[z][n][p]  = static_cast<double>(sepEnergy);
				} // end of particle type
				excEnergies[z][n][e] = static_cast<double>(excEnergy1);
				binW[z][n][e] = static_cast<double>(excBinWidth1);
				decayHists[z][n][e]->Write(); // save histograms
			} // end of excited states
		} // end of neutron number
	} // end of proton number
		
 	talysFile->Close();
 	return;
}

// Follows the decay chain of a given excited state (from Event), returns number of decay particles,
// adds decay particles to given eventi
//
// This is not the whole story, because TALYS stops when reaching ground state of unstable
// nuclei...
//
// Because of this I still need to add decays from ground states: TALYS database has some decay
// channels so I should work on accessing that as well
//
// For now I handle particle emission in the last part of this function manually, basically if
// separation energy is negative
//
//
unsigned int TalysData::decayParticles(Event& event) {
	if(xscn->decayFinal0==false || event.xscnExState==0) return 0;
	unsigned int pCtr = 0;	
	int zCur = 0; 
	int zReal = xscn->zFinalNuc;
	int nCur = 0; 
	int nReal = (xscn->aFinalNuc - xscn->zFinalNuc); 
	int exCur;
	double eDiff, eDiffMin = 100;
	// Find starting TALYS bin
	for (int eCtr = 0; eCtr < nExc; eCtr++) {
		eDiff = abs(excEnergies[zCur][nCur][eCtr] - xscn->excLevels[event.xscnExState].energyGnd);
		if( eDiff < eDiffMin ) {
			exCur = eCtr;
			eDiffMin = eDiff; 
		}
	}
	TFile* talysFile = new TFile(talysFileName,"READ");
	talysFile->cd();
	while(exCur>0 && zReal>0 && nReal>0) { // If nuclei is excited keep going
		double pTypeD, exOutD;
		TH2D* decayer = (TH2D*) talysFile->Get(getHName(zCur,nCur,exCur));
		decayer->GetRandom2(pTypeD,exOutD); // get decay mode
		//decayHists[zCur][nCur][exCur]->Print();//GetRandom2(pTypeD,exOutD); // Random decay
		int pType = round(pTypeD), exNext = round(exOutD); // Round to integer
		if(pType==0 && exNext==exCur) continue; // Do not allow decay to same bin...
		int zNext = zCur + ZRED[pType], nNext = nCur + NRED[pType];
		// Get Maximum and Minimum Possible Energies from Bin Widths and Get a random energy within
		double minEnergy = (excEnergies[zCur][nCur][exCur] - binW[zCur][nCur][exCur]) - 
			(excEnergies[zNext][nNext][exNext] + binW[zNext][nNext][exNext]) - 
			sepEnergies[zCur][nCur][pType];
		if( minEnergy<0 ) minEnergy = 0; // Cannot have negative KE
		double maxEnergy = (excEnergies[zCur][nCur][exCur] + binW[zCur][nCur][exCur]) - 
			(excEnergies[zNext][nNext][exNext] - binW[zNext][nNext][exNext]) - 
			sepEnergies[zCur][nCur][pType];
		double pEnergy = decayRand.Rndm()*(maxEnergy-minEnergy) + minEnergy;
		
		if ( pEnergy<=0 ) { // Warning if KE< 0...
			cout << " ***** WARNING: KINETIC ENERGY FOUND TO BE NEGATIVE ****************** " << endl;
			cout << "Z: " << zCur << " N: " << nCur << " Exc: " << exCur << endl; 
			cout << "Actual E: "  << xscn->excLevels[event.xscnExState].energyGnd <<
					" Bin E: " << excEnergies[zCur][nCur][exCur] << " Bin Width: " << 
					binW[zCur][nCur][exCur] << endl;
			cout << "Decay Bin E: " << excEnergies[zNext][nNext][exNext] << " Decay Bin Width: " <<
					binW[zNext][nNext][exNext] << endl;	
			cout << "Sel. Particle: " << pType << " Sel. Ex. State: " << exNext << endl;
			cout << " KE sel: " << pEnergy  << " KE max: " << maxEnergy << 
				" KE min: " << minEnergy << endl;
		}
		event.addParticle(PDGID[pType],pEnergy,0,1); // Add particle to the event at t=0, as KE
		// Update Current State
		zCur = zNext; nCur = nNext; exCur = exNext;
		zReal = zReal - ZRED[pType]; nReal = nReal - NRED[pType];
		pCtr++; // increase decay particle count
	}
	talysFile->Close();

	// TALYS stops when it reaches ground state of an unstable nuclei, manual particle decay
	
	while(true) {
		int zNext,nNext;
		if(zReal<=0 || nReal<= 0) break;	
		double maxSep = 0;
		int pType = 0;
		for (int pCtr = 0; pCtr < 7; pCtr++) {
			if(maxSep > sepEnergies[zCur][nCur][pCtr]) {
				maxSep = sepEnergies[zCur][nCur][pCtr];
				pType = pCtr;
			}
		}
		if (maxSep >= 0) break;
		double pEnergy = -1*maxSep;
		// Output to notify added decay
		cout << "Z: " << zCur << " N: " << nCur << endl; 
		cout << "Type: " << pType << " Energy: " << pEnergy << endl; 
		event.addParticle(PDGID[pType],pEnergy,0,1); // Add particle to the event at t=0, as KE
		// Update Current State
		zCur = zNext; nCur = nNext;
		zReal = zReal - ZRED[pType]; nReal = nReal - NRED[pType];
		pCtr++; // increase decay particle count	
	}	
	return pCtr;	
}

void TalysData::textOutput() { // Just to test accessibility of TALYS  results
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

TString TalysData::getHName(int z, int n, int exc) { // get a string for given values
	int zCur = xscn->zFinalNuc - z; 
	int nCur = (xscn->aFinalNuc - xscn->zFinalNuc) - n; 
	TString name = "h"+(TString)xscn->strName+Form("_z%d_n%d_exc%d",zCur,nCur,exc);
	return name;
}


