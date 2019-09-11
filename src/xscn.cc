/********
*Author: Baran Bodur
*Date: 2019-09-10
*Description: Class That Reads Xscn Data In and Allows Its Use 
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
#include "xscn.h"

using namespace std;

Xscn::Xscn(string cardName, string material) { // Constructor
	ifstream xscnCard(cardName); // Read xscn card file
	string line; // dummy to read lines in xscn card
	string name, value;
	unsigned int marker;

	while(!xscnCard.eof()){ // Until File Ends
		getline(xscnCard,line); // Get Line
		if(xscnCard.eof()) break; // Break if File Ended
		if(line[0]=='#') continue; // skip comments
		marker = line.find(' '); // Find Space, Where Value Starts
		value = line.substr(marker+1); // Get Value
		name = line.substr(0,marker); // Get Parameter Name
		cout << name << " - " << value <<  endl;
		if(name=="XSCN_NAME") strName = value;
		if(name=="WATERPERMOLECULE" && material=="WATER") targetPerMolecule = stod(value); 	
		if(name=="HEAVYWATERPERMOLECULE" && material=="HEAVYWATER") targetPerMolecule = stod(value); 	
		if(name=="FINALSTATES") {
			nFinalStates = stoi(value); 	
			strXscnVsEnergyAndAngle.reserve(nFinalStates);
		}
		if(name=="PDGT") pdgTarget = stoi(value); 	
		if(name=="PDGF") pdgFinal0 = stoi(value); 	
		if(name=="PDGL") pdgLepton = stoi(value); 	
		if(name=="MT") mTarget = stod(value); 	
		if(name=="MF") mFinal0 = stod(value); 	
		if(name=="ML") mLepton = stod(value); 	
		if(name=="WRITE_F") writeFinal0 = stoi(value); 	
		if(name=="DECAY_F") decayFinal0 = stoi(value); 		
		if(name=="NUE") intNu[0] = stoi(value); 	
		if(name=="NUEBAR") intNu[1] = stoi(value); 	
		if(name=="NUMU") intNu[2] = stoi(value); 	
		if(name=="NUMUBAR") intNu[3] = stoi(value); 	
		if(name=="NUTAU") intNu[4] = stoi(value); 	
		if(name=="NUTAUBAR") intNu[5] = stoi(value); 	
		if(name=="ANGBINS") angBins = stoi(value); 		
		if(name=="ENERGYBINS") energyBins = stoi(value); 		
		if(name=="DOUBLEDIFF") doubleDiff = stoi(value); 		
		if(name=="EXCDATA") excData = stoi(value); 		
		if(name=="XSCNDATA") strXscnVsEnergy = value; 		
		if(name=="EXCLEVELDATA") strExcStates = value; 		
		if(name=="EXCPROBDATA") strExcProb = value; 		
		if(name.substr(4)=="ANGXSCNDATA") {
			strXscnVsEnergyAndAngle.push_back = value; 		
		}
	}

	if(nFinalStates>1 && excData!=0) {
		
	}
}

TH1D* Xscn::readXscn() { // Read Single Diff Xscn
	ifstream inFile();
}

vector<TH2D*> Xscn::readXscnDouble() { // Read Double Diff Xscn Per Exc Energy

}

TH2D* Xscn::readExcProb() { // Read Exc Prob vs Energy

}

vector<exStates> Xscn::readExcStates() { // Read Exc State Info

}

