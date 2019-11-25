/********
*Author: Baran Bodur
*Date: 2019-11-25
*Description: Class that is responsible for overall organization of the generator 
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
#include "TFile.h"

// Project Headers
#include "organizer.h"

using namespace std;
using namespace TMath;

Organizer::Organizer(string cardName) {
	ifstream orgCard(cardName); // Read organizer card file
	string line; // dummy to read lines in xscn card
	string name, value;
	unsigned int marker;
	while(!orgCard.eof()){ // Until File Ends
		getline(orgCard,line); // Get Line
		if(orgCard.eof()) break; // Break if File Ended
		if(line[0]=='#') continue; // skip comments
		marker = line.find(' '); // Find Space, Where Value Starts
		value = line.substr(marker+1); // Get Value
		name = line.substr(0,marker); // Get Parameter Name
		cout << name << " - " << value <<  endl;
		if(name=="ORG_NAME") strName = value;
		if(name=="ROOTFILENAME") rootFileName = value;	
		if(name=="VECTORFILENAME") vectorFileName = value;
		if(name=="DETCARD") detCard = value;
		if(name.length()>9 && name.substr(3)=="FLUXCARD") {
			fluxCards.reserve(1);
			fluxCards.push_back(value); 		
		}	
		if(name.length()>9 && name.substr(3)=="XSCNCARD") {
			xscnCards.reserve(1);
			xscnCards.push_back(value); 		
		}
	} // End of file read
	
	// Construct root file, detector, xscns and fluxes based on card file input
	rootOutFile = new TFile(rootFileName,"RECREATE"); 		
	detector = new Detector(detCard,rootOutFile);	
	for (unsigned int i = 0; i < fluxCards.size(); i++) { // construct fluxes
		fluxes.reserve(1);
		Flux* tempFlux = new Flux(fluxCards[i],rootOutFile);	
		fluxes.push_back(tempFlux);
	}
	for (unsigned int i = 0; i < xscnCards.size(); i++) { // construct xscns
		Xscns.reserve(1);
		Xscn* tempXscn = new Xscn(xscnCards[i],detector->material,rootOutFile);	
		xscns.push_back(tempXscn);
	}
	orgCard.close();	
	orgRand.SetSeed(0);
}


Organizer::~Organizer() {
	rootOutFile->Close();
	delete rootOutFile;
};

// Function that generates events for input xscn and fluxes
void Organizer::generateEvents() {
	for (int fCtr = 0; fCtr < fluxes.size(); fCtr++) {
		vector<int> tempVecEv, tempVecExp; 	
		cout << "******* Events From Flux: " << fluxes[fCtr]->strName << " **********" << endl;
		for (int xCtr = 0; xCtr < xscns.size(); xCtr++) {
			// Get number of events generated (and expected) for each energy bin
			int tempEv = xscns[xCtr]->generateEventCountPerEnergy(*fluxes[fCtr], *detector);
			int tempExp = xscns[xCtr]->expectedVsEnergy->Integral();	
			tempVecEv.push_back(tempEv); // event counts per flux per xscn
			tempVecExp.push_back(tempExp); // expected counts per flux per xscn
			// Print counts for user to see
			cout << "Total " << xscns[xCtr]->strName << " Events: " << temp << endl;
			cout << "Expected " << xscns[xCtr]->strName << " Events: " << tempExp << endl;
			

		} // end of loop over xscns
		eventCounts.push_back(tempVecEv); // save event counts per flux
		expectedCounts.push_back(tempVecExp); // save expected counts per flux
	} // end of loop over fluxes
}

// Function saves events to tree and kin files with truth info
void Organizer::saveEvents() {

}

// Function saves diagnostic histograms into root file
void Organizer::plotHists() {
}

