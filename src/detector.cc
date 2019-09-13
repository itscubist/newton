/********
*Author: Baran Bodur
*Date: 2019-09-12
*Description: Class That Reads Flux Info and Organizes Access To It 
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
#include "TVector3.h"
#include "TF1.h"
#include "TRandom3.h"

// Project Headers
#include "detector.h"

using namespace std;
using namespace TMath;

	// Constructor
Detector::Detector(std::string cardName, TFile* outFile) {
	ifstream detCard(cardName); // Read flux card file
	string line; // dummy to read lines in xscn card
	string name, value;
	unsigned int marker;
	while(!detCard.eof()){ // Until File Ends
		getline(detCard,line); // Get Line
		if(detCard.eof()) break; // Break if File Ended
		if(line[0]=='#') continue; // skip comments
		marker = line.find(' '); // Find Space, Where Value Starts
		value = line.substr(marker+1); // Get Value
		name = line.substr(0,marker); // Get Parameter Name
		cout << name << " - " << value <<  endl;
		if(name=="DET_NAME") strName = value;
		if(name=="NAVOGADRO") nAvogadro = stod(value); 	
		if(name=="MATERIAL") material = value; 	
		if(name=="MOLARMASS") molarMass = stod(value); 	
		if(name=="DENSITY") density = stod(value); 	
		if(name=="DETMASS") detMass = stod(value); 	
		if(name=="DETYEARS") detYears = stod(value); 	
		if(name=="SECONDSINYEAR") yearsToSeconds  = stod(value); 	
		if(name=="UNITFIXER") unitFixer = stod(value); 	
		if(name=="DETGEOM") detGeom = value; 		
		if(name=="DETX") detX = stod(value); 		
		if(name=="DETY") detY = stod(value); 		
		if(name=="DETZ") detZ = stod(value); 		
		if(name=="DETR") detR = stod(value); 		
	} // end of reading card file
	detCard.close(); 
	// Unit fixing * Time in Seconds * Number of Atoms in the Detector
	overallCoeff = unitFixer*(detYears*yearsToSeconds)*(nAvogadro*detMass/molarMass);
	randCyl = new TF1("randCyl","x*x",0,detR);
	randSph = new TF1("randSph","x",0,detR);
	detRand.SetSeed(0);
}

// Function To Give a Random Point in the Detector Considering Detector Geometry
// Assumes Detector Center is the Center of Coordinate System
TVector3 Detector::randomInteractionVertex() {
	TVector3 vertexVector;
	if(detGeom=="RECTANGULAR") {
		vertexVector.SetX(detRand.Uniform((-0.5)*detX,0.5*detX));
		vertexVector.SetY(detRand.Uniform((-0.5)*detY,0.5*detY));
		vertexVector.SetZ(detRand.Uniform((-0.5)*detZ,0.5*detZ));
	}
	else if(detGeom=="CYLINDRICAL") {
		double randR = randCyl->GetRandom();
		double randAzi = detRand.Uniform(0,2*Pi());
		vertexVector.SetX(randR*cos(randAzi));
		vertexVector.SetY(randR*sin(randAzi));
		vertexVector.SetZ(detRand.Uniform((-0.5)*detZ,0.5*detZ));
		//vertexVector.SetPhi(randAzi);
		//vertexVector.SetPerp(randR);
	}
	else if(detGeom=="SPHERICAL") {
		double randR = randSph->GetRandom();
		double randAzi = detRand.Uniform(0,2*Pi());
		double randCosZen = detRand.Uniform(-1.0,1.0);	
		double randSinZen = sqrt(1-randCosZen*randCosZen);	
		vertexVector.SetX(randR*randSinZen*cos(randAzi));
		vertexVector.SetX(randR*randSinZen*sin(randAzi));
		vertexVector.SetX(randR*randCosZen);	
		//vertexVector.SetMag(randR);	
		//vertexVector.SetTheta(acos(randCosZen));
		//vertexVector.SetPhi(randAzi);
	}
	else cout << "Geometry Not Found, Check Detector Card File" << endl;
	return vertexVector;
}
