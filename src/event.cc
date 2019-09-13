/********
*Author: Baran Bodur
*Date: 2019-09-13
*Description: Class That 
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
#include "event.h"

using namespace std;
using namespace TMath;

// Constructor
Event::Event(double inEnergy, Xscn &inXscn, Flux &inFlux, Detector &inDet) { 
	xscnName = inXscn.strName;
	fluxName = inFlux.strName;
	nuEnergy = inEnergy;
	// Get Random Vertex According to Detector Geometry
	vertex = inDet.randomInteractionVertex(); 
	time = 0.0; // time is always 0 for now...
	// Should fix this part later, for selecting a neutrino flavor
	for (int fCtr = 0; fCtr < 6; fCtr++) {
		nuFlv = fCtr; if(inXscn.intNu[fCtr]==true) break;
	}
	// Get Random angles, excited states, directions ...
	fluxCosZen = inFlux.randomZenithAtEnergy(nuEnergy,nuFlv); 
	fluxAzi = inFlux.randomAzimuthAtEnergyAndZenith(nuEnergy,fluxCosZen,nuFlv); 
	xscnExState = inXscn.randomExLevelAtEnergy(nuEnergy); 
	initialTargetMass = inXscn.mTarget;
	finalHadMass = inXscn.mFinal0 + inXscn.excLevels[0].energyGnd;
	xscnCosAngle = inXscn.randomAngleAtEnergy(nuEnergy, xscnExState);
	xscnAzi = inXscn.randomAzimuth();
	
	// Lepton From Scattering
	Particle chLepton;
	chLepton.pdgId = inXscn.pdgLepton;
	chLepton.mass = inXscn.mLepton;
	chLepton.isSimulated = true;
	particles.push_back(chLepton);
	
	// Hadron From Scattering
	Particle finalHadron;
	chLepton.pdgId = inXscn.pdgFinal0;
	chLepton.mass = finalHadMass; // with excited state
	chLepton.isSimulated = inXscn.writeFinal0;
	particles.push_back(finalHadron);

	// Function to fill energy and direction of the charged lepton
	fillLeptonDirAndEnergy();
}

// Functions To Write Event to Kin File
// void writeEvent(ofstream &outKinFile);

// Calculate Lepton Direction and Energy
void Event::fillLeptonDirAndEnergy() {
	
	// For Energy
	double A = -2*initialTargetMass+2*nuEnergy*(xscnCosAngle-1);
	double B = 2*initialTargetMass*nuEnergy+pow(initialTargetMass,2)-
		pow(finalHadMass,2)-pow(particles[0].mass,2);
	double delta = sqrt(pow(B,2)+4*pow(particles[0].mass,2)*A);
	double ePlus = (-B + delta)/(2*A);
	double eMinus = (-B - delta)/(2*A);
	double elEnergy = eMinus;
	if(elEnergy<particles[0].mass) particles[0].energy = particles[0].mass;
	else particles[0].energy = eMinus;

	// For direction
	particles[0].direction = combineZenAzi(fluxCosZen,fluxAzi,xscnCosAngle,xscnAzi);
}

//void Event::fillHadronDirAndEnergy() {
//
//}

// A handy kinematic function to add zenith and azimuth
TVector3 Event::combineZenAzi(double cosZen1,double azi1, double cosZen2, double azi2) {
	// Get sines and cosines
	double sinZen1 = sqrt(1-pow(cosZen1,2));	
	double sinZen2 = sqrt(1-pow(cosZen2,2));	
	double sinAzi1 = sin(DegToRad()*azi1);	
	double sinAzi2 = sin(DegToRad()*azi2);	
	double cosAzi1 = cos(DegToRad()*azi1);	
	double cosAzi2 = cos(DegToRad()*azi2);	
	// Construct coordinate vectors with z axis in neutrino direction
	TVector3 xPrime,yPrime,zPrime;
	zPrime.SetXYZ(sinZen1*cosAzi1,sinZen1*sinAzi1,cosZen1);
	if(zPrime(2)>=1) xPrime.SetXYZ(1,0,0);
	else if(zPrime(2)<=-1) xPrime.SetXYZ(-1,0,0);
	else xPrime.SetXYZ(cosZen1,-1*sinZen1*cosAzi1,0);
	xPrime.SetMag(1.0);
	yPrime = zPrime.Cross(xPrime);
	yPrime.SetMag(1.0);
	// Combine to get electron direction in detector frame
	TVector3 combinedDir;
	combinedDir = sinZen2*cosAzi2*xPrime + sinZen2*sinAzi2*yPrime + cosZen2*zPrime;
	combinedDir.SetMag(1.0);
	return combinedDir;
}
