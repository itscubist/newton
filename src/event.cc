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
	
	// Lepton From Scattering - particle[0]
	Particle chLepton;
	chLepton.pdgId = inXscn.pdgLepton;
	chLepton.mass = inXscn.mLepton;
	chLepton.simDegree = 1; // primary particle
	chLepton.time = 0; // Created at vertex
	chLepton.isSimulated = true;
	particles.push_back(chLepton);
	
	// Hadron From Scattering - this is particle[1]
	Particle finalHadron;
	finalHadron.pdgId = inXscn.pdgFinal0;
	finalHadron.mass = finalHadMass; // with excited state
	finalHadron.simDegree = 1;
	finalHadron.time = 0;
	finalHadron.isSimulated = inXscn.writeFinal0; // read whether simulated from card file
	particles.push_back(finalHadron);
	
	// Function to fill energy and direction of the charged lepton
	fillLeptonDirAndEnergy();
	// Fill hadron energy and direction if the xscn card asks for hadron output
	fillHadronDirAndEnergy();
}

// Functions To Write Event to Kin File
// void writeEvent(ofstream &outKinFile);

// Calculate Lepton Direction and Energy
void Event::fillLeptonDirAndEnergy() {
	if (initialTargetMass==0.511 && particles[1].pdgId==0) { // If the target is electron 
		// Then energy is from nue-e elastic scattering
		double kappa = pow((particles[0].mass+nuEnergy),2) - pow((xscnCosAngle*nuEnergy),2); 
		double top = 2*particles[0].mass*pow((nuEnergy*xscnCosAngle),2);
		particles[0].energy = particles[0].mass + top/kappa;
	}
	else if(particles[1].mass!=0 && particles[1].pdgId!=0) { // If the target is hadron
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
	}

	// Momentum
	particles[0].momentum = sqrt(pow(particles[0].energy,2)-pow(particles[0].mass,2));
	// For direction
	particles[0].direction = combineZenAzi(fluxCosZen,fluxAzi,xscnCosAngle,xscnAzi);
	return;
}

void Event::fillHadronDirAndEnergy() {
	// If no hadron, or not simulated
	if(particles[1].mass==0 || particles[1].isSimulated==false) return; 	
	double pHz = nuEnergy - xscnCosAngle*particles[0].momentum;	
	double pHx = -sqrt(1-pow(xscnCosAngle,2))*particles[0].momentum;	
	particles[1].momentum = sqrt(pow(pHz,2)+pow(pHx,2));	
	particles[1].energy = sqrt(pow(pHz,2)+pow(pHx,2)+pow(particles[1].mass,2));	
	double cosHadron = pHz/particles[1].momentum; // cos(ThetaNuHadron)
	double aziHadron = xscnAzi+180; // Azimuth is 180 degree other way
	// Check whether azimuth is within 0-360
	if(aziHadron>=360) aziHadron-=360; 
	else if(aziHadron<=0) aziHadron+=360;
	// Get Direction
	particles[1].direction = combineZenAzi(fluxCosZen,fluxAzi,cosHadron,aziHadron); 
	return;
}

// Add gamma to a specific direction, at a specific time if so desired
void Event::addParticle(int pdgId, double energy, TVector3 direction, double time) {
	Particle decayP;
	decayP.pdgId = pdgId;
	// Get masses from pdg id in MeV
	if(pdgId==22) decayP.mass = 0;
	if(pdgId==11 || pdgId==-11) decayP.mass = 0.511; // electron or positron
	if(pdgId==2112) decayP.mass = 939.565; // neutron
	if(pdgId==2212) decayP.mass = 938.272; // proton
	if(pdgId==1000010020) decayP.mass = 1876.122; // deuteron
	if(pdgId==1000010030) decayP.mass = 2809.432; // triton
	if(pdgId==1000020030) decayP.mass = 2809.413; // He-3
	if(pdgId==1000020040) decayP.mass = 3728.401; // alpha
	decayP.energy = energy;
	decayP.momentum = sqrt(pow(decayP.energy,2)-pow(decayP.mass,2));
	decayP.time = time;
	decayP.isSimulated = true;
	decayP.direction = direction; 
	decayP.simDegree = 2; 
	particles.push_back(decayP);
	return;
}

// Adds gamma with given energy isotropically
void Event::addParticle(int pdgId, double energy, double time) {
	TVector3 dir = isotropicDirection();
	addParticle(pdgId,energy,dir,time); // call the function with direction setting
	return;
}

// Writes event to kinFile in NUANCE format
void Event::writeEvent(ofstream &outFile) {
	double curTime = particles[0].time;
	outFile << "$ begin" << endl;
	// Initial vertex
	outFile << "$ vertex " << vertex.X() <<" "<< vertex.Y() <<" "<< vertex.Z() <<" "<< curTime << endl;	
	for (int pCtr = 0; pCtr < particles.size(); pCtr++) { // loop over particles
		if(particles[pCtr].isSimulated==false) continue; // if particle is not simulated skip
		if(particles[pCtr].time!=curTime) { // Then need to add another vertex with new time
			curTime = particles[pCtr].time;
			// Rewrite vertex with new time
			outFile<<"$ vertex "<<vertex.X()<<" "<< vertex.Y() <<" "<< vertex.Z() <<" "<< curTime << endl;	
		}
		// Add particle properties
		outFile << "$ track " << particles[pCtr].pdgId <<" " << particles[pCtr].energy <<" " <<  
			particles[pCtr].direction.X() <<" "<< particles[pCtr].direction.Y() <<" "<< 
			particles[pCtr].direction.Z() << " 0." << endl;	
	}
	outFile << "$ end" << endl;
	return;
}

// Gives a unit magnitude TVector in a random direction
TVector3 Event::isotropicDirection() {
	TVector3 dir(0,0,1.0);
	TF1 fCosT("fCosT","1",-1.0,1.0);
	TF1 fAzi("fAzi","1",0.0,360.0);
	double randCosTheta = fCosT.GetRandom();
	double randAzi = fAzi.GetRandom();
	dir.SetMagThetaPhi(1.0,ACos(randCosTheta),randAzi);
	return dir;
}

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
