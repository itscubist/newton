/********
*Author: Baran Bodur
*Date: 2019-09-13
*Description: Class handling how an event looks? 
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
	nuanceCode = inXscn.nuanceCode;
	// Get Random Vertex According to Detector Geometry
	vertex = inDet.getInteractionVertex(); 
	// For selecting a neutrino flavor (should make this a function in future!)
	selectNuFlv(inXscn, inFlux);

	// Get Random angles, excited states, directions ...
	fluxCosZen = inFlux.randomZenithAtEnergy(nuEnergy,nuFlv); 
	fluxAzi = inFlux.randomAzimuthAtEnergyAndZenith(nuEnergy,fluxCosZen,nuFlv); 
	nuDir = zenAziToDir(fluxCosZen, fluxAzi, true);
	xscnExState = inXscn.randomExLevelAtEnergy(nuEnergy); 
	initialTargetMass = inXscn.mTarget;
	initialTargetPdg = inXscn.pdgTarget; // get pdg of target
	// Get mass of the final state (need to calculate energy available to the electron)
	finalHadMass = inXscn.mFinal0 + inXscn.excLevels[xscnExState].energyGnd;
	//cout << "Ex State: " << xscnExState << " Mass of Final Hadron: " << finalHadMass << endl; 
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
	
	//Initialize event totals
	totalNeutrons = 0;
	// If final hadron is a neutron increase neutron number
	if(finalHadron.pdgId == 2112) totalNeutrons++;
	totalGammaEnergy = 0;

}

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
	particles[0].direction = combineZenAziDot(fluxCosZen,fluxAzi,xscnCosAngle,xscnAzi);
	return;
}

void Event::fillHadronDirAndEnergy() {
	// If no hadron, or not simulated
	if(particles[1].mass==0) return; 	
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
	particles[1].direction = combineZenAziDot(fluxCosZen,fluxAzi,cosHadron,aziHadron); 
	return;
}

// Add gamma to a specific direction, at a specific time if so desired
// if isKE = 1 then energy input is considered kinetic energy
void Event::addParticle(int pdgId, double energy, TVector3 direction, double time, bool isKE) {
	
	Particle decayP;
	decayP.pdgId = pdgId;
	// Get masses from pdg id in MeV
	if(pdgId==22) decayP.mass = 0;
	else if(pdgId==11 || pdgId==-11) decayP.mass = 0.511; // electron or positron
	else if(pdgId==2112) decayP.mass = 939.565; // neutron
	else if(pdgId==2212) decayP.mass = 938.272; // proton
	else if(pdgId==1000010020) decayP.mass = 1876.122; // deuteron
	else if(pdgId==1000010030) decayP.mass = 2809.432; // triton
	else if(pdgId==1000020030) decayP.mass = 2809.413; // He-3
	else if(pdgId==1000020040) decayP.mass = 3728.401; // alpha
	else if(pdgId==1000080150) decayP.mass = 13975,266; // 15O (unstable so decays might matter)
	else if(pdgId==1000070160) decayP.mass = 14909.588; // 16N (unstable so decays might matter)	
	else decayP.mass = 1000000; // other nuclei, (no need to write the mass for every possibility) 
	
	decayP.energy = energy; // if input is total energy
	if(isKE) decayP.energy += decayP.mass; // if input is kinetic energy
	decayP.momentum = sqrt(pow(decayP.energy,2)-pow(decayP.mass,2));
	decayP.time = time;
	decayP.isSimulated = true; // Decide on simulating or not...
	//Particles cannot be simulated (and will not create Cherenkov light)
	if( pdgId > 1000000000 ) decayP.isSimulated = false;
	decayP.direction = direction; 
	decayP.simDegree = 2; 
	
	
	// To store total gamma energy and total neutron number, observables
	if(pdgId==2112) totalNeutrons++; // increase total neutron number
	if(pdgId==22) totalGammaEnergy+=decayP.energy; // neutron

	particles.push_back(decayP); // add particle
	return;
}

// Selects neutrino flavor (required for interactions that look the same with two different
// flavors. Ex: nutau-e, numu-e scattering. This function assigns those interactions a neutrino
// flavor proportional to the incoming flux (so if no nu-tau's in the flux then it would be all
// nu-mu...) For interactions with single nu flavor this will return that flavor obviously
void Event::selectNuFlv(Xscn &inXscn, Flux &inFlux) {
	vector<int> nonZeroNu; 
	vector<double> nonZeroFlux;
	double totalFlux=0;
	nuFlv = -1;
	for (int fCtr = 0; fCtr < 6; fCtr++) { // loop over potential flavours
		if(inXscn.intNu[fCtr]==true) { // if event can be generated from that type of neutrino
			nonZeroNu.push_back(fCtr); // record which index
			nonZeroFlux.push_back(inFlux.fluxAtEnergy(nuEnergy,fCtr)); // record flux
			totalFlux += inFlux.fluxAtEnergy(nuEnergy,fCtr); // record total flux
		}
	}
	if(nonZeroNu.size()==1) { // If only 1 type of neutrino can do this interaction no need to worry
		nuFlv = nonZeroNu[0];
		nuPdg = NU_PDG[nuFlv]; 
	}
	else if(nonZeroNu.size()>1) { // If multiple types, ratio according to flux
		int tempRand = inXscn.xscnRand.Rndm(); // random to select from fluxes
		double cumFlux = 0;
		for (int fCtr = 0; fCtr < nonZeroNu.size(); fCtr++) {
			cumFlux+=nonZeroFlux[fCtr]; // add flux of next flavor to try again
			if(tempRand < cumFlux/totalFlux) {nuFlv = nonZeroNu[fCtr]; break;} // if lesser, then that flv
		}
		if (nuFlv==-1) nuFlv = nonZeroNu[nonZeroNu.size()-1]; // in case of numerical issues
		nuPdg = NU_PDG[nuFlv]; // pdg id
	} // end of neutrino flavor selection
	return;
}

// Adds gamma with given energy isotropically
void Event::addParticle(int pdgId, double energy, double time, bool isKE) {
	TVector3 dir = isotropicDirection();
	addParticle(pdgId,energy,dir,time,isKE); // call the function with direction setting
	return;
}

// Writes event to kinFile in NUANCE format
void Event::writeEvent(ofstream &outFile, bool nuanceType) {
	double curTime = particles[0].time;
	outFile << "$ begin" << endl;
	if(nuanceType==true) // if reaction type
		outFile << "$ nuance " << nuanceCode << endl; 
	// Initial vertex
	outFile<<"$ vertex "<< vertex.X() <<" "<< vertex.Y() <<" "<< vertex.Z() <<" "<< curTime << endl;	
	// True nu info
	if(nuanceType==true) {// if reaction type then also add initial neutrino/target info
		outFile << "$ track " << nuPdg << " " << nuEnergy << " " << nuDir.X() << " " 
			<< nuDir.Y() << " " << nuDir.Z() << " -1" << endl;
		outFile << "$ track " << initialTargetPdg << " " << initialTargetMass << " " << 0 << " " 
			<< 0 << " " << 0 << " -1" << endl;
	}	
	// All outgoing particles for MC truth information
	// first and second entries will be the outgoing lepton and target as above
	// third is the lepton, fourth is the nuclei after interaction
	// from there on it is products of the decay of final nuclei
	for (int pCtr = 0; pCtr < particles.size(); pCtr++) { // loop over particles
		// Add particle properties
		outFile << "$ track " << particles[pCtr].pdgId <<" " << particles[pCtr].energy <<" " <<  
			particles[pCtr].direction.X() <<" "<< particles[pCtr].direction.Y() <<" "<< 
			particles[pCtr].direction.Z() << " -1" << endl;	
	}

	// Only particles that are simulated 
	for (int pCtr = 0; pCtr < particles.size(); pCtr++) { // loop over particles
		if(particles[pCtr].isSimulated==false) continue; // if particle is not simulated skip
		if(particles[pCtr].time!=curTime) { // Then need to add another vertex with new time
			curTime = particles[pCtr].time;
			// Rewrite vertex with new time
			outFile<<"$ vertex "<<vertex.X()<<" "<< vertex.Y()<<" "<< vertex.Z() <<" "<<curTime << endl;	
		}
		// Add particle properties
		outFile << "$ track " << particles[pCtr].pdgId <<" " << particles[pCtr].energy <<" " <<  
			particles[pCtr].direction.X() <<" "<< particles[pCtr].direction.Y() <<" "<< 
			particles[pCtr].direction.Z() << " 0" << endl;	
	}
	outFile << "$ end" << endl;
	return;
}

// Writes true info: Neutrino Energy, Neutrino Direction
void Event::writeNuInfo(ofstream &outFile) {
	outFile << "$ begin" << endl;
	outFile << "$ track " << 
		nuEnergy << " " << nuDir.X() << " " << nuDir.Y() << " " << nuDir.Z() << endl;
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

// A handy kinematic function to add zenith and azimuth, and older version and not very clear
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

// A handy kinematic function to add zenith and azimuth, current version with dot products
// USE THIS ONE!!!
TVector3 Event::combineZenAziDot(double cosZen1,double azi1, double cosZen2, double azi2) {
	// Get sines and cosines (should not need them anymore)
	
	//double sinZen1 = sqrt(1-pow(cosZen1,2));	
	//double sinZen2 = sqrt(1-pow(cosZen2,2));	
	//double sinAzi1 = sin(DegToRad()*azi1);	
	//double sinAzi2 = sin(DegToRad()*azi2);	
	//double cosAzi1 = cos(DegToRad()*azi1);	
	//double cosAzi2 = cos(DegToRad()*azi2);	
	
	// Construct coordinate vectors with z axis in neutrino direction
	TVector3 xPrime,yPrime,zPrime;
	TVector3 xOri(1,0,0),yOri(0,1,0),zOri(0,0,1);
	TVector3 x,y,z;
	//zPrime.SetXYZ(sinZen1*cosAzi1,sinZen1*sinAzi1,cosZen1); // use ready function below instead
	zPrime = zenAziToDir(cosZen1,azi1,true);
	zPrime.SetMag(1.0);
	if(zPrime(0)==0) xPrime.SetXYZ(1,0,0); // by default x axis is prependicular...
	else xPrime = yOri.Cross(zPrime);
	xPrime.SetMag(1.0);
	yPrime = zPrime.Cross(xPrime);
	yPrime.SetMag(1.0);
	// Now invert since we need original system in terms of neutrino coordinates
	x.SetXYZ(xOri.Dot(xPrime),xOri.Dot(yPrime),xOri.Dot(zPrime));	
	y.SetXYZ(yOri.Dot(xPrime),yOri.Dot(yPrime),yOri.Dot(zPrime));	
	z.SetXYZ(zOri.Dot(xPrime),zOri.Dot(yPrime),zOri.Dot(zPrime));	
	x.SetMag(1.0); y.SetMag(1.0); z.SetMag(1.0);
	// Combine to get electron direction in detector frame
	TVector3 combinedDir, inNuDir;
	//inNuDir.SetXYZ(sinZen2*cosAzi2,sinZen2*sinAzi2,cosZen2); //use ready function below instead
	inNuDir = zenAziToDir(cosZen2,azi2,false);
	combinedDir.SetXYZ(x.Dot(inNuDir),y.Dot(inNuDir),z.Dot(inNuDir));
	combinedDir.SetMag(1.0);
	return combinedDir;
}

// A handy kinematic function to combine zenith and azimuth
TVector3 Event::combineZenAziRodrigues(double cosZen1,double azi1, double cosZen2, double azi2) {
	// Get sines and cosines
	double sinZen1 = sqrt(1-pow(cosZen1,2));	
	double sinZen2 = sqrt(1-pow(cosZen2,2));	
	double sinAzi1 = sin(DegToRad()*azi1);	
	double sinAzi2 = sin(DegToRad()*azi2);	
	double cosAzi1 = cos(DegToRad()*azi1);	
	double cosAzi2 = cos(DegToRad()*azi2);	
	// Construct coordinate vectors with z axis in neutrino direction
	TVector3 z, zPrime, rotAxis, oriDir, combinedDir;
	z.SetXYZ(0,0,1);
	zPrime.SetXYZ(sinZen1*cosAzi1,sinZen1*sinAzi1,cosZen1);
	if(zPrime(2)>=1) combinedDir.SetXYZ(sinZen2*cosAzi2,sinZen2*sinAzi2,cosZen2);
	else if(zPrime(2) <= -1) combinedDir.SetXYZ(-1*sinZen2*cosAzi2,-1*sinZen2*sinAzi2,-1*cosZen2);
	else {
		rotAxis = zPrime.Cross(z);
		rotAxis.SetMag(1.0);
		oriDir.SetXYZ(sinZen2*cosAzi2,sinZen2*sinAzi2,cosZen2);
		oriDir.SetMag(1.0);
		// Combine to get electron direction in detector frame via Rodrigues' formula
		combinedDir=oriDir*cosZen1+rotAxis.Cross(oriDir)*sinZen1+rotAxis*rotAxis.Dot(oriDir)*(1-cosZen1);
		combinedDir.SetMag(1.0);
	}
	return combinedDir;
}

// Convert zenith and azimuth to direction vector 
// In case of flux instead of neutrino direction, the sky direction from which neutrino is
// coming might be specified (which is negative of the neutrino direction), to do that still
// input the same zenith and azimuth values, but set isFlux=true (default false)
TVector3 Event::zenAziToDir(double cosZen, double azi, bool isFlux) {
	// Get sines and cosines
	double sinZen = sqrt(1-pow(cosZen,2));	
	double sinAzi = sin(DegToRad()*azi);	
	double cosAzi = cos(DegToRad()*azi);	
	TVector3 dir;
	dir.SetXYZ(sinZen*cosAzi,sinZen*sinAzi,cosZen);
	if(isFlux==true) dir=-1*dir;
	return dir;
}

// Destructor (commented, use default)
//Event::~Event() {
//	cout << "An event is destroyed" << endl;
//}
