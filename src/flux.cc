/********
*Author: Baran Bodur
*Date: 2019-09-10
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

// Project Headers
#include "flux.h"

using namespace std;
using namespace TMath;

Flux::Flux(string cardName, TFile* outFile) {
	ifstream fluxCard(cardName); // Read flux card file
	string line; // dummy to read lines in xscn card
	string name, value;
	unsigned int marker;
	while(!fluxCard.eof()){ // Until File Ends
		getline(fluxCard,line); // Get Line
		if(fluxCard.eof()) break; // Break if File Ended
		if(line[0]=='#') continue; // skip comments
		marker = line.find(' '); // Find Space, Where Value Starts
		value = line.substr(marker+1); // Get Value
		name = line.substr(0,marker); // Get Parameter Name
		cout << name << " - " << value <<  endl;
		if(name=="FLUX_NAME") strName = value;
		if(name=="SCALE") scale = stod(value); 	
		if(name=="ESCALE") energyScale = stod(value); 	
		if(name=="STORETYPE") storageType = value; 	
		if(name=="DIRTYPE") dirType = value; 	
		if(name=="NUE") nuFluxes[0] = stod(value); 	
		if(name=="NUEBAR") nuFluxes[1] = stod(value); 	
		if(name=="NUMU") nuFluxes[2] = stod(value); 	
		if(name=="NUMUBAR") nuFluxes[3] = stod(value); 	
		if(name=="NUTAU") nuFluxes[4] = stod(value); 	
		if(name=="NUTAUBAR") nuFluxes[5] = stod(value); 	
		if(name=="COSTHETA") cosAngleRead = stod(value); 	
		if(name=="AZIMUTH") azimuthRead = stod(value); 	
		if(name=="ENERGY") energyRange[0] = stod(value); 	
		if(name=="ENERGYMAX") energyRange[1] = stod(value); 	
		if(name=="ANGBINS") angBins = stoi(value); 		
		if(name=="AZIBINS") aziBins = stoi(value); 		
		if(name=="ENERGYBINS") energyBins = stoi(value); 		
		if(name=="ROOTFILENAME") strFluxFileName = value;
		if(name=="H_NUE") strFluxHistName.push_back(value); 	
		if(name=="H_NUEBAR") strFluxHistName.push_back(value); 	
		if(name=="H_NUMU") strFluxHistName.push_back(value);	
		if(name=="H_NUMUBAR") strFluxHistName.push_back(value);
		if(name=="H_NUTAU") strFluxHistName.push_back(value);
		if(name=="H_NUTAUBAR") strFluxHistName.push_back(value);
	} // End of file read
	fluxCard.close();	
	// Store Information In Histograms and Structs Accordingly
	if(storageType=="rootFile") { // If Info in Root File
		TFile* inRootFile = new TFile((TString)strFluxFileName,"READ");
		//if(dirType=="rootFile") { // If Angular Info Also in Root File
			// Then Read TH3D Histograms with given names in card file
			for (int flv = 0; flv < 6; flv++) {
				// Read TH3D
				TString histName = (TString)strFluxHistName[flv];
				TH3D* tempRead = (TH3D*) inRootFile->Get(histName);
				outFile->cd(); // To put temp copy in the outfile rather than inFile
				TH3D* tempCopy = (TH3D*) tempRead->Clone();
				fluxVsEnergyDirectional.push_back(tempCopy);
				fluxVsEnergyDirectional[flv]->Scale(scale); // Scale for units to get nu/cm2/s/MeV/sr
				//Project TH1D for people who care about energy only
				TH1D* temp1D = (TH1D*) fluxVsEnergyDirectional[flv]->ProjectionZ(
						(TString)strFluxHistName[flv] + "_1D",1,angBins,1,aziBins);
				fluxVsEnergy.push_back(temp1D);
				// Scale projection to units of nu/cm2/s/MeV
				fluxVsEnergy[flv]->Scale(4*Pi()/(angBins*aziBins));
			} // end of loop over flavors	
		//} // end of "if" angular dist also given in root files
		inRootFile->Close();
	} // end of root file reading
	fluxRand.SetSeed(0);
}

// Builds single energy flux from card file
vector<TH3D*> Flux::buildSingleEnergyFlux() { 
	return fluxVsEnergyDirectional;
}

// Builds constant flux accross energies from card file
vector<TH3D*> Flux::buildConstantAlongEnergyFlux() {
	return fluxVsEnergyDirectional;
}

// total flux at given energy for a flavor
double Flux::fluxAtEnergy(double energy, unsigned int flavor) { //total flux at energy
	energy/=energyScale;
	return fluxVsEnergy[flavor]->Interpolate(energy);
}

// total flux at a given energy, zenith and azimuth for a flavor
double Flux::fluxAtDirection(double energy, double cosAngle, double azimuth, unsigned int flavor) { 
	//energy/=energyScale;
	return fluxVsEnergyDirectional[flavor]->Interpolate(cosAngle,azimuth,energy);
}

// random zenith angle at a given energy for a flavor
double Flux::randomZenithAtEnergy(double energy, unsigned int flavor) {
	double cosAngle;
	if(dirType=="rootFile") {
		energy/=energyScale;
		unsigned int selBinZ = fluxVsEnergyDirectional[flavor]->GetZaxis()->FindBin(energy);
		TH1D* zenithProbAtEnergy = fluxVsEnergyDirectional[flavor]->ProjectionX("zenP",1,aziBins,
				selBinZ,selBinZ);
		//zenithProbAtEnergy->Print("all");
		cosAngle = zenithProbAtEnergy->GetRandom();
		delete zenithProbAtEnergy;
	}
	else if(dirType == "singleAngle") cosAngle = cosAngleRead;
	else if(dirType == "isotropic") cosAngle = fluxRand.Uniform(-1.0,1.0);
	//cout << cosAngle << endl;
	return cosAngle;
}

// random azimuth angle at a given energy and zenith for a flavor
double Flux::randomAzimuthAtEnergyAndZenith(double energy, double cosAngle, unsigned int flavor) {
	double azimuth;
	if(dirType=="rootFile") {
		energy/=energyScale;
		unsigned int selBinZ = fluxVsEnergyDirectional[flavor]->GetZaxis()->FindBin(energy);
		unsigned int selBinX = fluxVsEnergyDirectional[flavor]->GetXaxis()->FindBin(cosAngle);
		TH1D* aziProbAtEnergyZen = fluxVsEnergyDirectional[flavor]->ProjectionY("aziP",selBinX,selBinX,
				selBinZ,selBinZ);
		azimuth = aziProbAtEnergyZen->GetRandom();
		delete aziProbAtEnergyZen;
	}
	else if(dirType == "singleAngle") azimuth = azimuthRead;
	else if(dirType == "isotropic") azimuth = fluxRand.Uniform(0,360);
	return azimuth;

}
