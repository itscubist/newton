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
			fluxCards.push_back(value); 		
		}	
		if(name.length()>9 && name.substr(3)=="XSCNCARD") {
			xscnCards.push_back(value); 		
		}
	} // End of file read
	
	// Construct root file, detector, xscns and fluxes based on card file input
	rootOutFile = new TFile((TString)rootFileName,"RECREATE"); 		
	detector = new Detector(detCard,rootOutFile);	
	cout << "# of Flux Cards: " << fluxCards.size() << " # of Xscn Cards: " << xscnCards.size() << endl;
	for (unsigned int i = 0; i < fluxCards.size(); i++) { // construct fluxes
		Flux* tempFlux = new Flux(fluxCards[i],rootOutFile);	
		fluxes.push_back(tempFlux);
		cout << "Flux Card: " << fluxCards[i] << " is created!" << endl;	
	}	
	for (unsigned int i = 0; i < xscnCards.size(); i++) { // construct xscns
		Xscn* tempXscn = new Xscn(xscnCards[i],detector->material,rootOutFile);	
		xscns.push_back(tempXscn);
		cout << "Xscn Card: " << xscnCards[i] << " is created!" << endl;	
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
	for (unsigned int fCtr = 0; fCtr < fluxes.size(); fCtr++) {
		vector<int> tempVecEv;
		vector<double> tempVecExp; 	
		cout << "******* Events From Flux: " << fluxes[fCtr]->strName << " **********" << endl;
		for (unsigned int xCtr = 0; xCtr < xscns.size(); xCtr++) {
			// Get number of events generated (and expected) for each energy bin
			int tempEv = xscns[xCtr]->generateEventCountPerEnergy(*fluxes[fCtr], *detector);
			double tempExp = xscns[xCtr]->expectedVsEnergy->Integral();	
			tempVecEv.push_back(tempEv); // event counts per flux per xscn
			tempVecExp.push_back(tempExp); // expected counts per flux per xscn
			// Print counts for user to see
			cout << "Total " << xscns[xCtr]->strName << " Events: " << tempEv << endl;
			cout << "Expected " << xscns[xCtr]->strName << " Events: " << tempExp << endl;
			for (unsigned int eCtr = 0; eCtr < xscns[xCtr]->genEnergies.size(); eCtr++) {
				Event tempEvent(xscns[xCtr]->genEnergies[eCtr],*xscns[xCtr],*fluxes[fCtr],*detector);
				// If can decay, and decay is requested
				if(xscns[xCtr]->nFinalStates>1 && xscns[xCtr]->decayFinal0==true) {
					xscns[xCtr]->talysDecayer.decayParticles(tempEvent);
				}
				events.push_back(tempEvent);
				double tempCosZen = tempEvent.particles[0].direction.CosTheta();
				double tempAzi = tempEvent.particles[0].direction.Phi()*RadToDeg();
				if(tempAzi<0) tempAzi+=360.0;
				double tempEne = tempEvent.particles[0].energy;
				// Fill TH3D of electron
				xscns[xCtr]->leptonDirEnergyDist->Fill(tempCosZen,tempAzi,tempEne);
				xscns[xCtr]->gammaEnergyHist->Fill(tempEvent.totalGammaEnergy);
				xscns[xCtr]->neutronNumberHist->Fill(tempEvent.totalNeutrons);
			}
		} // end of loop over xscns
		eventCounts.push_back(tempVecEv); // save event counts per flux
		expectedCounts.push_back(tempVecExp); // save expected counts per flux
	} // end of loop over fluxes
}

// Function saves events to tree and kin files with truth info
void Organizer::saveEvents() {
	// out Kinematic File
	ofstream outText(vectorFileName);
	ofstream outTrue(vectorFileName+".trueNu");
	rootOutFile->cd();
	//genTree = new TTree("genTree","tree containing all created event info");
  //genTree->Branch("xscnName",tempEvent,"eventNo/C");
	for (int eCtr = 0; eCtr < events.size(); eCtr++) {
		events[eCtr].writeEvent(outText);
		events[eCtr].writeNuInfo(outTrue);
	}
	outText.close();
	outTrue.close();
}

// Function saves diagnostic histograms into root file
void Organizer::plotHists() {
	TH1D* hAzi, *hZen, *hE;
	TH1D *aziTestH, *zenTestH;
	TString aziName = "hAzi_", zenName = "hZen_", eName = "hE_ ";
	rootOutFile->cd();
	// Write flux graphs
	for (unsigned int fCtr = 0; fCtr < fluxes.size(); fCtr++) {
		for (int i = 0; i < 6; i++) {
			fluxes[fCtr]->fluxVsEnergy[i]->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
			fluxes[fCtr]->fluxVsEnergy[i]->Write();
		}
		aziTestH = new TH1D("aziTestH_"+(TString)fluxes[fCtr]->strName,"Azimuth",12,0,360);
		zenTestH = new TH1D("zenTestH_"+(TString)fluxes[fCtr]->strName,"Cos Zenith",20,-1,1);
		for (int i = 0; i < 10000; i++) {
			zenTestH->Fill(fluxes[fCtr]->randomZenithAtEnergy(100,0));
			aziTestH->Fill(fluxes[fCtr]->randomAzimuthAtEnergyAndZenith(100,0.0,0));
		}
		zenTestH->Write();
		aziTestH->Write();
		delete zenTestH; delete aziTestH;
	}
	// Detector vertex test
	TH1D* zTestH = new TH1D("zTestH","z",100,-detector->detZ,detector->detZ);
	TH1D* xTestH = new TH1D("xTestH","x",100,-detector->detX,detector->detX);
	TH1D* yTestH = new TH1D("yTestH","y",100,-detector->detY,detector->detY);
	TH1D* rTestH = new TH1D("rTestH","r",100,0,detector->detR);
	TH1D* rSphTestH = new TH1D("rSphTestH","r",100,0,detector->detR);
	for (int i = 0; i < 10000; i++) {
		TVector3 vertex = detector->randomInteractionVertex();
		zTestH->Fill(vertex.z());
		xTestH->Fill(vertex.x());
		yTestH->Fill(vertex.y());
		rTestH->Fill(vertex.Perp());
		rSphTestH->Fill(vertex.Mag());
	}
	zTestH->Write();
	xTestH->Write();
	yTestH->Write();
	rTestH->Write();
	rSphTestH->Write();
	
	// Write xscn graphs
	for (unsigned int xCtr = 0; xCtr < xscns.size(); xCtr++) {
		if(xscns[xCtr]->nFinalStates>1) xscns[xCtr]->excProbVsEnergy->Write();
		xscns[xCtr]->xscnVsEnergy->SetLineWidth(2);
		xscns[xCtr]->xscnVsEnergy->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
		xscns[xCtr]->xscnVsEnergy->GetYaxis()
			->SetTitle("Cross Section (#10^{-38} #cm^{2})");
		xscns[xCtr]->xscnVsEnergy->Write();
		xscns[xCtr]->expectedVsEnergy->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
		xscns[xCtr]->expectedVsEnergy->Write();
		xscns[xCtr]->eventsVsEnergy->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
		xscns[xCtr]->eventsVsEnergy->Write();
		xscns[xCtr]->leptonDirEnergyDist->Write();
		xscns[xCtr]->gammaEnergyHist->Write();
		xscns[xCtr]->neutronNumberHist->Write();
		hZen = xscns[xCtr]->leptonDirEnergyDist->
				ProjectionX(zenName + (TString)xscns[xCtr]->strName,1,12,30,100,"e");
		//hZen->GetXaxis()->SetTitle("Lepton Energy (MeV)");
		//hZen->GetYaxis()->SetTitle("Events/22.5kTon/20years/cosBin");
		hZen->SetLineWidth(2);
		hZen->Write();		
		hAzi = xscns[xCtr]->leptonDirEnergyDist->
				ProjectionY(aziName+(TString)xscns[xCtr]->strName,5,15,30,100,"e");
		//hAzi->GetXaxis()->SetTitle("Lepton Energy (MeV)");
		//hAzi->GetYaxis()->SetTitle("Events/22.5kTon/20years/30 Degrees");
		hAzi->SetLineWidth(2);
		hAzi->Write();
		hE = xscns[xCtr]->leptonDirEnergyDist->
				ProjectionZ(eName+(TString)xscns[xCtr]->strName,0,20,0,12,"e");
		hE->GetXaxis()->SetTitle("Lepton Energy (MeV)");
		//hE->GetYaxis()->SetTitle("Events/22.5kTon/20years/MeV");
		hE->SetLineWidth(2);
		cout << "30-100 MeV Lepton Count in: " << xscns[xCtr]->strName << " is: " << 
			hE->Integral(30,100) << endl; 
		hE->Write();
		// Delete
		delete hE; delete hZen; delete hAzi;
	}
}

