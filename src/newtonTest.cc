/********
*Author: Baran Bodur
*Date:2019-09-10
*Description: Test file for earky features of NEWTON
*
********/

// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <queue>
#include <bitset>
// ROOT headers
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TF1.h"

// Project Headers
#include "xscn.h"
#include "flux.h"
#include "detector.h"
#include "event.h"
#include "exstates.h"

using namespace std;
using namespace TMath;


//Main Function
int main(int argc, const char *argv[]) {
	//****** Obtain Inputs
	if (argc<0) {
		cout <<"Usage: ./newtonTest "<< endl;
		return 0;
	} // end of usage printing statement
	
	// Outfile and detector
	TFile* outFile = new TFile("newtonTestOut.root","RECREATE");
	Detector superK("detData/cardSuperK.txt",outFile);
	// Xscns
	unsigned int nXscn = 7;
	Xscn ibdXscn("xscnData/cardIbd.txt",superK.material,outFile);
	Xscn nueOXscn("xscnData/cardNueO16.txt",superK.material,outFile);
	Xscn nuebarOXscn("xscnData/cardNuebarO16.txt",superK.material,outFile);
	Xscn nueEs("xscnData/cardNueEs.txt",superK.material,outFile);
	Xscn nuebarEs("xscnData/cardNuebarEs.txt",superK.material,outFile);
	Xscn nuxEs("xscnData/cardNuxEs.txt",superK.material,outFile);
	Xscn nuxbarEs("xscnData/cardNuxbarEs.txt",superK.material,outFile);
	
	Xscn *activeXscns[nXscn] = {&ibdXscn, &nueOXscn, &nuebarOXscn, 
		&nueEs,&nuebarEs,&nuxEs,&nuxbarEs};
	TString xscnNames[nXscn] = {"ibd","nueO16","nuebarO16","nueEs","nuebarEs","nuxEs","nuxbarEs"};

	// FLux
	Flux atmFlux("fluxData/cardAtmospheric.txt",outFile);
	// Convolve Flux and Xscns
	int xscnCounts[nXscn];
	for (int xCtr = 0; xCtr < nXscn; xCtr ++) {
		xscnCounts[xCtr] = 	activeXscns[xCtr]->generateEventCountPerEnergy(atmFlux, superK);
		cout << "Total " << xscnNames[xCtr] << " Events: " << xscnCounts[xCtr] << endl;
		activeXscns[xCtr]->xscnVsEnergy->Write();	
		activeXscns[xCtr]->xscnVsEnergyAngle[0]->Write();
		activeXscns[xCtr]->eventsVsEnergy->Write();
	}
	// Some Testing
	for (int i = 0; i < 6; i++) atmFlux.fluxVsEnergy[i]->Write();
	nueOXscn.excProbVsEnergy->Write();
	nuebarOXscn.excProbVsEnergy->Write();

	/*Event ev2(50,nuebarOXscn,atmFlux,superK);
	Event ev(50,nuebarOXscn,atmFlux,superK);
	ev.particles[0].direction.Print();
	cout << "Electron Energy: " << ev.particles[0].energy << endl;
*/
	vector<Event> genEvents;
	TString aziName = "hAzi_";
	TString zenName = "hZen_";
	TString ezName = "hEZ_ ";
	TH1D* hAzi[nXscn], *hZen[nXscn];
	for (unsigned int xCtr = 0; xCtr < nXscn; xCtr++) {
		cout << "Xscn: " << xscnNames[xCtr] << endl; // Print Xscn
		for (unsigned int eCtr = 0; eCtr < activeXscns[xCtr]->genEnergies.size(); eCtr++) {
			Event tempEvent(activeXscns[xCtr]->genEnergies[eCtr],*activeXscns[xCtr],atmFlux,superK);
			if(tempEvent.xscnExState>0) { // Add gamma
				tempEvent.addParticle(22,activeXscns[xCtr]->excLevels[tempEvent.xscnExState].energyGnd,0);	
			}
			genEvents.push_back(tempEvent);
			double tempCosZen = tempEvent.particles[0].direction.CosTheta();
			double tempAzi = tempEvent.particles[0].direction.Phi()*RadToDeg();
			if(tempAzi<0) tempAzi+=360.0;
			double tempEne = tempEvent.particles[0].energy;
			// Fill TH3D of electron
			activeXscns[xCtr]->leptonDirEnergyDist->Fill(tempCosZen,tempAzi,tempEne);
		}
		activeXscns[xCtr]->leptonDirEnergyDist->Write();
		hZen[xCtr] = activeXscns[xCtr]->leptonDirEnergyDist->
				ProjectionX(zenName + xscnNames[xCtr],1,12,30,100,"e");
		hZen[xCtr]->Write();		
		hAzi[xCtr] = activeXscns[xCtr]->leptonDirEnergyDist->
				ProjectionY(aziName+xscnNames[xCtr],10,11,30,100,"e");
		hAzi[xCtr]->Write();

	}
	// out Kinematic File
	ofstream outText("newtonTest.kin");
	for (int eCtr = 0; eCtr < genEvents.size(); eCtr++) genEvents[eCtr].writeEvent(outText);
	outText.close();
	
	// Testing
	//cout << "Overall Coefficient: " << superK.overallCoeff << endl;
	double cosAngle[10000], aveCos=0;
	TH1D* aziTestH = new TH1D("aziTestH","Azimuth",12,0,360);
	TH1D* zenTestH = new TH1D("zenTestH","Azimuth",20,-1,1);
	TH1D* zTestH = new TH1D("zTestH","z",100,-superK.detZ,superK.detZ);
	TH1D* xTestH = new TH1D("xTestH","x",100,-superK.detZ,superK.detZ);
	TH1D* yTestH = new TH1D("yTestH","y",100,-superK.detZ,superK.detZ);
	TH1D* rTestH = new TH1D("rTestH","r",100,0,superK.detR);
	for (int i = 0; i <10000; i++) {
		unsigned int randEx = nueOXscn.randomExLevelAtEnergy(100); 
		double exEnergy = nueOXscn.excLevels[randEx].energyInitial;
		int exSpin = nueOXscn.excLevels[randEx].spin;
		int exParity = nueOXscn.excLevels[randEx].parity;
		//cout << "Level: "<< randEx << " Energy: " << exEnergy << " Spin: " << exSpin << endl;

		TVector3 vertex = superK.randomInteractionVertex();
		//cout << "In loop: " << i << endl;
		zTestH->Fill(vertex.z());
		xTestH->Fill(vertex.x());
		yTestH->Fill(vertex.y());
		rTestH->Fill(vertex.Perp());

		zenTestH->Fill(atmFlux.randomZenithAtEnergy(100,0));
		aziTestH->Fill(atmFlux.randomAzimuthAtEnergyAndZenith(100,0.0,0));
		cosAngle[i] = ibdXscn.randomAngleAtEnergy(26.4,0);
		aveCos+=cosAngle[i];
	}
	aveCos/=10000.0;
	//cout << aveCos << endl;
	aziTestH->Write();
	zenTestH->Write();
	zTestH->Write();
	xTestH->Write();
	yTestH->Write();
	rTestH->Write();
	
	outFile->Close();
	return 0;
}
