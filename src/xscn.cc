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

Xscn::Xscn(string cardName, string material, TFile* outFile) { // Constructor
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
			//strXscnVsEnergyAndAngle.reserve(nFinalStates);
		}
		if(name=="NUANCECODE") nuanceCode = stoi(value);
		if(name=="PDGT") pdgTarget = stoi(value); 	
		if(name=="PDGF") {
			pdgFinal0 = stoi(value);
			if(value.length() == 10) { // not nucleus size
				zFinalNuc = stoi(value.substr(3,3)); // get Z from pdg id
				aFinalNuc = stoi(value.substr(6,3)); // get A from pdg id
				cout << "Z nuc: " << zFinalNuc << " - A nuc: " << aFinalNuc << endl;
			}
			else {
				zFinalNuc = 0; aFinalNuc = 0;
			}
		} 	
		if(name=="PDGL") pdgLepton = stoi(value); 	
		if(name=="MT") mTarget = stod(value); 	
		if(name=="MF") mFinal0 = stod(value); 	
		if(name=="ML") mLepton = stod(value); 	
		if(name=="WRITE_F") writeFinal0 = stoi(value); 	
		if(name=="DECAY_F") decayFinal0 = stoi(value); 		
		// Read the next 6 in SNOWGLOBES order
		if(name=="NUE") intNu[0] = stoi(value); 	
		if(name=="NUEBAR") intNu[3] = stoi(value); 	
		if(name=="NUMU") intNu[1] = stoi(value); 	
		if(name=="NUMUBAR") intNu[4] = stoi(value); 	
		if(name=="NUTAU") intNu[2] = stoi(value); 	
		if(name=="NUTAUBAR") intNu[5] = stoi(value); 	
		if(name=="ANGBINS") angBins = stoi(value); 		
		if(name=="ENERGYBINS") {
			energyBins = stoi(value); 		
			energyBins2 = energyBins; // assume energyBins2 = energyBins if not overwritten later
		}
		if(name=="ENERGYBINS2") energyBins2 = stoi(value);
		if(name=="MAXNUENERGY") maxNuEnergy = stoi(value);
		if(name=="DOUBLEDIFF") doubleDiff = stoi(value); 		
		if(name=="EXCDATA") excData = stoi(value); 		
		if(name=="XSCNDATA") strXscnVsEnergy = value; 		
		if(name=="EXCLEVELDATA") strExcStates = value; 		
		if(name=="EXCPROBDATA") strExcProb = value; 		
		if(name.length()>10 && name.substr(3)=="ANGXSCNDATA") {
			strXscnVsEnergyAndAngle.reserve(1);
			strXscnVsEnergyAndAngle.push_back(value); 		
		}
	} // End of file read
	xscnCard.close();


	// Store Information In Histograms and Structs Accordingly
	outFile->cd();
	cout << "Reading Xscn Vs Energy" << endl;
	xscnVsEnergy = readXscn();
	cout << "Reading Double Differential Xscn For Each Excited State" << endl;
	xscnVsEnergyAngle = readXscnDouble();
	
	talysDecayer.xscn = this;
	if(nFinalStates>1 && excData!=0) {
		cout << "Reading Excited Levels" << endl;
		readExcStates();
		cout << "Reading Probabilities of Each Excited State As a Function of Energy" << endl;
		excProbVsEnergy = readExcProb();
	  if(decayFinal0 == true) { // If decay is simulated, then run TALYS to get decay probs.
			talysDecayer.runTalys(); // run TALYS
			//talysDecayer.textOutput(); // for testing access
			talysDecayer.setFileName("talysInterface/talysOutputs/"+(TString)strName+"_decay.root");
			talysDecayer.initHists(); // init histograms based on TALYS output size
		}
	}
	if(nFinalStates==1) { // Fill the ground state with the given mass in xscn card
		Exstates excLevel0;
		excLevel0.energyInitial = mFinal0-mTarget;
		excLevel0.energyGnd = 0;
		excLevels.push_back(excLevel0);
	}

	xscnRand.SetSeed(0);
	// Lepton Energy and Direction Distribution Filled
	leptonDirEnergyDist = new TH3D((TString)strName+"_leptonDir","Lepton Dir & Energy",
			20,-1,1,12,0,360,100,0,100);
	gammaEnergyHist = new TH1D((TString)strName+"_gammaE","Total Gamma Energy",50,0,50);
	neutronNumberHist = new TH1D((TString)strName+"_neutronN","Total Neutron Number",10,-0.5,9.5);
}

TH1D* Xscn::readXscn() { // Read Single Diff Xscn
	ifstream inFile(strXscnVsEnergy);
	string line;
	getline(inFile,line);
	TH1D* xscnH = new TH1D((TString)(strName)+"_xscnH",(TString)strXscnVsEnergy,
			energyBins,0.5,(double)(energyBins)+0.5);
	while(!inFile.eof()){
		double tEne,tXscn,tDdXscn;
		inFile >> tEne >> tXscn;
		if(inFile.eof()) break;
		xscnH->SetBinContent(tEne,tXscn);
	}
	inFile.close();
	return xscnH;
}

vector<TH2D*> Xscn::readXscnDouble() { // Read Double Diff Xscn Per Exc Energy
	vector<TH2D*> ddXscnVector;
	ifstream inFile;
	string line;
	for (int ex = 0; ex < nFinalStates; ex++) {
		inFile.open(strXscnVsEnergyAndAngle[ex]);
		getline(inFile,line);
		TH2D* ddXscnH = new TH2D(Form(((TString)strName)+"_angXscFsH_%d",ex), // Declare Hist
			(TString)strXscnVsEnergyAndAngle[ex],energyBins2,0.5,(double)(energyBins2)+0.5,angBins,-1,1);	
		while(!inFile.eof()){ // Loop to read file
			double tEne;
			inFile >> tEne;
			if(inFile.eof()) break;
			for (int ang = 0; ang < angBins; ang++) { // loop over angles
				double tDdXscn;
				inFile >> tDdXscn;
				ddXscnH->SetBinContent(tEne,ang+1,tDdXscn);
			} // end of reading 1 row
		} // end of reading 1 file
		ddXscnVector.reserve(1);
		ddXscnVector.push_back(ddXscnH);
		inFile.close();
	} // end of loop over different excited states
	return ddXscnVector;
}

TH2D* Xscn::readExcProb() { // Read Exc Prob vs Energy, Call only if nFinalStates>1
	ifstream inFile(strExcProb);
	string line;
	getline(inFile,line);
	TH2D* excProbH = new TH2D((TString)(strName)+"_excProbH",(TString)strExcProb, // Declare Hist
			energyBins2,0.5,(double)(energyBins2)+0.5,nFinalStates,-0.5,(double)(nFinalStates)-0.5); 
	while(!inFile.eof()){
		double tEne;
		inFile >> tEne;
		if(inFile.eof()) break;
		for (int ex = 0; ex < nFinalStates; ex++) { // loop over excited states
			double tExcProb;
			inFile >> tExcProb;
			excProbH->SetBinContent(tEne,ex+1,tExcProb);
		} // end of reading 1 row
	} // end of file
	inFile.close();
	return excProbH;
}

void Xscn::readExcStates() { // Read Exc State Info
	ifstream inFile(strExcStates);
	string line;
	for (int i = 0; i < 6; i++) getline(inFile,line); // Skip 6 Lines
	int exCtr = 0;
	while(!inFile.eof()){
		double tEx, dummy[3];
		int tSpin, tParity;
		inFile >> tEx >> tSpin >> tParity >> dummy[0] >> dummy[1] >> dummy[2];
		if(inFile.eof()) break;
		// Fill excited state object array
		Exstates tempState;
		excLevels.push_back(tempState);
		excLevels[exCtr].energyInitial = tEx;
		excLevels[exCtr].energyGnd = tEx - excLevels[0].energyInitial;
		excLevels[exCtr].spin = tSpin;
		excLevels[exCtr].parity = (tSpin==1) ? -1:1;
		exCtr++;
	} // Finish Reading File
	inFile.close();
	return;
}

// Functions To Get Xscn At Asked Values
double Xscn::xscnAtEnergy(double energy) { // Xscn at Eenrgy
	return xscnVsEnergy->Interpolate(energy);	
}

// This may not work in general, xscn at given angle, energy and excited state
double Xscn::xscnAtEnergyAngle(double energy, double cosAngle, unsigned int exLevel) { //Xscn At 
	return xscnVsEnergyAngle[exLevel]->Interpolate(energy,cosAngle);	
}

unsigned int Xscn::randomExLevelAtEnergy(double energy) { // Random Excited State At Sel. Energy 
	if(nFinalStates<2) return 0; // only excited state is the zeroth
	unsigned int selBin = excProbVsEnergy->GetXaxis()->FindBin(energy);
	// Not great but if energy is too high, so that no excitation data is available,
	// use the highest energy data avilable
	if(selBin>excProbVsEnergy->GetXaxis()->GetNbins()) 
		selBin = excProbVsEnergy->GetXaxis()->GetNbins();

	TH1D* excProbAtEnergy = excProbVsEnergy->ProjectionY("exP",selBin,selBin);
	unsigned int exLevel = round(excProbAtEnergy->GetRandom());
	delete excProbAtEnergy;
	if(exLevel==0) cout << "Ground State Interaction!" << endl; 
	return exLevel;
}

double Xscn::randomAngleAtEnergy(double energy, unsigned int exLevel) { // Select A Random Angle
	unsigned int selBin = xscnVsEnergyAngle[exLevel]->GetXaxis()->FindBin(energy);
	// Not great but if energy is too high, so that no excitation data is available,
	// use the highest energy data available
	if(selBin>xscnVsEnergyAngle[exLevel]->GetXaxis()->GetNbins()) 
		selBin = xscnVsEnergyAngle[exLevel]->GetXaxis()->GetNbins();
	TH1D* angleProbAtEnergy = xscnVsEnergyAngle[exLevel]->ProjectionY("exP",selBin,selBin);
	double cosAngle = angleProbAtEnergy->GetRandom();
	delete angleProbAtEnergy;
	return cosAngle;
}
	
double Xscn::randomAzimuth() { // Select A Random Angle, Uniform Accross Azimuth
	return xscnRand.Uniform(0,360);
}

// Function To Generate Event Counts with Given Flux
int Xscn::generateEventCountPerEnergy(Flux inFlux, Detector inDet) {
	int nBins = static_cast<int>(maxNuEnergy);
	eventsVsEnergy = new TH1D((TString)strName+"_evCnt","Poisson",nBins,0,maxNuEnergy);
	expectedVsEnergy = new TH1D((TString)strName+"_expCnt","Poisson",nBins,0,maxNuEnergy);
	//Loop over bins and fill them with poisson distributed event counts
	for (unsigned int bCtr = 1; bCtr <= eventsVsEnergy->GetXaxis()->GetNbins(); bCtr++) {
		double energy = eventsVsEnergy->GetBinCenter(bCtr);				
		// Flux Times Cross Section
		double poissonPar = 0;
		for (unsigned int fCtr = 0; fCtr <6; fCtr++) { // Sum over all Flavors That Do This
			if(intNu[fCtr]==true) poissonPar += inFlux.fluxAtEnergy(energy,fCtr)*xscnAtEnergy(energy);
			//cout <<"Energy Bin: " << bCtr <<" Flavor: " << fCtr << " Flux: " << 
			//	inFlux.fluxAtEnergy(energy, fCtr) << endl;
		}
		// Weigh By How Many Target and How Much Time
		poissonPar*=inDet.overallCoeff*targetPerMolecule;
		// Fill expected event count per energy histogram (aka. Poisson means)
		expectedVsEnergy->SetBinContent(bCtr,poissonPar);			
	}
		
	if(inDet.fixedEventCount<=0) { // If using # of events expected in given time/det volume...
		for (unsigned int bCtr = 1; bCtr <= eventsVsEnergy->GetXaxis()->GetNbins(); bCtr++) {
		// Get Event Count, Fill Histogram
		double poissonPar = expectedVsEnergy->GetBinContent(bCtr); // get poisson mean
		int evCount;
		if(inDet.fixedEventCount<0) // If smaller than zero add poisson fluctuations
			evCount = xscnRand.Poisson(poissonPar); // find poisson fluctuation
		else if(inDet.fixedEventCount==0) // If equal to zero then use closest integer
			evCount = round(poissonPar); // round to nearest int
		eventsVsEnergy->SetBinContent(bCtr,evCount); // set bin content to that fluctuation
		// Create Events
		if(evCount<1) continue;
		for (int eCtr = 0; eCtr < evCount; eCtr++) { //Loop over bins, and create events with that energy
			double evEnergy = xscnRand.Uniform(eventsVsEnergy->GetXaxis()->GetBinLowEdge(bCtr),
					eventsVsEnergy->GetXaxis()->GetBinLowEdge(bCtr+1));
			genEnergies.push_back(evEnergy);
			}
		}
	}

	else if(inDet.fixedEventCount>0) { // If creating fixed number of events
		for (unsigned int eCtr = 0; eCtr < inDet.fixedEventCount ; eCtr++) {
			double evEnergy = expectedVsEnergy->GetRandom(); // get random energy for each event 
			eventsVsEnergy->Fill(evEnergy); // fill event histogram with the obtained random energy
			genEnergies.push_back(evEnergy); // to later create an event with that energy
		}
	}
	
	return eventsVsEnergy->Integral();
}

/*
// Destructor
Xscn::~Xscn() {
	delete xscnVsEnergy;
	delete excProbVsEnergy;
	delete eventsVsEnergy;
	delete leptonDirEnergyDist;
	delete gammaEnergyHist;
	delete neutronNumberHist;
	for(unsigned int v=0;v<xscnVsEnergyAngle.size();v++) {
		delete xscnVsEnergyAngle[v];
	}
}
*/
