/********
*Author: Baran Bodur
*Date: 2019-12-06
*Description: Main program to run NEWTON
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
#include "organizer.h"

using namespace std;
using namespace TMath;

//Main Function
int main(int argc, const char *argv[]) {
	//****** Obtain Inputs
	if (argc<2) {
		cout <<"Usage: ./newtonTest [main card]"<< endl;
		return 0;
	} // end of usage printing statement


	Organizer orgMain((string)argv[1]);	
	cout << "Read Input From Card Files, Initialized" << endl; 
	orgMain.generateEvents(); // generates events from input det,xscn and flux info
	cout << "Generated Events" << endl;
	if( !orgMain.batchMode) { // if batch mode is on events are saved on the go
		orgMain.saveEvents(); // saves events to ttree in root file and a kinematics file
		cout << "Saved Events" << endl;
	}
	orgMain.plotHists(); // plot histograms
	cout << "Plotted Histograms to Check" << endl;
/*
	// To test rotator
	double zen[5] = {0.9,0.45,0,-0.45,-0.9};
	double azi[4] = {0,90,180,270};
	for (int z1 = 0; z1 <5; z1++) {
		for (int z2 = 0; z2 <5; z2++) {
			for (int a1 = 0; a1 <4; a1++) {
				for (int a2 = 0; a2 <4; a2++) {
					cout << "Zen1: " << zen[z1] << " Zen2: " << zen[z2] 
						<< " Azi1: " << azi[a1] << " Azi2: " << azi[a2] << endl;
					TVector3 resNormal = orgMain.events[0].combineZenAzi(zen[z1],azi[a1],zen[z2],azi[a2]);
					TVector3 resNew = 
						orgMain.events[0].combineZenAziDot(zen[z1],azi[a1],zen[z2],azi[a2]);
					cout << "My old method" << endl;
					resNormal.Print(); 
					cout << "With new robust formula" << endl;
					resNew.Print();
					//cout << "Normal: " <<  << " Rodrigues': " << zen[z2] << endl;
				}
			}
		}
	}
*/
	return 0;
}
