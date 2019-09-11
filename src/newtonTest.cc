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
#include "TF1.h"

// Project Headers
//#include "xscn.h"

using namespace std;
using namespace TMath;

//Main Function
int main(int argc, const char *argv[]) {
	//****** Obtain Inputs
	if (argc<0) {
		cout <<"Usage: ./newtonTest "<< endl;
		return 0;
	} // end of usage printing statement

	TFile* outFile = new TFile("newtonTestOut.root","RECREATE");
	ifstream testIbd("xscnData/ibdxscn.txt");
	ifstream testIbdDd("xscnData/ibdxscnDoubleDiff.txt");
	string line;
	vector<double> energy;
	vector<double> xscn;
	vector<vector<double>> ddXscn;
	TH1D* xscnH = new TH1D("xscnH","xscnH",200,0.5,200.5);
	TH2D* ddXscnH = new TH2D("ddXscnH","ddXscnH",200,0.5,200.5,20,-1,1);

	getline(testIbd,line);
	getline(testIbdDd,line);
	while(!testIbd.eof()){
		double tEne,tXscn,tDdXscn;
		testIbd >> tEne >> tXscn;
		if(testIbd.eof()) break;
		xscnH->Fill(tEne,tXscn);
		testIbdDd >> tEne;
		cout << tEne;
		//vector<double> tempV(20);
		for (int i = 0; i < 20; i++) {
			testIbdDd >> tDdXscn;
			//tempV.push_back(tDdXscn);
			ddXscnH->Fill(tEne,-0.95+0.1*i,tDdXscn);
			cout << " " << tDdXscn;
		}
		cout << endl;
		//ddXscn.push_back(tempV);
			
	}
	cout << ddXscnH->GetBinContent(199,19) << endl;
	cout << ddXscnH->Interpolate(1.5,0.85) << endl;
	xscnH->Write();
	ddXscnH->GetZaxis()->SetRangeUser(1e-8,1e-1);
	ddXscnH->Write();
	outFile->Close();
	return 0;
}
