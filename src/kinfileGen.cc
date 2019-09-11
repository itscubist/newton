/********
*Author: Baran Bodur
*Date: 2019-02-26
*Description: To write events to kinFiles, and read them back
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

using namespace std;

void kinAddTrack(ofstream &kinFile,int pdgId, float energy, float* pos, float* dir, float time) {
	kinFile << "$ begin" << endl;
	kinFile << "$ vertex " << pos[0] << " " << pos[1] << " " << pos[2] << " " << time << endl;
	kinFile << "$ track " <<pdgId<<" "<<energy<<" "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<" 0"<< endl;
	kinFile << "$ end" << endl;
}

int kinNextTrack(ifstream &kinFile, int &pdgId, float &energy, float* pos, float* dir, float &time) {
	string dummy[2];
	int dumInt;
	int ctr = 0;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1];
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1] >> pos[0] >> pos[1] >> pos[2] >> time;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1] >> pdgId >> energy >> dir[0] >> dir[1] >> dir[2] >> dumInt;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1];
	return 0;
}

int kinGetTrack(ifstream &kinFile, int trackNo, int &pdgId, float &energy, 
		float* pos, float* dir, float &time) {
	kinFile.seekg(0, kinFile.beg); // Go to beginning
	int skippedLines = 4*(trackNo-1); // Skip all the previous lines
	string dummy[2];
	for (int i = 0; i < skippedLines; i++) {
		getline(kinFile,dummy[0]);
		if(kinFile.eof()) return -1;
	}
	int dumInt;
	int ctr = 0;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1];
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1] >> pos[0] >> pos[1] >> pos[2] >> time;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1] >> pdgId >> energy >> dir[0] >> dir[1] >> dir[2] >> dumInt;
	if(kinFile.eof()) return -1;
	kinFile >> dummy[0] >> dummy[1];
	return 0;
}

