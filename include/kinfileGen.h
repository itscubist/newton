/********
*Author: Baran Bodur
*Date: 	2019-02-26
*Description: To write events to kinFiles, and read them back
*
********/
#ifndef KINFILEGEN_INCLUDED
#define KINFILEGEN_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TString.h"


// Write atrack with specific properties to a kinfile
void kinAddTrack(ofstream &kinFile, int pdgId, float energy, float* pos, float* dir, float time);
// Read next track from a kinfile (Assuming you read tracks sqeuentially)
int kinNextTrack(ifstream &kinFile, int &pdgId, float &energy, float* pos, float* dir, float &time);
// Read the Nth (trackNo) track from a kinfile
int kinGetTrack(ifstream &kinFile, int trackNo, int &pdgId, float &energy, 
		float* pos, float* dir, float &time);
#endif
