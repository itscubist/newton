/********
*Author: Baran Bodur
*Date: 2019-11-25
*Description: To organize all generator operations
* like which flux, xscn, detector info to use,
* which channels are active, where to save the output...
* which outputs to save into root or txt files
********/
#ifndef ORGANIZER_INCLUDED
#define ORGANIZER_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TString.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// Project headers
#include "flux.h"
#include "detector.h"
#include "xscn.h"
#include "event.h"


class Organizer{
public:
	// Constructor
	Organizer(std::string cardName);
	~Organizer();

	// Functions
	void generateEvents(); // generates events from input det,xscn and flux info
	void saveEvents(); // saves events to ttree in root file and a kinematics file
	void plotHists(); // plot histograms

	// Variables Filled From Card File
	std::string strName; // name of organizer
	std::string detCard; //
	std::vector<std::string> xscnCards; // 
	std::vector<std::string> fluxCards; //
	std::string rootFileName; // output Root file name
	std::string vectorFileName; // output Vector file name

	// Created objects to run the generator as inputed
	Detector* detector;
	std::vector<Xscn*> xscns;
	std::vector<Flux*> fluxes;
	// Created events
	std::vector<Event> events;
	// Created root output file
	TFile* rootOutFile;
	// Created root tree (1 entry per event)
	TTree* genTree;
	// Stores event counts for flux and 	
	std::vector<std::vector<int> > eventCounts;
	std::vector<std::vector<double> > expectedCounts;


	// Created debug histograms
	//

private:
	// RNG
	TRandom3 orgRand;
};


#endif
