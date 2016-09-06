#ifndef ROOTPLOTTER_BARS_H
#define ROOTPLOTTER_BARS_H

#include <vector>
#include <string>

#include <TH2.h>
#include <TH1.h>
#include <TH1F.h>
#include <TStyle.h>
#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TVirtualHistPainter.h"
#include "TMath.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"

using namespace std;

//TCanvas *c1 = new TCanvas("c1","transparent pad",200,10,700,500);
//TPad *pad1 = new TPad("pad1","",0,0,1,1);

class DATASET
{
public:
	string name;
	vector<double> numbers;
	int size;

};

class SORTED
{
public:
	vector<double> X, Y, Z, E, L, B; 
	int size;
	double cut;
	string name;
};

class QSTATS
{
public:
	double E10_40, E10_30, E10_20, E20_40, E20_30, E30_40;
	int ROI;
	double cut;
};



#endif