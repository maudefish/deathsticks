#define rootplotter_c
#include <string>
#include <iostream>		//cin & cout
#include <sstream>		//two-way string stream
#include <fstream>		//file input & output
#include <vector>
#include <cmath>
#include <TH2.h>
#include <TH1.h>
#include <TH1F.h>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF2.h"
#include "TLegend.h"
#include "TMultiGraph.h"

#include "extract.h"

#include "rootplotter_bars.h"

using namespace std;

int rootplotter_bars()
{	
	vector<double> q1, q2, q3, q4, q5, q6, delta1, delta2, delta3, delta4, delta5, delta6;
	q1 = extract("output/q_e10_40");
	q2 = extract("output/q_e20_40");
	q3 = extract("output/q_e30_40");
	q4 = extract("output/q_e10_20");
	q5 = extract("output/q_e10_30");
	q6 = extract("output/q_e20_30");

	delta1 = extract("output/delta_e10_40");
	delta2 = extract("output/delta_e20_40");
	delta3 = extract("output/delta_e30_40");
	delta4 = extract("output/delta_e10_20");
	delta5 = extract("output/delta_e10_30");
	delta6 = extract("output/delta_e20_30");



	Color_t color1 = kRed+4;
	Color_t color2 = kRed+2;
	Color_t color3 = kViolet+9;
	Color_t color4 = kRed;
	Color_t color5 = kGreen;

	Color_t color6 = kMagenta;

	Size_t size1 = .75;
	Size_t size2 = 1;
	Size_t size3 = 1;
	Size_t size4 = 1;
	Size_t size5 = 1.75;

	Size_t size6 = 1.75;

	Style_t style1 = 5;
	Style_t style2 = 22;
	Style_t style3 = 8;
	Style_t style4 = 21;
	Style_t style5 = 29;

	Style_t style6 = 29;

	char title1[7] = "10 GeV";
	char title2[7] = "20 GeV";
	char title3[7] = "30 GeV";
	char title4[7] = "40 GeV";
	char title5[7] = "50 GeV";
	char title6[8] = "60+ GeV";

	//TCanvas *c3 = new TCanvas("c3","transparent pad",150,10,900,700);
	//c3->SetGrid();
	//TPad *pad3 = new TPad("pad3","",0,0,1,1);
	//TMultiGraph *mg = new TMultiGraph("mg","Statistical Prevalence of Extragalactic Magnetic Field Helicity");

	Double_t E10_40 [30];
	Double_t E10_30 [30];
	Double_t E10_20 [30];
	Double_t E20_40 [30];
	Double_t E20_30 [30];
	Double_t E30_40 [30];

	Double_t DEL10_40 [30];
	Double_t DEL10_30 [30];
	Double_t DEL10_20 [30];
	Double_t DEL20_40 [30];
	Double_t DEL20_30 [30];
	Double_t DEL30_40 [30];
	Double_t counter [30];

	for (int i=0; i<30; i++)
	{
		E10_40[i] = q1[i];
		E20_40[i] = q2[i];
		E30_40[i] = q3[i];
		E10_20[i] = q4[i];
		E10_30[i] = q5[i];
		E20_30[i] = q6[i];
		

		DEL10_40[i] = 45*delta1[i];
		DEL20_40[i] = 45*delta2[i];
		DEL30_40[i] = 45*delta3[i];
		DEL10_20[i] = 45*delta4[i];
		DEL10_30[i] = 45*delta5[i];
		DEL20_30[i] = 45*delta6[i];
		
		counter[i] = i+1;
	
	}

	TCanvas *can = new TCanvas("can", "ok", 150,100,1200,1000);
	//can->Draw();
	can->SetGrid();
	can->Divide(3,2,0.0001,0.0001,0);
	TPad* can_1 = (TPad*)(can->GetPrimitive("can_1"));
	TPad* can_2 = (TPad*)(can->GetPrimitive("can_2"));
	TPad* can_3 = (TPad*)(can->GetPrimitive("can_3"));
	TPad* can_4 = (TPad*)(can->GetPrimitive("can_4"));
	TPad* can_5 = (TPad*)(can->GetPrimitive("can_5"));
	TPad* can_6 = (TPad*)(can->GetPrimitive("can_6"));

	can_1->cd();
	TGraphErrors *g1 = new TGraphErrors(30,counter,E10_40, 0, DEL10_40);
	g1->SetMarkerColor(color1);
   	g1->SetMarkerSize(size1);
	g1->SetMarkerStyle(style1);
	g1->SetTitle("{10,40} GeV");
	g1->GetXaxis()->SetTitle("Radius (deg)");
  	g1->GetYaxis()->SetTitle("Q(R)");
  	g1->GetYaxis()->SetTitleOffset(1.5);
  	g1->GetXaxis()->CenterTitle();
   	g1->GetYaxis()->CenterTitle();
	g1->Draw();


	can_2->cd();
	TGraphErrors *g2 = new TGraphErrors(30,counter,E20_40, 0, DEL20_40);
	g2->SetMarkerColor(color2);
   	g2->SetMarkerSize(size2);
	g2->SetMarkerStyle(style2);
	g2->SetTitle("{20,40} GeV");
	g2->GetXaxis()->SetTitle("Radius (deg)");
  	g2->GetYaxis()->SetTitle("Q(R)");
  	g2->GetYaxis()->SetTitleOffset(1.5);
  	g2->GetXaxis()->CenterTitle();
   	g2->GetYaxis()->CenterTitle();
   	g2->Draw();

	can_3->cd();
	TGraphErrors *g3 = new TGraphErrors(30,counter,E30_40, 0, DEL30_40);
	g3->SetMarkerColor(color3);
   	g3->SetMarkerSize(size3);
	g3->SetMarkerStyle(style3);
	g3->SetTitle("{30,40} GeV");
	g3->GetXaxis()->SetTitle("Radius (deg)");
  	g3->GetYaxis()->SetTitle("Q(R)");
  	g3->GetYaxis()->SetTitleOffset(1.5);
  	g3->GetXaxis()->CenterTitle();
   	g3->GetYaxis()->CenterTitle();
	g3->Draw();

	can_4->cd();
	TGraphErrors *g4 = new TGraphErrors(30,counter,E10_20, 0, DEL10_20);
	g4->SetMarkerColor(color4);
   	g4->SetMarkerSize(size4);
	g4->SetMarkerStyle(style4);
	g4->SetTitle("{10,20} GeV");
	g4->GetXaxis()->SetTitle("Radius (deg)");
  	g4->GetYaxis()->SetTitle("Q(R)");
  	g4->GetYaxis()->SetTitleOffset(1.5);
  	g4->GetXaxis()->CenterTitle();
   	g4->GetYaxis()->CenterTitle();
	g4->Draw();

	can_5->cd();
	TGraphErrors *g5 = new TGraphErrors(30,counter,E10_30, 0, DEL10_30);
	g5->SetMarkerColor(color5);
   	g5->SetMarkerSize(size5);
	g5->SetMarkerStyle(style5);
	g5->SetTitle("{10,30} GeV");
	g5->GetXaxis()->SetTitle("Radius (deg)");
  	g5->GetYaxis()->SetTitle("Q(R)");
  	g5->GetYaxis()->SetTitleOffset(1.5);
  	g5->GetXaxis()->CenterTitle();
   	g5->GetYaxis()->CenterTitle();
	g5->Draw();

	can_6->cd();
	TGraphErrors *g6 = new TGraphErrors(30,counter,E20_30, 0, DEL20_30);
	g6->SetMarkerColor(color6);
   	g6->SetMarkerSize(size6);
	g6->SetMarkerStyle(style6);
	g6->SetTitle("{20,30} GeV");
	g6->GetXaxis()->SetTitle("Radius (deg)");
  	g6->GetYaxis()->SetTitle("Q(R)");
  	g6->GetYaxis()->SetTitleOffset(1.5);
  	g6->GetXaxis()->CenterTitle();
   	g6->GetYaxis()->CenterTitle();
	g6->Draw();

	TLegend *leg3 = new TLegend(0.75,0.4,0.97,0.1);
   	leg3->SetHeader("Photon Energy Parameters");
    leg3->AddEntry(g1,"(E1, E2) = 10, 40","P");
    leg3->AddEntry(g2,"(E1, E2) = 20, 40","P");
    leg3->AddEntry(g3,"(E1, E2) = 30, 40","P");
    leg3->AddEntry(g4,"(E1, E2) = 10, 20","P");
    leg3->AddEntry(g5,"(E1, E2) = 10, 30","P");
    leg3->AddEntry(g6,"(E1, E2) = 20, 30","P");

   

	//leg3->Draw();
	//c3->Update();


   	cout << endl << "should be done..." << endl;

	cout << "\n===============================================================================\n";
	string filename = "rootplotter.txt";

	cout << "\n===============================================================================\n\n";
	rename("CP_project_v4_TEMP.txt","CP_project_v4.txt");			//write over old txt file for future use
	remove("rootplotter_bars.txt");
	remove("rootplotter_bars_C_ACLiC_dict_rdict.pcm");
	remove("rootplotter_bars_C.d");
	remove("rootplotter_bars_C.so");


return 0;
}
