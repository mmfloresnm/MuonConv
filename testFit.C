// testFit.C
// Marcos Flores
// 2019 July 16;
//
// Testing bed for fitting date from ATLAS paper
//
// -------------------------------------------------------------------

#include <string>
#include <Riostream.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>


Double_t fitf(Double_t *x, Double_t *par)
{
	Double_t arg = par[3] + par[2]*TMath::Log(x[0]);
	Double_t fitval = par[0]*pow(1+x[0],par[1])/pow(x[0],arg);
	return fitval;
}

void testFit(){

	TFile *f = 0;

	f = new TFile("~/Desktop/HEPData-ins1419652-v1-Table_12.root");
	if(!f) throw "Can't open file";

	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",0,100,600,500);

	Canvas1->SetLogx();
	Canvas1->SetLogy();

	TGraphAsymmErrors *g1 =(TGraphAsymmErrors*)f->Get("Table 12/Graph1D_y1");

	TF1 *func = new TF1("fitf",fitf,0.5,50,4);
	func->SetParameters(1,1,2,-1);

    g1->Fit(func,"MRE");
    std::cout << "\nChiSquare / NDoF: " << func->GetChisquare() << " / " << func->GetNDF() << "\n\n";

	g1->Draw();

	f->Close();
}
