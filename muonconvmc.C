// muonconvmc.C
// Marcos Flores
// 2019 July 19;
//
// Monte-Carlo for \pi_0\to 2\gamma
// 
// All units are GeV unless otherwise stated
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
#include "TLorentzVector.h"

#define PIZERO_MASS 0.134977 // GeV
#define NDECAY 5e6 // Number of Decays computed

Double_t fitf(Double_t *x, Double_t *par)
// Fit Function for the ATLAS Data
{
	Double_t arg = par[3] + par[2]*TMath::Log(x[0]);
	Double_t fitval = par[0]*pow(1+x[0],par[1])/pow(x[0],arg);
	return fitval;
}

vector<Double_t> gammaprod(Double_t pt,Double_t eta,Double_t phi){
// For a given pi-zero pT calculates two photon energies and outputs 
// them as a vector

	vector<Double_t> outv;

	Double_t Epi, ctheta;
	Double_t gamma, beta;
	Double_t Egamma1, Egamma2;

	TLorentzVector pizero;
	pizero.SetPtEtaPhiM(pt,eta,phi,PIZERO_MASS);

	Epi = pizero.E();
	ctheta = gRandom->Uniform(-1,1);
	gamma = Epi/PIZERO_MASS;
	beta = sqrt(1-1/pow(gamma,2));

	Egamma1 = gamma*PIZERO_MASS/2*(1 + beta*ctheta);
	Egamma2 = gamma*PIZERO_MASS/2*(1 - beta*ctheta);

	outv.push_back(Egamma1);
	outv.push_back(Egamma2);

	return outv;
}

void muonconvmc(){

//--------------------------------------------------------------------
// Graphics
//--------------------------------------------------------------------

Int_t nbins = 100;

TCanvas *Canvas1 = new TCanvas("Canvas1","Photon Energy Dist.",0,100,600,500);
// TCanvas *Canvas2 = new TCanvas("Canvas2","",0,100,600,500);

TH1F *PhotonEHist = new TH1F("PhotonEHist","Photon Energy Dist.; Photon Energy [GeV]; Count [#]",nbins,0,25);

TH1F *AtlasH = new TH1F("AtlasH","MC Atlas Histogram ; pT [GeV]; Count [#]",nbins,0.5,50);

//--------------------------------------------------------------------
// Fitting ATLAS Function
//--------------------------------------------------------------------

// Import data from ALTAS
 
// Open ATLAS File
TFile *f = 0;
 
f = new TFile("~/Desktop/HEPData-ins1419652-v1-Table_12.root");
if(!f) throw "Can't open file";
 
// Get graph from data
TGraphAsymmErrors *g1 =(TGraphAsymmErrors*)f->Get("Table 12/Graph1D_y1");
 
// Fit the graph
//	NOTE: The starting value of the function is chosen due to 
//		  the divergent nature at the origin

TF1 *func = new TF1("fitf",fitf,1e-3,60,4);
func->SetParameters(1,1,2,-1);
g1->Fit(func,"QN");

Double_t *params = func->GetParameters();
func->SetParameters(params);

// Normalize function
Double_t norm = func->Integral(1e-3,60);
params[0] = params[0]/norm;
func->SetParameters(params);

//--------------------------------------------------------------------
// Enter main generation

Int_t num = NDECAY; // Number of Decays

Double_t pizero_pt; // Random pi zero pT based on ATLAS distribution
Double_t eta;    	// random eta of pion (assumed flat)
Double_t phi; 	 	// random phi of pion (assumed flat)

TLorentzVector pizero;
vector<Double_t> photovect;

for(UInt_t i = 0; i < NDECAY; ++i){

	// Get a pT for the neutral pion
	pizero_pt = func->GetRandom();

	// Generate random eta & phi (these are assumed flat for now)
	eta = gRandom->Uniform(-1,1);
	phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());

	// Generate two photons provided the above parameters
	photovect = gammaprod(pizero_pt,eta,phi);	

	PhotonEHist->Fill(photovect.at(0));
	PhotonEHist->Fill(photovect.at(1));

}

Canvas1->SetLogy();

Canvas1->cd();
PhotonEHist->Draw();


f->Close();

}