// muonconvmc.C
// Marcos Flores
// 2019 July 25;
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
#define MUON_MASS 0.10565837 // GeV
#define e_MASS 5.109989e-4 // GeV

#define NDECAY 1e6 // Number of Decays computed

#define ZSi 14.0
#define ASi 28.0855 

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

Double_t calcW(UInt_t  Z, Double_t A, Double_t Egamma, Double_t xplus){
// Calculates the W parameter described in "Monte Carlo Generator for Muon
// Pair Production"

	Double_t B;
	Double_t Dn;

	if(Z == 1){
		B = 202.4;
		Dn = 1.49;
	}else{
		B = 183;
		Dn = 1.54*pow(A,0.27);
	}

	Double_t Winf = B*pow(Z,-1.0/3.0)/Dn * MUON_MASS/e_MASS;

	if(xplus == 1){
		return Winf;
	}else{
		Double_t delta = pow(MUON_MASS,2)/(2*Egamma*xplus*(1-xplus));
		Double_t W = Winf*(1 + (Dn*sqrt(TMath::E()) - 2)*delta/MUON_MASS)/(1 + B*pow(Z,-1.0/3.0)*sqrt(TMath::E())*delta/e_MASS);
		return W;
	}

}

vector<Double_t> muonprod(Double_t Egamma){

	Double_t xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
	Double_t xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

	Double_t xplus = xmin + gRandom->Uniform(0,1)*(xmax-xmin);

	Double_t Eplus = xplus*Egamma;
	Double_t Eminus = Egamma - Eplus;


	Double_t W = calcW(ZSi, ASi, Egamma, xplus);
	Double_t Winf = calcW(ZSi, ASi, Egamma, 1);

	Double_t nDiffCross = (1-(4.0/3.0)*xplus*(1-xplus))*TMath::Log(W)/TMath::Log(Winf);

	if(nDiffCross < 0){
		nDiffCross = 0;
	}

	vector<Double_t> outv;

	outv.push_back(nDiffCross);
	outv.push_back(Eplus);
	outv.push_back(Eminus);

	return outv;

}

bool muontest(Double_t Egamma, Double_t Eplus){

	Double_t xplus = Eplus/Egamma;

	Double_t W = calcW(ZSi, ASi, Egamma, xplus);
	Double_t Wmax = calcW(ZSi, ASi, Egamma, 0.5);

	Double_t tDiffCross = (1-(4.0/3.0)*xplus*(1-xplus))*TMath::Log(W)/TMath::Log(Wmax);
	
	if(tDiffCross > gRandom->Uniform(0,1)){
		return true;
	}else{
		return false;
	}

}

void muonconvmc(){

//--------------------------------------------------------------------
// Graphics
//--------------------------------------------------------------------

Int_t nbins = 100;

TCanvas *Canvas1 = new TCanvas("Canvas1","Photon Energy Dist.",0,100,600,500);
TCanvas *Canvas2 = new TCanvas("Canvas2","Pi Zero Hist",0,100,600,500);
TCanvas *Canvas3 = new TCanvas("Canvas3","Generated Cross Section",0,100,600,500);
TCanvas *Canvas4 = new TCanvas("Canvas4","Generated Cross Section",0,100,600,500);
//TCanvas *Canvas5 = new TCanvas("Canvas5","X+ Spectra",0,100,600,500);

TH1F *PhotonEHist = new TH1F("PhotonEHist","Photon Energy Dist.; Photon Energy [GeV]; Count [#]",nbins,0,25);

TH1F *AtlasH = new TH1F("AtlasH","MC Atlas Histogram ; pT [GeV]; Count [#]",nbins,0.5,50);

TH1F *CrossHist = new TH1F("CrossHist","Generated Cross Section; Value; Count [#]",nbins,0,1);

TH1F *TCrossHist = new TH1F("TCrossHist","Test Generated Cross Section; x+; Count [#]",nbins,0,1);

TH1F *XPSpecHist = new TH1F("XPSpecHist","x+ spectra; x+; Count [#]",nbins,0,1);

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

TF1 *func = new TF1("fitf",fitf,1e-3,100,4);
func->SetParameters(1,1,2,-1);
g1->Fit(func,"QN");

Double_t *params = func->GetParameters();
func->SetParameters(params);

// Normalize function
Double_t norm = func->Integral(1e-3,100);
params[0] = params[0]/norm;
func->SetParameters(params);

//--------------------------------------------------------------------
// Enter main generation
//--------------------------------------------------------------------

Int_t num = NDECAY; // Number of Decays

Double_t pizero_pt; // Random pi zero pT based on ATLAS distribution
Double_t eta;    	// random eta of pion (assumed flat)
Double_t phi; 	 	// random phi of pion (assumed flat)
bool difftest;

TLorentzVector pizero;
vector<Double_t> photovect;
vector<Double_t> outv;

Double_t nDiffCross;

for(UInt_t i = 0; i < NDECAY; ++i){

	// Get a pT for the neutral pion
	pizero_pt = func->GetRandom(6,100);

	AtlasH->Fill(pizero_pt);

	// Generate random eta & phi (these are assumed flat for now)
	eta = gRandom->Uniform(-1,1);
	phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());

	// Generate two photons provided the above parameters
	photovect = gammaprod(pizero_pt,eta,phi);

	// photovect.push_back(6.1);
	// photovect.push_back(6.1);

	for(UInt_t j = 0; j < 2; ++j){
		if(photovect.at(j) > 6){
			PhotonEHist->Fill(photovect.at(j));
			outv = muonprod(photovect.at(j));
			// std::cout << "\t Egamma: " << photovect.at(j) << "\tnDiffCross: " << nDiffCross << std::endl;
			CrossHist->Fill(outv.at(0));

			difftest = muontest(photovect.at(j),outv.at(1));
			XPSpecHist->Fill(outv.at(1)/photovect.at(j));

			if(difftest == true){

				TCrossHist->Fill(outv.at(1)/photovect.at(j));
			}
		}

	}

}

f->Close();

Canvas1->cd();
PhotonEHist->Draw();

Canvas2->cd();
AtlasH->Draw();

Canvas3->cd();
CrossHist->Draw();

Canvas4->cd();
TCrossHist->Draw();

// Canvas5->cd();
// XPSpecHist->Draw();

}