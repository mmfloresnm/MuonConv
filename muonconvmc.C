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

#define NDECAY 1e7 // Number of Decays computed

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

	Double_t Epi, ctheta, gtheta, gphi;
	Double_t gamma, beta;
	Double_t Egamma1, Egamma2;

	TLorentzVector pizero;
	pizero.SetPtEtaPhiM(pt,eta,phi,PIZERO_MASS);

	TVector3 pizero3 = pizero.Vect();
	TVector3 pi_unitv = pizero3.Unit();
	TVector3 boostv =  -pizero.Beta()*pizero3.Unit();

	Epi = pizero.E();
	gamma = pizero.Gamma();
	beta = pizero.Beta();

	TLorentzVector gfourv1, gfourv2;

	ctheta = gRandom->Uniform(-1,1);
	gtheta = gRandom->Uniform(0,TMath::Pi());
	gphi = gRandom->Uniform(0,2*TMath::Pi());


	// Create forward photon in rest frame
	Double_t pg = PIZERO_MASS/2;

	gfourv1.SetE(pg);
	gfourv1.SetPx(pg*TMath::Sin(gtheta)*TMath::Cos(gphi));
	gfourv1.SetPy(pg*TMath::Sin(gtheta)*TMath::Sin(gphi));
	gfourv1.SetPz(pg*TMath::Cos(gtheta));

	// Create backward photon in rest frame

	gfourv2.SetE(pg);
	gfourv2.SetPx(-pg*TMath::Sin(gtheta)*TMath::Cos(gphi));
	gfourv2.SetPy(-pg*TMath::Sin(gtheta)*TMath::Sin(gphi));
	gfourv2.SetPz(-pg*TMath::Cos(gtheta));

	// Stuff for consisency tests -- should remove later

	TVector3 test1 = gfourv1.Vect();
	TVector3 test2 = gfourv2.Vect();

	Egamma1 = gamma*(PIZERO_MASS/2 - boostv*test1);
	Egamma2 = gamma*(PIZERO_MASS/2 - boostv*test2);

	// Boost back into the lab frame
	gfourv1.Boost(-boostv);
	gfourv2.Boost(-boostv);

	// std::cout  << "TLorentz: " << gfourv1.E() << "\tBy Hand E1: " << Egamma1 << std::endl;

	outv.push_back(gfourv1.E());
	outv.push_back(gfourv2.E());

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
// Generates plus/minus muon energies and the value of the differential cross section

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


Double_t cprefactor(Double_t Z, Double_t mass){
// Calculate the classical radius of a particle

	Double_t fstruc = 1.0/137.035999084;
	Double_t spdlgt = 299792458;
	Double_t echge = 1.602e-19;
	Double_t k = 8.99e9;

	Double_t kgmass = mass*pow(10,9)*echge/(spdlgt*spdlgt);

	Double_t rc = k*echge*echge/(kgmass*spdlgt*spdlgt);

	Double_t prefactor = 4.0 * fstruc * ZSi * ZSi * rc * rc;

	return prefactor;

}

Double_t cexpr(Double_t *x, Double_t *par){
// Expression for the differetial cross section

	Double_t xx = x[0];
	Double_t prefactor = cprefactor(ZSi,MUON_MASS);
	Double_t nCross = prefactor*(1 - 4.0/3.0*xx*(1-xx))*TMath::Log(calcW(ZSi,ASi,par[0],xx));
	return nCross;
}

bool convprob(Double_t Egamma){
// Determine whether a photo will pair produce

	Double_t sigmainf;
	Double_t crssec;

	Double_t xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
	Double_t xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

	TF1 *func = new TF1("cexpr",cexpr,xmin,xmax,1);
	func->SetParameter(0,Egamma);
	sigmainf = 7.0/9.0*cprefactor(ZSi,MUON_MASS)*TMath::Log(calcW(ZSi,ASi,Egamma,1));

	crssec = func->Integral(xmin,xmax)/sigmainf;

	if(gRandom->Uniform(0,1)< crssec){
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

TH1F *PhotonEHist = new TH1F("PhotonEHist","Photon Energy Dist.; Photon Energy [GeV]; Count [#]",nbins,6,25);

TH1F *AtlasH = new TH1F("AtlasH","MC Atlas Histogram ; pT [GeV]; Count [#]",nbins,0.5,50);

TH1F *CrossHist = new TH1F("CrossHist","Generated Cross Section; Value; Count [#]",nbins,0,1);

TH1F *TCrossHist = new TH1F("TCrossHist","Test Generated Cross Section; x+; Count [#]",nbins,0,1);

TH1F *XPSpecHist = new TH1F("XPSpecHist","x+ spectra; E(mu+); Count [#]",nbins,6,25);

TH1F *XPConvHist = new TH1F("XPConvHist","x+ spectra (with conversion calc); E(mu+); Count [#]",nbins,6,25);

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

	for(UInt_t j = 0; j < 2; ++j){
		if(photovect.at(j) > 6){
			PhotonEHist->Fill(photovect.at(j));
			outv = muonprod(photovect.at(j));
			XPSpecHist->Fill(outv.at(1));
			if(convprob(photovect.at(j)) == true){
				XPConvHist->Fill(outv.at(1));
			}

		}

	}

}

f->Close();

Canvas1->cd();
AtlasH->Draw();

Canvas2->cd();
PhotonEHist->Draw();

Canvas3->cd();
XPSpecHist->Draw();

Canvas4->cd();
XPConvHist->Draw();

// Canvas5->cd();
// XPSpecHist->Draw();

}