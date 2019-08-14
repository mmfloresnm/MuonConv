// muonconvmc.C
// Marcos Flores
// 2019 August 12;
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

#define ZEl 14.0
#define AEl 28.0855 

Double_t fitf(Double_t *x, Double_t *par)
// Fit Function for the ATLAS Data
{
	Double_t arg = par[3] + par[2]*TMath::Log(x[0]);
	Double_t fitval = par[0]*pow(1+x[0],par[1])/pow(x[0],arg);
	return fitval;
}

vector<TLorentzVector> gammaprod(Double_t pt,Double_t eta,Double_t phi){
// For a given pi-zero pT calculates two photon four vectors and outputs 
// them as a vector

	vector<TLorentzVector> outv;

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

	// TVector3 test1 = gfourv1.Vect();
	// TVector3 test2 = gfourv2.Vect();

	// Egamma1 = gamma*(PIZERO_MASS/2 - boostv*test1);
	// Egamma2 = gamma*(PIZERO_MASS/2 - boostv*test2);

	// Boost back into the lab frame
	gfourv1.Boost(-boostv);
	gfourv2.Boost(-boostv);

	// std::cout  << "TLorentz: " << gfourv1.E() << "\tBy Hand E1: " << Egamma1 << std::endl;

	outv.push_back(gfourv1);
	outv.push_back(gfourv2);

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

Double_t cprefactor(Double_t Z, Double_t mass){
// Calculate the classical radius of a particle

	Double_t fstruc = 1.0/137.035999084;
	Double_t spdlgt = 299792458;
	Double_t echge = 1.602e-19;
	Double_t k = 8.99e9;

	Double_t kgmass = mass*pow(10,9)*echge/(spdlgt*spdlgt);

	Double_t rc = k*echge*echge/(kgmass*spdlgt*spdlgt);

	Double_t prefactor = 4.0 * fstruc * ZEl * ZEl * rc * rc;

	return prefactor;

}

Double_t cexpr(Double_t *x, Double_t *par){
// Expression for the differetial cross section

	Double_t xx = x[0];
	Double_t prefactor = cprefactor(ZEl,MUON_MASS);
	Double_t nCross = prefactor*(1 - 4.0/3.0*xx*(1-xx))*TMath::Log(calcW(ZEl,AEl,par[0],xx));
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
	sigmainf = 7.0/9.0*cprefactor(ZEl,MUON_MASS)*TMath::Log(calcW(ZEl,AEl,Egamma,1));

	crssec = func->Integral(xmin,xmax)/sigmainf;

	if(gRandom->Uniform(0,1)< crssec){
		return true;
	}else{
		return false;
	}

}

Double_t calcC1(Double_t xplus, Double_t Egamma, Double_t A){

	Double_t C1 = pow(0.35*pow(A,0.27),2)/(xplus*(1-xplus)*Egamma/MUON_MASS);

	return C1;

}

Double_t calcC2(Double_t xplus, Double_t Egamma, Double_t t, Double_t Z){

	Double_t xminus = 1-xplus;
	Double_t Z13 = pow(Z,-1.0/3.0);

	Double_t leadcoff = 4/sqrt(xplus*xminus);
	Double_t term1 = pow(MUON_MASS/(2*Egamma*xplus*xminus*t),2);
	Double_t term2 = pow(e_MASS/(183*Z13*MUON_MASS),2);

	Double_t C2 = leadcoff*pow(term1 + term2, 2);

	return C2;

}

Double_t calcrhomax2(Double_t t, Double_t A){

	return 1.9/pow(A,0.27)*(1/t - 1);
}

Double_t calcbeta(Double_t C2, Double_t rhomax2){

	return TMath::Log((C2 + pow(rhomax2,2))/C2);
}

Double_t f1(Double_t *t, Double_t *par){

	Double_t tt = t[0];
	Double_t Egamma = par[0];
	Double_t xplus = par[1];

	Double_t C1 = calcC1(xplus,Egamma,ZEl);

	Double_t num = 1 - 2*xplus*(1-xplus) + 4*xplus*(1-xplus)*tt*(1 - tt);
	Double_t denom = 1 + C1/(tt*tt);

	return num/denom;

}

Double_t f2(Double_t *psi, Double_t *par){

	Double_t Psi = psi[0];
	Double_t xplus = par[0];
	Double_t xminus = 1-par[0];
	Double_t t = par[1];

	Double_t val;

	val = 1 - 2*xplus*xminus + 4*xplus*xminus*t*(1 - t)*(1 + TMath::Cos(2*Psi));

	return val;
}

vector<TLorentzVector> muonprod(TLorentzVector photo4vect){
// Generates 4 vectors for muon +- given a photon four vector

	Double_t Egamma = photo4vect.E();

	Double_t xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
	Double_t xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

	Double_t xplus = xmin + gRandom->Uniform(0,1)*(xmax-xmin);

	Double_t Eplus = xplus*Egamma;
	Double_t Eminus = Egamma - Eplus;

	Double_t pplus = sqrt(pow(Eplus,2) - pow(MUON_MASS,2));
	Double_t pminus = sqrt(pow(Eminus,2) - pow(MUON_MASS,2));

	// Variables needed for angular generation
	Double_t t, u;
	Double_t Psi;
	Double_t thetam, thetap, phi;
	Double_t rho, rhomax2, beta;
	Double_t C2;
	Double_t gammap, gammam;
	Double_t R1, R2, R3;

	TVector3 nplus, nminus;

	TF1 *func1 = new TF1("f1",f1,0,1,2);
	TF1 *func2 = new TF1("f2",f2,0,2*TMath::Pi(),2);

	func1->SetParameters(Egamma,xplus);

	bool valtest = true;

	// Run until satisfactory angles are generated
	while(valtest == true){

		t = gRandom->Uniform(0,1);
		R1 = gRandom->Uniform(0,1);
		R2 = gRandom->Uniform(0,1);
		R3 = gRandom->Uniform(0,1);

		if(func1->Eval(t) < R1) continue;

		func2->SetParameters(xplus,t);
		Psi = gRandom->Uniform(0,2*TMath::Pi());

		if(func2->Eval(Psi) < R2) continue;

		gammap = xplus*Egamma/MUON_MASS;
		gammam = (1-xplus)*Egamma/MUON_MASS;
		u = sqrt(1/t - 1);

		thetap = (u + rho/2*TMath::Cos(Psi))/gammap;
		thetam = (u - rho/2*TMath::Cos(Psi))/gammam;
		phi = rho/u*TMath::Sin(Psi);

		valtest = false;

	}

	Double_t phi0 = gRandom->Uniform(0,2*TMath::Pi());

	nplus.SetX(TMath::Sin(thetap)*TMath::Cos(phi0 + phi/2));
	nplus.SetY(TMath::Sin(thetap)*TMath::Sin(phi0 + phi/2));
	nplus.SetZ(TMath::Cos(phi0 + phi/2));

	nminus.SetX(-TMath::Sin(thetap)*TMath::Cos(phi0 - phi/2));
	nminus.SetY(-TMath::Sin(thetap)*TMath::Sin(phi0 - phi/2));
	nminus.SetZ(TMath::Cos(phi0 + phi/2));

	TVector3 zhat(0.0,0.0,1.0);
	TVector3 vec = photo4vect.Vect();
	vec = vec.Unit();

	TVector3 crossprod = vec.Cross(zhat);
	Double_t theta = vec.Theta();

	nplus.Rotate(-theta,crossprod.Unit());
	nminus.Rotate(-theta,crossprod.Unit());

	TLorentzVector muplus4vec(pplus*nplus,Eplus);
	TLorentzVector muminus4vec(pminus*nminus,Eminus);

	vector<TLorentzVector> outv;

	outv.push_back(muplus4vec);
	outv.push_back(muminus4vec);

	return outv;

}

//--------------------------------------------------------------------
// main()
//--------------------------------------------------------------------

void muonconvmc(){

//--------------------------------------------------------------------
// Graphics
//--------------------------------------------------------------------

Int_t nbins = 100;

// TCanvas *Canvas1 = new TCanvas("Canvas1","Photon Energy Dist.",0,100,600,500);
TCanvas *Canvas2 = new TCanvas("Canvas2","",0,100,600,500);
TCanvas *Canvas3 = new TCanvas("Canvas3","",0,100,600,500);
TCanvas *Canvas4 = new TCanvas("Canvas4","",0,100,600,500);
//TCanvas *Canvas5 = new TCanvas("Canvas5","X+ Spectra",0,100,600,500);

TH1F *PhotonEHist = new TH1F("PhotonEHist","Photon Pt Dist.; Photon Energy [GeV]; Count [#]",nbins,6,25);

TH1F *AtlasH = new TH1F("AtlasH","MC Atlas Histogram ; pT [GeV]; Count [#]",nbins,0.5,50);

TH1F *MuPlusHist = new TH1F("MuPlusHist","Positive Muon pT Spectra; pT [GeV]; Count [#]",nbins,6,30);

TH1F *MuMinusHist = new TH1F("MuMinusHist","Negative Muon pT Spectra; pT [GeV]; Count [#]",nbins,6,30);

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

TVector3 photo3vect;

vector<TLorentzVector> photo4vect;
Double_t photoarr[2];
vector<TLorentzVector> outv;

Double_t nDiffCross;

for(UInt_t i = 0; i < NDECAY; ++i){

	// Get a pT for the neutral pion
	pizero_pt = func->GetRandom(6,100);

	AtlasH->Fill(pizero_pt);

	// Generate random eta & phi (these are assumed flat for now)
	eta = gRandom->Uniform(-1,1);
	phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());

	// Generate two photons provided the above parameters
	photo4vect = gammaprod(pizero_pt,eta,phi);

	photoarr[0] = photo4vect.at(0).E();
	photoarr[1] = photo4vect.at(1).E();

	// std::cout << "Photon 1 Energy: " << photoarr[0] << 
	// "\tPhoton 2 Energy: " << photoarr[1] << std::endl;

	for(UInt_t j = 0; j < 2; ++j){
		if(photoarr[j] > 6){
			PhotonEHist->Fill(photo4vect.at(j).Pt());
			if(convprob(photoarr[j]) == true){
				outv = muonprod(photo4vect.at(j));
				if(outv.at(0).Pt() > 6 && outv.at(1).Pt() > 6){
					MuPlusHist->Fill(outv.at(0).Pt());
					MuMinusHist->Fill(outv.at(1).Pt());
				}
			}
		}
	}
}

f->Close();

// Canvas1->cd();
// AtlasH->Draw();

Canvas2->cd();
PhotonEHist->Draw();

Canvas3->cd();
MuPlusHist->Draw();

Canvas4->cd();
MuMinusHist->Draw();

}