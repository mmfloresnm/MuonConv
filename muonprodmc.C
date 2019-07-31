#include <string>
#include <Riostream.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>

#define MUON_MASS 0.10565837 // GeV
#define e_MASS 5.109989e-4 // GeV

#define ZSi 14.0
#define ASi 28.0855 

Double_t calcW(UInt_t  Z, Double_t A, Double_t Egamma, Double_t xplus){

	Double_t B;
	Double_t Dn;
	Double_t third = 1.0/3.0;

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

	if(Eplus > 0.9*Egamma){

		//std::cout << "Egamma: " << Egamma << "\txplus: " << xplus << "\t nDiffCross: " << nDiffCross << "\n";
	}

	vector<Double_t> outv;

	outv.push_back(nDiffCross);
	outv.push_back(Eplus);
	outv.push_back(Eminus);

	return outv;

}

void muonprodmc(){

UInt_t nbins = 50;

TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",0,100,1200,500);
TCanvas *Canvas2 = new TCanvas("Canvas2","Canvas2",0,100,1200,500);
TCanvas *Canvas3 = new TCanvas("Canvas3","Canvas3",0,100,1200,500);

THStack *hs = new THStack("hs","");
TH1F *CrossSecHist_6 = new TH1F("CrossSecHist_6","Cross Section - 6 GeV Photons; Cross Section; Count [#]",nbins,0,1);
TH1F *CrossSecHist_10 = new TH1F("CrossSecHist_10","Cross Section - 10 GeV Photons; Cross Section; Count [#]",nbins,0,1);
TH1F *CrossSecHist_15 = new TH1F("CrossSecHist_15","Cross Section - 15 GeV Photons; Cross Section; Count [#]",nbins,0,1);

TH1F *MuPlus_6 = new TH1F("MuPlus_6","Muon (+) Spectra - 6 GeV Photons; Energy [GeV]; Count [#]",nbins,0,20);
TH1F *MuPlus_10 = new TH1F("MuPlus_10","Muon (+) Spectra - 10 GeV Photons; Energy [GeV]; Count [#]",nbins,0,25);
TH1F *MuPlus_15 = new TH1F("MuPlus_15","Muon (+) Spectra - 15 GeV Photons; Energy [GeV]; Count [#]",nbins,0,25);

TH1F *Tests = new TH1F("Tests","Tests",nbins,0,1);

vector<TH1F* > CSHists;
vector<TH1F* > ESHists;

CSHists.push_back(CrossSecHist_6);	
CSHists.push_back(CrossSecHist_10);	
CSHists.push_back(CrossSecHist_15);	

ESHists.push_back(MuPlus_6);	
ESHists.push_back(MuPlus_10);	
ESHists.push_back(MuPlus_15);

UInt_t nev = 1e6;

UInt_t Z = 4;
Double_t A = 9.01218;

Double_t Egamma; 
Double_t xplus;
Double_t xmin;
Double_t xmax;
Double_t R[2];

Double_t W;
Double_t Wmax;
Double_t param;

Double_t count = 0;

vector<Double_t> outv;

Double_t earr[3] = {6, 10, 20};

// Calculation of actual cross section

for(UInt_t i = 0; i < 3; ++i){

	Egamma = earr[i];

	for(UInt_t j = 0; j < nev; ++j){

		outv = muonprod(Egamma);

		CSHists.at(i)->Fill(outv.at(0));
		ESHists.at(i)->Fill(outv.at(1));
	}

}

Canvas1->cd();
Canvas1->Divide(2, 1);
Canvas1->cd(1);
CSHists.at(0)->Draw();
Canvas1->cd(2);
ESHists.at(0)->Draw();

Canvas2->cd();
Canvas2->Divide(2, 1);
Canvas2->cd(1);
CSHists.at(1)->Draw();
Canvas2->cd(2);
ESHists.at(1)->Draw();

Canvas3->cd();
Canvas3->Divide(2, 1);
Canvas3->cd(1);
CSHists.at(2)->Draw();
Canvas3->cd(2);
ESHists.at(2)->Draw();
// Canvas2->cd();
// CSHists.at(1)->Draw();
// Canvas3->cd();
// CSHists.at(2)->Draw();

}
