#include <string>
#include <Riostream.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>


#define MUON_MASS 0.10565837 // GeV
#define e_MASS 5.109989e-4 // GeV

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
	Double_t delta = pow(MUON_MASS,2)/(2*Egamma*xplus*(1-xplus));
	Double_t W = Winf*(1 + (Dn*sqrt(TMath::E()) - 2)*delta/MUON_MASS)/(1 + B*pow(Z,-1.0/3.0)*sqrt(TMath::E())*delta/e_MASS);

	return W;

}

void pairprodmc(){

UInt_t nbins = 50;

TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",0,100,600,500);

THStack *hs = new THStack("hs","");
TH1F *MuonSpHist_10 = new TH1F("MuonSpHist_10","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);
TH1F *MuonSpHist_100 = new TH1F("MuonSpHist_100","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);
TH1F *MuonSpHist_1000 = new TH1F("MuonSpHist_1000","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);

TH1F *Tests = new TH1F("Tests","Tests",nbins,0,1);

vector<TH1F* > test;

test.push_back(MuonSpHist_10);	
test.push_back(MuonSpHist_100);	
test.push_back(MuonSpHist_1000);	

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

for(UInt_t i = 0; i < 3; ++i){

	Egamma = pow(10.0,i + 1);
	Wmax = calcW(Z,A,Egamma,0.5); // See definition in paper

	// Calculate xmin/xmax

	xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
	xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

	UInt_t j = 0;

	while(j < nev){
		
		for(UInt_t k = 0; k < 2; ++k){
		R[k] = gRandom->Uniform(0,1);
		}

		xplus = xmin + R[0]*(xmax-xmin);
		W = calcW(Z,A,Egamma,xplus);

		param = (1-(4.0/3.0)*xplus*(1-xplus))*TMath::Log(W)/TMath::Log(Wmax);

		if(param > R[1]){
			test.at(i)->Fill(xplus);
			++j;
		}

	}

}

test.at(0)->SetLineStyle(9);
test.at(1)->SetLineStyle(2);

Canvas1->cd();
test.at(0)->Draw();
test.at(1)->Draw("SAME");
test.at(2)->Draw("SAME");

}
