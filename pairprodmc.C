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

UInt_t nbins = 100;

TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",0,100,600,500);

THStack *hs = new THStack("hs","");
TH1F *MuonSpHist_10 = new TH1F("MuonSpHist_10","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);
TH1F *MuonSpHist_100 = new TH1F("MuonSpHist_100","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);
TH1F *MuonSpHist_1000 = new TH1F("MuonSpHist_1000","Muon + Energy Dist.; x+; Count [#]",nbins,0,1);

vector<TH1F* > test;

std::cout << "Size of test: " << test.size() << "\n";

test.push_back(MuonSpHist_10);	
test.push_back(MuonSpHist_100);	
test.push_back(MuonSpHist_1000);	

std::cout << "Size of test: " << test.size() << "\n";


UInt_t nev = 1e6;

UInt_t Z = 4;
Double_t A = 9.01218;


Double_t xplus;
Double_t xmin;
Double_t xmax;
Double_t r[6];

Double_t count = 0;

for(UInt_t i = 0; i < 3; ++i){

	Double_t Egamma = pow(10,i + 1);

	// Calculate xmin/xmax

	xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
	xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

	for(UInt_t j=0; j < nev; ++j){
		
		for(UInt_t k = 0; k < 6; ++k){
		r[k] = gRandom->Uniform();
		}

		xplus = xmin + r[0]*(xmax-xmin);
		Double_t W = calcW(Z,A,Egamma,xplus);
		Double_t Wmax = calcW(Z,A,Egamma,0.5); // See definition in paper

		Double_t param = (1-4.0/3.0*xplus*(1-xplus))*TMath::Log(W)/TMath::Log(Wmax);

		if(param > r[1]){
			test.at(i)->Fill(xplus);
			++count;
		}

	}

}

Canvas1->cd();
test.at(0)->Draw();
test.at(1)->Draw("SAME");
test.at(2)->Draw("SAME");

}
