// muonconvmc.C
// Marcos Flores
// 2019 July 9;
//
// Monte-Carlo for \pi_0\to 2\gamma
// 
//
// -------------------------------------------------------------------

#include <string>
#include <Riostream.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>

#define PIZERO_MASS 0.134977 // GeV
#define minE 0.5 // GeV
#define maxE 50// GeV
#define NEVENT 5e5

Double_t fitf(Double_t *x, Double_t *par)
// Fit Function for the ATLAS Data
{
	Double_t arg = par[3] + par[2]*TMath::Log(x[0]);
	Double_t fitval = par[0]*pow(1+x[0],par[1])/pow(x[0],arg);
	return fitval;
}

void muonconvmc(){

//--------------------------------------------------------------------

Int_t nbins = 250;

TCanvas *Canvas1 = new TCanvas("Canvas1","Photon Energy Dist.",0,100,600,500);
//TCanvas *Canvas2 = new TCanvas("Canvas2","Tests",0,100,600,500);


TH1F *PhotonEHist = new TH1F("PhotonEHist","Photon Energy Dist.; Photon Energy [GeV]; Count [#]",nbins,minE,maxE);

TH1F *TestsH = new TH1F("TestsH","Tests; pT [GeV]; Count [#]",nbins,minE,maxE);

//--------------------------------------------------------------------

// Import data from ALTAS
 
// Open ATLAS File
TFile *f = 0;
 
f = new TFile("~/Desktop/HEPData-ins1419652-v1-Table_12.root");
if(!f) throw "Can't open file";
 
TH1F *AHist =(TH1F*)f->Get("Table 12/Hist1D_y1");
 
 
// Fit the Histogram
TF1 *func = new TF1("fitf",fitf,0.5,50,4);
func->SetParameters(1,1,1,-1);
AHist->Fit(func,"QN");
 
Double_t *params = func->GetParameters();
Double_t *x;

//--------------------------------------------------------------------

Int_t num = NEVENT; // Number of Decays

Double_t ctheta; // Random cos(\theta)
Double_t Epi;	 // Random neutral pion energy
Double_t gamma;  // Gamma of pion
Double_t beta;	 // Beta of pion

Double_t Egamma1;// Forward photon energy [GeV]
Double_t Egamma2;// Backward photon energy [GeV]

// Creation of Photon Energy Distributions

for(UInt_t i = 0; i < num; ++i){

	ctheta = gRandom->Uniform(-1,1);
	Epi = gRandom->Uniform(minE,maxE);

	x = &Epi;

	gamma = Epi/PIZERO_MASS;
	beta = sqrt(1-1/pow(gamma,2));

	Egamma1 = gamma*PIZERO_MASS/2*(1 + beta*ctheta);
	Egamma2 = gamma*PIZERO_MASS/2*(1 - beta*ctheta);

	PhotonEHist->Fill(Egamma1);
	PhotonEHist->Fill(Egamma2);
}

Canvas1->cd();
PhotonEHist->Draw();

f->Close();

}