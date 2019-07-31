// crssecint.C
// Marcos Flores
// 2019 July 29;
//
// Numerically integrate to calculate the cross section for \gamma  
// to \mu^+
//
// -------------------------------------------------------------------

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
#define ASi 28.085

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

Double_t cprefactor(Double_t Z, Double_t mass){

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

	Double_t xx = x[0];
	Double_t prefactor = cprefactor(ZSi,MUON_MASS);
	Double_t nCross = prefactor*(1 - 4.0/3.0*xx*(1-xx))*TMath::Log(calcW(ZSi,ASi,par[0],xx));
	return nCross;
}


void crssecint(){

	Double_t xmin;
	Double_t xmax;
	Double_t Egamma;

	Double_t sigmainf;

	UInt_t num = 8;

	Double_t engs[num]; 
	Double_t crscs[num];

	for(UInt_t i = 0; i < num; ++i){
		Egamma = pow(10,i);
		engs[i] = Egamma;

		xmin = 1.0/2.0 - sqrt(1.0/4.0 - MUON_MASS/Egamma);
		xmax = 1.0/2.0 + sqrt(1.0/4.0 - MUON_MASS/Egamma);

		TF1 *func = new TF1("cexpr",cexpr,xmin,xmax,1);
		func->SetParameter(0,Egamma);

		sigmainf = 7.0/9.0*cprefactor(ZSi,MUON_MASS)*TMath::Log(calcW(ZSi,ASi,Egamma,1));

		crscs[i] = func->Integral(xmin,xmax)/sigmainf;

	}


TGraph *g = new TGraph(num,engs,crscs);

g->Draw("AC*");

}