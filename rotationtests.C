
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

#define NDECAY 5e6 // Number of Decays computed

#define ZEl 4
#define AEl 9.01218

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

void rotationtests(){

	Double_t xplus = 0.5;
	Double_t Egamma = 10;

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

	// nplus.SetX(TMath::Sin(thetap)*TMath::Cos(phi0 + phi/2));
	// nplus.SetY(TMath::Sin(thetap)*TMath::Sin(phi0 + phi/2));
	// nplus.SetZ(TMath::Cos(phi0 + phi/2));

	// nminus.SetX(-TMath::Sin(thetap)*TMath::Cos(phi0 - phi/2));
	// nminus.SetY(-TMath::Sin(thetap)*TMath::Sin(phi0 - phi/2));
	// nminus.SetZ(TMath::Cos(phi0 + phi/2));

	nplus.SetX(-1.0);
	nplus.SetY(0.0);
	nplus.SetZ(0.0);

	nminus.SetX(0.0);
	nminus.SetY(1.0);
	nminus.SetZ(0.0);

	TVector3 zhat(0.0,0.0,1.0);

	TVector3 vec(1.0,1.0,0.0);	// Photon vector
	vec = vec.Unit();

	TVector3 crossprod = vec.Cross(zhat);
	Double_t theta = vec.Theta();

	vec.Rotate(-theta,crossprod.Unit());

	nplus.Rotate(-theta,crossprod.Unit());
	nminus.Rotate(-theta,crossprod.Unit());

	std::cout << "nplus: " << "(" << nplus.X() << ", "  << nplus.Y() << ", " << nplus.Z() << ")" << std::endl;
	
	std::cout << "nminus: " << "(" << nminus.X() << ", "  << nminus.Y() << ", " << nminus.Z() << ")" << std::endl;


	// std::cout << "Vec Theta: " << vec.Theta() << std::endl;

	// std::cout << "Theta: " << vec.Theta() << "\tMag: " << vec.Mag() << std::endl;

}





