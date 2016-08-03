#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>

#ifndef __CINT__

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"

#include "RooFitResult.h"
#include "RooMyPdf.h"

#endif

// functions needed :
// Jacobi polynomials
// factorial function
// combination function


using namespace std;
using namespace RooFit;

check()
{


//RooRealVar x("x","x",-10,10) ;
//RooRealVar alpha("alpha","alpha",1.0,0.,10.) ;


//RooRealVar rooB0_mass("B0_mass","B0_mass",5.28,5.27,5.3); //5.28,5.22,5.34
RooRealVar rooKPi_mass("KPi_mass","KPi_mass",1.0,0.6,2.2);//0.89,0.5,3.0
//RooRealVar rooJpsi_mass("Jpsi_mass","Jpsi_mass",3.09,2.97,3.21); //3.0,2.0,4.0
//RooRealVar rooB0_3mom("B0_3mom","B0_3mom",10.0,10.0,10.0); //24.0,10.0,30.0
//cout <<rooB0_3mom->getVal();
RooRealVar rooTheta_Kstar("Theta_Kstar","Theta_Kstar",0.5*TMath::Pi(),0.0,TMath::Pi()); //0.0,3.15
RooRealVar rooTheta_Jpsi("Theta_Jpsi","Theta_Jpsi",0.5*TMath::Pi(),0.0,TMath::Pi()); //0.0,3.15
RooRealVar rooPhi("Phi","Phi",0.25,-TMath::Pi(),TMath::Pi()); //-3.15,3.15

/*

complex<double> hK0_800(1.12*TMath::Cos(2.3),1.12*TMath::Sin(2.3));
complex<double> hK892(1*TMath::Cos(0),1*TMath::Sin(0));
complex<double> hK1410(0.12*TMath::Cos(0.81),0.12*TMath::Sin(0.81));
complex<double> hK0_1430(0.89*TMath::Cos(-2.17),0.89*TMath::Sin(-2.17));
complex<double> hK2_1430(4.66*TMath::Cos(-0.32),4.66*TMath::Sin(-0.32));
complex<double> hK1680(0.14*TMath::Cos(-2.46),0.14*TMath::Sin(-2.46));
complex<double> hK3_1780(16.8*TMath::Cos(-1.43),16.8*TMath::Sin(-1.43));
complex<double> hK0_1950(0.24*TMath::Cos(-2.39),0.24*TMath::Sin(-2.39));
complex<double> hK2_1980(4.53*TMath::Cos(-0.26),4.53*TMath::Sin(-0.26));
complex<double> hK4_2045(590*TMath::Cos(-2.66),590*TMath::Sin(-2.66));

RooRealVar rooAmp_0_K0_800("rooAmp_0_K0_800","rooAmp_0_K0_800",0.0,-10.0,10.0);
RooRealVar rooPhase_0_K0_800("rooPhase_0_K0_800","rooPhase_0_K0_800",0.0,-3.15,3.15);

RooRealVar rooAmp_0_K892("rooAmp_0_K892","rooAmp_0_K0_800",0.0,-10.0,10.0); //edit
RooRealVar rooPhase_0_K892("rooPhase_0_K892","rooPhase_0_K0_800",0.0,-3.15,3.15); //edit
*/


//RooGenericPdf g("g","sqrt(abs(alpha*x))+0.1",RooArgSet(x,alpha)) ;
//RooMyPdf g("g","compiled class g",x,alpha) ;
//RooMyPdf g("g","compiled class g",rooB0_mass,rooKPi_mass,rooJpsi_mass,rooB0_3mom,rooTheta_Kstar,rooPhi,rooTheta_Jpsi) ;
RooMyPdf g("g","compiled class g",rooKPi_mass,rooTheta_Kstar,rooPhi,rooTheta_Jpsi) ;

//cout << "PDF val :" << g.getVal() << endl;
//return;

//RooDataSet* data = g.generate(RooArgSet(rooKPi_mass,rooTheta_Kstar,rooPhi,rooTheta_Jpsi),2000) ;
RooPlot* frame = rooKPi_mass.frame() ;
//TFile* f_data = new TFile ("pdf_gen.root");

//data->plotOn(frame);
//gData->plotOn(frame);
//g.plotOn(frame) ;
//alpha=1e-4 ;
g.plotOn(frame,LineColor(kRed)) ;
frame->Draw() ;

//TFile f("pdf_gen.root","RECREATE") ;
//data->Write();
//f.Close();


}