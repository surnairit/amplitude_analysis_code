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
#include "Dalitz_contour.h"

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
    RooRealVar rooKPi_mass("KPi_mass","KPi_mass",TMath::Sqrt(0.7),0.6,2.2);//1.0,0.6,2.2
    //RooRealVar rooJpsi_mass("Jpsi_mass","Jpsi_mass",3.09,2.97,3.21); //3.0,2.0,4.0
    RooRealVar rooJpsiPi_mass("JpsiPi_mass","JpsiPi_mass",TMath::Sqrt(23),3.2,4.9); //3.0,2.0,4.0
    //RooRealVar rooB0_3mom("B0_3mom","B0_3mom",10.0,10.0,10.0); //24.0,10.0,30.0
    //cout <<rooB0_3mom->getVal();
    //RooRealVar rooTheta_Kstar("Theta_Kstar","Theta_Kstar",0.5*TMath::Pi(),0.0,TMath::Pi()); //0.0,3.15
    RooRealVar rooTheta_Jpsi("Theta_Jpsi","Theta_Jpsi",0.5*TMath::Pi(),0.0,TMath::Pi()); //0.0,3.15
    RooRealVar rooPhi("Phi","Phi",0.25,-TMath::Pi(),TMath::Pi()); //-3.15,3.15
    
    
    RooFormulaVar rooKPi_mass2For("rooKPi_mass2For","m^{2}(K^{-}#pi^{+}) [GeV^{2}]","pow(rooKPi_mass,2)",rooKPi_mass);
    RooRealVar rooKPi_mass2("rooKPi_mass2","m^{2}(K^{-}#pi^{+}) [GeV^{2}]",TMath::Power(rooKPi_mass.getVal(),2),TMath::Power(rooKPi_mass.getMin(),2),TMath::Power(rooKPi_mass.getMax(),2));
    
    RooFormulaVar rooJpsiPi_mass2For("rooJpsiPi_mass2For","m^{2}(J/#psi #pi^{+}) [GeV^{2}]","pow(rooJpsiPi_mass,2)",rooJpsiPi_mass);
    RooRealVar rooJpsiPi_mass2("rooJpsiPi_mass2","m^{2}(J/#psi #pi^{+}) [GeV^{2}]",TMath::Power(rooJpsiPi_mass.getVal(),2),TMath::Power(rooJpsiPi_mass.getMin(),2),TMath::Power(rooJpsiPi_mass.getMax(),2));
    
    
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
    //RooMyPdf g("g","compiled class g",rooKPi_mass,rooTheta_Kstar,rooPhi,rooTheta_Jpsi) ;
    RooMyPdf g("g","compiled class g",rooKPi_mass,rooJpsiPi_mass,rooPhi,rooTheta_Jpsi) ;
    
    //cout << "PDF val :" << g.getVal() << endl;
    //return;
    
    //RooDataSet* data = g.generate(RooArgSet(rooKPi_mass,rooTheta_Kstar,rooPhi,rooTheta_Jpsi),2000) ;
    RooDataSet* data = g.generate(RooArgSet(rooKPi_mass,rooJpsiPi_mass,rooPhi,rooTheta_Jpsi),1000) ; //60000
    
    RooRealVar* var;
    RooRealVar* var2;
    RooArgSet* set;
    
    TH2F* hDalitz = new TH2F("hDalitz","Generated Dalitz Plot;m^{2}_{(K#pi)};m^{2}_{(J/#psi#pi)}",600,0,6,1700,9,26);
    
    for (Int_t iEvent=0; iEvent<data->numEntries(); ++iEvent) {
        //    kinematicVars = *(data->get(iEvent)) ; // this get method will propagate the RooRealVars values of the event to all the corresponding RooRealVars
        set = data->get(iEvent);
        var = (RooRealVar*)set->find("KPi_mass");
        var2 = (RooRealVar*)set->find("JpsiPi_mass");
        double mKPi_temp = var->getVal();
        double mJpsiPi_temp = var2->getVal();
        hDalitz->Fill(mKPi_temp*mKPi_temp,mJpsiPi_temp*mJpsiPi_temp);
        //cout << var->getVal() << endl;
        
    }
    
    Int_t m = 19800;
    Double_t x[19800], m23_max[19800], m23_min[19800];
    Double_t E2[19800], E3[19800];
    
    Double_t m_mother = mB;
    Double_t m_dau1 = mK, m_dau2 = mPi, m_dau3 = mJpsi ;
    
    Double_t m12_min = (m_dau1+m_dau2)*(m_dau1+m_dau2);
    Double_t m12_max = (m_mother-m_dau3)*(m_mother-m_dau3);
    Double_t step = (m12_max - m12_min)/(m-1);
    
    x[0] = m12_min + 0.0001;
    
    for (Int_t k=1; k<m; k++ )
        x[k] = x[k-1] + step;
    Int_t n = 19799;
    for (Int_t i=0; i<n; i++) {
        E2[i] = (x[i] - m_dau1*m_dau1 + m_dau2*m_dau2)/(2*sqrt(x[i]));
        E3[i] = (m_mother*m_mother - x[i] - m_dau3*m_dau3)/(2*sqrt(x[i]));
        m23_min[i] = (E2[i]+E3[i])*(E2[i]+E3[i]) - TMath::Power((sqrt(E2[i]*E2[i] - m_dau2*m_dau2) + sqrt(E3[i]*E3[i] - m_dau3*m_dau3)),2);
        m23_max[i] = (E2[i]+E3[i])*(E2[i]+E3[i]) - TMath::Power((sqrt(E2[i]*E2[i] - m_dau2*m_dau2) - sqrt(E3[i]*E3[i] - m_dau3*m_dau3)),2);
    }
    
    TGraph *cont_up = new TGraph(n,x,m23_min); cont_up->SetLineWidth(3);
    TGraph *cont_down = new TGraph(n,x,m23_max); cont_down->SetLineWidth(3);
    
    TCanvas *cDal = new TCanvas("cDal","cDal",800,600);
    
    cDal->cd();
    hDalitz->Draw("colz");
    cont_down->Draw("lsame");
    cont_up->Draw("lsame");
    //cDal->SaveAs("gen_dalitz.pdf");
    cDal->SaveAs("gen_dalitz_Zmodel.pdf");
    
    
    RooPlot* frame = rooKPi_mass.frame() ;
    data->plotOn(frame);
    g.plotOn(frame,LineColor(kRed)) ;
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    frame->Draw() ;
    //c1->SaveAs("gen_kpi.pdf");
    c1->SaveAs("gen_kpi_Zmodel.pdf");
    
    RooPlot* frame2 = rooJpsiPi_mass.frame() ;
    data->plotOn(frame2);
    g.plotOn(frame2,LineColor(kBlue)) ;
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    c2->cd();
    frame2->Draw() ;
    c2->SaveAs("gen_jpsipi.pdf");
    c2->SaveAs("gen_jpsipi_Zmodel.pdf");
    
    
    
    
    //TFile* f_data = new TFile ("pdf_gen.root");
    
    //data->plotOn(frame);
    //gData->plotOn(frame);
    //g.plotOn(frame) ;
    //alpha=1e-4 ;
    //g.plotOn(frame,LineColor(kRed)) ;
    //frame->Draw() ;
    
    
    
    //TFile f("pdf_gen.root","RECREATE") ;
    //data->Write();
    //f.Close();
    
    
}
