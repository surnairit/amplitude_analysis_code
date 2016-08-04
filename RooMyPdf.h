/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOMYPDF
#define ROOMYPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>


#include <math.h> 
#include "TMath.h" 

const double mB = 5.27961;        // [GeV] B0 mass
const double mK = 0.493677;     // [GeV] Kaon mass
const double mPi = 0.13957;     // [GeV] Pion mass
const double mJpsi = 3.096916;    // [GeV] Jpsi mass

const double rB = 3.00; //5.0;          // [GeV^-1] meson radial parameter for B
const double rR = 3.00; //1.5;          // [GeV^-1] meson radial parameter for intermediate resonances

const double mK0_800 = 0.682;  // [GeV^-1] K0(800) mass
const double wK0_800 = 0.547;   // [GeV^-1] K0(800) width

const double mK892 = 0.89581;   // [GeV^-1] K*(892) mass - neutral only -> see pdg
const double wK892 = 0.0474;    // [GeV^-1] K*(892) width - neutral only -> see pdg

const double mK1410 = 1.414;   // [GeV^-1] K*(1410) mass
const double wK1410 = 0.232;    // [GeV^-1] K*(1410) width

const double mK0_1430 = 1.425;  // [GeV^-1] K0(1430) mass //1.425
const double wK0_1430 = 0.27;   // [GeV^-1] K0(1430) width

const double mK2_1430 = 1.4324; // [GeV^-1] K2(1430) mass - neutral only -> see pdg
const double wK2_1430 = 0.109;  // [GeV^-1] K2(1430) width - neutral only -> see pdg

const double mK1680 = 1.717;   // [GeV^-1] K*(1680) mass
const double wK1680 = 0.322;    // [GeV^-1] K*(1680) width

const double mK3_1780 = 1.776;  // [GeV^-1] K3(1780) mass
const double wK3_1780 = 0.159;   // [GeV^-1] K3(1780) width

const double mK0_1950 = 1.945;  // [GeV^-1] K0(1950) mass
const double wK0_1950 = 0.201;   // [GeV^-1] K0(1950) width

const double mK2_1980 = 1.973; // [GeV^-1] K2(1980) mass - neutral only -> see pdg
const double wK2_1980 = 0.373;  // [GeV^-1] K2(1980) width - neutral only -> see pdg

const double mK4_2045 = 2.045;  // [GeV^-1] K4(2045) mass
const double wK4_2045 = 0.198;   // [GeV^-1] K4(2045) width


 
class RooMyPdf : public RooAbsPdf {
public:
  RooMyPdf() {} ; 
  RooMyPdf(const char *name, const char *title,
//	      RooAbsReal& _x,
//	      RooAbsReal& _alpha);
//        RooAbsReal& _rooB0_mass,
        RooAbsReal& _rooKPi_mass,
        RooAbsReal& _rooJpsiPi_mass,
//        RooAbsReal& _rooB0_3mom,
//        RooAbsReal& _rooTheta_Kstar,
        RooAbsReal& _rooPhi,
        RooAbsReal& _rooTheta_Jpsi
//        RooAbsReal& _rooPhi
        );

  RooMyPdf(const RooMyPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMyPdf(*this,newname); }
  inline virtual ~RooMyPdf() { }
  
//================ Decay Momentum ====================
double sq_calc(double x,double y, double z) const;
double dec2mm (double m0, double m1, double m2) const;

//================ Blatt-Weisskopf Form Factors ======
double bwff(int l, double q, double q0, double r) const;

//================ Breit-Wigner Amplitude ============
//complex<double> bwamp(double m0,double w0,double m,double m_d1,double m_d2,int l,double f,double q0,double q) const;
complex<double> bwamp(double m0,double w0,double m,int l,double f,double q0,double q) const;
 
//================ Jacobi Polynomial =================
double jacobi_Pn (int n, double a, double b, double x) const;

//================ factorial =========================
int Factorial(int x) const;

//================ combination =======================
int Combination(int n, int r) const;

//================ wigner d calculations =============
double wigner_d (int j, int m1, int m2, double theta ) const;

//================ phase space =======================
double PHSP(double mKPicalc) const;

//================ Signal Density Calculation ========
//double get_signal_density (double mBcalc, double mKPicalc, double mJpsicalc, double pB, double theta_k, double phi, double theta_jpsi ) const;
//double get_signal_density (double mKPicalc, double theta_k, double phi, double theta_jpsi ) const;
double get_signal_density (double mKPicalc, double mJpsiPicalc , double phi, double theta_jpsi ) const;
protected:

//  RooRealProxy x ;
//  RooRealProxy alpha ;
  
//        RooRealProxy rooB0_mass;
        RooRealProxy rooKPi_mass;
        RooRealProxy rooJpsiPi_mass;
//        RooRealProxy rooB0_3mom;
//        RooRealProxy rooTheta_Kstar;
        RooRealProxy rooPhi;
        RooRealProxy rooTheta_Jpsi;
    

  
  Double_t evaluate() const ;

private:

  ClassDef(RooMyPdf,1) // Your description goes here...
};
 
#endif
