//==========================================================================
// This file has been automatically generated for C++
// MadGraph 5 v. 2.0.0.beta4, 2013-06-22
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_HiggsCharac_UFO_H
#define Parameters_HiggsCharac_UFO_H

#include <complex> 

#include "read_slha.h"
using namespace std; 

class Parameters_HiggsCharac_UFO
{
  public:

    static Parameters_HiggsCharac_UFO * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double WX2, WX1, WX0, WW, WZ, WT, ymtau, ymt, ymb, aS, Gf, aEWM1, MX2, MX1,
        MX0, MZ, MTA, MT, MB, kw, kz, ka, kg, kl, kq, kz5, kz3, kz1, kw5, kw4,
        kw3, kw2, kw1, kfb, kfa, kHdw, kHdz, kHda, kAww, kHww, kAzz, kHzz,
        kAgg, kHgg, kAza, kHza, kAaa, kHaa, kAll, kHll, kAbb, kHbb, kAtt, kHtt,
        kSM, ca, Lambda, cabi, cos__cabi, sin__cabi, ca__exp__2, sa, g1, gw,
        MZ__exp__2, MZ__exp__4, sqrt__2, nb__2__exp__0_75, MX0__exp__2, aEW,
        MW, sqrt__aEW, ee, MW__exp__2, sw2, cw, sqrt__sw2, sw, ad, al, an, au,
        bd, bl, bn, bu, gwwz, vev, ee__exp__2, gAaa, cw__exp__2, gAza, gHaa,
        gHza, vev__exp__2, lam, yb, yt, ytau, muH, sw__exp__2;
    std::complex<double> CKM11, CKM12, CKM21, CKM22, complexi, conjg__CKM11,
        conjg__CKM12, conjg__CKM21, conjg__CKM22;
    // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2, gAgg, gHgg; 
    // Model couplings independent of aS
    std::complex<double> GC_61, GC_89, GC_35, GC_37, GC_56; 
    // Model couplings dependent on aS
    std::complex<double> GC_53, GC_18; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_HiggsCharac_UFO * instance; 
}; 

#endif  // Parameters_HiggsCharac_UFO_H

