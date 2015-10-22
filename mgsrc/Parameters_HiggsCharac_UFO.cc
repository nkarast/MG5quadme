//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 2.0.0.beta4, 2013-06-22
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_HiggsCharac_UFO.h"

// Initialize static instance
Parameters_HiggsCharac_UFO * Parameters_HiggsCharac_UFO::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_HiggsCharac_UFO * Parameters_HiggsCharac_UFO::getInstance()
{
  if (instance == 0)
    instance = new Parameters_HiggsCharac_UFO(); 

  return instance; 
}

void Parameters_HiggsCharac_UFO::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  WX2 = slha.get_block_entry("decay", 5000002, 7.567860e-02); 
  WX1 = slha.get_block_entry("decay", 5000001, 7.567860e-02); 
  WX0 = slha.get_block_entry("decay", 5000000, 7.567860e-02); 
  WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  MX2 = slha.get_block_entry("mass", 5000002, 1.250000e+02); 
  MX1 = slha.get_block_entry("mass", 5000001, 1.250000e+02); 
  MX0 = slha.get_block_entry("mass", 5000000, 1.250000e+02); 
  MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  kw = slha.get_block_entry("frblock", 38, 1.000000e+00); 
  kz = slha.get_block_entry("frblock", 37, 1.000000e+00); 
  ka = slha.get_block_entry("frblock", 36, 1.000000e+00); 
  kg = slha.get_block_entry("frblock", 35, 1.000000e+00); 
  kl = slha.get_block_entry("frblock", 34, 1.000000e+00); 
  kq = slha.get_block_entry("frblock", 33, 1.000000e+00); 
  kz5 = slha.get_block_entry("frblock", 32, 0.000000e+00); 
  kz3 = slha.get_block_entry("frblock", 31, 1.000000e+00); 
  kz1 = slha.get_block_entry("frblock", 30, 0.000000e+00); 
  kw5 = slha.get_block_entry("frblock", 29, 0.000000e+00); 
  kw4 = slha.get_block_entry("frblock", 28, 0.000000e+00); 
  kw3 = slha.get_block_entry("frblock", 27, 0.000000e+00); 
  kw2 = slha.get_block_entry("frblock", 26, 1.000000e+00); 
  kw1 = slha.get_block_entry("frblock", 25, 1.000000e+00); 
  kfb = slha.get_block_entry("frblock", 24, 1.000000e+00); 
  kfa = slha.get_block_entry("frblock", 23, 1.000000e+00); 
  kHdw = slha.get_block_entry("frblock", 22, 0.000000e+00); 
  kHdz = slha.get_block_entry("frblock", 21, 0.000000e+00); 
  kHda = slha.get_block_entry("frblock", 20, 0.000000e+00); 
  kAww = slha.get_block_entry("frblock", 19, 0.000000e+00); 
  kHww = slha.get_block_entry("frblock", 18, 0.000000e+00); 
  kAzz = slha.get_block_entry("frblock", 17, 0.000000e+00); 
  kHzz = slha.get_block_entry("frblock", 16, 0.000000e+00); 
  kAgg = slha.get_block_entry("frblock", 15, 1.000000e+00); 
  kHgg = slha.get_block_entry("frblock", 14, 1.000000e+00); 
  kAza = slha.get_block_entry("frblock", 13, 1.000000e+00); 
  kHza = slha.get_block_entry("frblock", 12, 1.000000e+00); 
  kAaa = slha.get_block_entry("frblock", 11, 1.000000e+00); 
  kHaa = slha.get_block_entry("frblock", 10, 1.000000e+00); 
  kAll = slha.get_block_entry("frblock", 9, 1.000000e+00); 
  kHll = slha.get_block_entry("frblock", 8, 1.000000e+00); 
  kAbb = slha.get_block_entry("frblock", 7, 1.000000e+00); 
  kHbb = slha.get_block_entry("frblock", 6, 1.000000e+00); 
  kAtt = slha.get_block_entry("frblock", 5, 1.000000e+00); 
  kHtt = slha.get_block_entry("frblock", 4, 1.000000e+00); 
  kSM = slha.get_block_entry("frblock", 3, 1.000000e+00); 
  ca = slha.get_block_entry("frblock", 2, 1.000000e+00); 
  Lambda = slha.get_block_entry("frblock", 1, 1.000000e+03); 
  cabi = slha.get_block_entry("ckmblock", 1, 2.277360e-01); 
  cos__cabi = cos(cabi); 
  CKM11 = cos__cabi; 
  sin__cabi = sin(cabi); 
  CKM12 = sin__cabi; 
  CKM21 = -sin__cabi; 
  CKM22 = cos__cabi; 
  ca__exp__2 = pow(ca, 2.); 
  sa = sqrt(1. - ca__exp__2); 
  g1 = 1.; 
  gw = 1.; 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  nb__2__exp__0_75 = pow(2., 0.75); 
  MX0__exp__2 = pow(MX0, 2.); 
  complexi = std::complex<double> (0., 1.); 
  conjg__CKM11 = conj(CKM11); 
  conjg__CKM12 = conj(CKM12); 
  conjg__CKM21 = conj(CKM21); 
  conjg__CKM22 = conj(CKM22); 
  aEW = 1./aEWM1; 
  MW = sqrt(MZ__exp__2/2. + sqrt(MZ__exp__4/4. - (aEW * M_PI * MZ__exp__2)/(Gf
      * sqrt__2)));
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.); 
  sw2 = 1. - MW__exp__2/MZ__exp__2; 
  cw = sqrt(1. - sw2); 
  sqrt__sw2 = sqrt(sw2); 
  sw = sqrt__sw2; 
  ad = (ee * (-0.5 + (2. * sw2)/3.))/(2. * cw * sw); 
  al = (ee * (-0.5 + 2. * sw2))/(2. * cw * sw); 
  an = ee/(4. * cw * sw); 
  au = (ee * (0.5 - (4. * sw2)/3.))/(2. * cw * sw); 
  bd = -ee/(4. * cw * sw); 
  bl = -ee/(4. * cw * sw); 
  bn = ee/(4. * cw * sw); 
  bu = ee/(4. * cw * sw); 
  gwwz = -((cw * ee)/sw); 
  vev = (2. * MW * sw)/ee; 
  ee__exp__2 = pow(ee, 2.); 
  gAaa = -ee__exp__2/(3. * pow(M_PI, 2.) * vev); 
  cw__exp__2 = pow(cw, 2.); 
  gAza = -((-5. + 8. * cw__exp__2) * sqrt(ee__exp__2 * Gf * MZ__exp__2))/(6. *
      nb__2__exp__0_75 * pow(M_PI, 2.) * vev);
  gHaa = (47. * ee__exp__2)/(72. * pow(M_PI, 2.) * vev); 
  gHza = ((-13. + 94. * cw__exp__2) * sqrt(ee__exp__2 * Gf * MZ__exp__2))/(36.
      * nb__2__exp__0_75 * pow(M_PI, 2.) * vev);
  vev__exp__2 = pow(vev, 2.); 
  lam = MX0__exp__2/(2. * vev__exp__2); 
  yb = (ymb * sqrt__2)/vev; 
  yt = (ymt * sqrt__2)/vev; 
  ytau = (ymtau * sqrt__2)/vev; 
  muH = sqrt(lam * vev__exp__2); 
  sw__exp__2 = pow(sw, 2.); 
}
void Parameters_HiggsCharac_UFO::setIndependentCouplings()
{
  GC_61 = (ee * complexi)/(sw * sqrt__2); 
  GC_89 = (ca * ee__exp__2 * complexi * kSM * vev)/(2. * sw__exp__2); 
  GC_35 = -((ca * complexi * kHdw)/Lambda); 
  GC_37 = -((ca * complexi * kHww)/Lambda); 
  GC_56 = (complexi * kAww * sa)/Lambda; 
}
void Parameters_HiggsCharac_UFO::setDependentParameters()
{
  sqrt__aS = sqrt(aS); 
  G = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
  gAgg = -G__exp__2/(8. * pow(M_PI, 2.) * vev); 
  gHgg = -G__exp__2/(12. * pow(M_PI, 2.) * vev); 
}
void Parameters_HiggsCharac_UFO::setDependentCouplings()
{
  GC_53 = (complexi * gAgg * kAgg * sa)/8.; 
  GC_18 = -(ca * complexi * gHgg * kHgg); 
}

// Routines for printing out parameters
void Parameters_HiggsCharac_UFO::printIndependentParameters()
{
  cout <<  "HiggsCharac_UFO model parameters independent of event kinematics:"
      << endl;
  cout << setw(20) <<  "WX2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WX2 << endl;
  cout << setw(20) <<  "WX1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WX1 << endl;
  cout << setw(20) <<  "WX0 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WX0 << endl;
  cout << setw(20) <<  "WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WW << endl;
  cout << setw(20) <<  "WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WZ << endl;
  cout << setw(20) <<  "WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WT << endl;
  cout << setw(20) <<  "ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ymtau << endl;
  cout << setw(20) <<  "ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ymt << endl;
  cout << setw(20) <<  "ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ymb << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "MX2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MX2 << endl;
  cout << setw(20) <<  "MX1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MX1 << endl;
  cout << setw(20) <<  "MX0 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MX0 << endl;
  cout << setw(20) <<  "MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MZ << endl;
  cout << setw(20) <<  "MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MTA << endl;
  cout << setw(20) <<  "MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MT << endl;
  cout << setw(20) <<  "MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MB << endl;
  cout << setw(20) <<  "kw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw << endl;
  cout << setw(20) <<  "kz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kz << endl;
  cout << setw(20) <<  "ka " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ka << endl;
  cout << setw(20) <<  "kg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kg << endl;
  cout << setw(20) <<  "kl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kl << endl;
  cout << setw(20) <<  "kq " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kq << endl;
  cout << setw(20) <<  "kz5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kz5 << endl;
  cout << setw(20) <<  "kz3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kz3 << endl;
  cout << setw(20) <<  "kz1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kz1 << endl;
  cout << setw(20) <<  "kw5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw5 << endl;
  cout << setw(20) <<  "kw4 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw4 << endl;
  cout << setw(20) <<  "kw3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw3 << endl;
  cout << setw(20) <<  "kw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw2 << endl;
  cout << setw(20) <<  "kw1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kw1 << endl;
  cout << setw(20) <<  "kfb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kfb << endl;
  cout << setw(20) <<  "kfa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kfa << endl;
  cout << setw(20) <<  "kHdw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHdw << endl;
  cout << setw(20) <<  "kHdz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHdz << endl;
  cout << setw(20) <<  "kHda " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHda << endl;
  cout << setw(20) <<  "kAww " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAww << endl;
  cout << setw(20) <<  "kHww " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHww << endl;
  cout << setw(20) <<  "kAzz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAzz << endl;
  cout << setw(20) <<  "kHzz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHzz << endl;
  cout << setw(20) <<  "kAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAgg << endl;
  cout << setw(20) <<  "kHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHgg << endl;
  cout << setw(20) <<  "kAza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAza << endl;
  cout << setw(20) <<  "kHza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHza << endl;
  cout << setw(20) <<  "kAaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAaa << endl;
  cout << setw(20) <<  "kHaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHaa << endl;
  cout << setw(20) <<  "kAll " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAll << endl;
  cout << setw(20) <<  "kHll " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHll << endl;
  cout << setw(20) <<  "kAbb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAbb << endl;
  cout << setw(20) <<  "kHbb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHbb << endl;
  cout << setw(20) <<  "kAtt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kAtt << endl;
  cout << setw(20) <<  "kHtt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kHtt << endl;
  cout << setw(20) <<  "kSM " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << kSM << endl;
  cout << setw(20) <<  "ca " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ca << endl;
  cout << setw(20) <<  "Lambda " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << Lambda << endl;
  cout << setw(20) <<  "cabi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << cabi << endl;
  cout << setw(20) <<  "cos__cabi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << cos__cabi << endl;
  cout << setw(20) <<  "CKM11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << CKM11 << endl;
  cout << setw(20) <<  "sin__cabi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sin__cabi << endl;
  cout << setw(20) <<  "CKM12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << CKM12 << endl;
  cout << setw(20) <<  "CKM21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << CKM21 << endl;
  cout << setw(20) <<  "CKM22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << CKM22 << endl;
  cout << setw(20) <<  "ca__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << ca__exp__2 << endl;
  cout << setw(20) <<  "sa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sa << endl;
  cout << setw(20) <<  "g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << g1 << endl;
  cout << setw(20) <<  "gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gw << endl;
  cout << setw(20) <<  "MZ__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MZ__exp__2 << endl;
  cout << setw(20) <<  "MZ__exp__4 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MZ__exp__4 << endl;
  cout << setw(20) <<  "sqrt__2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__2 << endl;
  cout << setw(20) <<  "nb__2__exp__0_75 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << nb__2__exp__0_75 << endl;
  cout << setw(20) <<  "MX0__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MX0__exp__2 << endl;
  cout << setw(20) <<  "complexi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << complexi << endl;
  cout << setw(20) <<  "conjg__CKM11 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << conjg__CKM11 << endl;
  cout << setw(20) <<  "conjg__CKM12 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << conjg__CKM12 << endl;
  cout << setw(20) <<  "conjg__CKM21 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << conjg__CKM21 << endl;
  cout << setw(20) <<  "conjg__CKM22 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << conjg__CKM22 << endl;
  cout << setw(20) <<  "aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEW << endl;
  cout << setw(20) <<  "MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MW << endl;
  cout << setw(20) <<  "sqrt__aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__aEW << endl;
  cout << setw(20) <<  "ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ee << endl;
  cout << setw(20) <<  "MW__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MW__exp__2 << endl;
  cout << setw(20) <<  "sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sw2 << endl;
  cout << setw(20) <<  "cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << cw << endl;
  cout << setw(20) <<  "sqrt__sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__sw2 << endl;
  cout << setw(20) <<  "sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sw << endl;
  cout << setw(20) <<  "ad " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ad << endl;
  cout << setw(20) <<  "al " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << al << endl;
  cout << setw(20) <<  "an " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << an << endl;
  cout << setw(20) <<  "au " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << au << endl;
  cout << setw(20) <<  "bd " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << bd << endl;
  cout << setw(20) <<  "bl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << bl << endl;
  cout << setw(20) <<  "bn " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << bn << endl;
  cout << setw(20) <<  "bu " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << bu << endl;
  cout << setw(20) <<  "gwwz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gwwz << endl;
  cout << setw(20) <<  "vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << vev << endl;
  cout << setw(20) <<  "ee__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << ee__exp__2 << endl;
  cout << setw(20) <<  "gAaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gAaa << endl;
  cout << setw(20) <<  "cw__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << cw__exp__2 << endl;
  cout << setw(20) <<  "gAza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gAza << endl;
  cout << setw(20) <<  "gHaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gHaa << endl;
  cout << setw(20) <<  "gHza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gHza << endl;
  cout << setw(20) <<  "vev__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << vev__exp__2 << endl;
  cout << setw(20) <<  "lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << lam << endl;
  cout << setw(20) <<  "yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << yb << endl;
  cout << setw(20) <<  "yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << yt << endl;
  cout << setw(20) <<  "ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ytau << endl;
  cout << setw(20) <<  "muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << muH << endl;
  cout << setw(20) <<  "sw__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << sw__exp__2 << endl;
}
void Parameters_HiggsCharac_UFO::printIndependentCouplings()
{
  cout <<  "HiggsCharac_UFO model couplings independent of event kinematics:"
      << endl;
  cout << setw(20) <<  "GC_61 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_61 << endl;
  cout << setw(20) <<  "GC_89 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_89 << endl;
  cout << setw(20) <<  "GC_35 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_35 << endl;
  cout << setw(20) <<  "GC_37 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_37 << endl;
  cout << setw(20) <<  "GC_56 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_56 << endl;
}
void Parameters_HiggsCharac_UFO::printDependentParameters()
{
  cout <<  "HiggsCharac_UFO model parameters dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "sqrt__aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "G__exp__2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G__exp__2 << endl;
  cout << setw(20) <<  "gAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gAgg << endl;
  cout << setw(20) <<  "gHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gHgg << endl;
}
void Parameters_HiggsCharac_UFO::printDependentCouplings()
{
  cout <<  "HiggsCharac_UFO model couplings dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "GC_53 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_53 << endl;
  cout << setw(20) <<  "GC_18 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_18 << endl;
}


