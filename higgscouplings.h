#ifndef _higgs_couplings_
#define _higgs_couplings_

#include "common.h"
#include "dcmplx.h"

struct couplings {
  // Model parameters independent of aS
  afloat_t WX2, WX1, WX0, WW, WZ, WT, ymtau, ymt, ymb, aS, Gf, aEWM1, MX2, MX1,
    MX0, MZ, MTA, MT, MB, kw, kz, ka, kg, kl, kq, kz5, kz3, kz1, kw5, kw4,
    kw3, kw2, kw1, kfb, kfa, kHdw, kHdz, kHda, kAww, kHww, kAzz, kHzz,
    kAgg, kHgg, kAza, kHza, kAaa, kHaa, kAll, kHll, kAbb, kHbb, kAtt, kHtt,
    kSM, ca, Lambda, cabi, cos__cabi, sin__cabi, ca__exp__2, sa, g1, gw,
    MZ__exp__2, MZ__exp__4, sqrt__2, nb__2__exp__0_75, MX0__exp__2, aEW,
    MW, sqrt__aEW, ee, MW__exp__2, sw2, cw, sqrt__sw2, sw, ad, al, an, au,
    bd, bl, bn, bu, gwwz, vev, ee__exp__2, gAaa, cw__exp__2, gAza, gHaa,
    gHza, vev__exp__2, lam, yb, yt, ytau, muH, sw__exp__2;
  dcmplx CKM11, CKM12, CKM21, CKM22, complexi, conjg__CKM11,
    conjg__CKM12, conjg__CKM21, conjg__CKM22;
  // Model parameters dependent on aS
  afloat_t sqrt__aS, G, G__exp__2, gAgg, gHgg; 
  // Model couplings independent of aS
  dcmplx GC_61, GC_89, GC_35, GC_37, GC_56; 
  // Model couplings dependent on aS
  dcmplx GC_53, GC_18; 
} ;  

_CUDA_DEVICE_ _CUDA_HOST_
void init(struct couplings* pars, afloat_t ca, afloat_t kSM, afloat_t kHWW, afloat_t kAWW, afloat_t kHdw, afloat_t lambda);

#endif
