#include "higgscouplings.h"

#ifndef __cplusplus
#define M_PI 3.14159265358979323846
#endif

_CUDA_DEVICE_ _CUDA_HOST_
void init(struct couplings* pars, afloat_t ca, afloat_t kSM, afloat_t kHWW, afloat_t kAWW, afloat_t lambda)
{
  /* 
     Note:

     default is ca=1, kSM=1, kHWW=kAWW=0 
  */

  pars->WX2 = 7.567860e-02;
  pars->WX1 = 7.567860e-02;
  pars->WX0 = 7.567860e-02;
  pars->WW = 2.085000e+00;
  pars->WZ = 2.495200e+00;
  pars->WT = 1.508336e+00;
  pars->ymtau = 1.777000e+00;
  pars->ymt = 1.720000e+02;
  pars->ymb = 4.700000e+00;
  pars->aS = 1.184000e-01;
  pars->Gf = 1.166370e-05;
  pars->aEWM1 = 1.279000e+02;
  pars->MX2 = 1.250000e+02;
  pars->MX1 = 1.250000e+02;
  pars->MX0 = 1.250000e+02;
  pars->MZ = 9.118760e+01;
  pars->MTA = 1.777000e+00;
  pars->MT = 1.720000e+02;
  pars->MB = 4.700000e+00;
  pars->kw = 1.000000e+00;
  pars->kz = 1.000000e+00;
  pars->ka = 1.000000e+00;
  pars->kg = 1.000000e+00;
  pars->kl = 1.000000e+00;
  pars->kq = 1.000000e+00;
  pars->kz5 = 0.000000e+00;
  pars->kz3 = 1.000000e+00;
  pars->kz1 = 0.000000e+00;
  pars->kw5 = 0.000000e+00;
  pars->kw4 = 0.000000e+00;
  pars->kw3 = 0.000000e+00;
  pars->kw2 = 1.000000e+00;
  pars->kw1 = 1.000000e+00;
  pars->kfb = 1.000000e+00;
  pars->kfa = 1.000000e+00;
  pars->kHdw = 0.000000e+00;
  pars->kHdz = 0.000000e+00;
  pars->kHda = 0.000000e+00;
  pars->kAww = kAWW;
  pars->kHww = kHWW;
  pars->kAzz = 0.000000e+00;
  pars->kHzz = 0.000000e+00;
  pars->kAgg = 1.000000e+00;
  pars->kHgg = 1.000000e+00;
  pars->kAza = 1.000000e+00;
  pars->kHza = 1.000000e+00;
  pars->kAaa = 1.000000e+00;
  pars->kHaa = 1.000000e+00;
  pars->kAll = 1.000000e+00;
  pars->kHll = 1.000000e+00;
  pars->kAbb = 1.000000e+00;
  pars->kHbb = 1.000000e+00;
  pars->kAtt = 1.000000e+00;
  pars->kHtt = 1.000000e+00;
  pars->kSM = kSM;
  pars->ca = ca;
  pars->Lambda = lambda; // 1.25e+2; // 1.000000e+03;
  pars->cabi = 2.277360e-01;
  pars->cos__cabi = cos(pars->cabi);
  pars->CKM11 = mkdcmplx(pars->cos__cabi,0); 
  pars->sin__cabi = sin(pars->cabi); 
  pars->CKM12 = mkdcmplx(pars->sin__cabi,0); 
  pars->CKM21 = mkdcmplx(-pars->sin__cabi,0);
  pars->CKM22 = mkdcmplx(pars->cos__cabi,0); 
  pars->ca__exp__2 = _POW_(pars->ca, 2.); 
  pars->sa = sqrt(1. - pars->ca__exp__2); 
  pars->g1 = 1.; 
  pars->gw = 1.; 
  pars->MZ__exp__2 = _POW_(pars->MZ, 2.); 
  pars->MZ__exp__4 = _POW_(pars->MZ, 4.); 
  pars->sqrt__2 = sqrt(2.); 
  pars->nb__2__exp__0_75 = _POW_(2., 0.75); 
  pars->MX0__exp__2 = _POW_(pars->MX0, 2.); 
  pars->complexi = mkdcmplx(0., 1.); 
  pars->conjg__CKM11 = conj(pars->CKM11); 
  pars->conjg__CKM12 = conj(pars->CKM12); 
  pars->conjg__CKM21 = conj(pars->CKM21); 
  pars->conjg__CKM22 = conj(pars->CKM22); 
  pars->aEW = 1./pars->aEWM1; 
  pars->MW = sqrt(pars->MZ__exp__2/2. + sqrt(pars->MZ__exp__4/4. - (pars->aEW * M_PI * pars->MZ__exp__2)/(pars->Gf * pars->sqrt__2)));
  pars->sqrt__aEW = sqrt(pars->aEW); 
  pars->ee = 2. * pars->sqrt__aEW * sqrt(M_PI); 
  pars->MW__exp__2 = _POW_(pars->MW, 2.); 
  pars->sw2 = 1. - pars->MW__exp__2/pars->MZ__exp__2; 
  pars->cw = sqrt(1. - pars->sw2); 
  pars->sqrt__sw2 = sqrt(pars->sw2); 
  pars->sw = pars->sqrt__sw2; 
  pars->ad = (pars->ee * (-0.5 + (2. * pars->sw2)/3.))/(2. * pars->cw * pars->sw); 
  pars->al = (pars->ee * (-0.5 + 2. * pars->sw2))/(2. * pars->cw * pars->sw); 
  pars->an = pars->ee/(4. * pars->cw * pars->sw); 
  pars->au = (pars->ee * (0.5 - (4. * pars->sw2)/3.))/(2. * pars->cw * pars->sw); 
  pars->bd = -pars->ee/(4. * pars->cw * pars->sw); 
  pars->bl = -pars->ee/(4. * pars->cw * pars->sw); 
  pars->bn = pars->ee/(4. * pars->cw * pars->sw); 
  pars->bu = pars->ee/(4. * pars->cw * pars->sw); 
  pars->gwwz = -((pars->cw * pars->ee)/pars->sw); 
  pars->vev = (2. * pars->MW * pars->sw)/pars->ee; 
  pars->ee__exp__2 = _POW_(pars->ee, 2.); 
  pars->gAaa = -pars->ee__exp__2/(3. * _POW_(M_PI, 2.) * pars->vev); 
  pars->cw__exp__2 = _POW_(pars->cw, 2.); 
  pars->gAza = -((-5. + 8. * pars->cw__exp__2) * sqrt(pars->ee__exp__2 * pars->Gf * pars->MZ__exp__2))/(6. *  pars->nb__2__exp__0_75 * _POW_(M_PI, 2.) * pars->vev);
  pars->gHaa = (47. * pars->ee__exp__2)/(72. * _POW_(M_PI, 2.) * pars->vev); 
  pars->gHza = ((-13. + 94. * pars->cw__exp__2) * sqrt(pars->ee__exp__2 * pars->Gf * pars->MZ__exp__2))/(36. * pars->nb__2__exp__0_75 * _POW_(M_PI, 2.) * pars->vev);
  pars->vev__exp__2 = _POW_(pars->vev, 2.); 
  pars->lam = pars->MX0__exp__2/(2. * pars->vev__exp__2); 
  pars->yb = (pars->ymb * pars->sqrt__2)/pars->vev; 
  pars->yt = (pars->ymt * pars->sqrt__2)/pars->vev; 
  pars->ytau = (pars->ymtau * pars->sqrt__2)/pars->vev; 
  pars->muH = sqrt(pars->lam * pars->vev__exp__2); 
  pars->sw__exp__2 = _POW_(pars->sw, 2.); 
  pars->GC_61 = (pars->ee * pars->complexi)/(pars->sw * pars->sqrt__2); 
  pars->GC_89 = (ca * pars->ee__exp__2 * pars->complexi * kSM * pars->vev)/(2. * pars->sw__exp__2); 
  pars->GC_35 = -((ca * pars->complexi * pars->kHdw)/lambda); 
  pars->GC_37 = -((ca * pars->complexi * pars->kHww)/lambda); 
  pars->GC_56 = (pars->complexi * pars->kAww * pars->sa)/lambda; 
  pars->sqrt__aS = sqrt(pars->aS); 
  pars->G = 2. * pars->sqrt__aS * sqrt(M_PI); 
  pars->G__exp__2 = _POW_(pars->G, 2.); 
  pars->gAgg = -pars->G__exp__2/(8. * _POW_(M_PI, 2.) * pars->vev); 
  pars->gHgg = -pars->G__exp__2/(12. * _POW_(M_PI, 2.) * pars->vev); 
  pars->GC_53 = (pars->complexi * pars->gAgg * pars->kAgg * pars->sa)/8.; 
  pars->GC_18 = -(ca * pars->complexi * pars->gHgg * pars->kHgg); 
}
