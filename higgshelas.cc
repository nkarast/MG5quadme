#include "higgshelas.h"

#ifndef __CINT__
#ifndef __cplusplus
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define abs(a) \
  ({ __typeof__ (a) _a = (a); \
     _a < 0 ? -_a : _a; })
#endif
#endif

_CUDA_DEVICE_ _CUDA_HOST_
afloat_t sgn(afloat_t a, afloat_t b)
{
  return(b < 0)? -fabs(a) : fabs(a); 
}

_CUDA_DEVICE_ _CUDA_HOST_
void vxxxxx(afloat_t p[4], afloat_t vmass, int nhel, int nsv, dcmplx vc[6])
{
  afloat_t hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = _POW_(0.5, 0.5); 
  hel = afloat_t(nhel); 
  nsvahl = nsv * fabs(hel); 
  pt2 = _POW_(p[1], 2) + _POW_(p[2], 2); 
  pp = min(p[0], _POW_(pt2 + _POW_(p[3], 2), 0.5)); 
  pt = min(pp, _POW_(pt2, 0.5)); 
  vc[0] = mkdcmplx(p[0] * nsv, p[3] * nsv); 
  vc[1] = mkdcmplx(p[1] * nsv, p[2] * nsv); 
  if(vmass != 0.0)
  {
    hel0 = 1.0 - fabs(hel); 
    if(pp == 0.0)
    {
      vc[2] = mkdcmplx(0.0, 0.0); 
      vc[3] = mkdcmplx(-hel * sqh, 0.0); 
      vc[4] = mkdcmplx(0.0, nsvahl * sqh); 
      vc[5] = mkdcmplx(hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = mkdcmplx(hel0 * pp/vmass, 0.0); 
      vc[5] = mkdcmplx(hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if(pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = mkdcmplx(hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = mkdcmplx(hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = mkdcmplx(-hel * sqh, 0.0); 
        vc[4] = mkdcmplx(0.0, nsvahl * sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = _POW_(_POW_(p[1], 2) + _POW_(p[2], 2), 0.5); 
    vc[2] = mkdcmplx(0.0, 0.0); 
    vc[5] = mkdcmplx(hel * pt/pp * sqh, 0.0); 
    if(pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = mkdcmplx(-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = mkdcmplx(-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = mkdcmplx(-hel * sqh, 0.0); 
      vc[4] = mkdcmplx(0.0, nsv * sgn(sqh, p[3])); 
    }
  }
  return; 
}

_CUDA_DEVICE_ _CUDA_HOST_
void oxxxxx(afloat_t p[4], afloat_t fmass, int nhel, int nsf, dcmplx fo[6])
{
  dcmplx chi[2]; 
  afloat_t sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = mkdcmplx(p[0] * nsf, p[3] * nsf); 
  fo[1] = mkdcmplx(p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if(fmass != 0.000)
  {
    pp = min(p[0], _POW_(_POW_(p[1], 2) + _POW_(p[2], 2) + _POW_(p[3], 2), 0.5)); 
    if(pp == 0.000)
    {
      sqm[0] = _POW_(fabs(fmass), 0.5); 
      sqm[1] = sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im =(1 + nh)/2 * nhel; 
      fo[2] = mkdcmplx(im * sqm[abs(ip)]); 
      fo[3] = mkdcmplx(ip * nsf * sqm[abs(ip)]); 
      fo[4] = mkdcmplx(im * nsf * sqm[abs(im)]); 
      fo[5] = mkdcmplx(ip * sqm[abs(im)]); 
    }
    else
    {
      pp = min(p[0], _POW_(_POW_(p[1], 2) + _POW_(p[2], 2) + _POW_(p[3], 2), 0.5)); 
      sf[0] = afloat_t(1 + nsf +(1 - nsf) * nh) * 0.5; 
      sf[1] = afloat_t(1 + nsf -(1 - nsf) * nh) * 0.5; 
      omega[0] = _POW_(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip =(1 + nh)/2; 
      im =(1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], (afloat_t)0.00); 
      chi[0] = mkdcmplx(_POW_(pp3 * 0.5/pp, 0.5), 0.00); 
      if(pp3 == 0.00)
      {
        chi[1] = mkdcmplx(-nh, 0.00); 
      }
      else
      {
        chi[1] = mkdcmplx(nh * p[1], -p[2])/_POW_(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and(p[2] == 0.00) and(p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = _POW_(max(p[0] + p[3], (afloat_t)0.00), 0.5) * nsf; 
    }
    chi[0] = mkdcmplx(sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = mkdcmplx(-nhel, 0.00) * _POW_(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = mkdcmplx(nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = mkdcmplx(0.00, 0.00); 
      fo[5] = mkdcmplx(0.00, 0.00); 
    }
    else
    {
      fo[2] = mkdcmplx(0.00, 0.00); 
      fo[3] = mkdcmplx(0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

_CUDA_DEVICE_ _CUDA_HOST_
void ixxxxx(afloat_t p[4], afloat_t fmass, int nhel, int nsf, dcmplx fi[6])
{
  dcmplx chi[2]; 
  afloat_t sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = mkdcmplx(-p[0] * nsf, -p[3] * nsf); 
  fi[1] = mkdcmplx(-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if(fmass != 0.0)
  {
    pp = min(p[0], _POW_((_POW_(p[1], 2) + _POW_(p[2], 2) + _POW_(p[3], 2)), 0.5)); 
    if(pp == 0.0)
    {
      sqm[0] = _POW_(fabs(fmass), 0.5); 
      sqm[1] = sgn(sqm[0], fmass); 
      ip =(1 + nh)/2; 
      im =(1 - nh)/2; 
      fi[2] = mkdcmplx(ip * sqm[ip]); 
      fi[3] = mkdcmplx(im * nsf * sqm[ip]); 
      fi[4] = mkdcmplx(ip * nsf * sqm[im]); 
      fi[5] = mkdcmplx(im * sqm[im]); 
    }
    else
    {
      sf[0] =(1 + nsf +(1 - nsf) * nh) * 0.5; 
      sf[1] =(1 + nsf -(1 - nsf) * nh) * 0.5; 
      omega[0] = _POW_(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip =(1 + nh)/2; 
      im =(1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], (afloat_t)0.0); 
      chi[0] = mkdcmplx(_POW_(pp3 * 0.5/pp, 0.5), 0); 
      if(pp3 == 0.0)
      {
        chi[1] = mkdcmplx(-nh, 0); 
      }
      else
      {
        chi[1] = mkdcmplx(nh * p[1], p[2])/_POW_(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if(p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = _POW_(max(p[0] + p[3], (afloat_t)0.0), (afloat_t)0.5) * nsf; 
    }
    chi[0] = mkdcmplx(sqp0p3, 0.0); 
    if(sqp0p3 == 0.0)
    {
      chi[1] = mkdcmplx(-nhel * _POW_(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = mkdcmplx(nh * p[1], p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fi[2] = mkdcmplx(0.0, 0.0); 
      fi[3] = mkdcmplx(0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = mkdcmplx(0.0, 0.0); 
      fi[5] = mkdcmplx(0.0, 0.0); 
    }
  }
  return; 
}

_CUDA_DEVICE_ _CUDA_HOST_
void sxxxxx(afloat_t p[4], int nss, dcmplx sc[3])
{
  sc[2] = mkdcmplx(1.00, 0.00); 
  sc[0] = mkdcmplx(p[0] * nss, p[3] * nss); 
  sc[1] = mkdcmplx(p[1] * nss, p[2] * nss); 
  return; 
}


_CUDA_DEVICE_ _CUDA_HOST_
void txxxxx(afloat_t p[4], afloat_t tmass, int nhel, int nst, dcmplx tc[18])
{
  dcmplx ft[6][4], ep[4], em[4], e0[4]; 
  afloat_t pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = _POW_(0.5, 0.5); 
  sqs = _POW_(0.5/3, 0.5); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], _POW_(pt2 + p[3] * p[3], 0.5)); 
  pt = min(pp, _POW_(pt2, 0.5)); 

  ft[4][0] = mkdcmplx(p[0] * nst, p[3] * nst); 
  ft[5][0] = mkdcmplx(p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = mkdcmplx(0, 0); 
      ep[1] = mkdcmplx(-sqh, 0); 
      ep[2] = mkdcmplx(0, nst * sqh); 
      ep[3] = mkdcmplx(0, 0); 
    }
    else
    {
      ep[0] = mkdcmplx(0, 0); 
      ep[3] = mkdcmplx(pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = mkdcmplx(-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = mkdcmplx(-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = mkdcmplx(-sqh, 0); 
        ep[2] = mkdcmplx(0, nst * sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = mkdcmplx(0, 0); 
      em[1] = mkdcmplx(sqh, 0); 
      em[2] = mkdcmplx(0, nst * sqh); 
      em[3] = mkdcmplx(0, 0); 
    }
    else
    {
      em[0] = mkdcmplx(0, 0); 
      em[3] = mkdcmplx(-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = mkdcmplx(-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = mkdcmplx(-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = mkdcmplx(sqh, 0); 
        em[2] = mkdcmplx(0, nst * sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(abs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = mkdcmplx(0, 0); 
      e0[1] = mkdcmplx(0, 0); 
      e0[2] = mkdcmplx(0, 0); 
      e0[3] = mkdcmplx(1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = mkdcmplx(pp/tmass, 0); 
      e0[3] = mkdcmplx(p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = mkdcmplx(p[1] * emp, 0); 
        e0[2] = mkdcmplx(p[2] * emp, 0); 
      }
      else
      {
        e0[1] = mkdcmplx(0, 0); 
        e0[2] = mkdcmplx(0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++)
    {
      for(i = 0; i < 4; i++)
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++)
    {
      for(i = 0; i < 4; i++)
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++)
    {
      for(i = 0; i < 4; i++)
        ft[i][j] = mkdcmplx(0,0); 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++)
      {
        for(i = 0; i < 4; i++)
          ft[i][j] = sqh *(ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++)
      {
        for(i = 0; i < 4; i++)
          ft[i][j] = sqs *(ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++)
      {
        for(i = 0; i < 4; i++)
          ft[i][j] = sqh *(em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      // std::cerr <<  "invalid helicity in txxxxx.\n"; 
      // std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++)
  {
    for(i = 0; i < 4; i++)
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}
_CUDA_DEVICE_ _CUDA_HOST_
void vvs2_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup, dcmplx* vertex)
{
  dcmplx ci = mkdcmplx(0., 1.); 
  dcmplx tmp2; 
  dcmplx tmp1; 
  afloat_t p1[4]; 
  afloat_t p2[4]; 
  p1[0] = v1[0].re; 
  p1[1] = v1[1].re; 
  p1[2] = v1[1].im; 
  p1[3] = v1[0].im; 
  p2[0] = v2[0].re; 
  p2[1] = v2[1].re; 
  p2[2] = v2[1].im; 
  p2[3] = v2[0].im; 
  tmp1 = -1. *(p1[0] *(p2[1] *(v2[5] * v1[4] - v2[4] * v1[5]) +(p2[2] *
     (v2[3] * v1[5] - v2[5] * v1[3]) + p2[3] *(v2[4] * v1[3] - v2[3] *
      v1[4]))) +(p1[1] *(p2[0] *(v2[4] * v1[5] - v2[5] * v1[4]) +(p2[2] *
     (v2[5] * v1[2] - v2[2] * v1[5]) + p2[3] *(v2[2] * v1[4] - v2[4] *
      v1[2]))) +(p1[2] *(p2[0] *(v2[5] * v1[3] - v2[3] * v1[5]) +(p2[1] *
     (v2[2] * v1[5] - v2[5] * v1[2]) + p2[3] *(v2[3] * v1[2] - v2[2] *
      v1[3]))) + p1[3] *(p2[0] *(v2[3] * v1[4] - v2[4] * v1[3]) +(p2[1] *
     (v2[4] * v1[2] - v2[2] * v1[4]) + p2[2] *(v2[2] * v1[3] - v2[3] *
      v1[2]))))));
  tmp2 = -1. *(p1[0] *(p2[1] *(v2[4] * v1[5] - v2[5] * v1[4]) +(p2[2] *
     (v2[5] * v1[3] - v2[3] * v1[5]) + p2[3] *(v2[3] * v1[4] - v2[4] *
      v1[3]))) +(p1[1] *(p2[0] *(v2[5] * v1[4] - v2[4] * v1[5]) +(p2[2] *
     (v2[2] * v1[5] - v2[5] * v1[2]) + p2[3] *(v2[4] * v1[2] - v2[2] *
      v1[4]))) +(p1[2] *(p2[0] *(v2[3] * v1[5] - v2[5] * v1[3]) +(p2[1] *
     (v2[5] * v1[2] - v2[2] * v1[5]) + p2[3] *(v2[2] * v1[3] - v2[3] *
      v1[2]))) + p1[3] *(p2[0] *(v2[4] * v1[3] - v2[3] * v1[4]) +(p2[1] *
     (v2[2] * v1[4] - v2[4] * v1[2]) + p2[2] *(v2[3] * v1[2] - v2[2] *
      v1[3]))))));
  (*vertex) = coup * 4. * s3[2] *(-ci *(tmp1) + ci *(tmp2)); 
}

_CUDA_DEVICE_ _CUDA_HOST_
void vvs2_6_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup1, dcmplx coup2, dcmplx* vertex)
{
  dcmplx ci = mkdcmplx(0., 1.); 
  afloat_t p1[4]; 
  afloat_t p2[4]; 
  dcmplx tmp; 
  vvs2_0(v1, v2, s3, coup1, vertex); 
  vvs6_0(v1, v2, s3, coup2, &tmp); 
  (*vertex) = (*vertex) + tmp; 
}

_CUDA_DEVICE_ _CUDA_HOST_
void vvs1_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  dcmplx tmp1; 
  afloat_t p1[4]; 
  afloat_t p2[4]; 
  afloat_t p3[4]; 
  dcmplx denom; 
  p1[0] = v1[0].re; 
  p1[1] = v1[1].re; 
  p1[2] = v1[1].im; 
  p1[3] = v1[0].im; 
  p2[0] = v2[0].re; 
  p2[1] = v2[1].re; 
  p2[2] = v2[1].im; 
  p2[3] = v2[0].im; 
  s3[0] = +v1[0] + v2[0]; 
  s3[1] = +v1[1] + v2[1]; 
  p3[0] = -s3[0].re; 
  p3[1] = -s3[1].re; 
  p3[2] = -s3[1].im; 
  p3[3] = -s3[0].im; 
  tmp1 = -1. *(p1[0] *(p2[1] *(v2[5] * v1[4] - v2[4] * v1[5]) +(p2[2] *
     (v2[3] * v1[5] - v2[5] * v1[3]) + p2[3] *(v2[4] * v1[3] - v2[3] *
      v1[4]))) +(p1[1] *(p2[0] *(v2[4] * v1[5] - v2[5] * v1[4]) +(p2[2] *
     (v2[5] * v1[2] - v2[2] * v1[5]) + p2[3] *(v2[2] * v1[4] - v2[4] *
      v1[2]))) +(p1[2] *(p2[0] *(v2[5] * v1[3] - v2[3] * v1[5]) +(p2[1] *
     (v2[2] * v1[5] - v2[5] * v1[2]) + p2[3] *(v2[3] * v1[2] - v2[2] *
      v1[3]))) + p1[3] *(p2[0] *(v2[3] * v1[4] - v2[4] * v1[3]) +(p2[1] *
     (v2[4] * v1[2] - v2[2] * v1[4]) + p2[2] *(v2[2] * v1[3] - v2[3] *
      v1[2]))))));
  denom = coup/(_POW_(p3[0], 2) - _POW_(p3[1], 2) - _POW_(p3[2], 2) - _POW_(p3[3], 2) -
		m3 *(mkdcmplx(m3) - ci * w3));
  s3[2] = denom * ci * tmp1; 
}

_CUDA_DEVICE_ _CUDA_HOST_
void vvs1_6_7_3(dcmplx v1[], dcmplx v2[], dcmplx coup1, dcmplx coup2, dcmplx coup3, afloat_t m3, afloat_t w3, dcmplx s3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  afloat_t p1[4]; 
  afloat_t p2[4]; 
  afloat_t p3[4]; 
  dcmplx denom; 
  dcmplx stmp[3]; 
  int i; 
  vvs1_3(v1, v2, coup1, m3, w3, s3); 
  vvs6_3(v1, v2, coup2, m3, w3, stmp); 
  i = 2; 
  while(i < 3)
  {
    s3[i] = s3[i] + stmp[i]; 
    i++; 
  }
  vvs7_3(v1, v2, coup3, m3, w3, stmp); 
  i = 2; 
  while(i < 3)
  {
    s3[i] = s3[i] + stmp[i]; 
    i++; 
  }
}

_CUDA_DEVICE_ _CUDA_HOST_
void vvs4_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  afloat_t p3[4]; 
  dcmplx denom; 
  dcmplx tmp3; 
  s3[0] = +v1[0] + v2[0]; 
  s3[1] = +v1[1] + v2[1]; 
  p3[0] = -s3[0].re; 
  p3[1] = -s3[1].re; 
  p3[2] = -s3[1].im; 
  p3[3] = -s3[0].im; 
  tmp3 =(v2[2] * v1[2] - v2[3] * v1[3] - v2[4] * v1[4] - v2[5] * v1[5]); 
  denom = coup/(_POW_(p3[0], 2) - _POW_(p3[1], 2) - _POW_(p3[2], 2) - _POW_(p3[3], 2) -
		m3 *(mkdcmplx(m3) - ci * w3));
  s3[2] = denom * ci * tmp3; 
}


_CUDA_DEVICE_ _CUDA_HOST_
void vvs7_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  dcmplx denom; 
  afloat_t p1[4]; 
  afloat_t p2[4]; 
  dcmplx tmp7; 
  afloat_t p3[4]; 
  dcmplx tmp6; 
  dcmplx tmp5; 
  dcmplx tmp4; 
  dcmplx tmp9; 
  dcmplx tmp3; 
  dcmplx tmp8; 
  p1[0] = v1[0].re; 
  p1[1] = v1[1].re; 
  p1[2] = v1[1].im; 
  p1[3] = v1[0].im; 
  p2[0] = v2[0].re; 
  p2[1] = v2[1].re; 
  p2[2] = v2[1].im; 
  p2[3] = v2[0].im; 
  s3[0] = +v1[0] + v2[0]; 
  s3[1] = +v1[1] + v2[1]; 
  p3[0] = -s3[0].re; 
  p3[1] = -s3[1].re; 
  p3[2] = -s3[1].im; 
  p3[3] = -s3[0].im; 
  tmp9 =mkdcmplx(p2[0] * p2[0] - p2[1] * p2[1] - p2[2] * p2[2] - p2[3] * p2[3],0); 
  tmp8 =mkdcmplx(p1[0] * p1[0] - p1[1] * p1[1] - p1[2] * p1[2] - p1[3] * p1[3],0); 
  tmp5 =(p1[0] * v1[2] - p1[1] * v1[3] - p1[2] * v1[4] - p1[3] * v1[5]); 
  tmp4 =(v2[2] * p1[0] - v2[3] * p1[1] - v2[4] * p1[2] - v2[5] * p1[3]); 
  tmp7 =(p2[0] * v1[2] - p2[1] * v1[3] - p2[2] * v1[4] - p2[3] * v1[5]); 
  tmp6 =(v2[2] * p2[0] - v2[3] * p2[1] - v2[4] * p2[2] - v2[5] * p2[3]); 
  tmp3 =(v2[2] * v1[2] - v2[3] * v1[3] - v2[4] * v1[4] - v2[5] * v1[5]); 
  denom = coup/(_POW_(p3[0], 2) - _POW_(p3[1], 2) - _POW_(p3[2], 2) - _POW_(p3[3], 2) -
      m3 *(m3 - ci * w3));
  s3[2] = denom *(tmp3 * - 1. *(+ci *(tmp8 + tmp9)) +(+ci *(tmp4 * tmp5 +
      tmp6 * tmp7)));
}


_CUDA_DEVICE_ _CUDA_HOST_
void vvs6_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  afloat_t p1[4]; 
  dcmplx tmp10; 
  afloat_t p2[4]; 
  dcmplx tmp7; 
  afloat_t p3[4]; 
  dcmplx denom; 
  dcmplx tmp4; 
  dcmplx tmp3; 
  p1[0] = v1[0].re; 
  p1[1] = v1[1].re; 
  p1[2] = v1[1].im; 
  p1[3] = v1[0].im; 
  p2[0] = v2[0].re; 
  p2[1] = v2[1].re; 
  p2[2] = v2[1].im; 
  p2[3] = v2[0].im; 
  s3[0] = +v1[0] + v2[0]; 
  s3[1] = +v1[1] + v2[1]; 
  p3[0] = -s3[0].re; 
  p3[1] = -s3[1].re; 
  p3[2] = -s3[1].im; 
  p3[3] = -s3[0].im; 
  tmp4 =(v2[2] * p1[0] - v2[3] * p1[1] - v2[4] * p1[2] - v2[5] * p1[3]); 
  tmp7 =(p2[0] * v1[2] - p2[1] * v1[3] - p2[2] * v1[4] - p2[3] * v1[5]); 
  tmp10 =mkdcmplx(p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2] - p1[3] * p2[3],0); 
  tmp3 =(v2[2] * v1[2] - v2[3] * v1[3] - v2[4] * v1[4] - v2[5] * v1[5]); 
  denom = coup/(_POW_(p3[0], 2) - _POW_(p3[1], 2) - _POW_(p3[2], 2) - _POW_(p3[3], 2) -
      m3 *(m3 - ci * w3));
  s3[2] = denom *(-ci *(tmp3 * tmp10) + ci *(tmp4 * tmp7)); 
}


_CUDA_DEVICE_ _CUDA_HOST_
void ffv3_3(dcmplx f1[], dcmplx f2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx v3[])
{
  dcmplx ci = mkdcmplx(0., 1.); 
  dcmplx denom; 
  dcmplx tmp0; 
  afloat_t p3[4]; 
  afloat_t om3; 
  om3 = 0.; 
  if(m3 != 0.)
    om3 = 1./_POW_(m3, 2); 
  v3[0] = +f1[0] + f2[0]; 
  v3[1] = +f1[1] + f2[1]; 
  p3[0] = -v3[0].re; 
  p3[1] = -v3[1].re; 
  p3[2] = -v3[1].im; 
  p3[3] = -v3[0].im; 
  tmp0 =(f1[2] *(f2[4] *(p3[0] + p3[3]) + f2[5] *(p3[1] + ci *(p3[2]))) +
      f1[3] *(f2[4] *(p3[1] - ci *(p3[2])) + f2[5] *(p3[0] - p3[3])));
  denom = coup/(_POW_(p3[0], 2) - _POW_(p3[1], 2) - _POW_(p3[2], 2) - _POW_(p3[3], 2) -
      m3 *(m3 - ci * w3));
  v3[2] = denom * - ci *(f2[4] * f1[2] + f2[5] * f1[3] - p3[0] * om3 * tmp0); 
  v3[3] = denom * - ci *(-f2[5] * f1[2] - f2[4] * f1[3] - p3[1] * om3 * tmp0); 
  v3[4] = denom * - ci *(-ci *(f2[5] * f1[2]) + ci *(f2[4] * f1[3]) - p3[2]
      * om3 * tmp0);
  v3[5] = denom * - ci *(f2[5] * f1[3] - f2[4] * f1[2] - p3[3] * om3 * tmp0); 
}


_CUDA_DEVICE_ _CUDA_HOST_
void vvs6_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup, dcmplx* vertex)
{
  dcmplx ci = mkdcmplx(0., 1.); 
  afloat_t p1[4]; 
  dcmplx tmp10; 
  afloat_t p2[4]; 
  dcmplx tmp7; 
  dcmplx tmp4; 
  dcmplx tmp3; 
  p1[0] = v1[0].re; 
  p1[1] = v1[1].re; 
  p1[2] = v1[1].im; 
  p1[3] = v1[0].im; 
  p2[0] = v2[0].re; 
  p2[1] = v2[1].re; 
  p2[2] = v2[1].im; 
  p2[3] = v2[0].im; 
  tmp4 =(v2[2] * p1[0] - v2[3] * p1[1] - v2[4] * p1[2] - v2[5] * p1[3]); 
  tmp7 =(p2[0] * v1[2] - p2[1] * v1[3] - p2[2] * v1[4] - p2[3] * v1[5]); 
  tmp10 =mkdcmplx(p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2] - p1[3] * p2[3],0); 
  tmp3 =(v2[2] * v1[2] - v2[3] * v1[3] - v2[4] * v1[4] - v2[5] * v1[5]); 
  (*vertex) = coup * s3[2] *(-ci *(tmp4 * tmp7) + ci *(tmp3 * tmp10)); 
}
