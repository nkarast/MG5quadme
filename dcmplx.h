#ifndef __cmplx_h__
#define __cmplx_h__

#include "common.h"

#ifdef _CUDA_READY_
typedef struct __align__(8)
{
  afloat_t re; afloat_t im;
} dcmplx;
#else
typedef struct {
  afloat_t re; afloat_t im;
} dcmplx;
#endif

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx mkdcmplx(afloat_t re, afloat_t im)
{
  dcmplx z; z.re=re; z.im=im;
  return z;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
afloat_t real(dcmplx a)
{
  return a.re;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
afloat_t imag(dcmplx a)
{
  return a.im;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx conj(dcmplx a)
{
  return mkdcmplx(a.re,-a.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
afloat_t fabsc(dcmplx a)
{
  return _SQRT_((a.re*a.re)+(a.im*a.im));
}

inline _CUDA_DEVICE_ _CUDA_HOST_
afloat_t fabsc_sqr(dcmplx a)
{
  return (a.re*a.re)+(a.im*a.im);
}

#ifdef __cplusplus // block operator and function overloading for ANSI C compliance

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx mkdcmplx(dcmplx z)
{
  return mkdcmplx(z.re,z.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx mkdcmplx(afloat_t a)
{
  return mkdcmplx(a, 0.);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator+(dcmplx a, dcmplx b)
{
  return mkdcmplx(a.re + b.re, a.im + b.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator+(afloat_t a, dcmplx b)
{
  return mkdcmplx(a + b.re, b.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
void operator+=(dcmplx &a, dcmplx b)
{
  a.re += b.re; a.im += b.im;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator+(dcmplx a)
{
  return mkdcmplx(+a.re, +a.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator-(dcmplx a, dcmplx b)
{
  return mkdcmplx(a.re - b.re, a.im - b.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator-(afloat_t a, dcmplx b)
{
  return mkdcmplx(a - b.re, -b.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
void operator-=(dcmplx &a, dcmplx b)
{
  a.re -= b.re; a.im -= b.im;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator-(dcmplx a)
{
  return mkdcmplx(-a.re, -a.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator*(dcmplx a, dcmplx b)
{
  return mkdcmplx((a.re * b.re) - (a.im * b.im),
		  (a.re * b.im) + (a.im * b.re));
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator*(dcmplx a, afloat_t s)
{
  return mkdcmplx(a.re * s, a.im * s);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator*(afloat_t s, dcmplx a)
{
  return mkdcmplx(a.re * s, a.im * s);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
void operator*=(dcmplx &a, afloat_t s)
{
  a.re *= s; a.im *= s;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator/(dcmplx a, dcmplx b)
{
  afloat_t tmpD=(1./(b.re*b.re+b.im*b.im));
  return mkdcmplx( ( (a.re * b.re) + (a.im * b.im))*tmpD,
		   (-(a.re * b.im) + (a.im * b.re))*tmpD );
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator/(dcmplx a, afloat_t s)
{
  afloat_t inv = 1. / s;
  return a * inv;
}

inline _CUDA_DEVICE_ _CUDA_HOST_
dcmplx operator/(afloat_t s, dcmplx a)
{
  afloat_t inv = s*(1./(a.re*a.re+a.im*a.im));
  return mkdcmplx(inv*a.re,-inv*a.im);
}

inline _CUDA_DEVICE_ _CUDA_HOST_
void operator/=(dcmplx &a, afloat_t s)
{
  afloat_t inv = 1. / s;
  a *= inv;
}

#endif

#endif
