
#ifndef _higgs_helas_
#define _higgs_helas_

// include <cmath> 

// #include <complex> 
// typedef std::complex<afloat_t> dcmplx;

#include "common.h"
#include "dcmplx.h"

_CUDA_DEVICE_ _CUDA_HOST_ afloat_t sgn(afloat_t e, afloat_t f); 

_CUDA_DEVICE_ _CUDA_HOST_ void oxxxxx(afloat_t p[4], afloat_t fmass, int nhel, int nsf, dcmplx fo[6]);
_CUDA_DEVICE_ _CUDA_HOST_ void txxxxx(afloat_t p[4], afloat_t tmass, int nhel, int nst, dcmplx fi[18]);
_CUDA_DEVICE_ _CUDA_HOST_ void ixxxxx(afloat_t p[4], afloat_t fmass, int nhel, int nsf, dcmplx fi[6]);
_CUDA_DEVICE_ _CUDA_HOST_ void sxxxxx(afloat_t p[4], int nss, dcmplx sc[3]); 
_CUDA_DEVICE_ _CUDA_HOST_ void vxxxxx(afloat_t p[4], afloat_t vmass, int nhel, int nsv, dcmplx v[6]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs2_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup, dcmplx* vertex);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs2_6_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup1, dcmplx coup2, dcmplx* vertex);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs1_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs1_6_7_3(dcmplx v1[], dcmplx v2[], dcmplx coup1, dcmplx coup2, dcmplx coup3, afloat_t m3, afloat_t w3, dcmplx s3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs4_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs7_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs6_3(dcmplx v1[], dcmplx v2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx s3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void ffv3_3(dcmplx f1[], dcmplx f2[], dcmplx coup, afloat_t m3, afloat_t w3, dcmplx v3[]);
_CUDA_DEVICE_ _CUDA_HOST_ void vvs6_0(dcmplx v1[], dcmplx v2[], dcmplx s3[], dcmplx coup, dcmplx* vertex);

#endif  
