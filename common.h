#ifndef _common_h_
#define _common_h_

#define _M_PI    3.141592653589793
#define _M_PI_2  9.869604401089358
#define _M_PI_3 31.006276680299816

#ifdef _CUDA_READY_

typedef float afloat_t;

#define _POW_ powf
#define _SQRT_ sqrtf
#define _EXP_ expf
#define _LOG_ logf
#define _ATANH_ atanhf
#define _TANH_ tanhf

#define _CUDA_HOST_ __host__
#define _CUDA_DEVICE_ __device__

#else

#include <cmath>
#include <stdlib.h>

typedef double afloat_t;

#define _POW_ ::pow
#define _SQRT_ std::sqrt
#define _EXP_ std::exp
#define _LOG_ std::log
#define _TANH_ std::tanh
#define _ATANH_ std::atanh

// afloat_t max(afloat_t x, afloat_t y); // { return  x >= y ? x : y; }
// afloat_t min(afloat_t x, afloat_t y); // { return  x <  y ? x : y; }
// int abs(int i); // { return i < 0 ? -i : i; } 

#define _CUDA_HOST_ 
#define _CUDA_DEVICE_ 

#endif

const afloat_t EBEAM = 4000.0;

#endif
