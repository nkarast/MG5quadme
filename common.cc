#include "common.h"

#ifndef _CUDA_READY_
float max(float x, float y) { return  x >= y ? x : y; }
float min(float x, float y) { return  x <  y ? x : y; }
#endif
