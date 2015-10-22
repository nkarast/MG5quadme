
#define _CUDA_READY_

#include "common.h"

#include "higgs.cc"
#include "higgscouplings.cc"
#include "higgshelas.cc"
#include "pdf.cc"

extern "C"
{
  __global__ void quadme(
			  afloat_t* result, afloat_t* input, afloat_t* grd, int* idxdim, int ndim, /* phase space data */
			  afloat_t* gpdf_data, int gpdf_nptsx, int gpdf_nptsq, afloat_t gpdf_lxmin, afloat_t gpdf_lxmax, afloat_t gpdf_lqmin, afloat_t gpdf_lqmax, /* gluon PDF data table */
			  afloat_t ca, afloat_t kSM, afloat_t kHWW, afloat_t kAWW, afloat_t lambda /* Higgs characterization model parameters */
			 )
  {
    // __shared__ struct couplings _couplings;
    // if(threadIdx.x==0)
    // {
    //   init(&_couplings, ca, kSM, kHWW, kAWW, lambda);
    // }
    // __syncthreads();

    struct couplings _couplings;
    init(&_couplings, ca, kSM, kHWW, kAWW, lambda);
    
    int idxthr = threadIdx.x + blockIdx.x * blockDim.x;
     
    const uint nkin = 4*nexternal;

    const uint IP0_E = 4*0; const uint IP0_X = 4*0+1; const uint IP0_Y = 4*0+2; const uint IP0_Z = 4*0+3;
    const uint IP1_E = 4*1; const uint IP1_X = 4*1+1; const uint IP1_Y = 4*1+2; const uint IP1_Z = 4*1+3;
    const uint IP2_E = 4*2; const uint IP2_X = 4*2+1; const uint IP2_Y = 4*2+2; const uint IP2_Z = 4*2+3;
    const uint IP3_E = 4*3; const uint IP3_X = 4*3+1; const uint IP3_Y = 4*3+2; const uint IP3_Z = 4*3+3;
    const uint IP4_E = 4*4; const uint IP4_X = 4*4+1; const uint IP4_Y = 4*4+2; const uint IP4_Z = 4*4+3;
    const uint IP5_E = 4*5; const uint IP5_X = 4*5+1; const uint IP5_Y = 4*5+2; const uint IP5_Z = 4*5+3;
    
    afloat_t p[nkin];
    
    if(ndim==0)
    {
#pragma unroll
      for(unsigned int i=0; i<nkin; ++i)
      {
	p[i] = input[nkin*idxthr+i]; // no integration
      }
    }    
    else
    {
#pragma unroll
      for(unsigned int i=0; i<nkin; ++i)
      {
	p[i] = input[i]; /* make copy of p a shared memory object */
      }
      for(unsigned int i=0; i<ndim; ++i)
      {   
	p[idxdim[i]] = grd[ndim*idxthr+i]; // set momenta from integration
      }

      // ensure E,p make sense for final state particles
      
      p[IP5_E] = _SQRT_(_POW_(p[IP5_X],2)+_POW_(p[IP5_Y],2)+_POW_(p[IP5_Z],2));
      p[IP4_E] = _SQRT_(_POW_(p[IP4_X],2)+_POW_(p[IP4_Y],2)+_POW_(p[IP4_Z],2));
      p[IP3_E] = _SQRT_(_POW_(p[IP3_X],2)+_POW_(p[IP3_Y],2)+_POW_(p[IP3_Z],2));
      p[IP2_E] = _SQRT_(_POW_(p[IP2_X],2)+_POW_(p[IP2_Y],2)+_POW_(p[IP2_Z],2));
    }
    
    // solve for the initial state of 2 -> 4 process

    afloat_t px_tot = p[IP2_X]+p[IP3_X]+p[IP4_X]+p[IP5_X];
    afloat_t py_tot = p[IP2_Y]+p[IP3_Y]+p[IP4_Y]+p[IP5_Y];
    afloat_t pz_tot = p[IP2_Z]+p[IP3_Z]+p[IP4_Z]+p[IP5_Z];
    afloat_t E_tot  = p[IP2_E]+p[IP3_E]+p[IP4_E]+p[IP5_E];

    afloat_t msqr_tot = _POW_(E_tot,2) - _POW_(pz_tot,2) - _POW_(py_tot,2) - _POW_(px_tot,2);

    afloat_t tmp = _SQRT_(_POW_(pz_tot,2) + msqr_tot);
    
    afloat_t e0 = (tmp + pz_tot)/2;
    afloat_t e1 = (tmp - pz_tot)/2;

    p[IP0_E] = e0; p[IP0_Z] =  e0; p[IP0_Y] = 0; p[IP0_X] = 0;
    p[IP1_E] = e1; p[IP1_Z] = -e1; p[IP1_Y] = 0; p[IP1_X] = 0;

    if(e0 > EBEAM || e1 > EBEAM)
    {
      result[idxthr] = 0;
      return;
    }         
   
    // afloat_t qsqr = _POW_(p[IP0_E]+p[IP1_E],2) - _POW_(p[IP0_Z]+p[IP1_Z],2);

    afloat_t scale = 1.0; // keep dLIPS^{-1} from blowing up ... FP precision issues

    afloat_t dLIPS = 1.0;
#pragma unroll
    for(unsigned int ip=2; ip<nexternal; ++ip)
    {
      dLIPS *= ( (2*p[ip*4]) / scale );
    }
    
    // calculate the matrix element and PDF

    result[idxthr] = ( (me(p, 0x0, 0, &_couplings)) *
		       (pdf(p[IP0_E]/EBEAM, _SQRT_(msqr_tot), gpdf_data, gpdf_nptsx, gpdf_nptsq, gpdf_lxmin, gpdf_lxmax, gpdf_lqmin, gpdf_lqmax)) *
		       (pdf(p[IP1_E]/EBEAM, _SQRT_(msqr_tot), gpdf_data, gpdf_nptsx, gpdf_nptsq, gpdf_lxmin, gpdf_lxmax, gpdf_lqmin, gpdf_lqmax)) / (dLIPS) );
  }
}
