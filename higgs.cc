#ifndef _higgs_cc_
#define _higgs_cc_

#include "common.h"
#include "dcmplx.h"
#include "higgscouplings.h"
#include "higgshelas.h"

#ifndef _CUDA_READY_
#include <iostream>
#endif

//
// tlv is a 6x4 unrolled array
//
// { {E,px,py,pz}_{0}, ..., {E,px,py,pz}_{5} }
//

const int ninitial	= 2; 
const int nexternal	= 6; 

_CUDA_HOST_ _CUDA_DEVICE_ afloat_t me(afloat_t tlv[], int goodhel[], int skiphel, couplings* pars)
{
  const int nprocesses	= 1; 
  const int ncomb	= 64; 
  const int nwavefuncs	= 10; 
  const int namplitudes = 2; 
  const int ngraphs     = 2; 
  const int ncolor      = 1; 
  
  // struct couplings default_pars;
  // init(&default_pars, 1, 1, 0, 0);

  // if(pars==0x0)
  // {
  //   return 0; // pars = &default_pars;
  // }

  afloat_t p[nexternal][4] = { {0,0,0,0},
			       {0,0,0,0},  
			       {0,0,0,0},
			       {0,0,0,0},
			       {0,0,0,0},
			       {0,0,0,0} };

#pragma unroll
  for(int i=0; i < nexternal; ++i)
  {
    p[i][0] = tlv[i*4+0]; // E 
    p[i][1] = tlv[i*4+1]; // px
    p[i][2] = tlv[i*4+2]; // py
    p[i][3] = tlv[i*4+3]; // pz
  }
  
  dcmplx w[nwavefuncs][18]; 
  dcmplx amp[namplitudes]; 
  
  const int helicities[ncomb][nexternal] = {
    {-1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1,  1},
    {-1, -1, -1, -1,  1, -1},
    {-1, -1, -1, -1,  1,  1},						 
    {-1, -1, -1,  1, -1, -1},
    {-1, -1, -1,  1, -1,  1},
    {-1, -1, -1,  1,  1, -1},
    {-1, -1, -1,  1,  1,  1},
    {-1, -1,  1, -1, -1, -1},
    {-1, -1,  1, -1, -1,  1},
    {-1, -1,  1, -1,  1, -1},
    {-1, -1,  1, -1,  1,  1},
    {-1, -1,  1,  1, -1, -1},
    {-1, -1,  1,  1, -1,  1},
    {-1, -1,  1,  1,  1, -1},
    {-1, -1,  1,  1,  1,  1},
    {-1,  1, -1, -1, -1, -1},
    {-1,  1, -1, -1, -1,  1},
    {-1,  1, -1, -1,  1, -1},
    {-1,  1, -1, -1,  1,  1},
    {-1,  1, -1,  1, -1, -1},
    {-1,  1, -1,  1, -1,  1},
    {-1,  1, -1,  1,  1, -1},
    {-1,  1, -1,  1,  1,  1},
    {-1,  1,  1, -1, -1, -1},
    {-1,  1,  1, -1, -1,  1},
    {-1,  1,  1, -1,  1, -1},
    {-1,  1,  1, -1,  1,  1},
    {-1,  1,  1,  1, -1, -1},
    {-1,  1,  1,  1, -1,  1},
    {-1,  1,  1,  1,  1, -1},
    {-1,  1,  1,  1,  1,  1},
    {1, -1, -1, -1, -1, -1},
    {1, -1, -1, -1, -1,  1},
    {1, -1, -1, -1,  1, -1},
    {1, -1, -1, -1,  1,  1},
    {1, -1, -1,  1, -1, -1},
    {1, -1, -1,  1, -1,  1},
    {1, -1, -1,  1,  1, -1},
    {1, -1, -1,  1,  1,  1},
    {1, -1,  1, -1, -1, -1},
    {1, -1,  1, -1, -1,  1},
    {1, -1,  1, -1,  1, -1},
    {1, -1,  1, -1,  1,  1},
    {1, -1,  1,  1, -1, -1},
    {1, -1,  1,  1, -1,  1},
    {1, -1,  1,  1,  1, -1},
    {1, -1,  1,  1,  1,  1},
    {1,  1, -1, -1, -1, -1},
    {1,  1, -1, -1, -1,  1},
    {1,  1, -1, -1,  1, -1},
    {1,  1, -1, -1,  1,  1},
    {1,  1, -1,  1, -1, -1},
    {1,  1, -1,  1, -1,  1},
    {1,  1, -1,  1,  1, -1},
    {1,  1, -1,  1,  1,  1},
    {1,  1,  1, -1, -1, -1},
    {1,  1,  1, -1, -1,  1},
    {1,  1,  1, -1,  1, -1},
    {1,  1,  1, -1,  1,  1},
    {1,  1,  1,  1, -1, -1},
    {1,  1,  1,  1, -1,  1},
    {1,  1,  1,  1,  1, -1},
    {1,  1,  1,  1,  1,  1} };
  
  afloat_t jampsqr[nprocesses][ncolor] = { { 0.0 } };
  dcmplx ztemp; 
  dcmplx jamp[ncolor]; 
  
  const afloat_t denom[ncolor] = {1}; 
  const afloat_t cf[ncolor][ncolor] = {{2}}; 
  const int denominators[nprocesses] = {256}; 
  
  afloat_t matrix_element=0;

#ifndef _CUDA_READY_
  // std::cout << "GC_61: " <<  pars->GC_61.re << "," << pars->GC_61.im << std::endl
  // 	    << "GC_89: " <<  pars->GC_89.re << "," << pars->GC_89.im << std::endl
  // 	    << "GC_35: " <<  pars->GC_35.re << "," << pars->GC_35.im << std::endl
  // 	    << "GC_37: " <<  pars->GC_37.re << "," << pars->GC_37.im << std::endl
  // 	    << "GC_56: " <<  pars->GC_56.re << "," << pars->GC_56.im << std::endl
  // 	    << "GC_53: " <<  pars->GC_53.re << "," << pars->GC_53.im << std::endl
  // 	    << "GC_18: " <<  pars->GC_18.re << "," << pars->GC_18.im << std::endl;
#endif
  
#pragma unroll
  for(int ihel = 0; ihel < ncomb; ihel++ )
  {
    if(skiphel != 0)
      if(!goodhel[ihel]) continue;
    
    // calculate all wavefunctions
    vxxxxx(p[0], 0.0, helicities[ihel][0], -1, w[0]); 
    vxxxxx(p[1], 0.0, helicities[ihel][1], -1, w[1]); 
    ixxxxx(p[2], 0.0, helicities[ihel][2], -1, w[2]); 
    oxxxxx(p[3], 0.0, helicities[ihel][3], +1, w[3]); 
    oxxxxx(p[4], 0.0, helicities[ihel][4], +1, w[4]); 
    ixxxxx(p[5], 0.0, helicities[ihel][5], -1, w[5]); 
    ffv3_3(w[2], w[3], pars->GC_61, pars->MW, pars->WW, w[6]);     
    ffv3_3(w[5], w[4], pars->GC_61, pars->MW, pars->WW, w[7]); 
    vvs4_3(w[7], w[6], pars->GC_89, pars->MX0, pars->WX0, w[8]); 
    vvs1_6_7_3(w[7], w[6], pars->GC_56, pars->GC_37, pars->GC_35, pars->MX0, pars->WX0, w[9]);
    
    // calculate all amplitudes
    vvs2_6_0(w[0], w[1], w[8], pars->GC_53, pars->GC_18, &(amp[0])); 
    vvs2_6_0(w[0], w[1], w[9], pars->GC_53, pars->GC_18, &(amp[1])); 
    
    jamp[0] = +2. * (-amp[0] - amp[1]); 

    // return fabsc(jamp[0]*(conj(jamp[0])));
    
    afloat_t matrix = 0; 
    for(int i = 0; i < ncolor; i++ )
    {
      ztemp = mkdcmplx(0.,0.); 
      for(int j = 0; j < ncolor; j++ )
	ztemp = ztemp + cf[i][j] * jamp[j]; 
      matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
    }
    
    if(matrix > 0 && skiphel != 0)
    {
      goodhel[ihel]=1;
    }
    
    // store the leading color flows for choice of color
    for(int i = 0; i < ncolor; i++ )
      jampsqr[0][i] += real(jamp[i] * conj(jamp[i])); 
    
    matrix_element += matrix / denominators[0];
  }
  
  return matrix_element;
}

#endif
