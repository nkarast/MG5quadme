#include "common.h"

#ifndef __CINT__
#include "higgs.cc"
#include "higgscouplings.cc"
#include "higgshelas.cc"
#endif

#include "pdf.h"

#include "TLorentzVector.h"

#include <iomanip> 
#include <iostream>


void quadme(
 	     afloat_t* result, afloat_t* input, unsigned int len, afloat_t* grd, int* idxdim, int ndim,  /* phase space data */
	     afloat_t* gpdf_data, int gpdf_nptsx, int gpdf_nptsq, afloat_t gpdf_lxmin, afloat_t gpdf_lxmax, afloat_t gpdf_lqmin, afloat_t gpdf_lqmax, /* gluon PDF data table */
	     afloat_t ca, afloat_t kSM, afloat_t kHWW, afloat_t kAWW, afloat_t kHdw, afloat_t lambda, /* Higgs characterization model parameters */
	     int sel_ind=-1 
	    )
{
  struct couplings default_couplings;
  init(&default_couplings, ca, kSM, kHWW, kAWW,kHdw, lambda);
  
  const uint nkin = 4*nexternal;
  
  const uint IP0_E = 4*0; const uint IP0_X = 4*0+1; const uint IP0_Y = 4*0+2; const uint IP0_Z = 4*0+3;
  const uint IP1_E = 4*1; const uint IP1_X = 4*1+1; const uint IP1_Y = 4*1+2; const uint IP1_Z = 4*1+3;
  const uint IP2_E = 4*2; const uint IP2_X = 4*2+1; const uint IP2_Y = 4*2+2; const uint IP2_Z = 4*2+3;
  const uint IP3_E = 4*3; const uint IP3_X = 4*3+1; const uint IP3_Y = 4*3+2; const uint IP3_Z = 4*3+3;
  const uint IP4_E = 4*4; const uint IP4_X = 4*4+1; const uint IP4_Y = 4*4+2; const uint IP4_Z = 4*4+3;
  const uint IP5_E = 4*5; const uint IP5_X = 4*5+1; const uint IP5_Y = 4*5+2; const uint IP5_Z = 4*5+3;
  
  // std::cout << "GC_61: " <<  default_couplings.GC_61.re << "," << default_couplings.GC_61.im << std::endl
  // 	    << "GC_89: " <<  default_couplings.GC_89.re << "," << default_couplings.GC_89.im << std::endl
  // 	    << "GC_35: " <<  default_couplings.GC_35.re << "," << default_couplings.GC_35.im << std::endl
  // 	    << "GC_37: " <<  default_couplings.GC_37.re << "," << default_couplings.GC_37.im << std::endl
  // 	    << "GC_56: " <<  default_couplings.GC_56.re << "," << default_couplings.GC_56.im << std::endl
  // 	    << "GC_53: " <<  default_couplings.GC_53.re << "," << default_couplings.GC_53.im << std::endl
  // 	    << "GC_18: " <<  default_couplings.GC_18.re << "," << default_couplings.GC_18.im << std::endl;
  
  afloat_t p[nkin];
  for(unsigned int ind=0; ind<(int)len; ++ind)
  {
    if(sel_ind >= 0 && ind != (unsigned)sel_ind)
      continue;

    if(ndim==0)
    {
      for(unsigned int i=0; i<nkin; ++i)
      {
	p[i] = input[nkin*ind+i];
      }
    }
    else
    {    
      for(unsigned int i=0; i<nkin; ++i)
      {
	p[i] = input[i];
      }
      for(unsigned int i=0; i<ndim; ++i)
      {   
	p[idxdim[i]] = grd[ndim*ind+i]; // set momenta from integration
      }

      // ensure E,p make sense for final state particles
      
      p[IP5_E] = _SQRT_(_POW_(p[IP5_X],2)+_POW_(p[IP5_Y],2)+_POW_(p[IP5_Z],2));
      p[IP4_E] = _SQRT_(_POW_(p[IP4_X],2)+_POW_(p[IP4_Y],2)+_POW_(p[IP4_Z],2));
      p[IP3_E] = _SQRT_(_POW_(p[IP3_X],2)+_POW_(p[IP3_Y],2)+_POW_(p[IP3_Z],2));
      p[IP2_E] = _SQRT_(_POW_(p[IP2_X],2)+_POW_(p[IP2_Y],2)+_POW_(p[IP2_Z],2));
    }
    
    // solve for the initial state of 2 -> 4 process
    
  /*  afloat_t px_tot = p[IP2_X]+p[IP3_X]+p[IP4_X]+p[IP5_X];
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
      result[ind] = 0;
      continue;
    }      
    */
      
    // afloat_t qsqr = _POW_(p[IP0_E]+p[IP1_E],2) - _POW_(p[IP0_Z]+p[IP1_Z],2);
    
    // if(ind % (len/15) == 0)
    //   std::cout << "sqrt(q**2) = " << _SQRT_(qsqr)
    // 		<< " x1 = " << p[IP0_E]/EBEAM 
    // 		<< " x2 = " << p[IP1_E]/EBEAM << std::endl;

    afloat_t offset = 1.0;

    afloat_t dLIPS = 1.0;

#pragma unroll
    for(unsigned int ip=2; ip<nexternal; ++ip)
    {
      dLIPS *= ( (2*p[ip*4]) / offset );
    }
    
    // calculate the matrix element

    // result[ind] = ( (me(p, 0x0, 0, &default_couplings)) *
    // 		    (pdf(p[IP0_E]/EBEAM, _SQRT_(msqr_tot), gpdf_data, gpdf_nptsx, gpdf_nptsq, gpdf_lxmin, gpdf_lxmax, gpdf_lqmin, gpdf_lqmax)) *
    // 		    (pdf(p[IP1_E]/EBEAM, _SQRT_(msqr_tot), gpdf_data, gpdf_nptsx, gpdf_nptsq, gpdf_lxmin, gpdf_lxmax, gpdf_lqmin, gpdf_lqmax)) / (dLIPS) );

    result[ind] = ( (me(p, 0x0, 0, &default_couplings)) ); //  / (dLIPS) );

    // if(ind % (len/15) == 0)
    //   std::cout << "ME = " << result[ind] << std::endl;
  }
}
