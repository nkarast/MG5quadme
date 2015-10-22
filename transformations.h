#ifndef _transformations_h_
#define _transformations_h_

#include <math.h>

#include "common.h"

// NOTE use units of TeV (or fraction of beam energy) if not double precision
// solutions with ISR need double precision

struct three_vec { afloat_t x, y, z; };

struct four_vec  { afloat_t x, y, z, E; };

struct variables
{
  afloat_t pax, pay, paz, paE;
  afloat_t pbx, pby, pbz, pbE;
  
  afloat_t lpx, lpy, lpz, lpE;  
  afloat_t lmx, lmy, lmz, lmE;
  afloat_t vx, vy, vz, vE;  
  afloat_t vbarx, vbary, vbarz, vbarE;
  
  afloat_t isrx, isry, isrz;
  
  afloat_t qwp, qwm, qh, gy, xa, xb;
  
  afloat_t jac;
};

int hwwy(struct variables* inputs, struct variables solutions[], int maxnsol, afloat_t scale) 
{  
  const afloat_t EBEAM = 4000.0/scale;
  const afloat_t RESOL = 1.0/scale; // momentum resolution parameter 
  
  const afloat_t PREC  = 1.0e-3; // precision
  
  /* 
   *  
   *  This function gets the two missing momenta for 
   *  the following topology:
   *  
   *                        ________ missing    p1
   *                       O________ visible    p3
   *                   r1 /
   *                     /
   *  ------------------O 
   *     r3              \ 
   *                   r2 \ ________ missing    p2
   *                       O________ visible    p4
   *  
   *  
   *  Following the solution found in Madweight
   *  
   */
  
  int nsol = 0;
  
  // std::vector<afloat_t> p2y_array;
  // std::vector<afloat_t> E2_array;
  
  afloat_t p2y_array[16];
  afloat_t E2_array[16];
  
  int ind_E2 = 0;
  int ind_p2y = 0;
  
  // basic variables
  
  afloat_t p1x, p1y, p1z, E1;
  afloat_t p2x, p2y, p2z, E2;
  
  afloat_t m1_sqr = 0;
  afloat_t m2_sqr = 0;
  
  afloat_t s13_sqr = inputs->qwm; // qwm_sqr;
  afloat_t s24_sqr = inputs->qwp; // qwp_sqr;
  afloat_t s_sqr   = inputs->qh; // qh_sqr;  
  
  // with pa = Ea, pb = Eb, have
  // (1) pa - pb = sinh(y) sqrt(s)
  // (2) pa + pb = cosh(y) sqrt(s)
  //
  // which yields
  //
  // => 2 pa = [sinh(y) + cosh(y)] sqrt(s) = exp(y) sqrt(s)
  // => 2 pb = [cosh(y) - sinh(y)] sqrt(s) = exp(-y) sqrt(s)
  
  afloat_t ea = sqrt(s_sqr)*exp(inputs->gy) / 2;
  afloat_t eb = sqrt(s_sqr)*exp(-inputs->gy) / 2;
  
  // TLorentzVector pa(0,0, ea,ea);
  // TLorentzVector pb(0,0,-eb,eb);
  
  struct four_vec pa = {0, 0, ea, ea};
  struct four_vec pb = {0, 0,-eb, eb};
  
  // std::cout << "system: " << (pa + pb).M() << ", " << (pa + pb).Px() << ", " << (pa + pb).Py()
  // 	    << ", " << (pa + pb).Rapidity() << ", " << (pa + pb).E() << std::endl;
  
  // TLorentzVector recoil;
  
  if(sqrt(pow(inputs->isrx,2) + pow(inputs->isry,2)) > 0) // (isr.Pt() > RESOL)
  {
    // recoil.SetPxPyPzE(inputs->isrx, inputs->isry, 0, sqrt(pow(isr.Pt(),2) + pow((pa + pb).Pz(),2) + pow((pa + pb).M(),2)));
    
    // afloat_t boostarry[3] = { recoil.BoostVector().X(),
    // 	                         recoil.BoostVector().Y(), 0 };
    
    struct four_vec recoil = {-inputs->isrx, 
			      -inputs->isry, 
			      pa.z + pb.z,
			      sqrt(pow(inputs->isrx,2) + pow(inputs->isry,2) + pow(pa.z + pb.z,2) + 4*pa.E*pb.E)};
    
    struct three_vec boostarray = {recoil.x/recoil.E, recoil.y/recoil.E, 0};
    
    // pa.Boost(boostarry);
    // pb.Boost(boostarry);
    
    afloat_t bp, buf, bpbuf;
    
    afloat_t bsqr = boostarray.x*boostarray.x + boostarray.y*boostarray.y + boostarray.z*boostarray.z;
    afloat_t gamma = 1.0 / sqrt(1.0 - bsqr);
    
    bp = boostarray.x*pa.x + boostarray.y*pa.y + boostarray.z*pa.z;
    buf = bsqr > 0 ? (gamma - 1.0)/bsqr : 0.0;
    
    // std::cerr << buf << " " << bp << " " << gamma << " " << boostarray.x << " " << boostarray.y << std::endl;
    
    bpbuf = bp * buf;
    pa.x = pa.x + bpbuf*boostarray.x + gamma*(boostarray.x*pa.E);
    pa.y = pa.y + bpbuf*boostarray.y + gamma*(boostarray.y*pa.E);
    pa.z = pa.z + bpbuf*boostarray.z + gamma*(boostarray.z*pa.E);
    pa.E = gamma*(pa.E + bp);
    
    bp = boostarray.x*pb.x + boostarray.y*pb.y + boostarray.z*pb.z;
    buf = bsqr > 0 ? (gamma - 1.0)/bsqr : 0.0;
    bpbuf = bp * buf;
    pb.x = pb.x + bpbuf*boostarray.x + gamma*(boostarray.x*pb.E);
    pb.y = pb.y + bpbuf*boostarray.y + gamma*(boostarray.y*pb.E);
    pb.z = pb.z + bpbuf*boostarray.z + gamma*(boostarray.z*pb.E);
    pb.E = gamma*(pb.E + bp);
    
    // std::cerr << pa.E << " " << pb.E << std::endl;
  }
  
  // ea = pa.E;
  // eb = pb.E;
  
  afloat_t xa = pa.E / EBEAM; // pa.E() / hepstd::beamEnergy;
  afloat_t xb = pb.E / EBEAM; // pb.E() / hepstd::beamEnergy;
  
  // std::cout << "system: " << (pa + pb).M() << ", " << (pa + pb).Px() << ", " << (pa + pb).Py()
  // 	    << ", " << (pa + pb).Rapidity() << ", " << (pa + pb).E() << std::endl;
  
  // TLorentzVector lm = inputs.lm;
  // TLorentzVector lp = inputs.lp;
  
  afloat_t p3x = inputs->lmx; // lm.Px();
  afloat_t p3y = inputs->lmy; // lm.Py();
  afloat_t p3z = inputs->lmz; // lm.Pz();
  afloat_t E3 = inputs->lmE;  // lm.E();
  afloat_t m3_sqr = pow(E3,2)-pow(p3x,2) - pow(p3y,2) - pow(p3z,2); // pow(lm.M(),2);
  
  afloat_t p4x = inputs->lpx; // lp.Px();
  afloat_t p4y = inputs->lpy; // lp.Py();
  afloat_t p4z = inputs->lpz; // lp.Pz();
  afloat_t E4 = inputs->lpE;  // lp.E();
  afloat_t m4_sqr = pow(E4,2)-pow(p4x,2) - pow(p4y,2) - pow(p4z,2); // pow(lp.M(),2);
  
  afloat_t px_miss = -(p4x+p3x+inputs->isrx);
  afloat_t py_miss = -(p4y+p3y+inputs->isry);
  
  // std::cout << "p3x= "<<Form("%0.8f",p3x) << std::endl;
  // std::cout << "p3y= "<<Form("%0.8f",p3y) << std::endl;
  // std::cout << "p3z= "<<Form("%0.8f",p3z) << std::endl;
  // std::cout << "E3 = "<<Form("%0.8f",E3) << std::endl;
  
  // std::cout << "p4x= "<<Form("%0.8f",p4x) << std::endl;
  // std::cout << "p4y= "<<Form("%0.8f",p4y) << std::endl;
  // std::cout << "p4z= "<<Form("%0.8f",p4z) << std::endl;
  // std::cout << "E4 = "<<Form("%0.8f",E4) << std::endl;
  
  // total measured E and pZ
  
  afloat_t e_out  =  E3 +  E4;
  afloat_t pz_out = p3z + p4z;
  
  //     start from the following expressions for p1z and p2x:
  //
  //     p1z=A1 + A2 E2 + A3 p2y + A4 p2x
  //     p2x=B1 + B2 E2 + B3 p2y + A4 p1z
  
  if(fabs(p3z) < RESOL) 
  {
    p3z = p3z < 0 ? p3z - RESOL : p3z + RESOL;
  }
  if(fabs(p4x) < RESOL) 
  {
    p4x = p4x < 0 ? p4x - RESOL : p4x + RESOL;
  }
  
  afloat_t A1=(m1_sqr + m3_sqr + 2.0*E3*(pa.E + pb.E - e_out) - 2.0*p3x*px_miss - 2.0*p3y*py_miss - s13_sqr)/(2.0*p3z);
  afloat_t A2=-E3/p3z;
  afloat_t A3=p3y/p3z;
  afloat_t A4=p3x/p3z;
  
  afloat_t B1=(m2_sqr + m4_sqr - 2.0*p4z*(pa.z + pb.z - pz_out) - s24_sqr)/(2.0*p4x);
  afloat_t B2=E4/p4x;
  afloat_t B3=-p4y/p4x;
  afloat_t B4=p4z/p4x;
  
  if(fabs(1 - A4*B4) < PREC)
  {
    // std::cout << "1 - A4 B4 ~ 0" << std::endl;
    return 0;
  }
  
  afloat_t p2x_ti=(B1+B4*A1)/(1.0-A4*B4);
  afloat_t p2x_E2=(B2+B4*A2)/(1.0-A4*B4);
  afloat_t p2x_p2y=(B3+B4*A3)/(1.0-A4*B4);
  
  afloat_t p1x_ti=px_miss-p2x_ti;
  afloat_t p1x_E2=-p2x_E2;
  afloat_t p1x_p2y=-p2x_p2y;
  
  afloat_t p1y_ti=py_miss;
  afloat_t p1y_p2y=-1.0;
  afloat_t p1y_E2=0.0;
  
  afloat_t p1z_ti=(A1+A4*B1)/(1.0-A4*B4);
  afloat_t p1z_E2=(A2+A4*B2)/(1.0-A4*B4);
  afloat_t p1z_p2y=(A3+A4*B3)/(1.0-A4*B4);
  
  afloat_t E1_ti=pa.E+pb.E-e_out;
  afloat_t E1_E2=-1.0;
  afloat_t E1_p2y=0.0;
  
  afloat_t p2z_ti=pa.z+pb.z-p1z_ti-pz_out;
  afloat_t p2z_E2=-p1z_E2;
  afloat_t p2z_p2y=-p1z_p2y;
  
  // std::cerr << A1 << " " << p1x_ti << " " << p1y_ti << std::endl;
  // std::cerr << p1z_ti << " " << p2x_ti << std::endl;
  
  //
  //     mass shell conditions:     
  //    
  //     E1^2-p1x^2-p1y^2-p1z^2=m1^2   (1)
  //     E2^2-p2x^2-p2y^2-p2z^2=m2^2   (2) 
  //
  //     equivalent to
  //
  //     g11*E2^2 + g22*p2y^2 + g12*E2*p2y + g10*E2 + g20*p2y + g00 = 0  (1)
  //     h11*E2^2 + h22*p2y^2 + h12*E2*p2y + h10*E2 + h20*p2y + h00 = 0  (2)
  //
  afloat_t g11=pow(E1_E2,2)-pow(p1x_E2,2)-pow(p1y_E2,2)-pow(p1z_E2,2);
  //
  afloat_t g22=pow(E1_p2y,2)-pow(p1x_p2y,2)-pow(p1y_p2y,2)-pow(p1z_p2y,2);
  //
  afloat_t g12=2.0*(E1_E2*E1_p2y-p1x_E2*p1x_p2y-p1y_E2*p1y_p2y-p1z_E2*p1z_p2y);
  //
  afloat_t g10=2.0*(E1_E2*E1_ti-p1x_E2*p1x_ti-p1y_E2*p1y_ti-p1z_E2*p1z_ti);
  // 
  afloat_t g20=2.0*(E1_ti*E1_p2y-p1x_ti*p1x_p2y-p1y_ti*p1y_p2y-p1z_ti*p1z_p2y);
  //
  afloat_t g00=pow(E1_ti,2)-pow(p1x_ti,2)-pow(p1y_ti,2)-pow(p1z_ti,2)-m1_sqr;
  
  //
  //afloat_t h11=1.0-pow(p2x_E2,2)-pow(p2z_E2,2);
  //
  //afloat_t h22=-1.0-pow(p2x_p2y,2)-pow(p2z_p2y,2);
  //
  //afloat_t h12=-2.0*(p2x_p2y*p2x_E2+p2z_p2y*p2z_E2);
  //
  afloat_t h10=-2.0*(p2x_ti*p2x_E2+p2z_ti*p2z_E2);
  //
  afloat_t h20=-2.0*(p2x_p2y*p2x_ti+p2z_p2y*p2z_ti);
  //
  afloat_t h00=-pow(p2x_ti,2)-pow(p2z_ti,2)-m2_sqr;
  
  // std::cout << h11 << ":" << g11 << std::endl;
  // std::cout << h22 << ":" << g22 << std::endl;
  // std::cout << h12 << ":" << g12 << std::endl;
  
  afloat_t alpha = -(g00-h00)/(g20-h20); 
  afloat_t beta  = -(g10-h10)/(g20-h20);
  
  // std::cerr << alpha << " " << beta << std::endl;
  
  int solved_E2  = 0;
  int solved_p2y = 0;
  
  // std::vector< std::pair<afloat_t, afloat_t> > E2_p2y_pairs;
  
  //     note that the two quartic equations are fake:
  //
  //     h11=g11
  //     h22=g22
  //     h12=g12
  //
  //     We are left with the equation
  //
  //     (g10-h10) E2 + (g20-h20) p2y + g00-h00 = 0
  //
  //      <=>   p2y = alpha + beta E2
  
  if(fabs(g20-h20) < PREC) // case #1: coefficient of p2y is ~ 0, then E2 = (h00-g00) / (g10-h10)
  {
    // // std::cout << "case #1" << std::endl;
    if(fabs(g10-h10) > PREC)
    {
      E2_array[ind_E2++] = -(g00-h00) / (g10-h10);
      solved_E2 = 1;
      if(E2_array[ind_E2-1] > 0) // use the quartic equation (1) to solve for p2y
      {       
	afloat_t dem = g22;
	afloat_t b = (g12*E2_array[ind_E2-1]+g20)/dem;
	afloat_t c = (g11*pow(E2_array[ind_E2-1],2)+g10*E2_array[ind_E2-1]+g00)/dem;
	afloat_t rho = pow(b,2)-4.0*c;
	if(rho == 0.0)  // maximum one solution for p2y
	{  
	  p2y_array[ind_p2y++] = -b/2.0;
	  solved_p2y = 1;
	  // goto L_SOLVED_P2Y;
	}
	if(rho > 0.0) // two solutions for p2y
	{
	  p2y_array[ind_p2y++] = (-b+sqrt(rho))/2.0;
	  p2y_array[ind_p2y++] = (-b-sqrt(rho))/2.0;
	  solved_p2y = 1;	  
	  // goto L_SOLVED_P2Y; 
	}
	if(rho < 0.0)
	{
	  // std::cout << "rho < 0" << std::endl;
	  return 0;
	}
      }
      if(solved_E2==0 || solved_p2y==0)
      {
	// std::cout << "no valid solutions" << std::endl;
	return 0; // no valid solutions for E2, p2y 
      }
    }
    else
    {
      // std::cout << "g10-h10 = 0" << std::endl;
      return 0;
    }
  }
  else  // case #2: coefficient of p2y is != 0, then p2y = alpha + beta E2
  {
    // std::cout << "case #2" << std::endl;
    
    // solve: g11*E2^2 + g22*p2y^2 + g12*E2*p2y + g10*E2 + g20*p2y + g00 = 0  (1)
    // 
    // using p2y = alpha + beta E2 ... i.e., a quadratic equation
    //
    // then g11*E2^2 + g22*(alpha + beta E2)^2 + g12*E2*(alpha + beta E2) + g10*E2 + g20*(alpha + beta E2) + g00 = 0
    // => (g11 + g22*beta^2 + g12*beta)*E2^2 + (2*g22*alpha*beta + g12*alpha + g10 + g20*beta)*E2 + g22*alpha^2 + g20*alpha + g00 = 0
    
    afloat_t dem = g11+g22*pow(beta,2)+g12*beta;
    
    if(fabs(dem) < PREC) // exactly one solution
    {
      afloat_t c0 = 2.0*g22*alpha*beta+g12*alpha+g10+g20*beta;
      afloat_t c1 = g22*pow(alpha,2)+g20*alpha+g00;
      if(fabs(c0) > PREC)
      {
	E2_array[ind_E2++] = -c1/c0;
	solved_E2 = 1;
      }
      else
      {
	// std::cout << "c0 = 0" << std::endl;
	return 0;
      }
      if(E2_array[ind_E2-1] <= 0.0) 
      {
	// std::cout << "E2 < 0" << std::endl;
	return 0;
      }
    }    
    else
    {
      afloat_t b=(2.0*g22*alpha*beta+g12*alpha+g10+g20*beta)/dem;
      afloat_t c=(g22*pow(alpha,2)+g20*alpha+g00)/dem;
      afloat_t rho=pow(b,2)-4.0*c;
      
      // dem*E2^2 + b*dem*E2 + c*dem = 0
      // rho = (b/dem)^2 - 4c/dem = 1/dem^2 * (b^2 - 4*dem*c)
      //
      // ==> E2 = -b +/- sqrt(rho)
      
      if(rho == 0.0) // maximum one solution 
      {
	if(b>=0.0)
	{
	  // std::cout << "E2 < 0" << std::endl;
	  return 0; // E2 must be > 0 of course
	}
	E2_array[ind_E2++] = -b/2.0;
	solved_E2 = 1;
      }
      
      if(rho > pow(b,2)) // only one physical solution (E2 > 0)
      {
	E2_array[ind_E2++] = (-b+sqrt(rho))/2.0;
	solved_E2 = 1;
      }
      
      if(rho < 0.0) 
      {
	// std::cout << "E2 complex: " << rho << std::endl;
	return 0; // no solutions
      }
      
      if(solved_E2==0)
      {
	E2_array[ind_E2++] = (-b+sqrt(rho))/2.0;
	E2_array[ind_E2++] = (-b-sqrt(rho))/2.0;   
	solved_E2 = 1;
      }
    }
  }
  
  if(solved_p2y==0)
  {
    if((ind_E2) > 1)
    {
      for(unsigned isol = 0; isol < (ind_E2); ++isol)
      {
	p2y_array[ind_p2y++] = alpha+beta*E2_array[isol];
      }
    }
    else
    {
      p2y_array[ind_p2y++] = alpha+beta*E2_array[0];
    }
  }
  
  // std::cout << ind_p2y << " " << ind_E2 << std::endl;
  
  while((ind_E2) < (ind_p2y))
  {
    E2_array[ind_E2++] = E2_array[0];
  }
  
  afloat_t s = 4.0 * pow(EBEAM,2);
  // afloat_t xxa, xxb;
  
  nsol = 0;
  for(unsigned isol = 0; isol < min(maxnsol,ind_p2y); ++isol)
  {
    // if(isol > maxnsol)
    // {
    //   // std::cout << "maxnsol" << std::endl;
    //   break;
    // }
    
    E2  = E2_array[isol];
    p2y = p2y_array[isol];
    
    // if((E2!=E2) || (p2y!=p2y)) // (std::isnan(E2) || std::isnan(p2y)) // NaN check ...
    // {
    //   // std::cout << "E2, p2y = " << E2 << ", " << p2y << std::endl;
    //   continue;
    // }
    
    // E1  = E1_ti  + E1_E2*E2  + E1_p2y*p2y;
    // p1z = p1z_ti + p1z_E2*E2 + p1z_p2y*p2y;
    // p2z = p2z_ti + p2z_E2*E2 + p2z_p2y*p2y;
    
    p1x=p1x_ti + p1x_E2*E2 + p1x_p2y*p2y;
    p1y=p1y_ti + p1y_E2*E2 + p1y_p2y*p2y;
    p1z=p1z_ti + p1z_E2*E2 + p1z_p2y*p2y;
    
    p2x=p2x_ti + p2x_E2*E2 + p2x_p2y*p2y;
    p2z=p2z_ti + p2z_E2*E2 + p2z_p2y*p2y;
    
    // check energy and momentum conservation for trial solutions ... 
    
    // if(fabs(pa.E+pb.E-e_out-E1-E2)    > RESOL ||
    //    fabs(pa.z+pb.z-pz_out-p1z-p2z) > RESOL)
    // {
    //   // std::cout << "E,pZ not conserved" << std::endl;
    //   continue;
    // }
    //
    // if(E1 < 0.0)
    // {
    //   // std::cout << "E1 < 0" << std::endl;
    //   continue;
    // }
    
    // if(fabs(pow(E1,2) - (pow(p1x,2)+pow(p1y,2)+pow(p1z,2)) /* norm(p1x, p1y, p1z),2) */ - m1_sqr) > RESOL) // check for massive neutrinos
    // {
    //   // std::cout << "m1 != 0 (" << fabs(pow(E1,2) - (pow(p1x,2)+pow(p1y,2)+pow(p1z,2)) - m1_sqr) /* pow(norm(p1x, p1y, p1z),2)) */  << ")" << std::endl;
    //   // continue;
    // }
    
    // if(fabs(pow(E2,2) - (pow(p2x,2)+pow(p2y,2)+pow(p2z,2)) /* pow(norm(p2x, p2y, p2z),2) */ - m2_sqr) > RESOL) // check for massive neutrinos
    // {
    //   // std::cout << "m2 != 0 (" << fabs(pow(E2,2) - (pow(p2x,2)+pow(p2y,2)+pow(p2z,2)) - m2_sqr) /* pow(norm(p2x, p2y, p2z),2)) */ << ")" << std::endl;
    //   // continue;
    // }
    
    afloat_t invjac = 16.0*(E4*(p1z*p2y*p3x - p1y*p2z*p3x - p1z*p2x*p3y + 
				p1x*p2z*p3y + p1y*p2x*p3z - p1x*p2y*p3z) + 
			    E2*p1z*p3y*p4x - E1*p2z*p3y*p4x - E2*p1y*p3z*p4x + 
			    E1*p2y*p3z*p4x - E2*p1z*p3x*p4y + E1*p2z*p3x*p4y + 
			    E2*p1x*p3z*p4y - E1*p2x*p3z*p4y + 
			    (E2*p1y*p3x - E1*p2y*p3x - E2*p1x*p3y + E1*p2x*p3y)*p4z + 
			    E3*(-(p1z*p2y*p4x) + p1y*p2z*p4x + p1z*p2x*p4y - 
				p1x*p2z*p4y - p1y*p2x*p4z + p1x*p2y*p4z));
    
    // if(fabs(invjac) < (PREC / 1.0e3)) 
    // {
    //   // std::cout << "1/J ~ 0" << std::endl;
    //   continue;
    // }
    
    // create parton container ...
    
    solutions[nsol].lpx = inputs->lpx;
    solutions[nsol].lpy = inputs->lpy;
    solutions[nsol].lpz = inputs->lpz;
    solutions[nsol].lpE = inputs->lpE;
    solutions[nsol].lmx = inputs->lmx;
    solutions[nsol].lmy = inputs->lmy;
    solutions[nsol].lmz = inputs->lmz;
    solutions[nsol].lmE = inputs->lmE;
    solutions[nsol].pax = pa.x;
    solutions[nsol].pay = pa.y;
    solutions[nsol].paz = pa.z;
    solutions[nsol].paE = pa.E;
    solutions[nsol].pbx = pb.x;
    solutions[nsol].pby = pb.y;
    solutions[nsol].pbz = pb.z;
    solutions[nsol].pbE = pb.E;
    solutions[nsol].vx = p2x;
    solutions[nsol].vy = p2y;
    solutions[nsol].vz = p2z;
    solutions[nsol].vbarx = p1x;
    solutions[nsol].vbary = p1y;
    solutions[nsol].vbarz = p1z;
    solutions[nsol].isrx = inputs->isrx;
    solutions[nsol].isry = inputs->isry;
    solutions[nsol].jac = 1.0 / fabs(invjac) / (2 * xa * xb * s) * 1.0 / s;
    
    // if( sqrt(pow(solutions[nsol].lpx + solutions[nsol].lmx + solutions[nsol].vx + solutions[nsol].vbarx + solutions[nsol].isrx,2) + 
    // 	     pow(solutions[nsol].lpy + solutions[nsol].lmy + solutions[nsol].vy + solutions[nsol].vbary + solutions[nsol].isry,2)) > RESOL )
    // {
    //   // std::cout << "Warning: transverse momentum balance failed" << std::endl;
    //   continue;
    // }
    
    nsol += 1;
  }
  
  return nsol; // sol.size();
}

#endif
