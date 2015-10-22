
_CUDA_HOST_ _CUDA_DEVICE_ afloat_t gg_fastpdf(afloat_t m, afloat_t Y, afloat_t sqrtsVal)
{
  //
  // copied shamelessly from JHU authors ...
  //
  
  afloat_t s0 = sqrtsVal*sqrtsVal;
  afloat_t s = m*m;
  afloat_t Q = m;
  afloat_t xa = _EXP_(Y)*_SQRT_(s/s0);
  afloat_t xb = _EXP_(-Y)*_SQRT_(s/s0);
  
  afloat_t weightg = (afloat_t)1.0;
  
  const afloat_t g0par0 = (afloat_t)0.2282;  const afloat_t g0par1 = (afloat_t)-0.0002252; const afloat_t g0par2 = (afloat_t)1.383e-07;
  const afloat_t g1par0 = (afloat_t)0.01968; const afloat_t g1par1 = (afloat_t)-0.0002993; const afloat_t g1par2 = (afloat_t)1.986e-07;
  const afloat_t g2par0 = (afloat_t)3.624;   const afloat_t g2par1 = (afloat_t)-0.003164;  const afloat_t g2par2 = (afloat_t)1.941e-06;
  const afloat_t g3par0 = (afloat_t)-0.578;  const afloat_t g3par1 = (afloat_t)-0.0003;    const afloat_t g3par2 = (afloat_t)1.828e-07;
  const afloat_t g4par0 = (afloat_t)-7.515;  const afloat_t g4par1 = (afloat_t)-0.001355;  const afloat_t g4par2 = (afloat_t)8.199e-07;
  
  afloat_t gluon0 = g0par0 + g0par1*Q + g0par2*Q*Q;
  afloat_t gluon1 = g1par0 + g1par1*Q + g1par2*Q*Q;
  afloat_t gluon2 = g2par0 + g2par1*Q + g2par2*Q*Q;
  afloat_t gluon3 = g3par0 + g3par1*Q + g3par2*Q*Q;
  afloat_t gluon4 = g4par0 + g4par1*Q + g4par2*Q*Q;
  
  afloat_t Funcga = (gluon0+gluon1*xa+gluon2*_POW_(xa,2))*_POW_((1-xa),4)*_POW_(xa,gluon3)*_EXP_(1.0+gluon4*xa);
  afloat_t Funcgb = (gluon0+gluon1*xb+gluon2*_POW_(xb,2))*_POW_((1-xb),4)*_POW_(xb,gluon3)*_EXP_(1.0+gluon4*xb);
  afloat_t FuncABg = Funcga*Funcgb/xa/xb;
  
  afloat_t _pdf = 2*m*((FuncABg)*weightg);
  
  if((m <= (afloat_t)600. && fabs(Y) > 20*_POW_(m,(afloat_t)-0.32)) || (m > (afloat_t)600. && fabs(Y) > 21*_POW_(m,(afloat_t)-0.34)))
  {
    afloat_t xa0 = _SQRT_(s/s0); // @ Y=0 xa=xb
    Funcga = (gluon0+gluon1*xa0+gluon2*_POW_(xa0,2))*_POW_((1-xa0),4)*_POW_(xa0,gluon3)*_EXP_((afloat_t)1.0+gluon4*xa0);
    FuncABg = Funcga*Funcga/xa0/xa0;
    afloat_t _tmp = 2*m*((FuncABg)*weightg);
    _pdf = ((afloat_t)1.e-5)*_tmp;
  }
  
  if (_pdf <= (afloat_t)0.) _pdf = (afloat_t)0.00001; 
  _pdf /= (afloat_t)10.;
  
  return _pdf*2;                                              
}

_CUDA_HOST_ _CUDA_DEVICE_ afloat_t qqbar_fastpdf(afloat_t m, afloat_t Y, afloat_t sqrts)
{
  
  //
  // copied shamelessly from JHU authors ...
  //

  afloat_t YVal = Y;
  afloat_t sqrtsVal = sqrts;
  afloat_t mVal = m;
  
  afloat_t weightu = 0.5;
  afloat_t weightd = 0.5;
  afloat_t weightc = 1.0;
  afloat_t weights = 1.0;
  afloat_t weightb = 1.0;
  
  afloat_t s0 = sqrtsVal*sqrtsVal;
  afloat_t s = mVal*mVal;
  afloat_t Q = mVal;
  afloat_t xa = _EXP_(YVal)*_SQRT_(s/s0);
  afloat_t xb = _EXP_(-YVal)*_SQRT_(s/s0);
  
  // PDF parameters
  // up params
  afloat_t u0par0 = 0.03134;	afloat_t u0par1 =-2.068e-05;	afloat_t u0par2 = 1.283e-08; 
  afloat_t u1par0 = 0.9;	afloat_t u1par1 =-0.0004307;	afloat_t u1par2 = 2.458e-07;
  afloat_t u2par0 =-0.1369;	afloat_t u2par1 = 0.003423;	afloat_t u2par2 =-2.155e-06;
  afloat_t u3par0 =-0.4013;	afloat_t u3par1 =-0.0002574;	afloat_t u3par2 = 1.561e-07;
  afloat_t u4par0 = 0.5782;	afloat_t u4par1 =-0.004728;	afloat_t u4par2 = 2.906e-06;
  afloat_t ubar0par0 = 0.02856;	afloat_t ubar0par1 =-2.112e-05;	afloat_t ubar0par2 = 1.272e-08;
  afloat_t ubar1par0 =-0.06822;	afloat_t ubar1par1 = 3.172e-05;	afloat_t ubar1par2 =-2.008e-08;
  afloat_t ubar2par0 = 0.1967;	afloat_t ubar2par1 =-0.000118;	afloat_t ubar2par2 = 6.871e-08;
  afloat_t ubar3par0 =-0.2251;	afloat_t ubar3par1 = 0.0001295;	afloat_t ubar3par2 =-7.181e-08;
  afloat_t ubar4par0 =-0.4068;	afloat_t ubar4par1 =-0.0002956;	afloat_t ubar4par2 = 1.783e-07;
  afloat_t ubar5par0 =-2.251;	afloat_t ubar5par1 =-0.0001699;	afloat_t ubar5par2 = 1.492e-07;
  // down params
  afloat_t d0par0 = 0.03278;	afloat_t d0par1 =-2.915e-05;	afloat_t d0par2 = 1.809e-08;
  afloat_t d1par0 = 0.479;	afloat_t d1par1 =-0.0002559;	afloat_t d1par2 = 1.557e-07;
  afloat_t d2par0 =-0.5972;	afloat_t d2par1 = 0.0003118;	afloat_t d2par2 =-1.905e-07;
  afloat_t d3par0 =-0.3892;	afloat_t d3par1 =-0.000317;	afloat_t d3par2 = 1.944e-07;
  afloat_t d4par0 = 0.5007;	afloat_t d4par1 =-0.001665;	afloat_t d4par2 = 9.895e-07;
  afloat_t dbar0par0 = 0.02328;	afloat_t dbar0par1 =-1.367e-05;	afloat_t dbar0par2 = 8.246e-09;
  afloat_t dbar1par0 = 0.09422;	afloat_t dbar1par1 =-0.0001019;	afloat_t dbar1par2 = 6.375e-08;
  afloat_t dbar2par0 =-0.5296;	afloat_t dbar2par1 = 0.000466;	afloat_t dbar2par2 =-2.896e-07;
  afloat_t dbar3par0 = 0.5354;	afloat_t dbar3par1 =-0.0004404;	afloat_t dbar3par2 = 2.728e-07;
  afloat_t dbar4par0 =-0.4386;	afloat_t dbar4par1 =-0.0002605;	afloat_t dbar4par2 = 1.582e-07;
  afloat_t dbar5par0 =-1.289;	afloat_t dbar5par1 =-0.001618;	afloat_t dbar5par2 = 9.601e-07;
  // charm, strange, bottom params
  afloat_t c0par0 = 0.01829;	afloat_t c0par1 =-6.93e-06;	afloat_t c0par2 = 3.796e-09;
  afloat_t c1par0 = 0.03081;	afloat_t c1par1 = 4.325e-05;	afloat_t c1par2 =-3.95e-08;
  afloat_t c2par0 = 0.5398;	afloat_t c2par1 =-4.284e-05;	afloat_t c2par2 =-1.362e-08;
  afloat_t c3par0 =-0.5986;	afloat_t c3par1 = 0.002565;	afloat_t c3par2 =-1.937e-06;
  afloat_t c4par0 =-0.4534;	afloat_t c4par1 =-0.0002329;	afloat_t c4par2 = 1.343e-07;
  afloat_t c5par0 =-8.657;	afloat_t c5par1 =-0.005157;	afloat_t c5par2 = 3.68e-06;
  afloat_t s0par0 = 0.01312;	afloat_t s0par1 =-3.743e-06;	afloat_t s0par2 = 2.076e-09;
  afloat_t s1par0 =-0.001416;	afloat_t s1par1 =-7.649e-06;	afloat_t s1par2 = 4.757e-09;
  afloat_t s2par0 = 0.2864;	afloat_t s2par1 =-6.693e-05;	afloat_t s2par2 = 3.566e-08;
  afloat_t s3par0 =-0.4857;	afloat_t s3par1 =-0.000253;	afloat_t s3par2 = 1.541e-07;
  afloat_t s4par0 =-10.33;	afloat_t s4par1 =-0.001601;	afloat_t s4par2 = 9.718e-07;
  afloat_t b0par0 = 0.005934;	afloat_t b0par1 = 2.516e-06;	afloat_t b0par2 =-1.828e-09;
  afloat_t b1par0 =-0.003063;	afloat_t b1par1 =-6.761e-06;	afloat_t b1par2 = 4.298e-09;
  afloat_t b2par0 = 0.1174;	afloat_t b2par1 = 3.752e-05;	afloat_t b2par2 =-2.863e-08;
  afloat_t b3par0 =-0.5549;	afloat_t b3par1 =-0.0002205;	afloat_t b3par2 = 1.334e-07;
  afloat_t b4par0 =-10.18;	afloat_t b4par1 =-0.001136;	afloat_t b4par2 = 6.931e-07;
  
  
  // PDF definition
  afloat_t up0 = u0par0 + u0par1*Q + u0par2*Q*Q;
  afloat_t up1 = u1par0 + u1par1*Q + u1par2*Q*Q;
  afloat_t up2 = u2par0 + u2par1*Q + u2par2*Q*Q;
  afloat_t up3 = u3par0 + u3par1*Q + u3par2*Q*Q;
  afloat_t up4 = u4par0 + u4par1*Q + u4par2*Q*Q;
  afloat_t antiup0 = ubar0par0 + ubar0par1*Q + ubar0par2*Q*Q;
  afloat_t antiup1 = ubar1par0 + ubar1par1*Q + ubar1par2*Q*Q;
  afloat_t antiup2 = ubar2par0 + ubar2par1*Q + ubar2par2*Q*Q;
  afloat_t antiup3 = ubar3par0 + ubar3par1*Q + ubar3par2*Q*Q;
  afloat_t antiup4 = ubar4par0 + ubar4par1*Q + ubar4par2*Q*Q;
  afloat_t antiup5 = ubar5par0 + ubar5par1*Q + ubar5par2*Q*Q;
  afloat_t down0 = d0par0 + d0par1*Q + d0par2*Q*Q;
  afloat_t down1 = d1par0 + d1par1*Q + d1par2*Q*Q;
  afloat_t down2 = d2par0 + d2par1*Q + d2par2*Q*Q;
  afloat_t down3 = d3par0 + d3par1*Q + d3par2*Q*Q;
  afloat_t down4 = d4par0 + d4par1*Q + d4par2*Q*Q;
  afloat_t antidown0 = dbar0par0 + dbar0par1*Q + dbar0par2*Q*Q;
  afloat_t antidown1 = dbar1par0 + dbar1par1*Q + dbar1par2*Q*Q;
  afloat_t antidown2 = dbar2par0 + dbar2par1*Q + dbar2par2*Q*Q;
  afloat_t antidown3 = dbar3par0 + dbar3par1*Q + dbar3par2*Q*Q;
  afloat_t antidown4 = dbar4par0 + dbar4par1*Q + dbar4par2*Q*Q;
  afloat_t antidown5 = dbar5par0 + dbar5par1*Q + dbar5par2*Q*Q;
  afloat_t charm0 = c0par0 + c0par1*Q + c0par2*Q*Q;
  afloat_t charm1 = c1par0 + c1par1*Q + c1par2*Q*Q;
  afloat_t charm2 = c2par0 + c2par1*Q + c2par2*Q*Q;
  afloat_t charm3 = c3par0 + c3par1*Q + c3par2*Q*Q;
  afloat_t charm4 = c4par0 + c4par1*Q + c4par2*Q*Q;
  afloat_t charm5 = c5par0 + c5par1*Q + c5par2*Q*Q;
  afloat_t strange0 = s0par0 + s0par1*Q + s0par2*Q*Q;
  afloat_t strange1 = s1par0 + s1par1*Q + s1par2*Q*Q;
  afloat_t strange2 = s2par0 + s2par1*Q + s2par2*Q*Q;
  afloat_t strange3 = s3par0 + s3par1*Q + s3par2*Q*Q;
  afloat_t strange4 = s4par0 + s4par1*Q + s4par2*Q*Q;
  afloat_t bottom0 = b0par0 + b0par1*Q + b0par2*Q*Q;
  afloat_t bottom1 = b1par0 + b1par1*Q + b1par2*Q*Q;
  afloat_t bottom2 = b2par0 + b2par1*Q + b2par2*Q*Q;
  afloat_t bottom3 = b3par0 + b3par1*Q + b3par2*Q*Q;
  afloat_t bottom4 = b4par0 + b4par1*Q + b4par2*Q*Q;
  
  afloat_t FuncAu1 = (up0+up1*xa+up2*_POW_(xa,2))*_POW_((1-xa),4)*_POW_(xa,up3)*_EXP_(1.0+up4*xa);
  afloat_t FuncBu1 = (antiup0+antiup1*xb+antiup2*_POW_(xb,2)+antiup3*_POW_(xb,3))*_POW_((1-xb),4)*_POW_(xb,antiup4)*_EXP_(1.0+antiup5*xb);
  afloat_t FuncAu2 = (up0+up1*xb+up2*_POW_(xb,2))*_POW_((1-xb),4)*_POW_(xb,up3)*_EXP_(1.0+up4*xb);
  afloat_t FuncBu2 = (antiup0+antiup1*xa+antiup2*_POW_(xa,2)+antiup3*_POW_(xa,3))*_POW_((1-xa),4)*_POW_(xa,antiup4)*_EXP_(1.0+antiup5*xa);
  afloat_t FuncABu = FuncAu1/xa*FuncBu1/xb+FuncAu2/xa*FuncBu2/xb;
  
  
  afloat_t FuncAd1 = (down0+down1*xa+down2*_POW_(xa,2))*_POW_((1-xa),4)*_POW_(xa,down3)*_EXP_(1.0+down4*xa);
  afloat_t FuncBd1 = (antidown0+antidown1*xb+antidown2*_POW_(xb,2)+antidown3*_POW_(xb,3))*_POW_((1-xb),4)*_POW_(xb,antidown4)*_EXP_(1.0+antidown5*xb);
  afloat_t FuncAd2 = (down0+down1*xb+down2*_POW_(xb,2))*_POW_((1-xb),4)*_POW_(xb,down3)*_EXP_(1.0+down4*xb);
  afloat_t FuncBd2 = (antidown0+antidown1*xa+antidown2*_POW_(xa,2)+antidown3*_POW_(xa,3))*_POW_((1-xa),4)*_POW_(xa,antidown4)*_EXP_(1.0+antidown5*xa);
  afloat_t FuncABd = FuncAd1/xa*FuncBd1/xb+FuncAd2/xa*FuncBd2/xb;
  
  afloat_t Funcca = (charm0+charm1*xa+charm2*_POW_(xa,2)+charm3*_POW_(xa,3))*_POW_((1-xa),4)*_POW_(xa,charm4)*_EXP_(1.0+charm5*xa);
  afloat_t Funccb = (charm0+charm1*xb+charm2*_POW_(xb,2)+charm3*_POW_(xb,3))*_POW_((1-xb),4)*_POW_(xb,charm4)*_EXP_(1.0+charm5*xb);
  afloat_t Funcsa = (strange0+strange1*xa+strange2*_POW_(xa,2))*_POW_((1-xa),4)*_POW_(xa,strange3)*_EXP_(1.0+strange4*xa);
  afloat_t Funcsb = (strange0+strange1*xb+strange2*_POW_(xb,2))*_POW_((1-xb),4)*_POW_(xb,strange3)*_EXP_(1.0+strange4*xb);
  afloat_t Funcba = (bottom0+bottom1*xa+bottom2*_POW_(xa,2))*_POW_((1-xa),4)*_POW_(xa,bottom3)*_EXP_(1.0+bottom4*xa);
  afloat_t Funcbb = (bottom0+bottom1*xb+bottom2*_POW_(xb,2))*_POW_((1-xb),4)*_POW_(xb,bottom3)*_EXP_(1.0+bottom4*xb);
  afloat_t FuncABc = Funcsa*Funcsb/xa/xb;
  afloat_t FuncABs = Funcca*Funccb/xa/xb;
  afloat_t FuncABb = Funcba*Funcbb/xa/xb;
  
  
  afloat_t totSec = 2*mVal*(
			  (FuncABu)*weightu
			  +(FuncABd)*weightd
			  +(FuncABc)*weightc
			  +(FuncABs)*weights
			  +(FuncABb)*weightb
			  );
  
  if(( mVal <= 600. && fabs(YVal) > 20*_POW_(afloat_t(mVal),afloat_t(-0.32))) || ( mVal > 600. && fabs(YVal) > 21*_POW_(afloat_t(mVal),afloat_t(-0.34))))
  {
    //Find totSec when mZZ, YVal=0
    afloat_t xa0 = _SQRT_(s/s0); //at YVal=0 xa=xb
    
    //up
    //if xa=xb then FuncAu1=FuncAu2 and FuncBu1=FuncBu2
    FuncAu1 = (up0+up1*xa0+up2*_POW_(xa0,2))*_POW_((1-xa0),4)*_POW_(xa0,up3)*_EXP_(1.0+up4*xa0);
    FuncBu1 = (antiup0+antiup1*xa0+antiup2*_POW_(xa0,2)+antiup3*_POW_(xa0,3))*_POW_((1-xa0),4)*_POW_(xa0,antiup4)*_EXP_(1.0+antiup5*xa0);
    FuncABu = 2*(FuncAu1/xa0*FuncBu1/xa0);
    
    //down
    //if xa=xb then FuncAd1=FuncAd2 and FuncBd1=FuncBd2
    FuncAd1 = (down0+down1*xa0+down2*_POW_(xa0,2))*_POW_((1-xa0),4)*_POW_(xa0,down3)*_EXP_(1.0+down4*xa0);
    FuncBd1 = (antidown0+antidown1*xa0+antidown2*_POW_(xa0,2)+antidown3*_POW_(xa0,3))*_POW_((1-xa0),4)*_POW_(xa0,antidown4)*_EXP_(1.0+antidown5*xa0);
    FuncABd = 2*(FuncAd1/xa0*FuncBd1/xa0);
    
    //sea
    Funcca = (charm0+charm1*xa0+charm2*_POW_(xa0,2)+charm3*_POW_(xa0,3))*_POW_((1-xa0),4)*_POW_(xa0,charm4)*_EXP_(1.0+charm5*xa0); //Funcca=Funccb
    Funcsa = (strange0+strange1*xa0+strange2*_POW_(xa0,2))*_POW_((1-xa0),4)*_POW_(xa0,strange3)*_EXP_(1.0+strange4*xa0); //Funcsa=Funcsb
    Funcba = (bottom0+bottom1*xa0+bottom2*_POW_(xa0,2))*_POW_((1-xa0),4)*_POW_(xa0,bottom3)*_EXP_(1.0+bottom4*xa0); //Funcba=Funcbb
    FuncABc = Funcsa*Funcsa/xa0/xa0;
    FuncABs = Funcca*Funcca/xa0/xa0;
    FuncABb = Funcba*Funcba/xa0/xa0;
    afloat_t totSec0 = 2*mVal*(
			     (FuncABu)*weightu
			     +(FuncABd)*weightd
			     +(FuncABc)*weightc
			     +(FuncABs)*weights
			     +(FuncABb)*weightb
			     );
    totSec = 1.e-5*totSec0;
  }
  
  if (totSec <= 0.) totSec = 0.00001;
  
  totSec*=2.; // correcting for the half weight factor
  return totSec;
}
