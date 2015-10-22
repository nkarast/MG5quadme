#ifndef _rambo_h_
#define _rambo_h_

#include <vector>

using namespace std;

class Random
{
  public:
    double ranmar();
    void rmarin(int ij, int kl);
    
  private:
    double ranu[98];
    double ranc, rancd, rancm;
    int iranmr, jranmr;
};


double rn(int idummy);

vector<vector<double> > get_momenta(int ninitial, double energy, 
				    vector<double> masses, double* wgt);

vector<vector<double> > rambo(double et, vector<double>& xm, double* wt);

#ifdef __CINT__
#pragma link C++ class std::vector<std::vector<double> >; 
#endif 

#endif
