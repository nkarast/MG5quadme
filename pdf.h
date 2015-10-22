#include "common.h"

#ifndef __CINT__
_CUDA_HOST_ _CUDA_DEVICE_ afloat_t pdf(afloat_t x, afloat_t Q, afloat_t pdf_data[], int pdf_nptsx, int pdf_nptsq, afloat_t pdf_lxmin, afloat_t pdf_lxmax, afloat_t pdf_lqmin, afloat_t pdf_lqmax);
#endif

#ifndef _CUDA_READY_

#ifndef _PDF_H
#define _PDF_H

#define USECTEQ // use CTEQ pdf routines directly (instead of LHAPDF) << faster

#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>

/*
 * a unified interface to CTEQ and LHAPDF PDF libraries and datasets
 */

#define NUMXSAMPLES 4000
#define NUMQSAMPLES 100

class PDFGrid
{
public:
  
  enum PartonType
    {
      kunknown = 100,
      kgluon = 0,
      kd = 1,
      ku = 2,
      ks = 3,
      kc = 4,
      kb = 5,
      kdbar = -1,
      kubar = -2,
      ksbar = -3,
      kcbar = -4,
      kbbar = -5
    };
  
  struct coord_t 
  {
  public:
    coord_t(afloat_t mlx, afloat_t mlQ)
    {
      lx = mlx;
      lQ = mlQ;
    }
    afloat_t lx;
    afloat_t lQ;
  };
  
  static PDFGrid* load(const std::string& pdfn, PartonType t=kunknown, unsigned nx = NUMXSAMPLES, unsigned nQ = NUMQSAMPLES, unsigned int isubset = 0);

  afloat_t operator()(PartonType t, afloat_t x, afloat_t Q) const;
  afloat_t operator()(int t, afloat_t x, afloat_t Q) const;
  
  afloat_t u(afloat_t x, afloat_t Q) const { return (*this)(ku, x, Q); }
  afloat_t d(afloat_t x, afloat_t Q) const { return (*this)(kd, x, Q); }
  afloat_t s(afloat_t x, afloat_t Q) const { return (*this)(ks, x, Q); }
  afloat_t c(afloat_t x, afloat_t Q) const { return (*this)(kc, x, Q); }
  
  afloat_t ubar(afloat_t x, afloat_t Q) const { return (*this)(kubar, x, Q); }
  afloat_t dbar(afloat_t x, afloat_t Q) const { return (*this)(kdbar, x, Q); }
  afloat_t sbar(afloat_t x, afloat_t Q) const { return (*this)(ksbar, x, Q); }
  afloat_t cbar(afloat_t x, afloat_t Q) const { return (*this)(kcbar, x, Q); }
  
  afloat_t g(afloat_t x, afloat_t Q) const { return (*this)(kgluon, x, Q); }
  
  afloat_t get_dx() const { return exp(m_dlx); }
  afloat_t get_dQ() const { return exp(m_dlQ); }
  
  afloat_t get_xmin() const { return exp(m_lxmin); }
  afloat_t get_xmax() const { return exp(m_lxmax); }
  afloat_t get_Qmin() const { return exp(m_lQmin); }
  afloat_t get_Qmax() const { return exp(m_lQmax); }
  
  unsigned int get_nx() const { return m_nx; }
  unsigned int get_nQ() const { return m_nQ; }
  
  std::string get_pdfn() const { return m_pdfn; }
  unsigned int get_isub() const { return m_isub; }

  afloat_t* get_grid_data(PartonType t) const;
  
private:
  
  typedef std::vector< std::vector<afloat_t> > grid_t;
  typedef std::vector<afloat_t> row_t;

  static PDFGrid* m_instance;
  
  PDFGrid(const std::string& pdfn, unsigned nx, unsigned nQ, unsigned int isub, PartonType t=kunknown);
  
  PDFGrid() { }
  ~PDFGrid() { }
  
  PDFGrid(const PDFGrid&);
  PDFGrid& operator=(const PDFGrid&);
  
  unsigned int m_nx;
  unsigned int m_nQ;
  
  afloat_t m_dlx;
  afloat_t m_dlQ;
  
  afloat_t m_lxmax;
  afloat_t m_lxmin;
  afloat_t m_lQmax;
  afloat_t m_lQmin;
  
  // grid_t m_u_data;
  // grid_t m_d_data;
  // grid_t m_s_data;
  // grid_t m_c_data;
  
  // grid_t m_ubar_data;
  // grid_t m_dbar_data;
  // grid_t m_sbar_data;
  // grid_t m_cbar_data;
  
  // grid_t m_g_data;

  std::vector<grid_t> m_data_vec;
  
  std::map<int, const grid_t*> m_grid_map;
  
  std::string m_pdfn;
  unsigned int m_isub;
  
  unsigned int get_ix(afloat_t x) const;
  unsigned int get_iQ(afloat_t Q) const;
  
  void init(PartonType t, grid_t* data);
  
  afloat_t lookup(unsigned int ix, unsigned int iQ, const grid_t*) const;
};

#endif 
#endif 
