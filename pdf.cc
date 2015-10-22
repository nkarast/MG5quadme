#include "common.h"

_CUDA_HOST_ _CUDA_DEVICE_ afloat_t pdf(afloat_t x, afloat_t Q, afloat_t pdf_data[], int pdf_nptsx, int pdf_nptsq, afloat_t pdf_lxmin, afloat_t pdf_lxmax, afloat_t pdf_lqmin, afloat_t pdf_lqmax)
{
  
  /*
   * bilinear interpolation of PDF grid provided as flattened m x n array: [q][x] -> [q_{1},x_{1}, ... q_{1},x_{n}, ... q_{m},x_{1}, ... q_{m},x_{n}]
   *
   * *** NOTE *** expectation is that grid is provided in log coordinates
   *
   */
  
  x = _LOG_(x);
  Q = _LOG_(Q);
  
  afloat_t dx = (pdf_lxmax - pdf_lxmin) / pdf_nptsx;
  unsigned ix = 0;
  if(pdf_lxmin < x) ix = static_cast<unsigned int>((x - pdf_lxmin) / dx);
  ix = ix < pdf_nptsx-1 ? ix : pdf_nptsx-2;
  
  afloat_t dq = (pdf_lqmax - pdf_lqmin) / pdf_nptsq;
  unsigned iq = 0;
  if(pdf_lqmin < Q) iq = static_cast<unsigned int>((Q - pdf_lqmin) / dq);
  iq = iq < pdf_nptsq-1 ? iq : pdf_nptsq-2;
  
  afloat_t c11x = pdf_lxmin + dx*ix;
  afloat_t c11q = pdf_lqmin + dq*iq;
  
  afloat_t c22x = pdf_lxmin + dx*(ix+1);
  afloat_t c22q = pdf_lqmin + dq*(iq+1);
  
  afloat_t norm = dx*dq;
  
  unsigned int i0j0 = iq*pdf_nptsx + ix;
  unsigned int i1j0 = (iq+1)*pdf_nptsx + ix;
  unsigned int i0j1 = iq*pdf_nptsx + ix + 1;
  unsigned int i1j1 = (iq+1)*pdf_nptsx + ix + 1;
  
  return (pdf_data[i0j0] / (norm) * (c22x - x)*(c22q - Q) +
	  pdf_data[i0j1] / (norm) * (x - c11x)*(c22q - Q) + 
	  pdf_data[i1j0] / (norm) * (c22x - x)*(Q - c11q) + 
	  pdf_data[i1j1] / (norm) * (x - c11x)*(Q - c11q));
}

// ------------------------- ======= ------------------------- ======= -------------------------

#ifndef _CUDA_READY_

#include "pdf.h"

// ------------------------- ======= ------------------------- ======= -------------------------
#ifdef USECTEQ
extern "C"
{
  void* setctq6_(int&);
  double ctq6pdf_(int&, double&, double&);
  void* setct10_(int&);
  double ct10pdf_(int&, double&, double&);
}

//   1  CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
//   2  CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
//   3  CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
//   4  CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
// 100  CT10     central NLO             0.118                  ct10.00.pds  

namespace 
{
  int CT10 = 100;
  int CTEQ6M  = 1;
  int CTEQ6D  = 2;
  int CTEQ6L  = 3;
  int CTEQ6L1 = 4;
  
  typedef double (*pdf_fun)(int&, double&, double&);
  
  pdf_fun ctpdf;
}
#else
#include "lhapdf/LHAPDF.h"
#endif
// ------------------------- ======= ------------------------- ======= -------------------------

PDFGrid* PDFGrid::m_instance = NULL;

PDFGrid* PDFGrid::load(const std::string& pdfn, PartonType t, unsigned nx, unsigned nQ, unsigned int isub)
{
  if(!m_instance)
  {
    m_instance = new PDFGrid(pdfn, nx, nQ, isub, t);
    return m_instance;
  }
  if(m_instance->get_pdfn() != pdfn ||
     m_instance->get_isub() != isub ||
     m_instance->get_nx() != nx || m_instance->get_nQ() != nQ)
  {
    delete m_instance;
    m_instance = new PDFGrid(pdfn, nx, nQ, isub, t);
    return m_instance;
  }
  return m_instance;
}

// ------------------------- ======= ------------------------- ======= -------------------------

PDFGrid::PDFGrid(const std::string& pdfn, unsigned nx, unsigned nQ, unsigned int isub, PartonType t) :
  m_nx(nx),
  m_nQ(nQ),
  m_pdfn(pdfn),
  m_isub(isub)
{
#ifndef USECTEQ
  LHAPDF::initPDFSet(pdfn.c_str(), LHAPDF::LHGRID, isubset);
  m_lxmin = std::log(LHAPDF::getXmin(isubset));
  m_lxmax = 0;
  m_lQmin = std::log(std::sqrt(LHAPDF::getQ2min(isubset)));
  m_lQmax = std::log(std::sqrt(LHAPDF::getQ2max(isubset)));
#else
  std::cout << "INFO :: initializing " << pdfn << " PDF table" << std::endl;
  if(pdfn == "CTEQ6L" || pdfn == "cteq6l") {
    setctq6_(CTEQ6L); ctpdf = ctq6pdf_; }
  if(pdfn == "CTEQ6D" || pdfn == "cteq6d") {
    setctq6_(CTEQ6D); ctpdf = ctq6pdf_; }
  if(pdfn == "CTEQ6M" || pdfn == "cteq6m") {
    setctq6_(CTEQ6M); ctpdf = ctq6pdf_; }
  if(pdfn == "CTEQ6L1" || pdfn == "cteq6l1") {
    setctq6_(CTEQ6L1); ctpdf = ctq6pdf_; }
  if(pdfn == "CT10" || pdfn == "ct10") {
    setct10_(CT10); ctpdf = ct10pdf_; }
  m_lxmin = std::log(1.0e-6);
  m_lxmax = 0;
  m_lQmin = std::log(1.3e00);
  m_lQmax = std::log(1.0e04); 
#endif
  
  m_dlQ = (m_lQmax - m_lQmin) / m_nQ;
  m_dlx = (m_lxmax - m_lxmin) / m_nx;
  
  if(t == kunknown)
  {
    m_data_vec.resize(10);
    unsigned int ipart=0;
    init(kd, &(m_data_vec[ipart++]));
    init(ku, &(m_data_vec[ipart++]));
    init(ks, &(m_data_vec[ipart++]));
    init(kc, &(m_data_vec[ipart++]));
    
    init(kdbar, &(m_data_vec[ipart++]));
    init(kubar, &(m_data_vec[ipart++]));
    init(ksbar, &(m_data_vec[ipart++]));
    init(kcbar, &(m_data_vec[ipart++]));
    
    init(kgluon, &(m_data_vec[ipart++]));
  }
  else
  {
    m_data_vec.push_back(grid_t());
    init(t, &(m_data_vec.back()));
  }
}

void PDFGrid::init(PartonType t, grid_t* data)
{
  // #ifdef USECTEQ
  //   return;
  // #else
  int id = static_cast<int>(t);
  afloat_t x;
  afloat_t Q;
  for(unsigned int iQ = 0; iQ < m_nQ; ++iQ)
  {
    afloat_t lQ = m_lQmin + iQ * m_dlQ;
    row_t buffer;
    for(unsigned int ix = 0; ix < m_nx; ++ix)
    {
      afloat_t lx = m_lxmin + ix * m_dlx;
#ifndef USECTEQ
      buffer.push_back(LHAPDF::xfx(exp(lx), exp(lQ), static_cast<int>(t)) / (exp(lx)));
#else
      x = exp(lx);
      Q = exp(lQ);
      double _xd(x); 
      double _Qd(Q);
      buffer.push_back(ctpdf(id, _xd, _Qd));
#endif
    }
    data->push_back(buffer);
  }
  m_grid_map[static_cast<int>(t)] = data;
  std::cout << "\tloaded PDF for " << id << " with " 
            << data->size() << "x" << data->at(0).size() << " grid" << std::endl;
  // #endif
}

afloat_t* PDFGrid::get_grid_data(PartonType t) const
{
  afloat_t* grd = new afloat_t[m_nQ*m_nx];
  for(unsigned int iQ=0; iQ < m_nQ; ++iQ)
    for(unsigned int ix=0; ix < m_nx; ++ix)
      grd[iQ*m_nx+ix] = (m_grid_map.find(t)->second)->at(iQ).at(ix);
  return grd;
}

unsigned int PDFGrid::get_ix(afloat_t lx) const
{
  if(lx <= m_lxmin) return 0;
  if(lx >= m_lxmax) return m_nx-2; 
  unsigned int ix = static_cast<unsigned int>((lx - m_lxmin) / m_dlx);
  return ix < m_nx-1 ? ix : m_nx-2;
}

unsigned int PDFGrid::get_iQ(afloat_t lQ) const
{
  if(lQ <= m_lQmin) return 0;
  if(lQ >= m_lQmax) return m_nQ-2;
  unsigned int iQ = static_cast<unsigned int>((lQ - m_lQmin) / m_dlQ);
  return iQ < m_nQ-1 ? iQ : m_nQ-2;
}

afloat_t PDFGrid::lookup(unsigned int ix, unsigned int iQ, const grid_t* data) const
{
  return (*data)[iQ][ix];
}  

afloat_t PDFGrid::operator()(PartonType t, afloat_t x, afloat_t Q) const
{
  return this->operator()(static_cast<int>(t), x, Q);
}

afloat_t PDFGrid::operator()(int t, afloat_t x, afloat_t Q) const
{
  // #ifdef USECTEQ
  //   return ctpdf(t, x, Q);
  // #else
  afloat_t lx = std::log(x);
  afloat_t lQ = std::log(Q);
  
  unsigned int ix = get_ix(lx);
  unsigned int iQ = get_iQ(lQ);
  
  afloat_t c11x = m_lxmin + m_dlx * ix;
  afloat_t c11q = m_lQmin + m_dlQ * iQ;
  
  afloat_t c22x = m_lxmin + m_dlx * (ix+1);
  afloat_t c22q = m_lQmin + m_dlQ * (iQ+1);
  
  afloat_t norm = m_dlx*m_dlQ;
  
  const grid_t* data = m_grid_map.find(t)->second;
  
  return ((*data)[iQ][ix]     / (norm) * (c22x - lx)*(c22q - lQ) +
	  (*data)[iQ][ix+1]   / (norm) * (lx - c11x)*(c22q - lQ) + 
	  (*data)[iQ+1][ix]   / (norm) * (c22x - lx)*(lQ - c11q) + 
	  (*data)[iQ+1][ix+1] / (norm) * (lx - c11x)*(lQ - c11q));
  // #endif
}

#endif
