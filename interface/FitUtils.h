#ifndef FITUTILS_H
#define FITUTILS_H

#include<iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooExtendPdf.h"
#include "RooRealConstant.h"

using namespace RooFit ; 
using namespace std;

namespace Option {
  enum Option { Binned=0, Unbinned=1 }; 
}
namespace Process{
  enum Process { Signal=0, Background=1, Transform=2 }; 
}

class FitUtils {

 public:

  FitUtils();
  ~FitUtils();

  void save_plots(const int&);
  void addPdf(TTree*, const Process::Process&, const double& =double(1.), const Option::Option& =Option::Unbinned);
  void run_test(const int&, const Option::Option&, const size_t& =2);

 private:

  RooRealVar* x;
  RooRealVar* x_transf;
  RooAbsPdf* pdf_sgn;
  RooAbsPdf* pdf_bkg;
  RooAbsReal* transformation;
  RooAbsPdf* pdf_sgn_transf;
  RooAbsPdf* pdf_bkg_transf;

  double norm_sgn;
  double norm_bkg;

  int savePlots;

  TFile* file;
  TTree* tree;

  float nsgn_gen;
  float nsgn_fit;
  float nsgn_err;
  float nbkg_gen;
  float nbkg_fit;
  float nbkg_err;
  float minNll;
  float edm;
  int status;

};

#endif
