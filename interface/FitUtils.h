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

class FitUtils {

 public:

  FitUtils();
  ~FitUtils();

  void addPdfSgn(TH1D*);
  void addPdfBkg(TH1D*);
  void addTransformPdf(TH1D*);

  void addPdf(TTree*, const string&, const bool& = bool(false));

  void run_test(const int&);
  void run_test_unbinned(const int&);

 private:

  RooRealVar* x;
  RooRealVar* x_transf;
  TH1D* h_pdf_sgn;
  TH1D* h_pdf_bkg;
  TH1D* cum_pdf;  
  RooAbsPdf* pdf_sgn;
  RooAbsPdf* pdf_bkg;
  RooAbsReal* transformation;
  RooAbsPdf* pdf_sgn_transf;
  RooAbsPdf* pdf_bkg_transf;

  TH1D* h_pdf_sgn_transf;
  TH1D* h_pdf_bkg_transf;

  TFile* f;
  bool savePlots;
};

#endif
