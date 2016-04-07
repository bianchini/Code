#include "TTH/MEIntegratorStandalone/interface/FitUtils.h"
#include<iostream>
#include "TH1D.h"
#include "TTree.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, char *argv[]){

  TRandom3* ran = new TRandom3();  

  FitUtils* fit = new FitUtils();

  /*
  TH1D h_sgn_pdf("h_pdf_sgn","sgn pdf", 500,-3,3);
  for(int i = 0 ; i < 100000; ++i)
    h_sgn_pdf.Fill( ran->Gaus(+2.,2.) );
  h_sgn_pdf.Scale(10./h_sgn_pdf.Integral());

  TH1D h_bkg_pdf("h_pdf_bkg","bkg pdf", 500,-3,3);
  for(int i = 0 ; i < 100000; ++i)
    h_bkg_pdf.Fill( ran->Gaus(-2.,2.) );
  h_bkg_pdf.Scale(100./h_bkg_pdf.Integral());

  TH1D h_transform_pdf("h_transform_bkg","pdf for the transformation", 1000,-3,3);
  for(int i = 0 ; i < 100000; ++i)
    h_transform_pdf.Fill( ran->Gaus(-1.5,2.) );

  fit->addTransformPdf(&h_transform_pdf);
  fit->addPdfBkg(&h_bkg_pdf);
  fit->addPdfSgn(&h_sgn_pdf);
  */

  TTree* tree_sgn = new TTree("tree_sgn", "signal tree");
  TTree* tree_bkg = new TTree("tree_bkg", "background tree");
  TTree* tree_tra = new TTree("tree_tra", "transform tree");
  double x;
  tree_sgn->Branch("x", &x, "x/D");
  tree_bkg->Branch("x", &x, "x/D");
  tree_tra->Branch("x", &x, "x/D");
  for(int i = 0 ; i < 5000; ++i){
    x = ran->Gaus(+2.,2.);
    tree_sgn->Fill();
    x = ran->Gaus(-2.,2.);
    tree_bkg->Fill();
    x = ran->Gaus(-2.,2.);
    tree_tra->Fill();
  }
  cout << "Make Transf Pdf..." << endl;
  fit->addPdf(tree_tra, "transf");
  cout << "Make Bkg Pdf..." << endl;
  fit->addPdf(tree_bkg, "bkg");
  cout << "Make Sgn Pdf..." << endl;
  fit->addPdf(tree_sgn, "sgn");
  
  fit->run_test(1);

  delete fit;
  delete ran;
}
