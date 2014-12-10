#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"

void CanvasAndLegend(TCanvas* c1, TLegend* leg, int logy=0){
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
}

TGraph* roc(TString fname="Test_toys50_smear0_btag0_gen1_test1.root", float step=0.1, float xMin=0., float xMax=20., 
	    int color=2 , TLegend* leg=0, TString leg_name="", int H0=0, int H1=1){

  TFile* f = TFile::Open( fname , "READ" );
  if(f==0 || f->IsZombie()) return 0;

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return 0;

  TH1F h("h", "", 400, xMin, xMax);
  t->Draw(Form("test__nll[%d]-test__nll[%d]>>h", H0, H1));
  float OF = h.GetBinContent( h.GetNbinsX() ) +  h.GetBinContent( h.GetNbinsX()+1) ;
  float UF = h.GetBinContent( 1 ) +  h.GetBinContent( 0 ) ;
  h.SetBinContent( h.GetNbinsX(), OF);
  h.SetBinContent( 1, UF );

  int total = h.GetEntries();
  const int steps = (xMax-xMin)/step;
  double x[steps];
  double y[steps];

  if(total==0) return 0;
  for(int i = 0 ; i < steps ; i++){
    x[i] = xMin + i*step;
    y[i] = h.Integral(1, h.FindBin(x[i]+step/2.) )/total;
  }

  TGraph* gROC = new TGraph(steps, x, y);
  gROC->SetLineWidth(3);
  gROC->SetLineColor(color);

  if(leg!=0){
    leg->AddEntry(gROC, leg_name, "L" );
  }

  f->Close();

  return gROC;
}

void roc_comp_ROC( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		   int H0=0, int H1=1,		  
                   float step=0.05, float xMin=-10., float xMax=50.,
                   TString title = "TopHad + TopLep",
                   TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
                   TString save_name = "tmp.png"
                   ){

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  TLegend* leg = new TLegend(0.48,0.13,0.88,0.33,NULL,"brNDC");
  CanvasAndLegend(c1, leg, 0);

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, 0, 1.0);
  hROC->SetMinimum(0.0);
  hROC->SetMaximum(1.0);

  TGraph* gROC1 = roc(fname1, step, xMin, xMax, 2, leg, leg_name1, H0, H1);
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, 3, leg, leg_name2, H0, H1);

  if(gROC1==0 || gROC2==0) return;

  TGraph* gROC = new TGraph( gROC1->GetN() , gROC1->GetY(), gROC2->GetY());
  gROC->SetLineWidth(3);
  gROC->SetLineColor(kRed);

  hROC->Draw();
  gROC->Draw("SAME");
  leg->Draw();
  c1->SaveAs( save_name );
}

void roc_comp_ROC_2( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		     int H0_1=0, int H1_1=1,
		     TString fname3="Test_toys50_smear0_btag0_gen1_test1.root", TString fname4="Test_toys50_smear1_btag0_gen1_test1.root",
		     int H0_2=0, int H1_2=2,
		     float step=0.05, float xMin=-10., float xMax=50.,
		     TString title = "TopHad + TopLep",
		     TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
		     TString save_name = "tmp.png"
		     ){

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  TLegend* leg = new TLegend(0.48,0.13,0.88,0.33,NULL,"brNDC");
  CanvasAndLegend(c1, leg, 0);

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, 0, 1.0);
  hROC->SetMinimum(0.0);
  hROC->SetMaximum(1.0);

  TGraph* gROC1 = roc(fname1, step, xMin, xMax, 2, leg, leg_name1, H0_1, H1_1);
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, 3, leg, leg_name2, H0_1, H1_1);
  TGraph* gROC3 = roc(fname3, step, xMin, xMax, 2, leg, leg_name1, H0_2, H1_2);
  TGraph* gROC4 = roc(fname4, step, xMin, xMax, 3, leg, leg_name2, H0_2, H1_2);

  if(gROC1==0 || gROC2==0) return;
  if(gROC3==0 || gROC4==0) return;

  TGraph* gROC_1 = new TGraph( gROC1->GetN() , gROC1->GetY(), gROC2->GetY());
  gROC_1->SetLineWidth(3);
  gROC_1->SetLineColor(kRed);

  TGraph* gROC_2 = new TGraph( gROC3->GetN() , gROC3->GetY(), gROC4->GetY());
  gROC_2->SetLineWidth(3);
  gROC_2->SetLineColor(kBlue);

  hROC->Draw();
  gROC_1->Draw("SAME");
  gROC_2->Draw("SAME");

  leg->Clear();
  leg->AddEntry(gROC_1, leg_name1 , "L");
  leg->AddEntry(gROC_2, leg_name2 , "L");

  leg->Draw();
  c1->SaveAs( save_name );

}

void doall(TString dir="./" ){

  ///////////////////////////////////////////////////////////////////////////////

  roc_comp_ROC_2( dir+"Test_tH_tL_Hh_smear1_btag1_gen2_test2.root", dir+"Test_tH_tL_Hh_smear1_btag1_gen12_test2.root",
		  0, 2,
		  dir+"Test_tH_tL_Hh_smear1_btag1_gen7_test2.root", dir+"Test_tH_tL_Hh_smear1_btag1_gen12_test2.root",
		  1, 2,
		  0.05, -10., 60.,
		  "lepton + six jets",
		  "tt+H vs tt+jj", "tt+bb vs tt+jj", dir+"TopHad_TopLep_Higgs_ROC.png");

  /////////////////////////////////////////////////////////////////



}
