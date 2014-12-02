#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"

TGraph* roc(TString fname="Test_toys50_smear0_btag0_gen1_test1.root", float step=0.1, float xMin=0., float xMax=20., 
	    int color=2 , TLegend* leg=0, TString leg_name=""){

  TFile* f = TFile::Open( fname , "READ" );
  if(f==0 || f->IsZombie()) return 0;

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return 0;

  TH1F h("h", "", 400, xMin, xMax);
  t->Draw("test__nll[0]-test__nll[1]>>h");
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

void roc_comp_two( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		   float step=0.1, float xMin=0., float xMax=20.,
		   TString title = "TopHad + TopLep",
		   TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
		   TString save_name = "tmp.png"
		   ){

  gStyle->SetOptStat(0);

  TLegend* leg = new TLegend(0.48,0.13,0.88,0.33,NULL,"brNDC");
  
  TH1F* hROC = new TH1F("hROC", title+";-2Log(L_{0}/L_{1}); P( -2Log(L_{0}/L_{1}) > x ) ", 100, xMin, xMax);

  TGraph* gROC1 = roc(fname1, step, xMin, xMax, 2, leg, leg_name1); 
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, 3, leg, leg_name2); 

  if(gROC1==0 || gROC2==0) return;

  hROC->Draw();
  gROC1->Draw("SAME");
  gROC2->Draw("SAME");
  leg->Draw();
  gPad->SaveAs( save_name );
}



void roc_comp_ROC( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
                   float step=0.1, float xMin=0., float xMax=20.,
                   TString title = "TopHad + TopLep",
                   TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
                   TString save_name = "tmp.png"
                   ){

  gStyle->SetOptStat(0);

  TLegend* leg = new TLegend(0.48,0.13,0.88,0.33,NULL,"brNDC");

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, 0, 1.0);
  hROC->SetMinimum(0.0);
  hROC->SetMaximum(1.0);

  TGraph* gROC1 = roc(fname1, step, xMin, xMax, 2, leg, leg_name1);
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, 3, leg, leg_name2);

  if(gROC1==0 || gROC2==0) return;

  TGraph* gROC = new TGraph( gROC1->GetN() , gROC1->GetY(), gROC2->GetY());
  gROC->SetLineWidth(3);
  gROC->SetLineColor(kRed);

  hROC->Draw();
  gROC->Draw("SAME");
  gPad->SaveAs( save_name );
}


void roc_comp_four( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root",
		     TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		     TString fname3="Test_toys50_smear1_btag1_gen1_test1.root",
		     TString fname4="Test_toys50_smear1_btag2_gen1_test1.root",
		     float step=0.1, float xMin=0., float xMax=10.,
		     TString title = "TopHad + TopLep",
		     TString leg_name1="Parton",  TString leg_name2="Smear", TString leg_name3="Smear+btag", TString leg_name4="Smear+btag+filter",
		     TString save_name = "tmp.png"
){

  gStyle->SetOptStat(0);
  TLegend* leg = new TLegend(0.48,0.13,0.88,0.33,NULL,"brNDC");
  TH1F* hROC   = new TH1F("hROC", title+";-2Log(L_{0}/L_{1}); P( -2Log(L_{0}/L_{1}) > x ) ", 100, xMin, xMax);
  TGraph* gROC1 = roc(fname1, step, xMin, xMax, 2, leg, leg_name1);
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, 3, leg, leg_name2);
  TGraph* gROC3 = roc(fname3, step, xMin, xMax, 4, leg, leg_name3);
  TGraph* gROC4 = roc(fname4, step, xMin, xMax, 5, leg, leg_name4);
  if(gROC1==0 || gROC2==0 || gROC3==0 || gROC4==0) return;

  hROC->Draw(); 
  gROC1->Draw("SAME");
  gROC2->Draw("SAME");
  gROC3->Draw("SAME");
  gROC4->Draw("SAME");
  leg->Draw();
  gPad->SaveAs( save_name );
}


void doall(TString dir="./" ){

  roc_comp_two( dir+"Test_toys50_smear0_btag0_gen1_test1.root", dir+"Test_toys50_smear0_btag0_gen6_test1.root",
		0.1, 0., 10.,
		"TopHad + TopLep (parton)",
		"H_{0}", "H_{1}", dir+"TopHad_TopLep_H0vsH1_parton.png");   
  
  roc_comp_two( dir+"Test_toys50_smear1_btag0_gen1_test1.root", dir+"Test_toys50_smear1_btag0_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_H0vsH1_smear.png");
  
  roc_comp_two( dir+"Test_toys50_smear1_btag1_gen1_test1.root", dir+"Test_toys50_smear1_btag1_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear+btag)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_H0vsH1_btag1.png");
  
  roc_comp_two( dir+"Test_toys50_smear1_btag2_gen1_test1.root", dir+"Test_toys50_smear1_btag2_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear+btag+filter)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_H0vsH1_btag2.png");


  roc_comp_two( dir+"Test_toys50_smear0_btag0_gen2_test2.root", dir+"Test_toys50_smear0_btag0_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (parton)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_H0vsH1_parton.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag0_gen2_test2.root", dir+"Test_toys50_smear1_btag0_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_H0vsH1_smear.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag1_gen2_test2.root", dir+"Test_toys50_smear1_btag1_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear+btag)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_H0vsH1_btag1.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag2_gen2_test2.root", dir+"Test_toys50_smear1_btag2_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear+btag+filter)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_H0vsH1_btag2.png");

  
  roc_comp_two( dir+"Test_toys50_smear0_btag0_gen3_test3.root", dir+"Test_toys50_smear0_btag0_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (parton)",
                "H_{0}", "H_{1}", "TopLep_TopLep_Higgs_H0vsH1_parton.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag0_gen3_test3.root", dir+"Test_toys50_smear1_btag0_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear)",
                "H_{0}", "H_{1}", "TopLep_TopLep_Higgs_H0vsH1_smear.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag1_gen3_test3.root", dir+"Test_toys50_smear1_btag1_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear+btag)",
                "H_{0}", "H_{1}", "TopLep_TopLep_Higgs_H0vsH1_btag1.png");

  roc_comp_two( dir+"Test_toys50_smear1_btag2_gen3_test3.root", dir+"Test_toys50_smear1_btag2_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear+btag+filter)",
                "H_{0}", "H_{1}", "TopLep_TopLep_Higgs_H0vsH1_btag2.png");

  
  ///////////////////////////////////////////////////////////////////////////////

  roc_comp_ROC( dir+"Test_toys50_smear0_btag0_gen1_test1.root", dir+"Test_toys50_smear0_btag0_gen6_test1.root",
		0.1, 0., 10.,
		"TopHad + TopLep (parton)",
		"H_{0}", "H_{1}", dir+"TopHad_TopLep_ROC_parton.png");   
  
  roc_comp_ROC( dir+"Test_toys50_smear1_btag0_gen1_test1.root", dir+"Test_toys50_smear1_btag0_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_ROC_smear.png");
  
  roc_comp_ROC( dir+"Test_toys50_smear1_btag1_gen1_test1.root", dir+"Test_toys50_smear1_btag1_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear+btag)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_ROC_btag1.png");
  
  roc_comp_ROC( dir+"Test_toys50_smear1_btag2_gen1_test1.root", dir+"Test_toys50_smear1_btag2_gen6_test1.root",
                0.1, 0., 10.,
                "TopHad + TopLep (smear+btag+filter)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_ROC_btag2.png");


  roc_comp_ROC( dir+"Test_toys50_smear0_btag0_gen2_test2.root", dir+"Test_toys50_smear0_btag0_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (parton)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_ROC_parton.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag0_gen2_test2.root", dir+"Test_toys50_smear1_btag0_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_ROC_smear.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag1_gen2_test2.root", dir+"Test_toys50_smear1_btag1_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear+btag)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_ROC_btag1.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag2_gen2_test2.root", dir+"Test_toys50_smear1_btag2_gen7_test2.root",
                0.1, 0., 10.,
                "TopHad + TopLep + Higgs (smear+btag+filter)",
                "H_{0}", "H_{1}", dir+"TopHad_TopLep_Higgs_ROC_btag2.png");

  
  roc_comp_ROC( dir+"Test_toys50_smear0_btag0_gen3_test3.root", dir+"Test_toys50_smear0_btag0_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (parton)",
                "H_{0}", "H_{1}", dir+"TopLep_TopLep_Higgs_ROC_parton.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag0_gen3_test3.root", dir+"Test_toys50_smear1_btag0_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear)",
                "H_{0}", "H_{1}", dir+"TopLep_TopLep_Higgs_ROC_smear.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag1_gen3_test3.root", dir+"Test_toys50_smear1_btag1_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear+btag)",
                "H_{0}", "H_{1}", dir+"TopLep_TopLep_Higgs_ROC_btag1.png");

  roc_comp_ROC( dir+"Test_toys50_smear1_btag2_gen3_test3.root", dir+"Test_toys50_smear1_btag2_gen8_test3.root",
                0.1, 0., 10.,
                "TopLep + TopLep + Higgs (smear+btag+filter)",
                "H_{0}", "H_{1}", dir+"TopLep_TopLep_Higgs_ROC_btag2.png");

  
  /////////////////////////////////////////////////////////////////




  roc_comp_four( dir+"Test_toys50_smear0_btag0_gen1_test1.root", dir+"Test_toys50_smear1_btag0_gen1_test1.root",
		  dir+"Test_toys50_smear1_btag1_gen1_test1.root", dir+"Test_toys50_smear1_btag2_gen1_test1.root",
		  0.1, 0., 5.,
		  "TopHad + TopLep",
		  "Parton", "Smear", "Smear+btag","Smear+btag+filter", dir+"TopHad_TopLep_TF.png");

  roc_comp_four( dir+"Test_toys50_smear0_btag0_gen2_test2.root", dir+"Test_toys50_smear1_btag0_gen2_test2.root",
                  dir+"Test_toys50_smear1_btag1_gen2_test2.root", dir+"Test_toys50_smear1_btag2_gen2_test2.root",
                  0.1, 0., 5.,
                  "TopHad + TopLep + Higgs",
                  "Parton", "Smear", "Smear+btag","Smear+btag+filter", dir+"TopHad_TopLep_Higgs_TF.png");

  roc_comp_four( dir+"Test_toys50_smear0_btag0_gen3_test3.root", dir+"Test_toys50_smear1_btag0_gen3_test3.root",
                  dir+"Test_toys50_smear1_btag1_gen3_test3.root", dir+"Test_toys50_smear1_btag2_gen3_test3.root",
                  0.1, 0., 5.,
                  "TopLep + TopLep + Higgs",
                  "Parton", "Smear", "Smear+btag","Smear+btag+filter", dir+"TopLep_TopLep_Higgs_TF.png");



}
