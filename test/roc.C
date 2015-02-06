#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCut.h"
#include "TLorentzVector.h"


void CanvasAndLegend(TCanvas* c1, TLegend* leg, int logy=0){
  c1->SetGrid(1,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy);
  c1->SetLogx(logy);
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

TGraphErrors* roc(TString fname="Test_toys50_smear0_btag0_gen1_test1.root", float step=0.1, float xMin=0., float xMax=20., TCut cut = "", 
		  int color=2 , TLegend* leg=0, TString leg_name="", TString var = "test__nll[1]-test__nll[0]"){

  TFile* f = TFile::Open( fname , "READ" );
  if(f==0 || f->IsZombie()) return 0;

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return 0;

  TH1F h("h", "", 800, xMin, xMax);
  //if(H0>=0 && H1>=0) 
  t->Draw(var+">>h", cut);
  //else
  //t->Draw("(-test__n_btag)>>h", cut);
  cout << "Histo filled" << endl;

  float OF = h.GetBinContent( h.GetNbinsX() ) +  h.GetBinContent( h.GetNbinsX()+1) ;
  float UF = h.GetBinContent( 1 ) +  h.GetBinContent( 0 ) ;
  h.SetBinContent( h.GetNbinsX(), OF);
  h.SetBinContent( 1, UF );

  int total = h.GetEntries();
  const int steps = (xMax-xMin)/step;
  double x[steps];
  double y[steps];
  double ex[steps];
  double ey[steps];

  if(total==0) return 0;
  cout << "Loop to fill eff. vs cut" << endl;
  for(int i = 0 ; i < steps ; i++){
    x[i] = xMin + i*step;
    ex[i] = 0.;
    y[i] = h.Integral(1, h.FindBin(x[i]+step/2.) )/total;
    ey[i] = sqrt(y[i]*(1-y[i])/total);
  }

  TGraphErrors* gROC = new TGraphErrors(steps, x, y, ex, ey);

  if(leg!=0){
    if( string(var.Data()).find("n_btag")!=string::npos ){
      gROC->SetLineWidth(3);
      gROC->SetLineColor(color);
      leg->AddEntry(gROC, leg_name, "L" );
    }
    else{
      gROC->SetMarkerSize(1);
      gROC->SetMarkerStyle(kFullCircle);
      leg->AddEntry(gROC, leg_name, "P" );
    }
  }

  f->Close();

  return gROC;
}


void roc_comp_ROC( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		   TString var = "test__nll[1]-test__nll[0]",		  
                   float step=0.05, float xMin=-10., float xMax=50., TCut cut="",
                   TString title = "TopHad + TopLep",
                   TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
                   TString save_name = "tmp.png",
		   float xLow = 1e-02, float yLow = 1e-04, int doLog=0
                   ){

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  TLegend* leg = new TLegend(0.12,0.72,0.47,0.89,NULL,"brNDC");
  CanvasAndLegend(c1, leg, doLog);

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, doLog ? xLow : 0., 1.0);
  hROC->SetMinimum(doLog ? yLow : 0.);
  hROC->SetMaximum(1.0);

  cout << "Doing ROC for mva" << endl;
  TGraphErrors* gROC1 = roc(fname1, step, xMin, xMax, cut, 2, leg, leg_name1+" vs "+leg_name2+" (nll)", var);
  TGraphErrors* gROC2 = roc(fname2, step, xMin, xMax, cut, 3, 0, leg_name2, var);

  cout << "Doing ROC for cut" << endl;
  TGraphErrors* gROC1cut = roc(fname1, step, xMin, xMax, cut, 2, leg, leg_name1+" vs "+leg_name2+" (cut)", "(-test__n_btag)");
  TGraphErrors* gROC2cut = roc(fname2, step, xMin, xMax, cut, 3, 0, leg_name2, "(-test__n_btag)");

  if(gROC1==0 || gROC2==0) return;

  TGraphErrors* gROC = new TGraphErrors( gROC1->GetN() , gROC1->GetY(), gROC2->GetY(), gROC1->GetEY(), gROC2->GetEY());
  gROC->SetLineWidth(3);
  gROC->SetLineColor(kRed);

  TGraphErrors* gROCcut = new TGraphErrors( gROC1cut->GetN() , gROC1cut->GetY(), gROC2cut->GetY(), gROC1cut->GetEY(), gROC2cut->GetEY());
  gROCcut->SetMarkerSize(1);
  gROCcut->SetMarkerStyle(kFullCircle);
  gROCcut->SetMarkerColor(kMagenta);

  hROC->Draw();
  gROC->Draw("SAME");
  gROCcut->Draw("EPSAME");

  TF1* diag = new TF1("diag", "x", 0.,1.);
  diag->SetLineColor(kBlack);
  diag->SetLineStyle(kDashed);
  diag->Draw("SAME");

  leg->Draw();
  c1->SaveAs( save_name );

  TFile* out = TFile::Open(save_name+".root", "RECREATE");
  out->cd();
  c1->Write();
  gROC1->Write("cut_vs_eff_nll_0");
  gROC2->Write("cut_vs_eff_nll_1");
  gROC1cut->Write("cut_vs_eff_cut_2");
  gROC2cut->Write("cut_vs_eff_cut_3");
  out->Write();
  out->Close();
}

/*
void roc_comp_ROC_2( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		     int H0_1=0, int H1_1=1,
		     TString fname3="Test_toys50_smear0_btag0_gen1_test1.root", TString fname4="Test_toys50_smear1_btag0_gen1_test1.root",
		     int H0_2=0, int H1_2=2,
		     float step=0.05, float xMin=-10., float xMax=50., TCut cut = "",
		     TString title = "TopHad + TopLep",
		     TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
		     TString save_name = "tmp.png", int doLog=0
		     ){

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  TLegend* leg = new TLegend(0.12,0.72,0.47,0.89,NULL,"brNDC");
  CanvasAndLegend(c1, leg, doLog);

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, 0, 1.0);
  hROC->SetMinimum(0.0);
  hROC->SetMaximum(1.0);

  TGraph* gROC1 = roc(fname1, step, xMin, xMax, cut, 2, leg, leg_name1, H0_1, H1_1);
  TGraph* gROC2 = roc(fname2, step, xMin, xMax, cut, 3, leg, leg_name2, H0_1, H1_1);
  TGraph* gROC3 = roc(fname3, step, xMin, xMax, cut, 2, leg, leg_name1, H0_2, H1_2);
  TGraph* gROC4 = roc(fname4, step, xMin, xMax, cut, 3, leg, leg_name2, H0_2, H1_2);

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

  TF1* diag = new TF1("diag", "x", 0.,1.);
  diag->SetLineColor(kBlack);
  diag->SetLineStyle(kDashed);
  diag->Draw("SAME");

  leg->Draw();
  c1->SaveAs( save_name );

  TFile* out = TFile::Open(save_name+".root", "RECREATE");
  out->cd();
  c1->Write();
  gROC1->Write("cut_vs_eff_0");
  gROC2->Write("cut_vs_eff_1");
  gROC3->Write("cut_vs_eff_2");
  gROC4->Write("cut_vs_eff_3");
  out->Write();
  out->Close();
}
*/

void roc_comp_ROC_4( TString fname1="Test_toys50_smear0_btag0_gen1_test1.root", TString fname2="Test_toys50_smear1_btag0_gen1_test1.root",
		     TString var1 = "test__nll[1]-test__nll[0]",

		     TString fname3="Test_toys50_smear0_btag0_gen1_test1.root", TString fname4="Test_toys50_smear1_btag0_gen1_test1.root",
		     TString var2 = "test__nll[1]-test__nll[0]",

		     TString fname5="Test_toys50_smear0_btag0_gen1_test1.root", TString fname6="Test_toys50_smear1_btag0_gen1_test1.root",
		     TString var3 = "test__nll[1]-test__nll[0]",

		     TString fname7="Test_toys50_smear0_btag0_gen1_test1.root", TString fname8="Test_toys50_smear1_btag0_gen1_test1.root",
		     TString var4 = "test__nll[1]-test__nll[0]",

		     float step=0.05, float xMin=-10., float xMax=50., TCut cut = "",
		     TString title = "TopHad + TopLep",
		     TString leg_name1="H_{0}",  TString leg_name2="H_{1}",
		     TString leg_name3="H_{2}",  TString leg_name4="H_{4}",
		     TString save_name = "tmp.png", int doLog=0		     
		     ){

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  TLegend* leg = new TLegend(0.12,0.72,0.47,0.89,NULL,"brNDC");
  CanvasAndLegend(c1, leg, doLog);

  TH1F* hROC = new TH1F("hROC", title+"; H_{0}; H_{1} ", 100, 0, 1.0);
  hROC->SetMinimum(0.0);
  hROC->SetMaximum(1.0);

  TGraphErrors* gROC1 = roc(fname1, step, xMin, xMax, cut, 2, leg, leg_name1, var1);
  TGraphErrors* gROC2 = roc(fname2, step, xMin, xMax, cut, 3, leg, leg_name2, var1);

  TGraphErrors* gROC3 = roc(fname3, step, xMin, xMax, cut, 2, leg, leg_name1, var2);
  TGraphErrors* gROC4 = roc(fname4, step, xMin, xMax, cut, 3, leg, leg_name2, var2);

  TGraphErrors* gROC5 = roc(fname5, step, xMin, xMax, cut, 2, leg, leg_name1, var3);
  TGraphErrors* gROC6 = roc(fname6, step, xMin, xMax, cut, 3, leg, leg_name2, var3);

  TGraphErrors* gROC7 = roc(fname7, step, xMin, xMax, cut, 2, leg, leg_name1, var4);
  TGraphErrors* gROC8 = roc(fname8, step, xMin, xMax, cut, 3, leg, leg_name2, var4);

  if(gROC1==0 || gROC2==0) return;
  if(gROC3==0 || gROC4==0) return;
  if(gROC5==0 || gROC6==0) return;
  if(gROC7==0 || gROC8==0) return;

  TGraphErrors* gROC_1 = new TGraphErrors( gROC1->GetN() , gROC1->GetY(), gROC2->GetY(), gROC1->GetEY(), gROC2->GetEY());
  gROC_1->SetLineWidth(3);
  gROC_1->SetLineColor(kRed);

  TGraphErrors* gROC_2 = new TGraphErrors( gROC3->GetN() , gROC3->GetY(), gROC4->GetY(), gROC3->GetEY(), gROC4->GetEY() );
  gROC_2->SetLineWidth(3);
  gROC_2->SetLineColor(kBlue);

  TGraphErrors* gROC_3 = new TGraphErrors( gROC5->GetN() , gROC5->GetY(), gROC6->GetY(), gROC5->GetEY(), gROC6->GetEY());
  gROC_3->SetLineWidth(3);
  gROC_3->SetLineColor(kMagenta);

  TGraphErrors* gROC_4 = new TGraphErrors( gROC7->GetN() , gROC7->GetY(), gROC8->GetY(), gROC7->GetEY(), gROC8->GetEY());
  gROC_4->SetLineWidth(3);
  gROC_4->SetLineColor(kGreen);

  hROC->Draw();
  gROC_1->Draw("SAME");
  gROC_2->Draw("SAME");
  gROC_3->Draw("SAME");
  gROC_4->Draw("SAME");

  leg->Clear();
  leg->AddEntry(gROC_1, leg_name1 , "L");
  leg->AddEntry(gROC_2, leg_name2 , "L");
  leg->AddEntry(gROC_3, leg_name3 , "L");
  leg->AddEntry(gROC_4, leg_name4 , "L");

  TF1* diag = new TF1("diag", "x", 0.,1.);
  diag->SetLineColor(kBlack);
  diag->SetLineStyle(kDashed);
  diag->Draw("SAME");

  leg->Draw();
  c1->SaveAs( save_name );

}



void doall(TString dir="./" ){

  if(false){

    // l+ddddbb vs l+dddddd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Ha2_smear1_btag1_gen15_test6.root", dir+"Test_tH_tL_Hh-Ha1_smear1_btag1_gen11_test6.root",
		  "test__nll[1]-test__nll[0]",
		  0.01, -10., 20.,"",
		  "W+bb discriminator (l+ddddbb vs l+dddddd)",
		  "H_{0} = l+ddddbb", "H_{1} = l+dddddd", 
		  dir+"TopHad_TopLep_Higgs_Ha2vsHa1_ROC.png",  0.1, 1e-03);
    
    // tt+X vs l+ddddbb
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hb1_smear1_btag1_gen12_test6.root", dir+"Test_tH_tL_Hh-Ha2_smear1_btag1_gen15_test6.root",
		  "test__nll[2]-test__nll[1]",
		  0.01, -10., 20.,"",
		  "tt+X discriminator (tt+dd vs l+ddddbb)",
		  "H_{0} = tt+X", "H_{1} = l+ddddbb", 
		  dir+"TopHad_TopLep_Higgs_Hb1vsHa2_ROC.png", 0.1, 1e-03);
    
    // tt+X vs l+dddddd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hb1_smear1_btag1_gen12_test6.root", dir+"Test_tH_tL_Hh-Ha1_smear1_btag1_gen11_test6.root",
		  "test__nll[2]-test__nll[0]",
		  0.01, -10., 20.,"",
		  "tt+X discriminator (tt+dd vs l+dddddd)",
		  "H_{0} = tt+X", "H_{1} = l+dddddd", 
		  dir+"TopHad_TopLep_Higgs_Hb1vsHa1_ROC.png",  0.1, 5e-04);
    
    // tt+bb vs tt+dd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hb3_smear1_btag1_gen7_test6.root", dir+"Test_tH_tL_Hh-Hb1_smear1_btag1_gen12_test6.root",
		  "test__nll[5]-test__nll[3]",
		  0.01, -10., 20.,"",
		  "tt+bb discriminator (tt+bb vs tt+dd)",
		  "H_{0} = tt+bb", "H_{1} = tt+dd", 
		  dir+"TopHad_TopLep_Higgs_Hb3vsHb1_ROC.png", 0.2, 5e-04,1);
    
    // tt+cc vs tt+dd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hb2_smear1_btag1_gen16_test6.root", dir+"Test_tH_tL_Hh-Hb1_smear1_btag1_gen12_test6.root",
		  "test__nll[4]-test__nll[3]",
		  0.01, -10., 20.,"",
		  "tt+cc discriminator (tt+cc vs tt+dd)",
		  "H_{0} = tt+cc", "H_{1} = tt+dd", 
		  dir+"TopHad_TopLep_Higgs_Hb2vsHb1_ROC.png", 0.1, 1e-03);
    
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(false){

    // tt+H vs l+dddddd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hs_smear1_btag1_gen2_test6.root", dir+"Test_tH_tL_Hh-Ha1_smear1_btag1_gen11_test6.root",
		  "test__nll[6]-test__nll[0]",		  
		  0.01, -10., 100.,"",
		  "tt+H discriminator (tt+H vs l+dddddd)",
		  "H_{0} = tt+H", "H_{1} = l+dddddd", 
		  dir+"TopHad_TopLep_Higgs_HsvsHa1_ROC.png", 0.9, 1e-03,1);
    
    // tt+H vs l+ddddbb
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hs_smear1_btag1_gen2_test6.root", dir+"Test_tH_tL_Hh-Ha2_smear1_btag1_gen15_test6.root",
		  "test__nll[6]-test__nll[1]",
		  0.01, -10., 30.,"",
		  "tt+H discriminator (tt+H vs l+ddddbb)",
		  "H_{0} = tt+H", "H_{1} = l+ddddbb", 
		  dir+"TopHad_TopLep_Higgs_HsvsHa2_ROC.png", 0.6, 1e-03,1);
    
    // tt+H vs tt+dd
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hs_smear1_btag1_gen2_test6.root", dir+"Test_tH_tL_Hh-Hb1_smear1_btag1_gen12_test6.root",
		  "test__nll[6]-test__nll[3]",
		  0.01, -10., 20.,"",
		  "tt+H discriminator (tt+H vs tt+dd)",
		  "H_{0} = tt+H", "H_{1} = tt+dd", 
		  dir+"TopHad_TopLep_Higgs_HsvsHb1_ROC.png", 0.2, 5e-04,1);

    // tt+H vs tt+cc
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hs_smear1_btag1_gen2_test6.root", dir+"Test_tH_tL_Hh-Hb2_smear1_btag1_gen16_test6.root",
		  "test__nll[6]-test__nll[4]",
		  0.01, -10., 20.,"",
		  "tt+H discriminator (tt+H vs tt+cc)",
		  "H_{0} = tt+H", "H_{1} = tt+cc", 
		  dir+"TopHad_TopLep_Higgs_HsvsHb2_ROC.png", 1e-02, 0.5e-03,0);
    
    // tt+H vs tt+bb
    roc_comp_ROC( dir+"Test_tH_tL_Hh-Hs_smear1_btag1_gen2_test6.root", dir+"Test_tH_tL_Hh-Hb3_smear1_btag1_gen7_test6.root",
		  "test__nll[6]-test__nll[5]",
		  0.01, -10., 20.,"",
		  "tt+H discriminator (tt+H vs tt+bb)",
		  "H_{0} = tt+H", "H_{1} = tt+bb", 
		  dir+"TopHad_TopLep_Higgs_HsvsHb3_ROC.png", 1e-02, 0.5e-02,0);

    // tt+H (full) bs tt+H (lost)
    roc_comp_ROC( dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		  "test__nll[5]-test__nll[3]",
		  0.01, -10., 20.,"",
		  "tt+H discriminator (tt+H full vs tt+H lost)",
		  "H_{0} = tt+H full", "H_{1} = tt+bb full", 
		  dir+"TopHad_TopLep_Higgs_Hs2vsHb3_ROC.png", 1e-02, 0.5e-02,0);
    
  }

  if(true){  
    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[5]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[5]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[5]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[5]-test__nll[3]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(full) vs tt+bb-(full)",
		    "tt+H-(lost) vs tt+bb-(full)", 
		    "tt+H-(lost) vs tt+bb-(lost)", 
		    "tt+H-(full) vs tt+bb-(full)", 
		    "tt+H-(full) vs tt+bb-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsFull-vs-Hb3Full_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[4]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[4]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[4]-test__nll[3]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[4]-test__nll[3]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(lost) vs tt+bb-(full)",
		    "tt+H-(lost) vs tt+bb-(full)", 
		    "tt+H-(lost) vs tt+bb-(lost)", 
		    "tt+H-(full) vs tt+bb-(full)", 
		    "tt+H-(full) vs tt+bb-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsLost-vs-Hb3Full_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[5]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[5]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[5]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[5]-test__nll[2]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(full) vs tt+bb-(lost)",
		    "tt+H-(lost) vs tt+bb-(full)", 
		    "tt+H-(lost) vs tt+bb-(lost)", 
		    "tt+H-(full) vs tt+bb-(full)", 
		    "tt+H-(full) vs tt+bb-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsFull-vs-Hb3Lost_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[4]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[4]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_bb-FullReco-6jets_smear1_btag1_gen7_test8.root",
		    "test__nll[4]-test__nll[2]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_bb-NotReco-6jets_smear1_btag1_gen19_test8.root",
		    "test__nll[4]-test__nll[2]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(lost) vs tt+bb-(lost)",
		    "tt+H-(lost) vs tt+bb-(full)", 
		    "tt+H-(lost) vs tt+bb-(lost)", 
		    "tt+H-(full) vs tt+bb-(full)", 
		    "tt+H-(full) vs tt+bb-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsLost-vs-Hb3Lost_ROC4.png");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[5]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[5]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[5]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[5]-test__nll[1]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(full) vs tt+dd-(full)",
		    "tt+H-(lost) vs tt+dd-(full)", 
		    "tt+H-(lost) vs tt+dd-(lost)", 
		    "tt+H-(full) vs tt+dd-(full)", 
		    "tt+H-(full) vs tt+dd-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsFull-vs-Hb1Full_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[4]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[4]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[4]-test__nll[1]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[4]-test__nll[1]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(lost) vs tt+dd-(full)",
		    "tt+H-(lost) vs tt+dd-(full)", 
		    "tt+H-(lost) vs tt+dd-(lost)", 
		    "tt+H-(full) vs tt+dd-(full)", 
		    "tt+H-(full) vs tt+dd-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsLost-vs-Hb1Full_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[5]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[5]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[5]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[5]-test__nll[0]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(full) vs tt+dd-(lost)",
		    "tt+H-(lost) vs tt+dd-(full)", 
		    "tt+H-(lost) vs tt+dd-(lost)", 
		    "tt+H-(full) vs tt+dd-(full)", 
		    "tt+H-(full) vs tt+dd-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsFull-vs-Hb1Lost_ROC4.png");

    roc_comp_ROC_4( dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[4]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-NotReco-6jets_smear1_btag1_gen17_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[4]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root", dir+"Test_tH_tL_dd-FullReco-6jets_smear1_btag1_gen12_test8.root",
		    "test__nll[4]-test__nll[0]",
		    dir+"Test_tH_tL_Hh-FullReco-6jets_smear1_btag1_gen2_test8.root",  dir+"Test_tH_tL_dd-NotReco-6jets_smear1_btag1_gen18_test8.root",
		    "test__nll[4]-test__nll[0]",
		    0.001, -10., 10., "", "Hypotheses: tt+H-(lost) vs tt+dd-(lost)",
		    "tt+H-(lost) vs tt+dd-(full)", 
		    "tt+H-(lost) vs tt+dd-(lost)", 
		    "tt+H-(full) vs tt+dd-(full)", 
		    "tt+H-(full) vs tt+dd-(lost)", 
		    dir+"TopHad_TopLep_Higgs_HsLost-vs-Hb1Lost_ROC4.png");

  }



  /*
  roc_comp_ROC_2( dir+"", dir+"",
		  , ,
		  dir+"",  dir+"",
		  , ,
		  0.001, -10., 10., "",
		  "",
		  "", "", 
		  dir+"TopHad_TopLep_Higgs__ROC.png");
  */
}


void get_eff(TString fname="Test_tH_tL_Hh_smear1_btag1_gen2_test2.root",
	     size_t hyp = 0){

  TFile* f = TFile::Open( fname , "READ" );
  if(f==0 || f->IsZombie()) return;

  TTree* tree = (TTree*)f->Get("tree");
  if(tree==0) return;

  double btag[99];
  int    dim [99]; // dim[0] = # of dimensions for hyp '0'
  tree->SetBranchAddress("test__obs_btag", btag);
  tree->SetBranchAddress("test__dim",      dim);

  size_t count_den{0};
  size_t count_num_1{0};
  size_t count_num_2{0};
  size_t count_num_3{0};
  size_t count_num_4{0};

  for(unsigned int i = 0 ; i < tree->GetEntries() ; ++i){
    tree->GetEntry(i);
    
    size_t pass{0};
    for(size_t n_jet = 0 ; n_jet < size_t(dim[hyp]) ; ++n_jet){
      //cout << n_jet << "th jet has " <<  btag[n_jet] << endl;
      if( btag[n_jet]>0.5 ) ++pass;
    }

    ++count_den;
    if( pass>=1 ) ++count_num_1;
    if( pass>=2 ) ++count_num_2;
    if( pass>=3 ) ++count_num_3;
    if( pass>=4 ) ++count_num_4;
  }

  cout << "Total: " << count_den << endl;
  cout << "Pass (>=1): " << float(count_num_1)/count_den << "%" << endl;
  cout << "Pass (>=2): " << float(count_num_2)/count_den << "%" << endl;
  cout << "Pass (>=3): " << float(count_num_3)/count_den << "%" << endl;
  cout << "Pass (>=4): " << float(count_num_4)/count_den << "%" << endl;

  f->Close();
}


void get_eff_W(TString fname="Test_tH_tL_Hh-FullReco_smear1_btag1_gen2_test7.root", float mLow=60., float mHigh=100.){

  TFile* f = TFile::Open( fname , "READ" );
  if(f==0 || f->IsZombie()) return;

  TTree* tree = (TTree*)f->Get("tree");
  if(tree==0) return;

  double btag[99];
  double pt  [99];
  double eta [99];
  double phi [99];
  int    dim [99]; // dim[0] = # of dimensions for hyp '0'

  tree->SetBranchAddress("test__obs_btag", btag);  
  tree->SetBranchAddress("test__obs_pt" , pt);
  tree->SetBranchAddress("test__obs_eta", eta);
  tree->SetBranchAddress("test__obs_phi", phi);
  tree->SetBranchAddress("test__dim",     dim);

  size_t count_den{0};
  size_t count_num_1{0};

  for(unsigned int i = 0 ; i < tree->GetEntries() ; ++i){
    tree->GetEntry(i);
    
    vector<TLorentzVector> vec;
    for(size_t n_jet = 0 ; n_jet < 6 ; ++n_jet){
      TLorentzVector v;
      if( pt[n_jet]>30. )
	v.SetPtEtaPhiM( pt[n_jet], eta[n_jet], phi[n_jet], 0.);
      if( btag[n_jet]>0.5 ) continue;
      vec.push_back( v );
    }

    ++count_den;   
    if(vec.size()<2) continue;

    int pass{0};

    for( size_t k = 0; k < vec.size()-1 && !pass; ++k){
      for( size_t j = k+1; j < vec.size() && !pass; ++j){
	double mass = (vec[k]+vec[j]).M();
	if(mass>mLow && mass<mHigh){
	  ++count_num_1; 
	  pass++;
	}
      }
    }

  }

  cout << "Total: " << count_den << endl;
  cout << "Pass: " << float(count_num_1)/count_den << endl;

  f->Close();
}
