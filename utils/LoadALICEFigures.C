#ifndef LoadALICEFigures_C
#define LoadALICEFigures_C
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLegend.h>


//methods for Preliminary Figure
void SetStyle(Bool_t graypalette=kFALSE);
void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, TString stat = "#int #it{L}d#it{t} =  70 nb^{-1}", Int_t font=42, TString system = "pp, #sqrt{#it{s}} = 8 TeV") ;
void DrawDmesonInJet(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines=kFALSE, Int_t font=42);
void DrawJet(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines=kFALSE, Int_t font=42);
void ZDef(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax,Int_t font=42);
void LogoPythia(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax,TString tune, Int_t font=42);
void DrawJetPart(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines=kFALSE, Int_t font=42);


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, TString stat, Int_t font, TString system) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"
  
  TPaveText *   tex = new TPaveText(xmin,ymin, xmax, ymax, "NDC");
  tex->SetFillStyle(0);
  tex->SetBorderSize(0);
  TString suffix= !logo ? "" : "Preliminary";
  if(logo==2) suffix="Performance";
  if(logo==3) suffix="Simulation";
  
  TString fulltext = TString::Format("ALICE %s",suffix.Data());
  if(!system.IsNull()) fulltext+=system;
  
  if(!stat.IsNull()) fulltext+=stat;
  tex->AddText(fulltext);
  tex->SetTextFont(font);
  
  tex->Draw();
  //TPaveText * text = new TPaveText (xmin, ymin, xmax, ymax, "NDC");
  
  //TPaveText * textL = new TPaveText (xmin,ymin, xmax, ymax, "NDC");
  
  
   
}

void DrawLogo(Int_t logo, Double_t xmin, Double_t ymin, TString stat, Int_t font, TString system) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"
  TLatex logotxt;
  logotxt.SetTextSize(0.07);
  logotxt.SetTextFont(font);
  
  TString suffix= !logo ? "" : "Preliminary";
  if(logo==2) suffix="Performance";
  if(logo==3) suffix="Simulation";
  //tex->AddText(TString::Format("ALICE %s",suffix.Data()));
  TString fulltext = TString::Format("ALICE %s",suffix.Data());
  if(!system.IsNull()) fulltext+=system;
  
  if(!stat.IsNull()) fulltext+=stat;
  
  logotxt.DrawLatex(xmin, ymin, fulltext);

}

void DrawJet(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines, Int_t font){
  TPaveText *textJet=new TPaveText(xmin, ymin,xmax, ymax, "NDC" );
  textJet->SetFillStyle(0);
  textJet->SetBorderSize(0);
  if(splitlines) {
     textJet->AddText("Anti-#it{k}_{T}, #it{R} = 0.4,");
     textJet->AddText("#it{p}_{T,track} > 0.15 GeV/#it{c},");
     textJet->AddText("#it{p}_{T,jet}^{ch,det} > 10 GeV/#it{c},");
     textJet->AddText("|#eta^{jet}|<0.5");
 }
  else textJet->AddText("#splitline{Anti-#it{k}_{T}, #it{R} = 0.4, |#eta^{jet}|<0.5,}{ #it{p}_{T,track} > 0.15 GeV/#it{c}, #it{p}_{T,jet}^{ch,det} > 10 GeV/#it{c}}");
  textJet->SetTextFont(font);
  textJet->Draw();

}

void DrawDmesonInJet(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines,Int_t font){
   TPaveText *textD=new TPaveText(xmin, ymin, xmax, ymax, "NDC");
   textD->SetFillStyle(0);
   textD->SetBorderSize(0);
   if(splitlines) {
      textD->AddText("D^{*+} #rightarrow D^{0}#pi^{+} #rightarrow K^{-}#pi^{+}#pi^{+}");
      textD->AddText("and charge conjugate,");

      textD->AddText("#it{p}_{T}(D^{*+}) > 500 MeV/#it{c},");

      textD->AddText("#it{R}(D^{*+}, jet) < 0.4");

   }
   else textD->AddText("#splitline{D^{*+} #rightarrow D^{0}#pi^{+} #rightarrow K^{-}#pi^{+}#pi^{+} and charge conjugate,}{#it{p}_{T}(D^{*+}) > 500 MeV/#it{c}, #it{R}(D^{*+}, jet) < 0.4}");
  
   textD->SetTextFont(font);
   textD->Draw();
  
}

void ZDef(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax,Int_t font){
   TPaveText *text=new TPaveText(xmin, ymin, xmax, ymax, "NDC");
   text->SetFillStyle(0);
   text->SetBorderSize(0);

   text->AddText("#it{z}^{obs} = #it{p}_{||,D} / #it{p}_{jet}^{ch,det}");
  
   text->SetTextFont(font);
   text->Draw();

}

void LogoPythia(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax,TString tune, Int_t font){
  TPaveText *   tex = new TPaveText(xmin,ymin, xmax, ymax, "NDC");
  tex->SetFillStyle(0);
  tex->SetBorderSize(0);
  tex->AddText(TString::Format("PYTHIA %s simulation",tune.Data()));
  tex->SetTextFont(font);
  tex->Draw();
}

void DrawJetPart(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Bool_t splitlines, Int_t font){
  TPaveText *textJet=new TPaveText(xmin, ymin,xmax, ymax, "NDC" );
  textJet->SetFillStyle(0);
  textJet->SetBorderSize(0);
  if(splitlines) {
     textJet->AddText("Anti-#it{k}_{T}, #it{R} = 0.4,");
     textJet->AddText("#it{p}_{T,jet}^{ch,part} > 10 GeV/#it{c},");
     textJet->AddText("|#eta^{part}|<0.8");
 }
  else textJet->AddText("#splitline{Anti-#it{k}_{T}, #it{R} = 0.4, |#eta^{part}|<0.8,}{#it{p}_{T,jet}^{ch,part} > 10 GeV/#it{c}}");
  textJet->SetTextFont(font);
  textJet->Draw();

}

void SetStyle(Bool_t graypalette) {
  Printf("Setting style!");
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  //gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

#endif

