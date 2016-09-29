#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>

void SetPlotStyle();
void style();
//void SetStyle(Bool_t graypalette);
void SetGrayScale();


void style()
{

/*
  static  int      myDarkRed     = TColor::GetColor(128,0,0);
  static  int      myDarkGreen   = TColor::GetColor(0,128,0);
  static  int      myDarkBlue    = TColor::GetColor(0,0,128);
*/
// Default white background for all plots
 gStyle->SetCanvasColor(10); // Witte achtergrond
 gStyle->SetPalette(1,0); // Fatsoenlijke kleuren
 gStyle->SetStatStyle(kBlack); // Statspad doorzichtig
 //gStyle->SetTitleColor(10);
 gStyle->SetPadColor(10);
 gStyle->SetMarkerStyle(20);
 gStyle->SetMarkerSize(1.2);
 gStyle->SetStatColor(kWhite); // Witte achtergrond stats
 gStyle->SetTitleFillColor(10);

// No canvas or pad borders in produced .eps files
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetPadBorderMode(0);

 gStyle->SetPadGridX(kFALSE);
 gStyle->SetPadGridY(kFALSE);

 gROOT->ForceStyle(); 

 gStyle->SetTextFont(4*10+2);
// gStyle->SetTextSize(gStyle->GetTextSize()*2.);
 gStyle->SetLabelSize(gStyle->GetLabelSize()*2.);

 gROOT->ForceStyle();

// SetPlotStyle();

 //nice size canvas >>>  TCanvas *c = new TCanvas("c","c",550,410)
}

/*
void SetStyle(Bool_t graypalette) {
  Printf("Setting style!");
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(1);//2
  gStyle->SetLabelSize(0.055,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  //  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");//0.2
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(30);//26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gROOT->ForceStyle();
}
*/
void SetPlotStyle() {
  const Int_t nRGBs = 5;
  const Int_t nCont = 255;
  
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}

void SetGrayScale() {
  //
  //source: http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=13395
  //
  
  // TCanvas *c  = new TCanvas("c","Contours",600,0,600,600);
  // TF2 *f1 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",1,3,1,3);  
  
  const UInt_t Number = 2;
  // Double_t Red[Number]   = { 0.1, 0.9};//1.00};
  // Double_t Green[Number] = { 0.1, 0.9};//1.00};
  // Double_t Blue[Number]  = { 0.1, 0.9};//1.00};
  // Double_t Stops[Number] = { 0.1, 0.9};//1.00};

  Double_t Red[Number]   = {1.,0.05};
  Double_t Green[Number] = {1.,0.05};
  Double_t Blue[Number]  = {1.,0.05};
  Double_t Stops[Number] = {0.,1.};
  
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
  //f2->SetContour(nb);
  //f2->Draw("surf1z");

}

