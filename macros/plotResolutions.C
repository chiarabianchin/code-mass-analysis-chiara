#include<THnSparse.h>
#include<TH1D.h>
#include<TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLegend.h>

//axes
//0 = (Mdet-Mpart), 1 = (Mdet-Mpart)/Mpart, 2 = Mpart, 3 = pTdet, 4 = pTpart
//list of functions
void SaveCv(TCanvas* c,TString suffix="", Int_t format=2);

TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Int_t *pTPBin, Int_t axisP, char* hprjname);

TPaveText* PrintInfoOnCanvas(TH1D* h, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString text, TVirtualPad* cvdest);

void DefineTrendingHistos(TString nameX, TString nameY, TString titleX, TString titleY, Int_t nbins, Double_t min, Double_t max, TH1D*& hMean, TH1D*& hWidth );

void FillTrendings(TH1D *hMeanDMvspT, TH1D* hWidthDMvspT, Int_t ipT,TH1D* hDM);

void ResetFullRanges(THnSparse *hns);

void CanvasDim(Int_t nptbins, Int_t& px, Int_t& py, Int_t& dimx, Int_t& dimy);

void SetStyle();
   
TList* ReadFile(TString strIn, TString strLst);

//main
void plotResolutions(Int_t axproj=0,TString strIn="JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC.root", TString strLst = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC") {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleXSize(0.04);
  
  TFile *f = new TFile(strIn.Data());
  if(!f->IsOpen()){
     Printf("File %s not found", strIn.Data());
     return;
  }
  TList *lst = static_cast<TList*>(f->Get(strLst.Data()));

  THnSparse *fhnDeltaMass = (THnSparse*)lst->FindObject("fhnDeltaMass");
  THnSparse *fhnDeltaMassCorr = (THnSparse*)lst->FindObject("fhnDeltaMassCorr");

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //pTpart bins

  const Int_t npTbins = 6;
  Float_t pTinterval=20., pTMin=20., pTMax=pTMin+npTbins*pTinterval; //GeV/c
  
  TCanvas *cDMvspTP=new TCanvas(Form("cDM%dvspTP",axproj), "Delta M vs pTPart", 1000,1000);
  cDMvspTP->Divide(3,2);
  Int_t axselrange=4;
  TH1D *hMeanDMvspTP;
  TH1D *hWidthDMvspTP;
  DefineTrendingHistos("pTP", Form("DM%d",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMvspTP, hWidthDMvspTP);
  if(!hMeanDMvspTP || !hWidthDMvspTP){
     Printf("Error in creating trending histograms");
     return;
  }
  TH1D *hMeanDMCorrvspTP;
  TH1D *hWidthDMCorrvspTP;
  DefineTrendingHistos("pTP", Form("DM%dCorr",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMCorrvspTP, hWidthDMCorrvspTP);
  
  if(!hMeanDMCorrvspTP || !hWidthDMCorrvspTP){
     Printf("Error in creating trending histograms");
     return;
  }
  hMeanDMCorrvspTP->SetMarkerColor(kRed+1);
  hWidthDMCorrvspTP->SetMarkerColor(kRed+1);
  for(Int_t ipT=0;ipT<npTbins;ipT++){ //loop on pT det bins
     TString ptstring=Form("%.0f<p_{T,jet}<%.0f GeV/c",pTMin, pTMin+pTinterval);
     
     TPaveText *pv=new TPaveText(0.1,0.4,0.45,0.55,"NDC");
     pv->SetFillStyle(0);
     pv->SetBorderSize(0);
     pv->AddText(ptstring);

     Int_t pTPBin[2] = {fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin), fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin+pTinterval)-1};
     Printf("Projections in pT bin range %d-%d",pTPBin[0],pTPBin[1]);
     
     TH1D* hDMpTP = (TH1D*)DoProjections(fhnDeltaMass, axselrange, pTPBin, axproj, Form("hDM%dpTP%d",axproj,ipT));
     
     TH1D* hDMCorrpTP = DoProjections(fhnDeltaMassCorr, axselrange, pTPBin, axproj, Form("hDM%dCorrpTP%d",axproj,ipT));
     hDMCorrpTP->SetLineStyle(2);

     cDMvspTP->cd(ipT+1);
     //gPad->SetLogy();
     hDMpTP->Draw();
     hDMCorrpTP->Draw("sames");
     pv->Draw();
     PrintInfoOnCanvas(hDMpTP, 0.15,0.25,0.9,0.4, "(det)", cDMvspTP->cd(ipT+1));
     PrintInfoOnCanvas(hDMCorrpTP, 0.15,0.10,0.9,0.35, "(det corr)", cDMvspTP->cd(ipT+1));
     
     FillTrendings(hMeanDMvspTP, hWidthDMvspTP, ipT, hDMpTP);
     FillTrendings(hMeanDMCorrvspTP, hWidthDMCorrvspTP, ipT, hDMCorrpTP);
     //     Int_t pTDBin[2] = {fhnDeltaMass->GetAxis(2)->FindBin(pTMin), fhnDeltaMass->GetAxis(2)->FindBin(pTMin+pTinterval)-1};
     
     pTMin+=pTinterval;
  }
  SaveCv(cDMvspTP);
  TCanvas *cSummaryDMpTP = new TCanvas(Form("cSummaryDM%dpTP",axproj), "Summary Delta Mass vs pT part", 500,900);
  cSummaryDMpTP->Divide(1,2);
  cSummaryDMpTP->cd(1);
  hMeanDMvspTP->SetMaximum(0);
  hMeanDMvspTP->Draw("P");
  hMeanDMCorrvspTP->Draw("Psames");
  cSummaryDMpTP->cd(2);
  hWidthDMvspTP->Draw("P");
  hWidthDMCorrvspTP->Draw("Psames");

  SaveCv(cSummaryDMpTP);
  
  ResetFullRanges(fhnDeltaMass);
  ResetFullRanges(fhnDeltaMassCorr);
 
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
  // pTdet bins (in bins of pT part to access feed-in)
  
  Int_t axselrange2=4; //pTpart
  
  const Int_t npTPbins=6;
  Double_t pTPinterval=20., pTPMin=0.; 
  Double_t pTPMax=pTPMin+npTPbins*pTPinterval; //GeV/c
  
  pTinterval=20.; pTMin=0.;  pTMax=pTMin+npTbins*pTinterval; //GeV/c
  Double_t pTMinCp=pTMin;
  
  TCanvas *cDMvspTDManyLines=new TCanvas(Form("cDMvspTDpTP"), Form("Delta M vs pTDet"), 1000,1000);
  
  // pTdet inclusive for comparison with pTdet in bins of pTpart
  Bool_t runInclusive=kFALSE;
  if(runInclusive) Printf("############### Run inclusive in pTPart ##########\n All canvas for different pTPart are the same\n Not sure if the low/high Feed out canvas make sense");
  
  //these are filled with the mean and widht of the delta mass in the case pTD bin == pTP bin
  TH1D *hMeanDMvspTDpTP;
  TH1D *hWidthDMvspTDpTP;
  DefineTrendingHistos(Form("pTDEqpTP"), Form("DM%d",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMvspTDpTP, hWidthDMvspTDpTP);
  if(!hMeanDMvspTDpTP || !hWidthDMvspTDpTP){
     Printf("Error in creating trending histograms");
     return;
  }
  TH1D *hMeanDMCorrvspTDpTP;
  TH1D *hWidthDMCorrvspTDpTP;
  DefineTrendingHistos(Form("pTDEqpTP"), Form("DM%dCorr",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMCorrvspTDpTP, hWidthDMCorrvspTDpTP);
  
  if(!hMeanDMCorrvspTDpTP || !hWidthDMCorrvspTDpTP){
     Printf("Error in creating trending histograms");
     return;
  }
  hMeanDMCorrvspTDpTP->SetMarkerColor(kRed+1);
  hWidthDMCorrvspTDpTP->SetMarkerColor(kRed+1);

  //ratios
  TH1D* hintegralDMrOrigSameCorrSame   = new TH1D("hintegralDMrOrigSameCorrSame","hintegralDMrOrigSameCorrSame; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //stay in the right bin
  //TH1D* hintegralDMrOrigLowerCorrSame  = new TH1D("hintegralDMrOrigLowerCorrSame","hintegralDMrOrigLowerCorrSame; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //go to the right bin
  //TH1D* hintegralDMrOrigHigherCorrSame = new TH1D("hintegralDMrOrigHigherCorrSame","hintegralDMrOrigHigherCorrSame; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //go to the right bin
  //TH1D* hintegralDMrOrigLowerCorrLower = new TH1D("hintegralDMrOrigLowerCorrLower","hintegralDMrOrigLowerCorrLower; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //stay lower  
  //TH1D* hintegralDMrOrigLowerCorrHigher= new TH1D("hintegralDMrOrigLowerCorrHigher","hintegralDMrOrigLowerCorrHigher; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //from lower to higher
  //TH1D* hintegralDMrOrigHigherCorrHigher= new TH1D("hintegralDMrOrigHigherCorrHigher","hintegralDMrOrigHigherCorrHigher; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);//stay higher  
  //TH1D* hintegralDMrOrigHigherCorrLower= new TH1D("hintegralDMrOrigHigherCorrLower","hintegralDMrOrigHigherCorrLower; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  //from higher to lower  
  //TH1D* hintegralDMrOrigSameCorrLower  = new TH1D("hintegralDMrOrigSameCorrLower","hintegralDMrOrigSameCorrLower; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax);  // from right to lower
  //TH1D* hintegralDMrOrigSameCorrHigher = new TH1D("hintegralDMrOrigSameCorrHigher","hintegralDMrOrigSameCorrHigher; #it{p}_{T,jet}^{Part} (GeV/#it{c}); Ratio Integrals",npTPbins,pTPMin, pTPMax); //from right to higher

  for(Int_t j=0;j<npTPbins;j++){
     //select pTpart bin and loop over all the pTdet bins
     
     Int_t pTPBin[2] = {fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPMin), fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPMin+pTPinterval)-1};
     
     TString ptpstring=Form("%.0f<p_{T,jet}^{Part}<%.0f GeV/c",pTPMin, pTPMin+pTPinterval);
     
     TPaveText *pvP=new TPaveText(0.1,0.6,0.45,0.75,"NDC");
     pvP->SetFillStyle(0);
     pvP->SetBorderSize(0);
     pvP->SetTextColor(kGreen+3);
     pvP->AddText(ptpstring);
     
     if(!runInclusive){
     	fhnDeltaMass->GetAxis(axselrange2)->SetRange(pTPBin[0], pTPBin[1]);
     	fhnDeltaMassCorr->GetAxis(axselrange2)->SetRange(pTPBin[0], pTPBin[1]);
     }
     TCanvas *cDMvspTD=new TCanvas(Form("cDM%dvspTDpTP%d",axproj,j), Form("Delta M vs pTDet (%.0f<pTPart<%.0f GeV/c)",pTPMin,pTPMin+pTPinterval), 1000,1000);
     cDMvspTD->Divide(3,2);
     TCanvas *cDMpTPartFeedOut = new  TCanvas(Form("cDM%dpTPartFeedOut%d",axproj,j), Form("Delta M in pTDetBin = pTPartBin (%.0f<pTPart<%.0f GeV/c) and outside",pTPMin,pTPMin+pTPinterval), 1000,500);
     cDMpTPartFeedOut->Divide(2,1);
     TLegend *legFeedOut = new TLegend(0.15,0.5,0.5,0.85);
     legFeedOut->SetFillStyle(0);
     legFeedOut->SetBorderSize(0);

     axselrange=3; //pTDet
          
     TH1D *hMeanDMvspTD;
     TH1D *hWidthDMvspTD;
     DefineTrendingHistos(Form("pTD-pTP%d",j), Form("DM%d",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMvspTD, hWidthDMvspTD);
     if(!hMeanDMvspTD || !hWidthDMvspTD){
     	Printf("Error in creating trending histograms");
     	return;
     }
     TH1D *hMeanDMCorrvspTD;
     TH1D *hWidthDMCorrvspTD;
     DefineTrendingHistos(Form("pTD-pTP%d",j), Form("DM%dCorr",axproj), "#it{p}_{T,jet} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",npTbins, pTMin, pTMax, hMeanDMCorrvspTD, hWidthDMCorrvspTD);
     
     if(!hMeanDMCorrvspTD || !hWidthDMCorrvspTD){
     	Printf("Error in creating trending histograms");
     	return;
     }
     hMeanDMCorrvspTD->SetMarkerColor(kRed+1);
     hWidthDMCorrvspTD->SetMarkerColor(kRed+1);
     
     //Delta Mass for the feed out bins
     TH1D* hDMFeedOut=0x0;
     TH1D* hDMFeedOutLow=0x0;
     TH1D* hDMFeedOutHigh=0x0;
     TH1D* hDMCorrFeedOut=0x0;
     TH1D* hDMCorrFeedOutLow=0x0;
     TH1D* hDMCorrFeedOutHigh=0x0;
     TH1D* hDMNoFeedOut=0x0;
     TH1D* hDMCorrNoFeedOut=0x0;
     
     //ratios
     TH1D* hDMrOrigSameCorrSame   = 0x0;  //stay in the right bin
     TH1D* hDMrOrigLowerCorrSame  = 0x0;  //go to the right bin
     TH1D* hDMrOrigHigherCorrSame = 0x0;  //go to the right bin
     TH1D* hDMrOrigLowerCorrLower = 0x0;  //stay lower  
     TH1D* hDMrOrigLowerCorrHigher= 0x0;  //from lower to higher
     TH1D* hDMrOrigHigherCorrHigher = 0x0;//stay higher  
     TH1D* hDMrOrigHigherCorrLower= 0x0;  //from higher to lower  
     TH1D* hDMrOrigSameCorrLower  = 0x0;  // from right to lower
     TH1D* hDMrOrigSameCorrHigher = 0x0; //from right to higher
     
     Double_t totalhDM, totalhDMCorr;
     
     for(Int_t ipT=0;ipT<npTbins;ipT++){ //loop on pT det bins
     	TString ptstring=Form("%.0f<p_{T,jet}<%.0f GeV/c",pTMin, pTMin+pTinterval);
     	
     	TPaveText *pv=new TPaveText(0.1,0.4,0.45,0.55,"NDC");
     	pv->SetFillStyle(0);
     	pv->SetBorderSize(0);
     	pv->AddText(ptstring);
     	
     	Int_t pTDBin[2] = {fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin), fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin+pTinterval)-1};
     	Printf("Projections in pT bin range %d-%d",pTDBin[0],pTDBin[1]);
     	
     	TH1D* hDMpTD = (TH1D*)DoProjections(fhnDeltaMass, axselrange, pTDBin, axproj, Form("hDM%dpTD%d",axproj,ipT));
     	totalhDM = hDMpTD->Integral();
     	TH1D* hDMCorrpTD = DoProjections(fhnDeltaMassCorr, axselrange, pTDBin, axproj, Form("hDM%dCorrpTD%d",axproj,ipT));
     	hDMCorrpTD->SetLineStyle(2);
     	totalhDMCorr = hDMCorrpTD->Integral();
     	
     	if(ipT != j){ //pTbin det different from pTbin Part
     	   if(!hDMFeedOut) { //define the histo by cloning
     	      hDMFeedOut = (TH1D*)hDMpTD->Clone("hDMFeedOut");   
     	      hDMCorrFeedOut = (TH1D*)hDMCorrpTD->Clone("hDMCorrFeedOut");
     	      //Printf("%s and %s defined.",hDMFeedOut->GetName(), hDMCorrFeedOut->GetName());
     	   } else { //Add the histos
     	      hDMFeedOut    ->Add(hDMpTD);
     	      hDMCorrFeedOut->Add(hDMCorrpTD);
     	   }
     	   if(ipT<j) {
     	      if(!hDMFeedOutLow) {
     	      	 hDMFeedOutLow = (TH1D*)hDMpTD->Clone("hDMFeedOutLow");
     	      	 hDMCorrFeedOutLow = (TH1D*)hDMCorrpTD->Clone("hDMCorrFeedOutLow");
     	      	 hDMFeedOutLow    ->SetLineColor(kGreen-4);
     	      	 hDMCorrFeedOutLow->SetLineColor(kGreen-4);
     	      	 //Printf("%s and %s defined.",hDMFeedOutLow->GetName(), hDMCorrFeedOutLow->GetName());
     	      	 
     	      } else {
     	      	 hDMFeedOutLow    ->Add(hDMpTD);
     	      	 hDMCorrFeedOutLow->Add(hDMCorrpTD);
     	      }
     	   } else{ //if(ipT>j)
     	      if(!hDMFeedOutHigh) {
     	      	 hDMFeedOutHigh = (TH1D*)hDMpTD->Clone("hDMFeedOutHigh");
     	      	 hDMCorrFeedOutHigh = (TH1D*)hDMCorrpTD->Clone("hDMCorrFeedOutHigh");
     	      	 hDMFeedOutHigh    ->SetLineColor(kRed+1);
     	      	 hDMCorrFeedOutHigh->SetLineColor(kRed+1);
     	      	 //Printf("%s and %s defined.",hDMFeedOutHigh->GetName(), hDMCorrFeedOutHigh->GetName());
     	      } else {
     	      	 hDMFeedOutHigh    ->Add(hDMpTD);
     	      	 hDMCorrFeedOutHigh->Add(hDMCorrpTD);
     	      }
     	      
     	   }
     	   
      	} else {
     	        	   
     	   hDMNoFeedOut = (TH1D*)hDMpTD->Clone("hDMNoFeedOut");
     	   hDMCorrNoFeedOut = (TH1D*)hDMCorrpTD->Clone("hDMCorrNoFeedOut");

     	   cDMpTPartFeedOut->cd(1);

     	   hDMNoFeedOut->GetXaxis()->SetRangeUser(-15,10);
     	   hDMNoFeedOut->Draw();
     	   hDMCorrNoFeedOut->Draw("sames");
     	   pvP->Draw();
     	}
     	
     	cDMvspTD->cd(ipT+1);
     	//gPad->SetLogy();
     	hDMpTD->Draw();
     	hDMCorrpTD->Draw("sames");
     	pv->Draw();
     	PrintInfoOnCanvas(hDMpTD, 0.15,0.25,0.9,0.4, "(det)", cDMvspTD->cd(ipT+1));
     	PrintInfoOnCanvas(hDMCorrpTD, 0.15,0.10,0.9,0.35, "(det corr)", cDMvspTD->cd(ipT+1));
     	if(ipT==j){
     	   if(!runInclusive) pvP->Draw();
     	   FillTrendings(hMeanDMvspTDpTP, hWidthDMvspTDpTP, ipT, hDMpTD);
     	   FillTrendings(hMeanDMCorrvspTDpTP, hWidthDMCorrvspTDpTP, ipT, hDMCorrpTD);
     	   
     	}
     	FillTrendings(hMeanDMvspTD, hWidthDMvspTD, ipT, hDMpTD);
     	FillTrendings(hMeanDMCorrvspTD, hWidthDMCorrvspTD, ipT, hDMCorrpTD);
     	//     Int_t pTDBin[2] = {fhnDeltaMass->GetAxis(2)->FindBin(pTMin), fhnDeltaMass->GetAxis(2)->FindBin(pTMin+pTinterval)-1};
     	
     	pTMin+=pTinterval;
     }
     pTMin=pTMinCp;
     
     SaveCv(cDMvspTD);

     //calculate ratios of feed out
     Printf("Filling the ratio plots");
     //stay in the right bin
     hDMrOrigSameCorrSame   = (TH1D*)hDMNoFeedOut->Clone("hDMrOrigSameCorrSame");
     hDMrOrigSameCorrSame->Divide(hDMCorrNoFeedOut);
     hDMrOrigSameCorrSame->SetLineWidth(2);
     hDMrOrigSameCorrSame->SetLineStyle(1);
     //hDMrOrigSameCorrSame->SetLineColor(kBlue);
     
     Double_t integralhDMrOrigSameCorrSame = hDMNoFeedOut->Integral()/totalhDM;
     hintegralDMrOrigSameCorrSame->SetBinContent(j+1, integralhDMrOrigSameCorrSame);
     Printf("NOW :: bin %d ratio in the right bin over total %f, after correction %f", j, integralhDMrOrigSameCorrSame, hDMCorrNoFeedOut->Integral()/totalhDM);


     Printf("hDMrOrigSameCorrSame");
     //go to the right bin
     if(hDMFeedOutLow) {
     	hDMrOrigLowerCorrSame  = (TH1D*)hDMFeedOutLow->Clone("hDMrOrigLowerCorrSame");
     	hDMrOrigLowerCorrSame->Divide(hDMCorrNoFeedOut);
     	hDMrOrigLowerCorrSame->SetLineWidth(2);
     	hDMrOrigLowerCorrSame->SetLineStyle(1);
     	hDMrOrigLowerCorrSame->SetLineColor(kGreen-4);
     	
     	Printf("hDMrOrigLowerCorrSame");
     	//stay lower
     	hDMrOrigLowerCorrLower =    (TH1D*)hDMFeedOutLow->Clone("hDMrOrigLowerCorrLower");
     	if(hDMCorrFeedOutLow){
     	hDMrOrigLowerCorrLower->Divide(hDMCorrFeedOutLow);
     	hDMrOrigLowerCorrLower->SetLineWidth(2);
     	hDMrOrigLowerCorrLower->SetLineStyle(2);
     	hDMrOrigLowerCorrLower->SetLineColor(kGreen-4);
     	} else hDMrOrigLowerCorrLower = 0x0;
     	Printf("hDMrOrigLowerCorrLower");
     	
     	//from lower to higher
     	hDMrOrigLowerCorrHigher =   (TH1D*)hDMFeedOutLow->Clone("hDMrOrigLowerCorrHigher");
     	if(hDMCorrFeedOutHigh) {
     	   hDMrOrigLowerCorrHigher->Divide(hDMCorrFeedOutHigh);
     	   hDMrOrigLowerCorrHigher->SetLineWidth(2);
     	   hDMrOrigLowerCorrHigher->SetLineStyle(3);
     	   hDMrOrigLowerCorrHigher->SetLineColor(kGreen-4);
     	} else hDMrOrigLowerCorrHigher = 0x0;
     	
     	Printf("hDMrOrigLowerCorrHigher");
     	
     }
     
     if(hDMFeedOutHigh) {
     	hDMrOrigHigherCorrSame = (TH1D*)hDMFeedOutHigh->Clone("hDMrOrigHigherCorrSame");
     	hDMrOrigHigherCorrSame->Divide(hDMCorrNoFeedOut);
     	hDMrOrigHigherCorrSame->SetLineWidth(2);
     	hDMrOrigHigherCorrSame->SetLineStyle(1);
     	hDMrOrigHigherCorrSame->SetLineColor(kRed+1);
     	Printf("hDMrOrigHigherCorrSame");
     	//stay higher
     	hDMrOrigHigherCorrHigher =  (TH1D*)hDMFeedOutHigh->Clone("hDMrOrigHigherCorrHigher");
     	if(hDMCorrFeedOutHigh){
     	   hDMrOrigHigherCorrHigher->Divide(hDMCorrFeedOutHigh);
     	   hDMrOrigHigherCorrHigher->SetLineWidth(2);
     	   hDMrOrigHigherCorrHigher->SetLineStyle(2);
     	   hDMrOrigHigherCorrHigher->SetLineColor(kRed+1);
     	} else hDMrOrigHigherCorrHigher=0x0;
     	Printf("hDMrOrigHigherCorrHigher");
     	//from higher to lower 
     	hDMrOrigHigherCorrLower =   (TH1D*)hDMFeedOutHigh->Clone("hDMrOrigHigherCorrLower"); 
     	hDMrOrigHigherCorrLower->Divide(hDMCorrFeedOutLow);
     	hDMrOrigHigherCorrLower->SetLineWidth(2);
     	hDMrOrigHigherCorrLower->SetLineStyle(3);
     	hDMrOrigHigherCorrLower->SetLineColor(kRed+1);
     	Printf("hDMrOrigHigherCorrLower");
     	
     	
     }

     // from right to lower
     hDMrOrigSameCorrLower = (TH1D*)hDMNoFeedOut->Clone("hDMrOrigSameCorrLower");
     if(hDMCorrFeedOutLow){
     hDMrOrigSameCorrLower->Divide(hDMCorrFeedOutLow);
     hDMrOrigSameCorrLower->SetLineWidth(2);
     hDMrOrigSameCorrLower->SetLineStyle(3);
     hDMrOrigSameCorrLower->SetLineColor(kCyan-5);
     }else hDMrOrigSameCorrLower = 0x0;
     Printf("hDMrOrigSameCorrLower");

     //from right to higher
     hDMrOrigSameCorrHigher = (TH1D*)hDMNoFeedOut->Clone("hDMrOrigSameCorrLower"); 
     hDMrOrigSameCorrHigher->Divide(hDMCorrFeedOutHigh);
     hDMrOrigSameCorrHigher->SetLineWidth(2);
     hDMrOrigSameCorrHigher->SetLineStyle(3);
     hDMrOrigSameCorrHigher->SetLineColor(kCyan-4);    
     Printf("hDMrOrigSameCorrHigher");
     
     TCanvas *cSummaryDMpTD = new TCanvas(Form("cSummaryDM%dpTDpTP%d",axproj,j), Form("Summary Delta Mass vs pT det (%.0f<pTPart<%.0f GeV/c)",pTPMin,pTPMin+pTPinterval), 500,900);
     cSummaryDMpTD->Divide(1,2);
     cSummaryDMpTD->cd(1);
     hMeanDMvspTD->SetMaximum(0);
     hMeanDMvspTD->Draw("P");
     hMeanDMCorrvspTD->Draw("Psames");
     cSummaryDMpTD->cd(2);
     hWidthDMvspTD->Draw("P");
     hWidthDMCorrvspTD->Draw("Psames");
     
     cDMvspTDManyLines->cd(2);
     hWidthDMvspTD->Draw("Psames");
     hWidthDMCorrvspTD->Draw("Psames");
     
     SaveCv(cSummaryDMpTD);
     
      
     cDMpTPartFeedOut->cd(2);
     if(hDMFeedOut)        {
     	hDMFeedOut->GetXaxis()->SetRangeUser(-15,10);
     	hDMFeedOut         ->Draw();
     	legFeedOut->AddEntry(hDMFeedOut, "p_{T,Det} != p_{T,Part} ","l");
     }
     if(hDMFeedOutLow)     {
     	hDMFeedOutLow      ->Draw("sames");
     	legFeedOut->AddEntry(hDMFeedOutLow, "p_{T,Det} < p_{T,Part} ","l");
     }
     if(hDMCorrFeedOutLow) hDMCorrFeedOutLow  ->Draw("sames");
     if(hDMFeedOutHigh)    {
     	hDMFeedOutHigh     ->Draw("sames");
     	legFeedOut->AddEntry(hDMFeedOutHigh, "p_{T,Det} > p_{T,Part} ","l");
     }
     if(hDMCorrFeedOutHigh)hDMCorrFeedOutHigh ->Draw("sames");
     if(hDMCorrFeedOut)    {
     	hDMCorrFeedOut     ->Draw("sames");
     	legFeedOut->AddEntry(hDMCorrFeedOut, "Corrected ","l");
     }
 
     legFeedOut->Draw();
     SaveCv(cDMpTPartFeedOut);
     
     TCanvas *cDMpTPartFeedOutRatios = new TCanvas(Form("cDMpTPartFeedOutRatios%d",j), Form("Ratios feed out (%.0f<p_{T,Part}<%.0f)",pTPMin,pTPMin+pTPinterval), 800,800);
     TLegend *legFeedOutR = new TLegend(0.15,0.5,0.5,0.85);
     legFeedOutR->SetFillStyle(0);
     legFeedOutR->SetBorderSize(0);
     if(hDMrOrigSameCorrSame)    legFeedOutR->AddEntry(hDMrOrigSameCorrSame    ,"Stay in the right bin","l");
     if(hDMrOrigLowerCorrSame)   legFeedOutR->AddEntry(hDMrOrigLowerCorrSame   ,"From lower to the right bin","l");
     if(hDMrOrigHigherCorrSame)  legFeedOutR->AddEntry(hDMrOrigHigherCorrSame  ,"From higher to the right bin","l");
     if(hDMrOrigLowerCorrLower)  legFeedOutR->AddEntry(hDMrOrigLowerCorrLower  ,"Stay lower","l");
     if(hDMrOrigLowerCorrHigher) legFeedOutR->AddEntry(hDMrOrigLowerCorrHigher ,"From lower to higher","l");
     if(hDMrOrigHigherCorrHigher)legFeedOutR->AddEntry(hDMrOrigHigherCorrHigher,"Stay higher","l");
     if(hDMrOrigHigherCorrLower) legFeedOutR->AddEntry(hDMrOrigHigherCorrLower ,"From higher to lower","l");
     if(hDMrOrigSameCorrLower)   legFeedOutR->AddEntry(hDMrOrigSameCorrLower   ,"From right bin to lower","l");
     if(hDMrOrigSameCorrHigher)  legFeedOutR->AddEntry(hDMrOrigSameCorrHigher  ,"From right bin to higher","l");

     cDMpTPartFeedOutRatios->cd();
     if(hDMrOrigSameCorrSame)    hDMrOrigSameCorrSame->Draw("");
     if(hDMrOrigLowerCorrSame)   hDMrOrigLowerCorrSame->Draw("sames");
     if(hDMrOrigHigherCorrSame)  hDMrOrigHigherCorrSame->Draw("sames");
     if(hDMrOrigLowerCorrLower)  hDMrOrigLowerCorrLower->Draw("sames");
     if(hDMrOrigLowerCorrHigher) hDMrOrigLowerCorrHigher->Draw("sames");
     if(hDMrOrigHigherCorrHigher)hDMrOrigHigherCorrHigher->Draw("sames");
     if(hDMrOrigHigherCorrLower) hDMrOrigHigherCorrLower->Draw("sames");
     if(hDMrOrigSameCorrLower)   hDMrOrigSameCorrLower->Draw("sames");
     if(hDMrOrigSameCorrHigher)  hDMrOrigSameCorrHigher->Draw("sames");
     legFeedOutR->Draw();
     //pvP->PaintPave(0.6,0.8,0.9,0.9, 0);
     pvP->Draw();
     
     SaveCv(cDMpTPartFeedOutRatios);
     
     pTPMin+=pTPinterval;
  }
  
  TCanvas *cSummaryDMpTDpTP = new TCanvas(Form("cSummaryDM%dpTDEqpTP",axproj), Form("Summary Delta Mass vs pT det = pTPart"), 500,900);
  cSummaryDMpTDpTP->Divide(1,2);
  cSummaryDMpTDpTP->cd(1);
  //hMeanDMvspTDpTP->SetMaximum(0);
  hMeanDMvspTDpTP->Draw("P");
  hMeanDMCorrvspTDpTP->Draw("Psames");
  cSummaryDMpTDpTP->cd(2);
  hWidthDMvspTDpTP->Draw("P");
  hWidthDMCorrvspTDpTP->Draw("Psames");
  SaveCv(cSummaryDMpTDpTP);
  
  TCanvas *cIntegral = new TCanvas("cIntegral", "Integral Ratios", 800,800);
  cIntegral->cd();
  hintegralDMrOrigSameCorrSame->Draw("P");
  
  ResetFullRanges(fhnDeltaMass);
  ResetFullRanges(fhnDeltaMassCorr);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  // MPart bins
  const Int_t nMbins = 10;
  Double_t MMin = 0, MMinCp=MMin, MMax = 25;
  Double_t Minterval=(MMax-MMin)/(Double_t)nMbins;
  Int_t px,py,dimx,dimy;
  CanvasDim(nMbins,px,py,dimx,dimy);
  
  Int_t npTPartLargeBins=2;
  Double_t pTPartLargeinterval=20., pTPartLargeMin = 40;
  /*
  Double_t  pTPartLargeMax = pTPartLargeinterval * npTPartLargeBins + pTPartLargeMin;
  TCanvas *cMassOutside = new TCanvas("cMassOutside", Form("Mass outside the selected pT ranges <%.0f and > %.0f", pTPartLargeMin, pTPartLargeMax),800,800 );
  Int_t pTPLargeBinTest[2] = {fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPartLargeMin), fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPartLargeMax)};
  fhnDeltaMass->GetAxis(axselrange2)->SetRange(pTPLargeBinTest[0],pTPLargeBinTest[1]);      fhnDeltaMassCorr->GetAxis(axselrange2)->SetRange(pTPLargeBinTest[0],pTPLargeBinTest[1]);

  Int_t MPBinTest[2] = {1, fhnDeltaMass->GetAxis(axselrange)->FindBin(MMin)};

  TH1D* hDMMPTest = (TH1D*)DoProjections(fhnDeltaMass, axselrange, MPBinTest, axproj, Form("hDM%dMPTest",axproj));
  
  TH1D* hDMCorrMPTest = DoProjections(fhnDeltaMassCorr, axselrange, MPBinTest, axproj, Form("hDM%dCorrMPTest",axproj));
  hDMCorrMPTest->SetLineStyle(2);
  cMassOutside->cd();
  hDMMPTest->Draw();
  hDMCorrMPTest->Draw("sames");
  ResetFullRanges(fhnDeltaMass);
  ResetFullRanges(fhnDeltaMassCorr);
  */
  for(Int_t j=0;j<npTPartLargeBins;j++){
     // range of pTPart
     Int_t pTPLargeBin[2] = {fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPartLargeMin), fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPartLargeMin+pTPartLargeinterval)-1};
     
     Printf("Projections in pTpart bin range %d-%d (%.0f-%.0f GeV/c)",pTPLargeBin[0],pTPLargeBin[1], pTPartLargeMin,pTPartLargeMin+pTPartLargeinterval);
     fhnDeltaMass->GetAxis(axselrange2)->SetRange(pTPLargeBin[0],pTPLargeBin[1]);      fhnDeltaMassCorr->GetAxis(axselrange2)->SetRange(pTPLargeBin[0],pTPLargeBin[1]); 
     
     TCanvas *cDMvsMP=new TCanvas(Form("cDM%dvsMP%d",axproj,j), Form("Delta M vs MPart %.0f<p_{T,Part}<%.0f GeV/c", pTPartLargeMin, pTPartLargeMin+pTPartLargeinterval), dimx,dimy);
     cDMvsMP->Divide(py,px);
     axselrange=2; //Mpart
     TH1D *hMeanDMvsMP;
     TH1D *hWidthDMvsMP;
     DefineTrendingHistos(Form("MP%d",j), Form("DM%d%d",axproj,j), "#it{M}_{part} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",nMbins, MMinCp, MMax, hMeanDMvsMP, hWidthDMvsMP);
     if(!hMeanDMvsMP || !hWidthDMvsMP){
     	Printf("Error in creating trending histograms");
     	return;
     }
     TH1D *hMeanDMCorrvsMP;
     TH1D *hWidthDMCorrvsMP;
     DefineTrendingHistos(Form("MP%d",j), Form("DM%dCorr%d",axproj,j), "#it{M}_{part} (GeV/#it{c})", "Delta Mass (GeV/#it{c^{2}})",nMbins, MMinCp, MMax, hMeanDMCorrvsMP, hWidthDMCorrvsMP);
     
     if(!hMeanDMCorrvsMP || !hWidthDMCorrvsMP){
     	Printf("Error in creating trending histograms");
     	return;
     }
     hMeanDMCorrvsMP->SetMarkerColor(kRed+1);
     hWidthDMCorrvsMP->SetMarkerColor(kRed+1);
     //Printf("X range = %.0f-%.0f, nbins = %d",hMeanDMvsMP->GetBinLowEdge(1), hMeanDMvsMP->GetBinLowEdge(nMbins+1), nMbins);
     //Printf("X range = %.0f-%.0f, nbins = %d",hMeanDMCorrvsMP->GetBinLowEdge(1), hMeanDMCorrvsMP->GetBinLowEdge(nMbins+1), nMbins);
     
     MMin=MMinCp;
     
     for(Int_t iM=0;iM<nMbins;iM++){ //loop on M det bins
     	TString ptstring=Form("%.0f<M_{part}<%.0f GeV/c",MMin, MMin+Minterval);
     	
     	TPaveText *pv=new TPaveText(0.1,0.4,0.45,0.55,"NDC");
     	pv->SetFillStyle(0);
     	pv->SetBorderSize(0);
     	pv->AddText(ptstring);
     	
     	Int_t MPBin[2] = {fhnDeltaMass->GetAxis(axselrange)->FindBin(MMin), fhnDeltaMass->GetAxis(axselrange)->FindBin(MMin+Minterval)-1};
     	Printf("Projections in M bin range %d-%d",MPBin[0],MPBin[1]);
     	
     	TH1D* hDMMP = (TH1D*)DoProjections(fhnDeltaMass, axselrange, MPBin, axproj, Form("hDM%dMP%d-%d",axproj,iM,j));
     	
     	TH1D* hDMCorrMP = DoProjections(fhnDeltaMassCorr, axselrange, MPBin, axproj, Form("hDM%dCorrMP%d-%d",axproj,iM,j));
     	hDMCorrMP->SetLineStyle(2);
     	
     	cDMvsMP->cd(iM+1);
     	//gPad->SetLogy();
     	hDMMP->Draw();
     	hDMCorrMP->Draw("sames");
     	pv->Draw();
     	PrintInfoOnCanvas(hDMMP, 0.15,0.25,0.9,0.4, "(det)", cDMvsMP->cd(iM+1));
     	PrintInfoOnCanvas(hDMCorrMP, 0.15,0.10,0.9,0.35, "(det corr)", cDMvsMP->cd(iM+1));
     	
     	FillTrendings(hMeanDMvsMP, hWidthDMvsMP, iM, hDMMP);
     	FillTrendings(hMeanDMCorrvsMP, hWidthDMCorrvsMP, iM, hDMCorrMP);
     	//     Int_t MPBin[2] = {fhnDeltaMass->GetAxis(2)->FindBin(MMin), fhnDeltaMass->GetAxis(2)->FindBin(MMin+Minterval)-1};
     	
     	MMin+=Minterval;
     }
     SaveCv(cDMvsMP);
     TCanvas *cSummaryDMMP = new TCanvas(Form("cSummaryDM%dMP%d",axproj,j), Form("Summary Delta Mass vs M part, %.0f<p_{T,part}<%.f",pTPartLargeMin, pTPartLargeMin+pTPartLargeinterval), 500,900);
     cSummaryDMMP->Divide(1,2);
     cSummaryDMMP->cd(1);
     hMeanDMvsMP->SetMaximum(0);
     hMeanDMvsMP->Draw("P");
     hMeanDMCorrvsMP->Draw("Psames");
     cSummaryDMMP->cd(2);
     hWidthDMvsMP->Draw("P");
     hWidthDMCorrvsMP->Draw("Psames");
     SaveCv(cSummaryDMMP);
     
     //Printf("Check pt bin = %.f",pTPartLargeMin);
     pTPartLargeMin+=pTPartLargeinterval;
     //Printf("Check pt bin after + %.0f => %.f",pTPartLargeinterval,pTPartLargeMin);
     
  }
  ResetFullRanges(fhnDeltaMass);
  ResetFullRanges(fhnDeltaMassCorr);

    //Double_t massBinW = fhnDeltaMass->GetAxis(1)->GetBinWidth(1);
  
  
} //end main

void FeedOut(Int_t axproj=0,TString strIn="JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC.root", TString strLst = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC"){
   //study the distribution of jets going to the right pT bin after correction
   
   SetStyle();
   
   TList *lst = ReadFile(strIn, strLst);
   if(!lst) return;
   
   THnSparse *fhnDeltaMass = (THnSparse*)lst->FindObject("fhnDeltaMass");
   THnSparse *fhnDeltaMassCorr = (THnSparse*)lst->FindObject("fhnDeltaMassCorr");
   
   Int_t axselrange2=4; //pTpart
   const Int_t npTPbins=6;
   Double_t pTPinterval=20., pTPMin=20.; 
   Double_t pTPMax=pTPMin+npTPbins*pTPinterval; //GeV/c
   
   Int_t axselrange=3; //pTdet
   const Int_t npTbins = 6;
   Double_t pTinterval=20., pTMin=20.,  pTMax=pTMin+npTbins*pTinterval; //GeV/c
   Double_t pTMinCp=pTMin;
   
   TLegend *legFeedOut = new TLegend(0.6,0.15,0.85,0.45);
   legFeedOut->SetFillStyle(0);
   legFeedOut->SetBorderSize(0);
   
   TH1F* hReq = new TH1F("hReq", "Ratio equivalent pT Part and Det;#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hReq->SetMarkerStyle(20);
   legFeedOut->AddEntry(hReq, "p_{T,Det} = p_{T,Part}","P");
   TH1F* hRCorreq = new TH1F("hRCorreq", "Ratio Corrected equivalent pT Part and Det;#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hRCorreq->SetMarkerStyle(20);
   hRCorreq->SetMarkerColor(kRed-3);
   
   
   TH1F* hRlo = new TH1F("hRlo", "Ratio pT Part lower than pT Det;#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hRlo->SetMarkerStyle(21);
   legFeedOut->AddEntry(hRlo, "p_{T,Det} > p_{T,Part}","P");
   TH1F* hRCorrlo = new TH1F("hRCorrlo", "Ratio Corrected pT Part lower than Det;#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hRCorrlo->SetMarkerStyle(21);
   hRCorrlo->SetMarkerColor(kRed-3);

   TH1F* hRhi = new TH1F("hRhi", "Ratio pT Part higher than pT Det;#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hRhi->SetMarkerStyle(33);
   legFeedOut->AddEntry(hRhi, "p_{T,Det} < p_{T,Part}","P");
   TH1F* hRCorrhi = new TH1F("hRCorrhi", "Ratio Corrected pT Part higher than Det'#it{p}_{T,jet} (GeV/#it{c});N_{jet}(p_{T,Det-sel})/N_{jet}", npTbins, pTMin, pTMax);
   hRCorrhi->SetMarkerStyle(33);
   hRCorrhi->SetMarkerColor(kRed-3);

   legFeedOut->AddEntry(hRCorreq, "Corrected","P");
   
   for(Int_t j=0;j<npTPbins;j++){
      //fixed pT Part bin
      Int_t pTPBin[2] = {fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPMin), fhnDeltaMass->GetAxis(axselrange2)->FindBin(pTPMin+pTPinterval)-1};
      
      TString ptpstring=Form("%.0f<p_{T,jet}^{Part}<%.0f GeV/c",pTPMin, pTPMin+pTPinterval);
      Printf("pT PART %s", ptpstring.Data());
      
      TPaveText *pvP=new TPaveText(0.1,0.6,0.45,0.75,"NDC");
      pvP->SetFillStyle(0);
      pvP->SetBorderSize(0);
      pvP->SetTextColor(kGreen+3);
      pvP->AddText(ptpstring);
      
      //Integrals
      Double_t totalhDM = 0, totalhDMCorr = 0,
      totalhDMjeqipT = 0, totalhDMCorrjeqipT = 0,
      totalhDMjloipT = 0, totalhDMCorrjloipT = 0,
      totalhDMjhiipT = 0, totalhDMCorrjhiipT = 0;
 
      TH1D* hDMpTD = (TH1D*)DoProjections(fhnDeltaMass, axselrange2, pTPBin, axproj, Form("hDM%dpTP%d",axproj,j));
      totalhDM = hDMpTD->Integral();
      TH1D* hDMCorrpTD = DoProjections(fhnDeltaMassCorr, axselrange2, pTPBin, axproj, Form("hDM%dCorrpTP%d",axproj,j));
      hDMCorrpTD->SetLineStyle(2);
      totalhDMCorr = hDMCorrpTD->Integral();

      for(Int_t ipT=0;ipT<npTbins;ipT++){
      	 TString ptstring=Form("%.0f<p_{T,jet}<%.0f GeV/c",pTMin, pTMin+pTinterval);
      	 
      	 TPaveText *pv=new TPaveText(0.1,0.4,0.45,0.55,"NDC");
      	 pv->SetFillStyle(0);
      	 pv->SetBorderSize(0);
      	 pv->AddText(ptstring);
      	 
      	 Int_t pTDBin[2] = {fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin), fhnDeltaMass->GetAxis(axselrange)->FindBin(pTMin+pTinterval)-1};
      	 //Printf("Projections in pT bin range %d-%d",pTDBin[0],pTDBin[1]);
      	 if(j == ipT){
      	    //pTD bin = pTP bin
      	    TH1D* hDMjeqpTD = (TH1D*)DoProjections(fhnDeltaMass, axselrange, pTDBin, axproj, Form("hDM%dpTD%d",axproj,ipT));
      	    totalhDMjeqipT = hDMjeqpTD->Integral();
      	    
      	    TH1D* hDMCorrjeqpTD = DoProjections(fhnDeltaMassCorr, axselrange, pTDBin, axproj, Form("hDM%dCorrpTD%d",axproj,ipT));
      	    hDMCorrjeqpTD->SetLineStyle(2);
      	    totalhDMCorrjeqipT = hDMCorrjeqpTD->Integral();
      	    
      	    hReq->SetBinContent(ipT+1,totalhDMjeqipT/totalhDM);
      	    
      	    hRCorreq->SetBinContent(ipT+1,totalhDMCorrjeqipT/totalhDMCorr);
      	 }
      	 
       	 if(j < ipT){
       	    
      	    //pTD bin > pTP bin
      	    TH1D* hDMjlopTD = (TH1D*)DoProjections(fhnDeltaMass, axselrange, pTDBin, axproj, Form("hDM%dpTD%d",axproj,ipT));
      	    totalhDMjloipT = hDMjlopTD->Integral();
      	    
      	    TH1D* hDMCorrjlopTD = DoProjections(fhnDeltaMassCorr, axselrange, pTDBin, axproj, Form("hDM%dCorrpTD%d",axproj,ipT));
      	    hDMCorrjlopTD->SetLineStyle(2);
      	    totalhDMCorrjloipT = hDMCorrjlopTD->Integral();
      	    
      	    hRlo->SetBinContent(ipT+1,totalhDMjloipT/totalhDM);
      	    
      	    hRCorrlo->SetBinContent(ipT+1,totalhDMCorrjloipT/totalhDMCorr);
      	 }
     	 
       	 if(j > ipT){
       	    
      	    //pTD bin < pTP bin
      	    TH1D* hDMjhipTD = (TH1D*)DoProjections(fhnDeltaMass, axselrange, pTDBin, axproj, Form("hDM%dpTD%d",axproj,ipT));
      	    totalhDMjhiipT = hDMjhipTD->Integral();
      	    
      	    TH1D* hDMCorrjhipTD = DoProjections(fhnDeltaMassCorr, axselrange, pTDBin, axproj, Form("hDM%dCorrpTD%d",axproj,ipT));
      	    hDMCorrjhipTD->SetLineStyle(2);
      	    totalhDMCorrjhiipT = hDMCorrjhipTD->Integral();
      	    
      	    hRhi->SetBinContent(ipT+1,totalhDMjhiipT/totalhDM);
      	    
      	    hRCorrhi->SetBinContent(ipT+1,totalhDMCorrjhiipT/totalhDMCorr);
      	 }

      	 pTMin+=pTinterval;
      }
      
      //reset needed
      pTMin=pTMinCp;
      pTPMin+=pTPinterval;
      ResetFullRanges(fhnDeltaMass);
      ResetFullRanges(fhnDeltaMassCorr);
   }
   
   TCanvas *cFeedOut = new TCanvas("cFeedOut","Feed Out", 800,800);
   cFeedOut->cd();
   gPad->SetLogy();
   hReq->SetMinimum(5e-4);
   hReq->Draw("P");
   hRCorreq->Draw("Psames");
   hRlo->Draw("Psames");
   hRCorrlo->Draw("Psames");
   hRhi->Draw("Phistsames");
   hRCorrhi->Draw("Phistsames");
   legFeedOut->Draw();
   SaveCv(cFeedOut);
}

//functions
TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Int_t* pTPBin, Int_t axisP, char* hprjname){
     fhnDeltaMass->    GetAxis(axisR)->SetRange(pTPBin[0],pTPBin[1]); //pT det
     TH1D* hDMpTP=(TH1D*)fhnDeltaMass->Projection(axisP);
     hDMpTP->SetName(hprjname);
     hDMpTP->SetLineWidth(2);
     return hDMpTP;
}

TPaveText* PrintInfoOnCanvas(TH1D* h, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString text, TVirtualPad* cvdest){
     TPaveText* pvinfo = new TPaveText(x1, y1, x2, y2, "NDC");
     pvinfo->SetFillStyle(0);
     pvinfo->SetBorderSize(0);
     Double_t meanDM = h->GetMean(), meanDMerr = h->GetMeanError(), widthDM = h->GetRMS();
     //print info on canvas
     pvinfo->AddText(Form("%s Mean = %.2f +- %.2f, Width %.2f ", text.Data(), meanDM, meanDMerr, widthDM));
     cvdest->cd();
     pvinfo->Draw();
     
     return pvinfo;
}

void DefineTrendingHistos(TString nameX, TString nameY, TString titleX, TString titleY, Int_t nbins, Double_t min, Double_t max, TH1D*& hMean, TH1D*& hWidth ){
  hMean = new TH1D(Form("hMean%svs%s", nameY.Data(), nameX.Data()), Form("Mean ; %s ; Mean %s", titleX.Data(), titleY.Data()), nbins, min, max);
  hWidth = new TH1D(Form("hSigma%svs%s", nameY.Data(), nameX.Data()), Form("Sigma; %s ; Width %s", titleX.Data(), titleY.Data()), nbins, min, max);
  hMean->SetMarkerStyle(20);
  hWidth->SetMarkerStyle(20);
}

void FillTrendings(TH1D *hMeanDMvspT, TH1D* hWidthDMvspT, Int_t ipT,TH1D* hDM){
   
   Double_t meanDM = hDM->GetMean(), meanDMerr = hDM->GetMeanError(), widthDM = hDM->GetRMS();
   //fill summay histograms
   //if(meanDM<1e-7) meanDM=1e-4;if(widthDM<1e-7) widthDM=1e-4;
   Printf("%d) Mean %f, Width = %f", ipT,meanDM, widthDM);
   hMeanDMvspT->SetBinContent(ipT+1, meanDM);
   //hMeanDMvspT->SetBinError(ipT+1, meanDMerr);
   hWidthDMvspT->SetBinContent(ipT+1, widthDM);
}

void ResetFullRanges(THnSparse *hns){
   for (Int_t i=0;i<hns->GetNdimensions();i++){
      hns->GetAxis(i)->SetRange(0,-1);
   }
}


void SaveCv(TCanvas* c,TString suffix, Int_t format){
   if(format > 0) c->SaveAs(Form("%s%s.png",c->GetName(),suffix.Data()));
   if(format > 1) c->SaveAs(Form("%s%s.pdf",c->GetName(),suffix.Data()));
   if(format > 2) c->SaveAs(Form("%s%s.eps",c->GetName(),suffix.Data()));
   if(format > 3 || format==0) c->SaveAs(Form("%s%s.root",c->GetName(),suffix.Data()));
}

void CanvasDim(Int_t nptbins, Int_t& px, Int_t& py, Int_t& dimx, Int_t& dimy){
   py=1;
   px=nptbins;
   dimy=300;

   if(nptbins>4){
      px=2;
      py=nptbins/px;
      if((nptbins % 2) !=0) px+=1; 
      
      dimy*=py;
   }
   dimx=400;
   dimx*=px;
   
   if(nptbins>10){
      dimy=200;
      px=3;
      py=nptbins/px;
      if((nptbins % px) !=0) px+=1; 
      
      dimy*=py;
   }
   
   dimx=200;
   dimx*=px;  
   return;
}

void SetStyle(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleXSize(0.04);

}

TList* ReadFile(TString strIn, TString strLst){
   
   TFile *f = new TFile(strIn.Data());
   if(!f->IsOpen()){
      Printf("File %s not found", strIn.Data());
      return 0x0;
   }
   TList *lst = static_cast<TList*>(f->Get(strLst.Data()));
   if(!lst){
      Printf("Error, list not found");
      return 0x0;
   }
   return lst;
}
