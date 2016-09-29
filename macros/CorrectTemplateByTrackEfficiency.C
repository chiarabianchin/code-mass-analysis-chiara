#include<THnSparse.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TH3F.h>
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

void SaveCv(TCanvas* c,TString suffix="", Int_t format=2);

void SaveTemplate(TString input = "/data/Work/jets/JetMass/DetectorCorrections/FastSimulation/UseTrResolution/AnalysisResultsCorrectTemplate.root", TString listName = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0000_E_scheme_TC", TString output = "/data/Work/jets/JetMass/DetectorCorrections/FastSimulation/UseTrResolution/TemplateTrReso.root"){
   //open
   TFile *fin = new TFile(input);
   if(!fin->IsOpen()){
      Printf("Input file %s not found", input.Data());
      return;
   }
   
   TString hTempName = "fh3JetPtDRTrackPt_0";
   TList *l = (TList*) fin->Get(listName);
   if(!l){
      Printf("List %s not found", listName.Data());
      fin->ls();
      return;
   }
   TH3F *h3Template = (TH3F*)l->FindObject(hTempName);
   if(!h3Template){
      Printf("File %s not found", hTempName.Data());
      l->ls();
      return;
   } else if(h3Template->GetEntries()==0){
      Printf("Histogram is empty");
      return;
   }
   TCanvas *c = new TCanvas("c", "Show template projections", 800,800);
   Printf("X axis = %s, Y axis = %s, Z axis = %s", h3Template->GetXaxis()->GetTitle(), h3Template->GetYaxis()->GetTitle(), h3Template->GetZaxis()->GetTitle());
   h3Template->GetXaxis()->SetRangeUser(60,80);
   TH2D* htpTR = (TH2D*)h3Template->Project3D("yz"); //r vs pttrack
   c->cd();
   htpTR->Draw();

   TFile *fout = new TFile(output, "recreate");
   fout->cd();
   h3Template->Write();
   fout->Close();
   
   Printf("SAVED");
}
void CorrectTemplateByTrackEfficiency(){
   
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   //gStyle->SetLineWidth(2);
   TString pathEff = "/data/Work/jets/JetMass/DetectorCorrections/LHC12a15e/Train399/EfficiencyCent10_LHC12a15e.root";
   TString hEffName = "hEffPosSum";
   
   TString pathTemplate = "/data/Work/jets/JetMass/DetectorCorrections/LHC12a15e/Train399/AnalysisResultsWeighted.root";
   TString listName = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC";
   TString hTempName = "fh3JetPtDRTrackPt_0";
   
   //open
   TFile *finE = new TFile(pathEff);
   if(!finE->IsOpen()){
      Printf("Efficiency file %s not found", pathEff.Data());
      return;
   }
   
   TH1D* hEff = (TH1D*)finE->Get(hEffName);
   if(!hEff){
      Printf("Efficiency histogram not found");
      finE->ls();
      return;
   }
   
   TCanvas *cEff = new TCanvas ("cEff","Efficiency positive tracks",800,800);
   cEff->cd();
   hEff->GetXaxis()->SetTitle("#it{p}_{T}");
   hEff->Draw();
   
   Int_t nBinsEf = hEff->GetXaxis()->GetNbins();
   Printf("Efficiency:\nX axis %s Nbins = %d", hEff->GetXaxis()->GetTitle(), nBinsEf);
   Double_t smallBinW = 0.9;
   Double_t pTFitRange[2]={6,100}; //reset in the loop, silent warning
   for(Int_t ibeff=0;ibeff<nBinsEf;ibeff++){
      Double_t bw = hEff->GetBinWidth(ibeff+1);
      if(bw>smallBinW){
      	 pTFitRange[0] = hEff->GetBinLowEdge(ibeff);
      	 
      	 break;
      
      }
   
   }
   pTFitRange[1] = hEff->GetBinLowEdge(nBinsEf+1);
   Printf("Fit efficiency from %.1f to %.1f", pTFitRange[0], pTFitRange[1]);
   TF1* fitfuncEff = new TF1("fitfuncEff", "[0]+x*[1]", pTFitRange[0], pTFitRange[1]);
   hEff->Fit("fitfuncEff", "R");
   
   TFile *finT = new TFile(pathTemplate);
   if(!finT->IsOpen()){
      Printf("Efficiency file %s not found", pathTemplate.Data());
      return;
   }
   
   TList *list = (TList*) finT->Get(listName);
   if(!list){
      Printf("List not found");
      finT->ls();
      return;
  
   }
   
   TH3F* hTemplate = (TH3F*) list->FindObject(hTempName);
   if(!hTemplate){
      Printf("Efficiency histogram %s not found", hTempName.Data());
      list->ls();
      return;
  
   }
   
   Int_t nBinsX = hTemplate->GetXaxis()->GetNbins(), nBinsY = hTemplate->GetYaxis()->GetNbins(),  nBinsZ = hTemplate->GetZaxis()->GetNbins();
   
   Printf("Template: \nX axis %s Nbins = %d,\nY axis %s Nbins = %d,\nZ axis %s Nbins = %d", hTemplate->GetXaxis()->GetTitle(), nBinsX, hTemplate->GetYaxis()->GetTitle(), nBinsY, hTemplate->GetZaxis()->GetTitle(), nBinsZ);
   //Double_t rangeXT[2] = {hTemplate->GetXaxis()->GetBinLowEdge(1), hTemplate->GetXaxis()->GetBinLowEdge(nBinsX+1) };
   //Double_t rangeYT[2] = {hTemplate->GetYaxis()->GetBinLowEdge(1), hTemplate->GetYaxis()->GetBinLowEdge(nBinsY+1) };
   //Double_t rangeZT[2] = {hTemplate->GetZaxis()->GetBinLowEdge(1), hTemplate->GetZaxis()->GetBinLowEdge(nBinsZ+1) };
   //
   TCanvas *cpTtrack = new TCanvas("cpTtrack", "p_{T,track} distribution Template",800,800);
   TH1D* hpTtrackT = (TH1D*)hTemplate->ProjectionZ("hpTtrackT");
   hpTtrackT->Sumw2();
   hpTtrackT->SetLineWidth(2);
   cpTtrack->cd();
   gPad->SetLogy();
   hpTtrackT->Draw("hist");
   
   //TH1D* hpTtrackTCorr = (TH1D*)hpTtrackT->Clone("hpTtrackTCorr");   
   //hpTtrackTCorr->Divide(hEff);
   
   //TCanvas *cpTtrackCorr = new TCanvas("cpTtrackCorr", "Corrected p_{T,track} distribution Template",800,800);
   
   //cpTtrackCorr->cd();
   
  
   TH1D* hpTtrackTCorr = (TH1D*)hpTtrackT->Clone("hpTtrackTCorr"); 
   hpTtrackTCorr->SetTitle("p_{T,track} Eff corrected;p_{T,track} (Gev/#it{c}");
   hpTtrackTCorr->Sumw2();
   hpTtrackTCorr->SetLineColor(kRed);
   hpTtrackTCorr->SetLineWidth(2);
   
   TH1D* hEffForCorr = (TH1D*)hpTtrackT->Clone("hEffForCorr");
   hEffForCorr->SetTitle("Efficiency with more bins ; p_{T,track} (Gev/#it{c}");
   hEffForCorr->Sumw2();
   hEffForCorr->SetMarkerStyle(22);
   hEffForCorr->SetMarkerColor(kCyan);
   
   TH1D* hEffInterpolated = (TH1D*)hpTtrackT->Clone("hEffInterpolated");
   hEffInterpolated->SetMarkerStyle(33);
   hEffInterpolated->SetMarkerColor(kRed-2);
   
   TH3F* hTemplateCorr = (TH3F*) hTemplate->Clone("hTemplateCorr"); 
   hTemplateCorr->Sumw2();
   Printf("Loop 1D");
   for(Int_t i=0;i<nBinsZ;i++){
      Double_t binContpTTrack = hpTtrackT->GetBinContent(i+1);
      Double_t binCentreT=hpTtrackT->GetBinCenter(i+1);
      Int_t bin = hEff->FindBin(binCentreT);
      Double_t binEff = hEff->GetBinContent(bin);
      hEffForCorr->SetBinContent(i+1,binEff);
      hEffForCorr->SetBinError(i+1,hEff->GetBinError(i+1));
      Double_t binContpTTrackCorr = binContpTTrack/binEff;
      //Printf("Bin %d, centre %f (binEff orig = %f at bin %d) content %f -> %f", i+1, binCentreT,binEff, bin, binContpTTrack, binContpTTrackCorr);
      hpTtrackTCorr->SetBinContent(i+1, binContpTTrackCorr);
      
      //by interpolation
      if(hEff->GetBinWidth(bin) > smallBinW){
     	 //Printf("Eval from Fit %d, %.2f", i+1, binCentreT);
     	 
      hEffInterpolated->SetBinContent(i+1,fitfuncEff->Eval(binCentreT));      } else {
       	 Double_t effInterp = hEff->Interpolate(binCentreT);
      	 hEffInterpolated ->SetBinContent(i+1,effInterp);
      	 hEffInterpolated ->SetBinError(i+1,0.);
      	 
      }
   }
   Printf("Loop ended, efficiency ready");
   
   cEff->cd();
   //hEffForCorr->Draw("samesP");
   hEffInterpolated->Draw("samesP");
   SaveCv(cEff);
   
   cpTtrack->cd();
   hpTtrackTCorr->Draw("sames");

   TH1D* hpTtrackTCorrOverpTtrackRatio = (TH1D*)hpTtrackTCorr->Clone("hpTtrackTRatio");
   hpTtrackTCorrOverpTtrackRatio->Divide(hpTtrackT);
   hpTtrackTCorrOverpTtrackRatio->SetMarkerStyle(20);
   hpTtrackTCorrOverpTtrackRatio->SetMarkerColor(kGreen+4);
  
   
   cpTtrack->cd();
   hpTtrackTCorrOverpTtrackRatio->Draw("samesP");

   TH1D* hpTtrackTDivEff = (TH1D*)hpTtrackT->Clone("hpTtrackTDivEff");
   Printf("Operation pTTrack/Eff right binning %d = %d",hpTtrackTDivEff->GetNbinsX(),hEffForCorr->GetNbinsX());
   hpTtrackTDivEff->Divide(hEffForCorr);
   Printf("Division done");
   hpTtrackTDivEff->SetLineStyle(2);
   hpTtrackTDivEff->SetLineColor(kCyan);
   cpTtrack->cd();
   hpTtrackTDivEff->Draw("sames");
   
   TH1D* hpTtrackCorrByHandHistoRatio = (TH1D*)hpTtrackTDivEff->Clone("hpTtrackCorrByHandHistoRatio");
   hpTtrackCorrByHandHistoRatio->SetYTitle("Ratio");
   hpTtrackCorrByHandHistoRatio->Divide(hpTtrackTCorr);
   hpTtrackCorrByHandHistoRatio->SetLineWidth(2);
   
   TH1D* hpTtrackTDivEffInterp = (TH1D*)hpTtrackT->Clone("hpTtrackTDivEffInterp");
   hpTtrackTDivEffInterp->Divide(hEffInterpolated);
   hpTtrackTDivEffInterp->SetLineStyle(2);
   hpTtrackTDivEffInterp->SetMarkerStyle(33);
   hpTtrackTDivEffInterp->SetLineColor(kRed-2);
   cpTtrack->cd();
   hpTtrackTDivEffInterp->Draw("sames");
   
   //draw the other projections to see the difference before and after correction
   TH1D* hr = (TH1D*)hTemplate->ProjectionY("hr");
   hr->SetLineWidth(2);
   hr->SetLineColor(kGreen+4);
   TCanvas *cr = new TCanvas("cr","Distance to the jet axis",800,800);
   hr->Draw("hist");
   
   TH1D* hpTjet = (TH1D*)hTemplate->ProjectionX("hpTjet");
   hpTjet->SetLineWidth(2);
   hpTjet->SetLineColor(kGreen+4);
   TCanvas *cpTjet = new TCanvas("cpTjet","Distance to the jet axis",800,800);
   hpTjet->Draw("hist");
 
   Printf("Start loop 3D");
   for(Int_t i=0;i<nBinsX;i++){ //
      //Printf("i = %d",i);
      for(Int_t j=0;j<nBinsY;j++){//
      	 for(Int_t k=0;k<nBinsZ;k++){
      	    
      	    Double_t bincontTemplate = hTemplate->GetBinContent(i+1, j+1, k+1);
      	    Double_t bincontTemplateErr = hTemplate->GetBinError(i+1, j+1, k+1);
      	    Double_t eff = hEffInterpolated->GetBinContent(k+1);
      	    //Double_t effErr = hEffInterpolated->GetBinError(k+1);
      	    
      	    Double_t binContTemplCorr = eff ? bincontTemplate/eff : 0;
      	    Double_t binContTemplCorrErr = bincontTemplateErr; //???????????? propagation, same error as before, what??
      	    
      	    hTemplateCorr->SetBinContent(i+1, j+1, k+1, binContTemplCorr);
      	    hTemplateCorr->SetBinError(i+1, j+1, k+1, binContTemplCorrErr);
      	 }
      }
   }
   
   TH1D* hpTtrackTCorproj = (TH1D*)hTemplateCorr->ProjectionZ("hpTtrackTCorproj");
   hpTtrackTCorproj->SetLineStyle(2);
   hpTtrackTCorproj->SetLineWidth(2);
   hpTtrackTCorproj->SetLineColor(kMagenta-4);
   cpTtrack->cd();
   hpTtrackTCorproj->Draw("samesP");
   
   TH1D* hpTCorXCheckRatio = (TH1D*)hpTtrackTCorproj->Clone("hpTCorXCheckRatio");
   hpTCorXCheckRatio->Divide(hpTtrackTCorr);
   hpTCorXCheckRatio->SetYTitle("Ratio");
   hpTCorXCheckRatio->SetLineColor(kRed+2);
   hpTCorXCheckRatio->SetLineWidth(2);
   
   TCanvas *ccheck= new TCanvas ("ccheck", "Checks",800,800);
   ccheck->cd();
   hpTCorXCheckRatio->Draw("hist");
   hpTtrackCorrByHandHistoRatio->Draw("histsames");
   
   TH1D* hrCorr = (TH1D*)hTemplateCorr->ProjectionY("hrCorr");
   hrCorr->SetLineWidth(2);
   hrCorr->SetLineStyle(2);
   hrCorr->SetLineColor(kGreen-4);
   cr->cd();
   hrCorr->Draw("histsames");
   
   
   TH1D* hpTjetCorr = (TH1D*)hTemplateCorr->ProjectionY("hpTjetCorr");
   hpTjetCorr->SetLineWidth(2);
   hpTjetCorr->SetLineStyle(2);
   hpTjetCorr->SetLineColor(kGreen-4);
   cpTjet->cd();
   gPad->SetLogy();
   hpTjetCorr->Draw("histsames");
   
   TCanvas *cRatios = new TCanvas("cRatios", "Ratios Corrected/Uncorrected",800,800);
   TH1D* hrR = (TH1D*)hrCorr->Clone("hrR");
   hrR->SetMarkerStyle(33);
   hrR->Divide(hr);
   TH1D* hpTtrackR = (TH1D*)hpTtrackTCorproj->Clone("hpTtrackR");
   hpTtrackR->Divide(hpTtrackT);
   hpTtrackR->SetMarkerStyle(33);
   TH1D* hpTjetR = (TH1D*)hpTjetCorr->Clone("hpTjetR");
   Printf("N bin ptjet corr %d and uncorr %d", hpTjetR->GetNbinsX(),hpTjet->GetNbinsX());
   hpTjetR->Divide(hpTjet);
   hpTjetR->SetLineStyle(3);
   
   cRatios->cd();
   hrR->Draw();
   hpTjetR->Draw("sames");
   hpTtrackR->Draw("sames");
   
   TCanvas *csummary = new TCanvas("csummary", "Correction steps", 2000,600);
   csummary->Divide(3,1);
   csummary->cd(1); //efficiency
   hEff->Draw();
   hEffInterpolated->Draw("sames");
   //hpTtrackR->Draw("Psames");
   
   csummary->cd(2); //pT track distributions and fraction
   gPad->SetLogy();
   hpTtrackT->Draw("hist");
   hpTtrackTCorproj->Draw("histsames");
   hpTtrackR->Draw("Psames");
   
   csummary->cd(3); //r distribution
   hr->Draw("hist");
   hrCorr->Draw("histsames");
   hrR->Draw("Psames");
   SaveCv(csummary);
}

void SaveCv(TCanvas* c,TString suffix, Int_t format){
   if(format > 0) c->SaveAs(Form("%s%s.png",c->GetName(),suffix.Data()));
   if(format > 1) c->SaveAs(Form("%s%s.pdf",c->GetName(),suffix.Data()));
   if(format > 2) c->SaveAs(Form("%s%s.eps",c->GetName(),suffix.Data()));
   if(format > 3 || format==0) c->SaveAs(Form("%s%s.root",c->GetName(),suffix.Data()));
}


