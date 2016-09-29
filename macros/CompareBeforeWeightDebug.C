#include <TCanvas.h>
#include <TString.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TLegend.h>
#include "/data/Work/MyCodeJetMass/classes/WeightClass.h"

void UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type);

void CompareBeforeWeightDebug(){
   
   const Int_t nSeries = 2;
   TString paths[nSeries] = {"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1096/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //paths[1] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/";
   TString listnameInp[nSeries] = {"PrepareInputForEmbedding", "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme"};
   //listnameInp[0] = "JetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme";
   //listnameInp[1] = "PrepareInputForEmbedding";
   
   TString respname[nSeries] = {"hResponse", "fhnMassResponse"};
   //respname[0] = "fhnMassResponse";
   //respname[1] = "hResponse";
   
   TString legF[nSeries] = {"Train1087", "Train576"};
   //legF[0] = "Train1096";
   //legF[1] = "Train1087";
   TString filename = "AnalysisResults.root";
   
   Int_t firstBin = 0;
   Int_t npthb =10;//2 ;//
   //
   THnSparseF **fhResponse = new THnSparseF*[npthb];
   TH1D* hRatio[4][npthb][nSeries];
   TH1D* hRespPj[4][npthb][nSeries];
   //canvas definition
   //TCanvas *cProjResp = new TCanvas(Form("cProjResp%s",legF[0].Data()), "TT");
   
   TCanvas *cProjResp = new TCanvas(Form("cProjResp%s%s", legF[0].Data(),  legF[1].Data()), Form("Response before weighting %s vs %s", legF[0].Data(),  legF[1].Data()), 1500, 1500);
   cProjResp->Divide(2,2);
   TCanvas *cProjRespRatio = new TCanvas(Form("cProjRespRatio%s%s", legF[0].Data(),  legF[1].Data()), Form("Response Ratio  %s vs %s", legF[0].Data(),  legF[1].Data()), 1500, 1500);
   cProjRespRatio->Divide(2,2);

   TLegend *legpTh = new TLegend(0.6, 0.6, 0.9, 0.9);
   legpTh->SetBorderSize(0);
   legpTh->SetFillStyle(0);
   
   TLegend *legFiles = new TLegend(0.3, 0.8, 0.6, 0.9);
   legFiles->SetBorderSize(0);
   legFiles->SetFillStyle(0);
   
   for(Int_t ifile = 0; ifile < nSeries; ifile++){
   //for(Int_t ifile = 0; ifile < 1; ifile++){
   
      WeightClass *doweight = new WeightClass(Form("weight%d", ifile), Form("Weight %d", ifile), npthb);
      doweight->Initialise(paths[ifile], listnameInp[ifile], filename, listnameInp[ifile]);
      
      //loop on pT hard bins
      for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
     	 //for(Int_t ipthb = firstBin; ipthb<1; ipthb++){
      	 
      	 Bool_t set = doweight->SetResponse(ipthb, respname[ifile]);
      	 if(!set) continue;
      	 Printf("Response name %s", doweight->GetResponsePtHbin(ipthb)->GetName());
      	 TH1D **hRespPjTmp = doweight->All1DResponseProj(ipthb, ipthb, Form("hPrjTmpFile%dPtHb%d", ifile, ipthb));
      	 
      	 
      	 for(Int_t ipj = 0; ipj<4; ipj++){
      	    
      	    hRespPj[ipj][ipthb][ifile] = (TH1D*)hRespPjTmp[ipj]->Clone(Form("hPrjFile%dPtHb%d_%d", ifile, ipthb, ipj));
      	    //hRespPj[ipj][ipthb][ifile]->SetName(Form("hPrjFile%dPtHb%d_%d", ifile, ipthb, ipj));
      	    
      	    Printf("%s bin %d -> number of bins %d", legF[ifile].Data(), ipthb, hRespPj[ipj][ipthb][ifile]->GetNbinsX());
      	    
      	    hRespPj[ipj][ipthb][ifile]->SetMarkerStyle(20+ifile);
      	    hRatio[ipj][ipthb][ifile] = (TH1D*)hRespPj[ipj][ipthb][ifile]->Clone(Form("hRatio%d%d_pTHb%d", ifile, ipj, ipthb));
      	      
       	       
      	    cProjResp->cd(ipj+1);
      	    gPad->SetLogy();
      	    
      	    if((ifile == 0) && (ipthb == firstBin)) hRespPj[ipj][ipthb][ifile]->Draw();
      	    else hRespPj[ipj][ipthb][ifile]->Draw("sames");
      	    Printf("Names: %s, %s, %s", hRespPj[ipj][ipthb][ifile]->GetName(), hRatio[ipj][ipthb][ifile]->GetName(), hRespPjTmp[ipj]->GetName());
      	   
      	 }
      	 
      	 if(ipthb == firstBin) legFiles->AddEntry(hRespPj[0][ipthb][ifile], legF[ifile], "P");
      	 if((ifile == 0)) legpTh->AddEntry(hRespPj[0][ipthb][ifile], Form("pThard %d", ipthb+1), "L");
      	 
      } //loop on pThard
      doweight->DrawCrossSec();
      doweight->DrawNTrials();
   }//loop on files
   
   cProjResp->cd(1);
   legpTh->Draw();
   cProjResp->cd(3);
   legFiles->Draw();
   //SaveCv(cProjResp);
   cProjResp->SaveAs("savethisplease.pdf");
 
   //TH1D*** hRatiopjRebin = new TH1D**[npthb];
    for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
      //for(Int_t ipthb = firstBin; ipthb<1; ipthb++){
      	 //hRatiopjRebin[ipthb] = new TH1D*[4];
      	 //loop on variables
       for(Int_t ipj = 0; ipj<4; ipj++){
       	  TH1* hRespPjDiv1 = 0x0;
       	  TH1* hRespPjDiv2 = 0x0;
       	  Printf("Before %d %d", hRatio[ipj][ipthb][0]->GetNbinsX(), hRatio[ipj][ipthb][1]->GetNbinsX());
       	  UniformTH1FForDivide(hRatio[ipj][ipthb][0], hRatio[ipj][ipthb][1], hRespPjDiv1, hRespPjDiv2, "TH1D");
       	  if(!hRespPjDiv1 || !hRespPjDiv2) {
       	     Printf("Error! %p %p", hRespPjDiv1, hRespPjDiv2); return;
       	  }
       	  Printf("%s %s", hRespPjDiv1->GetName(), hRespPjDiv2->GetName());
       	  Printf("%d %d", hRespPjDiv1->GetNbinsX(), hRespPjDiv2->GetNbinsX());
       	  //hRatiopjRebin[ipj][npthb] = (TH1D*)hRespPjDiv1->Clone(Form("hRatioRebin_pthb%d_pj%d", ipthb, ipj));
       	  //Printf("Now it is %s / %s", hRatiopjRebin[ipj][npthb]->GetName(), hRespPjDiv2->GetName());
       	  hRespPjDiv1->Divide(hRespPjDiv2);
      	  //hRatiopjRebin[ipj][ipthb]->Divide(hRespPjDiv2);
      	  Printf("Division successful");
      	  cProjRespRatio->cd(ipj+1);
      	  //gPad->SetLogy();
       	     
       	  //hRatiopjRebin[ipj][ipthb]->SetLineStyle(2);
       	  if(ipthb == firstBin) hRespPjDiv1->Draw(); //hRatiopjRebin[ipj][ipthb]->Draw();
       	  else hRespPjDiv1->Draw("sames"); //hRatiopjRebin[ipj][ipthb]->Draw("sames");
       	  
       } //end loop on variables
    } //end loop on pT hard bins
    
}

void UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type){
   if(!h1 || !h2){
      Printf("Need valid input!");
      return;
   }
   
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   if((nBins1 == nBins2) && ((binW1-binW2) < 1e-6) &&  ((low1-low2) < 1e-6) && ((up1-up2) < 1e-6)) {
      Printf("Histograms can be divided");
      h1bis = (TH1*)h1->Clone(Form("%sb", h1->GetName()));
      h2bis = (TH1*)h2->Clone(Form("%sb", h2->GetName()));
      return;
   }
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2);
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Common) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2);
   
   Int_t rebin = (Int_t)(commonW/minWidth);
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin);
   else h2->Rebin(rebin);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   if(type == "TH1F"){
      h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
   } else {
      if(type == "TH1D"){
      	 h1bis = new TH1D(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      	 h2bis = new TH1D(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      } else {
      	 Printf("%s Not defined!", type.Data());
      	 return;
      }
   }
   h1bis->Sumw2();
   h1bis->SetLineColor(h1->GetLineColor());
   h1bis->SetLineWidth(h1->GetLineWidth());
   h1bis->SetMarkerColor(h1->GetMarkerColor());
   h1bis->SetMarkerStyle(h1->GetMarkerStyle());
   
   h2bis->Sumw2();
   h2bis->SetLineColor(h2->GetLineColor());
   h2bis->SetLineWidth(h2->GetLineWidth());
   h2bis->SetMarkerColor(h2->GetMarkerColor());
   h2bis->SetMarkerStyle(h2->GetMarkerStyle());
   
   //refill the original histograms using the same binning
   for(Int_t i = 0; i < commonNBins; i++){
      
      //First histogram
      Float_t bincentre = h1bis ->GetBinCenter(i+1);
      Int_t fillingbin = h1->FindBin(bincentre);
      Double_t fillWith = 0, error = 0;
      if(fillingbin>0 && fillingbin <= nBins1) {
      	 Printf"Here, debug %d", nBins1);
      	 fillWith = (Double_t)h1->GetBinContent(fillingbin);
      	 error = (Double_t)h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h1->GetName(), bincentre, i+1, fillingbin, fillWith);
      
      //second histogram
      bincentre = h2bis ->GetBinCenter(i+1);
      fillingbin = h2->FindBin(bincentre);
      fillWith = 0; error = 0;
      if(fillingbin>0 && fillingbin <= nBins2) {
      	 fillWith = h2->GetBinContent(fillingbin);
      	 error = h2->GetBinError(fillingbin);
      }
      h2bis ->SetBinContent(i+1, fillWith);
      h2bis ->SetBinError(i+1, error);
      Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h2->GetName(), bincentre, i+1, fillingbin, fillWith);
   }
   
   
}
