#include <TProfile.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TF1.h>

#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"
#include "/data/Work/MyCodeJetMass/classes/WeightClass.h"

//global variables
const Int_t nptbins = 4;
Double_t ptlims[nptbins+1] = {40., 60., 80., 100., 120};

//________________________________________________________________________________________________

//methods
TH1D** All1DProjections(THnSparseF *hnsp, Int_t color = 0, TString nameh = "");
TTree* ReadTree (TString strIn, TString nameTree);
void ComparisonFinalResponse(const Int_t ninputs, TString files[], TString legs[], Int_t axSel, Double_t range[]);
void StudyDistributionShapeExcludingpTHbins(const Int_t n, Int_t iPtHb[], Double_t pTcut, TString inputFile, TString inputList, TString hname);
TH1D* GetProjection(THnSparseF *hresp, Int_t iax, TString basenamepj, Int_t colid, Int_t marker);
void CombineRhoBins(Int_t nrhobins, Double_t *rhobinlims, TString *inputFiles, Int_t *composition, Int_t *rhoindexinfile, TString description = "", TString filenameDataRho = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/RhoFromDataRatiosNPtDetbins2Deriv.root");
void ComparisonResponses(Int_t nresp, TString *files, TString *hspname, TString *legs, Int_t baseRatios);
//________________________________________________________________________________________________

//implementation
TH1D** All1DProjections(THnSparseF *hnsp, Int_t color, TString nameh){
   
   const Int_t ndim = hnsp->GetNdimensions();
   TH1D **hproj = new TH1D*[ndim];
   for(Int_t iax = 0; iax < ndim; iax++){
      hproj[iax] = hnsp->Projection(iax);
      if(!nameh.IsNull()) hproj[iax]->SetName(Form("%s_%d", nameh.Data(), iax));
      hproj[iax]->SetTitle( nameh.Data());
      hproj[iax]->SetLineColor(colors[color]);
      hproj[iax]->SetLineWidth(2);
      hproj[iax]->SetMarkerColor(colors[color]);
   }

   return hproj;
}
//________________________________________________________________________________________________

TTree* ReadTree (TString strIn, TString nameTree){
   TFile *fin = new TFile (strIn);
   if(!fin->IsOpen()){
      Printf("%s not found", strIn.Data());
      return 0x0;
   }
   TTree *treeArea = (TTree*) fin->Get(nameTree);
   if(!treeArea) {
      Printf("Tree %s not found %p", nameTree.Data(), treeArea);
      fin->ls();
      return 0x0;
   }
   
   return treeArea;
}

//________________________________________________________________________________________________

Double_t GetWeight(TString pathPythiaFile, Int_t ipthb, TString listPythiaFile = "PrepareInputForEmbedding"){
   
 
}
//________________________________________________________________________________________________
void ComparisonFinalResponse(const Int_t ninputs, TString files[], TString legs[], Int_t axSel, Double_t range[]){
   
   TString histo = "fhResponseFinal";
   Int_t npthb = 10;
   THnSparseF *hsp[ninputs];
   
   WeightClass *onlyproj = new WeightClass("onlyproj", "Weight class to project", npthb);
   onlyproj->SetOrderAxesRespose(0, 1, 2, 3); // leave the order as it is, it's just plotting
   
   TPaveText *pvt = new TPaveText(0.3, 0.7, 0.7, 0.8, "NDC");
   pvt->SetFillStyle(0);
   pvt->SetBorderSize(0);
   
   TCanvas *cres = 0x0;
   
   Int_t countNfound = 0;
   for(Int_t ifile = 0; ifile<ninputs; ifile++){
      
      TFile *fin = new TFile(files[ifile]);
      if(!fin->IsOpen()) {
      	 Printf("%s not found, skip", files[ifile].Data());
      	 continue;
      }
      Printf("Reading file %s", files[ifile].Data());
      hsp[countNfound] = (THnSparseF*)fin->Get(histo);
      if(!hsp[countNfound]) continue;
      
      hsp[countNfound]->SetName(Form("%s_%d", histo.Data(), ifile));
      
      Printf("** Performing selection %.0f - %.0f on the axis %d, %s", range[0], range[1], axSel, hsp[countNfound]->GetAxis(axSel)->GetTitle());
      Int_t rangeB[2] = {hsp[countNfound]->GetAxis(axSel)->FindBin(range[0]), hsp[countNfound]->GetAxis(axSel)->FindBin(range[1]) };
      hsp[countNfound]->GetAxis(axSel)->SetRange(rangeB[0], rangeB[1]);
      if(countNfound > 0){ //there are two, compare
      	 Printf("Compare %s and %s", hsp[countNfound-1]->GetName(), hsp[countNfound]->GetName());
      	 cres = onlyproj->CompareResponses(hsp[countNfound-1], hsp[countNfound], countNfound-1, countNfound, legs[countNfound-1], legs[countNfound]);
      	 
      } else {
      	 //use this if to fill the pavetext
      	 pvt->AddText(Form("%.0f < %s < %.0f", range[0], hsp[countNfound]->GetAxis(axSel)->GetTitle(), range[1]));
      }
      if(cres) {
      	 cres->cd(3);
      	 pvt->Draw();
      }
      //last thing: increase counter
      countNfound++;
   }
}

//_____________________________________________________________________________

void StudyDistributionShapeExcludingpTHbins(const Int_t n, Int_t iPtHb[], Double_t pTcut, TString inputFile, TString inputList, TString hname){
   
   /// Compare the pT (M) shape with a given pT cut using the distribution coming from iPtHb up to 10 wrt to the full one
   /// n is the number of trial of excluded PtH bins
   /// pTcut is the chosen pT cut (to be used in embedding)
   
   Int_t npthb = 10;
   Int_t firstBin = 0;
   Int_t axisPtcut = -1; //defined below from the weighted response
   
   Int_t MpMdptpptd[4] = {1, 0, 3, 2} ; //"fhnMassResponse"
   if(hname == "hResponse") {
      MpMdptpptd[0] = 0; MpMdptpptd[1] = 1; MpMdptpptd[2] = 2; MpMdptpptd[3] = 3;
   }
   if(hname == "fhnMassResponse") {
      MpMdptpptd[0] = 1; MpMdptpptd[1] = 0; MpMdptpptd[2] = 3; MpMdptpptd[3] = 2;
   }
   
   //histograms
   THnSparseF *hFinRespW[n+1];
   
   //Drawing histograms
   //only pTcut
   TCanvas *cDistributionsW = new TCanvas(Form("cDistributionsWpTCut%.0f", pTcut), Form("Distributions Weighted (#it{p}_{T} > %.0f GeV/#it{c}) ", pTcut), 1000, 1000);
   cDistributionsW->Divide(2,2);
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   //pt ranges of the results
   
   TCanvas **cDistrWPtBins = new TCanvas*[4];
   TString cvname[4] = {"cMassDetWPtBins", "cMassParWPtBins", "cPtDetWPtBins", "cPtParWPtBins"};
   TString cvtitle[4] = {"Mass Detector level Weighted", "Mass Particle level Weighted", "Pt Detector level Weighted", "Pt Particle level Weighted"};
  
   //legend
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.9);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   for(Int_t id = 0; id < 4; id++){
      
      cDistrWPtBins[id] = new TCanvas(Form("%spTCut%.0f", cvname[id].Data(), pTcut), Form("%s (#it{p}_{T} > %.0f GeV/#it{c})", cvtitle[id].Data(), pTcut), dx, dy);
      cDistrWPtBins[id]->Divide(nx,ny);
   }
   
   WeightClass *doweight = new WeightClass(Form("weight" ), Form("Weight "), npthb);
   
   
   doweight->SetOrderAxesRespose(MpMdptpptd[1], MpMdptpptd[0], MpMdptpptd[3], MpMdptpptd[2]);
   
   
   doweight->Initialise(inputFile, inputList, inputFile, inputList);
   
   //loop over pT-hard bins
   for(Int_t iipthb = firstBin; iipthb<npthb; iipthb++){
      Bool_t weighttaken = doweight->ReadWeight(iipthb);
      if(!weighttaken) {
      	 Printf("Failed reading weight, continue");
      	 continue;
      }
      doweight->SetResponse(iipthb, hname);
      Bool_t success = doweight->WeightResponse(iipthb);
      if(!success) Printf("Weight response returned %d", success);
      
      
   } // end loop over pT-hard bins
   for(Int_t i = 0; i < n+1; i++){
      //clean response within doweight to do the sum starting from a given pT hard bin
      doweight->CleanPointerResponseFinal();
      //starting loop when showing results
      if(i>0) doweight->SetFirstPtHard(iPtHb[i-1]);
      
      Bool_t res = doweight->Result();
      if(!res){
      	 Printf("Error! return");
      	 return;
      }
      THnSparseF* respTmp = doweight->GetResponseWFinal();
      axisPtcut = doweight->GetPtParAxis();
      Printf("Cut %s (axis n %d) > %f applied", respTmp->GetAxis(axisPtcut)->GetTitle(), axisPtcut, pTcut);
      Int_t binPtDmin = respTmp->GetAxis(axisPtcut)->FindBin(pTcut); // pT par cut
      respTmp->GetAxis(3)->SetRange(binPtDmin, -1);
      
      hFinRespW[i] = (THnSparseF*) respTmp->Clone(Form("hFinRespW_PtHB%d", iPtHb[i]));
      
      //if(i>0) hFinRespW[i]->SetName(Form("hFinRespW_PtHB%d", iPtHb[i]));
      
      TH1D** hpjRespFin = doweight->All1DProjections(respTmp, i, Form("hpjRespFin%d", i));
      
      if(i>0) leg->AddEntry(hpjRespFin[0], Form("p_{T, HARD} #geq %d ", iPtHb[i-1]), "L");
      else leg->AddEntry(hpjRespFin[0], Form("All p_{T, HARD}"), "L");
      
      for(Int_t id = 0; id < 4; id++){
      	 if(!hpjRespFin[id]){
      	    Printf("Projection %d not found, continue", id);
      	 }
      	 cDistributionsW->cd(id+1);
      	 
      	 if(i == 0) {
      	    gPad->SetLogy();
      	    hpjRespFin[id]->Draw();
      	 }
      	 else  hpjRespFin[id]->Draw("sames");
      }
      
   }// end loop on minimum pT hard bin used
   cDistributionsW->cd(1);
   leg->Draw();
   
   //plot in bins of the final analysis
   for(Int_t i = 0; i < n+1; i++){ //first is the unbiased
      for(Int_t ipt = 0; ipt < nptbins; ipt++){
      	 TPaveText *pvt = new TPaveText(0.3, 0.8, 0.7, 0.9);
      	 pvt->SetBorderSize(0);
      	 pvt->SetFillStyle(0);
      	 pvt->AddText(Form("%.0f < #it{p}_{T} < %.0f ", ptlims[ipt], ptlims[ipt+1]));
      	 
      	 Int_t binPtDRange[2] = {hFinRespW[i]->GetAxis(axisPtcut)->FindBin(ptlims[ipt]), hFinRespW[i]->GetAxis(axisPtcut)->FindBin(ptlims[ipt+1]-1)};
      	 
      	 hFinRespW[i]->GetAxis(axisPtcut)->SetRange(binPtDRange[0], binPtDRange[1]);
      	 WeightClass *docomparisononly = new WeightClass("weight", "Weight", 10);
      	 
      	 TH1D** hpjRespFinSel = docomparisononly->All1DProjections(hFinRespW[i], i, Form("hpjRespFinPTH%d_pTbin%d", i, ipt));
      	 
      	 for(Int_t id = 0; id < 4; id++){
      	    if(!hpjRespFinSel[id]){
      	       Printf("Projection %d not found, continue", id);
      	    }
      	    cDistrWPtBins[id]->cd(ipt+1);
      	    if(i == 0) {
      	       gPad->SetLogy();
      	       hpjRespFinSel[id]->Draw();
      	       pvt->Draw();
      	    }
      	    else  hpjRespFinSel[id]->Draw("sames");
      	    
      	    if(ipt == 0) {
      	       cDistrWPtBins[id]->cd(ipt+1);
      	       leg->Draw();
      	    }
      	 }
      }
      
   }
   for(Int_t id = 0; id < 4; id++) SaveCv(cDistrWPtBins[id]);
   
   SaveCv(cDistributionsW);
}

//________________________________________________________________________________________________
//________________________________________________________________________________________________

//main

//obsolete
void WeightResponse(TString inputFile, TString pathPythiaFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1067/outputpTHardBins/mergeRuns/", TString listPythiaFile = "PrepareInputForEmbedding"){
   
   Int_t firstBin = 0;
   Int_t offset = 1;
   const Int_t npthB = 10;
   //names of the histograms needed for the weights
   TString namehtrials = "fHistTrials", namehxsec = "fHistXsection", namehevents = "fHistEvents";
   
   TFile *fin = new TFile(inputFile);
   
   //variables describind the response
   Int_t        axisMp = 1, axisMd = 0, axisPtp = 3, axisPtd = 2;
   Int_t        axesResp[4] = {axisMd, axisMp, axisPtd, axisPtp};
   
   //Histograms
   
   //raw
   TH2D         *hptMpar[npthB];
   THnSparseF   *hResponse[npthB];
   THnSparseF   *hResponseN[npthB];
   //weighted
   
   TH2D         *hptMparW[npthB];
   THnSparseF   *hResponseW[npthB];
   
   //final
   THnSparseF   *hResponseWFinal = 0x0;
   TH2D         *hptMparWFinal = 0x0;
   TH2D         *hptMdetWFinal = 0x0;
   TH1D         *hFracPoverD = 0x0;
   TH1D         *hMCorrFactorP = 0x0;
   
   TString namehresponse = "fhnMassResponse_0";
   
   //QA
   
   TH1F *hCrossSec = new TH1F("hCrossSec", "Cross sections; p_{T}-hard bin; cross section", npthB, 0.5, npthB+0.5);
   TH1F *hTrials   = new TH1F("hTrials", "Number of trials; p_{T}-hard bin; N of trials"  , npthB, 0.5, npthB+0.5);
   
   //canvas
   // - per pt hard bin
   Int_t nx, ny, dx, dy;
   CalculatePads(npthB, nx, ny, dx, dy);
   
   TCanvas *cptMpart = new TCanvas("cptMpart", "pt vs M part", dx, dy);
   cptMpart->Divide(nx, ny);
   TCanvas *cptMpartW = new TCanvas("cptMpartW", "weighted pt vs M part", dx, dy);
   cptMpartW->Divide(nx, ny);
   
   // - per variable
   TCanvas *cProjResp[4];
   TCanvas *cProjRespWFinal[4];
   for(Int_t i=0; i<4; i++){
      cProjResp[i] = new TCanvas(Form("cProjResp%d", i), Form(" Response Projection %d", i));//, dx, dy);
      //cProjResp[i]->Divide(nx, ny);
      cProjRespWFinal[i] = new TCanvas(Form("cProjRespWFinal%d", i), Form("Response Projection Scaled and weighted dashed %d", i));
      
   }
   
   TLegend *legpTh = new TLegend(0.1, 0.4, 0.4, 0.7);
   legpTh->SetBorderSize(0);
   legpTh->SetFillStyle(0);
   
   // - summed up
   TCanvas *cptMWFinal = new TCanvas("cptMWFinal", "weighted Sum pt vs M ", 900, 500);
   cptMWFinal->Divide(2,1);
   
   Int_t nnx, nny, ddx, ddy;
   CalculatePads(nptbins, nnx, nny, ddx, ddy);
   
   TCanvas *cMWFinalpTjetDbins = new TCanvas("cMWFinalpTjetDbins", "Mass particle and detector level in detector level pT bins of the jet", ddx, ddy);
   cMWFinalpTjetDbins->Divide(nnx, nny);
   
   TCanvas *cMWFinalpTjetPbins = new TCanvas("cMWFinalpTjetPbins", "Mass particle and detector level in particle level pT bins of the jet", ddx, ddy);
   cMWFinalpTjetPbins->Divide(nnx, nny);
   
   TCanvas *cFracPoverD = new TCanvas("cFracPoverD", "Mass particle/detector", ddx, ddy);
   cFracPoverD->Divide(nnx, nny);
   
   // - QA
   TCanvas *cWeights = new TCanvas("cWeights", "Weights", 900, 500);
   cWeights->Divide(2,1);
   TCanvas *cWeightsInput = new TCanvas("cWeightsInput", "Weights Input", dx, dy);
   cWeightsInput->Divide(nx, ny);
   TCanvas *cTreeEntries = new TCanvas("cTreeEntries", "Entries of tree used per pT-hard", dx, dy);
   cTreeEntries->Divide(nx, ny);
   
   // Output file
   TFile *fout = new TFile("FullResponseOutput.root", "recreate");
   // Saving response, mass projections and correction factors
   
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   //loop over pT-hard bins
   for(Int_t ipthb = firstBin; ipthb<npthB; ipthb++){
      
      
      //retrieve histograms with weights for this pT hard bin
      TList* listPythia = ReadFile(Form("%s%02d/AnalysisResults.root", pathPythiaFile.Data(), ipthb+1), listPythiaFile);
      if(!listPythia){
      	 Printf("Cannot proceed with pT-hard %d, weight list not found", ipthb+1);
      	 //I don't even want to continue, because I miss one bin
      	 return;
      }
      
      TH1F*     htrials = (TH1F*)listPythia->FindObject(namehtrials);
      TProfile* hxsec = (TProfile*)listPythia->FindObject(namehxsec);
      TH1F*     hevents = (TH1F*)listPythia->FindObject(namehevents); //not sure if needed
      if(!htrials || !hxsec || !hevents){
      	 Printf("Histograms for weights not found, please check file in %s%02d/, I exit", pathPythiaFile.Data(), ipthb+1);
      	 listPythia->ls();
      }
      cWeightsInput->cd(ipthb+1);
      hxsec->Draw();
      
      Double_t xsec   = hxsec->GetBinContent(ipthb+1+offset);
      Double_t trials = htrials->GetBinContent(ipthb+1+offset);
      Double_t events = hevents->GetBinContent(ipthb+1+offset);
      Printf("Pt-hard %d xsec = %e, trials = %f", ipthb+1, xsec, trials);
      
      hCrossSec->Fill(ipthb+1, xsec);
      hTrials  ->Fill(ipthb+1, trials);
      
      Double_t scale = xsec/trials;
      
      //------------------------------------------------------------------------------------------
      // read list Embedding
      
      TString listnameEmb = "SingleTrackEmbedding";
      TList *lemb = (TList*) fin->Get(Form("%s%d",listnameEmb.Data(), ipthb+1));
      if(lemb){
      	 
      	 TH1F* htreeEntries = (TH1F*)lemb->FindObject("fhTreeEntriesUsed");
      	 cTreeEntries->cd(ipthb+1);
      	 
      	 if(htreeEntries) htreeEntries->Draw();
      	 
      }
      
      //------------------------------------------------------------------------------------------
      //read list of response to be weighted
      TString listnameInp = Form("JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH%d_pT0150_E_scheme_TC", ipthb+1);
      TList *list = (TList*) fin->Get(listnameInp);
      if(!list){
      	 Printf("Cannot proceed with pT-hard %d, response list %s not found", ipthb+1, listnameInp.Data());
      	 //I don't even want to continue, because I miss one bin
      	 fin->ls();
      
	 return;
      }
      THnSparseF *hnsp = (THnSparseF*)list->FindObject(namehresponse);
      if(!hnsp){
      	 Printf("Cannot proceed with pT-hard %d, %s not found", ipthb+1, namehresponse.Data());
      	 //I don't even want to continue, because I miss one bin
      	 list->ls();
      	 return;
      }
      //TH1D *htest = hnsp->Projection(axisMp);
      //htest->SetName(Form("htest%d", ipthb));
      
      hResponse[ipthb] = (THnSparseF*)hnsp->ProjectionND(4, axesResp, "E");
      hResponse[ipthb]->SetName(Form("hResponse_pTh%d", ipthb+1));
      hResponse[ipthb]->Sumw2();
      
      hptMpar[ipthb]   = hResponse[ipthb]->Projection(axisMp, axisPtp);
      hptMpar[ipthb]->SetName(Form("hptMpar_pTh%d", ipthb+1));
      TH1D *hpT = hptMpar[ipthb]->ProjectionY(Form("hpT%d", ipthb+1));
      //hptMpar[ipthb]->SetLineColor(colors[ipthb]);
      
      cptMpart->cd(ipthb+1);
      hptMpar[ipthb]->Draw("colz");
      //htest->Draw("E");
      //hpT->Draw("Esames");
      
      Double_t integralptMpar = hptMpar[ipthb]->Integral();

      hResponseN[ipthb] = (THnSparseF*)hResponse[ipthb]->Clone(Form("hResponseN_pTh%d", ipthb+1));
      hResponseN[ipthb]->Scale(1./integralptMpar);
      
      //draw each projection in each pt hard bin (Pad)
      TH1D** hRespProj = All1DProjections(hResponse[ipthb], ipthb);
      TH1D** hRespNProj = All1DProjections(hResponseN[ipthb], ipthb);
      Bool_t bfirst = kFALSE;
      for(Int_t i=0; i<4; i++){
      	 //Printf("Drawing projection %d of pt hard %d", i, ipthb+1);
      	 
      	 if(i == 0) legpTh->AddEntry(hRespProj[i], Form("pThard %d", ipthb+1), "L");
      	 if(ipthb == firstBin) {
      	    cProjResp[i]->cd();
      	    gPad->SetLogy();
      	    hRespProj[i]->Draw("E");
      	    
      	    cProjRespWFinal[i]->cd();
      	    gPad->SetLogy();
      	    hRespNProj[i]->Draw("E");
      	 }
      	 else {
      	    cProjResp[i]->cd();
      	    hRespProj[i]->Draw("Esames");
      	    cProjRespWFinal[i]->cd();
      	    hRespNProj[i]->Draw("Esames");
      	 }

      }
      
      
     //------------------------------------------------------------------------------------------ 
     
     //Apply weights
     
     hptMparW[ipthb] = (TH2D*)hptMpar[ipthb]->Clone(Form("hptMparW_pTh%d", ipthb+1));
     hptMparW[ipthb]->Sumw2();
     hptMparW[ipthb]->Scale(scale);
     Double_t integralptMparW = hptMparW[ipthb]->Integral();
     cptMpartW->cd(ipthb+1);
     hptMparW[ipthb]->Draw("colz");
      
     hResponseW[ipthb] = (THnSparseF*)hResponse[ipthb]->Clone(Form("hResponseW_pTh%d", ipthb+1));
     hResponseW[ipthb]->Sumw2();
     hResponseW[ipthb]->Scale(scale/integralptMparW);
     
     if(!hResponseWFinal) {
     	hResponseWFinal = (THnSparseF*)hResponseW[ipthb]->Clone("hResponseW");
     	hResponseWFinal->Sumw2();
     } else {
     	hResponseWFinal->Add(hResponseW[ipthb]);
     }
     
     
   } // end loop on pt-hard bins
   
   //save final response THnSparse
   
   fout->cd();
   hResponseWFinal->Write();
   
   //projections of the final thnsparse
   
   //draw each projection in each pt hard bin (Pad) apply a cut that selects only the region of pT considered later
   Int_t binDAR[2] = {hResponseWFinal->GetAxis(axisPtd)->FindBin(ptlims[0]), hResponseWFinal->GetAxis(axisPtd)->FindBin(ptlims[nptbins])};
      
   hResponseWFinal->GetAxis(axisPtd)->SetRange(binDAR[0], binDAR[1]);
   TH1D** hRespWFinalProj = All1DProjections(hResponseWFinal, 0, "hRespWProj");
   
   for(Int_t i=0; i<4; i++){
      //Printf("Drawing projection %d of pt hard %d", i, ipthb+1);
      cProjRespWFinal[i]->cd();
      hRespWFinalProj[i]->SetLineStyle(2);
      hRespWFinalProj[i]->Draw("Esames");
      
   }
   
   fout->cd();
   hRespWFinalProj[axisPtd]->Write();
   hRespWFinalProj[axisPtp]->Write();
   //2D projections
   
   hptMparWFinal = hResponseWFinal->Projection(axisMp,axisPtp);
   
   hptMdetWFinal = hResponseWFinal->Projection(axisMd,axisPtd);
   
   // ptjet bins
   
   //TPaveText **pvpt = GetPavePtBins(nptbins, ptlims);
   
   //TPaveText** pvpt = new TPaveText*[nptbins];
   
   Int_t x1 = 0.3, y1 = 0.8, x2 = 0.8, y2 = 0.9;
   
   //for(Int_t ip = 0; ip<nptbins; ip++){
   //   pvpt[ip] = new TPaveText(x1, y1, x2, y2, "NDC");
   //   pvpt[ip]->SetBorderSize(0);
   //   pvpt[ip]->SetFillStyle(0);
   //   pvpt[ip]->AddText(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", ptlims[ip], ptlims[ip+1]));
   //   Printf("pave %p", pvpt[ip]);
   //}
   TLegend *legM = new TLegend(0.6, 0.6, 0.9, 0.8);
   legM->SetFillStyle(0);
   legM->SetBorderSize(0);
   
   Double_t errrelmax = 0.08, maxrange = 1, maxrangecut = 4;
   for(Int_t ipt = 0; ipt < nptbins; ipt++){
      TPaveText *pvpt = new TPaveText(x1, y1, x2, y2, "NDC");
      pvpt->SetName(Form("pvpt%d", ipt));
      pvpt->AddText(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", ptlims[ipt], ptlims[ipt+1]));
      // detector level pt selection
      Int_t binD[2] = {hResponseWFinal->GetAxis(axisPtd)->FindBin(ptlims[ipt]), hResponseWFinal->GetAxis(axisPtd)->FindBin(ptlims[ipt+1]) -1};
      
      hResponseWFinal->GetAxis(axisPtd)->SetRange(binD[0], binD[1]);
      TH1D** hRespWFinalProjPtD = All1DProjections(hResponseWFinal, 0, Form("hpTD%d", ipt));
      cMWFinalpTjetDbins->cd(ipt+1);
      hFracPoverD = (TH1D*)hRespWFinalProjPtD[axisMp]->Clone(Form("hFracPoverD_pt%d", ipt));
      hFracPoverD->Divide(hRespWFinalProjPtD[axisMd]);
      
      hRespWFinalProjPtD[axisMd]->SetLineColor(colors[2]); //detector = red
      hRespWFinalProjPtD[axisMp]->SetLineColor(colors[1]); //particle = blue
      if(ipt == 0){
      	 legM->AddEntry(hRespWFinalProjPtD[axisMd], "Detector", "L");
      	 legM->AddEntry(hRespWFinalProjPtD[axisMp], "Particle", "L");
      }
      hRespWFinalProjPtD[axisMd]->GetXaxis()->SetRangeUser(0., 25.);
      
      hRespWFinalProjPtD[axisMd]->Draw("E");
      hRespWFinalProjPtD[axisMp]->Draw("Esames");
      
      //pvpt[ipt]->DrawClone();
      pvpt->Draw();
      legM->Draw();
      for(Int_t ibin = 0; ibin < hFracPoverD->GetNbinsX(); ibin++){
      	 Double_t cont = hFracPoverD->GetBinContent(ibin+1), err = hFracPoverD->GetBinError(ibin+1);
      	 maxrange = hFracPoverD->GetBinCenter(ibin+1);
      	 Printf("Max range = %f, err = %f, cont = %f, err/cont = %f", maxrange, err, cont, err/cont);
      	 if(err/cont > errrelmax && maxrange > maxrangecut) { //the second condition is required for a specific case where the first bins have very large errors
      	    maxrange = hFracPoverD->GetBinCenter(ibin);
      	    break;
      	 }
      }
      Printf("Set range 0 - %f", maxrange);
      hFracPoverD->GetXaxis()->SetRangeUser(0., maxrange);
      cFracPoverD->cd(ipt+1);
      gPad->SetLogy();
      hFracPoverD->Draw();
      
      //save mass projections binned at detector level and their ratio
      fout->cd();
      hFracPoverD->Write();
      hRespWFinalProjPtD[axisMd]->Write();
      hRespWFinalProjPtD[axisMp]->Write();
      hRespWFinalProjPtD[axisPtd]->Write();
      hRespWFinalProjPtD[axisPtp]->Write();
      //reset
      hResponseWFinal->GetAxis(axisPtd)->SetRange(0, -1);
      
      // particle level pt selection
      Int_t binP[2] = {hResponseWFinal->GetAxis(axisPtp)->FindBin(ptlims[ipt]), hResponseWFinal->GetAxis(axisPtp)->FindBin(ptlims[ipt+1]) -1};
      
      hResponseWFinal->GetAxis(axisPtp)->SetRange(binP[0], binP[1]);
      TH1D** hRespWFinalProjPtP = All1DProjections(hResponseWFinal, 0, Form("hpTP%d", ipt));
      cMWFinalpTjetPbins->cd(ipt+1);
      hRespWFinalProjPtP[axisMd]->SetLineColor(colors[2]); //detector = red
      hRespWFinalProjPtP[axisMp]->SetLineColor(colors[1]); //particle = blue
      hRespWFinalProjPtP[axisMd]->GetXaxis()->SetRangeUser(0., 25.);
      hRespWFinalProjPtP[axisMd]->Draw("E");
      hRespWFinalProjPtP[axisMp]->Draw("Esames");
      //pvpt[ipt]->Draw();
      pvpt->Draw();
      legM->Draw();
      
      //reset
      hResponseWFinal->GetAxis(axisPtp)->SetRange(0, -1);
   }
   
   //drawing the results
   cptMWFinal->cd(1);
   hptMparWFinal->Draw("colz");
   cptMWFinal->cd(2);
   hptMdetWFinal->Draw("colz");
   
   //drawing QA
   cWeights->cd(1);
   hCrossSec->Draw();
   cWeights->cd(2);
   hTrials->Draw();
   
   //refining canvas
   for(Int_t i=0; i<4; i++) {
      cProjResp[i]->cd();
      legpTh->Draw();
      //save Canvas
      SaveCv(cProjRespWFinal[i]);
      SaveCv(cProjResp[i]);
   }
   
   //save Canvas
   SaveCv(cptMWFinal);
   SaveCv(cWeights);
   SaveCv(cWeightsInput);
   SaveCv(cptMpart);
   SaveCv(cptMpartW);
   SaveCv(cMWFinalpTjetDbins);
   SaveCv(cMWFinalpTjetPbins);
   SaveCv(cFracPoverD);
}

//________________________________________________________________________________________________

// read the tree and prepare response
void PrepareResponseFromTree(TString pathPythiaFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/", TString treeFileName = "AnalysisResults.root", TString treename = "fTreeJet", Double_t ptDcut = -1, Double_t ptPcut = 20, Int_t firstBin = 0){
   
   /// ptDcut, ptPcut :  if -1 not applied, if > -1 applied at det/par level
   
   //define THnSparse to be filled
   const Int_t dim = 4;
   const Int_t nBinsfineM     = 200;
   const Int_t nBinsfinePt    = 200;
   const Double_t minPt       = 0.;
   const Double_t maxPt       = 200.;
   const Double_t minM        = 0.;
   const Double_t maxM        = 50.;
   Int_t npthb =10;
   Int_t nBins[dim] = {nBinsfineM, nBinsfineM, nBinsfinePt, nBinsfinePt};
   Double_t xmin[dim]  = {minM, minM, minPt, minPt};
   Double_t xmax[dim]  = {maxM, maxM, maxPt, maxPt};
   TString hsptitle = "Mass-pT response; #it{M}_{part}; #it{M}_{det};#it{p}_{T,part}; #it{p}_{T,det}";
   TString h2dtitle = "pT response;  #it{p}_{T,det};#it{p}_{T,part}";
   
   THnSparseF **fhResponse = new THnSparseF*[npthb];
   TH2D       **fh2Response = new TH2D*[npthb];
   
   TString suff = "";
   if(ptPcut > -1) suff = Form("CutP%.0f", ptPcut);
   if(ptDcut > -1) suff = Form("CutD%.0f", ptDcut);
   
   TFile *fout = new TFile(Form("ResponseFromTree%s.root", suff.Data()), "recreate");
   
   
   
   for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
      //   for(Int_t ipthb = firstBin; ipthb<2; ipthb++){
      fhResponse[ipthb] = new THnSparseF(Form("hResponse%d", ipthb+1), hsptitle.Data(), dim, nBins, xmin, xmax);
      fh2Response[ipthb] = new TH2D(Form("h2Response%d", ipthb+1), h2dtitle.Data(), nBinsfinePt, minPt, maxPt, nBinsfinePt, minPt, maxPt);
      
      TString inputfile = "";
      if(treeFileName.Contains(".root")) inputfile = Form("%s%02d/%s", pathPythiaFile.Data(), ipthb+1, treeFileName.Data());
      else inputfile = Form("%s%02d/%s%d.root", pathPythiaFile.Data(), ipthb+1, treeFileName.Data(),  ipthb+1);
      
      TTree* tree = ReadTree (inputfile, treename);
      if(!tree) return;
      TString branches[2] = {"fJetDet.", "fJetPart."};
      TLorentzVector *vecsP = 0x0, *vecsD = 0x0;
      
      tree->SetBranchAddress(branches[0], &vecsD);
      tree->SetBranchAddress(branches[1], &vecsP);
      Printf(" -- ipthb = %d", ipthb);
      Int_t nentries = tree->GetEntries();
      Printf("N entries tree %d", nentries);
      
      for(Int_t i = 0; i<nentries; i++){
      	 
      	 Int_t ient = tree->GetEntry(i);
      	 //if (i % 100000 == 0) Printf("Entry %d: This is pt hard %d", i, ipthb);
      	 Double_t ptgen, ptrec, mrec, mgen;
      	 
      	 ptrec = vecsD->Pt();
      	 mrec  = vecsD->M();
      	 ptgen = vecsP->Pt();
      	 mgen  = vecsP->M();
      	 
      	 if(ptPcut > -1 && ptgen <= ptPcut) continue;
      	 if(ptDcut > -1 && ptrec <= ptDcut) continue;
      	 //if(applyptcut && ptrec <= ptlims[0] && ptrec >= ptlims[nptbins]) continue; //obsolete. Used in /data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/ResponseEntriesUsed
      	 
      	 Double_t response[dim] = {mgen, mrec, ptgen, ptrec};
      	 
      	 fhResponse[ipthb]->Fill(response);
      	 fh2Response[ipthb]->Fill(ptrec, ptgen);
      }
      fout->cd();
      fhResponse[ipthb]->Write();
      fh2Response[ipthb]->Write();
   }
}

//________________________________________________________________________________________________

//-> draw and save weighted histograms
void WeightedTree(Bool_t donorm = kFALSE, TString pathPythiaFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/", TString listnameInp = "PrepareInputForEmbedding", TString treename = "", TString filename = "", TString hname = "fhnMassResponse", TString pathDataFile = "", Int_t trueFrom4vecEmb = 0, Int_t ptHardFirst = 1, Int_t firstBin = 0, Bool_t fineBin = kFALSE, Bool_t weight4Rho = kFALSE, TString pathDataRho = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/AnalysisResults.root", TString listDataRho = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv", TString histoDataRho = "fhnRhoVsRhoMVsLeadJetPtVsMassVsCent"){
   /// Four uses : 
   /// 1) read the thnsparses filled with the tree entries in PrepareResponseFromTree
   ///    Input:
   ///          - pathPythiaFile: path pointing to the base dir for the bin-by-bin output of the task that contains the xsec and ntrials histograms
   ///          - listnameInp :  name of the list where to read xsec and ntrials
   ///          - treename: it's the suffix of the file "ResponseFromTree<suffix>.root"
   ///          - filename: ""
   /// 2) read the THnsparses with the response from a TList 
   ///    Input:
   ///          - pathPythiaFile: path pointing to the base dir for the bin-by-bin output of the task that contains the xsec and ntrials histograms
   ///          - listnameInp :  name of the list where to read xsec and ntrials and the input response
   ///          - treename: ""
   ///          - filename: name of the file where the input response is
   ///          - hname : name of the THnSparse. Either "fhnMassResponse" or "hResponse"
   ///          - pathDataFile = ""
   /// 3) weights and input data are in different files (this is the case for the embedding)
   ///          - pathPythiaFile: path pointing to the base dir for the bin-by-bin output of the task that contains the xsec and ntrials histograms (the PrepareInputForEmbedding used for the embedding)
   ///          - listnameInp :  firt part of the name of the list where to read the input response (before the string "PtH", e.g. "JetShapeConst_JetRhosub_AKTChargedR040_PicoTracks"). The list name for the weights is hard coded to "PrepareInputForEmbedding" when !pathDataFile.IsNull()
   ///          - treename: ""
   ///          - filename: name of the file where the input response is
   ///          - hname : name of the THnSparse. Should be "fhnMassResponse_0"
   ///          - pathDataFile : path of the output file for the response that contains all the pT hard bins
   /// 4) change the axes of the detector level quantities to the detector level embedded. This is available only for JetShape task tag>= 20160227. It's a cross check of the detector level only
   ///          - trueFrom4vecEmb = 1	
   /// For later tags >=20160730 (check) available also the non-bkg-subtracted M, pT 
   ///          - trueFrom4vecEmb = 2
   ///
   /// ptHardFirst = <number> : start merging the output from pThard = <number> . Counter starts from 1. In must be > than firstBin
   /// firstBin: first pT hard bin available in the output. Counter starts from 0
   /// The histogram "fhnDeltaMass_0" contains dM (axis 0) dpT (axis 1). Axis 2 and 4 contain M and pT at det (or det+fluc) level, Axis 3 and 5 contain M and pT at particle level  
   /// The histogram "fhnDeltaMassAndBkgInfo" contains dM, dpT, msub, munsub, ptsub, ptunsub, rho, rhom
   /// weight4Rho = when true the weighting according to rho it's performed. It requires additional input:
   /// - pathDataRho = path and filename of the data output containing the THnSparse with the rho vs rho_m vs pt vs M vs centrality
   /// - listDataRho = listname within the above file
   /// - histoDataRho = name of the thnsparse, see default
   
   TStopwatch watch;
   watch.Start();
   
   Int_t npthb = 10;//2 ;//3;
   
   Printf("INFO : first pT hard bin read = %d", firstBin);
   if(trueFrom4vecEmb == 1) Printf("INFO : Detector effect only response matrix, taken from embedding job");
   if(trueFrom4vecEmb == 2) Printf("INFO : No background subtracted response matrix, taken from embedding job");
   
   THnSparseF **fhResponse = new THnSparseF*[npthb];
   TFile *fin = 0x0;
   
   if(!treename.IsNull()){
      fin = new TFile(Form("ResponseFromTree%s.root", treename.Data()));
      if(!fin->IsOpen()){
      	 Printf("Input file not found");
      	 return;
      }
   }
   
   Int_t MpMdptpptd[4] = {1, 0, 3, 2} ; //"fhnMassResponse", "fhnMassResponse_0"
   Int_t axesEmbMPt[2] = {-1, -1};
   Int_t rhos[2] = {6, 7};
   if(0) {
   	   Printf("Using rho and rho m");
   MpMdptpptd[0] = rhos[0];
   MpMdptpptd[1] = rhos[1];
   }
   if(!weight4Rho) { // this needed in case of older output where rho and rhom were not present
   	   rhos[0] = -1; 
   	   rhos[1] = -1;
   }
   
   if(hname == "hResponse" || filename.IsNull()) {
      MpMdptpptd[0] = 0; MpMdptpptd[1] = 1; MpMdptpptd[2] = 2; MpMdptpptd[3] = 3;
   }
   if(trueFrom4vecEmb == 1 ) {
      MpMdptpptd[0] = 1; MpMdptpptd[1] = 4; MpMdptpptd[2] = 3; MpMdptpptd[3] = 5;
   }
   
   if(trueFrom4vecEmb == 2) {
      MpMdptpptd[0] = 1; MpMdptpptd[1] = 8; MpMdptpptd[2] = 3; MpMdptpptd[3] = 9;
   }
   
   if(hname == "fhnMassResponse_0" && !trueFrom4vecEmb){
      //Add also the det level embedded axes, in older version it's not available and will crash
      //axesEmbMPt[0] = 4;
      //axesEmbMPt[1] = 5;
   	 
   }
   if(hname == "fhnDeltaMass_0"){
   	   MpMdptpptd[0] = 1; MpMdptpptd[1] = 0; MpMdptpptd[2] = 5; MpMdptpptd[3] = 3;
   }
   
   if(hname == "fhnDeltaMassAndBkgInfo"){
   	   MpMdptpptd[0] = 6; MpMdptpptd[1] = 7; MpMdptpptd[2] = 2; MpMdptpptd[3] = 4;
   	   // 6, 7 rho, rhoM
   	   // 2, 4 M, pt sub
   	   // 3, 5 M, pt unsub
   }
   WeightClass *doweight = new WeightClass("weight", "Weight", npthb);
   
   if(pathDataFile.IsNull()) doweight->Initialise(pathPythiaFile, listnameInp, filename, listnameInp);
   else doweight->Initialise(pathPythiaFile, "PrepareInputForEmbedding", Form("%s", pathDataFile.Data()), listnameInp);
   
   TCanvas *cProjResp = new TCanvas(Form("cProjResp%s", trueFrom4vecEmb ? "4VecDet" : ""), Form("Response before weighting "), 1000, 1000);
   cProjResp->Divide(2,2);
   TCanvas *cProjRespW = new TCanvas(Form("cProjRespW%s", trueFrom4vecEmb ? "4VecDet" : ""), Form("Response after weighting "), 1000, 1000);
   cProjRespW->Divide(2,2);
   
   TCanvas *cProjRespWN = new TCanvas(Form("cProjRespWN%s", trueFrom4vecEmb ? "4VecDet" : ""), Form("Response after weighting and norm"), 1000, 1000);
   cProjRespWN->Divide(2,2);
   
   TCanvas *cProjPartForNorm = new TCanvas(Form("cProjPartForNorm%s", trueFrom4vecEmb ? "4VecDet" : ""), "Particle level mass and pt for normalization", 1000, 600);
   cProjPartForNorm->Divide(3,1);
   
   TCanvas *cProjPartNormed = new TCanvas(Form("cProjPartNormed%s", trueFrom4vecEmb ? "4VecDet" : ""), "Particle level mass and pt normalized", 1000, 600);
   cProjPartNormed->Divide(3,1);
   
   doweight->SetDoNormalisation(donorm);
   doweight->SetFineFinning(fineBin);
   
   if(!pathDataFile.IsNull()) doweight->ReadNumberOfJetPerPtHBin(Form("%s/NJetsPerPtHBin.root", pathPythiaFile.Data()));
   
   // settings for the rho weighting
   if(weight4Rho){
   	   doweight->SetDoWeightForBkg(kTRUE);
   	   doweight->SetDataPathForBkd(pathDataRho);
   	   doweight->SetDataListForBkdName(listDataRho);
   	   doweight->SetDataHistForBkdName(histoDataRho);
   }
   
   TLegend *legpTh = new TLegend(0.1, 0.4, 0.4, 0.7);
   legpTh->SetBorderSize(0);
   legpTh->SetFillStyle(0);
   
   //loop on pT hard bins
   for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
   //     for(Int_t ipthb = firstBin; ipthb<2; ipthb++){
      //TString listnameInp = Form("JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH%d_pT0150_E_scheme_TC", ipthb+1);
      
      // setorder has not to go in each pT hard bin loop because the order gets canged each time!!
      doweight->SetOrderAxesRespose(MpMdptpptd[0], MpMdptpptd[1], MpMdptpptd[2], MpMdptpptd[3]);
      doweight->SetOrderAxesMPt3D(MpMdptpptd[0], MpMdptpptd[1], MpMdptpptd[2], MpMdptpptd[3], axesEmbMPt[0], axesEmbMPt[1]);
      if(weight4Rho) doweight->SetOrderAxesResposeWithRho(MpMdptpptd[0], MpMdptpptd[1], MpMdptpptd[2], MpMdptpptd[3], rhos[0], rhos[1]);
      
      if(fin){ //case 1
      	 fhResponse[ipthb] = (THnSparseF*)fin->Get(Form("hResponse%d", ipthb+1));
      	 doweight->ReadInputData(ipthb, fhResponse[ipthb]);
      } else { //case 2
      	 Bool_t set = doweight->SetResponse(ipthb, hname);
      	 if(!set) return;
      	 else Printf("Response successfully set");
      }
      
      Bool_t weighttaken = doweight->ReadWeight(ipthb);
      if(!weighttaken) {
      	 Printf("Failed reading weight, continue");
      	 continue;
      } else Printf("Weight successfully set");
      
      Bool_t success = doweight->WeightResponse(ipthb);
      if(!success) Printf("Weight response returned %d", success);
       else Printf("Weighting procedure was successful");
      
      
      //drawing input
      TH1D**  hprojResp  = doweight->All1DResponseProj(ipthb, ipthb);
      TH1D**  hprojRespW = doweight->All1DResponseWProj(ipthb, ipthb);
      TH1D**  hprojRespWN = doweight->All1DResponseWNormProj(ipthb, ipthb);
      
      for(Int_t ipj = 0; ipj<4; ipj++){
      	 
      	 if(ipj == 0) legpTh->AddEntry(hprojResp[ipj], Form("pThard %d", ipthb+1), "L");
      	 
      	 cProjResp->cd(ipj+1);
      	 gPad->SetLogy();
      	 
      	 if(ipthb == firstBin) hprojResp[ipj]->Draw();
      	 else hprojResp[ipj]->Draw("sames");
      	 
      	 cProjRespW->cd(ipj+1);
      	 gPad->SetLogy();
      	 
      	 hprojRespW[ipj]->SetLineStyle(2);
      	 if(ipthb == firstBin) hprojRespW[ipj]->Draw();
      	 else hprojRespW[ipj]->Draw("sames");
      	 if(donorm){
      	    cProjRespWN->cd(ipj+1);
      	    gPad->SetLogy();
      	    hprojRespWN[ipj]->SetLineStyle(3);
      	    if(ipthb == firstBin) hprojRespWN[ipj]->Draw();
      	    else hprojRespWN[ipj]->Draw("sames");
      	 }
      }
      
      TH2D* h2DPartProj  = (TH2D*)doweight->Get2DPartProjPtHbin(ipthb);
      TH1D* hPtPartProj  = h2DPartProj->ProjectionX();
      hPtPartProj->SetLineColor(colors[ipthb]);
      TH1D* hMPartProj   = h2DPartProj->ProjectionY();
      hMPartProj->SetLineColor(colors[ipthb]);
   

      TH2D* h2DPartProjW = (TH2D*)doweight->Get2DPartProjWPtHbin(ipthb);
      TH1D* hPtPartProjW = h2DPartProjW->ProjectionX();
      hPtPartProjW->SetLineColor(colors[ipthb]);
      hPtPartProjW->SetLineStyle(2);
      TH1D* hMPartProjW  = h2DPartProjW->ProjectionY();
      hMPartProjW->SetLineColor(colors[ipthb]);
      hMPartProjW->SetLineStyle(2);
      
      //draw the weighted histograms before dividing by N
      cProjPartNormed->cd(2);
      gPad->SetLogy();
      if(ipthb == firstBin) hPtPartProjW->Draw();
      else hPtPartProjW->Draw("sames");
      cProjPartNormed->cd(3);
      gPad->SetLogy();
      if(ipthb == firstBin) hMPartProjW->Draw();
      else hMPartProjW->Draw("sames");
      
      // draw the original histograms before any weighting
      cProjPartForNorm->cd(2);
      gPad->SetLogy();
      if(ipthb == firstBin) hPtPartProj->Draw();
      else hPtPartProj->Draw("sames");
      cProjPartForNorm->cd(3);
      gPad->SetLogy();
      if(ipthb == firstBin) hMPartProj->Draw();
      else hMPartProj->Draw("sames");
      
   } // end loop on pT hard bins
  
   cProjResp->cd(1);
   legpTh->Draw();
   cProjPartNormed->cd(2);
   legpTh->Draw();
   cProjPartForNorm->cd(2);
   legpTh->Draw();
   
   if(!donorm && !filename.Contains(".root")){
      doweight->SaveNumberOfJetPerPtHBin(Form("%s/NJetsPerPtHBin.root", pathPythiaFile.Data()));
   }
   if(ptHardFirst > 0) doweight->SetFirstPtHard(ptHardFirst);
   
   Bool_t results = doweight->Result();
   if(!results) Printf("The sum of the pT hard bins failed");
   else Printf("Final response successfully obtained");
   
   THnSparseF*  fResponsefinal = doweight->GetResponseWFinal();
   THnSparseF*  fMPt3Dfinal = doweight->GetMPt3DWFinal();
   
   TH1D**  hprojRespFin = doweight->All1DResponseWFinalProj(1);
   
   TCanvas *cRho = new TCanvas("crho", "Rho, several bins", 600, 600);
   TCanvas *cRhom = new TCanvas("crhom", "Rho_m, several bins", 600, 600);
   
   if(weight4Rho){
   	   Bool_t res = doweight->ReadRhoFromData();
   	   if(res){
   	   	   TList *listRhoH = doweight->GetListOfRhoAndRhoMProjections();
   	   	   if(listRhoH){
   	   	   	   for(Int_t il = 0; il<listRhoH->GetEntries(); il+=2){
   	   	   	   	   TH1D *hrho = (TH1D*)listRhoH->At(il);
   	   	   	   	   TH1D *hrhom= (TH1D*)listRhoH->At(il);
   	   	   	   	   
   	   	   	   	   cRho->cd();
   	   	   	   	   hrho->Draw("sames");
   	   	   	   	   cRhom->cd();
   	   	   	   	   hrhom->Draw("sames");
   	   	   	   }
   	   	   }
   	   	   
   	   	 doweight->  WeightForRho();
   	   }
   }
   
   //doweight->NormalizePerIntegral();

   TCanvas *cProjRespFinal = new TCanvas(Form("cProjRespFinal%s", trueFrom4vecEmb ? "4VecDet" : ""), "Final response", 1000, 1000);
   cProjRespFinal->Divide(2,2);
   for(Int_t ipj = 0; ipj<4; ipj++){
   	   if(!hprojRespFin) break;
      cProjRespFinal->cd(ipj+1);
      gPad->SetLogy();
      if(hprojRespFin[ipj]) hprojRespFin[ipj]->Draw();
      else Printf("Error! response projections not produced");
   }
  
   SaveCv(doweight->DrawScaleF() );
   SaveCv(doweight->DrawNormInt());
   SaveCv(doweight->DrawCrossSec());
   SaveCv(doweight->DrawNTrials());
   
   SaveCv(cProjRespFinal);
   SaveCv(cProjResp);
   SaveCv(cProjRespW);
   SaveCv(cProjPartForNorm);
   SaveCv(cProjPartNormed);

   Int_t rebinall[4] = {4,4,10,10};
   THnSparseF *fResponsefinalReb = fResponsefinal;
   //if(!donorm){
   //  fResponsefinalReb = (THnSparseF*)fResponsefinal->Rebin(rebinall);
   //}
   
   TH1::AddDirectory(kFALSE);
   
   /*
   Int_t rebinall1[5] = {4,4,10,10,1};
   TList *l576 = ReadFile("/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged");
   
   THnSparseF *hsp576 = (THnSparseF*)l576->FindObject("fhnMassResponse");
   hsp576->SetName("hResp576Merged");
   //
   //hsp576 = (THnSparseF*)hsp576->Rebin(rebinall1);
   SaveCv(doweight->CompareResponse(hsp576));
   
   TString filename1 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/AnalysisResults.root", list1 = "PrepareInputForEmbedding";
   TList *lPrep = ReadFile(filename1, list1);
   if(!lPrep) return;
   THnSparseF *hPrep = (THnSparseF*)lPrep->FindObject("hResponse");
   hPrep->SetName("hRespPrepMerged");
   //rebinall[0] = 10;
   //rebinall[1] = 10;
   //hPrep = (THnSparseF*)hPrep->Rebin(rebinall);
   SaveCv(doweight->CompareResponses(fResponsefinal, hPrep, 1 , 0, Form("4D Weighted Response %s", listnameInp.Data()), "Merged with Weights 1087"));
   
   // compare with 576 weighted by WeightClass
   TString fname576W = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root";
   
   TFile* f576W = new TFile(fname576W);
   if(f576W->IsOpen()){
      THnSparseF *hResp576W = (THnSparseF*)f576W->Get("fhResponseFinal");
      if(!hResp576W){
      	 Printf("Crashing....");
      	 f576W->ls();
      }
      hResp576W->SetName("hResp576W");
      SaveCv(doweight->CompareResponses(fResponsefinal, hResp576W, 1 , 3, Form("4D Weighted Response %s", listnameInp.Data()), "576 Weighted Response"));
   
   } else Printf("%s not found", fname576W.Data());
   */
   TString suff = ""; 
   if(!filename.IsNull()) suff = listnameInp;
  TFile *fout = new TFile(TString::Format("ResponseW%s.root", suff.Data()), "recreate");
  fout->cd();
  if(!fResponsefinal) Printf("ERROR: Response final doesn't exist and can't be saved.");
  else fResponsefinal->Write();
  if(fMPt3Dfinal) fMPt3Dfinal->Write();
  fout->Close();
  //fResponsefinN->Write();
   
  watch.Stop();
  watch.Print();
}

//________________________________________________________________________________________________
// comparison of weighted responses

void Comparison(TString filename1 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/AnalysisResults.root", TString list1 = "PrepareInputForEmbedding", TString filename2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString list2 = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString filename3 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1096/outputpTHardBins/mergeRuns/AnalysisResults.root", TString list3 = "JetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme"){

   //Train with response from task TreePreparation (pTdet > 2 GeV/c)
   TList *lPrep = ReadFile(filename1, list1);
   if(!lPrep) return;
   THnSparseF *hPrep = (THnSparseF*)lPrep->FindObject("hResponse");
   
   // Original response from Marta
   TList *l576 = ReadFile(filename2, list2);
   if(!l576) return;
   THnSparseF *hsp576 = (THnSparseF*)l576->FindObject("fhnMassResponse");
   
   // Train 1096 Rerun train 576 (no bkg sub)
   TList *l1096 = ReadFile(filename3, list3);
   if(!l1096) return;
   THnSparseF *hsp1096 = (THnSparseF*)l1096->FindObject("fhnMassResponse");
   hsp1096->SetName("fhnMassResponse1096");
   
   
   //rebin
   Int_t rebinall[5] = {10,10,10,10,1};
   //hPrep = (THnSparseF*)hPrep->Rebin(rebinall);
   
   Int_t rebinall1[5] = {4,4,10,10,1};
   //hsp576 = (THnSparseF*)hsp576->Rebin(rebinall1);
   //hsp1096 = (THnSparseF*)hsp1096->Rebin(rebinall);
   Int_t select[2] = {hsp576->GetAxis(3)->FindBin(2.), hsp576->GetAxis(3)->GetNbins()};
   
   hsp576->GetAxis(3)->SetRange(select[0], select[1]);
   hsp1096->GetAxis(3)->SetRange(select[0], select[1]);
   
   WeightClass *docomparisononly = new WeightClass("weight", "Weight", 10);
   SaveCv(docomparisononly->CompareResponses(hPrep, hsp576, 1, 2, "From Tree", "Train576"));
   SaveCv(docomparisononly->CompareResponses(hsp1096, hsp576, 3, 2, "Train1096", "Train576"));
   SaveCv(docomparisononly->CompareResponses(hPrep, hsp1096, 1, 3, "From Tree", "Train1096"));
   
}

//________________________________________________________________________________________________
// comparison of weighted responses

void ComparisonDavideFilip(TString filename1 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/AnalysisResults.root", TString list1 = "PrepareInputForEmbedding", TString filename4 = "/data/Work/jets/JetMass/DetectorCorrections/Davide2DRespLHC11a1ppPass2Weighted/h2_pT_DetPart.root", TString filename5 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Filip/filip/JetRespMatrixPP5_LHC13b4_plus_train756_AKT04_Zero_NoOutliers_2c.root", TString file4DWeighted = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/analysis/ResponseWPrepareInputForEmbedding.root", TString file2DWeighted = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/analysis/Output2DResponsePrepareInputForEmbedding.root", TString file4DWeighted2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", 
TString file2DWeighted2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/Output2DResponseJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", 
TString file576 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString list576 = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged"){
   
   //Float_t minptDet = 10, maxptDet = 120;
   Float_t minptDet = 0, maxptDet = 10000;
   Int_t arrayPartDetPartDet[4] = {0, 1, 2, 3}; // this is usually good for Preparation of for Final Response!!!
   Int_t arrayDetPartDetPart[4] = {1, 0, 3, 2}; //this is usually good for JetMassResponseDet
   
   //Train with response from task TreePreparation (pTdet > 2 GeV/c)
   TList *lPrep = ReadFile(filename1, list1);
   if(!lPrep) return;
   THnSparseF *hPrep = (THnSparseF*)lPrep->FindObject("hResponse");
   hPrep->GetAxis(2)->SetRange(hPrep->GetAxis(2)->FindBin(minptDet), hPrep->GetAxis(2)->FindBin(maxptDet));
   Bool_t drawprep = kFALSE;
   
   // response from already weighted file
   TList *l576 = ReadFile(file576, list576);
   if(!l576) return;
   THnSparseF *hsp576 = (THnSparseF*)l576->FindObject("fhnMassResponse");
   if(!hsp576) return;
   
   
   //Response from LHC11a1 (pp pass2, from Davide)
   TFile *file4= new TFile(filename4);
   if(!file4) return;
   TH2F* hptpardetD = (TH2F*)file4->Get("hpTDetPart");
   Int_t rangePjD[2] = {hptpardetD->GetXaxis()->FindBin(minptDet), hptpardetD->GetXaxis()->FindBin(maxptDet)};
   TH1D* hpTpD = hptpardetD->ProjectionY("hpTpD", rangePjD[0], rangePjD[1]);
   hpTpD->SetLineColor(colors[3]);
   hpTpD->SetLineWidth(2);
   TH1D* hpTdD = hptpardetD->ProjectionX("hpTdD");
   hpTdD->GetXaxis()->SetRangeUser(minptDet, maxptDet);
   hpTdD->SetLineColor(colors[3]);
   hpTdD->SetLineWidth(2);
   
   //Response from LHC13b4_plus from Filip
   TFile *file5= new TFile(filename5);
   if(!file5) return;
   TH2D* hptpardetF = (TH2D*)file5->Get("hResponseMatrix");
   Int_t rangePjF[2] = {hptpardetF->GetXaxis()->FindBin(minptDet), hptpardetF->GetXaxis()->FindBin(maxptDet)};
   TH1D* hpTpF = hptpardetF->ProjectionY("hpTpF", rangePjF[0], rangePjF[1]);
   hpTpF->SetLineColor(colors[4]);
   TH1D* hpTdF = hptpardetF->ProjectionX("hpTdF");
   hpTdF->GetXaxis()->SetRangeUser(minptDet, maxptDet);
   hpTpF->SetLineWidth(2);
   hpTdF->SetLineColor(colors[4]);
   hpTdF->SetLineWidth(2);
   
   //draw projections from 4D (the histogram is already weighted, not done with WeightClass)
   WeightClass *docomparisononly = new WeightClass("weight", "Weight", 10);
   docomparisononly->SetOrderAxesRespose(arrayPartDetPartDet[0], arrayPartDetPartDet[1], arrayPartDetPartDet[2], arrayPartDetPartDet[3]);
   
   TH1D** hPrepPj = docomparisononly->All1DProjections(hPrep, 0);
   
   //for(Int_t i = 0; i< 4; i++) {
   //   hPrepPj[i]->SetMarkerStyle(25);
   //   hPrepPj[i]->SetMarkerColor(hPrepPj[i]->GetLineColor());
   //}
   
   //draw projections from 4D (the histogram is already weighted, not done with WeightClass)
   TH1D** h576Pj = docomparisononly->All1DProjections(hsp576, 8);
   
   // weighted 4D response
   TFile *file4DW= new TFile(file4DWeighted);
   if(!file4DW) return;
   THnSparseF *h4DRespW = (THnSparseF*)file4DW->Get("fhResponseFinal");
   h4DRespW->GetAxis(2)->SetRange(h4DRespW->GetAxis(2)->FindBin(minptDet), h4DRespW->GetAxis(2)->FindBin(maxptDet));
   TH1D** h4DRespWPj = docomparisononly->All1DProjections(h4DRespW, 7, "hResponseFinal4D");
   //for(Int_t i = 0; i< 4; i++) {
   //   h4DRespWPj[i]->SetMarkerStyle(25);
   //   h4DRespWPj[i]->SetMarkerColor(h4DRespWPj[i]->GetLineColor());
   //}
   
   // weighted 2D response
   TString name = "fh2DResponseFinal";
   TFile *file2DW= new TFile(file2DWeighted);
   if(!file2DW) return;
   file2DW->ls();
   TH2D* hptpardetW = (TH2D*)file2DW->Get(Form("%s", name.Data()));
   TH1D* hpTpW = (TH1D*)file2DW->Get(Form("%s_pj1Y", name.Data()));
   TH1D* hpTdW = (TH1D*)file2DW->Get(Form("%s_pj0X", name.Data()));
   hpTpW->SetLineStyle(2);
   hpTpW->SetName(Form("%s1087_Par", name.Data()));
   hpTdW->SetLineStyle(2);
   hpTdW->SetLineColor(hpTpW->GetLineColor());
   hpTdW->SetName(Form("%s1087_Det", name.Data()));
   Int_t rangePjW[2] = {hptpardetW->GetXaxis()->FindBin(minptDet), hptpardetW->GetXaxis()->FindBin(maxptDet)};
   //TH1D* hpTpWpj = hptpardetW->ProjectionY("hpTpW", rangePjW[0], rangePjW[1]);
   TH1D* hpTpWpj = hptpardetW->ProjectionY("hpTpW", 0, -1, "E");
   //TH1D* hpTpWpj = (TH1D*)docomparisononly->GetProj(hptpardetW, 0, Form("%sGetProj_0", name.Data()));
   hpTpWpj->SetLineColor(hpTpW->GetLineColor());
   //hpTpWpj->SetMarkerStyle(24);
   //hpTpWpj->SetMarkerColor(hpTpWpj->GetLineColor());
   TH1D* hpTdWpj = hptpardetW->ProjectionX("hpTdW", 0, -1, "E");
   //TH1D* hpTdWpj = (TH1D*)docomparisononly->GetProj(hptpardetW, 1, Form("%sGetProj_1", name.Data()));
   hpTdWpj->GetXaxis()->SetRangeUser(minptDet, maxptDet);
   hpTpWpj->SetLineWidth(2);
   hpTdWpj->SetLineColor(hpTpW->GetLineColor());
   hpTdWpj->SetLineWidth(2);
   //hpTdWpj->SetMarkerStyle(24);
   //hpTdWpj->SetMarkerColor(hpTdWpj->GetLineColor());
   
   // weighted 4D response another train (default 576)
   TFile *file4DW2= new TFile(file4DWeighted2);
   if(!file4DW2) return;
   THnSparseF *h4DRespW2 = (THnSparseF*)file4DW2->Get("fhResponseFinal");
   h4DRespW2->GetAxis(2)->SetRange(h4DRespW2->GetAxis(2)->FindBin(minptDet), h4DRespW2->GetAxis(2)->FindBin(maxptDet));
   TH1D** h4DRespW2Pj = docomparisononly->All1DProjections(h4DRespW2, 9, "hResponseFinal4D2");
   
   // weighted 2D response another train (default 576)
   //TString name = "fh2DResponseFinal";
   TFile *file2DW2= new TFile(file2DWeighted2);
   if(!file2DW) return;
   TH2D* hptpardetW2 = (TH2D*)file2DW2->Get(Form("%s", name.Data()));
   TH1D* hpTpW2 = (TH1D*)file2DW2->Get(Form("%s_pj1Y", name.Data()));
   TH1D* hpTdW2 = (TH1D*)file2DW2->Get(Form("%s_pj0X", name.Data()));
   hpTpW2->SetLineStyle(2);
   hpTpW2->SetName(Form("%s576_Par", name.Data()));
   hpTpW2->SetLineColor(hpTdW2->GetLineColor());
   hpTdW2->SetLineStyle(2);
   hpTdW2->SetName(Form("%s576_Det", name.Data()));
   rangePjW[0] = hptpardetW->GetXaxis()->FindBin(minptDet); 
   rangePjW[1] = hptpardetW->GetXaxis()->FindBin(maxptDet);
   //TH1D* hpTpW2pj = hptpardetW2->ProjectionY("hpTpW2", rangePjW[0], rangePjW[1]);
   TH1D* hpTpW2pj = hptpardetW2->ProjectionY("hpTpW2", 0, -1, "E");
   //TH1D* hpTpW2pj = (TH1D*)docomparisononly->GetProj(hptpardetW2, 0, Form("%s576GetProj_0", name.Data()));
   hpTpW2pj->SetLineColor(hpTdW2->GetLineColor());
   //hpTpWpj->SetMarkerStyle(24);
   //hpTpWpj->SetMarkerColor(hpTpWpj->GetLineColor());
   TH1D* hpTdW2pj = hptpardetW2->ProjectionX("hpTdW2", 0, -1, "E");
   //TH1D* hpTdW2pj = (TH1D*)docomparisononly->GetProj(hptpardetW2, 1, Form("%s576GetProj_1", name.Data()));
   hpTdW2pj->GetXaxis()->SetRangeUser(minptDet, maxptDet);
   hpTpW2pj->SetLineWidth(2);
   hpTdW2pj->SetLineColor(hpTdW2->GetLineColor());
   hpTdW2pj->SetLineWidth(2);

   Bool_t readyproj = kTRUE;
   TLegend *leg = new TLegend(0.4, 0.5, 0.8, 0.9);
   leg->SetFillStyle(0);
   leg->AddEntry(hpTpD, "Davide LHC11a1", "l");
   leg->AddEntry(hpTpF, "Filip LHC13b4", "l");
   if(drawprep) leg->AddEntry(hPrepPj[0], "Preparation LHC13b4", "l");
   if(readyproj){
      leg->AddEntry(hpTpW, "Weighted 2D LHC13b4 - Train1087", "l");
      leg->AddEntry(hpTpW2, "Weighted 2D LHC13b4 - Train576", "l");
   } else{
      leg->AddEntry(hpTpWpj, "Weighted 2D LHC13b4", "l");
      leg->AddEntry(hpTpW2pj, "Weighted 2D LHC13b4 Train576", "l");
   }
   leg->AddEntry(h4DRespWPj[0], "Weighted 4D LHC13b4 - 1087", "l");
   leg->AddEntry(h4DRespW2Pj[0], "Weighted 4D LHC13b4 - 576", "l");
   leg->AddEntry(h576Pj[0], "Original LHC13b4 - 576", "l");
   TPaveText *pv = new TPaveText(0.1, 0.2, 0.5, 0.3, "NDC");
   pv->SetFillStyle(0);
   pv->SetBorderSize(0);
   
   pv->AddText(Form("%.1f < p_{T, det} < %.0f", minptDet, maxptDet));
   
   TCanvas *cPrepD = new TCanvas("cPrepD", "Comparison Prep vs Davide", 600, 800);
   cPrepD->Divide(1,2);
   TCanvas *cPrepF = new TCanvas("cPrepF", "Comparison Prep vs Filip", 600, 800);
   cPrepF->Divide(1,2);
   TCanvas *cFD = new TCanvas("cFD", "Comparison Filip vs Davide", 600, 800);
   cFD->Divide(1,2);
   
   TCanvas *cCmpProj2D = new TCanvas("cCmpProj2D", "Comparison projections 2D", 600, 800);
   cCmpProj2D->Divide(1,2);
   
   //comparison with Davide
   cPrepD->cd(1);
   gPad->SetLogy();
   
   hpTpD->Draw("");
   if(readyproj){
      hpTpW->Draw("sames");
      hpTpW2->Draw("sames");
   } else {
      hpTpWpj->Draw("sames");
      hpTpW2pj->Draw("sames");
   }
   h4DRespWPj[3]->Draw("sames");
   h4DRespW2Pj[3]->Draw("sames");
   h576Pj[3]->Draw("sames");
   if(drawprep) hPrepPj[3]->Draw("sames");
   pv->Draw();
   cPrepD->cd(2);
   gPad->SetLogy();
   
   hpTdD->Draw("");
   if(readyproj){
      hpTdW->Draw("sames");
      hpTdW2->Draw("sames");
   } else {
      hpTdWpj->Draw("sames");
      hpTdW2pj->Draw("sames");
   }
   h4DRespWPj[2]->Draw("sames");
   h4DRespW2Pj[2]->Draw("sames");
   h576Pj[2]->Draw("sames");
   if(drawprep) hPrepPj[2]->Draw("sames");
   leg->Draw();
   
   
   //Comparison with Filip
   cPrepF->cd(1);
   gPad->SetLogy();
  
   hpTpF->Draw("");
   hpTpWpj->Draw("sames");
   if(drawprep) hPrepPj[3]->Draw("sames");
   cPrepF->cd(2);
   pv->Draw();
   gPad->SetLogy();
   hpTdF->Draw("");
   hpTdWpj->Draw("sames");
   if(drawprep) hPrepPj[2]->Draw("sames");
   leg->Draw();
   
   //Comparison Davide Filip
   cFD->cd(1);
   gPad->SetLogy();
   hpTpF->Draw();
   hpTpD->Draw("sames");
   pv->Draw();
   cFD->cd(2);
   gPad->SetLogy();
   hpTdF->Draw();
   hpTdD->Draw("sames");
   leg->Draw();
   
   // comparison of the projection saved to File and those redone here (temporary for cross-check)
   cCmpProj2D->cd(1);
   //particle level
   hpTpW->Draw("");
   hpTpWpj->Draw(""); //Y
   //hpTpW2->Draw("sames");
   //hpTpW2pj->Draw("sames");
   
   cCmpProj2D->cd(2);
   //detector level
   hpTdW->Draw("");
   hpTdWpj->Draw("sames"); //X
   //hpTdW2->Draw("sames");
   //hpTdW2pj->Draw("sames");
   
   SaveCv(cPrepD);
   SaveCv(cPrepF);
   SaveCv(cFD);
   
   // do ratios
   
   // 4D over 2D of the same train 1087 and 576

   TCanvas *cR4over2ptP = new TCanvas("cR4over2ptP", "Ratio 4D/2D weighted ptP distributions");
   TCanvas *cR4over2ptD = new TCanvas("cR4over2ptD", "Ratio 4D/2D weighted ptD distributions");
   
   //particle level
   TH1 *hR4over2ptP576  = 0x0;
   TH1 *hR4over2ptP1087 = 0x0;
   
   TH1 *hpTpW2copy = 0x0;
   TH1 *hpTpWcopy  = 0x0;
   
   UniformTH1FForDivide(h4DRespW2Pj[3], hpTpW2, hR4over2ptP576, hpTpW2copy);
   
   UniformTH1FForDivide(h4DRespWPj[3], hpTpW, hR4over2ptP1087, hpTpWcopy);
   
   hR4over2ptP576 ->SetName("hR4over2ptP576"); 
   hR4over2ptP1087->SetName("hR4over2ptP1087");
   hR4over2ptP576 ->Divide(hpTpW2copy); 
   hR4over2ptP1087->Divide(hpTpWcopy);
   
   //detector level
   TH1 *hR4over2ptD576  = 0x0;
   TH1 *hR4over2ptD1087 = 0x0;
   
   TH1 *hpTdW2copy = 0x0;
   TH1 *hpTdWcopy  = 0x0;
   
   UniformTH1FForDivide(h4DRespW2Pj[2], hpTdW2, hR4over2ptD576, hpTdW2copy);
   
   UniformTH1FForDivide(h4DRespWPj[2], hpTdW, hR4over2ptD1087, hpTdWcopy);

   hR4over2ptD576 ->SetName("hR4over2ptD576"); 
   hR4over2ptD1087->SetName("hR4over2ptD1087");
   hR4over2ptD576 ->Divide(hpTdW2copy);
   hR4over2ptD1087->Divide(hpTdWcopy);

   cR4over2ptP->cd();
   hR4over2ptP576 ->Draw("");
   hR4over2ptP1087->Draw("sames");
   
   cR4over2ptD->cd();
   hR4over2ptD576 ->Draw("");
   hR4over2ptD1087->Draw("sames");
   
   // compare train 576 and 1087 4D response
   TCanvas *cR1087over576ptP = new TCanvas("cR1087over576ptP", "Ratio 4D 1087/576 weighted ptP distributions");
   TCanvas *cR1087over576ptD = new TCanvas("cR1087over576ptD", "Ratio 4D 1087/576 weighted ptD distributions");

   //particle level
   
   TH1 *hR1087over576ptP  = 0x0;
   TH1* h4DRespW2PjPcopy = 0x0;
   UniformTH1FForDivide(h4DRespWPj[3], h4DRespW2Pj[3], hR1087over576ptP, h4DRespW2PjPcopy);
   hR1087over576ptP->SetName("hR1087over576ptP");
   hR1087over576ptP->Divide(h4DRespW2PjPcopy);
   
   //detector level
   
   TH1 *hR1087over576ptD  = 0x0;
   TH1* h4DRespW2PjDcopy = 0x0;
   UniformTH1FForDivide(h4DRespWPj[2], h4DRespW2Pj[2], hR1087over576ptD, h4DRespW2PjDcopy);
   hR1087over576ptD->SetName("hR1087over576ptD");
   hR1087over576ptD->Divide(h4DRespW2PjDcopy);
   
   cR1087over576ptP->cd();
   hR1087over576ptP->Draw();
   
   cR1087over576ptD->cd();
   hR1087over576ptD->Draw();
}

//________________________________________________________________________________________________
// Define input files for comparison of N generic responses weighted
void ComparisonFinalResponseFixedInput(){
   
   //number of files to be comared
   const Int_t ninputs = 2;
   
   //file paths
   TString files[ninputs] =
   
   // local test output
   //{"ResponseWJetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme.root", "ResponseWPrepareInputForEmbedding.root"};
   
   // 576 before set axis order, after set axis order
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/MPMDpPpD/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root"};

   // 1087 before set axis order, after set axis order
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/analysis/ResponseWPrepareInputForEmbedding.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/analysis/MPMDpPpD/ResponseWPrepareInputForEmbedding.root"};
   
   // 576 weighted and 1087 weighted correct axis order
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/MPMDpPpD/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/analysis/MPMDpPpD/ResponseWPrepareInputForEmbedding.root"};
   
   // 576 weighted and 1096 weighted correct axis order
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/MPMDpPpD/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1096/outputpTHardBins/mergeRuns/analysis/MPMDpPpD/ResponseWJetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme.root"};

   // 1119 contains both tasks
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/analysis/ResponseWJetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/analysis/ResponseWPrepareInputForEmbedding.root"};
   
   // 576 weighted and embedding (different trials)
   {"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/BinByBinCorrDetPtParPt/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160222/Input1134/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root"};
   // 1134
   files[1] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root";
   //files[1] = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160226/outputAfterRuedBugFix/3Runs/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root";
   // before and after Ruediger bug fix
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160226/output195351BeforeBugFix/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160226/outputAfterRuedBugFix/195351/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root"};
   
   //compare embedding in different runs
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160226/outputAfterRuedBugFix/195389/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160226/outputAfterRuedBugFix/195351/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root"};

   // legends
   TString legs[ninputs] = 
   // local test output
   //{"ResponseTask", "PrepareTask"};
   // 576 weighted and 1087 weighted correct axis order
   //{"576", "1087"};
   
   // 576 before set axis order, after set axis order
   //{"576 Before bug fix", "576 After bug fix"};
   
   // 1087 before set axis order, after set axis order
   //{"1087 Before bug fix", "1087 After bug fix"};
   
   // 576 weighted and 1096 weighted correct axis order
   //{"576", "1096"};
   {"576", "1134"};
   
   // 576 weighted and embedding
   //{"576", "Embedding"};
   //legs[0] = "1134";
   
   //compare embedding in different runs
   //{"195389", "195351"};
   // before and after Ruediger bug fix
   //{"Before", "After"};
   Int_t axSel = 2; // pt particle level
   axSel = 3; // pt detector level
   Double_t ptrange[2] = {40, 120};
   //start the comparison
   ComparisonFinalResponse(ninputs, files, legs, axSel, ptrange);
   
}



//________________________________________________________________________________________________
// perform weighting in 2D instead of 4 using class WeightClass

void Weight2D(TString pathPythiaFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/", TString listnameInp = "PrepareInputForEmbedding", TString filename = "", Int_t parax = 1, Int_t detax = 0, TString hname = "fhnMassResponse"){
   /// Two uses : 
   /// 1) read the thnsparses filled with the tree entries in PrepareResponseFromTree
   ///    Input:
   ///          - pathPythiaFile: path pointing to the base dir for the bin-by-bin output of the task that contains the xsec and ntrials histograms
   ///          - listnameInp :  name of the list where to read xsec and ntrials
   ///          - filename: must be ""
   /// 2) read the THnsparses with the response from a TList 
   ///    Input:
   ///          - pathPythiaFile: path pointing to the base dir for the bin-by-bin output of the task that contains the xsec and ntrials histograms
   ///          - listnameInp :  name of the list where to read xsec and ntrials and the input response
   ///          - filename: name of the file where the input response is
   ///          - hname : name of the THnSparse. Either "fhnMassResponse" or "hResponse"

   Int_t firstBin = 0;
   Int_t npthb =10;//2 ;//
   
   TH2 **fhResponse = new TH2*[npthb];
   TFile *fin = 0x0;
   
   if(filename.IsNull()){
      
      fin = new TFile("ResponseFromTree.root");
      if(!fin->IsOpen()){
      	 Printf("Input file not found");
      	 return;
      }
   }

   WeightClass *doweight = new WeightClass("weight", "Weight", npthb);
   doweight->SetIs2D(kTRUE, parax, detax);
   doweight->Initialise(pathPythiaFile, listnameInp, filename, listnameInp);
   
   TCanvas *cProjResp = new TCanvas(Form("cProjResp"), Form("Response before weighting "), 600, 800);
   cProjResp->Divide(1,2);
   TCanvas *cProjRespW = new TCanvas(Form("cProjRespW"), Form("Response after weighting "), 600, 800);
   cProjRespW->Divide(1,2);
   
   TLegend *legpTh = new TLegend(0.1, 0.4, 0.4, 0.7);
   legpTh->SetBorderSize(0);
   legpTh->SetFillStyle(0);
   
   //loop on pT hard bins
   for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
      //     for(Int_t ipthb = firstBin; ipthb<2; ipthb++){
      //TString listnameInp = Form("JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH%d_pT0150_E_scheme_TC", ipthb+1);
      
      if(fin){ //case 1
      	 fhResponse[ipthb] = (TH2*)fin->Get(Form("h2Response%d", ipthb+1));
      	 doweight->ReadInputData2D(ipthb, fhResponse[ipthb]);
      } else { //case 2
      	 Bool_t set = doweight->SetResponse(ipthb, hname);
      	 if(!set) return;
      }
      
      doweight->ReadWeight(ipthb);

      cProjResp->cd(1);
      gPad->SetLogy();
      TH1D* hpTpar = (TH1D*)doweight->GetResponse2DProjPar(ipthb);
      Printf(" entries %f", hpTpar->GetEntries());
      if(hpTpar) {
      	 if(ipthb == firstBin) hpTpar->Draw();
      	 else hpTpar->Draw("sames");
      }
      cProjResp->cd(2);
      gPad->SetLogy();
      TH1* hpTdet = doweight->GetResponse2DProjDet(ipthb);
      if(hpTdet) {
      	 if(ipthb == firstBin) hpTdet->Draw();
      	 else hpTdet->Draw("sames");
      }
      Bool_t success = doweight->WeightResponse(ipthb);
      if(!success) Printf("Weight response returned %d", success);
      
      cProjRespW->cd(1);
      gPad->SetLogy();
      TH1D* hpTparW = (TH1D*)doweight->GetResponse2DWProjPar(ipthb);
      Printf(" entries %f", hpTparW->GetEntries());
      if(hpTparW) {
      	 if(ipthb == firstBin) hpTparW->Draw();
      	 else hpTparW->Draw("sames"); 
      }
      cProjRespW->cd(2);
      gPad->SetLogy();
      TH1* hpTdetW = doweight->GetResponse2DWProjDet(ipthb);
      if(hpTdetW) {
      	 if(ipthb == firstBin) hpTdetW->Draw();
      	 else hpTdetW->Draw("sames");
      }
   }
   
   if(!doweight->Result()) {
      Printf("Something went wrong");
      return;
   }

   TH2* hFinalResp = doweight->GetResponse2DWFinal();
   TH1* hFinalPar = doweight->GetResponse2DWFinalProjPar();
   TH1* hFinalDet = doweight->GetResponse2DWFinalProjDet();
   
   TCanvas *cProjRespFinal = new TCanvas("cProjRespFinal", "Final response", 600, 800);
   cProjRespFinal->Divide(1,2);
   cProjRespFinal->cd(1);
   gPad->SetLogy();
   hFinalPar->Draw();
   cProjRespFinal->cd(2);
   gPad->SetLogy();
   hFinalDet->Draw();
   
   TString suff = ""; 
   if(!filename.IsNull()) suff = listnameInp;
   TFile *fout = new TFile(TString::Format("Output2DResponse%s.root",  suff.Data()), "recreate");
   fout->cd();
   hFinalResp->Write();
   hFinalPar ->Write();
   hFinalDet ->Write();
   
}

//________________________________________________________________________________________________
// comparison of distributions before weighting

void CompareBeforeWeight(){
   
   const Int_t nSeries = 2;
   TString paths[nSeries] = //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1096/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/StudyDiscrepancyJetMassResponseDetVsPrepareInputForEmbedding/", "/data/Work/jets/JetMass/DetectorCorrections/StudyDiscrepancyJetMassResponseDetVsPrepareInputForEmbedding/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1127/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1127/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1119/outputpTHardBins/mergeRuns/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1123/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1132/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1133/outputpTHardBins/mergeRuns/"};
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1123/outputpTHardBins/mergeRuns/"};
   //paths[1] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1087/outputpTHardBins/mergeRuns/";
   //paths[0] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1131/outputpTHardBins/mergeRuns/";
   //paths[0] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1132/outputpTHardBins/mergeRuns/";
   //paths[0] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1133/outputpTHardBins/mergeRuns/";
   {"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160222/Input1134/analysis/"};
   
   TString listnameInp[nSeries] = 
   //{"PrepareInputForEmbedding", "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme"};
   //listnameInp[0] = "JetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme";
   //listnameInp[0] = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme";
   //listnameInp[1] = "PrepareInputForEmbedding";
   //{"PrepareInputForEmbedding", "JetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme"};
   //listnameInp[1] = "JetMassResponseDet_Jet_AKTChargedR040_tracks_pT0150_E_scheme";
   {"JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme", "JetShapeConst_JetRhosub_AKTChargedR040_PicoTracks"};
   
   TString respname[nSeries] = {"hResponse", "fhnMassResponse"};
   respname[0] = "fhnMassResponse";
   //respname[1] = "hResponse";
   respname[1] = "fhnMassResponse_0";
   
   TString legF[nSeries] = 
   {"Train1087", "Train576"};
   //{"TaskPrepare", "TaskResponse"};
   //legF[0] = "Train1096";
   //legF[0] = "Train1119";
   //legF[0] = "Train1123";
   //legF[0] = "Train1127";
   //legF[0] = "Train1131";
   //legF[0] = "Train1132";
   //legF[0] = "Train1133";
   legF[0] = "Train1134";
   //legF[1] = "Train1087";
   //legF[1] = "Train1119";
   //legF[1] = "Train1123";
   //legF[1] = "Train1134";
   legF[1] = "Embedding";
   TString filename = "AnalysisResults.root";
   
   Int_t firstBin = 0;
   const Int_t npthb =10; //2 ;//
   
   THnSparseF **fhResponse = new THnSparseF*[npthb];
   TH1D* hRatio[4][npthb][nSeries];
   TH1D* hRespPj[4][npthb][nSeries];
   //canvas definition
   
   TCanvas *cProjResp = new TCanvas(Form("cProjResp%s%s", legF[0].Data(),  legF[1].Data()), Form("Response before weighting %s vs %s", legF[0].Data(),  legF[1].Data()), 1000, 1000);
   cProjResp->Divide(2,2);
   TCanvas *cProjRespRatio = new TCanvas(Form("cProjRespRatio%s%s", legF[0].Data(),  legF[1].Data()), Form("Response Ratio  %s vs %s", legF[0].Data(),  legF[1].Data()), 1000, 1000);
   cProjRespRatio->Divide(2,2);

   TLegend *legpTh = new TLegend(0.6, 0.6, 0.9, 0.9);
   legpTh->SetBorderSize(0);
   legpTh->SetFillStyle(0);
   
   TLegend *legFiles = new TLegend(0.3, 0.8, 0.6, 0.9);
   legFiles->SetBorderSize(0);
   legFiles->SetFillStyle(0);

   Int_t MpMdptpptd[4] = {1, 0, 3, 2} ; //"fhnMassResponse"
      
   for(Int_t ifile = 0; ifile < nSeries; ifile++){
   //for(Int_t ifile = 0; ifile < 1; ifile++){
   
      WeightClass *doweight = new WeightClass(Form("weight%d", ifile), Form("Weight %d", ifile), npthb);
      if(respname[ifile] == "hResponse") {
      	 MpMdptpptd[0] = 0; MpMdptpptd[1] = 1; MpMdptpptd[2] = 2; MpMdptpptd[3] = 3;
      }
      if(respname[ifile] == "fhnMassResponse") {
      	 MpMdptpptd[0] = 1; MpMdptpptd[1] = 0; MpMdptpptd[2] = 3; MpMdptpptd[3] = 2;
      }
      
      doweight->SetOrderAxesRespose(MpMdptpptd[0], MpMdptpptd[1], MpMdptpptd[2], MpMdptpptd[3]);
      if(legF[ifile] != "Embedding") doweight->Initialise(paths[ifile], listnameInp[ifile], paths[ifile], listnameInp[ifile]);
      else doweight->Initialise(paths[0], listnameInp[0], paths[1], listnameInp[1]);
      //loop on pT hard bins
      for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
      	 //for(Int_t ipthb = firstBin; ipthb<1; ipthb++){
      	 //doweight->SetCut(MpMdptpptd[3], 5.);
      	 Bool_t set = doweight->SetResponse(ipthb, respname[ifile]);
      	 if(!set) {
      	    Printf("Error File %d, pThard bin %d", ifile, ipthb);
      	    return;
      	 }
      	 if(ipthb == 0) {
      	    THnSparse *htempsSp = doweight->GetResponsePtHbin(0);
      	    Printf("++++++++++++++++ %s order %d = %s, %d = %s, %d = %s, %d = %s", respname[ifile].Data(), MpMdptpptd[0], htempsSp->GetAxis(MpMdptpptd[0])->GetTitle(),  MpMdptpptd[1], htempsSp->GetAxis(MpMdptpptd[1])->GetTitle(), MpMdptpptd[2], htempsSp->GetAxis(MpMdptpptd[2])->GetTitle(), MpMdptpptd[3], htempsSp->GetAxis(MpMdptpptd[3])->GetTitle());
      	  
      	 }
      	 
      	 Printf("Response name %s", doweight->GetResponsePtHbin(ipthb)->GetName());
      	 TH1D **hRespPjTmp = doweight->All1DResponseProj(ipthb, ipthb, Form("hPrjTmpFile%dPtHb%d", ifile, ipthb));
      	    
      	 for(Int_t ipj = 0; ipj<4; ipj++){
      	    
      	    hRespPj[ipj][ipthb][ifile] = (TH1D*)hRespPjTmp[ipj]->Clone(Form("hPrjFile%dPtHb%d_%d", ifile, ipthb, ipj));
      	    //hRespPj[ipj][ipthb][ifile]->SetName(Form("hPrjFile%dPtHb%d_%d", ifile, ipthb, ipj));
      	    
      	    Printf("%s bin %d -> number of bins %d", legF[ifile].Data(), ipthb, hRespPj[ipj][ipthb][ifile]->GetNbinsX());
      	    
      	    hRespPj[ipj][ipthb][ifile]->SetMarkerStyle(markers[ifile]);
      	    hRatio[ipj][ipthb][ifile] = (TH1D*)hRespPj[ipj][ipthb][ifile]->Clone(Form("hRatio%d%d_pTHb%d", ifile, ipj, ipthb));
      	      
       	    Printf("File %d, Pt hard bin %d, Distr %d, Content Bin 10 %f", ifile, ipthb, ipj, hRespPj[ipj][ipthb][ifile]->GetBinContent(10));
       	    
       	    //uncomment this to draw with the original binning (and comment in the loop below)
      	    //cProjResp->cd(ipj+1);
      	    //gPad->SetLogy();
      	    //
      	    //if((ifile == 0) && (ipthb == firstBin)) hRespPj[ipj][ipthb][ifile]->Draw();
      	    //else hRespPj[ipj][ipthb][ifile]->Draw("sames");
      	    //Printf("Names: %s, %s, %s", hRespPj[ipj][ipthb][ifile]->GetName(), hRatio[ipj][ipthb][ifile]->GetName(), hRespPjTmp[ipj]->GetName());
      	   
      	 }
      	 
      	 if(ipthb == firstBin) legFiles->AddEntry(hRespPj[0][ipthb][ifile], legF[ifile], "P");
      	 if((ifile == 0)) legpTh->AddEntry(hRespPj[0][ipthb][ifile], Form("pThard %d", ipthb+1), "L");
      	 
      	 //apply weights
      	 
      	 
      } //loop on pThard
   }//loop on files
   
   
   
    //ratios
    //loop on pT hard bins
    TH1D* hRatiopjRebin[4][npthb];
    for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
      //for(Int_t ipthb = firstBin; ipthb<1; ipthb++){
      	 //loop on variables
       for(Int_t ipj = 0; ipj<4; ipj++){
       	  TH1* hRespPjDiv1 = 0x0;
       	  TH1* hRespPjDiv2 = 0x0;
       	  //Printf("Before %d %d", hRatio[ipj][ipthb][0]->GetNbinsX(), hRatio[ipj][ipthb][1]->GetNbinsX());
       	  Printf("Ratios -- Pt hard bin %d, Distr %d, Content Bin 10: File1 = %.1f File1 = %.1f , Ratio = %f", ipthb, ipj, hRespPj[ipj][ipthb][0]->GetBinContent(10), hRespPj[ipj][ipthb][1]->GetBinContent(10), hRespPj[ipj][ipthb][0]->GetBinContent(10)/hRespPj[ipj][ipthb][1]->GetBinContent(10) );
       	  
       	  UniformTH1FForDivide(hRatio[ipj][ipthb][0], hRatio[ipj][ipthb][1], hRespPjDiv1, hRespPjDiv2, "TH1D");
       	  if(!hRespPjDiv1 || !hRespPjDiv2) {
       	     Printf("Error! %p %p", hRespPjDiv1, hRespPjDiv2); return;
       	  }
       	  Printf("%s %s", hRespPjDiv1->GetName(), hRespPjDiv2->GetName());
       	  Printf("%d %d", hRespPjDiv1->GetNbinsX(), hRespPjDiv2->GetNbinsX());
       	  
       	  hRespPjDiv1->SetName(Form("hRatioRebin_pthb%d_pj%d", ipthb, ipj));
       	  hRespPjDiv1->SetTitle(Form("hRatioRebin_pthb%d_pj%d;%s;%s/%s", ipthb, ipj, hRespPjDiv1->GetXaxis()->GetTitle(), legF[0].Data(), legF[1].Data()));
       	  //uncomment this to draw with the common binning used for dividing (and comment in the above loop)
       	  cProjResp->cd(ipj+1);
       	  gPad->SetLogy();
       	  if(ipthb == firstBin) {
       	     hRespPjDiv1->DrawClone();
       	     hRespPjDiv2->DrawClone("sames");
       	  }
       	  else {
       	     hRespPjDiv1->DrawClone("sames");
       	     hRespPjDiv2->DrawClone("sames");
       	  }
       	  hRespPjDiv1->Divide(hRespPjDiv2);
      	  Printf("Division successful -> Bin 10 gives %f", hRespPjDiv1->GetBinContent(10));
      	  cProjRespRatio->cd(ipj+1);
      	  //gPad->SetLogy();
       	  hRespPjDiv1->SetLineStyle(2);
       	  if(ipthb == firstBin) hRespPjDiv1->Draw();
       	  else hRespPjDiv1->Draw("sames");
       	  
       } //end loop on variables
    } //end loop on pT hard bins
    
    cProjResp->cd(1);
    legpTh->Draw();
    cProjResp->cd(3);
    legFiles->Draw();
    SaveCv(cProjResp);
    cProjRespRatio->cd(1);
    legpTh->Draw();
    SaveCv(cProjRespRatio);
}

//________________________________________________________________________________________________

void RunStudyDistributionShapeExcludingpTHbins(){
   
   
   /// Conclusion: I can use the PtHard bins >2 or 3, it's good enough 
   const Int_t n = 5;
   Int_t iPtHb[n] = {2, 3, 4, 5, 6};
   Double_t pTcut = 10; //20;// 30; // 80; //
   TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/";
   TString inputList = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme";
   TString hname = "fhnMassResponse";
   
   StudyDistributionShapeExcludingpTHbins(n, iPtHb, pTcut, inputFile, inputList, hname);

}

//_________________________________________________________________________________________________

void PerformRhoWeighting(TString pathResponseFinal = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/forRhoW/analysis/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", TString pathDataRho = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1200-1202/AnalysisResults.root", TString listDataRho = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv", TString histoDataRho = "fhnRhoVsRhoMVsLeadJetPtVsMassVsCent"){
	
	// never used, never finished probably
	
	TFile *finResp = new TFile(pathResponseFinal);
	if(!finResp){
		Printf("File containing reponse not found in %s", pathResponseFinal.Data());
		return;
	}
	
	WeightClass *doweight = new WeightClass("weight", "Weight", 1); // the pthbin is not relevant here
	
	doweight->SetDoWeightForBkg(kTRUE);
	doweight->SetDataPathForBkd(pathDataRho);
	doweight->SetDataListForBkdName(listDataRho);
	doweight->SetDataHistForBkdName(histoDataRho);
	
	THnSparseF*  fResponsefinal = (THnSparseF*) finResp->Get("fnhResponseFinal");
	
	doweight->SetWeightedResponseFinal(fResponsefinal);
	
	
	Bool_t res = doweight->ReadRhoFromData();
	if(!res){
		
		Printf("No data input available");
		return;
	}
	
	TList *listRhoH = doweight->GetListOfRhoAndRhoMProjections();
	if(!listRhoH){
		Printf("List of rho histograms not found");
		return;
	}
	
	TCanvas *cRho = new TCanvas("crho", "Rho, several bins", 600, 600);
	TCanvas *cRhom = new TCanvas("crhom", "Rho_m, several bins", 600, 600);

	for(Int_t il = 0; il<listRhoH->GetEntries(); il+=2){
		TH1D *hrho = (TH1D*)listRhoH->At(il);
		TH1D *hrhom= (TH1D*)listRhoH->At(il);
		
		cRho->cd();
		hrho->Draw("sames");
		cRhom->cd();
		hrhom->Draw("sames");
	}
	
	
}

void RunWeightForRhoBins(Int_t numberrhobins = 4, Bool_t donorm = kTRUE, TString pathPythiaFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/", Int_t ninputrespfiles = 1, TString pathDataFile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/forRhoW/output/", TString listnameInp1 = "JetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks", Bool_t detFrom4vecEmb = kFALSE, Bool_t fineBin = kFALSE){
	// run the normal pt-hard bin weighting for nrhobins. Each of them have to be weighted a posteriori and summed
	// ninputrespfiles  = if each rho bin is in a different file, ninputrespfiles is the number of such files. Assumed path pathDataFile/00, ... pathDataFile/0(ninputrespfiles-1)

	TStopwatch watch;
	watch.Start();
   
	
	Int_t nrhobins = 0;
	Double_t drho = 0;
	Double_t *rhobinlims = 0x0;
	if(numberrhobins > 0){
		if(numberrhobins == 5) {
			nrhobins = 5; 
			drho = 1.2;
		}
		if(numberrhobins == 4) {
			nrhobins = 4; 
			drho = 1.5;
		}
		
		rhobinlims = new Double_t[nrhobins+1];
		for(Int_t ibrho = 0; ibrho < nrhobins+1; ibrho++){
			rhobinlims[ibrho] = drho*ibrho;
		}
		//include "all" in the last bin
		rhobinlims[nrhobins] = 10.;
	} else {
		if(numberrhobins == -2) {
			nrhobins = 2; 
			rhobinlims = new Double_t[nrhobins+1];
			rhobinlims[0] = 4.5;
			rhobinlims[1] = 7;
			rhobinlims[2] = 10;
		}
		if(numberrhobins == -5){
			//negative means that bins have different width, setting the array by hand
			nrhobins = 5;
			rhobinlims = new Double_t[nrhobins+1];
			rhobinlims[0] = 0  ;
			rhobinlims[1] = 1.5;
			rhobinlims[2] = 3;
			rhobinlims[3] = 4.5;
			rhobinlims[4] = 7;
			rhobinlims[5] = 10;
			
		}
		
	}
	TString filename = "AnalysisResults.root";
	TString hname = "fhnMassResponse_0";
	Int_t npthb = 10;
	Int_t ptHardFirst = 2;
	Int_t firstBin = 1;
	if(ptHardFirst < firstBin) ptHardFirst = firstBin+1; 
	
	THnSparseF*  hresponsefinal[nrhobins];
	THnSparseF*  hresponsefinalRhoW = 0x0;
	/*
	TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv";
	if(listnameInp1.Contains("Const")) listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub_TCConst";
	TString hdataname = "fhnRhoVsRhoMVsLeadJetPtVsMassVsCent";
	
	TString datarhofile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/output/AnalysisResults.root";
	
	TList *listdata = ReadFile(datarhofile, listname);
	if(!listdata) {
		Printf("Error, no data found");
		return;
	}
	THnSparseF *hspdata = (THnSparseF*)listdata->FindObject(hdataname);
	TH1D *hRhoData = hspdata->Projection(0);
	
	hRhoData->Scale(1./hRhoData->Integral());
	*/
	
	//read the ratio directly from the prepared file (if not present run PrepareRhoWeightingFactors)
	TString datarhofile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/RhoFromDataRatiosNPtDetbins2Deriv.root";
	if(listnameInp1.Contains("Const")) datarhofile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/RhoFromDataRatiosNPtDetbins2Const.root";
	TString hdataname = "hRhoDataR_Bin1";
	
	TFile *finRhoData = new TFile(datarhofile);
	if(!finRhoData->IsOpen()) {
		Printf(" %s not found, check path or run PrepareRhoWeightingFactors first", datarhofile.Data());
		return;
	}
	
	TH1D* hRhoData = (TH1D*)finRhoData->Get(hdataname);
	if(!hRhoData) {
		Printf("Rho factors %s not found in %s", hdataname.Data(), datarhofile.Data());
		finRhoData->ls();
		return;
	}
	
	TCanvas *cProjRespRhoBins = new TCanvas("cProjRespRhoBins", "cProjRespRhoBins", 800, 800);
	cProjRespRhoBins->Divide(2,2);
	
	TCanvas *cRhoVsPtH = new TCanvas("cRhoVsPtH", "QA : Rho vs pT hard bin", 600, 600);
	
	TFile *fout = new TFile(Form("Respose%sRhoW.root", listnameInp1.Data()), "recreate");
	
	TString basepathDataFile = pathDataFile;
	
	for(Int_t ibrho = 0; ibrho < nrhobins; ibrho++){
		//set rho factor
		Int_t binlimsRhoFactor[2] = { hRhoData->FindBin(rhobinlims[ibrho]), hRhoData->FindBin(rhobinlims[ibrho+1]- 0.001)};
		Int_t binRhoFactor = (binlimsRhoFactor[0] + binlimsRhoFactor[1]) * 0.5;  
		Double_t respRhoFact = hRhoData->GetBinContent(binRhoFactor);
		
		Printf("Scaling the rho bin %d by %f", ibrho, respRhoFact);
		
		WeightClass *doweight = new WeightClass(Form("weight%d", ibrho), "Weight", npthb);
		
		doweight->SetRhoSelBin(rhobinlims[ibrho], rhobinlims[ibrho+1]);
		
		//THnSparseF **fhResponse = new THnSparseF*[npthb];
		
		Int_t MpMdptpptd[4] = {1, 0, 3, 2};
		
		//Set strings for initialisation
		TString path = basepathDataFile;
		TString list = listnameInp1;
		if(ninputrespfiles > 1) {
			path+=Form("/0%d/", ibrho);
			list += Form("Rho%d", ibrho);
		}
		Printf("Input response in %s", path.Data());
		
		doweight->Initialise(pathPythiaFile, "PrepareInputForEmbedding", Form("%s", path.Data()), list);
		
		doweight->SetDoNormalisation(donorm);
		doweight->SetFineFinning(fineBin);
		doweight->ReadNumberOfJetPerPtHBin(Form("%s/NJetsPerPtHBin.root", pathPythiaFile.Data()));
		
		
		/*
		// Get the factor for this rho bin
		//assuming same binning in rho for the resp and the data
		Double_t datarhoContent = hRhoData->Integral(hRhoData->FindBin(rhobinlims[ibrho]), hRhoData->FindBin(rhobinlims[ibrho+1]-0.001));
		*/
		//loop on pT hard bins
		for(Int_t ipthb = firstBin; ipthb<npthb; ipthb++){
			doweight->SetOrderAxesRespose(MpMdptpptd[0], MpMdptpptd[1], MpMdptpptd[2], MpMdptpptd[3]);
			
			Bool_t set = doweight->SetResponse(ipthb, hname);
			if(!set) continue;
			else Printf("Response successfully set");
			
			Bool_t outliersremoved = doweight->CleanOutliers(ipthb);
			if(!outliersremoved){
				Printf("Ouliers could not be removed");
				continue;
			} else Printf("Successfully removed outliers");
			
			THnSparseF* hRespOutlRem = doweight->GetResponsePtHbin(ipthb);
			hRespOutlRem->SetName(Form("hRespRho%dPtHb%d", ibrho, ipthb));
			fout->cd();
			hRespOutlRem->Write();
			
			Bool_t weighttaken = doweight->ReadWeight(ipthb);
			if(!weighttaken) {
				Printf("Failed reading weight, continue");
				continue;
			} else Printf("Weight successfully set");
			
			Bool_t success = doweight->WeightResponse(ipthb);
			if(!success) Printf("Weight response returned %d", success);
			else Printf("Weighting procedure was successful");
			
			TH1D**  hprojRespWN = doweight->All1DResponseWNormProj(ipthb, ipthb);
			
		}
		
		if(ptHardFirst > 0) doweight->SetFirstPtHard(ptHardFirst);
		Bool_t results = doweight->Result();
		if(!results) Printf("The sum of the pT hard bins failed");
		else Printf("Final response successfully obtained");
		
		//TH1D *hRhoVsPtH = doweight->GetRhoVsPtH();
		//hRhoVsPtH->SetName(Form("hRhoVsPtH%d", ibrho));
		//hRhoVsPtH->SetLineColor(colors[ibrho]);
		//cRhoVsPtH->cd();
		//if(ibrho == 0) hRhoVsPtH->Draw();
		//else hRhoVsPtH->Draw("sames");
		
		//Double_t respRhoFact = doweight->GetRhoFactor();
		
		hresponsefinal[ibrho] = doweight->GetResponseWFinal();
		if(!hresponsefinal[ibrho]){
			Printf("Error! response is null");
			continue;
		}
		hresponsefinal[ibrho]->SetName(Form("%sRho%d", hresponsefinal[ibrho]->GetName(), ibrho));
		//Double_t scale = 1.;
		//if(respRhoFact > 1.e-6) {
		//	scale =  datarhoContent/respRhoFact;
		//	hresponsefinal[ibrho]->Scale(scale);
		//}
		TH1D** hpj = All1DProjections(hresponsefinal[ibrho], colors[ibrho], Form("hPjAll%d", ibrho));
		
		for(Int_t i = 0; i<4; i++){
			cProjRespRhoBins->cd(i+1);
			if(ibrho == 0) hpj[i]->Draw();
			else  hpj[i]->Draw("sames");
			
		}
		if(!hresponsefinalRhoW) hresponsefinalRhoW = (THnSparseF*) hresponsefinal[ibrho]->Clone("hresponsefinalRhoW");
		else hresponsefinalRhoW->Add(hresponsefinal[ibrho], respRhoFact); 
		fout->cd();
		hresponsefinal[ibrho]->Write();
		Printf("clearing mem bin %d ", ibrho);
		doweight->ClearMem();
	}

	fout->cd();
	if(hresponsefinalRhoW) hresponsefinalRhoW->Write();
	
	watch.Stop();
	watch.Print();

}

void ReAnalyseRespRhoBins(Int_t numberrhobins = 4,Int_t NotUselastbin = 1, TString filenameResp = "ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoW.root", TString filenameDataRho = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/RhoFromDataRatiosNPtDetbins2Deriv.root"){
	
	// filenameResp: note that also the output will be written in that directory, adding "Out" to the name
	
	//rho bins
	Int_t nrhobins = 0;
	Double_t drho = 0;
	if(numberrhobins == 5) {
		nrhobins = 5; 
		drho = 1.2;
	}
	if(numberrhobins == 4) {
		nrhobins = 4; 
		drho = 1.5;
	}
	if(numberrhobins == 2) {
		nrhobins = 2; 
		drho = 2.5;
	}
	if(numberrhobins == 1){
		nrhobins = 1;
		drho = 1.5;
	}
	Double_t rhobinlims[nrhobins+1];
	for(Int_t ibrho = 0; ibrho < nrhobins+1; ibrho++){
		rhobinlims[ibrho] = drho*ibrho;
	}
	//include "all" in the last bin
	rhobinlims[nrhobins] = 10.;

	if(numberrhobins == 1){
		rhobinlims[0] = 1.5;
		rhobinlims[1] = 3;
	}
	TFile *finR = new TFile(filenameResp);
	if(!finR->IsOpen()){
		Printf("File %s not found", filenameResp.Data());
		return;
	}
	
	TString hRespname = "fhResponseFinalRho";
	
	
	TFile *finD = new TFile(filenameDataRho);
	if(!finD->IsOpen()){
		Printf("File %s not found", filenameDataRho.Data());
		return;
	}
	TString hnameRhoDRatio = "hRhoDataR_Bin1";
	TH1D *hRhoFactor = (TH1D*)finD->Get(hnameRhoDRatio);
	if(!hRhoFactor){
		Printf("%s not found, might crash", hnameRhoDRatio.Data());
	}
	//fit in the range 0-6 with a pol1
	TF1 *fpol1 = new TF1("fpol1", "pol1", 0., 6.);
	hRhoFactor->Fit(fpol1, "RL0");
	fpol1->SetLineColor(hRhoFactor->GetLineColor());
	
	TCanvas *cRhoDataRatio = new TCanvas("cRhoDataRatio", "Rho ratio, factor to be applied to the response");
	cRhoDataRatio->cd();
	hRhoFactor->Draw();
	fpol1->Draw("sames");
	
	TCanvas *cDistrib = new TCanvas("cDistrib", "Distributions in rho bins", 800, 800);
	cDistrib->Divide(2,2);
	
	TCanvas *cDistribFinal = new TCanvas("cDistribFinal", "Distributions final response", 800, 800);
	cDistribFinal->Divide(2,2);
	
	TCanvas *cRatioDistribF = new TCanvas("cRatioDistribF", "Ratios of the final distributions wrt no weights", 800, 800);
	cRatioDistribF->Divide(2,2);
	
	TH1D *hRhoFactorApplied = new TH1D("hRhoFactorApplied", "Rho factor applied; #rho; Factor", nrhobins, rhobinlims);
	hRhoFactorApplied->SetLineColor(colors[2]);
	hRhoFactorApplied->SetMarkerColor(colors[2]);
	hRhoFactorApplied->SetMarkerStyle(25);
	
	TH1D *hRhoFactorFromFit = new TH1D("hRhoFactorFromFit", "Rho factor from fit; #rho; Factor", nrhobins, rhobinlims);
	hRhoFactorFromFit->SetLineColor(colors[3]);
	hRhoFactorFromFit->SetMarkerColor(colors[3]);
	hRhoFactorFromFit->SetMarkerStyle(24);
	
	
	THnSparseF *hRespInFile = (THnSparseF*)finR->Get("hresponsefinalRhoW");
	if(!hRespInFile) Printf("Resp in file named hresponsefinalRhoW not found");
	
	THnSparseF *hRespFacCalc1 = 0x0;
	THnSparseF *hRespFacFit1  = 0x0;
	THnSparseF *hRespXCheck   = 0x0;

	TLegend *legRho = new TLegend(0.4, 0.4, 0.7, 0.9);
	legRho->SetBorderSize(0);
	legRho->SetFillStyle(0);
	
	TLegend *legCorF = new TLegend(0.5, 0.6, 0.8, 0.9);
	legCorF->SetBorderSize(0);
	legCorF->SetFillStyle(0);
	
	
	for(Int_t ibrho = 0; ibrho < nrhobins-NotUselastbin; ibrho++){
		Double_t centralRhoValue = (rhobinlims[ibrho] + rhobinlims[ibrho+1]) * 0.5;
		
		Int_t binlimsRhoFactor[2] = { hRhoFactor->FindBin(rhobinlims[ibrho]), hRhoFactor->FindBin(rhobinlims[ibrho+1]- 0.001)};
		Int_t binRhoFactor = (binlimsRhoFactor[0] + binlimsRhoFactor[1]) * 0.5;
		Double_t respRhoFact = hRhoFactor->GetBinContent(binRhoFactor);
		
		hRhoFactorApplied->SetBinContent(ibrho+1, respRhoFact);
		Printf("Factor %d, %f", ibrho+1, centralRhoValue);
		Printf("Histo bin limits %f - %f", hRhoFactorApplied->GetBinLowEdge(ibrho+1), hRhoFactorApplied->GetBinLowEdge(ibrho+2));
		Printf("Bin limits: %f -> %d, %f -> %d\n Avg bin %d, content %f", rhobinlims[ibrho], binlimsRhoFactor[0], rhobinlims[ibrho+1]- 0.001, binlimsRhoFactor[1], binRhoFactor, respRhoFact);
		
		Double_t binCenterF = hRhoFactor->GetBinCenter(binRhoFactor);
		Double_t factorFromFit = fpol1->Eval(centralRhoValue);
		hRhoFactorFromFit->SetBinContent(ibrho+1, factorFromFit);
		Printf("BinCenter %f -> %d, fPol1 = %f",  centralRhoValue, ibrho+1,  factorFromFit);
		
		THnSparseF *hResp = (THnSparseF*)finR->Get(Form("%s%d", hRespname.Data(), ibrho));
		if(!hResp) {
			Printf("%s%d not found, skipping", hRespname.Data(), ibrho);	
			continue;
		}
		if(!hRespFacCalc1) {
			hRespFacCalc1 = (THnSparseF*)hResp->Clone(Form("hRespFacCalc1"));
			hRespFacCalc1->Scale(respRhoFact);
		} else {
			hRespFacCalc1->Add(hResp, respRhoFact);
		}
		if(!hRespFacFit1) {
			hRespFacFit1 = (THnSparseF*)hResp->Clone(Form("hRespFacFit1"));
			hRespFacFit1->Scale(respRhoFact);
		} else {
			hRespFacFit1->Add(hResp, factorFromFit);
		}
		
		if(!hRespXCheck){
			hRespXCheck = (THnSparseF*)hResp->Clone("hRespNoRhoW");
		} else {
			hRespXCheck->Add(hResp);
		}
		for(Int_t iax = 0; iax < 4; iax++){
			
			TH1D* hPj = hResp->Projection(iax);
			hPj->SetName(Form("hPjAx_%d_Rho_%d", iax, ibrho));
			hPj->SetLineColor(colors[ibrho]);
			hPj->SetLineWidth(2);
			cDistrib->cd(iax+1);
			gPad->SetLogy();
			if(ibrho == 0) hPj->Draw();
			else hPj->Draw("sames");
			
			if(iax == 0) {
				legRho->AddEntry(hPj, Form("%.1f< #rho < %.1f", rhobinlims[ibrho], rhobinlims[ibrho+1]), "l");
			}
		}
	}
	
	cRhoDataRatio->cd();
	hRhoFactorApplied->Draw("Psames");
	hRhoFactorFromFit->Draw("Psames");
	
	for(Int_t iax = 0; iax < 4; iax++){
		
		TH1D* hPj0 = 0x0;
		if(hRespInFile) {
			hPj0 = hRespInFile->Projection(iax);
			hPj0->SetName(Form("hPjInFileAx_%d", iax));
			hPj0->SetLineColor(colors[0]);
			hPj0->SetMarkerColor(colors[0]);
			hPj0->SetMarkerStyle(20);
		}
		
		TH1D* hPj1 = 0x0;
		if(hRespFacCalc1) { 
			hPj1 = hRespFacCalc1->Projection(iax);
			hPj1->SetName(Form("hPjFacCalc1Ax_%d", iax));
			hPj1->SetLineColor(colors[2]);
			hPj1->SetMarkerColor(colors[2]);
			hPj1->SetMarkerStyle(25);
		}
		
		TH1D* hPj2 = 0x0;
		if(hRespFacFit1){
			hPj2 = hRespFacFit1->Projection(iax);
			hPj2->SetName(Form("hPjFacFit1Ax_%d", iax));
			hPj2->SetLineColor(colors[3]);
			hPj2->SetMarkerColor(colors[3]);
			hPj2->SetMarkerStyle(24);
		}
		
		TH1D* hPj3 = 0x0;
		if(hRespXCheck){
			hPj3 = hRespXCheck->Projection(iax);
			hPj3->SetName(Form("hPjNoFacAx_%d", iax));
			hPj3->SetLineColor(colors[4]);
			hPj3->SetMarkerColor(colors[4]);
			hPj3->SetMarkerStyle(26);
		}

		cDistribFinal->cd(iax+1);
		gPad->SetLogy();
		if(hPj0) hPj0->Draw("P");
		if(hPj1) hPj1->Draw("Psames");
		if(hPj2) hPj2->Draw("Psames");
		if(hPj3) hPj3->Draw("Psames");
		
		if(iax == 0 ){
			//if(hPj0) legCorF->AddEntry(hPj0, "From file", "LP");
			if(hPj1) legCorF->AddEntry(hPj1, "Mean-bin content", "LP");
			if(hPj2) legCorF->AddEntry(hPj2, "Eval from pol1", "LP");
			if(hPj3) legCorF->AddEntry(hPj3, "No Factor", "LP");
		}
		
		//ratios:
		
		TH1D* hRatioPj0 = 0x0;
		if(hPj0)  hRatioPj0 = (TH1D*)hPj0->Clone(Form("hPjInFileAx_%d", iax));

		TH1D* hRatioPj1 = 0x0;
		if(hPj1)  hRatioPj1 = (TH1D*)hPj1->Clone(Form("hRatioPjFacCalc1Ax_%d", iax));
		
		TH1D* hRatioPj2 = 0x0;
		if(hPj2)  hRatioPj2 = (TH1D*)hPj2->Clone(Form("hPjFacFit1Ax_%d", iax));
		if(hPj3){
		if(hRatioPj0) hRatioPj0->Divide(hPj3);
		if(hRatioPj1) hRatioPj1->Divide(hPj3);
		if(hRatioPj2) hRatioPj2->Divide(hPj3);
		}
		
		if(hRatioPj0) hRatioPj0->GetYaxis()->SetRangeUser(0, 10);
		if(hRatioPj1) hRatioPj1->GetYaxis()->SetRangeUser(0, 10);
		if(hRatioPj2) hRatioPj2->GetYaxis()->SetRangeUser(0, 10);
		
		cRatioDistribF->cd(iax+1);
		gPad->SetGridy();
		if(hRatioPj0) hRatioPj0->Draw("P");
		if(hRatioPj1) hRatioPj1->Draw("Psames");
		if(hRatioPj2) hRatioPj2->Draw("Psames");
		
	}
	cDistribFinal->cd(1);
	legCorF->Draw();
	cDistrib->cd(1);
	legRho->Draw();
	cRhoDataRatio->cd();
	legCorF->Draw();
	cRatioDistribF->cd(1);
	legCorF->Draw();
	
	SaveCv(cDistribFinal);
	SaveCv(cDistrib);
	SaveCv(cRhoDataRatio);
	SaveCv(cRatioDistribF);
	
	TString newname = filenameResp(0, filenameResp.Sizeof()-6);
	
	TFile *fout = new TFile(Form("%sOut.root", newname.Data()), "recreate");
	fout->cd();
	if(hRespFacCalc1) hRespFacCalc1->Write();
	if(hRespFacFit1 ) hRespFacFit1 ->Write();
	if(hRespXCheck  ) hRespXCheck  ->Write();
	
	fout->Close();
	Printf("Written file %sOut.root", newname.Data());
	
}

//_____________________________________________________________________________
void CombineRhoBinsOvlExcluAndHighRho(){

	Int_t nrhobins = 4;
	Double_t rhobinlims[nrhobins+1] = {0., 1.5, 3., 4.5, 7};
	TString inputFiles[nrhobins] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160707/OvlExclu/ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoW.root", "", "", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160707/OvlExclu/HighRho/ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoW.root"};
	
	Int_t composition[nrhobins] = {0, 0, 0, 3};
	
	Int_t rhoindexinfile[nrhobins] = {0, 1, 2, 0};
	
	TString description = "OvlExcluAndHighRho";
	
	CombineRhoBins(nrhobins, rhobinlims, inputFiles, composition, rhoindexinfile, description);
}

//_____________________________________________________________________________

void CombineRhoBins(Int_t nrhobins, Double_t *rhobinlims, TString *inputFiles, Int_t *composition, Int_t *rhoindexinfile, TString description, TString filenameDataRho){
	
	// inputFiles is an array of nrhobins. Only indeces where the input file changes are filled
	// composition: contains the number of the filled inputFiles propagated for all the rho bins that get input from that entry. E.g. nrhobins = 5; inputFiles[5] = {"file1", "", "", "file2", "file3"}; composition[5] = {0, 0, 0, 3, 4}
	// rhoindexinfile: which rho bin has to be taken relativly to the file that contains the non-weighted responses. Delicate, requires manually opening and looking into the files
	
	//read ratio rho pT bin regions
	
	TFile *finD = new TFile(filenameDataRho);
	if(!finD->IsOpen()){
		Printf("File %s not found", filenameDataRho.Data());
		return;
	}
	TString hnameRhoDRatio = "hRhoDataR_Bin1";
	TH1D *hRhoFactor = (TH1D*)finD->Get(hnameRhoDRatio);
	if(!hRhoFactor){
		Printf("%s not found, might crash", hnameRhoDRatio.Data());
	}
	//fit in the range 0-6 with a pol1
	TF1 *fpol1 = new TF1("fpol1", "pol1", 0., 6.);
	hRhoFactor->Fit(fpol1, "RL0");
	fpol1->SetLineColor(hRhoFactor->GetLineColor());
	
	TCanvas *cRhoDataRatio = new TCanvas("cRhoDataRatio", "Rho ratio, factor to be applied to the response");
	cRhoDataRatio->cd();
	hRhoFactor->Draw();
	fpol1->Draw("sames");
	
	TH1D *hRhoFactorFromFit = new TH1D("hRhoFactorFromFit", "Rho factor from fit; #rho; Factor", nrhobins, rhobinlims);
	hRhoFactorFromFit->SetLineColor(colors[3]);
	hRhoFactorFromFit->SetMarkerColor(colors[3]);
	hRhoFactorFromFit->SetMarkerStyle(24);
	
	TCanvas *cDistribFinal = new TCanvas("cDistribFinal", "Distributions final response", 800, 800);
	cDistribFinal->Divide(2,2);
	
	TFile *finR = 0x0;
	TString hRespname = "fhResponseFinalRho";
	
	THnSparseF *hResponseFinal = 0x0;
	
	TString filename = TString(inputFiles[0], inputFiles[0].Sizeof() - 6);
	filename+=description;
	TFile *fout = new TFile(Form("%s.root", filename.Data()), "recreate");
	
	for(Int_t ibrho = 0; ibrho < nrhobins; ibrho++){
		
		//determination of the factor
		Int_t binlimsRhoFactor[2] = { hRhoFactor->FindBin(rhobinlims[ibrho]), hRhoFactor->FindBin(rhobinlims[ibrho+1]- 0.001)};
		Int_t binRhoFactor = (binlimsRhoFactor[0] + binlimsRhoFactor[1]) * 0.5;  
		Double_t respRhoFact = hRhoFactor->GetBinContent(binRhoFactor);
		
		//using the fit
		Double_t centralRhoValue = (rhobinlims[ibrho] + rhobinlims[ibrho+1]) * 0.5;
		Double_t factorFromFit = fpol1->Eval(centralRhoValue);
		hRhoFactorFromFit->SetBinContent(ibrho+1, factorFromFit);
		Printf("BinCenter %f -> %d, fPol1 = %f",  centralRhoValue, ibrho+1,  factorFromFit);
		
		
		if(!finR || (composition[ibrho] != composition[ibrho-1])){
			delete finR;
			Printf("Reading file %s", inputFiles[composition[ibrho]].Data());
			finR = new TFile(inputFiles[composition[ibrho]]);
			if(!finR->IsOpen()){
				Printf("File not found, can't proceed, but go to next");
				continue;
			}
			
		}
		//read the response
		THnSparseF *hResp = (THnSparseF*)finR->Get(Form("%s%d", hRespname.Data(), rhoindexinfile[ibrho]));
		if(!hResp) {
			Printf("%s%d not found, skipping", hRespname.Data(), ibrho);	
			continue;
		}
		hResp->SetName(Form("%s_InF%d", hResp->GetName(), composition[ibrho]));
		if(!hResponseFinal) {
			hResponseFinal = (THnSparseF*)hResp->Clone(Form("hResponseFinal"));
			hResponseFinal->Scale(respRhoFact);
		} else {
			hResponseFinal->Add(hResp, factorFromFit);
		}
		
		fout->cd();
		hResp->Write();
		
	}// loop on rho bins
	
	cRhoDataRatio->cd();
	hRhoFactorFromFit->Draw("Psames");
	
	SaveCv(cRhoDataRatio);
	
	TCanvas *cAxPtBins[4];
	Int_t axisRange = 3; //pt par 
	TLegend *legPt = new TLegend(0.6, 0.2, 0.9, 0.4);
	legPt->SetBorderSize(0);
	legPt->SetFillStyle(0);
   
	for(Int_t iax = 0; iax < 4; iax++){
		if(iax == axisRange) continue;
		cAxPtBins[iax] = new TCanvas(Form("cAx%dPtBins", iax), Form("Projections in pT bins, Axis %s", hResponseFinal->GetAxis(iax)->GetTitle()), 600, 500);
		
		TH1D* hpj = GetProjection(hResponseFinal, iax, Form("hFinal%s", description.Data()), 0, 20);
		
		cDistribFinal->cd(iax + 1);
		gPad->SetLogy();
		hpj->Draw();
		
		for(Int_t ipt = 0; ipt < nptbins; ipt++){
			Int_t bin[2] = {hResponseFinal->GetAxis(axisRange)->FindBin(ptlims[ipt]), hResponseFinal->GetAxis(axisRange)->FindBin(ptlims[ipt+1]-0.01)};
			hResponseFinal->GetAxis(axisRange)->SetRange(bin[0], bin[1]);
			
			TH1D* hpjRange = GetProjection(hResponseFinal, iax, Form("hFinal%sPtBin%d", description.Data(), ipt), ipt+1, 20);
			
			if(iax == 0){
				legPt->AddEntry(hpjRange, Form("%.0f < %s < %.0f", ptlims[ipt],  hResponseFinal->GetAxis(axisRange)->GetTitle(), ptlims[ipt+1]));
			}
			cAxPtBins[iax]->cd();
			gPad->SetLogy();
			if(iax == axisRange) hpjRange->GetXaxis()->SetRange(0, -1);
			if(ipt == 0) hpjRange->Draw();
			else hpjRange->Draw("sames");
			
			hResponseFinal->GetAxis(axisRange)->SetRange(0, -1);
		}
		cAxPtBins[iax]->cd();
		legPt->Draw();
		SaveCv(cAxPtBins[iax]);
	}
	
	fout->cd();
	hResponseFinal->Write();
	hRhoFactorFromFit->Write();
	
	

}

//_____________________________________________________________________________
void CompareRespRhoWNoRhoW(){
	Int_t nresp = 2; 
	TString files[nresp] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160707/OvlExclu/ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoWOvlExcluAndHighRho.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	TString hspname[nresp] = {"hResponseFinal", "fhResponseFinal"};
	TString legs[nresp] = {"DetFlucDerivRhoW", "DetFlucDeriv"};
	Int_t baseRatios = 1;
	ComparisonResponses(nresp, files, hspname, legs, baseRatios);
	return;
}

//_____________________________________________________________________________

void ComparisonResponses(Int_t nresp, TString *files, TString *hspname, TString *legs, Int_t baseRatios){
	
	TCanvas *cAxPtBins[4];
	TCanvas *cAxPtBinsRatios[4];
	Int_t axisRange = 3; //pt par 
	
	TLegend *leg = new TLegend(0.7, 0.3, 0.9, 0.4);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nptbins, nx, ny, dx, dy);
   
	for(Int_t iax = 0; iax < 4; iax++) {
		cAxPtBins[iax] = new TCanvas(Form("cAx%dPtBins", iax), "Projections", 900, 900);
		cAxPtBins[iax]->Divide(nx, ny);
		
		cAxPtBinsRatios[iax] = new TCanvas(Form("cAx%dPtBinsRatios", iax), "Ratios", 600, 500);
		//cAxPtBinsRatios[iax]->Divide(nx, ny);
	}
	
	
	TH1D *hRatios[4][nresp][nptbins];
	for(Int_t iresp = 0; iresp < nresp; iresp++){
		for(Int_t iax = 0; iax < 4; iax++){
			for(Int_t ipt = 0; ipt < nptbins; ipt++){
				hRatios[iax][iresp][ipt] = 0x0;
			}
		}
	}
	
	for(Int_t iresp = 0; iresp < nresp; iresp++){
		
		TFile *fin = new TFile(files[iresp]);
		
		if(!fin->IsOpen()){
			Printf("%s not found", files[iresp].Data());
			continue;
		}
		
		THnSparseF* hResp = (THnSparseF*)fin->Get(hspname[iresp]);
		if(!hResp){
			Printf("%s not found", hspname[iresp].Data());
			continue;
		}
		
		const char* axRangetitle = hResp->GetAxis(axisRange)->GetTitle();
		
		TCanvas *cDistribFinal = new TCanvas(Form("cDistribFinal%s", legs[iresp].Data()), Form("Distributions final response %s", legs[iresp].Data()), 800, 800);
		cDistribFinal->Divide(2,2);
		//projecting and plotting
		
		for(Int_t iax = 0; iax < 4; iax++){
			
			
			TH1D* hpj = GetProjection(hResp, iax, Form("hFinal%s", legs[iresp].Data()), 0, 20);
			
			cDistribFinal->cd(iax + 1);
			gPad->SetLogy();
			hpj->Draw();
			SaveCv(cDistribFinal);
			
			//if(iax == 0) leg->AddEntry(hpj, legs[iresp], "P");
			
			cAxPtBins[iax]->SetTitle(Form("Projection Axis %s", hResp->GetAxis(iax)->GetTitle()));
			cAxPtBinsRatios[iax]->SetTitle(Form("Ratios Axis %s", hResp->GetAxis(iax)->GetTitle()));
			if(iax == axisRange) {
				hRatios[iax][iresp][0] = (TH1D*)hpj->Clone(Form("hRatioResp%sPj%dAllPt", legs[iresp].Data(), iax));
				continue;
			}
			for(Int_t ipt = 0; ipt < nptbins; ipt++){
				
				TPaveText *pv = new TPaveText(0.3, 0.8, 0.7, 0.9, "NDC");
				pv->SetFillStyle(0);
				pv->SetBorderSize(0);
				pv->AddText(Form("%.0f < %s < %.0f", ptlims[ipt],  axRangetitle, ptlims[ipt+1]));
				
				Int_t bin[2] = {hResp->GetAxis(axisRange)->FindBin(ptlims[ipt]), hResp->GetAxis(axisRange)->FindBin(ptlims[ipt+1]-0.01)};
				hResp->GetAxis(axisRange)->SetRange(bin[0], bin[1]);
				
				TH1D* hpjRange = GetProjection(hResp, iax, Form("hFinal%sPtBin%d", legs[iresp].Data(), ipt), ipt+1, 20+iresp);
				
				//copy the histogram in another object used later to perform the ratios
				hRatios[iax][iresp][ipt] = (TH1D*)hpjRange->Clone(Form("hRatioResp%sPj%dPtB%d", legs[iresp].Data(), iax, ipt));
				
				
				cAxPtBins[iax]->cd(ipt+1);
				gPad->SetLogy();
				if(iresp == 0) hpjRange->Draw();
				else hpjRange->Draw("sames");
				pv->Draw();
				
				hResp->GetAxis(axisRange)->SetRange(0, -1);
			}
			cAxPtBins[iax]->cd();
			//leg->Draw();
			SaveCv(cAxPtBins[iax]);
		}
		
		//delete fin;
	}
	
	//loop again over axes and pt bins to perform the ratios
	
	
	for(Int_t iax = 0; iax < 4; iax++){
		for(Int_t iresp = 0; iresp < nresp; iresp++){
			
			if(iresp == baseRatios) continue;
			
			for(Int_t ipt = 0; ipt < nptbins; ipt++){
				if(!hRatios[iax][iresp][ipt]) continue;
				hRatios[iax][iresp][ipt]->Divide(hRatios[iax][baseRatios][ipt]);
				cAxPtBinsRatios[iax]->cd();
				if(ipt == 0) hRatios[iax][iresp][ipt]->Draw();
				else hRatios[iax][iresp][ipt]->Draw("sames");
			}
			
		}
		SaveCv(cAxPtBinsRatios[iax]);
	}
	return;
}

//_____________________________________________________________________________

TH1D* GetProjection(THnSparseF *hresp, Int_t iax, TString basenamepj, Int_t colid, Int_t marker) {
	if(!hresp) {
		Printf("Input is null");
		return 0x0;
	}
	
	TH1D* hPj0 = hresp->Projection(iax);
	hPj0->SetName(Form("%s_%d", basenamepj.Data(), iax));
	hPj0->SetLineColor(colors[colid]);
	hPj0->SetMarkerColor(colors[colid]);
	hPj0->SetMarkerStyle(marker);
	
	return hPj0;
	
}
