#include <iostream>
#include <fstream>
using namespace std;
#include <TF1.h>
#include <TList.h>
#include <TFile.h>
#include <TKey.h>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

#include </data/Work/CodeJetMass/PtHardUtil/PtHardBinUtilities.h>

TGraph* DefineFactorsForMerging(TString fileDataNtracks = "/data/Work/jets/JetMass/testpTTaskSample/AnalysisResults.root", TString listname="AliAnalysisTaskEmcalJetSample_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub_Rho_EMCAL_histos");
void MergeToyModelNchBins(TString pathToytxt = "pathToyModel.txt", TString liststxt = "lists.txt");

void MergeToyModelNchBins(TString pathToytxt, TString liststxt){
   //init weights
   TGraph *gFactors = DefineFactorsForMerging();
   Int_t nbins = gFactors->GetN();
   Double_t weights[nbins];
   Printf("Wights:");
   for(Int_t i = 0; i<nbins; i++){
      Double_t x,y;
      gFactors->GetPoint(i, x, y);
      weights[i] = y;
      Printf("%d) %e", i, weights[i]);
   }
   
   TString linefile = "";
   Bool_t masstask = kFALSE;
   TString linel="JetByJetCorrectionOutput";
   
   ifstream read(liststxt.Data());
   while (read>>linefile){
      Printf("line: %s",linefile.Data());
      
      PtHardBinUtilities *ptHardUtil = new PtHardBinUtilities();
      ptHardUtil->SetNameHistosDir("");
      ptHardUtil->SetNameHistosList(linefile);
      if(masstask){
      	 ptHardUtil->SetNameHistosDir(linefile.Data());
      	 ptHardUtil->SetNameHistosList(linel.Data()); //for double linked lists
      }
      
      ptHardUtil->SetArrWeights(nbins, weights);
      ptHardUtil->InitFileList(pathToytxt);
      ptHardUtil->DoWeighting();
      
      Printf("Weighting done, store output in root file");
      
      TFile *fout = new TFile(Form("%s.root",linefile.Data()),"RECREATE");
      TList *listOut = ptHardUtil->GetOutputList();
      listOut->Print();
      fout->WriteTObject(listOut,listOut->GetName());
      //listOut->Write();
      fout->Write();
      fout->Close();
      
      ptHardUtil->CleanMemory();
   }
}

TGraph* DefineFactorsForMerging(TString fileDataNtracks, TString listname){
   Printf("Define weights...");
   TList *list = ReadFile(fileDataNtracks, listname);
   if(!list) return 0x0;
   TString hname = "fHistNTracks_0";
   TH1F *hNtracks = (TH1F*)list->FindObject(hname);
   if(!hNtracks){
      Printf("%s not found", hname.Data());
      list->ls();
      return 0x0;
   }
   Int_t nbins = hNtracks->GetNbinsX();
   Double_t maxtr = hNtracks->GetBinLowEdge(nbins+1);
   Printf("Limit %f", maxtr);
   TH1F *hNtracksg0 = new TH1F("hNtracksg0", "Number of tracks (>0); NTracks; Entries",nbins - 1 , 1., maxtr);
   for(Int_t i=0; i<nbins-1 ; i++){
      //Printf("Bin %d : %f error %f", i+1, hNtracks->GetBinContent(i+2), hNtracks->GetBinError(i+2));
      hNtracksg0->SetBinContent(i+1, hNtracks->GetBinContent(i+2));
      hNtracksg0->SetBinError(i+1, hNtracks->GetBinError(i+2));
   }
   
   hNtracksg0->Scale(1./hNtracksg0->Integral());
   TCanvas *cNch = new TCanvas("cNch", "Track distribution (normalized to 1)", 800, 800);
   hNtracksg0->Draw();
   SaveCv(cNch);
   const Int_t nMultBins = 7;
   Int_t Nch[nMultBins] = {45, 36.2, 30.5, 23.2, 16.1, 9.8, 4.4};
   TGraph* gFactors = new TGraph(nMultBins);
   gFactors->SetMarkerStyle(20);
   gFactors->SetTitle("Fraction of event with N tracks; N_{tracks}; Fraction");
   gFactors->SetName("gFactors");
   for(Int_t ibin = 0; ibin<nMultBins; ibin++){
      Int_t b = hNtracksg0->FindBin(Nch[ibin]);
      Double_t factor = hNtracksg0->GetBinContent(b);
      gFactors->SetPoint(ibin, Nch[ibin], factor);
   
   }
   
   TCanvas *cFactors = new TCanvas("cFactors", "Factors per multiplicity bin", 800, 800);
   cFactors->cd();
   gFactors->Draw("AP");
   SaveCv(cFactors);
   
   Printf("Done!");
   return gFactors;
   
}

void WritepathToytxt(){
   TString base =
   "/data/Work/jets/JetMass/BkgFluctStudies/EmbedSinglePart/ExcludeOverlap/";
   //"/data/Work/jets/JetMass/BkgFluctStudies/EmbedSinglePart/ExcludeOverlap/JetWithSingleTrack/";
   //"/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/";
   Int_t n = 7;
   TString path="";
   ofstream out("pathToyModel.txt", ios::out);
   for(Int_t i = 0; i< n; i++){
      path = Form("%s%02d/output/AnalysisResults.root\n", base.Data(), i+1);
      out<<path;
   }
}

void Writeliststxt(TString input){
   TString path="";
   ofstream out("lists.txt", ios::out);

   TFile *f = new TFile(input);
   if(!f->IsOpen()){
      Printf("File with lists not found");
      return;
   }
   TList *lobj = f->GetListOfKeys();
   Int_t n = lobj->GetEntries();
   Printf("%d lists", n);
   for(Int_t i=0 ; i< n; i++){
      TKey *l = (TKey*)lobj->At(i);
      out<<l->GetName()<<endl;
   }
}
