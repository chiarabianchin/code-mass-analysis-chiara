#include "Riostream.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include "TGraph.h"
#include <TSystem.h>
#include <TObjString.h>

void MergeEmbeddingThrmModelBins(TString listname = "JetShapeConst_Jet_AKTChargedR040_ThrmTracksEmbSingle_pT0150_E_scheme_TCMCMatch", TString inputWeights = "/data/Work/jets/JetMass/BkgFluctStudies/CleanEmbedding/WeightFactors.root", TString graphWname = "gfractionsvsNch", TString inputfiletxt = "/data/Work/jets/JetMass/BkgFluctStudies/CleanEmbedding/SingleTrackInPYTHIAPlusThermal/inputfiles.txt", TString listwithEventsForSure = "BackFlucRandomCone", TString heventname = "fHistEventCount"){
   
   TFile *fw = new TFile(inputWeights);
   if(!fw->IsOpen()){
      Printf("File weights not found");
      return;
   
   }
   TGraph *gweight = (TGraph*) fw->Get(graphWname);
   if(!gweight){
      Printf("Graph with weights not found");
      return;
   
   }
   Printf("The area of the wright graph is %f", gweight->Integral());
   static const Int_t nBins = gweight->GetN();
   Double_t factors[nBins];
   for(Int_t i = 0; i<nBins; i++){
      Double_t x;
      gweight->GetPoint(i, x, factors[i]);
      Printf("Weight %d = %f", i, factors[i]);
   
   }
  
   ifstream input(inputfiletxt);
   TString path;
   Int_t count = 0;
    //whatch out, the name might vary from list to list, find a better way to do it
   Double_t nevents = 0;
   Double_t weight = 1.;

   TList* finalList = 0x0;
   
   while(input){
      input>>path;
      if(path.IsNull()) break;
      Printf("Reading %s", path.Data());
      
      TFile *f = new TFile(path);
      if(!f->IsOpen()){
      	 Printf("File not found");
      	 return;
      
      }
      
      TList *list = (TList*)f->Get(listname); // this is the TList considered for merging
      if(!list) {
      	 Printf("List %s not found, exit", listname.Data());
      	 f->ls();
      	 return;
      }
      Int_t nentrieslist = list->GetEntries();
      
      TH1F *hevents = (TH1F*)list->FindObject(heventname);
      if(!hevents) {
      	 Printf("Number of events histograms %s not found, take it from %s", heventname.Data(), listwithEventsForSure.Data());
      	 TList *listNevents = (TList*)f->Get(listwithEventsForSure);
      	 hevents = (TH1F*)listNevents->FindObject(heventname);
      	 if(!hevents) {
      	    Printf("Histogram with number of event NOT FOUND, check the last 2 inputs!!!");
      	    continue;
      	 }
      }
      nevents = hevents->GetBinContent(1);
      if(nevents > 0){
      	 weight = factors[count]/nevents;
      	 Printf("Weight for %s is %f", path.Data(), weight);
      }
      
      if(count == 0){
      	 finalList = static_cast<TList*>(list->Clone(Form("%s",list->GetName())));
      	 TH1 *h1 = 0x0;
      	 TH2 *h2 = 0x0;
      	 TH3 *h3 = 0x0;
      	 THnSparse *hn = 0x0; 
      	 TObjArray *fObjArrayHistoNames = new TObjArray(); //is this needed???? try to remove it
      	 fObjArrayHistoNames->SetOwner(kTRUE);
      	 
      	 for (Int_t j=0; j<finalList->GetEntries(); ++j) {
      	    TObjString *ostr = 0x0;
      	    if(finalList->At(j)->InheritsFrom("THnSparse")) {
      	       hn = dynamic_cast<THnSparse*>(finalList->At(j));
      	       if(hn) {
      	       	  hn->Sumw2();
      	       	  hn->Reset();
      	       	  ostr = new TObjString(hn->GetName());
      	       	  fObjArrayHistoNames->Add(ostr);
      	       	  continue;
      	       }
      	    }
      	    else if(finalList->At(j)->InheritsFrom("TH3")) {
      	       h3 = dynamic_cast<TH3*>(finalList->At(j));
      	       if(h3) {
      	       	  h3->Sumw2();
      	       	  h3->Reset();
      	       	  ostr = new TObjString(h3->GetName());
      	       	  fObjArrayHistoNames->Add(ostr);
      	       	  continue;
      	       }
      	    }
      	    else if(finalList->At(j)->InheritsFrom("TH2")) {
      	       h2 = dynamic_cast<TH2*>(finalList->At(j));
      	       if(h2) {
      	       	  h2->Sumw2();
      	       	  h2->Reset();
      	       	  ostr = new TObjString(h2->GetName());
      	       	  fObjArrayHistoNames->Add(ostr);
      	       	  continue;
      	       }
      	    }
      	    else if(finalList->At(j)->InheritsFrom("TH1")) {
      	       h1 = dynamic_cast<TH1*>(finalList->At(j));
      	       if(h1) { 
      	       	  h1->Sumw2();
      	       	  h1->Reset();
      	       	  ostr = new TObjString(h1->GetName());
      	       	  fObjArrayHistoNames->Add(ostr);
      	       	  continue;
      	       }
      	    }
      	 }
      	 
      }
      
      TH1F * histo1d = 0x0;
      TH2F * histo2d = 0x0;
      TH3F * histo3d = 0x0; 
      THnSparse *histosparse = 0x0;
      Printf("Now adding...");

      for(Int_t j = 0; j<nentrieslist; j++){

      	 if(finalList->At(j)->InheritsFrom("TH1F")){
      	    histo1d = (TH1F*)finalList->At(j);
      	    TH1F* htmp = (TH1F*)list->At(j);
      	    htmp->Sumw2();
      	    htmp->Scale(weight);
      	    histo1d->Add(htmp);
      	 }
      	 
      	 if(finalList->At(j)->InheritsFrom("TH2F")){
      	    histo2d = (TH2F*)finalList->At(j);
      	    TH2F* htmp = (TH2F*)list->At(j);
      	    htmp->Sumw2();
      	    htmp->Scale(weight);
      	    histo2d->Add(htmp);
      	    
      	 }
      	 if(finalList->At(j)->InheritsFrom("TH3F")){
      	    histo3d = (TH3F*)finalList->At(j);
      	    TH3F* htmp = (TH3F*)list->At(j);
      	    htmp->Sumw2();
      	    htmp->Scale(weight);
      	    histo3d->Add(htmp);
      	    
      	 }

      	 if(finalList->At(j)->InheritsFrom("THnSparse")){
      	    histosparse = (THnSparse*)finalList->At(j);
      	    THnSparse* htmp = (THnSparse*)list->At(j);
      	    htmp->Sumw2();
      	    htmp->Scale(weight);
      	    histosparse->Add(htmp);
      	 }
      }
      
      
      count++;
   }

   Printf("Writing file");
   gSystem->MakeDirectory("merge");
   
   TFile *fout = new TFile(Form("merge/%s.root", listname.Data()), "recreate");
   fout->WriteTObject(finalList,finalList->GetName()); //use to get the list written, not only the histograms
   fout->Write();
   fout->Close();
   
}

