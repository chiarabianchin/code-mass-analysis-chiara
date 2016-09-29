#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>

#include </data/Work/MyCodeJetMass/macros/Convolution.C>

Bool_t ReadHistograms(TH2F*& hpTMrec, TH2F*& hdMdpT, TString file1 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/UnfoldingMatrix.root", TString file2 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/analysis/UnfoldingThnSparse.root", TString hnameRec = "fhnMassResponse_proj_0_2", TString hnameFluc = "fhnDeltaMass_0_proj_1_0");

Bool_t ReadHistograms(TH2F*& hpTMrec, TH2F*& hdMdpT, TString file1, TString file2, TString hnameRec, TString hnameFluc){
   
   TFile *fin1 = new TFile(file1);
   if(!fin1->IsOpen()){
      Printf("Input file %s not found", file1.Data());
      return kFALSE;
   }
   
   hpTMrec = (TH2F*)fin1->Get(hnameRec);
   if(!hpTMrec) {
      Printf("W: %s not found!", hnameRec.Data());
      return kFALSE;
   }

   TFile *fin2 = new TFile(file2);
   if(!fin2->IsOpen()){
      Printf("Input file %s not found", file2.Data());
      return kFALSE;
   }
   
   hdMdpT = (TH2F*)fin2->Get(hnameFluc);
   if(!hdMdpT) {
      Printf("W: %s not found!", hnameFluc.Data());
      return kFALSE;
   }
   
   return kTRUE;
}

void PractiseMultiplyMatrices(){
   TH2F *hpTMrec = 0x0;
   TH2F *hdMdpT  = 0x0;
   if(!ReadHistograms(hpTMrec, hdMdpT)) return;
   
   TH2F* hwhat = Convolution2D(hpTMrec, hdMdpT);
   TCanvas *c = new TCanvas("c", "?");
   c->cd();
   hwhat->Draw("colz");
}
