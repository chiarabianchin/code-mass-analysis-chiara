#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

void DrawJetSampleTask(TString listname = "SampleJet_JetRhosub_AKTChargedR040_PicoTracksTree_pT0150_E_scheme_TpcRhoSp_EXLJ_TPC_histos", TString filename = "AnalysisResults.root"){
   TList *list = ReadFile(filename, listname);
   if(!list) return;
   const Int_t nh2 = 2;
   TString h2name[nh2] = {"fHistJetsPhiEta_0", "fHistJetsPtLeadHad_0"};
   
   TCanvas *c2[nh2];
   
   for(Int_t i = 0; i< nh2; i++){
      TH2F *h2 = dynamic_cast<TH2F*>(list->FindObject(h2name[i]));
      if(!h2){
      	 Printf("%s not found", h2name[i].Data());
      	 continue;
      }
      c2[i] = new TCanvas(Form("c2_%d",i), Form("%s", h2name[i].Data()), 1000, 500);
      c2[i] ->Divide(2,1);
      c2[i] ->cd(1);
      gPad->SetLogy();
      h2->ProjectionX()->Draw("E");
      c2[i] ->cd(2);
      gPad->SetLogy();
      h2->ProjectionY()->Draw("E");
      
   }
   return;
}
