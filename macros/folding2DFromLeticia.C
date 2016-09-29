#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TKey.h>
//#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>


void folding2DFromLeticia(TString pathRaw = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/", TString pathMtx = "", TString pathoutput = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train921/analysis1/"){
   TString filenameRaw = "MassvspT.root";
   filenameRaw.Prepend(pathRaw);
   TFile *fRaw = new TFile(filenameRaw);
   if(!fRaw->IsOpen()){
      Printf("File %s not found", filenameRaw.Data());
      return;
   }
   TString hRawname = "hMasspT";
   TH2D* htrue = (TH2D*)fRaw->Get(hRawname);
   if(!htrue){
      Printf("%s not found in %s", hRawname.Data(), filenameRaw.Data());
      fRaw->ls();
      return;
   }
   TH2D *hfold = (TH2D*)htrue->Clone("hfold"); 
   
   TString filenameMtx = "UnfoldingThnSparse.root";
   filenameMtx.Prepend(pathMtx);
   TFile *fMtx = new TFile(filenameMtx);
   if(!fMtx->IsOpen()){
      Printf("File %s not found", filenameMtx.Data());
      return;
   }
   TString hMtxname = "fhnDeltaMass_0_proj_2_3_4_5";
   THnSparse *hresponse = (THnSparse*) fMtx->Get(hMtxname);
   if(!hresponse){
      Printf("%s not found in %s", hMtxname.Data(), filenameMtx.Data());
      fMtx->ls();
      return;
   }
   
   const Int_t n = 4;
   Int_t nbins[n] = {hresponse->GetAxis(0)->GetNbins(), hresponse->GetAxis(1)->GetNbins(), hresponse->GetAxis(2)->GetNbins(), hresponse->GetAxis(3)->GetNbins()};
   
   Printf("Response matrix:");
   for(int i=0;i<n;i++){
      Printf("\n- axis %s (%d) %d bins", hresponse->GetAxis(i)->GetTitle(), i, nbins[i]);
   }
   
   for(int i=0;i<nbins[0];i++){
      for(int j=0;j<nbins[2];j++){
      	 double effects=0;
      	 for(int k=0;k<nbins[1];k++){
      	    for(int l=0;l<nbins[3];l++){
      	       Int_t binresp[n] = {i+1,k+1,j+1,l+1};
      	       
      	       effects=effects+htrue->GetBinContent(k+1,l+1)*hresponse->GetBinContent(binresp);
      	       
      	    }
      	 }
      	 hfold->SetBinContent(i+1,j+1,effects);
      }
   }
   
   TCanvas *cOutputCheck = new TCanvas("cOutputCheck", "Project Raw/Unfolded distributions", 400, 800);
   cOutputCheck->Divide(1,2);
   
   TH1D* hpTRaw = (TH1D*) htrue->ProjectionX("hpTRaw");
   TH1D* hpTUnf = (TH1D*) hfold->ProjectionX("hpTUnf");
   TH1D* hMRaw  = (TH1D*) htrue->ProjectionY("hMRaw");
   TH1D* hMUnf  = (TH1D*) hfold->ProjectionY("hMUnf");
   
   hpTRaw->SetLineColor(kRed+2);
   hpTRaw->SetLineWidth(2);
   hMRaw ->SetLineColor(kRed+2);
   hMRaw ->SetLineWidth(2);
   
   hpTUnf->SetLineColor(kGreen+2);
   hpTUnf->SetLineWidth(2);
   hMUnf ->SetLineColor(kGreen+2);
   hMUnf ->SetLineWidth(2);
   
   cOutputCheck->cd(1);
   gPad->SetLogy();
   hpTRaw->Draw();
   hpTUnf->Draw("sames");
   
   cOutputCheck->cd(2);
   hMRaw->Draw();
   hMUnf->Draw("sames");
   
}
