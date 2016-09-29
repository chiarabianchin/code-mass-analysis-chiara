#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
#include <TF1.h>

void DefineTH1FInputEmbedding(TString inputMeanpTVsNch = "/data/Work/jets/JetMass/BkgFluctStudies/AvgpTVsNCh/GraphAvgpTvsMult.root" , TString graphname= "gAvgpTvsMult"){
   const Int_t nNchBins = 7;
   Double_t Nch[nNchBins] = {45, 36.2, 30.5, 23.2, 16.1, 9.8, 4.4};
   
   TFile* fin = new TFile(inputMeanpTVsNch);
   if(!fin->IsOpen()){
      Printf("Check input file name");
      return;
   }
   
   TGraphAsymmErrors *gpT = (TGraphAsymmErrors*) fin->Get(graphname);
   
   if(!gpT){
      Printf("Check input graph name");
      return;
   }

   printf("Double_t AvgpTArray[%d] = {", nNchBins);
   for(Int_t i = 0; i<nNchBins; i++){
      Double_t x, y;
      Int_t point = gpT->GetXaxis()->FindBin(Nch[i]);
      gpT->GetPoint(point, x, y);
      if(i<nNchBins-1) printf("%f, ", y);
      else Printf("%f}", y);   
      
   }
   Printf("Remember that this is in eta<0.3, pT 0.150-10 GeV/c");
   //printf in a form of array
   

}

void FitpTDistribution(TString inputpTDistr = "/data/Work/jets/JetMass/BkgFluctStudies/pTdistribution/GraphpTChpPb.root", TString graphname= "gpTDistrpPb"){

   TFile* fin = new TFile(inputpTDistr);
   if(!fin->IsOpen()){
      Printf("Check input file name");
      return;
   }
   
   TGraphAsymmErrors *gpT = (TGraphAsymmErrors*) fin->Get(graphname);
   
   if(!gpT){
      Printf("Check input graph name");
      return;
   }

   TCanvas *cpT = new TCanvas("cpT", "pT distribution");
   gpT->Draw();
   TF1* func = new TF1("func", "[0]*TMath::Power([1], 2)*x*TMath::Exp(-[1]*x)", 0.150, 15);
   func->SetParameter(1, 0.696);
   
   gpT->Fit(func);
   
   Int_t nbins = gpT->GetN();
   
   Double_t xvalues[nbins+1];
   
   for(Int_t i=0;i<nbins;i++){
      Double_t y;
      gpT->GetPoint(i, xvalues[i],y);
      Double_t dx = gpT->GetErrorXlow(i);
//      Printf("point %d x = %f - %f", i, xvalues[i], dx);
      xvalues[i] -= dx;
      Printf("point %d x = %f", i, xvalues[i]);
      if(i==nbins-1) {
      	 dx = gpT->GetErrorXhigh(i);
      	 xvalues[i+1] = xvalues[i]+ dx;
      	 Printf("point %d x = %f", i+1, xvalues[i+1]);
      }
   }
   
   TH1F* hpT = new TH1F("hpT", ";p_{T};dN/dp_{T}", nbins, xvalues);
   
   hpT->SetLineColor(2);
   
   Printf("N = %d ", nbins);
   for(Int_t i=0;i<nbins;i++){
      Double_t x,y;
      gpT->GetPoint(i,x,y);
      Int_t bin = i+1;//hpT->FindBin(x);
      hpT->SetBinContent(bin, y);
      hpT->SetBinError(bin, gpT->GetErrorY(i));
   }
   
   cpT->cd();
   hpT->Draw("sames");
   TFile *fout = new TFile(Form("%s.root", hpT->GetName()), "recreate");
   
   hpT->Write();
   
}
