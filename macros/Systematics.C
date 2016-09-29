#include </data/Work/MyCodeJetMass/macros/plotJetMasspPb.C>

//#include </data/macros/LoadALICEFigures.C>
//#include <TF1.h>
//#include <TH1F.h>
//#include <TCanvas.h>
//#include <TString.h>
//#include <TFile.h>

//Max pol1 (60.000000 GeV/c) = 8.750615, Min pol1 (120.000000 GeV/c) = 7.811229, Compare to pol0 param 0 = 8.359860
void SystematicTrigTurnOnFactor(TString input = "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output794-795/TurnOnCurves.root"){

   const Int_t nmethods = 3; //
   TString methodname[nmethods] = {"ConstSub", "Area-MRaw", "Area-MDeriv"};
   const Int_t ntags = 2;
   TString tagsname[ntags] = {"All", "Tagged"};//, "TaggedMatch"
   PlotUtilEmcalJetMass::JetSelection jetSel[ntags] = {PlotUtilEmcalJetMass::kJetAll, PlotUtilEmcalJetMass::kJetTagged};
   Int_t minbias[nmethods]    = {0,1,2};
   Int_t triggersJ1[nmethods] = {3,4,5};
   TString lists[ntag*nmethods] = {"ConstSub_TC", "_TCRaw", "_TCDeriv", "ConstSub_TCConstJ1", "_TCRawJ1", "_TCDerivJ1"};

   const Int_t nvariation = 3;
   Double_t factors[nvariation][nmethods][ntags];
   
   TFile *fin = new TFile(input);
   if(!fin->IsOpen()){
      Printf("File %s not found", input.Data());
      return;
   }
   
   //determination of the factor
   TString hname = "hpTRJet";
   TString triggername = "J1";
   Double_t range[2] = {60.,120.};
   TF1 *fpol0 = new TF1("fpol0", "pol0",range[0], range[1]);
   TF1 *fpol1 = new TF1("fpol1", "pol1",range[0], range[1]);
   TGraph *gTrFactor[ntags][nvariation];
   TCanvas *cFactors = new TCanvas("cFactors", "Factors");

   for(Int_t j = 0; j< ntags; j++){
      for(Int_t n = 0; n< nvariation; n++){
      	 gTrFactor[j][n] = new TGraph(nmethods);
      	 gTrFactor[j][n] ->SetName(Form("gFactor%s", tagsname[j].Data()));
      	 gTrFactor[j][n] ->SetTitle(Form("Trigger factor (Variation %d) %s Jets; Bkg Method; Factor", n, tagsname[j].Data()));
      	 gTrFactor[j][n] ->SetMarkerColor(colors[n]);
      	 gTrFactor[j][n] ->SetMarkerStyle(20+j);
 
      }
   }
   for(Int_t j = 0; j< ntags; j++){
      for(Int_t i = 0; i < nmethods; i++){
      	 
      	 
      	 TString thisname = Form("%s%s%s%d", hname.Data(), tagsname[j].Data(), triggername.Data(), i);
      	 TH1D *hturnon = (TH1D*)fin->Get(thisname);
      	 if(!hturnon){
      	    Printf("%s not found", thisname.Data());
      	    continue;
      	 }
      	 
      	 //fit
      	 hturnon->Fit(fpol0,"0RL+");
      	 hturnon->Fit(fpol1,"0RL+");
      	 
      	 factors[0][i][j] = fpol1->Eval(range[0]);
      	 factors[1][i][j] = fpol0->GetParameter(0);
      	 factors[2][i][j] = fpol1->Eval(range[1]);
      	 
      	 gTrFactor[j][0]->SetPoint(i, i, factors[0][i][j]);
      	 gTrFactor[j][1]->SetPoint(i, i, factors[1][i][j]);
      	 gTrFactor[j][2]->SetPoint(i, i, factors[2][i][j]);
      	 
      	 Printf("Max pol1 (%f GeV/c) = %f, Min pol1 (%f GeV/c) = %f, Compare to pol0 param 0 = %f", range[0], fpol1->Eval(range[0]), range[1], fpol1->Eval(range[1]), fpol0->GetParameter(0));
      }
      cFactors->cd();
      if(j==0) {
      	 gTrFactor[j][0]->Draw("AP");
      	  gTrFactor[j][0]->GetYaxis()->SetRangeUser(7,10);
      }
      else gTrFactor[j][0]->Draw("P");
      gTrFactor[j][1]->Draw("P");
      gTrFactor[j][2]->Draw("P");
   }
   
   
   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   Bool_t bEmb = kFALSE;
   
   PlotUtilEmcalJetMass *util=0x0;
   InitUtilFiles(util, ntag, lists, nfi, pathfile, R, pTjet[0], pTjet[npTjbins], type, centBin, bEmb);
   
   TString suffix = "TrFactorVariation";
   
   TFile *fout = new TFile(Form("MassOutput%s.root",suffix.Data()), "recreate");
   /*
   TGraphErrors* hMassDiffToNominal[nmethods][ntags][nvariation-1];

   TCanvas *cjmptj[ntags][nmethods];
   
   for(Int_t j = 0; j<ntags; j++) {
      for(Int_t i = 0; i<nmethods; i++) {
      	 cjmptj[i][j] = new TCanvas(Form("cjmptj%s%d",tagsname[j].Data(), i), Form("%s Jet mass %s in bins of jet pT",tagsname[j].Data(), methodname[i].Data()), 1000,1000);
      	 cjmptj[i][j]->Divide(3,2);
      	 Printf("Defined %s", cjmptj[i][j]->GetName());
      	 
      }
   }
   */
   for(Int_t iptj=0 ; iptj<npTjbins; iptj++){
      
      TPaveText *pvpT = new TPaveText(0.6, 0.7, 0.8, 0.9, "NDC");
      pvpT->SetFillStyle(0);
      pvpT->SetBorderSize(0);
      pvpT->AddText(Form("p_{T,jet} %.0f - %.0f", pTjet[iptj],pTjet[iptj+1]));
      
      //pT leading track integrated plots
      util->SetMinLeadTrackPt(0);
      util->SetJetPtRange(pTjet[iptj], pTjet[iptj+1]);
      
      TH1D *hMass[nmethods][ntags];
      TGraphErrors* hMeanMass[nmethods][ntags];
      
      
      for(Int_t j = 0; j<ntags; j++) {
      	 
      	 for(Int_t i = 0; i<nmethods; i++) {
      	    hMass[i][j]=0;
      	    hMeanMass[i][j]=0;
      	    Printf("Filling pT lead integrated pTj %.0f - %.0f, Jet %s, %s", pTjet[iptj], pTjet[iptj+1], tagsname[j].Data(), methodname[i].Data());
      	    
      	    //min bias
      	    TH1D* hMassMB = util->GetJetMassDistribution(jetSel[j],minbias[i]);
      	    
      	    //if(i==0) printf("Filling pT lead integrated pTj %.0f - %.0f: min bias", pTjet[iptj], pTjet[iptj+1]);
      	    hMassMB->SetName(Form("hMassMB%s%dpTj%d",tagsname[j].Data(), i, iptj));
      	    Printf("Number of entries mb %.0f , %f, %f", hMassMB->GetEntries(), hMassMB->GetBinLowEdge(1), hMassMB->GetBinLowEdge(hMassMB->GetNbinsX()+1));
      	    
      	    if(pTjet[iptj]>59){
      	       //J1
      	       for(Int_t sys=0; sys<nvariation ; sys++){
      	       	  util->SetScaleFactor(factors[sys][i][j]);
      	       	  hMass[i][j] = util->GetJetMassDistribution(jetSel[j],triggersJ1[i]);
      	       	  hMass[i][j]->SetName(Form("hMass%s%dpTj%dSys%d",tagsname[j].Data(), i, iptj,sys));
      	       	  hMass[i][j]->SetTitle(Form("Mass %s %dpTj%dSys%d",tagsname[j].Data(), i, iptj,sys));
      	       	  Printf("Integral %f ", hMass[i][j]->Integral());
      	       	  hMass[i][j]->Add(hMassMB);
      	       	  
      	       	  Printf("Number of entries mb+j1 %.0f", hMass[i][j]->GetEntries());   
      	       	  hMass[i][j]->SetLineColor(colors[sys]);
      	       	  hMass[i][j]->SetMarkerColor(hMass[i][j]->GetLineColor());
      	       	  hMass[i][j]->SetMarkerStyle(20);
      	       	  fout->cd();
      	       	  hMass[i][j]->Write();
      	       }
      	    } else hMass[i][j] = (TH1D*)hMassMB->Clone(Form("hMass%s%dpTj%d",tagsname[j].Data(), i, iptj));
      	    if(!hMass[i][j]) continue;
      	    //hMAllpTleadint[i]->Scale(1./norm);
      	    
      	    //hMAllpTleadint[i]->GetYaxis()->SetRangeUser(0.,0.4);
      	    //hMass[i][j]->GetXaxis()->SetRangeUser(-5,20.);
      	    
      	    //cjmptj[i][j]->cd(iptj+1);
      	    //hMass[i][j]->Draw("samesP");
      	    
      	 }
      }
   }
}

void CorrectionFluctuations(TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/UnfoldingMatrix.root"){
   
   TFile *fin = new TFile(Form("%s", inputFile.Data()));
   if(!fin->IsOpen()){
      Printf("Input file for correction not found");
      return;
   }

   TString hFactorname = "hMassPartToDetc";
   TString append = "Less40";
   
   TF1* funcpol = new TF1("funcpol", "[0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x", 0, 50.);
   
   TCanvas *cFits = new TCanvas("cFits", "Fits", 1000,800);
   cFits->Divide(3,2);
   for(Int_t iptj = 0; iptj < npTjbins ; iptj++){
      TH1D *hFactor = (TH1D*) fin->Get(Form("%s%s_pTj%d", hFactorname.Data(), append.Data(), iptj));
      
      funcpol->SetRange(hFactor->GetBinLowEdge(1), hFactor->GetBinLowEdge(hFactor->GetNbinsX()));
      
      //Fit (bad)
      hFactor->Fit(funcpol, "R+0");
      cFits->cd(iptj+1);
      hFactor->Draw();
      hFactor->GetFunction("funcpol")->Draw("sames");
   
      //Uncertainty Band? cosa serve? l'errore statistico gia' ne tiene conto
   }
   

}
