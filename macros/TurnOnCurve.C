#include "/data/Work/MyCodeJetMass/macros/plotJetMassLeadTrkDataUtil.C"

const Int_t nfi = 2;//1;
TString pathfile[nfi] = //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output771-772/AnalysisResults.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output775-774/AnalysisResults.root"}; //MB, EJE
//{"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output792-793/AnalysisResults.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output794-795/AnalysisResults.root"}; //MB, EJE
{"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output810-811/AnalysisResults.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output806-807/AnalysisResults.root"}; //MB, EJE with rho and rhom
//{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train833/outputpTHardBins/mergeRuns/merge1to9/AnalysisResults.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train833/outputpTHardBins/mergeRuns/merge1to9/AnalysisResults.root"};
//MC LHC13b4_plus trigger and MB
//{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root"};
//MB only (modify also plotJetMasspPb)
//{"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output834-835/AnalysisResults.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output836-837/AnalysisResults.root"}; //MB, EJE with rho and rhom and centrality bins ZNA
//{"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output838-839/AnalysisResults.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output840-841/AnalysisResults.root"}; //MB, EJE with rho and rhom and centrality bins V0A
//{"/data/Work/jets/JetMass/ppJetMassAnalysis/merge333-336/Flat/AnalysisResults.root" };
//2010 pp MB (uncomment ntag = 1, comment ntag = 9)
//const Int_t ntag = 1;
//TString turnontags[ntag] = {"_TCRaw"};


const Int_t ntag = 9; // 1; TString turnontags[ntag] = {"_TCMerged"};
TString turnontags[ntag] = {"ConstSub_TC", "_TCRaw", "_TCDeriv", "ConstSub_TCConstJ1", "_TCRawJ1", "_TCDerivJ1", "ConstSub_TCConstJ2", "_TCRawJ2", "_TCDerivJ2"};

Int_t availableto[ntag];

//rebin
const Int_t nentries = 25;
//Double_t newarray[nentries] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120};
Double_t newarray[nentries] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120};

void DefNewBins(TH1D* hTrig, Double_t* rebinpt, Int_t newnptbins);
Bool_t InitUtilFiles(PlotUtilEmcalJetMass *&util, Int_t nt, TString *tags, Int_t nf, TString *strf, Double_t R, Double_t ptmin, Double_t ptmax, TString type, Int_t centBin, Bool_t bEmb);

Bool_t InitUtilFiles(PlotUtilEmcalJetMass *&util, Int_t nt, TString *tags, Int_t nf, TString *strf, Double_t R, Double_t ptmin, Double_t ptmax, TString type, Int_t centBin, Bool_t bEmb){
//  gStyle->SetOptStat(0);
//  gStyle->SetOptTitle(0);
  //gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
  //gROOT->LoadMacro("$gitJetMass/utils/style.C");
  SetStyle(0);
 
  Double_t centMin[4] = {0.,10.,30.,50.};
  Double_t centMax[4] = {10.,30.,50.,80.};

  TString outName = Form("JetMassCent%dR%03dPtMin%.0fPtMax%.0fType%s.root",centBin,(Int_t)(R*100.),ptmin,ptmax,type.Data());

  //gROOT->LoadMacro("/data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.cxx+");
  util = new PlotUtilEmcalJetMass();
  for(Int_t ifi=0; ifi<nf;ifi++){
     util->SetInputFileName(strf[ifi]);
     util->SetJetRadius(R);
     util->SetJetType(type);
     //util->SetCentBin(centBin);
     util->SetJetPtRange(ptmin,ptmax);
     util->LoadFile();
     
     for(Int_t j=0;j<nt;j++){
      	
     	util->SetTag(tags[j]); //"ConstSub_TC"
     	util->SetConstTag("");
     	if(util->LoadList()) {
     	   availableto[j] = 1;
     	   Printf("Loaded lists %d (%s) --> %d", j, tags[j].Data(), availableto[j]);  
     	}
      	else continue;
 
     }
  }
  Printf("Load done.");
  return kTRUE;
}

void TurnOnCurve(Double_t R = 0.4, Double_t ptmin = 0., Double_t ptmax = 120., TString type = "Charged", Int_t centBin = 0, Bool_t bEmb = kFALSE){
   PlotUtilEmcalJetMass *util=0x0;
   InitUtilFiles(util, ntag, turnontags, nfi, pathfile, R, ptmin, ptmax, type, centBin, bEmb);
   util->SetMinLeadTrackPt(0);
   const Int_t nmethods = 3;
   Int_t minbias[nmethods]    = {0,1,2};
   Int_t triggersJ1[nmethods] = {3,4,5};
   Int_t triggersJ2[nmethods] = {6,7,8};
   
   TCanvas *ctrTurn[nmethods];
   for (Int_t i=0;i<nmethods;i++) ctrTurn[i] = new TCanvas(Form("ctrTurn%d", i), Form("Turn on curves J1 and J2 - %s", turnontags[minbias[i]].Data()));
   const Int_t ntrig = 2;
   TString trigname[ntrig] = {"J1", "J2"} ; //for pPb 2013: J1 high threshold, J2 low threshold 
   TCanvas *ctrTurn2[ntrig];
   for (Int_t i=0;i<ntrig;i++)   ctrTurn2[i] = new TCanvas(Form("ctrTurnTr%d", i), Form("Turn on curves %s", trigname[i].Data()));
   TCanvas *cpTExample = new TCanvas("cpTExample", "pT distribution", 700, 700);
   
   TH1D *hpTJetAll[ntrig+1][nmethods];
   TH1D *hpTJetTagged[ntrig+1][nmethods];
   TH1D *hpTTaggedMatch[ntrig+1][nmethods];
   TH1D *hpTRJetAll[ntrig][nmethods];
   TH1D *hpTRJetTagged[ntrig][nmethods];
   TH1D *hpTRJetTaggedMatch[ntrig][nmethods];

   Int_t ktag = 0;
   const Int_t nnewbins = 100;
   Double_t newbinlims[nnewbins];

   for(Int_t i=0; i < ntrig+1; i++){ // loop on triggers
      TString thistrigname = (i==0) ? "MB" : trigname[i-1];
      for(Int_t j=0; j < nmethods; j++){ //loop on sub method
      	 Printf("Reading from list %d (%s)", ktag, turnontags[ktag].Data());
      	 hpTJetAll[i][j] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetAll,ktag);
      	 hpTJetAll[i][j]->SetName(Form("hJpTAll%s%d", thistrigname.Data(), j));
      	 hpTJetAll[i][j]->SetLineWidth(2);    
      	 hpTJetAll[i][j]->SetLineColor(colors[i]);
      	 hpTJetAll[i][j]->SetMarkerColor(colors[i]);
      	 hpTJetAll[i][j]->SetMarkerStyle(20+j);
      	 hpTJetTagged[i][j] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTagged,ktag);
      	 hpTJetTagged[i][j]->SetName(Form("hJpTTagged%s%d", thistrigname.Data(), j));      	 
      	 hpTJetTagged[i][j]->SetLineColor(colors[i]);
      	 hpTJetTagged[i][j]->SetLineWidth(2);
      	 hpTJetTagged[i][j]->SetMarkerColor(colors[i]);
      	 hpTJetTagged[i][j]->SetMarkerStyle(20+j);
      	 
      	 cpTExample->cd();
      	 gPad->SetLogy();
      	 hpTJetTagged[i][j]->Draw("sames");
      	 
      	 hpTTaggedMatch[i][j] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,ktag);
      	 hpTTaggedMatch[i][j]->SetName(Form("hJpTTaggedMatch%s%d", thistrigname.Data(), j));
      	 hpTTaggedMatch[i][j]->SetLineColor(colors[i]);
      	 hpTTaggedMatch[i][j]->SetLineWidth(2);
      	 hpTTaggedMatch[i][j]->SetMarkerColor(colors[i]);
      	 hpTTaggedMatch[i][j]->SetMarkerStyle(20+j);
      	 
      	 hpTJetAll[i][j]      = (TH1D*)hpTJetAll[i][j]->Rebin(nentries-1, Form("%srb", hpTJetAll[i][j]->GetName()), newarray);
      	 hpTJetTagged[i][j]   = (TH1D*)hpTJetTagged[i][j]->Rebin(nentries-1, Form("%srb", hpTJetTagged[i][j]->GetName()), newarray);
      	 hpTTaggedMatch[i][j] = (TH1D*)hpTTaggedMatch[i][j]->Rebin(nentries-1, Form("%srb", hpTTaggedMatch[i][j]->GetName()), newarray);
      	 Printf("Histo name %s ", hpTJetAll[i][j]->GetName());
      	 
      	 if(i>0){
      	    hpTRJetAll[i-1][j] = (TH1D*)hpTJetAll[i][j]->Clone(Form("hpTRJetAll%s%d", thistrigname.Data(), j));
      	    hpTRJetAll[i-1][j]->Divide(hpTJetAll[0][j]);
      	    
      	    hpTRJetTagged[i-1][j] = (TH1D*)hpTJetTagged[i][j]->Clone(Form("hpTRJetTagged%s%d", thistrigname.Data(), j));
      	    hpTRJetTagged[i-1][j]->Divide(hpTJetTagged[0][j]);
      	    
      	    hpTRJetTaggedMatch[i-1][j] = (TH1D*)hpTTaggedMatch[i][j]->Clone(Form("hpTRJetTaggedMatch%s%d", thistrigname.Data(), j));
      	    hpTRJetTaggedMatch[i-1][j]->Divide(hpTTaggedMatch[0][j]);
      	    Printf("Ratio Histo name %s ", hpTRJetAll[i-1][j]->GetName());
      	 }
      	 ktag++;
      	 
      	 
      }
      
   
   }
   
   cpTExample->SaveAs(Form("%s.pdf", cpTExample->GetName()));
   
   //loops for drawing
      
   TFile *fout = new TFile("TurnOnCurves.root", "recreate");

   for(Int_t i=0;i<ntrig;i++){
      for(Int_t j=0; j < nmethods; j++){
      	 ctrTurn2[i]->cd(); //canvas with the different triggers turn on for tagged jets
      	 if(j==0) hpTRJetTagged[i][j]->Draw("P");
      	 else hpTRJetTagged[i][j]->Draw("samesP");
      	 ctrTurn2[i]->SaveAs(Form("%s.pdf",ctrTurn2[i]->GetName()));
      	 
      	 ctrTurn[j]->cd();
      	 if(i==0) hpTRJetTagged[i][j]->Draw("P");
      	 else hpTRJetTagged[i][j]->Draw("samesP");
      	 ctrTurn[j]->SaveAs(Form("%s.pdf",ctrTurn[j]->GetName()));
      	 
      	 fout->cd();
      	 hpTRJetAll[i][j]->Write();
      	 hpTRJetTagged[i][j]->Write();
      	 hpTRJetTaggedMatch[i][j]->Write();
      	 
      	 
      }
      
      
      
   }
   
   
     	 
      	 //
      	 //gPad->SetLogy();
      	 //if(j == 0) hpTJetTagged[i][j]->Draw();
      	 //else hpTJetTagged[i][j]->Draw("sames");
      	 //
      	 //gPad->SetLogy();
   
   

}

void DefNewBins(TH1D* hTrig, Double_t* rebinpt, Int_t newnptbins){
   const Int_t nsteps=4;
   //Double_t steps[nsteps]={28,35,50,70}; //(in GeV/c)
   //Int_t deltapt[nsteps]={2,6,15,45};//rebinning at each step 
   Double_t steps[nsteps]={20,50,70,90}; //(in GeV/c)
   Int_t deltapt[nsteps]={2,6,15,45};//rebinning at each step
   
      Int_t ii=0,jj=0;
      for(Int_t k=0;k<newnptbins;k++){
      	 if(hTrig->GetBinLowEdge(k)<0) {continue;}
      	 if(hTrig->GetBinLowEdge(k)<steps[0]) rebinpt[k]=hTrig->GetBinLowEdge(k+1);
      	
      	 for(Int_t s=0;s<nsteps-1;s++){
      	    Printf("LOOP on steps corresponding to bin edge %f",hTrig->GetBinLowEdge(k));
      	    if(hTrig->GetBinLowEdge(k)>= steps[s] && hTrig->GetBinLowEdge(k)<steps[s+1]){
      	       Printf("S = %d , ii= %d, deltapt = %d, k = %d",s,ii,deltapt[s],k);
      	       rebinpt[k]=rebinpt[k-1]+deltapt[s];
      	       ii++;
      	       break;
      	    }
      	 }
      	 
      	 if(hTrig->GetBinLowEdge(k)>=steps[nsteps-1]){
      	    Printf("S = %d , jj= %d, deltapt= %d, k = %d",nsteps-1,jj,deltapt[nsteps-1],k);
      	    rebinpt[k]=rebinpt[k-1]+deltapt[nsteps-1];;
      	    jj++;
      	 }
      	 Printf("%d => %f",k,rebinpt[k]);
      }
      Printf("Array for rebin, %d entries. New numb bins %d ",newnptbins,newnptbins-1);
   
   
}

