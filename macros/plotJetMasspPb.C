#include "/data/Work/MyCodeJetMass/macros/TurnOnCurve.C"
#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"

const Int_t npTjbins = 5;
Double_t pTjet[npTjbins+1] = {20.,40.,60.,80.,100.,120.};
//const Int_t npTjbins = 4;
//Double_t pTjet[npTjbins+1] = {40.,60.,80.,100.,140.};
Double_t GetFactorMBJetTrigger(TString type, TString trig, Int_t methodnumb);
TH1D* GetTurnOnTrigger(TString type, TString trig, Int_t methodnumb);
//void SaveCv(TCanvas *c, TString suff = "");
Double_t GetQuantile(TH1 *h1, Double_t prob=0.5);
void Superimpose(Int_t nfiles, TString paths[], TString legs[], Int_t methodIndex, TString tag);
void SumHistogramstoMassOutput(Int_t nfiles, TString filenames[], TString meaningfulflags[], Int_t MBEJETrigComb = 1 /*0 = only MB, 1 = MB for pT < 80, EJE pT>80, 2 = sum MB and EJE for pT> 80, 3 = include J1 */, Int_t centBin = -1);
void MassVsCentrality(TString tag = "_TCDeriv", TString path = "");
void plotJetMasspPbSimple(Int_t ntag = 1, TString turnontags[] = 0, Bool_t bEmb = kFALSE, TString pathfile = "/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/merge/weights/analysis/AnalysisResults.root" /*Embedding PYTHIA thermal*/);
void RhopPbvsPYTHIA(TString dataPath, TString MCpath, TString trkname[] = 0, Int_t listnametype[] = 0);
void RhopPbvsPYTHIA(Int_t ninput, TString inputF[], TString trkname[] = 0, Int_t listnametype[] = 0, Bool_t cent = kTRUE);
void RhopPbvsPYTHIApTLeading(TString dataPath, TString MCpath, TString tag = "Deriv");

//-------------------------------------------------------------------------------------------

void plotJetMasspPb(TString suffix="", Int_t MBEJETrigComb = 0 /*0 = only MB, 1 = MB for pT < 80, EJE pT>80, 2 = sum MB and EJE for pT> 80, 3 = include J1 */, Bool_t bEmb = kFALSE){
   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   
   if (bEmb) type = "ThrmTracksEmb";
   
   Printf("READING...\n MB %s\n EJE %s", pathfile[0].Data(), pathfile[1].Data());
   
   PlotUtilEmcalJetMass *util=0x0;
   InitUtilFiles(util, ntag, turnontags, nfi, pathfile, R, pTjet[0], pTjet[npTjbins], type, centBin, bEmb);
   
   //lists, as in TurnOnCurve where InitUtilFiles is taken from
   const Int_t nmethods = 3;
   
   Int_t minbias[nmethods]    = {0,1,2};
   Int_t triggersJ1[nmethods] = {3,4,5};
   Int_t triggersJ2[nmethods] = {6,7,8};
   
   const Int_t ntrig = 2;
   TString trigname[ntrig] = {"J1", "J2"} ; //for pPb 2013: J1 high threshold, J2 low threshold 
  //TString trigname[ntrig] = {"J2", "J1"} ; //MC LHC13b4_plus swapped!!!
   //Printf("Using MC trigger names!!!!!!!!!!!! Thresholds swapped wrt data");
   
   //suffix="MB";
   //suffix=trigname[0];
   //suffix=trigname[1];
   const Int_t nLeadTrkBins = 6;
   Double_t leadTrkPt[nLeadTrkBins] = {0.,2.,6.,10.,15.,19.};
   
   //Double_t factorJ1 = GetFactorMBJetTrigger("All", trigname[0], 1);
   
   //Double_t factorJ2 = GetFactorMBJetTrigger("All", trigname[1], 1);
   //Printf("Factor J1 = %f, Factor J2 = %f", factorJ1, factorJ2);

   // 2D distributions, M vs pT for MB+triggers according to MBEJETrigComb:  = 0 only MB, 1 = MB for pT < 80, EJE pT>80, 2 = sum MB and EJE for pT> 80, 3 = include J1
   
   //pT bins
   
   for(Int_t iptj=0 ; iptj<npTjbins; iptj++){
      
      util->SetJetPtRange(pTjet[iptj], pTjet[iptj+1]);
      TH1D *hMAll[nmethods][nLeadTrkBins];
      TH1D *hMTagged[nmethods][nLeadTrkBins];
      TH1D *hMTaggedMatch[nmethods][nLeadTrkBins];
      
      TGraphErrors *grMeanPtLeadTrAll[nmethods];
      TGraphErrors *grMeanPtLeadTrTagged[nmethods];
      TGraphErrors *grMeanPtLeadTrTaggedMatch[nmethods];
      
      TLegend *leg[nmethods];

      for(Int_t i = 0; i<nmethods; i++) {
      	 grMeanPtLeadTrAll[i] = new TGraphErrors();
      	 grMeanPtLeadTrAll[i]->SetTitle(";#it{p}_{T,lead track};#LT#it{M}_{jet}#GT");
      	 grMeanPtLeadTrAll[i]->SetMarkerStyle(20);
      	 grMeanPtLeadTrAll[i]->SetMarkerColor(colors[i]);
      	 
      	 grMeanPtLeadTrTagged[i] = new TGraphErrors();
      	 grMeanPtLeadTrTagged[i]->SetTitle(";#it{p}_{T,lead track};#LT#it{M}_{jet}#GT");      	 grMeanPtLeadTrTagged[i]->SetMarkerStyle(24);      	 grMeanPtLeadTrTagged[i]->SetMarkerColor(colors[i]);
      	 
      	 grMeanPtLeadTrTaggedMatch[i] = new TGraphErrors();
      	 grMeanPtLeadTrTaggedMatch[i]->SetTitle(";#it{p}_{T,lead track};#LT#it{M}_{jet}#GT");
      	 grMeanPtLeadTrTaggedMatch[i]->SetMarkerStyle(27);
      	 grMeanPtLeadTrTaggedMatch[i]->SetMarkerColor(colors[i]);
      	 
      	 
      	 leg[i] = CreateLegend(0.25,0.6,0.52,0.94,strLeg[i].Data());
      	 leg[i]->SetTextSize(leg[i]->GetTextSize()*0.8);
      }
      TCanvas *cmassptjbin = new TCanvas(Form("cmassptjbin%d", iptj), Form("Mass for pTj %.1f - %.1f", pTjet[iptj], pTjet[iptj+1]), 800, 800);
      cmassptjbin->Divide(3,2);
      //cmassptjbin->Divide(2,2);
      
      Double_t norm = 0.;
      
      for(Int_t j = 0; j<nLeadTrkBins; j++) {
      	 util->SetMinLeadTrackPt(leadTrkPt[j]);
      	 for(Int_t i = 0; i<nmethods; i++) {
      	    //for(Int_t i = 0; i<1; i++) { //when 1 method only, MB only
      	    //min bias
      	    hMAll[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,minbias[i]);
      	    if(!hMAll[i][j]) {
      	       Printf("Error, mass histogram is null");
      	       continue;  
      	    }
      	    hMTagged[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,minbias[i]);
      	    if(!hMTagged[i][j]) {
      	       Printf("Error, mass tagged histogram is null");
      	       continue;  
      	    }
      	    hMTaggedMatch[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,minbias[i]);
      	    if(!hMTaggedMatch[i][j]) {
      	       Printf("Error, mass tagged match histogram is null");
      	       continue;  
      	    }
      	    if(i==0 && j==0) printf("pTj %.0f - %.0f: min bias", pTjet[iptj], pTjet[iptj+1]);
      	    hMAll[i][j]->SetName(Form("%spTj%d", hMAll[i][j]->GetName(), iptj));
      	    hMTagged[i][j]->SetName(Form("%spTj%d", hMTagged[i][j]->GetName(), iptj));
      	    hMTaggedMatch[i][j]->SetName(Form("%spTj%d", hMTaggedMatch[i][j]->GetName(), iptj));
      	    if(MBEJETrigComb){
      	       if(pTjet[iptj]>59){
      	       	  //J1
      	       	  if(MBEJETrigComb == 2){
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMAll[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggersJ1[i]));
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMTagged[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggersJ1[i]));
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMTaggedMatch[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggersJ1[i]));
      	       	  }
      	       	  if(MBEJETrigComb == 1){
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMAll[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggersJ1[i]);
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMTagged[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggersJ1[i]);
      	       	     //util->SetScaleFactor(factorJ1);
      	       	     hMTaggedMatch[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggersJ1[i]);
      	       	     
      	       	     
      	       	  }
      	       	  if(i==0 && j==0) printf(" + J1 (ht)");
      	       }
      	       
      	       if(pTjet[iptj]>39){
      	       	  //J2
      	       	  if(MBEJETrigComb == 3){
      	       	     //util->SetScaleFactor(factorJ2);
      	       	     hMAll[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggersJ2[i]));
      	       	     //util->SetScaleFactor(factorJ2);
      	       	     hMTagged[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggersJ2[i]));
      	       	     //util->SetScaleFactor(factorJ2);
      	       	     hMTaggedMatch[i][j]->Add( util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggersJ2[i]));
      	       	  }
      	       	  if(i==0 && j==0) printf(" + J2 (lt)");
      	       }
      	    }
      	    if(i==0 && j==0) printf("\n");
      	    grMeanPtLeadTrAll[i]->SetPoint(grMeanPtLeadTrAll[i]->GetN(),leadTrkPt[j],hMAll[i][j]->GetMean());
      	    grMeanPtLeadTrTagged[i]->SetPoint(grMeanPtLeadTrTagged[i]->GetN(),leadTrkPt[j],hMTagged[i][j]->GetMean());
      	    grMeanPtLeadTrTaggedMatch[i]->SetPoint(grMeanPtLeadTrTaggedMatch[i]->GetN(),leadTrkPt[j],hMTaggedMatch[i][j]->GetMean());
      	    
      	    grMeanPtLeadTrAll[i]->SetPointError(grMeanPtLeadTrAll[i]->GetN()-1,0.,hMAll[i][j]->GetMeanError());
      	    grMeanPtLeadTrTagged[i]->SetPointError(grMeanPtLeadTrTagged[i]->GetN()-1,0.,hMTagged[i][j]->GetMeanError());
      	    grMeanPtLeadTrTaggedMatch[i]->SetPointError(grMeanPtLeadTrTaggedMatch[i]->GetN()-1,0.,hMTaggedMatch[i][j]->GetMeanError());
      	    
      	    leg[i]->AddEntry(hMTagged[i][j],Form("#it{p}_{T,lead trk}=%.0f #LT#it{M}_{jet}#GT=%.1f RMS=%.1f",leadTrkPt[j],hMTagged[i][j]->GetMean(),hMTagged[i][j]->GetRMS()),"p");
      	    
      	    if(j==0) norm = hMAll[i][j]->Integral();
      	    if(norm>0.) {
      	       //Printf("norm: %f",norm);
      	       hMAll[i][j]->Scale(1./norm);
      	       
      	       hMAll[i][j]->SetLineColor(colors[j]);
      	       hMAll[i][j]->SetMarkerColor(hMAll[i][j]->GetLineColor());
      	       hMAll[i][j]->SetMarkerStyle(20);
      	       hMAll[i][j]->GetYaxis()->SetRangeUser(0.,0.4);
      	       hMAll[i][j]->GetXaxis()->SetRangeUser(-5,15.);
      	       hMTagged[i][j]->Scale(1./norm);
      	       hMTagged[i][j]->SetLineColor(colors[j]);
      	       hMTagged[i][j]->SetMarkerColor(hMTagged[i][j]->GetLineColor());
      	       hMTagged[i][j]->SetMarkerStyle(24);
      	       hMTaggedMatch[i][j]->Scale(1./norm);
      	       hMTaggedMatch[i][j]->SetLineColor(colors[j]);
      	       hMTaggedMatch[i][j]->SetMarkerColor(hMTaggedMatch[i][j]->GetLineColor());
      	       hMTaggedMatch[i][j]->SetMarkerStyle(27);
      	       
      	    }
      	    cmassptjbin->cd(i+1);
      	    hMAll[i][j]->Draw("Psames");
      	    hMTagged[i][j]->Draw("Psames");
      	    hMTaggedMatch[i][j]->Draw("Psames");
      	    if(j == nLeadTrkBins-1) {
      	       leg[i]->Draw();
      	    }
      	    
      	    
      	 }
      }
      
      
      for(Int_t i=0; i<nmethods;i++){
      	 
      	 //draw mean mass distribution
      	 
      	 cmassptjbin->cd(4);
      	 //if(i==0) grMeanPtLeadTrAll[i]->Draw("AP");
      	 //else  grMeanPtLeadTrAll[i]->Draw("P");
      	 grMeanPtLeadTrAll[i]->GetYaxis()->SetRangeUser(1.8,12.5);
      	 grMeanPtLeadTrTagged[i]->GetYaxis()->SetRangeUser(1.8,12.5);
      	 grMeanPtLeadTrTaggedMatch[i]->GetYaxis()->SetRangeUser(1.8,12.5);
      	 if(i==0) grMeanPtLeadTrTagged[i]->Draw("ALP");
      	 else  grMeanPtLeadTrTagged[i]->Draw("LP");
      	 //grMeanPtLeadTrTaggedMatch[i]->Draw("P");
      	      	 
       }
      SaveCv(cmassptjbin, suffix.Data());
      
   }

   
}

//-------------------------------------------------------------------------------------------

void plotMassCompareTriggers(Int_t MBEJETrigComb = 1 /*0 = only MB, 1 = MB for pT < 80, EJE pT>80, 2 = sum MB and EJE for pT> 80, 3 = include J1 */, Int_t var = 1 /*1 = M, 2 = pT*/, Int_t centBin = -1){
   
   Double_t R = 0.4;
   TString type = "Charged";

   Bool_t bEmb = kFALSE;
   
   Printf("READING...\n MB %s\n EJE %s", pathfile[0].Data(), pathfile[1].Data());
   PlotUtilEmcalJetMass *util=0x0;
   InitUtilFiles(util, ntag, turnontags, nfi, pathfile, R, pTjet[0], pTjet[npTjbins], type, centBin, bEmb);
   
   
   //lists, as in TurnOnCurve where InitUtilFiles is taken from
   const Int_t nmethodsAll = 3;
   Int_t minbias[nmethodsAll]    = {0,1,2}; //Const sub, Raw, Deriv
   const Int_t ntrigc = 2;
   Int_t triggers[ntrigc][nmethodsAll];
   TString methnames[nmethodsAll] = {"Const sub", "No sub", "Deriv sub"};
   Int_t high = 0, low = 1;    
   TString trigname[ntrigc] = {"J1", "J2"} ; //for pPb 2013: J1 high threshold, J2 low threshold
   Double_t thresh[ntrigc]; thresh[high] = 79; thresh[low] = 59;
   triggers[high][0] = 3; triggers[high][1] = 4; triggers[high][2] = 5;
   triggers[low][0] = 6; triggers[low][1] = 7; triggers[low][2] = 8;
   
   Int_t ntrig = ntrigc; //change name
   Int_t nmethods = nmethodsAll;
   if(MBEJETrigComb == 0) {
      ntrig=0; //or put it to zero
      //nmethods =1;
   }
   
   
   TString suffix="";
   
   Printf("High threshold trigger set to %s", trigname[high].Data());
   Printf("Low threshold trigger set to %s", trigname[low].Data());
   Double_t factorJ1 = GetFactorMBJetTrigger("All", trigname[0], 1);
   
   Double_t factorJ2 = GetFactorMBJetTrigger("All", trigname[1], 1);
   Printf("Factor J1 = %f, Factor J2 = %f", factorJ1, factorJ2);
   
   TCanvas *cjmptj = new TCanvas("cjmptj", "Jet mass in bins of jet pT", 1000,1000);
   cjmptj->Divide(3,2);
   TCanvas *cStatUncjmptj = new TCanvas("cStatUncjmptj", "Statistical uncertainty on jet mass in bins of jet pT", 1000,1000);
   cStatUncjmptj->Divide(3,2);
   TLegend *legSub = new TLegend(0.2, 0.4, 0.5, 0.7);
   legSub->SetFillStyle(0);
   legSub->SetBorderSize(0);
   
   TCanvas *cjmptjTrigcmp[nmethods];
   TCanvas *cjmptjTrigcmpR[nmethods];
   TCanvas *cTrigRwidepT[nmethods];
   
   for(Int_t i=0;i<nmethods;i++){
      cjmptjTrigcmp[i]= new TCanvas(Form("cjmptjTrigcmp%d",i), Form("Jet mass in bins of jet pT Comparison of triggers - Tagged jets - %s", strLeg[i].Data()), 1000,1000);
      cjmptjTrigcmp[i]->Divide(3,2);
      cjmptjTrigcmpR[i]= new TCanvas(Form("cjmptjTrigcmpR%d",i), Form("RATIOS Jet mass in bins of jet pT Comparison of triggers - Tagged jets - %s", strLeg[i].Data()), 1000,1000);
      cjmptjTrigcmpR[i]->Divide(3,2);
      cTrigRwidepT[i]= new TCanvas(Form("cTrigRwidepT%d",i), Form("RATIOS Jet mass plateau jet pT Comparison of triggers - Tagged jets - %s", strLeg[i].Data()), 600,600);
   }
   
   TLegend *legTr=new TLegend(0.15, 0.6, 0.55, 0.9, "Trigger");
   legTr->SetFillStyle(0);
   legTr->SetBorderSize(0);
   
   TFile *fout = new TFile(Form("MassOutput%s.root",suffix.Data()), "recreate");
   const Int_t nLeadTrkBins = 6;
   Double_t leadTrkPt[nLeadTrkBins] = {0.,2.,6.,10.,15.,19.};
   Int_t loopend;
   if(var == 1 ) loopend = npTjbins;
   if(var == 2 ) loopend = nLeadTrkBins;
   
   
   TH2D *hMassPt[nmethods];
   TH1D *hTurnOnJ1[nmethods];
   TH1D *hTurnOnJ2[nmethods];
   Float_t triggerNorm1 = 8.2; // this value was obtained from a pol0 fit of hTurnOnJ1
   Float_t triggerNorm2 = 0.12;// this value was obtained from a pol0 fit of hTurnOnJ2
   
   TCanvas *cTurnOn = new TCanvas("cTurnOn", "Turn On", 600, 600);
   TCanvas *cMPt = new TCanvas("cMPt", "Mass vs Pt and Projections", 600, 800);
   cMPt->Divide(1,3);
   
   for(Int_t i = 0; i< nmethods; i++){
      
      hMassPt[i] = 0x0;
      // min bias
      TH2D *hMassPtMB = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetTagged, minbias[i]);
      hMassPtMB->SetName(Form("hMassPtMBTagged_BkgSub%d", i));
      // trig high threshold
      TH2D *hMassPtJ1 = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetTagged, triggers[0][i]);
      hMassPtJ1->SetName(Form("hMassPtJ1Tagged_BkgSub%d", i));
      // trig low threshold
      TH2D *hMassPtJ2 = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetTagged, triggers[1][i]);
      hMassPtJ2->SetName(Form("hMassPtJ2Tagged_BkgSub%d", i));
      fout->cd();
      hMassPtMB->Write();
      hMassPtJ1->Write();
      hMassPtJ2->Write();
      
      hTurnOnJ1[i] = hMassPtJ1->ProjectionX(Form("hTurnOnJ1_BkgSub%d", i));
      hTurnOnJ1[i]->SetLineStyle(1);
      hTurnOnJ1[i]->SetLineColor(colors[i]);
      
      hTurnOnJ2[i] = hMassPtJ2->ProjectionX(Form("hTurnOnJ2_BkgSub%d", i));
      hTurnOnJ2[i]->SetLineStyle(2);
      hTurnOnJ2[i]->SetLineColor(colors[i]);
      TH1D *hPtMB = hMassPtMB->ProjectionX(Form("hPtMB%d",i));
      hTurnOnJ1[i]->Divide(hPtMB);
      hTurnOnJ2[i]->Divide(hPtMB);
      
      cTurnOn->cd();
      if(i == 0) {
      	 hTurnOnJ1[i]->Draw("h");
      	 //hTurnOnJ2[i]->Draw("h");
      }
      else{
      	 hTurnOnJ1[i]->Draw("hsames");
      	 //hTurnOnJ2[i]->Draw("hsames");
      }
      
      Int_t nmaxbinspt = hMassPtMB->GetXaxis()->GetNbins();
      Int_t nmaxbinsm = hMassPtMB->GetYaxis()->GetNbins();
      if(MBEJETrigComb != 0){
      	 Int_t nminbin = hMassPtMB->GetXaxis()->FindBin(thresh[high]);
      	 Printf("Bin threshold %d + 1 = %d", nminbin, nminbin+1);
      	 //remove J1 entries below thresh[high]
      	 //for(Int_t ib = 1; ib<=nminbin; ib++){
      	 for(Int_t ib = 1; ib<=nmaxbinspt; ib++){
      	    for(Int_t ibm = 1; ibm<=nmaxbinsm; ibm++){
      	       if(ib<=nminbin) hMassPtJ1->SetBinContent(ib,ibm, 0);
      	       else {
      	       	  hMassPtJ1->SetBinContent(ib,ibm, hMassPtJ1->GetBinContent(ib,ibm)/triggerNorm1);
      	       	  hMassPtJ1->SetBinError(ib,ibm, hMassPtJ1->GetBinError(ib,ibm)/triggerNorm1);
      	       }
      	    }
      	 }
      	 if(MBEJETrigComb == 1){
      	   //remove MB entries above thresh[high]
      	    for(Int_t ib = nminbin+1; ib<nmaxbinspt; ib++){
      	       for(Int_t ibm = 1; ibm<=nmaxbinsm; ibm++){
      	       	  hMassPtMB->SetBinContent(ib,ibm, 0);
      	       }
      	    }
      	    // build the final histogram
      	    hMassPt[i] = hMassPtMB;
      	    hMassPt[i] ->Add(hMassPtJ1);
      	 }
      	 if(MBEJETrigComb == 2){
      	     // build the final histogram
      	    hMassPt[i] = hMassPtMB;
      	    hMassPt[i] ->Add(hMassPtJ1);
      	 }
      	 if(MBEJETrigComb == 3){
      	    Int_t nminbin2 = hMassPtJ2->GetXaxis()->FindBin(thresh[low]);
      	    //remove J1 entries below thresh[high]
      	    for(Int_t ib = 1; ib<=nminbin2; ib++){
      	       for(Int_t ibm = 1; ibm<=nmaxbinsm; ibm++){
      	       	  hMassPtJ2->SetBinContent(ib,ibm, 0);
      	       }
      	    }
      	    hMassPt[i] = hMassPtMB;
      	    hMassPt[i] ->Add(hMassPtJ1);
      	    hMassPt[i] ->Add(hMassPtJ2);
      	    
      	 }
      } else hMassPt[i] = hMassPtMB;
      
      hMassPt[i]->SetName(Form("h2MPtTagged_BkgSub%d_TrgCmb%d", i, MBEJETrigComb));
      
      cMPt->cd(1);
      if(i == 0 ) hMassPt[i]->Draw("colz");
      
      TH1D *hPtpj = hMassPt[i]->ProjectionX(Form("%s_X", hMassPt[i]->GetName()));
      hPtpj->SetMarkerStyle(20);
      hPtpj->SetMarkerColor(colors[i]);
      
      TH1D *hMpj = hMassPt[i]->ProjectionY(Form("%s_Y", hMassPt[i]->GetName()));
      hMpj->SetMarkerStyle(20);
      hMpj->SetMarkerColor(colors[i]);
      
      cMPt->cd(2);
      if(i == 0 ) hMpj->Draw("P");
      else hMpj->Draw("Psames");
      
      cMPt->cd(3);
      gPad->SetLogy();
      if(i == 0 ) hPtpj->Draw("P");
      else hPtpj->Draw("Psames");
      
      fout->cd();
      hMassPt[i]->Write();
   }
   
   SaveCv(cMPt);
   
   TH1D *hMTaggedTrigRwidepT[nmethods][ntrig];
   TH1D *hMTaggedMBRwidepT[nmethods][ntrig];
   for(Int_t i = 0; i<nmethods; i++) {
      for(Int_t k = 0; k<ntrig; k++){
      	 hMTaggedTrigRwidepT[i][k]=0;
      	 hMTaggedMBRwidepT[i][k]=0;
      }
   }
   
   for(Int_t iptj=0 ; iptj<loopend; iptj++){
      //for(Int_t iptj=1 ; iptj<2; iptj++){
      TH1D *hMAllpTleadint[nmethods];
      TH1D *hMTaggedpTleadint[nmethods];
      TH1D *hMTaggedMatchpTleadint[nmethods];
      
      TH1D *hStatUncMTaggedpTleadint[nmethods];
      
      TH1D *hMAllPerTrig[nmethods][ntrig+1];
      TH1D *hMTaggedPerTrig[nmethods][ntrig+1];
      TH1D *hMTaggedMatchPerTrig[nmethods][ntrig+1];
      TH1D *hMAllTrigR[nmethods][ntrig];
      TH1D *hMTaggedTrigR[nmethods][ntrig];
      TH1D *hMTaggedMatchTrigR[nmethods][ntrig];
      
      
      
      TPaveText *pvpT = new TPaveText(0.4, 0.7, 0.8, 0.9, "NDC");
      pvpT->SetFillStyle(0);
      pvpT->SetBorderSize(0);
      if(var == 1) pvpT->AddText(Form("p_{T,jet} %.0f - %.0f GeV/#it{c}", pTjet[iptj],pTjet[iptj+1]));
      if(var == 2) pvpT->AddText(Form("p_{T,lead} > %.0f GeV/#it{c}", leadTrkPt[iptj]));
      
      //pT leading track integrated plots
      util->SetMinLeadTrackPt(0);
      if(var == 1) util->SetJetPtRange(pTjet[iptj], pTjet[iptj+1]);
      if(var == 2) util->SetMinLeadTrackPt(leadTrkPt[iptj]);
      Double_t norm=1;
      for(Int_t i = 0; i<nmethods; i++) {
      	 hMAllpTleadint[i]=0;
      	 hMTaggedpTleadint[i]=0;
      	 hMTaggedMatchpTleadint[i]=0;
      	 hStatUncMTaggedpTleadint[i]=0;
      	 for(Int_t k=0; k<ntrig+1; k++){
      	    
      	    hMAllPerTrig[i][k]=0;
      	    hMTaggedPerTrig[i][k]=0;
      	    hMTaggedMatchPerTrig[i][k]=0;
      	    if(k<ntrig){
      	       hMAllTrigR[i][k]=0;
      	       hMTaggedTrigR[i][k]=0;
      	       hMTaggedMatchTrigR[i][k]=0;
      	    }
      	 }
      	 //min bias
      	 if(var == 1 ){
      	    hMAllPerTrig[i][0] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,minbias[i]);
      	    hMTaggedPerTrig[i][0] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,minbias[i]);
      	    hMTaggedMatchPerTrig[i][0] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,minbias[i]);
      	    hMAllPerTrig[i][0]->SetName(Form("hMassAll%dpTj%dMB", i, iptj));
      	    hMTaggedPerTrig[i][0]->SetName(Form("hMassTagged%dpTj%dMB", i, iptj));
      	    hMTaggedMatchPerTrig[i][0]->SetName(Form("hMassTaggedMatch%dpTj%dMB", i, iptj));
      	 }
      	 if(var == 2 ){
      	    hMAllPerTrig[i][0] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetAll,minbias[i]);
      	    hMTaggedPerTrig[i][0] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTagged,minbias[i]);
      	    hMTaggedMatchPerTrig[i][0] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,minbias[i]);
      	    hMAllPerTrig[i][0]->SetName(Form("hMassAll%dpTj%dMB", i, iptj));
      	    hMTaggedPerTrig[i][0]->SetName(Form("hMassTagged%dpTj%dMB", i, iptj));
      	    hMTaggedMatchPerTrig[i][0]->SetName(Form("hMassTaggedMatch%dpTj%dMB", i, iptj));
      	    
      	    
      	 }
      	 if(i==0) Printf("Filling pT lead integrated pTj %.0f - %.0f: min bias", pTjet[iptj], pTjet[iptj+1]);
      	 
      	 
      	 //Printf("Number of entries mb %.0f ()", hMAllpTleadint[i]->GetEntries());
      	 Printf("Number of entries mb %.0f method %d", hMTaggedPerTrig[i][0]->GetEntries(), i);
      	 
      	 //at pT less than threshold 1 or for case 2 (MB+J2 all over the pT spectrum)
      	 if(MBEJETrigComb == 0 || MBEJETrigComb >= 2 || pTjet[iptj]<=thresh[high]){ //start with the minimum bias
      	    hMAllpTleadint[i]=(TH1D*)hMAllPerTrig[i][0]->Clone(Form("hMassAll%dpTj%d", i, iptj));
      	    hMTaggedpTleadint[i]=(TH1D*)hMTaggedPerTrig[i][0]->Clone(Form("hMassTagged%dpTj%d", i, iptj));
      	    hMTaggedMatchpTleadint[i]=(TH1D*)hMTaggedMatchPerTrig[i][0]->Clone(Form("hMassTaggedMatch%dpTj%d", i, iptj));
      	    
      	 }
      	 
      	 if(MBEJETrigComb){
      	    
      	    Printf("TRIGGER %s", trigname[high].Data());
      	    //J1
      	    //util->SetScaleFactor(factorJ1);
      	    if(var == 1){
      	       hMAllPerTrig[i][high+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggers[high][i]);
      	       
      	       hMTaggedPerTrig[i][high+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggers[high][i]);
      	       //util->SetScaleFactor(factorJ1);
      	       hMTaggedMatchPerTrig[i][high+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggers[high][i]);
      	       hMAllPerTrig[i][high+1]->SetName(Form("hMassAll%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       hMTaggedPerTrig[i][high+1]->SetName(Form("hMassTagged%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       hMTaggedMatchPerTrig[i][high+1]->SetName(Form("hMassTaggedMatched%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       	  
      	       
      	    }
      	    if(var == 2) {
      	       hMAllPerTrig[i][high+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetAll,triggers[high][i]);
      	       
      	       hMTaggedPerTrig[i][high+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTagged,triggers[high][i]);
      	       //util->SetScaleFactor(factorJ1);
      	       hMTaggedMatchPerTrig[i][high+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggers[high][i]);
      	       hMAllPerTrig[i][high+1]->SetName(Form("hMassAll%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       hMTaggedPerTrig[i][high+1]->SetName(Form("hMassTagged%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       hMTaggedMatchPerTrig[i][high+1]->SetName(Form("hMassTaggedMatched%dpTj%d%s", i, iptj, trigname[high].Data()));
      	       
      	    }
       	   
      	    if(var == 1){
      	       hMAllPerTrig[i][low+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggers[low][i]);
      	     
      	       hMTaggedPerTrig[i][low+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggers[low][i]);
      	       //util->SetScaleFactor(factorJ1);
      	       hMTaggedMatchPerTrig[i][low+1] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggers[low][i]);
      	       hMAllPerTrig[i][low+1]->SetName(Form("hMassAll%dpTj%d%s", i, iptj, trigname[low].Data()));
      	       hMTaggedPerTrig[i][low+1]->SetName(Form("hMassTagged%dpTj%d%s", i, iptj, trigname[low].Data()));
      	       hMTaggedMatchPerTrig[i][low+1]->SetName(Form("hMassTaggedMatched%dpTj%d%s", i, iptj, trigname[low].Data()));
      	    }
      	    if(var == 2 ){
      	       hMAllPerTrig[i][low+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetAll,triggers[low][i]);
      	       
      	       hMTaggedPerTrig[i][low+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTagged,triggers[low][i]);
      	       //util->SetScaleFactor(factorJ1);
      	       hMTaggedMatchPerTrig[i][low+1] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggers[low][i]);
      	       
      	       hMAllPerTrig[i][low+1]->SetName(Form("hMassAll%dpTj%d%s", i, iptj, trigname[low].Data()));
      	       hMTaggedPerTrig[i][low+1]->SetName(Form("hMassTagged%dpTj%d%s", i, iptj, trigname[low].Data()));
      	       hMTaggedMatchPerTrig[i][low+1]->SetName(Form("hMassTaggedMatched%dpTj%d%s", i, iptj, trigname[low].Data()));
      	       
      	       
      	    }
      	    if(pTjet[iptj]>thresh[low]){
      	       if(MBEJETrigComb == 3){
      	       	  Printf("INCLUDE TRIGGER %s in the total", trigname[low].Data());
      	       	  hMAllpTleadint[i]->Add(hMAllPerTrig[i][low+1]);
      	       	  Printf("Number of entries mb+j1 %.0f", hMAllpTleadint[i]->GetEntries());   
      	       	  //util->SetScaleFactor(factorJ1);
      	       	  hMTaggedpTleadint[i]->Add(hMTaggedPerTrig[i][low+1]);
      	       	  hMTaggedMatchpTleadint[i]->Add(hMTaggedMatchPerTrig[i][low+1]);
      	       	  
      	       	  Printf("Number of entries mb+j1 %.0f", hMTaggedpTleadint[i]->GetEntries());
      	       }
      	    } else Printf("Threshold < %f", thresh[low]);
      	    

      	    if(pTjet[iptj]>thresh[high]){   
      	       if(MBEJETrigComb >= 2) {
      	       	  Printf("INCLUDE TRIGGER %s in the total", trigname[high].Data());
      	       	  hMAllpTleadint[i]->Add(hMAllPerTrig[i][high+1]);
      	       	  Printf("Number of entries mb+j2 %.0f", hMAllpTleadint[i]->GetEntries());   
      	       	  //util->SetScaleFactor(factorJ1);
      	       	  
      	       	  hMTaggedpTleadint[i]->Add(hMTaggedPerTrig[i][high+1]);
      	       	  Printf("Number of entries mb+j2 %.0f", hMTaggedpTleadint[i]->GetEntries());
      	       	   hMTaggedMatchpTleadint[i]->Add(hMTaggedMatchPerTrig[i][high+1]);
      	       	  
      	       }
      	       if(MBEJETrigComb == 1) {
      	       	  if(hMTaggedpTleadint[i]) Printf("YOU SHOULD BE NULL");
      	       	  hMAllpTleadint[i]=(TH1D*)hMAllPerTrig[i][high+1]->Clone(Form("hMassAll%dpTj%d", i, iptj));
      	       	  hMTaggedpTleadint[i]=(TH1D*)hMTaggedPerTrig[i][high+1]->Clone(Form("hMassTagged%dpTj%d", i, iptj));
      	       	  hMTaggedMatchpTleadint[i]=(TH1D*)hMTaggedMatchPerTrig[i][high+1]->Clone(Form("hMassTaggedMatch%dpTj%d", i, iptj));
      	       }
      	       
      	    } else Printf("Threshold < %f", thresh[high]);
      	    
      	 }
      	 /*
      	 if(pTjet[iptj]>39){
      	 //J2
      	 //util->SetScaleFactor(factorJ2);
      	 hMAllPerTrig[i][2] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,triggersJ2[i]);
      	 hMAllPerTrig[i][2]->SetName(Form("%s%s", hMAllpTleadint[i]->GetName(), trigname[1].Data()));
      	 // we don't use the J2 trigger
      	 //hMAllpTleadint[i]->Add(hMAllPerTrig[i][2] );
      	 
      	 //util->SetScaleFactor(factorJ2);
      	 hMTaggedPerTrig[i][2] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,triggersJ2[i]);
      	 hMTaggedPerTrig[i][2]->SetName(Form("%s%s", hMTaggedpTleadint[i]->GetName(), trigname[1].Data()));
      	 
      	 //hMTaggedpTleadint[i]->Add(hMTaggedPerTrig[i][2] );
      	 //util->SetScaleFactor(factorJ2);
      	 hMTaggedMatchPerTrig[i][2] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,triggersJ2[i]);
      	 hMTaggedMatchPerTrig[i][2]->SetName(Form("%s%s", hMTaggedMatchpTleadint[i]->GetName(), trigname[1].Data()));
      	 if(MBEJETrigComb == 3){
      	 hMAllpTleadint[i]->Add(hMAllPerTrig[i][2]);
      	 hMTaggedpTleadint[i]->Add(hMTaggedPerTrig[i][2]);
      	 hMTaggedMatchpTleadint[i]->Add(hMTaggedMatchPerTrig[i][2]);
      	 }
      	 //if(i==0) printf(" + J2 (lt)");
      	 }
      	 */
      	 if(i==0) printf("\n");
      	 
     
      //statistical uncertainty
      hStatUncMTaggedpTleadint[i] = (TH1D*)hMTaggedpTleadint[i]->Clone(Form("hStatUncMTagged%dpTj%d", i, iptj));
      
      Double_t finalWidth = 2.;
      Printf("Bin width now %f", hStatUncMTaggedpTleadint[i]->GetBinWidth(3));
      Int_t finalrebin = (Int_t)(finalWidth/hStatUncMTaggedpTleadint[i]->GetBinWidth(3));
      hStatUncMTaggedpTleadint[i]->Rebin(finalrebin);
      Printf("Rebin %d", finalrebin);
      for(Int_t ibin = 0; ibin < hStatUncMTaggedpTleadint[i]->GetNbinsX(); ibin++){
      	 Double_t cont = hStatUncMTaggedpTleadint[i]->GetBinContent(ibin+1), err = hStatUncMTaggedpTleadint[i]->GetBinError(ibin+1);
      	 //Printf("Filling bin %d with %f", ibin+1, err/cont);
      	 if (cont>0.) {
      	    hStatUncMTaggedpTleadint[i]->SetBinContent(ibin+1, err/cont);
      	    hStatUncMTaggedpTleadint[i]->SetBinError(ibin+1, 0);
      	 }
      	 
      }
      
      //if(i==0) norm = hMAllpTleadint[i]->Integral();
      //hMAllpTleadint[i]->Scale(1./norm);
      hMAllpTleadint[i]->SetLineColor(colors[i]);
      hMAllpTleadint[i]->SetMarkerColor(hMAllpTleadint[i]->GetLineColor());
      hMAllpTleadint[i]->SetMarkerStyle(20);
      //hMAllpTleadint[i]->GetYaxis()->SetRangeUser(0.,0.4);
      if(var == 1) hMAllpTleadint[i]->GetXaxis()->SetRangeUser(-5,20.);
      
      //hMTaggedpTleadint[i]->Scale(1./norm);
      hMTaggedpTleadint[i]->SetLineColor(colors[i]);
      hMTaggedpTleadint[i]->SetMarkerColor(hMTaggedpTleadint[i]->GetLineColor());
      hMTaggedpTleadint[i]->SetMarkerStyle(24);
      
      //hMTaggedMatchpTleadint[i]->Scale(1./norm);
      hMTaggedMatchpTleadint[i]->SetLineColor(colors[i]);
      hMTaggedMatchpTleadint[i]->SetMarkerColor(hMTaggedMatchpTleadint[i]->GetLineColor());
      hMTaggedMatchpTleadint[i]->SetMarkerStyle(27);
      
      hStatUncMTaggedpTleadint[i]->SetLineColor(colors[i]);
      hStatUncMTaggedpTleadint[i]->SetFillColor(colors[i]);
      hStatUncMTaggedpTleadint[i]->SetMarkerColor(hStatUncMTaggedpTleadint[i]->GetLineColor());
      //hStatUncMTaggedpTleadint[i]->SetMarkerStyle(20);
      if(var == 1){
      	 hStatUncMTaggedpTleadint[i]->GetXaxis()->SetRangeUser(-5,24.);
      	 hStatUncMTaggedpTleadint[i]->GetYaxis()->SetRangeUser(0.,1.);
      }
      fout->cd();
      hMAllpTleadint[i]            ->Write();
      hMTaggedpTleadint[i]         ->Write();  
      hMTaggedMatchpTleadint[i]    ->Write();
      hStatUncMTaggedpTleadint[i]  ->Write();
   
      Double_t normTr = 1;
      for(Int_t k=0; k<ntrig+1; k++){
      	 
      	 Printf("++++++++++++ K=%d", k);
      	 if(hMAllPerTrig[i][k]){
      	    //if(k==0) normTr = hMTaggedPerTrig[i][k]->Integral();
      	    //hMAllPerTrig[i][k]->Scale(1./normTr/hMAllPerTrig[i][k]->Integral());
      	    hMAllPerTrig[i][k]->Scale(1./hMAllPerTrig[i][k]->Integral());
      	    hMAllPerTrig[i][k]->Rebin(2);
      	    hMAllPerTrig[i][k]->SetLineColor(colors[k+nmethods]);
      	    hMAllPerTrig[i][k]->SetMarkerColor(hMAllPerTrig[i][k]->GetLineColor());
      	    hMAllPerTrig[i][k]->SetMarkerStyle(20);
      	    //hMAllPerTrig[i]->GetYaxis()->SetRangeUser(0.,0.4);
      	    if(var == 1) hMAllPerTrig[i][k]->GetXaxis()->SetRangeUser(-5,20.);
      	    if(var == 2) hMAllPerTrig[i][k]->GetXaxis()->SetRangeUser(0,140.);
      	    
      	    if(i==0 && iptj==4) legTr->AddEntry(hMTaggedPerTrig[i][k], Form("%s", (k==0) ? "MB" : trigname[k-1].Data()), "PL");
      	    //hMTaggedPerTrig[i][k]->Scale(1./normTr/hMTaggedPerTrig[i][k]->Integral());
      	    hMTaggedPerTrig[i][k]->Scale(1./hMTaggedPerTrig[i][k]->Integral());
      	    hMTaggedPerTrig[i][k]->Rebin(2);
      	    hMTaggedPerTrig[i][k]->SetLineColor(colors[k+nmethods]);
      	    hMTaggedPerTrig[i][k]->SetMarkerColor(hMTaggedPerTrig[i][k]->GetLineColor());
      	    hMTaggedPerTrig[i][k]->SetMarkerStyle(24);
      	    if(var == 1) {
      	       hMTaggedPerTrig[i][k]->GetXaxis()->SetRangeUser(-5,20.);
      	       hMTaggedPerTrig[i][k]->GetYaxis()->SetRangeUser(0,0.6);
      	    }
      	    if(var == 2) {
      	       hMTaggedPerTrig[i][k]->GetXaxis()->SetRangeUser(0,140.);
      	       hMTaggedPerTrig[i][k]->GetYaxis()->SetRangeUser(1e-7,1.);
      	    }
      	    //hMTaggedMatchPerTrig[i][k]->Scale(1./normTr/hMTaggedMatchPerTrig[i][k]->Integral());
      	    hMTaggedMatchPerTrig[i][k]->Scale(1./hMTaggedMatchPerTrig[i][k]->Integral());
      	    hMTaggedMatchPerTrig[i][k]->Rebin(2);
      	    hMTaggedMatchPerTrig[i][k]->SetLineColor(colors[k+nmethods]);
      	    hMTaggedMatchPerTrig[i][k]->SetMarkerColor(hMTaggedMatchPerTrig[i][k]->GetLineColor());
      	    hMTaggedMatchPerTrig[i][k]->SetMarkerStyle(27);
      	    if(var == 1) hMTaggedMatchPerTrig[i][k]->GetXaxis()->SetRangeUser(-5,20.);
      	    if(var == 2) hMTaggedMatchPerTrig[i][k]->GetXaxis()->SetRangeUser(0,140.);
      	    
      	    if(k>0){
      	       //Printf("Ratio ipt %d, k %d, i %d", iptj, k, i);
      	       hMAllTrigR[i][k-1] = (TH1D*)hMAllPerTrig[i][k]->Clone(Form("h%sR%s", hMAllPerTrig[i][k]->GetName(), trigname[k-1].Data()));
      	       if(hMAllTrigR[i][k-1]) hMAllTrigR[i][k-1]->Divide(hMAllPerTrig[i][0]);
      	       
      	       hMTaggedTrigR[i][k-1] = (TH1D*)hMTaggedPerTrig[i][k]->Clone(Form("h%sR%s", hMTaggedPerTrig[i][k]->GetName(), trigname[k-1].Data()));
      	       if(hMTaggedTrigR[i][k-1]) hMTaggedTrigR[i][k-1]->Divide(hMTaggedPerTrig[i][0]);
      	       
      	       
      	       if(var == 1) {
      	       	  hMTaggedTrigR[i][k-1]->GetXaxis()->SetRangeUser(-5,20.);
      	       	  hMTaggedTrigR[i][k-1]->GetYaxis()->SetRangeUser(0.5,1.5);
      	       }
      	       if(var == 2) {
      	       	  hMTaggedTrigR[i][k-1]->GetXaxis()->SetRangeUser(0,140.);
      	       	  hMTaggedTrigR[i][k-1]->GetYaxis()->SetRangeUser(0,10.);
      	       }
      	       
      	       
      	       Printf("000000000000000 Clone %s to %s and divide by %s", hMTaggedPerTrig[i][k]->GetName(), hMTaggedTrigR[i][k-1]->GetName(), hMTaggedPerTrig[i][0]->GetName());
      	       
      	       
      	       hMTaggedMatchTrigR[i][k-1] = (TH1D*)hMTaggedMatchPerTrig[i][k]->Clone(Form("h%sR%s", hMTaggedMatchPerTrig[i][k]->GetName(), trigname[k-1].Data()));
      	       if(hMTaggedMatchTrigR[i][k-1]) hMTaggedMatchTrigR[i][k-1]->Divide(hMTaggedMatchPerTrig[i][0]);
      	       
      	       
      	       if(pTjet[iptj] >= thresh[k-1]) {
      	       	  if(!hMTaggedTrigRwidepT[i][k-1]) {
      	       	     hMTaggedTrigRwidepT[i][k-1] = (TH1D*)hMTaggedPerTrig[i][k]->Clone(Form("hRlowpThMTaggedTrigRwidepT%d%.0f", i, thresh[k-1]));
      	       	     hMTaggedTrigRwidepT[i][k-1]->SetMarkerColor(hMTaggedPerTrig[i][k]->GetMarkerColor());
      	       	     hMTaggedTrigRwidepT[i][k-1]->SetLineColor(hMTaggedPerTrig[i][k]->GetMarkerColor());
      	       	     Printf(">>>>>>>>>> Created %s", hMTaggedTrigRwidepT[i][k-1]->GetName());
      	       	  }
      	       	  else hMTaggedTrigRwidepT[i][k-1]->Add(hMTaggedPerTrig[i][k]);
      	       }
      	    } else { //k = 0 MB
      	       if(MBEJETrigComb){
      	       	  if(pTjet[iptj] >= thresh[low]) {
      	       	     if(!hMTaggedMBRwidepT[i][low]) {
      	       	     	hMTaggedMBRwidepT[i][low] = (TH1D*)hMTaggedPerTrig[i][k]->Clone(Form("hMTaggedMBRwidepTlowpT%d%.0f", i, thresh[low]));
      	       	     	Printf(">>>>>>>>>> Created %s", hMTaggedMBRwidepT[i][low]->GetName());
      	       	     }
      	       	     else{
      	       	     	Printf("Pointers %p, %p -- i = %d, k = %d, low = %d", hMTaggedMBRwidepT[i][low], hMTaggedPerTrig[i][k], i, k, low);
      	       	     	hMTaggedMBRwidepT[i][low]->Add(hMTaggedPerTrig[i][k]);
      	       	     }
      	       	  } 
      	       	  if(pTjet[iptj] >= thresh[high]) {
      	       	     if(!hMTaggedMBRwidepT[i][high]) {
      	       	     	hMTaggedMBRwidepT[i][high] = (TH1D*)hMTaggedPerTrig[i][k]->Clone(Form("hMTaggedMBRwidepTlowpT%d%.0f", i, thresh[high]));
      	       	     	Printf(">>>>>>>>>> Created %s", hMTaggedMBRwidepT[i][high]->GetName());
      	       	     } else hMTaggedMBRwidepT[i][high]->Add(hMTaggedPerTrig[i][k]);
      	       	     
      	       	     
      	       	  }
      	       } //else Printf("Not found ipt %d, k %d, i %d", iptj, k, i);
      	       
      	    }
      	 }
      }
      
      //rebin raw counts
      hMAllpTleadint[i]->Rebin(2);
      hMTaggedpTleadint[i]->Rebin(2);
      cjmptj->cd(iptj+1);
      if(var == 2 ) gPad->SetLogy();
      if(var == 1 ) {
      	 hMAllpTleadint[i]->GetXaxis()->SetRangeUser(-5, 20);
      }
      if(i==0) hMAllpTleadint[i]->Draw("P");
      else hMAllpTleadint[i]->Draw("Psames");
      hMTaggedpTleadint[i]->Draw("Psames");
      hMTaggedMatchpTleadint[i]->Draw("Psames");
      if(iptj == 0) legSub->AddEntry(hMTaggedpTleadint[i], Form("%s", methnames[i].Data()), "P");
      
      cStatUncjmptj->cd(iptj+1);
      gPad->SetGridx(); gPad->SetGridy();
      hStatUncMTaggedpTleadint[i]->Draw("histsames");
      
      
      
      for(Int_t k=0; k<ntrig; k++){
      	 //hMAllPerTrig[i][0]->Draw("sames");
      	 /*
      	 */
      	 //Printf("i = %d, k = %d, pointers %p, %p", i, k, hMTaggedPerTrig[i][k], hMTaggedTrigR[i][k]);
      	 //cjmptjTrigcmp[i]->cd(iptj+1);
      	 //if(hMTaggedPerTrig[i][k]) hMTaggedPerTrig[i][k]->Draw("sames");
      	 cjmptjTrigcmp[i]->cd(iptj+1);
      	 if(var == 2 ) gPad->SetLogy();
      	 if(hMTaggedPerTrig[i][k]){
      	    hMTaggedPerTrig[i][k]->Draw("sames");
      	    legTr->Draw();
      	    pvpT->Draw();
      	 }
      
      	 cjmptjTrigcmpR[i]->cd(iptj+1);
      	 if(hMTaggedTrigR[i][k]) {//k<ntrig+1 && 
      	    hMTaggedTrigR[i][k]->Draw("sames");
      	    legTr->Draw();
      	    pvpT->Draw();
      	 }
      }
      
      SaveCv(cjmptjTrigcmp[i]);
      SaveCv(cjmptjTrigcmpR[i]);
   }
   cjmptj->cd(iptj+1);
   pvpT->Draw();
   
   cStatUncjmptj->cd(iptj+1);
   pvpT->Draw();
   }
   cjmptj->cd(1);
   legSub->Draw();
   
   SaveCv(cjmptj,suffix.Data());
   SaveCv(cStatUncjmptj);
   for(Int_t i = 0 ; i < nmethods; i++){
      for(Int_t k=0; k<ntrig; k++){
      	 cTrigRwidepT[i]->cd();
      	 if(hMTaggedMBRwidepT[i][k] && hMTaggedMBRwidepT[i][k]->Integral()>0){
      	    //hMTaggedMBRwidepT[i][k]->Rebin(2);
      	    hMTaggedMBRwidepT[i][k]->Scale(1./hMTaggedMBRwidepT[i][k]->Integral());
      	 
      	  //  hMTaggedMBRwidepT[i][k]->Draw("sames");
      	 }
      	 if(hMTaggedTrigRwidepT[i][k] && hMTaggedTrigRwidepT[i][k]->Integral()>0) {
      	    //hMTaggedTrigRwidepT[i][k]->Rebin(2);
      	    hMTaggedTrigRwidepT[i][k]->Scale(1./hMTaggedTrigRwidepT[i][k]->Integral());
      	    hMTaggedTrigRwidepT[i][k]->Divide(hMTaggedMBRwidepT[i][k]);
      	    hMTaggedTrigRwidepT[i][k]->GetYaxis()->SetRangeUser(0,2);
      	    hMTaggedTrigRwidepT[i][k]->Draw("sames");
      	    TF1 *fitpol0 = new TF1("fitpol0", "pol0", -5, 20);
      	    hMTaggedTrigRwidepT[i][k]->Fit("fitpol0","R");
      	    TPaveText *pstat = new TPaveText(0.3, 0.7, 0.9, 0.8, "NDC");
      	    pstat->SetBorderSize(0);
      	    pstat->SetFillStyle(0);
      	    pstat->AddText(Form("Par = %.2f +- %.2f; Prob = %f", fitpol0->GetParameter(0), fitpol0->GetParError(0), fitpol0->GetProb()));
      	    pstat->Draw();
      	    /*GetPtob :
      	     Computation of the probability for a certain Chi-squared (chi2) and number of degrees of freedom (ndf).

      	     Calculations are based on the incomplete gamma function P(a,x), where a=ndf/2 and x=chi2/2.
      	     P(a,x) represents the probability that the observed Chi-squared for a correct model should be less than the value chi2.
      	    
      	     The returned probability corresponds to 1-P(a,x),
      	     which denotes the probability that an observed Chi-squared exceeds the value chi2 by chance, even for a correct model.
      	    */
      	    
      	 }
      	 
      }
      SaveCv(cTrigRwidepT[i]);
   }
}

void OverlapPtMProjectionsFrom2D(TString f1name, TString f2name, Int_t bkgcode1, Int_t bkgcode2, Int_t trigger1, Int_t trigger2){
   
   Int_t nfiles = 2;
   TString names[nfiles] = {f1name, f2name};
   TString basehname = "h2MPtTagged_BkgSub";
   Int_t bkgcode[nfiles] = {bkgcode1, bkgcode2};
   Int_t trigger[nfiles] = {trigger1, trigger2};
   
   TCanvas *cAll = new TCanvas(Form("cCmpBkg%d_%d_Trg%d_%d", bkgcode1, bkgcode2, trigger1, trigger2), Form("Comparison Bkg%d_%d_Trg%d_%d", bkgcode1, bkgcode2, trigger1, trigger2), 800, 800);
   cAll->Divide(2,2);
   
   for(Int_t ifi = 0; ifi< nfiles; ifi++){
      TFile *fin = new TFile(names[ifi]);
      if(!fin->IsOpen()){
      	 Printf("File %s not found", names[ifi].Data());
      	 return;
      }
      
      TString histoname = Form("%s%d_TrgCmb%d", basehname.Data(), bkgcode[ifi], trigger[ifi]);
      
      TH2D *h2 =(TH2D*)fin->Get(histoname);
      if(!h2){
      	 Printf("Histogram %s not found, check settings", histoname.Data());
      	 return;
      }
      cAll->cd(1+ifi);
      h2->Draw("colz");
      
      TH1D *hPt = h2->ProjectionX(Form("%s_Pt", histoname.Data()));
      TH1D *hM = h2->ProjectionY(Form("%s_M", histoname.Data()));
      hPt->SetLineColor(colors[ifi]);
      hM ->SetLineColor(colors[ifi]);
      cAll->cd(3);
      gPad->SetLogy();
      if(ifi==0) hPt->Draw();
      else hPt->Draw("sames");
      cAll->cd(4);
      gPad->SetLogy();
      if(ifi==0) hM->Draw();
      else hM->Draw("sames");
      
      
   }
   SaveCv(cAll);
}

//------------------------------------------------------------------------------------------

void CorrectMasspTjbinbybin(TString suffixout = "", TString mcCorrFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/UnfoldingMatrix.root", TString cutMaxUnc = "Less40", TString dataFilePath = ".", TString pathCmpUnfold = "/data/cernbox/JetMassAnalysis/Data/pPb/150501/", Bool_t ispp = kFALSE){
   // pathCmpUnfold = "/data/cernbox/JetMassAnalysis/Data/pPb/Train792-793MB/unfold/Area" match with data sample output794-795 (no rhom)
   
   SetStyle();
   TString dataFileName = "MassOutput.root";
   TString fileUnfoldMeanM = "MeanJetMassAllSyst_Area.root", fileUnfoldMassDistr = "DistJetMassAllSyst_Area.root";
   TString fileUnfoldMedianM = "MedianJetMassAllSyst_Area.root";
   
   TFile *fdata = new TFile(Form("%s/%s", dataFilePath.Data(), dataFileName.Data()));
   TFile *fMC   = new TFile(Form("%s", mcCorrFile.Data()));
   
   if(!fdata->IsOpen() || !fMC->IsOpen()){
      Printf("Input not available");
      return;
   }
   
   TFile *fUnfMean = new TFile(Form("%s/%s", pathCmpUnfold.Data(), fileUnfoldMeanM.Data()));
   TFile *fUnfDistr = new TFile(Form("%s/%s", pathCmpUnfold.Data(), fileUnfoldMassDistr.Data()));
   TFile *fUnfMedian = new TFile(Form("%s/%s", pathCmpUnfold.Data(), fileUnfoldMedianM.Data()));
   TGraphErrors *gMeanMassUnf=0x0;
   TGraphErrors *gMeanMassUnfSyst=0x0;
   TGraphErrors *gMedianMassUnf=0x0;
   TGraphErrors *gMedianMassMC=0x0;
   if(fUnfMean->IsOpen()){
      gMeanMassUnf     = (TGraphErrors*)fUnfMean->Get("grUnfMean");
      gMeanMassUnf     ->SetTitle(";#it{p}_{T, jet};#LT#it{M}_{jet}#GT");
      gMeanMassUnfSyst = (TGraphErrors*)fUnfMean->Get("grUnfMeanSyst"); //"grSyst2"
      gMeanMassUnfSyst ->SetTitle(";#it{p}_{T, jet};#LT#it{M}_{jet}#GT");
      if(fUnfMedian) {
      	 gMedianMassUnf   = (TGraphErrors*)fUnfMedian->Get("grUnfMedian");
      	 gMedianMassMC    = (TGraphErrors*)fUnfMedian->Get("grMedianPyt");
      }
   }
   const Int_t nmethodsAll = 3;
   TString methodname[nmethodsAll] = {"ConstSub", "Area-MRaw", "Area-MDeriv"};
   Int_t nmethods = nmethodsAll;
   if(ispp) {
      nmethods = 1;
      methodname[0] = "Area-MRaw";
   }
   const Int_t ntags = 2;
   TString tagsname[ntags] = {"All", "Tagged"};//, "TaggedMatch"
   TString hDataName = "hMass";
   TString hCorrName = Form("hMassPartToDetc%s", cutMaxUnc.Data());
   Int_t bigcvX = 1000, bigcvY = 600;
   TCanvas *cCorrection = new TCanvas(Form("cCorrection"), Form("Correction factor"), bigcvX, bigcvY);
   cCorrection->Divide(3,2);
   //add on top the ratio of the corrected/uncorrected to check that it match with the initial factor
   //cCorrection[i][j]->Divide(2,2);
   
   TLegend *legMethMC = new TLegend(0.6, 0.4, 0.9, 0.8);
   legMethMC->SetBorderSize(0);
   legMethMC->SetFillStyle(0);
   
   TCanvas *cMassCorr[nmethods][ntags];
   TCanvas *cMass[nmethods][ntags];
   TCanvas *cMassCheck[nmethods][ntags];
   TCanvas *cMassRatios[nmethods][ntags];
   TCanvas *cMeanMass[nmethods][ntags];
   TCanvas *cMassClosure[nmethods][ntags];
   TCanvas *cMassRelUnc[ntags];
   TCanvas *cMassCorrMethMC[ntags];
   
   for(Int_t j = 0; j< ntags; j++){
      cMassRelUnc[j] = new TCanvas(Form("cMassRelUnc%s", tagsname[j].Data()), Form("Relative uncertainty (%s)",  tagsname[j].Data()), bigcvX, bigcvY);
      cMassRelUnc[j]->Divide(3,2);
      
      cMassCorrMethMC[j] = new TCanvas(Form("cMassCorrMethMC%s", tagsname[j].Data()), Form("Corrected AreaDeriv and Const plus MC part (%s)",  tagsname[j].Data()), bigcvX, bigcvY);
      cMassCorrMethMC[j]->Divide(3,2);
   }
   TGraphErrors *gMeanMass[nmethods][ntags];
   TGraphErrors *gMeanMassCor[nmethods][ntags];
   TGraphErrors *gMedianMass[nmethods][ntags];
   TGraphErrors *gMedianMassCor[nmethods][ntags];
   
   TGraphErrors *gMeanMassDMC = (TGraphErrors*)fMC->Get("gMeanMassDetc");
   gMeanMassDMC->SetTitle("; #it{p}_{T, jet};#LT#it{M}_{jet}#GT");
   TGraphErrors *gMeanMassPMC = (TGraphErrors*)fMC->Get("gMeanMassPart");
   gMeanMassPMC->SetTitle("; #it{p}_{T, jet};#LT#it{M}_{jet}#GT");
   
   for(Int_t i = 0; i < nmethods; i++){
      for(Int_t j = 0; j< ntags; j++){

      	 cMassCorr[i][j] = new TCanvas(Form("cMassCorr%s%s", tagsname[j].Data(), methodname[i].Data()), Form("Mass Corrected %s (%s)", methodname[i].Data(), tagsname[j].Data()), bigcvX, bigcvY);
      	 cMassCorr[i][j]->Divide(3,2);
      	 //cMassCorr[i][j]->Divide(2,2);
      	 cMass[i][j] = new TCanvas(Form("cMass%s%s", tagsname[j].Data(), methodname[i].Data()), Form("Mass Uncorrected %s (%s)", methodname[i].Data(), tagsname[j].Data()), bigcvX, bigcvY);
      	 cMass[i][j]->Divide(3,2);
      	 cMassCheck[i][j] = new TCanvas(Form("cMassCheck%s%s", tagsname[j].Data(), methodname[i].Data()), Form("Mass Check %s (%s)", methodname[i].Data(), tagsname[j].Data()), bigcvX, bigcvY);
      	 cMassCheck[i][j]->Divide(3,2);
      	 
      	 cMassRatios[i][j] = new TCanvas(Form("cMassRatios%s%s", tagsname[j].Data(), methodname[i].Data()), Form("Ratios Mass distributions %s (%s)", methodname[i].Data(), tagsname[j].Data()), bigcvX, bigcvY);
      	 cMassRatios[i][j]->Divide(3,2);
      	 
      	 cMassClosure[i][j] = new TCanvas(Form("cMassClosure%s%s", tagsname[j].Data(), methodname[i].Data()), Form("MC closure %s (%s)", methodname[i].Data(), tagsname[j].Data()), bigcvX, bigcvY);
      	 cMassClosure[i][j]->Divide(3,2);
      	 
      	 
      	 
      	 
      	 cMeanMass[i][j] = new TCanvas(Form("cMeanMass%s%s", tagsname[j].Data(), methodname[i].Data()), Form("Mean Mass %s (%s)", methodname[i].Data(), tagsname[j].Data()));
      	 gMeanMass[i][j] = new TGraphErrors(npTjbins);
      	 gMeanMass[i][j]->SetName(Form("gMeanMass%s%s",tagsname[j].Data(), methodname[i].Data()));
      	 gMeanMass[i][j]->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
      	 gMeanMass[i][j]->SetMarkerStyle(20);
       	 gMeanMassCor[i][j] = new TGraphErrors(npTjbins);
      	 gMeanMassCor[i][j]->SetName(Form("gMeanMassCor%s%s",tagsname[j].Data(), methodname[i].Data()));
      	 gMeanMassCor[i][j]->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
      	 gMeanMassCor[i][j]->SetMarkerStyle(21);
      	 
      	 gMedianMass[i][j] = new TGraphErrors(npTjbins);
      	 gMedianMass[i][j]->SetName(Form("gMedianMass%s%s",tagsname[j].Data(), methodname[i].Data()));
      	 gMedianMass[i][j]->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
      	 gMedianMass[i][j]->SetMarkerStyle(24);
       	 gMedianMassCor[i][j] = new TGraphErrors(npTjbins);
      	 gMedianMassCor[i][j]->SetName(Form("gMedianMassCor%s%s",tagsname[j].Data(), methodname[i].Data()));
      	 gMedianMassCor[i][j]->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
      	 gMedianMassCor[i][j]->SetMarkerStyle(25);
      	 
      	 
     }
   }
   
   TLegend *legCorrUncorr = new TLegend(0.15, 0.6, 0.45, 0.8);
   legCorrUncorr->SetFillStyle(0);
   legCorrUncorr->SetBorderSize(0);
   legCorrUncorr->AddEntry(gMeanMass[0][0], "Uncorrected", "P");
   legCorrUncorr->AddEntry(gMeanMassCor[0][0], "Corrected", "P");
   
   TCanvas *cMeanMassUncoMethodCmp = new TCanvas("cMeanMassUncoMethodCmp", "Uncorrected Mean mass Tagged Jets for different Bkg sub");
   TCanvas *cMeanMassCorrMethodCmp = new TCanvas("cMeanMassCorrMethodCmp", "Corrected Mean mass Tagged Jets for different Bkg sub");
   
   TLegend *legDistr = new TLegend(0.6,0.4, 0.9, 0.8);
   legDistr->SetBorderSize(0);
   legDistr->SetFillStyle(0);
   TLegend *legDistrCor = (TLegend*)legDistr->Clone();
   TFile *fres = new TFile(Form("Results%s.root", suffixout.Data()), "recreate");
   gDirectory->cd(0);
   
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      
      TPaveText *pvpT = new TPaveText(0.2, 0.7, 0.55, 0.9, "NDC");
      pvpT->SetFillStyle(0);
      pvpT->SetBorderSize(0);
      pvpT->AddText(Form("%.0f < p_{T,jet} < %.0f GeV/#it{c}", pTjet[iptj],pTjet[iptj+1]));

      TGraphErrors* gMassDistrUnf = (TGraphErrors*)fUnfDistr->Get(Form("grUnfMass_PtBin%d",iptj)); //pTbin 1 is 40-60 GeV/c
      TH1D* hMassDistrUnf=0x0;
      if(gMassDistrUnf) {
      	 gMassDistrUnf->SetLineColor(kOrange+7);
      	 gMassDistrUnf->SetMarkerColor(kOrange+7);
      	 hMassDistrUnf = (TH1D*)gMassDistrUnf->GetHistogram();
      	 Printf("npoints = %d, N bins %d", gMassDistrUnf->GetN(), hMassDistrUnf->GetNbinsX());
      	 if (iptj == 2) {
      	    legCorrUncorr->AddEntry(gMassDistrUnf, "Unfolded", "P");
      	    legDistrCor->AddEntry(gMassDistrUnf, "Unfolded", "P");
      	 }
      }

      TH1D* hMassDetc = (TH1D*)fMC->Get(Form("hMassDetc_pTjDetc%d",iptj));
      TH1D* hMassPart = (TH1D*)fMC->Get(Form("hMassPart_pTjPart%d",iptj));
      hMassDetc->Rebin(2);
      hMassPart->Rebin(2);
      hMassDetc->Scale(1./hMassDetc->Integral("width"));
      hMassPart->Scale(1./hMassPart->Integral("width"));
      
      
      
      if(iptj==0) {
      	 legDistrCor->AddEntry(hMassPart, "MC Part", "P");
      	 legDistr->AddEntry(hMassDetc, "MC Detc", "P");
      }
      
      TString hcorrnamept = Form("%s_pTj%d", hCorrName.Data(), iptj);
      Printf("Ratio name %s ", hcorrnamept.Data());
      	    
      TH1D* hcorr = (TH1D*)fMC->Get(hcorrnamept);
      if(!hcorr){
      	 Printf("Not found!");
      	 continue;
      }
      	    
      //hcorr->Scale(1./hcorr->Integral());
      cCorrection->cd(iptj+1);
      hcorr->Draw();
      
      Printf("%s, nbins %d, width = %f", hcorr->GetName(), hcorr->GetNbinsX(), hcorr->GetBinWidth(3));
      
      hcorr->Rebin(2);
      Int_t nbinCorr = hcorr->GetNbinsX();
      Printf("Warning! hcorr was rebinned!, now %f", hcorr->GetBinWidth(3));
      Double_t minCorr = hcorr->GetBinLowEdge(1);
      
      TH1D* hMassDetcRangeMC = new TH1D(Form("%sRangeMC",hMassDetc->GetName()), Form("%s ; Mass; ", hMassDetc->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1)); 
      hMassDetcRangeMC->SetMarkerStyle(hMassDetc->GetMarkerStyle());
      hMassDetcRangeMC->SetMarkerColor(hMassDetc->GetMarkerColor());
      TH1D* hMassPartRangeMC = new TH1D(Form("%sRangeMC",hMassPart->GetName()), Form("%s ; Mass; ", hMassPart->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1));
      hMassPartRangeMC->SetMarkerStyle(hMassPart->GetMarkerStyle());
      //hMassPartRangeMC->SetMarkerColor(hMassPart->GetMarkerColor());
      hMassPartRangeMC->SetMarkerColor(colors[nmethods+1]);
      //TH1D* hMassUnfoldRangeMC = new TH1D(Form("%sRangeMC",gMassDistrUnf->GetName()), gMassDistrUnf->GetTitle(), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1));
      Bool_t bRerangeMC = kTRUE;
      /* //moved below
      for(Int_t k=0;k<nbinCorr;k++){
      	 hMassDetcRangeMC->SetBinContent(k+1, hMassDetc->GetBinContent(k+1));
      	 hMassDetcRangeMC->SetBinError(k+1, hMassDetc->GetBinError(k+1));
      	 hMassPartRangeMC->SetBinContent(k+1, hMassPart->GetBinContent(k+1));
      	 hMassPartRangeMC->SetBinError(k+1, hMassPart->GetBinError(k+1));
      	 //if(hMassDistrUnf){
      	 //   hMassDistrUnf->SetBinContent(k+1, hMassDistrUnf->GetBinContent(k+1));
      	 //   hMassDistrUnf->SetBinError(k+1, hMassDistrUnf->GetBinError(k+1));
      	 //}

      }
      */
      
      for(Int_t i = 0; i < nmethods; i++){
      	 for(Int_t j = 0; j< ntags; j++){
      	    //for(Int_t i = 0; i < 1; i++){
      	 //for(Int_t j = 0; j< 1; j++){
      	    
      	    TString hdataname = Form("%s%s%dpTj%d", hDataName.Data(), tagsname[j].Data(), i, iptj);
      	    Printf("Histo name %s ", hdataname.Data());
      	    
      	    TH1D* hmass = (TH1D*)fdata->Get(hdataname);
      	    hmass->GetXaxis()->SetRangeUser(-5,20.);
      	    if(!hmass){
      	       Printf("Not found!");
      	       continue;
      	    }
      	    hmass->Sumw2();
      	    hmass->Scale(1./hmass->Integral("width"));
      	    hmass->SetMarkerStyle(33);
      	    
      	    cMassCheck[i][j]->cd(iptj+1);
      	    hmass->DrawClone("Psames");
      	    
      	    hMassDetc->DrawClone("Psames");
      	    hMassDetcRangeMC->SetMarkerStyle(25);
      	    hMassDetcRangeMC->DrawClone("Psames");
      	    hMassDetcRangeMC->SetMarkerStyle(hMassDetc->GetMarkerStyle());
      	    hMassDetcRangeMC->SetMarkerColor(hMassDetc->GetMarkerColor());
      	    pvpT->Draw();
      	    
      	    Printf("hMASS bin width %f", hmass->GetBinWidth(3));
      	    

      	    //TH1D *hmassReUncorrected = (TH1D*) hcorr->Clone(Form("%sReUCor",hdataname.Data()));
      	    //hmassReUncorrected->SetMarkerStyle(22);
      	    //Double_t norm = 1.; 
      	    //if(gMassDistrUnf) norm = gMassDistrUnf->Integral();
      	    //else norm = hmass->Integral();
      	    if(!TMath::AreEqualRel(hmass->GetBinWidth(3), hcorr->GetBinWidth(3), 1.E-12)){
      	       Double_t binWm = hmass->GetBinWidth(3), binWc = hcorr->GetBinWidth(3);
      	       Printf("ERRORRRR!! different bin width! hmass = %.3f, hcorr = %.3f ", binWm, binWc);
      	       Double_t rebD = binWm/binWc;
      	       Int_t reb = (Int_t)rebD;
      	       if(rebD-reb > 1.E-12){
      	       	  Printf("Fuck... not possible to rebin properly");
      	       } else {
      	       	  
      	       	  hcorr->Rebin(reb);
      	       	  hcorr->Scale(1./hcorr->Integral("width"));
      	       	  if(bRerangeMC){ //a bit dangerous, incremented below
      	       	     hMassDetc->Rebin(reb);
      	       	     hMassDetc->Scale(1./hcorr->Integral("width"));
      	       	     hMassDetcRangeMC->Rebin(reb);
      	       	     hMassDetcRangeMC->Scale(1./hcorr->Integral("width"));
      	       	     hMassPart->Rebin(reb);
      	       	     hMassPart->Scale(1./hcorr->Integral("width"));
      	       	     hMassPartRangeMC->Rebin(reb);
      	       	     hMassPartRangeMC->Scale(1./hcorr->Integral("width"));
     	       	  }
      	       	  if(TMath::AreEqualRel(hmass->GetBinWidth(3), hcorr->GetBinWidth(3), 1.E-12)){
      	       	     Printf("Rebinned hcorr, now all fine");
      	       	  }
      	       }
      	    }
      	    
      	    TH1D *hmassRangeMC = new TH1D(Form("%sRangeMC",hdataname.Data()), Form("%s ; Mass; ",hmass->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1));
      	    hmassRangeMC->SetMarkerStyle(hmass->GetMarkerStyle());
      	    hmassRangeMC->SetMarkerColor(hmass->GetMarkerColor());
      	    hmassRangeMC->SetLineColor(hmass->GetMarkerColor());
      	    hmassRangeMC->Sumw2();
      	    if(iptj==0 && i==0 && j==0) legDistr->AddEntry(hmassRangeMC, "Uncorrected", "P");
      	    TH1D *hmassCorrected = new TH1D(Form("%sCor",hdataname.Data()), Form("%s ; Mass; ", hmass->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1)); //(TH1D*) hmass->Clone(Form("%sCor",hdataname.Data()));
      	    hmassCorrected->SetMarkerStyle(21);
      	    hmassCorrected->SetMarkerColor(colors[i]);
      	    
      	    TH1D *hStatUncMassCorrected = new TH1D(Form("%sCorStatUnc",hdataname.Data()), Form("%s ; Mass; Relative Stat Uncertainty", hmass->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1));
      	    hStatUncMassCorrected->SetMarkerStyle(24);
      	    hStatUncMassCorrected->SetMarkerColor(colors[i]);
      	    hStatUncMassCorrected->SetLineColor(colors[i]);
      	    TH1D *hStatUncExactFactMassCorrected = new TH1D(Form("%sCorStatUncExactFact",hdataname.Data()), Form("%s (Exact Correction Factor); Mass; Relative Stat Uncertainty", hmass->GetTitle()), nbinCorr, minCorr, hcorr->GetBinLowEdge(nbinCorr+1));
      	    hStatUncExactFactMassCorrected->SetMarkerStyle(25);
      	    hStatUncExactFactMassCorrected->SetMarkerColor(colors[i]);
      	    hStatUncExactFactMassCorrected->SetLineColor(colors[i]);
      	    
      	    TH1D *hmassClosureRangeMC=0x0;
      	    
      	    Int_t binminmass = hmass->FindBin(minCorr);
      	    if(nbinCorr>hmass->GetNbinsX()){
      	    Printf("PROBLEEEMMM");
      	    }
      	    for(Int_t k=0; k<nbinCorr;k++){
      	       Double_t content = hmass->GetBinContent(binminmass + k);
      	       Double_t factor =  hcorr->GetBinContent(k+1);
      	       
      	       Double_t error = hmass->GetBinError(binminmass + k);
      	       Double_t errorFactor = hcorr->GetBinError(k+1); 
      	       Double_t errorRec = TMath::Sqrt((error*error* factor*factor)+(content*content* errorFactor*errorFactor));
      	       
      	       Double_t corrCont = content * factor;
      	       hmassCorrected->SetBinContent(k+1, corrCont);
      	       hmassCorrected->SetBinError(k+1, errorRec);
      	       if(corrCont > 0.){
      	       	  hStatUncMassCorrected->SetBinContent(k+1, errorRec/corrCont);
      	       	  hStatUncMassCorrected->SetBinError(k+1, 0);
      	       	  hStatUncExactFactMassCorrected->SetBinContent(k+1, error/corrCont);
      	       	  hStatUncExactFactMassCorrected->SetBinError(k+1, 0);
      	       }
      	       hmassRangeMC->SetBinContent(k+1, content);
      	       hmassRangeMC->SetBinError(k+1, error);
      	       if(bRerangeMC){
      	       	  hMassDetcRangeMC->SetBinContent(k+1, hMassDetc->GetBinContent(k+1));
      	       	  hMassDetcRangeMC->SetBinError(k+1, hMassDetc->GetBinError(k+1));
      	       	  hMassPartRangeMC->SetBinContent(k+1, hMassPart->GetBinContent(k+1));
      	       	  hMassPartRangeMC->SetBinError(k+1, hMassPart->GetBinError(k+1));
      	       	  bRerangeMC=kTRUE;
      	       }
      	       //Printf("hMass Centre %.2f, content %f", hmass->GetBinCenter(binminmass + k), hmass->GetBinContent(binminmass + k));
      	       //Printf("hmassRangeMC Centre %.2f, content %f", hmassRangeMC->GetBinCenter(k+1), hmassRangeMC->GetBinContent(k+1));
      	       
      	    }
      	    
      	    hMassDetcRangeMC->Scale(1./hMassDetcRangeMC->Integral("width"));
      	    hMassPartRangeMC->Scale(1./hMassPartRangeMC->Integral("width"));
      	    
      	    Printf("Integral mass CORRECTED = %f, entries %f", hmassCorrected->Integral(), hmassCorrected->GetEntries());
      	    //hmassCorrected->Rebin(2);
      	    hmassCorrected->Scale(1./hmassCorrected->Integral("width"));
      	    //hmassRangeMC->Rebin(2);
      	    hmassRangeMC->Scale(1./hmassRangeMC->Integral("width"));
      	    
      	    hmassClosureRangeMC = (TH1D*)hmassCorrected->Clone("hmassClosureRangeMC");
      	    hmassClosureRangeMC->SetMarkerStyle(25);
      	    
      	    hmassClosureRangeMC->Divide(hcorr);
      	    hmassClosureRangeMC->Scale(1./hmassClosureRangeMC->Integral("width"));
      	    
      	    cMassCheck[i][j]->cd(iptj+1);
      	    hmassRangeMC->SetMarkerStyle(27);
      	    hmassRangeMC->DrawClone("Psames");
      	    pvpT->Draw();
      	    
      	    TH1D* hmassRDetcMC = (TH1D*)hmassRangeMC->Clone(Form("hmassRDetcMC%d%d",i,j));
      	    hmassRDetcMC->SetMarkerColor(hMassDetc->GetMarkerColor());
      	    hmassRDetcMC->SetMarkerStyle(hMassDetc->GetMarkerStyle());
      	    hmassRDetcMC->SetTitle(";Mass ; Data / MC");
      	    hmassRDetcMC->Divide(hMassDetcRangeMC);
      	    
      	    TH1D* hmassCorrectedRPartMC = (TH1D*) hmassCorrected->Clone(Form("hmassCorrectedRPartMC%d%d",i,j));
      	    hmassCorrectedRPartMC->SetTitle(";Mass ; Data / MC");
      	    hmassCorrectedRPartMC->SetMarkerColor(hMassPart->GetMarkerColor());
      	    hmassCorrectedRPartMC->SetMarkerStyle(hMassPart->GetMarkerStyle());
      	    Printf("Division masscorrected (n = %d, w = %.1f) /mass part (n = %d, w = %.1f)", hmassCorrected->GetNbinsX(), hmassCorrected->GetBinWidth(3),hMassPartRangeMC->GetNbinsX(), hMassPartRangeMC->GetBinWidth(3));
      	    hmassCorrectedRPartMC->Divide(hMassPartRangeMC);
      	    
      	    TH1D* hmassCorrectedRUnfold =  (TH1D*)hmassCorrected->Clone(Form("hmassCorrectedRUnfold%d%d",i,j));
      	    
      	    if(hMassDistrUnf) {
      	       //hMassDistrUnf->Scale(1./hMassDistrUnf->Integral("width"));
      	       //hmassCorrectedRUnfold->SetMarkerColor(gMassDistrUnf->GetMarkerColor());
      	       //hmassCorrectedRUnfold->SetMarkerStyle(gMassDistrUnf->GetMarkerStyle());
      	       //Printf("Division masscorrected (n = %d, w = %.1f) /mass unfolded (n = %d, w = %.1f)", hmassCorrected->GetNbinsX(), hmassCorrected->GetBinWidth(3) ,hMassDistrUnf->GetNbinsX(), hMassDistrUnf->GetBinWidth(3));
      	       //hmassCorrectedRUnfold->Divide(hMassDistrUnf);
      	    }
      	    Printf("Divisions done");
      	    //hmassCorrected->Scale(1./hmassCorrected->Integral());
      	    //hmassReUncorrected->Divide(hmassCorrected);
      	    //hmassCorrected->Scale(1./hmassCorrected->Integral());
      	    if(iptj==0 && i==0 && j==0) legDistrCor->AddEntry(hmassCorrected, "corrected", "P");
      	    
      	    cMassCorr[i][j]->cd(iptj+1);
      	    //gPad->SetLogy();
      	    hMassPartRangeMC->GetYaxis()->SetRangeUser(0.,0.4);
      	    hMassPartRangeMC->GetXaxis()->SetRangeUser(0.,20);
      	    hMassPartRangeMC->Draw("Psames");
      	    if(gMassDistrUnf) gMassDistrUnf->Draw("P");
      	    hmassCorrected->Draw("Psames");
      	    legDistrCor->Draw();
      	    pvpT->Draw();
      	    
      	    fres->cd();
      	    hmassCorrected->Write();
      	    
      	    cMass[i][j]->cd(iptj+1);
      	    hMassDetcRangeMC->GetXaxis()->SetRangeUser(0.,20);
      	    hMassDetcRangeMC->GetYaxis()->SetRangeUser(0.,0.4);
      	    hMassDetcRangeMC->Draw("Psames");
      	    hmassRangeMC->Draw("Psames");
      	    legDistr->Draw();
      	    pvpT->Draw();
      	    
      	    cMassRatios[i][j]->cd(iptj+1);
      	    TLegend *legR = new TLegend(0.5, 0.5, 0.9, 0.7);
      	    legR->SetFillStyle(0);
      	    legR->SetBorderSize(0);
      	    legR->AddEntry(hmassCorrectedRPartMC, "Corrected", "P");
      	    legR->AddEntry(hmassRDetcMC, "Uncorrected", "P");
      	    hmassCorrectedRPartMC->GetYaxis()->SetRangeUser(0.,4.);
      	    //hmassCorrectedRUnfold->Draw("Psames");
      	    hmassCorrectedRPartMC->Draw("Psames");
      	    hmassRDetcMC->Draw("Psames");
      	    pvpT->Draw();
      	    legR->Draw();
      	    
      	    cMassClosure[i][j]->cd(iptj+1);
      	    hmassClosureRangeMC->Draw("samesP");
      	    hmassRangeMC->Draw("samesP");
      	    pvpT->Draw();
      	    
      	    cMassRelUnc[j]->cd(iptj+1);
      	    gPad->SetGridx();
      	    gPad->SetGridy();
      	    hStatUncMassCorrected->GetYaxis()->SetRangeUser(0.,0.4);
      	    hStatUncMassCorrected->Draw("samesP");
      	    hStatUncExactFactMassCorrected->Draw("samesP");
      	    pvpT->Draw();
      	    
      	    if(j==0 && iptj==0){
      	       if(i==0) legMethMC->AddEntry(hMassPartRangeMC, "MC, particle level", "P");
      	       legMethMC->AddEntry(hmassCorrected, Form("%s", methodname[i].Data()), "P");
      	    }
      	    cMassCorrMethMC[j]->cd(iptj+1);
      	    hMassPartRangeMC->GetXaxis()->SetRangeUser(0.,20);
      	    hMassPartRangeMC->GetYaxis()->SetRangeUser(0.,0.4);
      	    hMassPartRangeMC->Draw("Psames");
      	    hmassCorrected->Draw("Psames");
      	    legMethMC->Draw();
      	    pvpT->Draw();
      	    
      	    SaveCv(cMassCorr[i][j]);
      	    SaveCv(cMass[i][j]);
      	    SaveCv(cMassRatios[i][j]);
      	    SaveCv(cMassCheck[i][j]);
      	    SaveCv(cMassClosure[i][j]);
      	    
      	    
      	    gMeanMass[i][j]->SetMarkerColor(hmass->GetMarkerColor());
      	    gMeanMass[i][j]->SetLineColor(hmass->GetMarkerColor());
      	    gMeanMass[i][j]->SetPoint(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5 + pTjet[iptj], hmass->GetMean());
      	    gMeanMass[i][j]->SetPointError(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5,hmass->GetMeanError());
      	    
      	    gMeanMassCor[i][j]->SetMarkerColor(hmass->GetMarkerColor());
      	    gMeanMassCor[i][j]->SetLineColor(hmass->GetMarkerColor());
      	    gMeanMassCor[i][j]->SetPoint(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5 + pTjet[iptj],hmassCorrected->GetMean());
      	    gMeanMassCor[i][j]->SetPointError(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5,hmassCorrected->GetMeanError());

      	    gMedianMass[i][j]->SetMarkerColor(hmass->GetMarkerColor());
      	    gMedianMass[i][j]->SetLineColor(hmass->GetMarkerColor());
      	    gMedianMass[i][j]->SetPoint(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5 + pTjet[iptj], GetQuantile(hmass));
      	    gMedianMass[i][j]->SetPointError(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5,0);
      	    
      	    gMedianMassCor[i][j]->SetMarkerColor(hmass->GetMarkerColor());
      	    gMedianMassCor[i][j]->SetLineColor(hmass->GetMarkerColor());
      	    gMedianMassCor[i][j]->SetPoint(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5 + pTjet[iptj],GetQuantile(hmassCorrected));
      	    gMedianMassCor[i][j]->SetPointError(iptj,(pTjet[iptj+1] - pTjet[iptj])*0.5,0);
      	 }
      }
   }
   for(Int_t j=0; j<ntags; j++ ){
      SaveCv(cMassRelUnc[j]);
      SaveCv(cMassCorrMethMC[j]);
   }
   SaveCv(cCorrection);
   
   TLegend *legCmp = new TLegend(0.6, 0.2, 0.9, 0.5);
   legCmp->SetBorderSize(0);
   legCmp->SetFillStyle(0);
   TLegend *legCmp2 = (TLegend*)legCmp->Clone();
   cMeanMassCorrMethodCmp->cd();
   if(gMeanMassUnf && gMeanMassUnfSyst){
      gMeanMassUnfSyst->SetFillColor(kOrange);
      gMeanMassUnf->SetLineColor(kOrange+7);
      gMeanMassUnf->SetMarkerColor(kOrange+7);
      gMeanMassUnfSyst->Draw("A2");
      gMeanMassUnf->Draw("P"); 
      legCmp->AddEntry(gMeanMassUnf, "From Unfolding", "PL");
   }
   legCmp->AddEntry(gMeanMassPMC, "LHC13b4_plus Part", "PL");
   if(gMeanMassPMC) gMeanMassPMC->Draw("P");
   legCmp->Draw();

   cMeanMassUncoMethodCmp->cd();
   legCmp2->AddEntry(gMeanMassDMC, "LHC13b4_plus Detc", "PL");
   if(gMeanMassDMC) gMeanMassDMC->Draw("AP");

   for(Int_t i = 0; i < nmethods; i++){
      for(Int_t j = 0; j< ntags; j++){
      	 cMeanMass[i][j]->cd();
      	 gMeanMass[i][j]->Draw("AP");
      	 gMeanMassCor[i][j]->Draw("P");
      	 gMedianMass[i][j]->Draw("P");
      	 gMedianMassCor[i][j]->Draw("P");
      	 legCorrUncorr->Draw();
      	 SaveCv(cMeanMass[i][j]);
      }
      cMeanMassUncoMethodCmp->cd();
      legCmp2->AddEntry(gMeanMass[i][1], Form("%s", methodname[i].Data()), "PL");
      //if(i==0) gMeanMass[i][1]->Draw("P");
      //else 
      gMeanMass[i][1]->Draw("P");
      //gMedianMass[i][1]->Draw("P");
      
      cMeanMassCorrMethodCmp->cd();
      if(i==1) continue; //do not draw the non bkg-sub mass (Area-Raw) 
      legCmp->AddEntry(gMeanMassCor[i][1], Form("%s", methodname[i].Data()), "PL");
      //if(i==0) gMeanMassCor[i][1]->Draw("AP");
      //else 
      gMeanMassCor[i][1]->Draw("P");
      //gMedianMassCor[i][1]->Draw("P");
      fres->cd();
      gMeanMassCor[i][1]->Write();
      
   }
   cMeanMassUncoMethodCmp->cd();
   legCmp2->Draw();
   SaveCv(cMeanMassCorrMethodCmp);
   SaveCv(cMeanMassUncoMethodCmp);

}
//-------------------------------------------------------------------------------------------

void SuperimposeFixedInput(){
   SetStyle(0);
   const Int_t nfiles = 2;
   TString inputF1[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output794-795/resultspT20-120/Results.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output794-795/Systematics/TrackingEff/ResultsTrEff.root"};
   
   TString legs1[nfiles] = {"Default", "Tracking Efficiency"};
   
   Superimpose(nfiles, inputF1, legs1, 2, "Tagged");
}

void Superimpose(Int_t nfiles, TString paths[], TString legs[], Int_t methodIndex, TString tag){
   
   TH1D* histograms[nfiles][npTjbins];
   
   TString namebase = "hMass";
   
   for(Int_t ifile = 0; ifile < nfiles; ifile++){
   
      TFile *fin = new TFile(paths[ifile]);
      if(!fin->IsOpen()){
      	 Printf("File %s not found ", paths[ifile].Data());
      	 continue;
      }
      fin->ls();
      for(Int_t iptj = 0; iptj<npTjbins; iptj++){
      	 TString hname = Form("%s%s%dpTj%dCor", namebase.Data(), tag.Data(), methodIndex, iptj);
      	 
      	 histograms[ifile][iptj] = 0x0;
      	 
      	 histograms[ifile][iptj] = (TH1D*)fin->Get(hname.Data());
      	 
      	 if(!histograms[ifile][iptj]){
      	    Printf("%s not found ", hname.Data());
      	    continue;
      	 }
      	 histograms[ifile][iptj]->SetName(Form("%sV%d", hname.Data(), ifile));
      	 histograms[ifile][iptj]->SetLineColor(colors[ifile]);
      	 histograms[ifile][iptj]->SetMarkerColor(histograms[ifile][iptj]->GetLineColor());
      	 
      }
   }
   Int_t nx, ny, dx, dy;
   CalculatePads(npTjbins, nx, ny, dx, dy);
   TCanvas *cMassvspTj = new TCanvas("cMassvspTj", "Comparison plot",  dx, dy);
   cMassvspTj->Divide(nx, ny);
   TCanvas *cMassRatiosvspTj = new TCanvas("cMassRatiosvspTj", "Ratio plot",  dx, dy);
   cMassRatiosvspTj->Divide(nx, ny);
    
   
   TLegend *legend = new TLegend(0.6, 0.5, 0.9, 0.8);
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   
   for(Int_t iptj = 0; iptj<npTjbins; iptj++){
      for(Int_t ifile = 0; ifile < nfiles-1; ifile++){
      	 if(!histograms[ifile][iptj] || !histograms[ifile+1][iptj]){
      	    Printf("Histogram missing, comparison not possible");
      	    continue;
      	 }
      	 if(iptj == 0){
      	    legend->AddEntry(histograms[ifile][iptj], legs[ifile], "P");
      	    legend->AddEntry(histograms[ifile+1][iptj], legs[ifile+1], "P");
      	 }
      	 TH1 *hRatio = CompareDistributions(cMassvspTj->cd(iptj+1), histograms[ifile][iptj], histograms[ifile+1][iptj]);
      	 legend->Draw();
      	 
      	 cMassRatiosvspTj->cd(iptj+1);
      	 hRatio->Draw();
      	 
      }
   }
}
//-------------------------------------------------------------------------------------------

void MassVsCentralityMB(){
   Int_t n = 3;
   TString path = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output834-835/";
   TString tags[n] = {"_TCRaw", "_TCDeriv", "ConstSub_TC"};
   
   for(Int_t i=0; i<n; i++){
      MassVsCentrality(tags[i], path);
   
   }


}

void MassVsCentralityJ1(){
   Int_t n = 3;
   TString path = "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output836-837/";
   TString tags[n] = {"_TCRawJ1", "_TCDerivJ1", "ConstSub_TCConstJ1"};
   
   for(Int_t i=0; i<n; i++){
      MassVsCentrality(tags[i], path);
   
   }


}

void MassVsCentrality(TString tag /*e.g. "_TCDeriv"*/, TString path){

   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   Bool_t bEmb = kFALSE;

   const Int_t nCentBins = 5;
   Double_t centBins[nCentBins+1] = {0,10,20,40,60,100};
   PlotUtilEmcalJetMass* util = new PlotUtilEmcalJetMass();
   util->SetInputFileName(Form("%s/AnalysisResults.root", path.Data()));
   util->SetJetRadius(R);
   util->SetJetType(type);
   
   //util->SetJetPtRange(ptmin,ptmax);
   util->LoadFile();
   
   util->SetTag(tag); //"ConstSub_TC"
   util->SetConstTag("");
   if(util->LoadList()) {
      Printf("Found! all good");
   } else {
      Printf("List not found");
      return;
   }
   Int_t nx, ny, dx, dy;
   CalculatePads(npTjbins, nx, ny, dx, dy);
   
   TPaveText *pvfixedpt = new TPaveText(0.5, 0.3, 0.85, 0.5, "NDC");
   pvfixedpt->SetFillStyle(0);
   pvfixedpt->SetBorderSize(0);
   pvfixedpt->AddText(Form("%.0f < p_{T,jet} < %.0f GeV/#it{c}", pTjet[0],pTjet[npTjbins]));
   
   TPaveText *pvpT[npTjbins];
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      pvpT[iptj]= new TPaveText(0.5, 0.3, 0.85, 0.5, "NDC");
      pvpT[iptj]->SetFillStyle(0);
      pvpT[iptj]->SetBorderSize(0);
      pvpT[iptj]->AddText(Form("%.0f < p_{T,jet} < %.0f GeV/#it{c}", pTjet[iptj],pTjet[iptj+1]));
   }
   
   TFile *fout = new TFile(Form("MassOutputvsCentrality%s.root", tag.Data()),"recreate");
   gDirectory->cd(0);
   TCanvas *cmassCP = new TCanvas(Form("cmassCP%s", tag.Data()), Form("Central/peripheral mass ratio (%s)", tag.Data()));
   
   TCanvas *cmassCPpT = new TCanvas(Form("cmassCPpT%s", tag.Data()), Form("Central/peripheral mass ratio (%s)", tag.Data()), dx, dy);
   cmassCPpT->Divide(nx, ny);
   TH1F *hmCP;
   TH1F *hmCPpT[npTjbins];
   
   for(Int_t ic = 0; ic < nCentBins; ic++){
      TCanvas *cmass = new TCanvas(Form("cmassC%s%d", tag.Data(), ic), Form("Mass %s in centrality %.0f - %.0f perc", tag.Data(), centBins[ic], centBins[ic+1]));
      
      TCanvas *cmasspT = new TCanvas(Form("cmasspTC%s%d", tag.Data(), ic), Form("pT bins, Mass %s in centrality %.0f - %.0f perc", tag.Data(), centBins[ic], centBins[ic+1]), dx, dy);
      cmasspT->Divide(nx, ny);

      TPaveText* pvCent =  new TPaveText(0.5, 0.7, 0.75, 0.8, "NDC");
      pvCent->SetFillStyle(0);
      pvCent->SetBorderSize(0);
      pvCent->AddText(Form("%.0f-%.0f%s", centBins[ic], centBins[ic+1], "%") );
      util->SetJetPtRange(pTjet[0], pTjet[npTjbins]);
      TH1D* hM = util->GetMassInCentrality(centBins[ic], centBins[ic+1], PlotUtilEmcalJetMass::kJetTagged, 0); //iList = 0, only one list present
      cmass->cd();
      hM->Sumw2();
      hM->Draw("PE");
      pvfixedpt->Draw();
      pvCent->Draw();
      
      SaveCv(cmass);
      fout->cd();
      hM->Write();
      if(ic==0) hmCP = (TH1F*) hM->Clone("hmCP");
      if(ic==nCentBins-1) hmCP->Divide(hM);
      
      for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      	 util->SetJetPtRange(pTjet[iptj], pTjet[iptj+1]);
      	 TH1D* hMpT = util->GetMassInCentrality(centBins[ic], centBins[ic+1], PlotUtilEmcalJetMass::kJetTagged, 0); //iList = 0, only one list present
      	 //hMpT->SetName(Form("hMassCent%d_pTj%d", ic, iptj));
      	 cmasspT->cd(iptj+1);
      	 hMpT->Sumw2();
      	 hMpT->Draw("PE");
      	 pvpT[iptj]->Draw();
      	 pvCent->Draw();
      	 SaveCv(cmasspT);
      	 fout->cd();
      	 hMpT->Write();
      	 
      	 if(ic==0) hmCPpT[iptj] = (TH1F*) hMpT->Clone(Form("hMassCP_pTj%d", iptj));
      	 if(ic==nCentBins-1) hmCPpT[iptj]->Divide(hMpT);

      }
      
   }
   cmassCP->cd();
   cmassCP->SetGridy();
   hmCP->GetYaxis()->SetRangeUser(0,3.5);
   hmCP->GetXaxis()->SetRangeUser(0,20);
   hmCP->Draw("PE");
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      cmassCPpT->cd(iptj+1);
      hmCPpT[iptj]->GetYaxis()->SetRangeUser(0,3.5);
      hmCPpT[iptj]->GetXaxis()->SetRangeUser(0,20);

      hmCPpT[iptj]->Draw("PE");
      pvpT[iptj]->Draw();
   }
   
   SaveCv(cmassCP);
   SaveCv(cmassCPpT);
}
//-------------------------------------------------------------------------------------------

void RunSumHistogramstoMassOutput(){
   Int_t nmethods = 3;
   //
   TString tags[nmethods] = {"_TCRawJ1", "_TCDerivJ1", "ConstSub_TCConstJ1"};
   
   const Int_t nfiles = 2;
   TString path[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output834-835/analysis/MassvsCent/", "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output836-837/analysis/MvsCent/"};
   
   TString basefilename = "MassOutputvsCentrality";
   TString meaningfulflags[nfiles] = {"MB", "J1"};
   Int_t MBEJETrigComb = 1;
   Int_t centBin = 4; //see the bins in MassVsCentrality and SumHistogramstoMassOutput

   for(Int_t j = 0; j<nfiles; j++){
      
      for(Int_t i=0; i<nmethods; i++){
      	 TString filenames[nfiles] = {Form("%s%s%s",path[0].Data(), tags[i].Data(), basefilename.Data()), Form("%s%s%s",path[1].Data(), tags[i].Data(), basefilename.Data())};
      	 
      	 SumHistogramstoMassOutput(nfiles, filenames, meaningfulflags, MBEJETrigComb, centBin); 
      	 
      }
      
   }
   
}

void SumHistogramstoMassOutput(Int_t nfiles, TString filenames[], TString meaningfulflags[], Int_t MBEJETrigComb /*0 = only MB, 1 = MB for pT < 80, EJE pT>80, 2 = sum MB and EJE for pT> 80, 3 = include J1 */, Int_t centBin){
   
   TString namehbase = "fh1MassTaggedLst0_";
   TH1F *hMasspT[npTjbins];
   Int_t high = 0, low = 1;
   Int_t ntrig = 2;
   TString trigname[ntrig] = {"J1", "J2"} ; //for pPb 2013: J1 high threshold, J2 low threshold
   Double_t thresh[ntrig]; thresh[high] = 79; thresh[low] = 59;
   const Int_t nCentBins = 5;
   Double_t centBins[nCentBins+1] = {0,10,20,40,60,100};
   
   Int_t nx, ny, dx, dy;
   CalculatePads(npTjbins, nx, ny, dx, dy);

   TCanvas *cMasspT = new TCanvas(Form("cMasspT_Cent%.0f-%.0f", centBins[centBin], centBins[centBin+1]), Form("Jet Mass in pT bins Cent%.0f-%.0f", centBins[centBin], centBins[centBin+1]), dx, dy);
   cMasspT->Divide(nx, ny);
   
   for(Int_t iptj = 0; iptj<npTjbins; iptj++){
      
      hMasspT[iptj] = 0;
      
      for(Int_t i=0; i< nfiles; i++){ //files
      	 TFile* file = new TFile(filenames[i]);
      	 if(!file->IsOpen()){
      	    Printf("File %s not found, skipping", filenames[i].Data());
      	    continue;
      	 }
      	 TString hname = Form("%sCent%.0f-%.0f_pTj%.0f-%.0f", namehbase.Data(), centBins[centBin], centBins[centBin+1], pTjet[iptj],pTjet[iptj+1]);
      	 TH1F* hMasspTtmp = (TH1F*)file->Get(hname);
      	 if(!hMasspTtmp){
      	    Printf("%s not found, skip", hname.Data());
      	    continue;
      	    
      	 }
      	 if(pTjet[iptj]<thresh[high]){
      	    
      	    if(meaningfulflags[i].Contains("MB")) {
      	       if(!hMasspT[iptj]){
      	       	  hMasspT[iptj] = (TH1F*)hMasspTtmp->Clone(Form("%sSum", hMasspTtmp->GetName()));
      	       	  
      	       } else Printf("Should not enter here");
      	       
      	    }
      	    
      	 }
      	 if(pTjet[iptj]>thresh[high]){
      	    if(meaningfulflags[i].Contains(trigname[high])) {
      	       if(MBEJETrigComb == 1){
      	       	  hMasspT[iptj] = (TH1F*)hMasspTtmp->Clone(Form("%sSum", hMasspTtmp->GetName()));
      	       } else { if(MBEJETrigComb>1) //2,3
      	       	  hMasspT[iptj]->Add(hMasspTtmp);
      	       }
      	       
      	    }
      	    
      	 }
      	 if(pTjet[iptj]>thresh[low]){
      	    if(meaningfulflags[i].Contains(trigname[low])) {
      	       if(MBEJETrigComb == 3){ 
      	       	  hMasspT[iptj]->Add(hMasspTtmp);
      	       }
      	    }
      	 }
      } //files
      hMasspT[iptj]->Sumw2();
      cMasspT->cd(iptj+1);
      hMasspT[iptj]->Draw();
      
   } //pT
}

//-------------------------------------------------------------------------------------------
void plotJetMasspPbSimple(Int_t ntag, TString turnontags[], Bool_t bEmb, TString pathfile){

   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   
   if (bEmb) type = "ThrmTracksEmb";
   
   Printf("READING... %s", pathfile.Data());
   
   PlotUtilEmcalJetMass *util=0x0;
   TString pathfilearr[1] = {pathfile};
   InitUtilFiles(util, ntag, turnontags, 1 , pathfilearr, R, pTjet[0], pTjet[npTjbins], type, centBin, bEmb);
   util->SetMinLeadTrackPt(0);
   TPaveText *pvpT[npTjbins];
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      pvpT[iptj]= new TPaveText(0.5, 0.3, 0.85, 0.5, "NDC");
      pvpT[iptj]->SetFillStyle(0);
      pvpT[iptj]->SetBorderSize(0);
      pvpT[iptj]->AddText(Form("%.0f < p_{T,jet} < %.0f GeV/#it{c}", pTjet[iptj],pTjet[iptj+1]));
   }
   
   Int_t nx, ny, dx, dy;
   CalculatePads(npTjbins, nx, ny, dx, dy);
   //TCanvas *cmassptjbin = new TCanvas(Form("cmassptjbin%d", iptj), Form("Mass for pTj %.1f - %.1f", pTjet[iptj], pTjet[iptj+1]), dx, dy);
   TCanvas *cmassptjbin[ntag];
   for(Int_t itag = 0; itag< ntag; itag++){
      cmassptjbin[itag] = new TCanvas(Form("cmassptjbin%d", itag), Form("Mass %s", turnontags[itag].Data()), dx, dy);
      cmassptjbin[itag]->Divide(nx, ny);
   }
   
   for(Int_t iptj=0 ; iptj < npTjbins; iptj++){ util->SetJetPtRange(pTjet[iptj], pTjet[iptj+1]);
      for(Int_t itag = 0; itag< ntag; itag++){
      	 TH1D *hMass = util-> GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,itag);
      	 cmassptjbin[itag]->cd(iptj+1);
      	 hMass->SetName(Form("hMass%s_pT%d", turnontags[itag].Data(), iptj));
      	 hMass->GetXaxis()->SetRangeUser(-5, 20);
      	 hMass->Draw();
      	 pvpT[iptj]->Draw();
      }
   }
   for(Int_t itag = 0; itag< ntag; itag++){
      SaveCv(cmassptjbin[itag]);
   }
}
//-------------------------------------------------------------------------------------------
void plotJetMasspPbSimpleEmbedding(){
   Int_t ntag = 1;
   TString turnontags[ntag] = {"_TC"};
   Bool_t bEmb = kTRUE;
   TString pathfile = "/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/merge/weights/analysis/AnalysisResults.root";
   plotJetMasspPbSimple(ntag, turnontags, bEmb, pathfile);
   
}
//-------------------------------------------------------------------------------------------

void RhopPbvsPYTHIAFixedPath(){
   //used for data-MC comparison or comparison of one specific bin of the toy model (change below)
   TString tracks[2] = {"PicoTracks", "PicoTracks"};
   Int_t listnametype[2] = {0,0};
   TString data = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output810-811/AnalysisResults.root";
   TString mc = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train833/outputpTHardBins/mergeRuns/merge1to9/AnalysisResults.root";
   Bool_t toy = kTRUE;
   if(toy){
      mc = //"/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/merge/factor1/analysis/AnalysisResults.root";
      //"/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/merge/weights/AnalysisResults.root";
      "/data/Work/jets/JetMass/BkgFluctStudies/EmbedSinglePart/ExcludeOverlap/JetWithSingleTrack/merge/weights/analysis/AnalysisResults.root";
      tracks[1] = "ThrmTracksEmb";
      listnametype[1] = 1;
   }
   RhopPbvsPYTHIA(data, mc, tracks, listnametype);
}

void RhoMany(){
   
   // used for toy model in bins of Nch
   
   Int_t nbins = 6;
   TString basedir = "/data/Work/jets/JetMass/BkgFluctStudies/EmbedpTDistr/MultipBins/";
   TString paths[nbins];
   TString lists[nbins];
   TString trkname[nbins];
   Int_t listnametype[nbins];
   
   for(Int_t i = 0; i < nbins; i++){
      paths[i] = Form("%s%02d/output/AnalysisResults.root", basedir.Data(), i+1);
      trkname[i] = "ThrmTracksEmb";
      listnametype[i] = 1;
   
   }
   
   RhopPbvsPYTHIA(nbins, paths, trkname, listnametype, kTRUE);
   
}

void RhopPbvsPYTHIA(TString dataPath, TString MCpath, TString trkname[], Int_t listnametype[]){
   Int_t ninput =2; 
   TString inputF[ninput]={dataPath, MCpath};
   RhopPbvsPYTHIA(ninput, inputF, trkname, listnametype);
}

void RhopPbvsPYTHIA(Int_t ninput, TString inputF[], TString trkname[], Int_t listnametype[], Bool_t cent){
   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   Bool_t bEmb = kFALSE;
   TString tag = "";
   
   Bool_t manycolorstyle = kFALSE;
   if(ninput>2) manycolorstyle = kTRUE;
   
   Int_t nCentBins = 5;
   Double_t centBins[nCentBins+1] = {0,10,20,40,60,100};
   if(!cent) nCentBins = 1;
   PlotUtilEmcalJetMass* util[ninput];
   TGraphErrors *gMeanRhovsC[ninput];
   TGraphErrors *gMeanRhoMvsC[ninput]; 
   TGraphErrors *gSigmaRhovsC[ninput];
   TGraphErrors *gSigmaRhoMvsC[ninput];
   
   TGraphErrors *gMeanRhoIncl[ninput];
   TGraphErrors *gSigmaRhoIncl[ninput] ;
   TGraphErrors *gMeanRhoMIncl[ninput] ;
   TGraphErrors *gSigmaRhoMIncl[ninput];
   
   TH1D *hRhoIncl[ninput] ;
   TH1D *hRhoMIncl[ninput];
   TH1D *hRhoAll = 0x0; //sum of all
   TH1D *hRhoMAll = 0x0; //sum of all
   
   TCanvas *crhoIncl = new TCanvas(Form("crhoIncl"), Form("Rho inclusive, no centrality"), 1000, 1000);
   crhoIncl->Divide(2,2);
   //TCanvas *crhoAllIncl = new TCanvas(Form("crhoAllIncl"), Form("Rho inclusive, no centrality"), 1000, 1000);
   //crhoAllIncl->Divide(1,2);
   
   TCanvas *crhovsCent[ninput];
   //save rho and rho_m inclusive
   TFile *foutRhos = new TFile("RhoDistributions.root", "recreate");
   
   for(Int_t f = 0; f<ninput; f++){
      
      util[f] = new PlotUtilEmcalJetMass();
      util[f]->SetInputFileName(inputF[f]);
      util[f]->SetJetRadius(R);
      util[f]->SetJetType(type);
      
      //util->SetJetPtRange(ptmin,ptmax);
      util[f]->LoadFile();
      
      util[f]->SetTag(tag); //"ConstSub_TC"
      util[f]->SetConstTag("");
      util[f]->SetListNameType(listnametype[f]);
      
      if(util[f]->LoadRhoLists(trkname[f])) {
      	 Printf("Data found! all good");
      } else {
      	 Printf("List for data not found");
      	 
      	 if(manycolorstyle) continue;
      	 else return;
      }
      if(manycolorstyle) 
      	 crhovsCent[f] = new TCanvas(Form("crhovsCent%d", f), Form("Rho vs Centrality (%d)", f), 1000, 600);
      else crhovsCent[f] = new TCanvas(Form("crhovsCent%s", f==0 ? "DATA" : "MC"), Form("Rho vs Centrality (%s)", f==0 ? "DATA" : "MC"), 1000, 600);
      crhovsCent[f]->Divide(2,1);
      crhovsCent[f]->cd(1);
      gPad->SetLogz();
      TH2D *h2drho = util[f]->GetRho(1);//default vs centrality
      if(h2drho) h2drho->Draw("colz"); 
      else Printf("shit");
      crhovsCent[f]->cd(2);
      gPad->SetLogz();
      TH2D *h2drhom = util[f]->GetRhoM(1);//default vs centrality
      if(h2drhom) h2drhom->Draw("colz"); 
      
      SaveCv(crhovsCent[f]);
      
      gMeanRhovsC[f] = new TGraphErrors(nCentBins);
      if(manycolorstyle){
      	 gMeanRhovsC[f]->SetMarkerStyle(20);
      	 gMeanRhovsC[f]->SetMarkerColor(colors[f]);
      } else{
      	 if(f == 0){
      	    gMeanRhovsC[f]->SetMarkerStyle(20);
      	 } else {
      	    gMeanRhovsC[f]->SetMarkerStyle(24);
      	 }
      }
      if(manycolorstyle) {
      	 gMeanRhovsC[f]->SetName(Form("gMeanRhovsC%d", f));
      	 
      } else gMeanRhovsC[f]->SetName(Form("gMeanRhovsC%s", f == 0 ? "data" : "MC"));
      gMeanRhovsC[f]->SetTitle("Mean #rho; Centrality bin");
      
      gMeanRhoMvsC[f] = (TGraphErrors*)gMeanRhovsC[f]->Clone(Form("gMeanRhoMvsC%s", f == 0 ? "data" : "MC"));
      gMeanRhoMvsC[f]->SetTitle("Mean #rho_{m}; Centrality bin");
      
      gSigmaRhovsC[f] = (TGraphErrors*)gMeanRhovsC[f]->Clone(Form("gSigmaRhovsC%s", f == 0 ? "data" : "MC"));
      gSigmaRhovsC[f]->SetTitle("Sigma #rho; Centrality bin");
       
      gSigmaRhoMvsC[f] = (TGraphErrors*)gMeanRhovsC[f]->Clone(Form("gsigmaRhoMvsC%s", f == 0 ? "data" : "MC"));
      gSigmaRhoMvsC[f]->SetTitle("Sigma #rho_{m}; Centrality bin");
      
      
      
      hRhoIncl[f] = util[f]->GetRhoProjection(-10, 110, 1);
      hRhoMIncl[f] = util[f]->GetRhoMProjection(-10, 110, 1);
      hRhoIncl[f]->Sumw2(); hRhoMIncl[f]->Sumw2();
      if(manycolorstyle) {
      	 hRhoIncl[f]->SetName(Form("hRhoIncl%d", f ));
      	 hRhoMIncl[f]->SetName(Form("hRhoMIncl%d", f ));
      
      } else {
      	 hRhoIncl[f]->SetName(Form("hRhoIncl%s", f == 0 ? "data" : "MC"));
      	 hRhoMIncl[f]->SetName(Form("hRhoMIncl%s", f == 0 ? "data" : "MC"));
      }
      hRhoIncl[f]->SetMarkerColor(colors[nCentBins]);
      hRhoMIncl[f]->SetMarkerColor(colors[nCentBins]);
      hRhoIncl[f]-> SetMarkerColor(colors[f]);
      hRhoMIncl[f]->SetMarkerColor(colors[f]);
      
      if(f==0 || manycolorstyle) {
      	 hRhoIncl[f]->SetMarkerStyle(20);
      	 hRhoMIncl[f]->SetMarkerStyle(20);
      }
      else {
      	 hRhoIncl[f]->SetMarkerStyle(24);
      	 hRhoMIncl[f]->SetMarkerStyle(24);
      	 hRhoIncl[f]-> SetMarkerColor(kOrange - 5);
      	 hRhoMIncl[f]->SetMarkerColor(kOrange - 5);
      }
      
      
      hRhoIncl[f]->Scale(1./hRhoIncl[f]->Integral());
      hRhoMIncl[f]->Scale(1./hRhoMIncl[f]->Integral());
      
      crhoIncl->cd(1);
      gPad->SetLogy();
      hRhoIncl[f]->GetXaxis()->SetRangeUser(0,14);
      hRhoIncl[f]->Draw("Psames");
      crhoIncl->cd(3);
      hRhoMIncl[f]->Draw("Psames");
      
      //write to root file rho and rho_m inclusive per each file
      foutRhos->cd();
      hRhoIncl[f]->Write();
      hRhoMIncl[f]->Write();
      if(manycolorstyle){
      	 if(!hRhoAll) {
      	    hRhoAll = (TH1D*)hRhoIncl[f]->Clone("hRhoAll");
      	    hRhoAll->SetMarkerStyle(25);
      	    hRhoAll->SetMarkerColor(kCyan);
      	    Printf("Init rho sum");
      	 }
      	 hRhoAll->Add(hRhoIncl[f]);
      	 Printf("rho sum, integral now %f", hRhoAll->Integral());
      	 if(!hRhoMAll) {
      	    hRhoMAll = (TH1D*)hRhoMIncl[f]->Clone("hRhoMAll");
      	    hRhoMAll->SetMarkerStyle(25);
      	    hRhoMAll->SetMarkerColor(kCyan);
      	 }
      	 hRhoMAll->Add(hRhoMIncl[f]);
      }
      if(!manycolorstyle && f == 1){
      	 TH1D *hRhoRatioDataoverMCInclus = (TH1D*)hRhoIncl[0]->Clone("hRhoRatioDataoverMCInclus");
      	 hRhoRatioDataoverMCInclus->Divide(hRhoIncl[1]);
      	 hRhoRatioDataoverMCInclus->GetYaxis()->SetRangeUser(0.5, 5);
      	 TH1D *hRhoMRatioDataoverMCInclus = (TH1D*)hRhoMIncl[0]->Clone("hRhoMRatioDataoverMCInclus");
	 hRhoMRatioDataoverMCInclus->Divide(hRhoMIncl[1]);
      	 hRhoMRatioDataoverMCInclus->GetYaxis()->SetRangeUser(0.5, 5);
	 crhoIncl->cd(2);
	 hRhoRatioDataoverMCInclus->GetXaxis()->SetRangeUser(0,14);
	 hRhoRatioDataoverMCInclus->Draw();
	 crhoIncl->cd(4);
	 hRhoMRatioDataoverMCInclus->Draw();
	 
      }
      	 
      gMeanRhoIncl[f] = new TGraphErrors(1);
      gMeanRhoIncl[f]->SetName(Form("gMeanRhoIncl%s", f == 0 ? "data" : "MC" ));
      if(f== 0) gMeanRhoIncl[f]  ->SetMarkerStyle(20);
      else gMeanRhoIncl[f]  ->SetMarkerStyle(24);
      gMeanRhoIncl[f]  ->SetMarkerColor(colors[nCentBins]);
      
      gSigmaRhoIncl[f]   = (TGraphErrors*)gMeanRhoIncl[f]->Clone(Form("gSigmaRhoIncl%s", f == 0 ? "data" : "MC" ));
      gMeanRhoMIncl[f]   = (TGraphErrors*)gMeanRhoIncl[f]->Clone(Form("gMeanRhoMIncl%s", f == 0 ? "data" : "MC") );
      gSigmaRhoMIncl[f]  = (TGraphErrors*)gMeanRhoIncl[f]->Clone(Form("gSigmaRhoMIncl%s", f == 0 ? "data" : "MC"));
      
      
      gMeanRhoIncl[f]   ->SetTitle("Mean #rho");
      gSigmaRhoIncl[f]  ->SetTitle("Sigma #rho");
      gMeanRhoMIncl[f]  ->SetTitle("Mean #rho_{m}");
      gSigmaRhoMIncl[f] ->SetTitle("Sigma #rho_{m}");
      
      gMeanRhoIncl[f]  ->SetPoint(0,2,hRhoIncl[f]->GetMean());
      gMeanRhoIncl[f]  ->SetPointError(0,2,hRhoIncl[f]->GetMeanError());
      gSigmaRhoIncl[f] ->SetPoint(0,2,hRhoIncl[f]->GetRMS());
      gSigmaRhoIncl[f] ->SetPointError(0,2,hRhoIncl[f]->GetRMSError());
      
      gMeanRhoMIncl[f] ->SetPoint(0,2,hRhoMIncl[f]->GetMean());
      gMeanRhoMIncl[f] ->SetPointError(0,2,hRhoMIncl[f]->GetMeanError());
      gSigmaRhoMIncl[f]->SetPoint(0,2,hRhoMIncl[f]->GetRMS());
      gSigmaRhoMIncl[f]->SetPointError(0,2,hRhoMIncl[f]->GetRMSError());
      
   }
   if(manycolorstyle){
      crhoIncl->cd(1);
      Printf ("The histogram %s is here, it has %f entries, %.2f integral, %f/%f minimum/maximum edge", hRhoAll->GetName(), hRhoAll->GetEntries(), hRhoAll->Integral(), hRhoAll->GetBinLowEdge(1), hRhoAll->GetBinLowEdge(hRhoAll->GetNbinsX()));
      hRhoAll->Scale(1./hRhoAll->Integral());
      hRhoAll->Draw("sames");
      crhoIncl->cd(3);
      hRhoMAll->Scale(1./hRhoMAll->Integral());
      hRhoMAll->Draw("sames");
      //write to root file the sum of rho and rho_m inclusive in each file (this is meanignful only for the toy model simulation that uses bins in Nch). It has to be compared with the inclusive distribution from data 
      foutRhos->cd();
      hRhoAll->Write();
      hRhoMAll->Write();
      SaveCv(crhoIncl);
   }
   TCanvas *crho = new TCanvas(Form("crho"), Form("Rho"), 1000, 1000);
   crho->Divide(2,2);
   TCanvas *crho2 = new TCanvas(Form("crho2"), Form("Rho - MC no centrality"), 1000, 1000);
   crho2->Divide(2,2);
   TCanvas *crhoparam = new TCanvas("crhoparam", "Rho mean and sigma", 1000, 1000);
   crhoparam->Divide(2,2);
   
   TLegend *legC = new TLegend(0.5, 0.2, 0.8, 0.7);
   legC->SetBorderSize(0);
   legC->SetFillStyle(0);
    
   for(Int_t ic = 0; ic < nCentBins; ic++){
      TH1D *hRhovsC[ninput];
      TH1D *hRhoMvsC[ninput];
      for(Int_t id = 0; id<ninput; id++){ //loop on the list of inputs
      	 
      	 //rho
      	 if(cent) hRhovsC[id] = util[id]->GetRhoProjection(centBins[ic], centBins[ic+1], 1);
      	 else hRhovsC[id] = util[id]->GetRhoProjection(0, 100, 1);
      	 //scale integral
      	 if(!hRhovsC[id]) continue;
      	 hRhovsC[id]->Scale(1./hRhovsC[id]->Integral());
      	 hRhovsC[id]->SetLineColor(colors[ic]);
      	 hRhovsC[id]->SetMarkerColor(colors[ic]);
      	 hRhovsC[id]->SetLineWidth(2);
      	 if(manycolorstyle){
      	    hRhovsC[id]->SetName(Form("%s%d", hRhovsC[id]->GetName(), id));
      	    hRhovsC[id]->SetMarkerStyle(20);
      	    hRhovsC[id]->SetMarkerColor(colors[id]);
      	    legC->AddEntry(hRhovsC[id], Form("%.0f-%.0f%s", centBins[ic], centBins[ic+1], "%"));
      	 }
      	 else{
      	 if(id==0) {
      	    hRhovsC[id]->SetName(Form("%sdata", hRhovsC[id]->GetName()));
      	    hRhovsC[id]->SetMarkerStyle(20);
      	    legC->AddEntry(hRhovsC[id], Form("%.0f-%.0f%s", centBins[ic], centBins[ic+1], "%"));
      	 }
      	 if(id==1) {
      	    hRhovsC[id]->SetName(Form("%sMC", hRhovsC[id]->GetName()));
      	    hRhovsC[id]->SetMarkerStyle(24); //MC
      	    TH1D *hRhoRatioDataoverMC = (TH1D*)hRhovsC[0]->Clone(Form("hRhoRatioDataoverMC_Cent%.0f-%.0f", centBins[ic], centBins[ic+1]));
      	    hRhoRatioDataoverMC->SetTitle("Ratio #rho Data/MC");
      	    hRhoRatioDataoverMC->Divide(hRhovsC[1]);
      	    crho->cd(2);
      	    hRhoRatioDataoverMC->GetYaxis()->SetRangeUser(-0.1, 5);
      	    hRhoRatioDataoverMC->GetXaxis()->SetRangeUser(0,14);
      	    hRhoRatioDataoverMC->Draw("sames");
      	    
      	    TH1D *hRhoRatioDataoverMCNoC = (TH1D*)hRhovsC[0]->Clone(Form("hRhoRatioDataoverMCNoC_Cent%.0f-%.0f", centBins[ic], centBins[ic+1]));
      	    hRhoRatioDataoverMCNoC->SetTitle("Ratio #rho Data/MC");
      	    hRhoRatioDataoverMCNoC->Divide(hRhoIncl[id]);
      	    crho2->cd(2);
      	    hRhoRatioDataoverMCNoC->GetXaxis()->SetRangeUser(0,14);
      	    hRhoRatioDataoverMCNoC->GetYaxis()->SetRangeUser(-0.1, 5);
      	    hRhoRatioDataoverMCNoC->Draw("sames");
      	    
      	 }
      	 }
      	 crho->cd(1);
      	 gPad->SetLogy();
      	 hRhovsC[id]->GetXaxis()->SetRangeUser(0,14);
      	 hRhovsC[id]->Draw("sames");
      	 
      	 crho2->cd(1);
      	 gPad->SetLogy();
      	 hRhovsC[id]->GetXaxis()->SetRangeUser(0,14);
      	 if(id == 0) hRhovsC[id]->Draw("sames");
      	 else if(ic == 0) hRhoIncl[id]->Draw("sames");
      	 
      	 gMeanRhovsC[id]->SetPoint(ic, ic, hRhovsC[id]->GetMean());
      	 gMeanRhovsC[id]->SetPointError(ic, 0.5, hRhovsC[id]->GetMeanError());
      	 gSigmaRhovsC[id]->SetPoint(ic, ic, hRhovsC[id]->GetRMS());
      	 gSigmaRhovsC[id]->SetPointError(ic, 0.5, hRhovsC[id]->GetRMSError());
      	 
      	 //rhom
      	 hRhoMvsC[id] = util[id]->GetRhoMProjection(centBins[ic], centBins[ic+1], 1);
      	 //scale integral
      	 hRhoMvsC[id]->Scale(1./hRhoMvsC[id]->Integral());
      	 
      	 hRhoMvsC[id]->SetLineColor(colors[ic]);
      	 hRhoMvsC[id]->SetMarkerColor(colors[ic]);
      	 hRhoMvsC[id]->SetLineWidth(2);
      	 if(id==0) {
      	    hRhoMvsC[id]->SetName(Form("%sdata", hRhoMvsC[id]->GetName()));
      	    hRhoMvsC[id]->SetMarkerStyle(20);

      	 }
      	 if(id==1) {
      	    hRhoMvsC[id]->SetName(Form("%sMC", hRhoMvsC[id]->GetName()));
      	    hRhoMvsC[id]->SetMarkerStyle(24); //MC
      	    TH1D *hRhoMRatioDataoverMC = (TH1D*)hRhoMvsC[0]->Clone(Form("hRhoMRatioDataoverMC_Cent%.0f-%.0f", centBins[ic], centBins[ic+1]));
      	    hRhoMRatioDataoverMC->SetTitle("Ratio #rho_{m} Data/MC");
      	    hRhoMRatioDataoverMC->Divide(hRhoMvsC[1]);
      	    crho->cd(4);
      	    hRhoMRatioDataoverMC->GetYaxis()->SetRangeUser(-0.1, 5);
      	    hRhoMRatioDataoverMC->Draw("sames");
      	    
       	    TH1D *hRhoMRatioDataoverMCNoC = (TH1D*)hRhoMvsC[0]->Clone(Form("hRhoMRatioDataoverMCNoC_Cent%.0f-%.0f", centBins[ic], centBins[ic+1]));
      	    hRhoMRatioDataoverMCNoC->SetTitle("Ratio #rho_{m} Data/MC");
      	    hRhoMRatioDataoverMCNoC->Divide(hRhoMIncl[id]);
      	    crho2->cd(4);
      	    hRhoMRatioDataoverMCNoC->GetYaxis()->SetRangeUser(-0.1, 5);
      	    hRhoMRatioDataoverMCNoC->Draw("sames");

     	 }
      	 
      	 crho->cd(3);
      	 gPad->SetLogy();
      	 hRhoMvsC[id]->Draw("sames");

      	 crho2->cd(3);
      	 gPad->SetLogy();
      	 if(id == 0) hRhoMvsC[id]->Draw("sames");
      	 else if(ic == 0) hRhoMIncl[id]->Draw("sames");

      	 gMeanRhoMvsC[id]->SetPoint(ic, ic, hRhoMvsC[id]->GetMean());
      	 gMeanRhoMvsC[id]->SetPointError(ic, 0.5, hRhoMvsC[id]->GetMeanError());
      	 gSigmaRhoMvsC[id]->SetPoint(ic, ic, hRhoMvsC[id]->GetRMS());
      	 gSigmaRhoMvsC[id]->SetPointError(ic, 0.5, hRhoMvsC[id]->GetRMSError());
      	 
      } //loop on data and MC
   } //centrality
   
   crho->cd(2);
   legC->Draw();
   crho2->cd(2);
   legC->Draw();
   
   SaveCv(crho);
   SaveCv(crho2);
   
   for(Int_t id = 0; id<ninput; id++){ //loop on data and MC
      crhoparam->cd(1);
      if(id == 0) gMeanRhovsC[id]->Draw("AP");
      else{
      	 gMeanRhovsC[id]->Draw("P");
      }
      gMeanRhoIncl[id]->Draw("P");
      
      crhoparam->cd(3);
      if(id == 0) gMeanRhoMvsC[id]->Draw("AP");
      else{
      	 gMeanRhoMvsC[id]->Draw("P");
      }
      gMeanRhoMIncl[id]->Draw("P");
      
      crhoparam->cd(2);
      if(id == 0) gSigmaRhovsC[id]->Draw("AP");
      else{
      	 gSigmaRhovsC[id]->Draw("P");
      }
      gSigmaRhoIncl[id]->Draw("P");
      
      crhoparam->cd(4);
      gSigmaRhoMvsC[id]->GetYaxis()->SetRangeUser(0.018, 0.034);
      if(id == 0) gSigmaRhoMvsC[id]->Draw("AP");
      else{
      	 gSigmaRhoMvsC[id]->Draw("P");
      }
      gSigmaRhoMIncl[id]->Draw("P");
   }
   
   SaveCv(crhoparam);
}
//------------------------------------------------------------------------------------------------
/*
void MatchDataModel(TString pathData, TString pathModel){
   
   TString filename = "RhoDistributions.root";
   
   // toy model, read file
   TFile *fModel = new TFile(Form("%s/%s",pathModel.Data(), filename.Data()));
   Int_t nBins = (fModel->GetEntries() - 2 ) / 2;
   TH1F* hRhoModel[nBins];
   TH1F* hRhomModel[nBins];
   
   //data read file
   TFile *fData = new TFile(Form("%s/%s",pathData.Data(), filename.Data()));
   
   TH1F* hRhoData = (TH1F*)fData->Get("hRhoIncldata");
   TH1F* hRhomData = (TH1F*)fData->Get("hRhoMIncldata");

   
   
}
*/
//------------------------------------------------------------------------------------------------
void RhopPbvsPYTHIApTLeadingFixedInput(){
   TString data = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train852/AnalysisResults.root";
   TString mc = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train851/outputpTHardBins/mergeRuns/AnalysisResults.root";

   RhopPbvsPYTHIApTLeading(data, mc, "Deriv");
}
//-------------------------------------------------------------------------------------
void RhopPbvsPYTHIApTLeading(TString dataPath, TString MCpath, TString tag){
   Double_t R = 0.4;
   TString type = "Charged";
   Int_t centBin = 0; 
   Bool_t bEmb = kFALSE;
   
   
   TString inputF[2]={dataPath, MCpath};
   const Int_t nCentBins = 5;
   Double_t centBins[nCentBins+1] = {0,10,20,40,60,100};
   TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC";
   TString hnamerho = "fh2RhoVsLeadJetPtVsCent", hnamerhom = "fh2RhoMVsLeadJetPtVsCent";
   Int_t nptl = 5;
   Double_t ptlead[nptl] = {0., 2., 5., 10., 15.};
   for(Int_t f = 0; f<1; f++){
      
      TList* list = ReadFile(inputF[f], Form("%s%s", listname.Data(), tag.Data()));
      if(list){
      	 Printf("Data found! all good");
      } else {
      	 Printf("List for data not found");
      	 return;
      }
      TH3F* h3rho = (TH3F*)list->FindObject(hnamerho);
      TH3F* h3rhom = (TH3F*)list->FindObject(hnamerhom);
      
      if(!h3rho || !h3rhom) {
      	 Printf("h3rho %p h3rhom %p", h3rho, h3rhom );
      	 return;
      }
      TCanvas *crhoCent[nCentBins];
      TCanvas *crhomCent[nCentBins];
      for(Int_t ic = 0; ic<nCentBins; ic++){
      	 crhoCent[ic] = new TCanvas(Form("crhoCent%d", ic), Form("Rho centrality %.0f-%.0f", centBins[ic], centBins[ic+1]), 800, 800);
      	 crhomCent[ic] = new TCanvas(Form("crhomCent%d", ic), Form("Rho %.0f-%.0f", centBins[ic], centBins[ic+1]), 800, 800);
      }
      
      TCanvas *crhovspTlead = new TCanvas("crhovspTlead", "Rho vs pT leading jet", 800);
      TH2F *h2rhopTlead = (TH2F*)h3rho->Project3D("xy");
      crhovspTlead->cd();
      h2rhopTlead->Draw("colz");
      TH1D *hrhoall = h3rho->ProjectionX(Form("hrho"));
      //hrhoall->Draw();
      
      TLegend *legptl = new TLegend(0.5, 0.5, 0.8, 0.9);
      legptl->SetBorderSize(0);
      
      for(Int_t iptl = 0; iptl < nptl; iptl++){
      	 
      	 Int_t binptlm  = h3rho->GetYaxis()->FindBin(ptlead[iptl]);
      	 Int_t binptlx  = h3rho->GetYaxis()->GetNbins() -1;
      	 //Int_t binptlmrhom = h3rhom->GetYaxis()->FindBin(ptlead[iptl]);
      	 Printf("Bin pt leading %d , %f (bin %d to %d)",iptl, ptlead[iptl], binptlm, binptlx);
      	       	 
      	 for(Int_t ic = 0; ic<nCentBins; ic++){
      	    Int_t bincent [2] = {h3rho->GetZaxis()->FindBin(centBins[ic]), h3rho->GetZaxis()->FindBin(centBins[ic+1])-1};
      	    Printf("Centrality bin %.0f - %.0f (%d to %d)", centBins[ic], centBins[ic+1], bincent[0], bincent[1]);
      	    TH1D* hrho  = h3rho->ProjectionX(Form("hrho%dCent_pTL%d", ic, iptl), binptlm, binptlx, bincent[0], bincent[1]);
      	    Printf("Rho integral = %.0f", hrho->Integral());
      	    
      	    TH1D* hrhom = h3rhom->ProjectionX(Form("hrhom%dCent_pTL%d", ic, iptl), binptlm, binptlx, bincent[0], bincent[1]);
      	    //Printf("Rhom integral = %.0f", hrhom->Integral());
      	    
      	    hrho->Scale(1./hrho->Integral());
      	    hrhom->Scale(1./hrhom->Integral());
      	    
      	    hrho->GetXaxis()->SetRangeUser(0, 8);
      	    hrhom->GetXaxis()->SetRangeUser(0, 0.4);
      	    
      	    hrho->SetLineColor(colors[iptl]);
      	    hrho->SetMarkerColor(colors[iptl]);
      	    hrho->SetLineWidth(2);
      	    hrhom->SetLineColor(colors[iptl]);
      	    hrhom->SetMarkerColor(colors[iptl]);
      	    hrhom->SetLineWidth(2);
      	    
      	    crhoCent[ic] ->cd();
      	    hrho->Draw("sames");
      	    crhomCent[ic] ->cd();
      	    hrhom->Draw("sames");
      	    
      	    if(ic==0) legptl->AddEntry(hrho, Form("#it{p}_{T, lead} > %.0f", ptlead[iptl]), "l");
      	 }
      	 
      }
      
      h3rho->GetZaxis()->SetRange(0, -1);
      h3rhom->GetYaxis()->SetRange(0, -1);
      
      for(Int_t ic = 0; ic<nCentBins; ic++){
      	 crhoCent[ic] ->cd();
      	 legptl->Draw();
      	 crhomCent[ic] ->cd();
      	 legptl->Draw();
      	 SaveCv(crhoCent[ic]);
      	 SaveCv(crhomCent[ic]);
      }
   }
}

//-------------------------------------------------------------------------------------------

Double_t GetFactorMBJetTrigger(TString type, TString trig, Int_t methodnumb){

   TH1D*  hturn = GetTurnOnTrigger(type, trig, methodnumb);
   if(!hturn){
      return -1;
   }
   Double_t threshold=0;
   if(trig.Contains("J1") ) threshold = 60;
   if(trig.Contains("J2") ) threshold = 40;
   
   Int_t bin1 = hturn->FindBin(threshold), bin2 = hturn->GetNbinsX();
   Double_t factor = 0;
   for(Int_t i=bin1; i<bin2;i++){
      factor+=hturn->GetBinContent(i);
   }

   factor/=(Double_t)(bin2-bin1);
   return factor;
   
}

//-------------------------------------------------------------------------------------------

TH1D* GetTurnOnTrigger(TString type, TString trig, Int_t methodnumb){
   TString pathToTurnOn = //"/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output775-774/TurnOnCurves.root";
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output794-795/TurnOnCurves.root";
   TFile *fto = new TFile(pathToTurnOn);
   
   if(!fto->IsOpen()){
      Printf("%s not found, check if it exist or run TurnOnCurve.C", pathToTurnOn.Data());
      return 0x0;
   }
   TString hname = "hpTRJet";
   hname+=type;
   hname+=trig;
   hname+=methodnumb;
   
   TH1D* hturn = (TH1D*)fto->Get(hname);
   if(!hturn){
      Printf("%s not found, check", hname.Data());
   }
   return hturn;
}


//void SaveCv(TCanvas *c, TString suff){
//   
//   c->SaveAs(Form("%s%s.png", c->GetName(),suff.Data()));
//   c->SaveAs(Form("%s%s.pdf", c->GetName(),suff.Data()));
//   c->SaveAs(Form("%s%s.eps", c->GetName(),suff.Data()));
//}

Double_t GetQuantile(TH1 *h1, Double_t prob) {
  double q[1];
  double probs[1] = {prob};
  h1->GetQuantiles( 1, q, probs );
  Printf("probs: %f  q = %f",prob,q[0]);
  return q[0];
}



