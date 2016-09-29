Int_t colors[7] = {1,kGreen+2,kBlue,kCyan,kRed,kMagenta,kYellow+1};
const Bool_t bNorm = kFALSE;//kTRUE;

void plotJetMassStructureUtil(TString strf = "AnalysisResults.root", Double_t R = 0.4, Double_t ptmin = 60., Double_t ptmax = 80., TString type = "Charged", Int_t centBin = 0, Bool_t bEmb = kFALSE) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
  gROOT->LoadMacro("$gitJetMass/utils/style.C");
  SetStyle(0);
 
  Double_t centMin[4] = {0.,10.,30.,50.};
  Double_t centMax[4] = {10.,30.,50.,80.};

  TString strLeg[3] = {"ConstSub","AreaSub-Raw","AreaSub-Deriv"};
  //TString strLeg[nLists] = {"Raw","Raw","Constituent"};
  TString strCent[4] = {"0-10%","10-30%","30-50%","50-80%"};

  const Int_t nLeadTrkBins = 1;
  //Double_t leadTrkPt[nLeadTrkBins] = {0.,2.,6.,10.,15.,19.};

  TString outName = Form("JetMassStructureCent%dR%03dPtMin%.0fPtMax%.0fType%s.root",centBin,(Int_t)(R*100.),ptmin,ptmax,type.Data());

  Int_t available[3] = {0,0,0};
  gROOT->LoadMacro("/data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.cxx+");
  PlotUtilEmcalJetMass *util = new PlotUtilEmcalJetMass();
  util->SetInputFileName(strf);
  util->SetTaskName("Structure");
  util->SetJetRadius(R);
  util->SetJetType(type);
  //util->SetCentBin(centBin);
  util->SetJetPtRange(ptmin,ptmax);
  util->SetJetPtPRange(ptmin,ptmax);
  util->LoadFile();

  util->SetTag("EfCorpTPois");
  util->SetConstTag("");
  if(util->LoadList()) available[0] = 1;

  /*
  //util->SetJetName();
  util->SetTag("Raw");
  util->SetConstTag("");
  if(util->LoadList()) available[1] = 1;
  
  util->SetTag("Deriv");
  //util->SetTag("Raw");
  //util->SetConstTag("ConstSub");
  if(util->LoadList()) available[2] = 1;
*/
  Printf("Loaded lists %d-%d-%d", available[0],available[1],available[2]);

  TH1D *hpTJetAll[3];
  TH1D *hpTJetTagged[3];
  TH1D *hpTaggedMatch[3];
  TH1D *hMAll[3];
  TH1D *hMTagged[3];
  for(Int_t i = 0; i<3; i++) {// lists
     if(available[i]==0) {
	Printf("Skipping list %d",i);
	continue;
     }
     hpTJetAll[i] = util->GetJetpTDetFromTHnSparse(PlotUtilEmcalJetMass::kJetAll,i);
     hpTJetAll[i]->SetLineWidth(2);     
     hpTJetAll[i]->SetLineColor(colors[i]);
     
  }
  
  
  //  TGraphErrors *grMeanPtLeadTrAll[3];
  //for(Int_t i = 0; i<3; i++) {
  //    grMeanPtLeadTrAll[i] = new TGraphErrors();
  //}
  for(Int_t i = 0; i<3; i++) {
     if(available[i]==0) {
	Printf("Skipping list %d",i);
	continue;
     }
     hMAll[i] = util->GetJetMassPartFromTHnSparse(PlotUtilEmcalJetMass::kJetAll,i);
     Printf("%p (%s) has %d entries",  hMAll[i] , hMAll[i]->GetName() ,  hMAll[i]->GetEntries());
     hMTagged[i] = util->GetJetMassDetFromTHnSparse(PlotUtilEmcalJetMass::kJetTagged,i);
     /*
     hMTaggedMatch[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,i);
     grMeanPtLeadTrAll[i]->SetPoint(grMeanPtLeadTrAll[i]->GetN(),leadTrkPt[j],hMAll[i][j]->GetMean());
     grMeanPtLeadTrTagged[i]->SetPoint(grMeanPtLeadTrTagged[i]->GetN(),leadTrkPt[j],hMTagged[i][j]->GetMean());
     grMeanPtLeadTrTaggedMatch[i]->SetPoint(grMeanPtLeadTrTaggedMatch[i]->GetN(),leadTrkPt[j],hMTaggedMatch[i][j]->GetMean());
     
     grMeanPtLeadTrAll[i]->SetPointError(grMeanPtLeadTrAll[i]->GetN()-1,0.,hMAll[i][j]->GetMeanError());
     grMeanPtLeadTrTagged[i]->SetPointError(grMeanPtLeadTrTagged[i]->GetN()-1,0.,hMTagged[i][j]->GetMeanError());
     grMeanPtLeadTrTaggedMatch[i]->SetPointError(grMeanPtLeadTrTaggedMatch[i]->GetN()-1,0.,hMTaggedMatch[i][j]->GetMeanError());
     */
  }

  Printf("Got all histos");

  Double_t mmax = 55.;
  if(R==0.2) mmax = 20.;
  if(R==0.4) mmax = 30.;

  TCanvas *c1 =new TCanvas("c1","c1: jet mass",900,700);
  c1->Divide(2,2);
  TH1F *frame[3];
  TLegend *leg1[3];
  for(Int_t i = 0; i<3; i++) {
     if(available[i]==0) continue;
     c1->cd(i+1);
     frame[i] = DrawFrame(-5.,mmax,0.,0.3,"#it{M}_{jet}","#frac{#it{N}_{jet}(#it{p}_{T,lead trk}>X)}{#it{N}_{jet}(all)}");
     leg1[i] = CreateLegend(0.25,0.6,0.52,0.94,strLeg[i].Data());
     leg1[i]->SetTextSize(leg1[i]->GetTextSize()*0.8);
     Double_t norm = 0.;
     
     
     hMAll[i]->SetLineColor(colors[i]);
     hMAll[i]->SetMarkerColor(hMAll[i]->GetLineColor());
     hMAll[i]->SetMarkerStyle(20);
     
     hMAll[i]->DrawCopy("same");
     
     hMTagged[i]->Scale(1./norm);
     hMTagged[i]->SetLineColor(colors[i]);
     hMTagged[i]->SetMarkerColor(hMTagged[i]->GetLineColor());
     hMTagged[i]->SetMarkerStyle(24);
     
     hMTagged[i]->DrawCopy("same");
     
     
     
     leg1[i]->AddEntry(hMTagged[i],Form("#LT#it{M}_{jet}#GT=%.1f RMS=%.1f",hMTagged[i]->GetMean(),hMTagged[i]->GetRMS()),"p");
  }
  
  TCanvas *cpTDistr = new TCanvas("cpTDistr","pT distributions", 800,800);
  cpTDistr->Divide(2,2);
  for(Int_t i = 0; i<3; i++) {
    if(available[i]==0) continue;
    c1->cd(i+1);
    leg1[i]->Draw();
    cpTDistr->cd(1);
    gPad->SetLogy();
    if(hpTJetAll[i] )   hpTJetAll[i]->Draw("sames");
    cpTDistr->cd(2);
    gPad->SetLogy();
    if(hpTJetTagged[i]) hpTJetTagged[i]->Draw("sames");
    //if(hpTJetTaggedMatch[i]) hpTJetTaggedMatch[i]->Draw("sames");
  }

  c1->cd(4);
  TH1F *frame4 = DrawFrame(-0.5,20.5,8.,18.,"#it{p}_{T,lead track}","#LT#it{M}_{jet}#GT");
  for(Int_t i = 0; i<3; i++) {
    if(available[i]==0) continue;
    Printf("Drawing list %d",i);
   //grMeanPtLeadTrAll[i]->SetMarkerColor(colors[i]);
   //grMeanPtLeadTrAll[i]->SetMarkerStyle(20);
   //grMeanPtLeadTrAll[i]->Draw("PL");
   //
   //grMeanPtLeadTrTagged[i]->SetMarkerColor(colors[i]);
   //grMeanPtLeadTrTagged[i]->SetMarkerStyle(24);
   //grMeanPtLeadTrTagged[i]->Draw("PL");
   //
   //grMeanPtLeadTrTaggedMatch[i]->SetMarkerColor(colors[i]);
   //grMeanPtLeadTrTaggedMatch[i]->SetMarkerStyle(27);
   //grMeanPtLeadTrTaggedMatch[i]->Draw("PL");
  }

  c1->SaveAs(Form("MassLeadTrkPtBinVsCent%d_%s.png",centBin,type.Data()));


  TFile *fout = new TFile(outName.Data(),"RECREATE");
  for(Int_t i = 0; i<3; i++) {
    if(available[i]==0) continue;
    //grMeanPtLeadTrAll[i]->Write(Form("grMeanPtLeadTrAll_List%d",i));
    //grMeanPtLeadTrTagged[i]->Write(Form("grMeanPtLeadTrTagged_List%d",i));
    //grMeanPtLeadTrTaggedMatch[i]->Write(Form("grMeanPtLeadTrTaggedMatch_List%d",i));

   
      hMAll[i]->Write();
      hMTagged[i]->Write();
      //hMTaggedMatch[i]->Write();
  }
}
