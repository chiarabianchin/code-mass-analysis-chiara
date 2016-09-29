#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraphErrors.h>

#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/style.C>
#include </data/Work/MyCodeJetMass/utils/plotUtils.C>
#include </data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.h>

//
// root[] .L /data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.cxx+
// root[] .L plotJetMassLeadTrkDataUtil.C+
//


const Bool_t bNorm = kFALSE;//kTRUE;
TString strLeg[3] = {"ConstSub","AreaSub-Raw","AreaSub-Deriv"};
TString specifictags[3] = {"ConstSub_TC", "_TCRaw", "_TCDeriv"};
Int_t available[3] = {0, 0, 0};
//Int_t colors[7] = {1,kGreen+2,kBlue,kCyan,kRed,kMagenta,kYellow+1};
   
Bool_t InitUtil(PlotUtilEmcalJetMass *&util, TString *tags, TString strf = "AnalysisResults.root", Double_t R = 0.4, Double_t ptmin = 60., Double_t ptmax = 80., TString type = "Charged", Int_t centBin = 0, Bool_t bEmb = kFALSE);

Bool_t InitUtil(PlotUtilEmcalJetMass *&util, TString *tags, TString strf, Double_t R, Double_t ptmin, Double_t ptmax, TString type, Int_t centBin, Bool_t bEmb){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
  //gROOT->LoadMacro("$gitJetMass/utils/style.C");
  SetStyle(0);
 
  Double_t centMin[4] = {0.,10.,30.,50.};
  Double_t centMax[4] = {10.,30.,50.,80.};

  TString outName = Form("JetMassCent%dR%03dPtMin%.0fPtMax%.0fType%s.root",centBin,(Int_t)(R*100.),ptmin,ptmax,type.Data());

  //gROOT->LoadMacro("/data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.cxx+");
  util = new PlotUtilEmcalJetMass();
  util->SetInputFileName(strf);
  util->SetJetRadius(R);
  util->SetJetType(type);
  //util->SetCentBin(centBin);
  util->SetJetPtRange(ptmin,ptmax);
  
  util->LoadFile();

  util->SetTag(tags[0]); //"ConstSub_TC"
  util->SetConstTag("");
  if(util->LoadList()) available[0] = 1;

  //util->SetJetName();
  util->SetTag(tags[1]); //"_TCRaw"
  util->SetConstTag("");
  if(util->LoadList()) available[1] = 1;
  
  util->SetTag(tags[2]); //"_TCDeriv"
  //util->SetTag("Raw");
  //util->SetConstTag("ConstSub");
  if(util->LoadList()) available[2] = 1;

  Printf("Loaded lists %d-%d-%d", available[0],available[1],available[2]);
  return kTRUE;

}
void plotJetMassLeadTrkDataUtil(TString strf = "AnalysisResults.root", Double_t R = 0.4, Double_t ptmin = 60., Double_t ptmax = 80., TString type = "Charged", Int_t centBin = 0, Bool_t bEmb = kFALSE) {

   PlotUtilEmcalJetMass *util=0x0;
   InitUtil(util, specifictags, strf, R, ptmin, ptmax, type, centBin, bEmb);
  //TString strLeg[nLists] = {"Raw","Raw","Constituent"};
  TString strCent[4] = {"0-10%","10-30%","30-50%","50-80%"};

  const Int_t nLeadTrkBins = 6;
  Double_t leadTrkPt[nLeadTrkBins] = {0.,2.,6.,10.,15.,19.};


  TH1D *hpTJetAll[3];
  TH1D *hpTJetTagged[3];
  TH1D *hpTaggedMatch[3];
  util->SetMinLeadTrackPt(leadTrkPt[2]);
  util->SetMaxLeadTrackPt(leadTrkPt[5]);
  for(Int_t i = 0; i<3; i++) {// lists
     if(available[i]==0) {
	Printf("Skipping list %d",i);
	continue;
     }
     hpTJetAll[i] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetAll,i);
     hpTJetAll[i]->SetLineWidth(2);     
     hpTJetAll[i]->SetLineColor(colors[i]);
     hpTJetTagged[i] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTagged,i);
     hpTJetTagged[i]->SetLineColor(colors[i]);
     if(!hpTJetTagged[i]) Printf("%d not found", i);
     //hpTJetTaggedMatch[i] = util->GetJetPtDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,i);
     //if(!hpTJetTaggedMatch[i]) Printf("ERROR!!!"); hpTJetTaggedMatch[i]->SetLineColor(colors[2]);
  
  }
  
  TH1D *hMAll[3][nLeadTrkBins];
  TH1D *hMTagged[3][nLeadTrkBins];
  TH1D *hMTaggedMatch[3][nLeadTrkBins];

  TGraphErrors *grMeanPtLeadTrAll[3];
  TGraphErrors *grMeanPtLeadTrTagged[3];
  TGraphErrors *grMeanPtLeadTrTaggedMatch[3];
  for(Int_t i = 0; i<3; i++) {
    grMeanPtLeadTrAll[i] = new TGraphErrors();
    grMeanPtLeadTrTagged[i] = new TGraphErrors();
    grMeanPtLeadTrTaggedMatch[i] = new TGraphErrors();
  }
  for(Int_t j = 0; j<nLeadTrkBins; j++) {
    util->SetMinLeadTrackPt(leadTrkPt[j]);
    for(Int_t i = 0; i<3; i++) {
      if(available[i]==0) {
	Printf("Skipping list %d",i);
	continue;
      }
      hMAll[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll,i);
      hMTagged[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged,i);
      hMTaggedMatch[i][j] = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTaggedMatch,i);
      grMeanPtLeadTrAll[i]->SetPoint(grMeanPtLeadTrAll[i]->GetN(),leadTrkPt[j],hMAll[i][j]->GetMean());
      grMeanPtLeadTrTagged[i]->SetPoint(grMeanPtLeadTrTagged[i]->GetN(),leadTrkPt[j],hMTagged[i][j]->GetMean());
      grMeanPtLeadTrTaggedMatch[i]->SetPoint(grMeanPtLeadTrTaggedMatch[i]->GetN(),leadTrkPt[j],hMTaggedMatch[i][j]->GetMean());

      grMeanPtLeadTrAll[i]->SetPointError(grMeanPtLeadTrAll[i]->GetN()-1,0.,hMAll[i][j]->GetMeanError());
      grMeanPtLeadTrTagged[i]->SetPointError(grMeanPtLeadTrTagged[i]->GetN()-1,0.,hMTagged[i][j]->GetMeanError());
      grMeanPtLeadTrTaggedMatch[i]->SetPointError(grMeanPtLeadTrTaggedMatch[i]->GetN()-1,0.,hMTaggedMatch[i][j]->GetMeanError());
      
    }
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
    for(Int_t nt = 0; nt<nLeadTrkBins; nt++) {
      if(nt==0) norm = hMAll[i][nt]->Integral();
      if(norm>0.) {
	//Printf("norm: %f",norm);
	hMAll[i][nt]->Scale(1./norm);

	hMAll[i][nt]->SetLineColor(colors[nt]);
	hMAll[i][nt]->SetMarkerColor(hMAll[i][nt]->GetLineColor());
	hMAll[i][nt]->SetMarkerStyle(20);
	
	hMAll[i][nt]->DrawCopy("same");

	hMTagged[i][nt]->Scale(1./norm);
	hMTagged[i][nt]->SetLineColor(colors[nt]);
	hMTagged[i][nt]->SetMarkerColor(hMTagged[i][nt]->GetLineColor());
	hMTagged[i][nt]->SetMarkerStyle(24);
	
	hMTagged[i][nt]->DrawCopy("same");

	hMTaggedMatch[i][nt]->Scale(1./norm);
	hMTaggedMatch[i][nt]->SetLineColor(colors[nt]);
	hMTaggedMatch[i][nt]->SetMarkerColor(hMTaggedMatch[i][nt]->GetLineColor());
	hMTaggedMatch[i][nt]->SetMarkerStyle(27);
	
	hMTaggedMatch[i][nt]->DrawCopy("same");
      
	leg1[i]->AddEntry(hMTagged[i][nt],Form("#it{p}_{T,lead trk}=%.0f #LT#it{M}_{jet}#GT=%.1f RMS=%.1f",leadTrkPt[nt],hMTagged[i][nt]->GetMean(),hMTagged[i][nt]->GetRMS()),"p");
      }
    }
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
  TH1F *frame4 = DrawFrame(-0.5,20.5,4.,12.,"#it{p}_{T,lead track}","#LT#it{M}_{jet}#GT");
  for(Int_t i = 0; i<3; i++) {
    if(available[i]==0) continue;
    Printf("Drawing list %d",i);
    grMeanPtLeadTrAll[i]->SetMarkerColor(colors[i]);
    grMeanPtLeadTrAll[i]->SetMarkerStyle(20);
    grMeanPtLeadTrAll[i]->Draw("PL");

    grMeanPtLeadTrTagged[i]->SetMarkerColor(colors[i]);
    grMeanPtLeadTrTagged[i]->SetMarkerStyle(24);
    grMeanPtLeadTrTagged[i]->Draw("PL");

    grMeanPtLeadTrTaggedMatch[i]->SetMarkerColor(colors[i]);
    grMeanPtLeadTrTaggedMatch[i]->SetMarkerStyle(27);
    grMeanPtLeadTrTaggedMatch[i]->Draw("PL");
  }

  c1->SaveAs(Form("MassLeadTrkPtBinVsCent%d_%sPtMin%.0fPtMax%.0f.png",centBin,type.Data(), ptmin, ptmax));

  TString outName = Form("JetMassCent%dR%03dPtMin%.0fPtMax%.0fType%s.root",centBin,(Int_t)(R*100.),ptmin,ptmax,type.Data());

  TFile *fout = new TFile(outName.Data(),"RECREATE");
  for(Int_t i = 0; i<3; i++) {
    if(available[i]==0) continue;
    grMeanPtLeadTrAll[i]->Write(Form("grMeanPtLeadTrAll_List%d",i));
    grMeanPtLeadTrTagged[i]->Write(Form("grMeanPtLeadTrTagged_List%d",i));
    grMeanPtLeadTrTaggedMatch[i]->Write(Form("grMeanPtLeadTrTaggedMatch_List%d",i));

    for(Int_t nt = 0; nt<nLeadTrkBins; nt++) {
      hMAll[i][nt]->Write();
      hMTagged[i][nt]->Write();
      hMTaggedMatch[i][nt]->Write();
    }
  }
}

void plotJetMasspTjetbins(TString strf = "AnalysisResults.root", Double_t R = 0.4, TString type = "Charged", Int_t centBin = 0, Bool_t bEmb = kFALSE){
   Int_t colors[7] = {1,kGreen+2,kBlue,kCyan,kRed,kMagenta,kYellow+1};
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   //gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
   //gROOT->LoadMacro("$gitJetMass/utils/style.C");
   SetStyle(0);
   
   PlotUtilEmcalJetMass *util=0x0;
   const Int_t npTj = 4;
   Double_t pTjlims[npTj+1]={20.,40.,60.,80.,100.};
   InitUtil(util, specifictags, strf, R, pTjlims[0], pTjlims[npTj-1], type, centBin, bEmb);
   util->SetMinLeadTrackPt(2.); //leadTrkPt[2]
   util->SetMaxLeadTrackPt(50.); //leadTrkPt[5]
  
   //TString outName = Form("JetMassCent%dR%03dPtMin%.0fType%s.root",centBin,(Int_t)(R*100.),type.Data());
   TString outName = Form("JetMassCent%dR%03dType%s.root",centBin,(Int_t)(R*100.),type.Data());
      
   Printf("Loaded lists %d-%d-%d", available[0],available[1],available[2]);
   
   TCanvas *cjmptj = new TCanvas("cjmptj", "Jet mass in bins of jet pT", 1000,1000);
   cjmptj->Divide(2,2);
   TCanvas *cjmoneptbin = new TCanvas("cjmoneptbin", "Jet mass in one bin of jet pT", 800,800);
  
   TLegend *legpTj = new TLegend(0.2, 0.25, 0.7, 0.75);
   legpTj->SetBorderSize(0);
   legpTj->SetFillStyle(0);
   Int_t norm=1;
   TString tag="Tagged Jets";
   TCanvas *cMeanMass = new TCanvas("cMeanMass", "Mean mass", 800,800);
   for(Int_t i = 0; i<3; i++) {
      if(!available[i]) continue;
      TPaveText *pvmean = new TPaveText(0.2, 0.4, 0.5, 0.6,"NDC");
      pvmean->SetFillStyle(0);
      TGraphErrors *gMeanMass = new TGraphErrors(npTj);
      gMeanMass->SetName(Form("gMeanMass%s",strLeg[i].Data()));
      gMeanMass->SetMarkerStyle(20+i);
      for(Int_t j = 0; j < npTj; j++){
      	 
      	 Printf("J ====== %d ,pTjet = %.0f - %.0f", j, pTjlims[j],pTjlims[j+1]);
      	 TPaveText *pvpT = new TPaveText(0.6, 0.7, 0.8, 0.9, "NDC");
      	 pvpT->SetFillStyle(0);
      	 pvpT->SetBorderSize(0);
      	 pvpT->AddText(Form("p_{T,jet} %.0f - %.0f", pTjlims[j],pTjlims[j+1]));
      	 util->SetJetPtRange(pTjlims[j],pTjlims[j+1]);
      	 TH1D* h = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetTagged, i);
      	 if(h->GetEntries()==0)  {
      	    h = util->GetJetMassDistribution(PlotUtilEmcalJetMass::kJetAll, i);
      	    Printf("WARNING! using kJetAll, not Tagged for %s", strLeg[i].Data());
      	    tag = "All Jets";
      	 } else tag="Tagged Jets";
      	 h->SetName(Form("hJetMassfile%dpTj%d", i, j));
      	 
      	 h->SetMarkerColor(colors[i]);
      	 h->SetMarkerStyle(20+i);
      	 h->SetLineColor(colors[i]);
      	 Double_t maximum= (Double_t)h->GetMaximum();
      	 Printf("Maximum is %f", maximum);
      	 if(i==0){
      	    //legpTj->AddEntry(h, Form("#it{p}_{T,jet} = %.0f-%.0f", pTjlims[j],pTjlims[j+1]), "PL");
      	    
      	 }
      	 if(j==0){
      	    legpTj->AddEntry(h, Form("%s (%s)", strLeg[i].Data(), tag.Data()), "P");
      	 }
      	 Printf("------> %f", pTjlims[j]+(pTjlims[j+1]-pTjlims[j])*0.5);
      	 gMeanMass->SetPoint(j, pTjlims[j]+(pTjlims[j+1]-pTjlims[j])*0.5, h->GetMean());
      	 gMeanMass->SetPointError(j, (pTjlims[j+1]-pTjlims[j])*0.5, h->GetRMS());
      	 
      	 norm=h->Integral();
      	 if(norm!=0)h->Scale(1./norm);
      	 //pvmean->AddText(Form("%d) Mean %.0f, sigma %.0f", j, h->GetMean(), 
      	 //h->GetRMS()));
      	 h->GetYaxis()->SetRangeUser(0.,0.4);
      	 h->GetXaxis()->SetRangeUser(-5,20.);
      	 
      	 cjmptj->cd(j+1);
      	 //if(j=0) h->Draw("P");
      	 //else 
      	 h->DrawClone("Psames");
      	 pvpT->Draw();
      	 if(j==1) {
      	    cjmoneptbin->cd();
      	    h->DrawClone("Psames");
      	    pvpT->Draw();

      	 }
      }
      cjmptj->cd(4);
      legpTj->Draw();
      cjmoneptbin->cd();
      legpTj->Draw();
      //cjmptj->cd(i+1);
      //pvmean->Draw();
      cMeanMass->cd();
      if(i==0) gMeanMass->Draw("AP");
      else gMeanMass->Draw("P");
  }
  
  cjmptj->SaveAs(Form("%s.pdf",cjmptj->GetName()));
}

