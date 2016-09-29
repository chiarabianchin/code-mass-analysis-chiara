#include<THnSparse.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TH3D.h>
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
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
void CompareCorrected(TString strIn[], TString strLst[], TString strLeg[], Int_t n=2);
void DrawTemplate(Double_t pTJMin, Double_t pTJMax, TString strIn, TString strLstIn, TString strOut, TString strLstOut);

TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Int_t* pTPBin, Int_t axisP, char* hprjname);
TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Double_t* pTP, Int_t axisP, char* hprjname);

void DefNewBins(TH1F* hTrig, Double_t* rebinpt, Int_t newnptbins);

Int_t GetJetPtBinTemp(Double_t ptj, Double_t delta, Double_t min = 0);

void plotDetectorResponseJetByJetCorr(TString strIn="JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC.root", TString strLst = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", Bool_t setaxrange = kTRUE) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->LoadMacro("$gitJetMass/utils/utilsMV.C");
  //gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
  //gROOT->LoadMacro("$gitJetMass/utils/style.C");
  SetStyle(0);
  //gStyle->SetOptStat(1111);
  //SetPlotStyle() ;
  TFile *f = new TFile(strIn.Data());
  if(!f->IsOpen()){
     Printf("File %s not found", strIn.Data());
     return;
  }
  TList *lst = static_cast<TList*>(f->Get(strLst.Data()));

  THnSparse *fhnMassResponse = (THnSparse*)lst->FindObject("fhnMassResponse");
  THnSparse *fhnMassResponseCorr = (THnSparse*)lst->FindObject("fhnMassResponseCorr");
  fhnMassResponse->GetAxis(3)->SetTitle("#it{p}_{T, part}");
  fhnMassResponseCorr->GetAxis(3)->SetTitle("#it{p}_{T, part}");
  TH2D *h2 = (TH2D*)fhnMassResponse->Projection(3,2);
  TH2D *h2Corr = (TH2D*)fhnMassResponseCorr->Projection(3,2);
  TH2D *hm = (TH2D*)fhnMassResponse->Projection(1,0);
  TH2D *hmCorr = (TH2D*)fhnMassResponseCorr->Projection(1,0);
  

  h2->SetAxisRange(10.,150.,"X");
  h2->SetAxisRange(10.,150.,"Y");
  if(setaxrange) h2->SetAxisRange(0.,1.e-5,"Z");
  h2Corr->SetAxisRange(10.,150.,"X");
  h2Corr->SetAxisRange(10.,150.,"Y");
  if(setaxrange) h2Corr->SetAxisRange(0.,1.e-5,"Z");
  TCanvas* cresp = new TCanvas("cresp", "Detector Response", 800, 800);
  cresp->cd();
  
  h2->Draw("colz");
  
  hm->SetAxisRange(0.,30.,"X");
  hm->SetAxisRange(0.,30.,"Y");
  if(setaxrange) hm->SetAxisRange(0.,1.e-5,"Z");
  hmCorr->SetAxisRange(0.,30.,"X");
  hmCorr->SetAxisRange(0.,30.,"Y");
  if(setaxrange) hmCorr->SetAxisRange(0.,1.e-5,"Z");
  TF1 *f1 = new TF1("f1","x*1.",0.,140.);
  f1->SetLineColor(2);

  TCanvas *c1 = new TCanvas("c1","c1",630,600);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogz();
  h2->Draw("colz");
  f1->Draw("same");
  c1->cd(2);
  gPad->SetLogz();
  h2Corr->Draw("colz");
  f1->Draw("same");
  c1->cd(3);
  gPad->SetLogz();
  hm->Draw("colz");
  f1->Draw("same");
  c1->cd(4);
  gPad->SetLogz();
  hmCorr->Draw("colz");
  f1->Draw("same");
  
  TLegend *leg=new TLegend(0.45,0.5,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  SaveCv(c1);
  
  //pT analysis
  TCanvas *cpT=new TCanvas("cpT","p_{T} distributions", 1200,400);
  cpT->Divide(3,1);
  
  
  //TH1D *hpTPartCorr=(TH1D*)fhnMassResponseCorr->Projection(3);
  //hpTPartCorr->GetXaxis()->SetTitle("p_{T}");
  //hpTPartCorr->SetLineColor(kBlue+1);
  //hpTPartCorr->SetLineWidth(2);
  //hpTPartCorr->SetLineStyle(2);
  //leg->AddEntry(hpTPartCorr, "Particle l. corrected","L");
  
  TH1D *hpTDet1=(TH1D*)fhnMassResponse->Projection(2);
  Int_t newNpTBins = 50;
  Double_t newpTbinranges[newNpTBins];
  Double_t ptval = hpTDet1->GetBinLowEdge(1);
  Double_t rightEdge =  hpTDet1->GetBinLowEdge(hpTDet1->GetNbinsX()+1);
  Printf("Going from %.3f to %.3f", ptval, rightEdge);
  
  for(Int_t i=0; i< newNpTBins+1; i++){
     
     if(ptval<20)
     	newpTbinranges[i] =  ptval+i*2;
     if(ptval<50)
     	newpTbinranges[i] =  ptval+i*5;
     else if(ptval<80)
     	newpTbinranges[i] =  ptval+i*10;
     else if(ptval<120)
     	newpTbinranges[i] =  ptval+i*15;
     else newpTbinranges[i] =  ptval+i*20;
     Printf("%d pT %f", i, newpTbinranges[i]);
     if(newpTbinranges[i] > rightEdge){
     	newpTbinranges[i] = rightEdge;
     	newNpTBins = i-1;
     	break;
     }
     
  }
  //DefNewBins((TH1F*)hpTDet1, newpTbinranges, newNpTBins);
  TH1D *hpTDet = (TH1D*)hpTDet1->Rebin(newNpTBins, Form("%sreb", hpTDet1->GetName()), newpTbinranges);
  hpTDet->SetLineWidth(2);
  hpTDet->SetLineColor(kRed);
  hpTDet->GetXaxis()->SetTitle("p_{T}");
  leg->AddEntry(hpTDet, "Detector level","L");
  
  TH1D *hpTDetCorr=(TH1D*)(fhnMassResponseCorr->Projection(2))->Rebin(newNpTBins, "hreb", newpTbinranges);
  hpTDetCorr->SetLineWidth(2);
  hpTDetCorr->SetLineStyle(2);
  hpTDetCorr->SetLineColor(kRed);
  hpTDetCorr->GetXaxis()->SetTitle("p_{T}");
  leg->AddEntry(hpTDetCorr, "Detector l. corrected","L");
  
  TH1D *hpTPart=(TH1D*)(fhnMassResponse->Projection(3))->Rebin(newNpTBins, "hreb", newpTbinranges);
  hpTPart->GetXaxis()->SetTitle("p_{T}");
  hpTPart->SetLineWidth(2);
  hpTPart->SetLineColor(kBlue);
  leg->AddEntry(hpTPart, "Particle level","L");

  cpT->cd(1);
  gPad->SetLogy();
  hpTPart->Draw();
  hpTDet->Draw("sames");
  leg->Draw();
  cpT->cd(2);
  gPad->SetLogy();
  hpTPart->Draw();
  hpTDetCorr->Draw("sames");
  
  cpT->cd(3);
  gPad->SetLogy();
  hpTDet->Draw();
  hpTDetCorr->Draw("sames");
  SaveCv(cpT);
  
  TCanvas *cpTratios=new TCanvas("cpTratios","p_{T} distribution ratios", 1200,400);
  cpTratios->Divide(3,1);
  
  TH1D *hrPartDet=(TH1D*)hpTPart->Clone("hrPartDet");
  hrPartDet->Divide(hpTDet);
    
  TH1D *hrPartDetC=(TH1D*)hpTPart->Clone("hrPartDetC");
  hrPartDetC->Divide(hpTDetCorr);
  
  TH1D* hrDetDetC=(TH1D*)hpTDet->Clone("hrDetDetC");
  hrDetDetC->Divide(hpTDetCorr);
   
  cpTratios->cd(1);
  hrPartDet->Draw();
  cpTratios->cd(2);
  hrPartDetC->Draw();
  cpTratios->cd(3);
  hrDetDetC->Draw();
  SaveCv(cpTratios);
  

  //mass analysis
  Float_t pTinterval=20.,pTMin=40.; //GeV/c
  const Int_t npTbins = 2;
  Double_t massBinW = fhnMassResponse->GetAxis(1)->GetBinWidth(1);
  const Int_t nMassBins = fhnMassResponse->GetAxis(1)->GetNbins();
  Int_t dmassMin = -30, dmassMax = massBinW * nMassBins + dmassMin; 
  
  TH1D* hMeanDMvspT = new TH1D("hMeanDMvspT", "Mean Delta Mass ; #it{p}_{T,jet} (GeV/#it{c}) ; Mean Delta Mass (GeV/#it{c^{2}})", npTbins, pTMin, pTMin+npTbins*pTinterval);
  TH1D* hMeanDMCorrvspT = new TH1D("hMeanDMCorrvspT", "Mean Delta Mass Corrected; #it{p}_{T,jet} (GeV/#it{c}) ; Mean Delta Mass (GeV/#it{c^{2}})", npTbins, pTMin, pTMin+npTbins*pTinterval);
  TH1D* hWidthDMvspT = new TH1D("hSigmaDMvspT", "Sigma Delta Mass; #it{p}_{T,jet} (GeV/#it{c}) ; Width Delta Mass (GeV/#it{c^{2}})", npTbins, pTMin, pTMin+npTbins*pTinterval);
  TH1D* hWidthDMCorrvspT = new TH1D("hSigmaDMCorrvspT", "Sigma Delta Mass Corrected; #it{p}_{T,jet} (GeV/#it{c}) ; Width Delta Mass (GeV/#it{c^{2}})", npTbins, pTMin, pTMin+npTbins*pTinterval);
  hMeanDMvspT    ->SetMarkerStyle(20);
  hMeanDMCorrvspT->SetMarkerStyle(21);
  hWidthDMvspT    ->SetMarkerStyle(20);
  hWidthDMCorrvspT->SetMarkerStyle(21);
       
  for(Int_t ipT=0;ipT<npTbins;ipT++){ //loop on pT det bins
     TString ptstring=Form("%.0f<p_{T,jet}<%.0f GeV/c",pTMin, pTMin+pTinterval);
     
     //TPaveText *pvpT=new TPaveText(0.5,0.7.0.8.0.8."NDC");
     //pvpT->SetBorderSize(0);
     //pvpT->SetFillStyle(0);
     //pvpT->AddText(ptstring);
     
     Int_t SelPartAndDetPT=-1; //-1 = pt Det only, 0 = pt part only, 1 = pt part and det, 2 = pt part and min pt det, 3 = pt part and max pt det SET THIS!
     TString suff=Form("pTMin%.0f",pTMin);
     
     
     Int_t pTPBin[2] = {fhnMassResponse->GetAxis(3)->FindBin(pTMin), fhnMassResponse->GetAxis(3)->FindBin(pTMin+pTinterval)-1};
     Int_t pTDBin[2] = {fhnMassResponse->GetAxis(2)->FindBin(pTMin), fhnMassResponse->GetAxis(2)->FindBin(pTMin+pTinterval)-1};
     Printf("Projections in pT bin range %d-%d",pTPBin[0],pTPBin[1]);
     if(SelPartAndDetPT > -1){
     	fhnMassResponse->    GetAxis(3)->SetRange(pTPBin[0],pTPBin[1]); //pT part
     	fhnMassResponseCorr->GetAxis(3)->SetRange(pTPBin[0],pTPBin[1]); //pT part
     	
     	if(SelPartAndDetPT) {
     	   suff+="PartAndDet";
     	   if(SelPartAndDetPT == 2) {
     	      pTDBin[1]=-1; //overflow bin
     	      suff+="Min";     	
     	   }
     	   if(SelPartAndDetPT == 3) {
     	      pTDBin[0]=0;  //underflow bin
     	      suff+="Max";
     	   }
     	   fhnMassResponse-> GetAxis(2)->SetRange(pTDBin[0],pTDBin[1]); //pT det; remove feed-in and feed-out
     	   fhnMassResponseCorr->GetAxis(2)->SetRange(pTDBin[0],pTDBin[1]); //pT det; remove feed-in and feed-out
     	   
     	}
     } else {
     	suff+="Det";
     	
     	fhnMassResponse-> GetAxis(2)->SetRange(pTDBin[0],pTDBin[1]);
     	fhnMassResponseCorr->GetAxis(2)->SetRange(pTDBin[0],pTDBin[1]);
     }
     
     
     TH1D* hMassPart =     (TH1D*)fhnMassResponse    ->Projection(1);  hMassPart     ->SetName("hMassPart"); hMassPart->SetTitle("hMassPart;M");
     TH1D* hMassDet  =     (TH1D*)fhnMassResponse    ->Projection(0); hMassDet      ->SetName("hMassDet" ); hMassDet->SetTitle("hMassDet;M");
     //TH1D* hMassPartCorr = (TH1D*)fhnMassResponseCorr->Projection(1); hMassPartCorr ->SetName("hMassPartCorr"); hMassPartCorr->SetTitle("hMassPartCorr;M");
     TH1D* hMassDetCorr  = (TH1D*)fhnMassResponseCorr->Projection(0); hMassDetCorr  ->SetName("hMassDetCorr" ); hMassDetCorr->SetTitle("hMassDetCorr;M");
     
     hMassPart->SetLineWidth(2);
     hMassPart->SetLineColor(kBlue);
     
     hMassDet ->SetLineWidth(2);
     hMassDet->SetLineColor(kRed);
     
     //hMassPartCorr->SetLineColor(kBlue+1);
     //hMassPartCorr->SetLineWidth(2);
     //hMassPartCorr->SetLineStyle(2);
     
     hMassDetCorr->SetLineWidth(2);
     hMassDetCorr->SetLineStyle(2);
     hMassDetCorr->SetLineColor(kRed);
     
     TPaveText *pv=new TPaveText(0.2,0.4,0.45,0.55,"NDC");
     pv->SetFillStyle(0);
     pv->SetBorderSize(0);
     pv->AddText(ptstring);
     
     TCanvas *cM = new TCanvas(Form("cM%d",ipT), Form("Mass distributions %d",ipT),1200,400);
     cM->Divide(3,1);
     
     cM->cd(1);
     gPad->SetLogy();
     hMassPart->Draw();
     hMassDet->Draw("sames");
     leg->Draw();
     cM->cd(2);
     gPad->SetLogy();
     hMassPart->Draw();
     hMassDetCorr->Draw("sames");
     pv->Draw();
     
     cM->cd(3);
     gPad->SetLogy();
     hMassDet->Draw();
     hMassDetCorr->Draw("sames");
     SaveCv(cM,suff);
     
     TCanvas *cMassratios=new TCanvas(Form("cMassratios%d",ipT),Form("Mass distribution ratios %d",ipT), 1200,400);
     cMassratios->Divide(3,1);
     
     TH1D *hrMassPartDet=(TH1D*)hMassPart->Clone("hrMassPartDet");
     hrMassPartDet->Divide(hMassDet);
     
     TH1D *hrMassPartDetC=(TH1D*)hMassPart->Clone("hrMassPartDetC");
     hrMassPartDetC->Divide(hMassDetCorr);
     
     TH1D* hrMassDetDetC=(TH1D*)hMassDet->Clone("hrMassDetDetC");
     hrMassDetDetC->Divide(hMassDetCorr);
     
     cMassratios->cd(1);
     hrMassPartDet->Draw();
     cMassratios->cd(2);
     hrMassPartDetC->Draw();
     cMassratios->cd(3);
     hrMassDetDetC->Draw();
     SaveCv(cMassratios,suff);
     
     
     //delta Mass
     
     TH2D* h2Mass = (TH2D*)fhnMassResponse->Projection(0,1); 
     h2Mass->SetName(Form("h2Mass%d",ipT));
     TH1D *hDeltaMass = new TH1D(Form("hDeltaMass%d", ipT), Form("Delta Mass %d ; #it{M}_{det}- #it{M}_{part} (GeV/#it{c^{2}})", ipT), nMassBins, dmassMin, dmassMax);
     hDeltaMass->SetLineWidth(2);
     hDeltaMass->SetLineColor(kRed);
     TH2D* h2MassCorr = (TH2D*)fhnMassResponseCorr->Projection(0,1); 
     h2MassCorr->SetName(Form("h2MassCorr%d",ipT));
     TH1D *hDeltaMassCorr = new TH1D(Form("hDeltaMassCorr%d", ipT), Form("Delta Mass Corr %d ; #it{M}_{det}- #it{M}_{part} (GeV/#it{c^{2}})", ipT), nMassBins, dmassMin, dmassMax);
     hDeltaMassCorr->SetLineWidth(2);
     hDeltaMassCorr->SetLineStyle(2);
     hDeltaMassCorr->SetLineColor(kRed);
     
     for(Int_t iX = 0; iX < nMassBins; iX++){
     	for(Int_t iY = 0; iY < nMassBins; iY++){
     	   
     	   Double_t deltaMass = h2Mass->GetXaxis()->GetBinCenter(iX+1) -  h2Mass->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContent = h2Mass->GetBinContent(iX+1,iY+1);
     	   hDeltaMass->Fill(deltaMass,binContent);
     	   
     	   Double_t deltaMassCorr = h2MassCorr->GetXaxis()->GetBinCenter(iX+1) -  h2MassCorr->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContentCorr = h2MassCorr->GetBinContent(iX+1,iY+1);
     	   hDeltaMassCorr->Fill(deltaMassCorr,binContentCorr);
     	   
     	}
     }
     
     TCanvas *cDeltaMass = new TCanvas(Form("cDeltaMass%d", ipT), Form("Delta Mass %d", ipT), 800,800);
     
     //cDeltaMass->Divide(1,2); cDeltaMass->cd(2); h2Mass->Draw("colz");
     cDeltaMass->cd(1);     
     hDeltaMass->Draw();
     hDeltaMassCorr->Draw("sames");
     
     TPaveText* pvinfo = new TPaveText(0.1,0.2,0.45,0.3,"NDC");
     pvinfo->SetFillStyle(0);
     pvinfo->SetBorderSize(0);
     Double_t meanDM = hDeltaMass->GetMean(), meanDMerr = hDeltaMass->GetMeanError(), widthDM = hDeltaMass->GetRMS();
     Double_t meanDMCorr = hDeltaMassCorr->GetMean(), meanDMCorrerr = hDeltaMassCorr->GetMeanError(), widthDMCorr = hDeltaMassCorr->GetRMS();
     //fill summay histograms
     hMeanDMvspT->SetBinContent(ipT+1, meanDM);
     hMeanDMvspT->SetBinError(ipT+1, meanDMerr);
     hMeanDMCorrvspT->SetBinContent(ipT+1, meanDMCorr);
     hMeanDMCorrvspT->SetBinError(ipT+1, meanDMCorrerr);
     hWidthDMvspT->SetBinContent(ipT+1, widthDM);
     hWidthDMCorrvspT->SetBinContent(ipT+1, widthDMCorr);
     
     //print info on canvas
     pvinfo->AddText(Form("(det) Mean = %.2f +- %.2f, Width %.2f ", meanDM, meanDMerr, widthDM));
     pvinfo->AddText(Form("(det corr) Mean = %.2f +- %.2f, Width %.2f ", meanDMCorr, meanDMCorrerr, widthDMCorr));
     
     cDeltaMass->cd(1);
     pvinfo->Draw();
     pv->Draw();
     SaveCv(cDeltaMass,suff);
     
     pTMin+=pTinterval;
  }
  
  TCanvas *cSummary = new TCanvas("cSummary", "Summary", 500,900);
  cSummary->Divide(1,2);
  //set axis range
  hMeanDMvspT->SetMinimum(0);
  hWidthDMvspT->SetMinimum(1.5);
  hWidthDMvspT->SetMaximum(2.5);  
  cSummary->cd(1);
  hMeanDMvspT->Draw("P");
  hMeanDMCorrvspT->Draw("Psames");
  
  cSummary->cd(2);
  hWidthDMvspT->Draw("P");
  hWidthDMCorrvspT->Draw("Psames");
  
  SaveCv(cSummary);


  //delta pT
  //reset the pT ranges
  fhnMassResponse->GetAxis(2)->SetRange(0,-1);
  fhnMassResponse->GetAxis(3)->SetRange(0,-1);
  fhnMassResponseCorr->GetAxis(2)->SetRange(0,-1);
  fhnMassResponseCorr->GetAxis(3)->SetRange(0,-1);
  
  //bins of particle level pT for the loop
  const Int_t npTPartBins=8;
  Double_t pTPartBinW = 20.; //GeV/c
  Double_t pTPartMin = 0, pTPartMax = pTPartBinW * npTPartBins + pTPartMin;
  //bins of detector level pT for the delta pT spectra
  Double_t pTBinW = fhnMassResponse->GetAxis(2)->GetBinWidth(1);
  const Int_t npTthnsBins = fhnMassResponse->GetAxis(2)->GetNbins();
  //Int_t dpTMin = -40, dpTMax = pTBinW * npTthnsBins + dpTMin; 
  Double_t dpTMin = -80, dpTMax = pTBinW * npTthnsBins + dpTMin;
  Double_t dpTRelMin = -5, dpTRelMax = 5;
  // delta pT vs pT Part
  TH1D* hMeanDpTvspTPart = new TH1D("hMeanDpTvspTP", "Mean Delta pT", npTPartBins, pTPartMin, pTPartMax);
  TH1D* hMeanDpTCorrvspTPart = new TH1D("hMeanDpTCorrvspTP", "Mean Delta pT Corrected", npTPartBins, pTPartMin, pTPartMax);
  TH1D* hWidthDpTvspTPart = new TH1D("hSigmaDpTvspTP", "Sigma Delta pT", npTPartBins, pTPartMin, pTPartMax);
  TH1D* hWidthDpTCorrvspTPart = new TH1D("hSigmaDpTCorrvspTP", "Sigma Delta pT Corrected", npTPartBins, pTPartMin, pTPartMax);
  hMeanDpTvspTPart    ->SetMarkerStyle(20);
  hMeanDpTCorrvspTPart->SetMarkerStyle(21);
  hWidthDpTvspTPart    ->SetMarkerStyle(20);
  hWidthDpTCorrvspTPart->SetMarkerStyle(21);
 
  TH2D* h2pT = (TH2D*)fhnMassResponse->Projection(3,2); 
  h2pT->SetName("h2pT");
  TH2D* h2pTCorr = (TH2D*)fhnMassResponseCorr->Projection(3,2); 
  h2pTCorr->SetName("h2pTCorr");
  Printf("Axis name X = %s, Y = %s",fhnMassResponse->GetAxis(2)->GetTitle(),fhnMassResponse->GetAxis(3)->GetTitle());
  TCanvas *cDeltapT = new TCanvas("cDeltapT", "Delta pT ", 1000,1000);
  cDeltapT->Divide(3,3);
  TCanvas *cDeltapTRel = new TCanvas("cDeltapTRel", "Delta pT Rel ", 1000,1000);
  cDeltapTRel->Divide(3,3);
  //draw 2D projections
  //cDeltapT->cd(npTPartBins+1);
  //h2pT->DrawClone("colz");
  //cDeltapT->cd(npTPartBins+2);
  //h2pTCorr->DrawClone("colz");
  
  //Printf("Range X %f, %f (bin max = %d)",h2pT->GetXaxis()->GetBinLowEdge(1),h2pT->GetXaxis()->GetBinLowEdge(npTthnsBins+1) , npTthnsBins);
  //Printf("Range Y %f, %f (bin max = %d)",h2pT->GetYaxis()->GetBinLowEdge(1),h2pT->GetYaxis()->GetBinLowEdge(npTthnsBins+1) , npTthnsBins);
 
  //pT integrated
  TH1D *hAllDeltapT = new TH1D("hAllDeltapT", "All Delta pT ; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) ", npTthnsBins, dpTMin, dpTMax);
  hAllDeltapT->SetLineWidth(2);
  hAllDeltapT->SetLineColor(kRed);
  TH1D *hAllDeltapTCorr = new TH1D("hAllDeltapTCorr", "All Delta pT Corr; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) ", npTthnsBins, dpTMin, dpTMax);
  hAllDeltapTCorr->SetLineWidth(2);
  hAllDeltapTCorr->SetLineStyle(2);
  hAllDeltapTCorr->SetLineColor(kGreen+3);
  TH1D *hAllDeltapTRel = new TH1D("hAllDeltapTRel", "All Delta pT ; (#it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part})/ #it{p}_{T,jet}^{part} (GeV/#it{c}) ", npTthnsBins, dpTRelMin, dpTRelMax);
  hAllDeltapTRel->SetLineWidth(2);
  hAllDeltapTRel->SetLineColor(kRed);
  TH1D *hAllDeltapTRelCorr = new TH1D("hAllDeltapTRelCorr", "All Delta pT Corr; (#it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part})/#it{p}_{T,jet}^{part} (GeV/#it{c}) ", npTthnsBins, dpTRelMin, dpTRelMax);
  hAllDeltapTRelCorr->SetLineWidth(2);
  hAllDeltapTRelCorr->SetLineStyle(2);
  hAllDeltapTRelCorr->SetLineColor(kGreen+3);
  
  TCanvas *cDeltapTIntegrated=new TCanvas("cDeltapTIntegrated", "Delta pT integrated in pTpart", 800,800);
  TCanvas *cDeltapTRelIntegrated=new TCanvas("cDeltapTRelIntegrated", "Delta pT Rel integrated in pTpart", 800,800);
  for(Int_t iX = 0; iX < npTthnsBins; iX++){
     for(Int_t iY = 0; iY < npTthnsBins; iY++){
     	
     	Double_t deltapT = h2pT->GetXaxis()->GetBinCenter(iX+1) -  h2pT->GetYaxis()->GetBinCenter(iY+1);
     	Double_t deltapTRel = deltapT/h2pT->GetYaxis()->GetBinCenter(iY+1);
     	
     	Double_t binContent = h2pT->GetBinContent(iX+1,iY+1);
     	hAllDeltapT->Fill(deltapT,binContent);
     	hAllDeltapTRel->Fill(deltapTRel,binContent);
     	//if((iX==0 || iX==10 || iX==50 || iX==100) && iY==10) {Printf("iX = %d , deltapT = %f, bincontent = %f",iX,deltapT,binContent);}
     	Double_t deltapTCorr = h2pTCorr->GetXaxis()->GetBinCenter(iX+1) -  h2pTCorr->GetYaxis()->GetBinCenter(iY+1);
     	Double_t deltapTRelCorr = deltapTCorr/h2pTCorr->GetYaxis()->GetBinCenter(iY+1);
     	Double_t binContentCorr = h2pTCorr->GetBinContent(iX+1,iY+1);
     	hAllDeltapTCorr->Fill(deltapTCorr,binContentCorr);
     	hAllDeltapTRelCorr->Fill(deltapTRelCorr,binContentCorr);
     }
  }
  cDeltapTIntegrated->cd();
  gPad->SetLogy();
  hAllDeltapT->Draw();
  hAllDeltapTCorr->Draw("sames");
  cDeltapTRelIntegrated->cd();
  gPad->SetLogy();
  hAllDeltapTRel->Draw();
  hAllDeltapTRelCorr->Draw("sames");

  TH1D* hSum=0x0;
  for(Int_t ipTPart=0;ipTPart<npTPartBins; ipTPart++){ //loop on pT particle level bins
     Double_t pTP[2]={pTPartMin + pTPartBinW*ipTPart, pTPartBinW * (ipTPart+1) + pTPartMin};
     Int_t bin0 = h2pTCorr->FindBin(0, pTP[0]);
     Int_t bin1 = h2pTCorr->FindBin(0, pTP[1]);
     Printf("Set Range user %f, %f", pTP[0], pTP[1]);
     //Printf("Set Range %d, %d", bin0, bin1);
     Int_t binx0,biny0, binx1,biny1, zdummy;
     h2pTCorr->GetBinXYZ(bin0, binx0, biny0, zdummy);
     h2pTCorr->GetBinXYZ(bin1, binx1, biny1, zdummy);
     //TH1D* hpjpTDet=h2pT->ProjectionX(Form("hprojpTDet%d",ipTPart),biny0, biny1);
     //h2pT->GetYaxis()->SetRange(biny0, biny1);
     //h2pTCorr->GetYaxis()->SetRange(biny0, biny1);
     
     TH1D *hDeltapT = new TH1D(Form("hDeltapT%d",ipTPart), Form("Delta pT %d; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) ",ipTPart), npTthnsBins, dpTMin, dpTMax);
     hDeltapT->SetLineWidth(2);
     hDeltapT->SetLineColor(kRed);
     TH1D *hDeltapTRel = new TH1D(Form("hDeltapTRel%d",ipTPart), Form("Delta pT %d; (#it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part})/#it{p}_{T,jet}^{part} (GeV/#it{c}) ",ipTPart), npTthnsBins, dpTRelMin, dpTRelMax);
     hDeltapTRel->SetLineWidth(2);
     hDeltapTRel->SetLineColor(kRed);
     TH1D *hDeltapTCorr = new TH1D(Form("hDeltapTCorr%d",ipTPart), Form("Delta pT Corr; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) %d",ipTPart), npTthnsBins, dpTMin, dpTMax);
     hDeltapTCorr->SetLineWidth(2);
     hDeltapTCorr->SetLineStyle(2);
     hDeltapTCorr->SetLineColor(kGreen+2);
     TH1D *hDeltapTCorrRel = new TH1D(Form("hDeltapTCorrRel%d",ipTPart), Form("Delta pT Corr; (#it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part})/ #it{p}_{T,jet}^{part}(GeV/#it{c}) %d",ipTPart), npTthnsBins, dpTRelMin, dpTRelMax);
     hDeltapTCorrRel->SetLineWidth(2);
     hDeltapTCorrRel->SetLineStyle(2);
     hDeltapTCorrRel->SetLineColor(kGreen+2);
    
     Printf("pTPart between bin %d and %d",biny0, biny1-1);
     
     for(Int_t iX = 0; iX < npTthnsBins; iX++){
     	for(Int_t iY = biny0; iY < biny1; iY++){
     	   
     	   Double_t deltapT =  h2pT->GetXaxis()->GetBinCenter(iX+1) -  h2pT->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t deltapTRel = deltapT/ h2pT->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContent = h2pT->GetBinContent(iX+1,iY+1);
     	   hDeltapT->Fill(deltapT,binContent);
     	   hDeltapTRel->Fill(deltapTRel,binContent);
     	   //if((iX==0 || iX==10 || iX==50 || iX==100) && iY==10) {Printf("iX = %d , deltapT = %f, bincontent = %f",iX,deltapT,binContent);}
     	   Double_t deltapTCorr = h2pTCorr->GetXaxis()->GetBinCenter(iX+1) -  h2pTCorr->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t deltapTCorrRel = deltapTCorr/ h2pTCorr->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContentCorr = h2pTCorr->GetBinContent(iX+1,iY+1);
     	   hDeltapTCorr->Fill(deltapTCorr,binContentCorr);
     	   hDeltapTCorrRel->Fill(deltapTCorrRel,binContentCorr);
     	   
     	}
     }

     if(!hSum) {
     	hSum = (TH1D*)hDeltapTCorrRel->Clone("hSum");
     	hSum->SetLineColor(kBlack);
     	hSum->SetLineStyle(5);
     	hSum->SetLineWidth(3);
     	
     }
     else hSum->Add(hDeltapTCorrRel);     
     
     
     //cDeltapT->Divide(1,2); cDeltapT->cd(2); h2pT->Draw("colz"); 
     cDeltapT->cd(ipTPart+1);
     
     //hpjpTDet->Draw();
     hDeltapT->Draw();
     hDeltapTCorr->Draw("sames");
     TString pTPstring=Form("%.0f< #it{p}_{T,Part} < %.0f",pTP[0], pTP[1]);
     TPaveText *pv=new TPaveText(0.6,0.4,0.9,0.6,"NDC");
     pv->SetFillStyle(0);
     pv->SetBorderSize(0);
     pv->AddText(pTPstring);
     pv->Draw();

     cDeltapTRel->cd(ipTPart+1);
     gPad->SetLogy();
     hDeltapTRel->Draw();
     hDeltapTCorrRel->Draw("sames");
     pv->Draw();
     
     TPaveText* pvinfo1 = new TPaveText(0.1,0.2,0.45,0.3,"NDC");
     pvinfo1->SetFillStyle(0);
     pvinfo1->SetBorderSize(0);
     Double_t meanpT = hDeltapT->GetMean(), meanpTerr = hDeltapT->GetMeanError(), widthpT = hDeltapT->GetRMS();
     Double_t meanpTCorr = hDeltapTCorr->GetMean(), meanpTCorrerr = hDeltapTCorr->GetMeanError(), widthpTCorr = hDeltapTCorr->GetRMS();
     //fill summay histograms
     hMeanDpTvspTPart->SetBinContent(ipTPart+1, meanpT);
     hMeanDpTvspTPart->SetBinError(ipTPart+1, meanpTerr);
     hMeanDpTCorrvspTPart->SetBinContent(ipTPart+1, meanpTCorr);
     hMeanDpTCorrvspTPart->SetBinError(ipTPart+1, meanpTCorrerr);
     hWidthDpTvspTPart->SetBinContent(ipTPart+1, widthpT);
     hWidthDpTCorrvspTPart->SetBinContent(ipTPart+1, widthpTCorr);
     
     //print info on canvas
     pvinfo1->AddText(Form("(det) Mean = %.2f +- %.2f, Width %.2f ", meanpT, meanpTerr, widthpT));
     pvinfo1->AddText(Form("(det corr) Mean = %.2f +- %.2f, Width %.2f ", meanpTCorr, meanpTCorrerr, widthpTCorr));
     
     cDeltapT->cd(ipTPart+1);
     gPad->SetLogy();
     pvinfo1->Draw();
     
  }
  //compare integrated and sum pT bins
  cDeltapTRelIntegrated->cd();
  if(hSum) {
     //hSum->Draw("same");
     Printf("hSum mean = %.2f, width = %.2f, integral = %.2f", hSum->GetMean(), hSum->GetRMS(), hSum->Integral());
  }
  
  SaveCv(cDeltapT);
  SaveCv(cDeltapTRel);
  
  TCanvas *cSummarydpT=new TCanvas("cSummarydpT", "Summary Delta pT", 500,900);
  cSummarydpT->Divide(1,2);
  //set axis range
  hMeanDpTvspTPart->SetMinimum(0);
  hWidthDpTvspTPart->SetMinimum(1.5);
  hWidthDpTvspTPart->SetMaximum(2.5);  
  cSummarydpT->cd(1);
  hMeanDpTvspTPart->Draw("P");
  hMeanDpTCorrvspTPart->Draw("Psames");
  
  cSummarydpT->cd(2);
  hWidthDpTvspTPart->Draw("P");
  hWidthDpTCorrvspTPart->Draw("Psames");
  
  SaveCv(cSummarydpT);
  
  //as a function of Mass

  Int_t nMPthnsBins=fhnMassResponse->GetAxis(1)->GetNbins();
  Printf("Bins of mass %d, Total Range %f-%f",nMPthnsBins, fhnMassResponse->GetAxis(1)->GetBinLowEdge(1), fhnMassResponse->GetAxis(1)->GetBinLowEdge(nMPthnsBins));
  
  //bins of particle level Mass for the loop
  const Int_t nMPartBins=6;
  Double_t MPartBinW = 10.; //GeV/c
  Double_t MPartMin = 0, MPartMax = MPartBinW * nMPartBins + MPartMin;
  
  // delta pT vs M Part

  TH1D* hMeanDpTvsMPart = new TH1D("hMeanDpTvsMP", "Mean Delta pT", nMPartBins, MPartMin, MPartMax);
  TH1D* hMeanDpTCorrvsMPart = new TH1D("hMeanDpTCorrvsMP", "Mean Delta pT Corrected", nMPartBins, MPartMin, MPartMax);
  TH1D* hWidthDpTvsMPart = new TH1D("hSigmaDpTvsMP", "Sigma Delta pT", nMPartBins, MPartMin, MPartMax);
  TH1D* hWidthDpTCorrvsMPart = new TH1D("hSigmaDpTCorrvsMP", "Sigma Delta pT Corrected", nMPartBins, MPartMin, MPartMax);
  hMeanDpTvsMPart    ->SetMarkerStyle(20);
  hMeanDpTCorrvsMPart->SetMarkerStyle(21);
  hWidthDpTvsMPart    ->SetMarkerStyle(20);
  hWidthDpTCorrvsMPart->SetMarkerStyle(21);
  
  
  TCanvas *cDeltapTMPBins = new TCanvas("cDeltapTMPBins", "Delta pT in MParticle Bins ", 1000,1000);
  cDeltapTMPBins->Divide(3,2);
  
  for(Int_t iMPart=0;iMPart<nMPartBins;iMPart++){//loop on Particle Level Mass bins
     Double_t MP[2] = {
     MPartMin+MPartBinW*iMPart, MPartMin+  MPartBinW*(iMPart+1)};
     
     Int_t binMP[2] = {
     fhnMassResponse->GetAxis(1)->FindBin(MP[0]), 
     fhnMassResponse->GetAxis(1)->FindBin(MP[0])};
     
     fhnMassResponse->    GetAxis(1)->SetRange(binMP[0], binMP[1]);
     fhnMassResponseCorr->GetAxis(1)->SetRange(binMP[0], binMP[1]);
     
     Printf("M Part range = %f-%f", MP[0], MP[1]);
     Printf("M Part range bin = %d-%d", binMP[0], binMP[1]);
     
     TH2D* h2pTMPSel = fhnMassResponse->Projection(3,2);
     h2pTMPSel->SetName("h2pTMPSel");
     TH2D* h2pTCorrMPSel = fhnMassResponseCorr->Projection(3,2);
     h2pTCorrMPSel->SetName("h2pTCorrMPSel");
     
     TH1D *hDeltapTMBins = new TH1D(Form("hDeltapTMBins%d",iMPart), Form("Delta pT %d; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) ",iMPart), npTthnsBins, dpTMin, dpTMax);
     hDeltapTMBins->SetLineWidth(2);
     hDeltapTMBins->SetLineColor(kRed);
     TH1D *hDeltapTCorrMBins = new TH1D(Form("hDeltapTCorrMBins%d",iMPart), Form("Delta pT Corr; #it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{part} (GeV/#it{c}) %d",iMPart), npTthnsBins, dpTMin, dpTMax);
     hDeltapTCorrMBins->SetLineWidth(2);
     hDeltapTCorrMBins->SetLineStyle(2);
     hDeltapTCorrMBins->SetLineColor(kGreen+2);

     for(Int_t iX = 0; iX < npTthnsBins; iX++){
     	for(Int_t iY = 0; iY < npTthnsBins; iY++){
     	   
     	   Double_t deltapT = h2pTMPSel->GetXaxis()->GetBinCenter(iX+1) -  h2pT->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContent = h2pTMPSel->GetBinContent(iX+1,iY+1);
     	   hDeltapTMBins->Fill(deltapT,binContent);
     	   //if((iX==0 || iX==10 || iX==50 || iX==100) && iY==10) {Printf("iX = %d , deltapT = %f, bincontent = %f",iX,deltapT,binContent);}
     	   Double_t deltapTCorr = h2pTCorrMPSel->GetXaxis()->GetBinCenter(iX+1) -  h2pTCorrMPSel->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContentCorr = h2pTCorrMPSel->GetBinContent(iX+1,iY+1);
     	   hDeltapTCorrMBins->Fill(deltapTCorr,binContentCorr);
     	   
     	}
     }
     cDeltapTMPBins->cd(iMPart+1);
     gPad->SetLogy();
     hDeltapTMBins->Draw();
     hDeltapTCorrMBins->Draw("sames");
     
     TString MPstring=Form("%.0f< #it{M}_{Part} < %.0f",MP[0], MP[1]);
     TPaveText *pv=new TPaveText(0.6,0.4,0.9,0.6,"NDC");
     pv->SetFillStyle(0);
     pv->SetBorderSize(0);
     pv->AddText(MPstring);
     pv->Draw();
  }
  
  SaveCv(cDeltapTMPBins);
}

void LoadLibs()
{
  // Load common libraries (better too many than too few)
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  //  gSystem->Load("libCDB");
  //  gSystem->Load("libRAWDatabase");
 

  //load CGAL, Fastjet and SISCone
  //  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");
  gSystem->Load("libfastjetcontribfragile");

  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libSTAT.so");
  // gSystem->Load("libPWGPP.so");
  // gSystem->Load("libPWGHFbase");
  // gSystem->Load("libPWGDQdielectron");
  // gSystem->Load("libPWGHFhfe");
  // gSystem->Load("libEMCALUtils");
  // gSystem->Load("libPHOSUtils");
  // gSystem->Load("libPWGCaloTrackCorrBase");
  // gSystem->Load("libEMCALraw");
  // gSystem->Load("libEMCALbase");
  // gSystem->Load("libEMCALrec");
  // gSystem->Load("libTRDbase");
  // gSystem->Load("libVZERObase");
  // gSystem->Load("libVZEROrec");
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");

  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGJEEMCALJetTasks");

  gSystem->Load("libPWGJE.so");

  // include paths
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
  gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");
  gSystem->AddIncludePath("-I$FASTJET_ROOT/include");
}

void CompareCorrected(TString strIn[], TString strLst[], TString strLeg[], Int_t n){
    
   SetStyle(0);
   TCanvas *cmassInt = new TCanvas("cmassInt", "Mass integrated over pT",1000,1000);
   cmassInt->Divide(2,2);
   
   const Int_t npTPartBins=6;
   Double_t pTPartBinW = 20.; //GeV/c
   Double_t pTPartMin = 20., pTPartMax = pTPartBinW * npTPartBins + pTPartMin;
   Int_t axispTP=3;
      
   Int_t axispTD=2;
   TCanvas *cpT = new TCanvas("cpT", "Corrected pT Det spectra in pTP bins",800,800);
   TCanvas *cpTr = new TCanvas("cpTr", "Corrected/Uncorr pT Det specra in pTP bins",800,800);
   Int_t axisMDet = 0;
   TCanvas *cM=new TCanvas("cM", "Corrected Mass Det specra in pTP bins",800,800);
   //TCanvas *cMr=new TCanvas("cMr", "Corrected/Uncorr Mass Det specra in pTP bins",800,800);
   TLegend *legpTPBins = new TLegend(0.6,0.65,0.9,0.85);
   legpTPBins->SetFillStyle(0);
   legpTPBins->SetBorderSize(0);
   TLegend *leg = new TLegend(0.6,0.45,0.85,0.55);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   Int_t axispTD2 = 3;
   Int_t axispTP2 = 4;
   Int_t axisDM   = 0;
   Int_t axisDMRel= 1;
   
   TCanvas* cDMrel=new TCanvas("cDMrel", "Corrected DeltaM/M in pTP bins",800,800);
   
   for (Int_t in=0;in<n;in++){
      TList *l=ReadFile(strIn[in], strLst[in]);
      if(!l) {
      	 Printf("List not found");
      	 continue;
      }
      THnSparse *fhnMassResponse = (THnSparse*)l->FindObject("fhnMassResponse");
      THnSparse *fhnMassResponseCorr = (THnSparse*)l->FindObject("fhnMassResponseCorr");
      Int_t naxes = fhnMassResponse->GetNdimensions();
      
      THnSparse *fhnDeltaMass = (THnSparse*)l->FindObject("fhnDeltaMass");
      THnSparse *fhnDeltaMassCorr = (THnSparse*)l->FindObject("fhnDeltaMassCorr");
      
      for (Int_t iaxis=0;iaxis<naxes;iaxis++){
      	 TH1D *hproj = (TH1D*)fhnMassResponse->Projection(iaxis);
      	 hproj->SetName(Form("hproj%d_%d",iaxis,in));
      	 hproj->SetLineColor(colors[in]);
      	 TH1D *hprojCorr = (TH1D*)fhnMassResponseCorr->Projection(iaxis);
      	 hprojCorr->SetLineColor(colors[in]);
      	 if(in>0 && in % 2 > 0) hprojCorr->SetLineStyle(2);
      	 hproj->SetName(Form("hCorrproj%d_%d",iaxis,in));
      	 hproj->SetTitle(Form(";%s",fhnMassResponseCorr->GetAxis(iaxis)->GetTitle()));
      	 cmassInt->cd(iaxis+1);
      	 gPad->SetLogy();
      	 hprojCorr->Draw("histsames");
      }
      
      for(Int_t ipTP=0; ipTP<npTPartBins; ipTP++){
      	 Double_t pTP[2]={pTPartMin + pTPartBinW*ipTP, pTPartBinW * (ipTP+1) + pTPartMin};
      	 TString pTPstring=Form("%.0f< #it{p}_{T,Part} < %.0f",pTP[0], pTP[1]);
      	 TPaveText *pv=new TPaveText(0.6,0.4,0.9,0.6,"NDC");
      	 pv->SetFillStyle(0);
      	 pv->SetBorderSize(0);
      	 pv->AddText(pTPstring);
      	 //pv->Draw();
      	 //pT detector level 
      	 TH1D* hpTUCorr=DoProjections(fhnMassResponse, axispTP, pTP, axispTD, Form("hpTDetU_pTP%d_%d",ipTP,in));
      	 //pT detector level corrected
      	 TH1D* hpTCorr=DoProjections(fhnMassResponseCorr, axispTP, pTP, axispTD, Form("hpTDetC_pTP%d_%d",ipTP,in));
      	 hpTCorr->SetLineColor(colors[ipTP]);
      	 if(in>0 && in % 2 > 0) hpTCorr->SetLineStyle(2);
      	 else {
      	    legpTPBins->AddEntry(hpTCorr, pTPstring, "L");
      	 }
      	 if(ipTP == 0){
      	    leg->AddEntry(hpTCorr, strLeg[in],"L");
      	 
      	 }
      	 
      	 cpT->cd();
      	 gPad->SetLogy();
      	 hpTCorr->Draw("sames");
      	 TH1F *hpTR = (TH1F*)hpTUCorr->Clone(Form("hpTRatio_pTP%d_%d",ipTP,in));
      	 hpTR->Divide(hpTCorr);
      	 hpTR->SetMarkerStyle(20);
      	 hpTR->SetLineColor(colors[ipTP]);
      	 hpTR->SetMarkerColor(colors[ipTP]);
      	 if(in>0 && in % 2 > 0) hpTCorr->SetLineStyle(2);
      	 
      	 cpTr->cd();
      	 hpTR->Draw("Psames");
      	 leg->Draw();
      	 //Mass detector level corrected
      	 TH1D* hMCorr=DoProjections(fhnMassResponseCorr, axispTP, pTP, axisMDet, Form("hMDet_pTP%d_%d",ipTP,in));
      	 hMCorr->SetLineColor(colors[ipTP]);
      	 if(in % 2 > 0) hMCorr->SetLineStyle(2);
      	 cM->cd();
      	 gPad->SetLogy();
      	 hMCorr->Draw("sames");
      	 
      	 //relative Mass difference (M_det - M_part)/M_part
      	 TH1D* hDMRelCorr=DoProjections(fhnDeltaMassCorr, axispTP2, pTP, axisDMRel, Form("hpTDet_pTP%d_%d",ipTP,in));
      	 hDMRelCorr->SetLineColor(colors[ipTP]);
      	 if(in % 2 > 0) hDMRelCorr->SetLineStyle(2);
      	 
      	 cDMrel->cd();
      	 gPad->SetLogy();
      	 hDMRelCorr->Draw("sames");
      	 
      	 //TH1D* hDMRelCorr=DoProjections(fhnDeltaMassCorr, axispTP, pTP, axisMDet, Form("hMDet_pTP%d",ipTP));
      	 //hMCorr->SetLineColor(colors[ipTP]);
      	 //if(in % 2 > 0) hMCorr->SetLineStyle(2);
      	 //cM->cd();
      	 //gPad->SetLogy();
      	 //hMCorr->Draw("sames");
      	 
      	 
      	 
      }
  }
  cpT->cd();
  legpTPBins->Draw();
  leg->Draw();
  cpTr->cd();
  legpTPBins->Draw();
  cM->cd();
  legpTPBins->Draw();
  leg->Draw();
  //cMr->cd();
  //legpTPBins->Draw();
  
  SaveCv(cpT);
  SaveCv(cM);
  
  cDMrel->cd();
  legpTPBins->Draw();
  SaveCv(cDMrel);
}

void RunCompareCorrected(){
   
   Bool_t sameListName=false;
   const Int_t n=2;
   TString path[n] = {"/data/Work/jets/JetMass/DetectorCorrections/FastSimulation/Eff08/output20k/", "/data/Work/jets/JetMass/DetectorCorrections/FastSimulation/Eff08/output20k/"};
   TString filename[n] = {"AnalysisResults.root", "AnalysisResults.root"};
   
   TString lname = "";
   TString suffixlname[n] = {"EpT", "P"};
   lname = "JetMassStructure_Jet_AKTChargedR040_PicoTracksFastSimu80_pT0000_E_scheme_TC"; //_pT0150
   TString nameLeg[n] = {"EpTdep", "Econst+Smearing"};
   
   TString strIn[n];
   TString strLst[n];
   for (Int_t j=0;j<n;j++){
      
      strIn[j]=path[j]+filename[j];
      if(sameListName) strLst[j]=lname;
      else strLst[j] = lname + suffixlname[j];
   }

   CompareCorrected(strIn, strLst, nameLeg, n);
   
   return;
}

void RunDrawTemplate(Double_t pTJMin=60, Double_t pTJMax=80){
   
   TString strIn = "/data/Work/jets/JetMass/DetectorCorrections/LHC12a15e/Train399/AnalysisResultsWeighted.root"; //need to figure out if this is at generated or reco level
   TString strLstIn = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC";
   TString strOut = "/data/Work/jets/JetMass/DetectorCorrections/LHC12a15e/Train410/outputpTHardBins/169838/JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC.root";
   TString strLstOut = "JetMassStructure_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC";
   DrawTemplate(pTJMin, pTJMax, strIn, strLstIn, strOut, strLstOut);

}
void DrawTemplate(Double_t pTJMin, Double_t pTJMax, TString strIn, TString strLstIn, TString strOut, TString strLstOut){ //"In" is before correction -input file, "out" is after the correction
   SetStyle(0);
   
   TPaveText *t=new TPaveText(0.3,0.8,0.8,0.9,"NDC");
   t->SetBorderSize(0);
   t->SetFillStyle(0);
   t->AddText("Scaled to GetEntries()");
   
   TString ptstring=Form("%.0f<p_{T,jet}<%.0f GeV/c",pTJMin,pTJMax);
     
   TPaveText *pvpTJ=new TPaveText(0.2,0.7,0.8,0.8,"NDC");
   pvpTJ->SetBorderSize(0);
   pvpTJ->SetFillStyle(0);
   pvpTJ->AddText(ptstring);
   
   TPaveText *ttitle=new TPaveText(0.3,0.9,0.8,0.99,"NDC");
   ttitle->SetBorderSize(0);
   ttitle->SetFillStyle(0);
   ttitle->AddText("Uncorrected");
   TPaveText *ttitleC=new TPaveText(0.3,0.9,0.8,0.99,"NDC");
   ttitleC->SetBorderSize(0);
   ttitleC->SetFillStyle(0);
   ttitleC->AddText("Corrected");
   
   TList *l=ReadFile(strIn, strLstIn);
   TH3D* fh3JetPtDRTrackPt = (TH3D*) l->FindObject(Form("fh3JetPtDRTrackPt_0")); //0 is the centrality
   fh3JetPtDRTrackPt->GetXaxis()->SetRangeUser(pTJMin,pTJMax);
   TH2D* hpTTR = (TH2D*) fh3JetPtDRTrackPt->Project3D("zy");
   hpTTR->SetName("hpTTR");
   hpTTR->Scale(hpTTR->GetEntries());
   //Double_t max=hpTTR->GetMaximum();
   TH1D* hpTT = (TH1D*) fh3JetPtDRTrackPt->Project3D("z");
   hpTT->SetName("hpTT");
   hpTT->Scale(hpTT->GetEntries());
   hpTT->SetMarkerStyle(20);

   TList *lc=ReadFile(strOut, strLstOut);
   TH3D* fh3JetPtDRTrackPtC = (TH3D*) lc->FindObject(Form("fh3JetPtDRTrackPt_0")); //0 is the centrality
   fh3JetPtDRTrackPtC->GetXaxis()->SetRangeUser(pTJMin,pTJMax);
   TH2D* hpTTRC = (TH2D*) fh3JetPtDRTrackPtC->Project3D("zy");
   hpTTRC->SetName("hpTTRC");
   hpTTRC->Scale(hpTTRC->GetEntries());
   //Double_t maxC=hpTTRC->GetMaximum();
   TH1D* hpTTC = (TH1D*) fh3JetPtDRTrackPtC->Project3D("z");
   hpTTC->SetName("hpTTC");
   hpTTC->Scale(hpTTC->GetEntries());
   hpTTC->SetMarkerStyle(24);
  
   //if(maxC > max) hpTTR->SetMaximum(maxC);
   //else hpTTRC->SetMaximum(max);
   TCanvas *cTemplate = new TCanvas("cTemplate", "Templates", 900,500);
   cTemplate->Divide(2,1);
   cTemplate->cd(1);
   //gPad->SetLogz();
   hpTTR->Draw("colz");
   t->Draw();
   pvpTJ->Draw();
   ttitle->Draw();
   cTemplate->cd(2);
   //gPad->SetLogz();
   hpTTRC->Draw("colz");
   t->Draw();
   pvpTJ->Draw();
   ttitleC->Draw();
   SaveCv(cTemplate);
   
   TCanvas *cpTT = new TCanvas("cpTT", "p_{T} distribution", 800, 800);
   cpTT->cd();
   gPad->SetLogy();
   hpTT->Draw("P");
   hpTTC->Draw("Psames");
}

void Nmissing(TString strIn, TString strLst1, TString strLst2 = "JetByJetCorrectionOutput"){

   TList* list1 = ReadFile(strIn, strLst1);
   TList* list=0x0;
   if(!list1){
      list = ReadFile(strIn, strLst2);
      if(!list){
      	 Printf("We have a problem");
      }
      list->ls();
   } else {
    list = (TList*)list1->FindObject(strLst2);
   }
   const Int_t npTjbins = 6;
   Double_t ptjrange[npTjbins+1] = {0,20.,40.,60.,80.,100.,120.};
   
   TString nameNadded = "fNmissing";//"fhNmissing"
   TH3F* hAvgNMissNaddptj = (TH3F*)list->FindObject(nameNadded);   
   
   
   if(!hAvgNMissNaddptj){
      Printf("%s not found",nameNadded.Data());
      return;
   }

   TH1D* hNadded = hAvgNMissNaddptj->ProjectionY("hNadded");
   hNadded->SetLineWidth(2);
   hNadded->SetLineColor(colors[1]);
   
   TH1D* hAvgNMiss = hAvgNMissNaddptj->ProjectionZ("hAvgNMiss");
   hAvgNMiss->SetLineWidth(2);
   hAvgNMiss->SetLineColor(colors[2]);
   
   TH2F* hNadded2d = (TH2F*)hAvgNMissNaddptj->Project3D("yx");

   TH2F* hAvgNMiss2d = (TH2F*)hAvgNMissNaddptj->Project3D("zx");

   TH2F* hAvgNMissNadded2d = (TH2F*)hAvgNMissNaddptj->Project3D("zy");
  
   Int_t nx, ny, dx, dy;
   CalculatePads(3, nx, ny, dx, dy);

   TCanvas *cNmissNadd= new TCanvas("cNmissNadd", "N missing and N added pT integrated",dx, dy);
   cNmissNadd->Divide(nx, ny);
   cNmissNadd->cd(1);
   hAvgNMiss2d->Draw();
   cNmissNadd->cd(2);
   hNadded2d->Draw();
   cNmissNadd->cd(3);
   hAvgNMiss->Draw();
   hNadded->Draw("sames");
   cNmissNadd->cd(4);
   hAvgNMissNadded2d->Draw("colz");
   SaveCv(cNmissNadd);
   
   
   CalculatePads(npTjbins, nx, ny, dx, dy);
   TCanvas *cNmissNaddpTj= new TCanvas("cNmissNaddpTj", "N missing and N added pT integrated",dx, dy);
   cNmissNaddpTj->Divide(nx,ny);
   
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      Printf("%.0f < pTj < %.0f GeV/c",ptjrange[iptj], ptjrange[iptj+1]);
      hNadded2d->GetXaxis()->SetRange(hNadded2d->GetXaxis()->FindBin(ptjrange[iptj]), hNadded2d->GetXaxis()->FindBin(ptjrange[iptj+1])-1);
      hAvgNMiss2d->GetXaxis()->SetRange(hAvgNMiss2d->GetXaxis()->FindBin(ptjrange[iptj]), hAvgNMiss2d->GetXaxis()->FindBin(ptjrange[iptj+1])-1);
      TH1D* hNaddedpTj = hNadded2d->ProjectionY(Form("hNaddedpTj%d",iptj));
      TH1D* hAvgNMisspTj = hAvgNMiss2d->ProjectionY(Form("hAvgNMisspTj%d",iptj));
      hAvgNMisspTj->SetLineWidth(2);
      hNaddedpTj->SetLineWidth(2);
      hNaddedpTj->SetLineColor(colors[1]);
      cNmissNaddpTj->cd(iptj+1);
      gPad->SetLogy();
      hNaddedpTj->Draw();
      hAvgNMisspTj->Draw("sames");
      	 
   }
   hNadded2d->GetXaxis()->SetRange(0,-1); //reset ptjet axis range
   hAvgNMiss2d->GetXaxis()->SetRange(0,-1); //reset ptjet axis range
   
   SaveCv(cNmissNaddpTj);

   TString nameNmisscorr = "fhCmpNmissStrategy";
   TH2F* hNMissAvgpTEffNTruth = (TH2F*)list->FindObject(nameNmisscorr);   
   
   
   if(!hNMissAvgpTEffNTruth || !hNMissAvgpTEffNTruth->GetEntries()){
      Printf("%s not found or empty",nameNmisscorr.Data());
      return;
   }
   
   TCanvas *cNmissCmp = new TCanvas("cNmissCmp", "N avg pT eff, N truth", 1300, 500);
   cNmissCmp->Divide(3,1);
   cNmissCmp->cd(1);
   hNMissAvgpTEffNTruth->Draw("colz");
   Printf("%s Axis X = %s, Axis Y = %s", hNMissAvgpTEffNTruth->GetName(), hNMissAvgpTEffNTruth->GetXaxis()->GetTitle(), hNMissAvgpTEffNTruth->GetYaxis()->GetTitle());
   
   Int_t nbinsX = hNMissAvgpTEffNTruth->GetXaxis()->GetNbins(), nbinsY = hNMissAvgpTEffNTruth->GetYaxis()->GetNbins();
   Int_t k=0;
   TLegend *legNmissAvg = new TLegend(0.4, 0.6, 0.9, 0.8, "N_{missed-Avg} = ");
   legNmissAvg->SetFillStyle(0);
   legNmissAvg->SetBorderSize(0);
   for(Int_t i=0; i < nbinsX; i++){
      //Printf("Bin %d content = %f", i+1, hNMissAvgpTEffNTruth->GetBinContent(i+1, 5));
      if(hNMissAvgpTEffNTruth->GetBinContent(i+1, 5) < 1e-6) continue;
      TH1D* hProjNY = (TH1D*)hNMissAvgpTEffNTruth->ProjectionY(Form("hProjNY%d", i), i+1, i+1);
      Printf("%d, maximum %.0f, mean = %.3f, sigma = %.3f", k,hProjNY->GetMaximum(),hProjNY->GetMean(),hProjNY->GetRMS());
      hProjNY->SetLineColor(colors[k]);
      hProjNY->SetFillColor(colors[k]);
      legNmissAvg->AddEntry(hProjNY, Form("%.0f, N_{Avg}: #mu = %.3f, #sigma = %.3f", hNMissAvgpTEffNTruth->GetXaxis()->GetBinLowEdge(i+1), hProjNY->GetMean(),hProjNY->GetRMS()), "LF");
      cNmissCmp->cd(2);
      gPad->SetLogy();
      if(k == 0 ) hProjNY->Draw();
      else hProjNY->Draw("sames");
      k++;   
   }
   cNmissCmp->cd(2); legNmissAvg->Draw();
   
   k=0;
   TLegend *legNmissTrth = new TLegend(0.4, 0.3, 0.9, 0.8, "N_{missed-Truth} = ");
   legNmissTrth->SetFillStyle(0);
   legNmissTrth->SetBorderSize(0);
   for(Int_t i=0; i < nbinsX; i++){
      //Printf("Bin %d content = %f", i+1, hNMissAvgpTEffNTruth->GetBinContent(i+1, 5));
      if(hNMissAvgpTEffNTruth->GetBinContent(1, i+1) < 1e-6) continue;
      TH1D* hProjNY = (TH1D*)hNMissAvgpTEffNTruth->ProjectionX(Form("hProjNX%d", i), i+1, i+1);
      Printf("%d, maximum %.0f, mean = %.3f, sigma = %.3f", k,hProjNY->GetMaximum(),hProjNY->GetMean(),hProjNY->GetRMS());
      hProjNY->SetLineColor(colors[k]);
      hProjNY->SetFillColor(colors[k]);
      legNmissTrth->AddEntry(hProjNY, Form("%.0f, N_{truth}: #mu = %.3f, #sigma = %.3f", hNMissAvgpTEffNTruth->GetXaxis()->GetBinLowEdge(i+1), hProjNY->GetMean(),hProjNY->GetRMS()), "LF");
      cNmissCmp->cd(3);
      gPad->SetLogy();
      if(k == 0 ) hProjNY->Draw();
      else hProjNY->Draw("sames");
      k++;   
   }
   cNmissCmp->cd(3); legNmissTrth->Draw();
   
   SaveCv(cNmissCmp);
}

void NmissingOld(TString strIn, TString strLst1, TString strLst2 = "JetByJetCorrectionOutput"){

   TList* list1 = ReadFile(strIn, strLst1);
   TList* list = (TList*)list1->FindObject(strLst2); 
   const Int_t npTjbins = 6;
   Double_t ptjrange[npTjbins+1] = {0,20.,40.,60.,80.,100.,120.};
   
   TString nameAvgNMiss = "fAvgNmiss", nameNadded = "fNmissing";
   TH2F* hAvgNMiss2d = (TH2F*)list->FindObject(nameAvgNMiss);   
   
   
   if(!hAvgNMiss2d){
      Printf("%s not found",nameAvgNMiss.Data());
      return;
   }
   
   TH1D* hAvgNMiss = hAvgNMiss2d->ProjectionY("hAvgNMiss");
   hAvgNMiss->SetLineWidth(2);
   
   TH2F* hNadded2d = (TH2F*)list->FindObject(nameNadded);

   if(!hNadded2d){
      Printf("%s not found",nameNadded.Data());
      return;
   }

   TH1D* hNadded = hNadded2d->ProjectionY("hNadded");
   hNadded->SetLineWidth(2);
   hNadded->SetLineColor(colors[1]);
   
   Int_t nx, ny, dx, dy;
   CalculatePads(3, nx, ny, dx, dy);

   TCanvas *cNmissNadd= new TCanvas("cNmissNadd", "N missing and N added pT integrated",dx, dy);
   cNmissNadd->Divide(nx, ny);
   cNmissNadd->cd(1);
   hAvgNMiss2d->Draw();
   cNmissNadd->cd(2);
   hNadded2d->Draw();
   cNmissNadd->cd(3);
   hAvgNMiss->Draw();
   hNadded->Draw("sames");
   SaveCv(cNmissNadd);
   
   
   CalculatePads(npTjbins, nx, ny, dx, dy);
   TCanvas *cNmissNaddpTj= new TCanvas("cNmissNaddpTj", "N missing and N added pT integrated",dx, dy);
   cNmissNaddpTj->Divide(nx,ny);
   
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      Printf("%.0f < pTj < %.0f GeV/c",ptjrange[iptj], ptjrange[iptj+1]);
      hNadded->GetXaxis()->SetRange(hNadded->GetXaxis()->FindBin(ptjrange[iptj]), hNadded->GetXaxis()->FindBin(ptjrange[iptj+1])-1);
      hAvgNMiss2d->GetXaxis()->SetRange(hAvgNMiss2d->GetXaxis()->FindBin(ptjrange[iptj]), hAvgNMiss2d->GetXaxis()->FindBin(ptjrange[iptj+1])-1);
      TH1D* hNaddedpTj = hNadded2d->ProjectionY(Form("hNaddedpTj%d",iptj));
      TH1D* hAvgNMisspTj = hAvgNMiss2d->ProjectionY(Form("hAvgNMisspTj%d",iptj));
      hAvgNMisspTj->SetLineWidth(2);
      hNaddedpTj->SetLineWidth(2);
      hNaddedpTj->SetLineColor(colors[1]);
      cNmissNaddpTj->cd(iptj+1);
      gPad->SetLogy();
      hNaddedpTj->Draw();
      hAvgNMisspTj->Draw("sames");
      
      
      	 
   }
   
   SaveCv(cNmissNaddpTj);

}

void ConstituentpTDistribution(TString strIn, TString strLst1){
   const Int_t pTjbin = 4;
   //Double_t pTjlims[pTjbin+1] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90, 100.};
   Double_t pTjlims[pTjbin+1] = {20.,40.,60.,80., 100.};
   
   const Int_t avgpTtbin = 4;
   Double_t avgpTt[avgpTtbin+1] = {0.,4., 6., 10., 40.};
   TList* list = ReadFile(strIn, strLst1);
   if(!list) {
      Printf("List not found");
      return;
   }
   
   TH2F *h2Rec = (TH2F*)list->FindObject("fhAllpTRec");
   TH2F *h2Gen = (TH2F*)list->FindObject("fhAllpTGen");
   TH2F *h2Cor = (TH2F*)list->FindObject("fhAllpTCor");
   if(!h2Rec || !h2Gen) {
      Printf("Histograms not found");
      return;
   }
   if(!h2Cor){
      Printf("Maybe older version of the macro without fhAllpTCor");
      return;
   }
   THnSparseD *hNconstDiff = (THnSparseD*)list->FindObject("fhConstRecGen");
   if(!hNconstDiff){
      Printf("THnSparse not found");
      return;
   }
   
   TCanvas *cpTconst = new TCanvas("cpTconst", "pT track distributions in pT jet bins (need to refine the def)", 1000,1000);
   //cpTconst->Divide(3,4); //good with 10 bins
   cpTconst->Divide(2,2); //4 bins
   TCanvas *cNconstDiff = new TCanvas ("cNconstDiff", "Constituent difference", 1000, 1200);
   //cNconstDiff->Divide(4,5); //good with 10 bins
   cNconstDiff->Divide(2,3); // 3 bins cut the last bin 80-100
   TCanvas *cNconstEff = new TCanvas ("cNconstEff", "Constituent efficiency", 1000, 1200);
   //cNconstEff->Divide(4,5); //good with 10 bins
   cNconstEff->Divide(2,3); // 3 bins cut the last bin 80-100
   TCanvas *cNconstEffproj = new TCanvas ("cNconstEffproj", "Constituent efficiency projection", 1000, 1200);
   //cNconstEff->Divide(4,5); //good with 10 bins
   cNconstEffproj->Divide(2,3); // 3 bins cut the last bin 80-100
   TCanvas *cNconstEffvspTt = new TCanvas ("cNconstEffvspTt", "Constituent efficiency vs pT track", 800, 800);
   
   Printf("X axis is %s, Y axis is %s", h2Rec->GetXaxis()->GetTitle(), h2Rec->GetYaxis()->GetTitle());
   
   TH1D* hpTjRec = (TH1D*)h2Rec->ProjectionY("hpTjRec");
   TH1D* hpTjGen = (TH1D*)h2Gen->ProjectionY("hpTjGen");
   
   TH1D* hpTGenOverRec = (TH1D*)hpTjGen->Clone("hpTGenOverRec");
   hpTGenOverRec->Divide(hpTjRec);
   
   TCanvas *cRatios = new TCanvas("cRatios", "Generated/Reconstructed -> Correction factor", 1000, 1000);
   cRatios->Divide(2,3);
   cRatios->cd(1);
   hpTGenOverRec->Draw();
   
   TPaveText *pvgen = new TPaveText(0.4,0.4, 0.55,0.6, "NDC");
   pvgen->SetFillStyle(0);
   TPaveText *pvrec = (TPaveText*)pvgen->Clone();
   pvrec->AddText("REC"); pvgen->AddText("GEN");
   TLegend *legpT = new TLegend(0.2, 0.7, 0.6, 0.9);
   legpT->SetBorderSize(0);
   legpT->SetFillStyle(0);
   
   TLegend *legpTj = new TLegend(0.2, 0.2, 0.6, 0.4);
   legpTj->SetBorderSize(0);
   legpTj->SetFillStyle(0);
   
   //read templates
   TString listJbJname = "JetByJetCorrectionOutput";
   TList* listJbJ = (TList*)list->FindObject(listJbJname);
   if(!listJbJ) {
      Printf("List JbJ not found");
      
   }
   Int_t nTemplHist = 30;
   TString hTname = "hPtR_";
   Double_t DeltaPtJ = 5;
   TH2D *hTemplate[nTemplHist];
   for(Int_t iT = 0; iT<nTemplHist; iT++){
      hTemplate[iT]=0x0;
      if(!listJbJ) continue;
      Double_t pTJmin = iT * DeltaPtJ, pTJmax = (iT+1) * DeltaPtJ;
      hTemplate[iT] = (TH2D*)listJbJ->FindObject(Form("%s%.0f_%.0f", hTname.Data(), pTJmin, pTJmax));
   
   }
   TH2D *hTemplate2[pTjbin];
   
   for(Int_t i = 0; i < pTjbin; i++){
      //sum Template in Template2
      Int_t bintosummin = GetJetPtBinTemp(pTjlims[i], DeltaPtJ), bintosummax = GetJetPtBinTemp(pTjlims[i+1], DeltaPtJ);
      hTemplate2[i] = (TH2D*)hTemplate[bintosummin]->Clone(Form("%s%.0f_%.0f", hTname.Data(), pTjlims[i], pTjlims[i+1]));
      for(Int_t iT = bintosummin+1 ; iT < bintosummax; iT++){
      	 hTemplate2[i]->Add(hTemplate[iT]);
      }
      
      TPaveText *pvpTj = new TPaveText(0.7, 0.3, 0.9, 0.4, "NDC");
      pvpTj->SetFillStyle(0);
      pvpTj->AddText(Form("pT jet %.0f - %.0f", pTjlims[i], pTjlims[i+1]));
      Printf("pT jet %.0f - %.0f", pTjlims[i], pTjlims[i+1]);
      Int_t b1 = h2Rec->GetYaxis()->FindBin(pTjlims[i]), b2 = h2Rec->GetYaxis()->FindBin(pTjlims[i+1])-1;
      Printf("Select bin %d, %d",  b1, b2);
      hpTjRec->GetXaxis()->SetRange(b1,b2);
      TPaveText pv(0.4,0.4, 0.65,0.6, "NDC");
      pv.SetBorderSize(0);
      pv.SetFillStyle(0);
      pv.AddText(Form("Mean Rec = %.2f", hpTjRec->GetMean()));
      hpTjGen->GetXaxis()->SetRange(b1,b2);
      pv.AddText(Form("Mean Gen = %.2f", hpTjGen->GetMean()));
      
      TH1D* hpTtrackRec = (TH1D*)h2Rec->ProjectionX(Form("hpTtrackRec%d",i), b1, b2);
      b1 = h2Gen->GetYaxis()->FindBin(pTjlims[i]); b2 = h2Gen->GetYaxis()->FindBin(pTjlims[i+1])-1;
      TH1D* hpTtrackGen = (TH1D*)h2Gen->ProjectionX(Form("hpTtrackGen%d",i), b1, b2);
      TH1D* hpTtrackCor = (TH1D*)h2Cor->ProjectionX(Form("hpTtrackCor%d",i), b1, b2);
      TH1D* hpTtrackGenOverRec = (TH1D*)hpTtrackGen->Clone(Form("hpTtrackGenOverRec%d",i));
      hpTtrackGenOverRec->Divide(hpTtrackRec);
            
      hpTtrackGen->SetMarkerStyle(20);
      hpTtrackRec->SetMarkerStyle(21);
      hpTtrackCor->SetMarkerStyle(25);
      hpTtrackGenOverRec->SetMarkerStyle(21);
      hpTtrackGen->SetMarkerColor(kBlue+3);
      hpTtrackRec->SetMarkerColor(kMagenta);
      hpTtrackCor->SetMarkerColor(kMagenta);
      hpTtrackGenOverRec->SetMarkerColor(kMagenta);
      hpTtrackGenOverRec->SetLineColor(kMagenta);
      
      
      TH1D* hpTtrackRecCor = (TH1D*)hpTtrackRec->Clone("hpTtrackRecCor");
      hpTtrackRecCor->Add(hpTtrackCor);
      hpTtrackRecCor->SetMarkerColor(kMagenta-2);
      
      TH1D* hpTtrackGenOverRecCor = (TH1D*)hpTtrackGen->Clone(Form("hpTtrackGenOverRecCor%d",i));
      hpTtrackGenOverRecCor->Divide(hpTtrackRecCor);

      hpTtrackGenOverRecCor->SetMarkerStyle(25);
      hpTtrackGenOverRecCor->SetMarkerColor(kMagenta-2);
      hpTtrackGenOverRecCor->SetLineColor(kMagenta-2);
      if(i==0) {
      	 legpT->AddEntry(hpTtrackGen, "Generated", "PL");
      	 legpT->AddEntry(hpTtrackRec, "Reconstructed", "PL");
      	 legpT->AddEntry(hpTtrackCor, "Added in correction", "PL");
      	 legpT->AddEntry(hpTtrackRecCor, "Reco+Added", "PL");
      }
      
      cpTconst->cd(i+1);
      gPad->SetLogy();
      hpTtrackGen->Draw();
      hpTtrackRec->Draw("sames");
      hpTtrackCor->Draw("Psames");
      hpTtrackRecCor->Draw("Psames");
      hTemplate2[i]->ProjectionY()->Draw("sames");
      legpT->Draw();
      pv.DrawClone();
      pvpTj->Draw();
      cRatios->cd(i+2);
      hpTtrackGenOverRec->Draw("Psames");
      hpTtrackGenOverRecCor->Draw("Psames");
      pv.DrawClone();
      pvpTj->Draw();
      
      if(pTjlims[i] > 60) {
      	 Printf("Exclude last bin for thsparse plot");
      	 break; 
      }
      Printf("Projection %s vs %s or %s", hNconstDiff->GetAxis(0)->GetTitle(), hNconstDiff->GetAxis(1)->GetTitle(), hNconstDiff->GetAxis(2)->GetTitle());
      Int_t pTjbinsel[2] = {hNconstDiff->GetAxis(3)->FindBin(pTjlims[i]), hNconstDiff->GetAxis(3)->FindBin( pTjlims[i+1])};
      hNconstDiff->GetAxis(3)->SetRange(pTjbinsel[0], pTjbinsel[1]);
      //projection Nconstdiff vs average pT,track rec 
      TH2D* hNconstdiffvspTRec = hNconstDiff->Projection(0,1);
      hNconstdiffvspTRec->SetName(Form("hNconstdiffvspTRec%d",i));
      //projection Nconstdiff vs average pT,track gen 
      TH2D* hNconstdiffvspTGen = hNconstDiff->Projection(0,2);
      hNconstdiffvspTGen->SetName(Form("hNconstdiffvspTGen%d",i));
      cNconstDiff->cd(2*i+1);
      hNconstdiffvspTRec->Draw("colz");
      pvrec->Draw();
      pvpTj->Draw();
      cNconstDiff->cd(2*(i+1));
      hNconstdiffvspTGen->Draw("colz");
      pvgen->Draw();
      pvpTj->Draw();
      
      Printf("Projection %s vs %s or %s", hNconstDiff->GetAxis(5)->GetTitle(), hNconstDiff->GetAxis(1)->GetTitle(), hNconstDiff->GetAxis(2)->GetTitle());
      //projection Efficiency vs average pT,track rec
      TH2D* hNconsteffvspTRec = hNconstDiff->Projection(5,1);
      hNconsteffvspTRec->SetName(Form("hNconsteffvspTRec%d",i));
      TPaveText *pvmeanEffRec = new TPaveText(0.2, 0.2, 0.4, 0.45, "NDC");
      pvmeanEffRec->SetFillStyle(0);
      pvmeanEffRec->SetBorderSize(0);
      TPaveText *pvmeanEffGen = (TPaveText*)pvmeanEffRec->Clone();
      pvmeanEffRec->AddText(Form("Mean Eff = %.2f", hNconsteffvspTRec->GetMean(2)));
      pvmeanEffRec->AddText(Form("Sigma Eff = %.2f", hNconsteffvspTRec->GetRMS(2)));
      //projection Efficiency vs average pT,track gen 
      TH2D* hNconsteffvspTGen = hNconstDiff->Projection(5,2);
      hNconsteffvspTGen->SetName(Form("hNconsteffvspTGen%d",i));
      pvmeanEffGen->AddText(Form("Mean Eff = %.2f", hNconsteffvspTGen->GetMean(2)));
      pvmeanEffGen->AddText(Form("Sigma Eff = %.2f", hNconsteffvspTGen->GetRMS(2)));
      cNconstEff->cd(2*i+1);
      hNconsteffvspTRec->Draw("colz");
      pvrec->Draw();
      pvpTj->Draw();
      pvmeanEffRec->Draw();
      cNconstEff->cd(2*(i+1));
      hNconsteffvspTGen->Draw("colz");
      pvgen->Draw();
      pvpTj->Draw();
      pvmeanEffGen->Draw();
      
      //per each bin in ptj slice in pttrack and draw the distribution of the efficiency and its mean value
      TGraphErrors *gmeanEffRecvspTt = new TGraphErrors(avgpTtbin);
      gmeanEffRecvspTt->SetName(Form("gmeanEffRecvspTt-pTj%d",i));
      gmeanEffRecvspTt->SetTitle(Form("Mean Eff Rec (pTj %d); #it{p}_{T,track};Mean Eff", i));
      gmeanEffRecvspTt->SetMarkerColor(colors[i]);
      gmeanEffRecvspTt->SetLineColor(colors[i]);
      gmeanEffRecvspTt->SetMarkerStyle(20);
      TGraphErrors *gmeanEffGenvspTt = (TGraphErrors*)gmeanEffRecvspTt->Clone("");
      gmeanEffGenvspTt->SetMarkerStyle(21);
      legpTj->AddEntry(gmeanEffRecvspTt, Form("%.0f < #it{p}_{T,jet} < %.0f ", pTjlims[i], pTjlims[i+1]), "P");
      	 
      for(Int_t j = 0; j<avgpTtbin ; j++){
      	 Int_t bt1 = hNconsteffvspTRec->GetXaxis()->FindBin(avgpTt[j]), bt2 = hNconsteffvspTRec->GetXaxis()->FindBin(avgpTt[j+1])-1;
      	 TH1D *hNconsteffRec = hNconsteffvspTRec->ProjectionY(Form("hNconsteffRec-ptj%d-ptt%d",i, j), bt1, bt2);
      	 hNconsteffRec->Rebin(3);
      	 TH1D *hNconsteffGen = hNconsteffvspTGen->ProjectionY(Form("hNconsteffGen-ptj%d-ptt%d",i, j), bt1, bt2);
      	 hNconsteffGen->Rebin(3);
      	 
      	 //TProfile *hNconsteffRec = hNconsteffvspTRec->ProfileY(Form("hNconsteffRec-ptj%d-ptt%d",i, j), bt1, bt2);
      	 //
      	 //TProfile *hNconsteffGen = hNconsteffvspTGen->ProfileY(Form("hNconsteffGen-ptj%d-ptt%d",i, j), bt1, bt2);
      	 //
      	 hNconsteffRec->SetMarkerColor(colors[j]);
      	 hNconsteffRec->SetMarkerStyle(21);
      	 hNconsteffGen->SetMarkerColor(colors[j]);
      	 hNconsteffGen->SetMarkerStyle(20);
      	 cNconstEffproj->cd(2*i+1);
      	 hNconsteffRec->Draw("Psames");
      	 pvrec->Draw();
      	 pvpTj->Draw();
      	 cNconstEffproj->cd(2*(i+1));
      	 hNconsteffGen->Draw("Psames");
      	 pvgen->Draw();
      	 pvpTj->Draw();
      	 
      	 gmeanEffRecvspTt->SetPoint(j,avgpTt[j] + (avgpTt[j+1]-avgpTt[j])*0.5, hNconsteffRec->GetMean());
      	 gmeanEffGenvspTt->SetPoint(j, avgpTt[j] + (avgpTt[j+1]-avgpTt[j])*0.5, hNconsteffGen->GetMean());
      	 gmeanEffRecvspTt->SetPointError(j, (avgpTt[j+1]-avgpTt[j])*0.5, hNconsteffRec->GetRMS());
      	 gmeanEffGenvspTt->SetPointError(j, (avgpTt[j+1]-avgpTt[j])*0.5, hNconsteffGen->GetRMS());
      }
      cNconstEffvspTt->cd();
      if(i==0) {
      	 gmeanEffRecvspTt->GetYaxis()->SetRangeUser(0.5,0.9);
      	 gmeanEffRecvspTt->Draw("AP");
      }
      else gmeanEffRecvspTt->Draw("P");
      gmeanEffGenvspTt->Draw("P");
      SaveCv(cNconstEffproj);
   }
   cNconstEffvspTt->cd();
   legpTj->Draw();
   
   
   TCanvas *cpT20const = new TCanvas("cpT20const", "pT track distributions for pTj > 20 GeV/c", 800, 800);
   
   Int_t b201 = h2Rec->GetYaxis()->FindBin(20), b202 = h2Rec->GetYaxis()->FindBin(pTjlims[pTjbin]) -1;
   Printf("Projection pTj 20 to %f, bins %d, %d", pTjlims[10], b201, b202);
   TH1F* hpTtrackRec20 = (TH1F*)h2Rec->ProjectionX("hpTtrackRec20", b201, b202);
   b201 = h2Gen->GetYaxis()->FindBin(20); b202 = h2Gen->GetYaxis()->FindBin(pTjlims[pTjbin]) -1;
   TH1F* hpTtrackGen20 = (TH1F*)h2Gen->ProjectionX("hpTtrackGen20", b201, b202);
   b201 = h2Cor->GetYaxis()->FindBin(20); b202 = h2Cor->GetYaxis()->FindBin(pTjlims[pTjbin]) -1;
   TH1F* hpTtrackCor20 = (TH1F*)h2Cor->ProjectionX("hpTtrackCor20", b201, b202);
   
   hpTtrackGen20->SetMarkerStyle(20);
   hpTtrackRec20->SetMarkerStyle(21);
   hpTtrackCor20->SetMarkerStyle(25);
   hpTtrackGen20->SetMarkerColor(kBlue+3);
   hpTtrackRec20->SetMarkerColor(kMagenta);
   hpTtrackCor20->SetMarkerColor(kMagenta);
   TH1F* hpTtrackRecCor20 = (TH1F*)hpTtrackRec20->Clone("hpTtrackRecCor20");
   hpTtrackRecCor20->Add(hpTtrackCor20);
   cpT20const->cd(1);
   gPad->SetLogy();
   hpTtrackGen20->Draw();
   hpTtrackRec20->Draw("sames");
   hpTtrackCor20->Draw("sames");
   hpTtrackRecCor20->Draw("sames");
   legpT->Draw();
   
   SaveCv(cpTconst);
   SaveCv(cRatios);
   SaveCv(cpT20const);
   SaveCv(cNconstDiff);
   SaveCv(cNconstEff);
   SaveCv(cNconstEffvspTt);

   Double_t pTjsel[2] = {60,80};
}

void SaveUnfoldingMatrix(TString strIn, TString strLstIn, TString suffixout, Int_t pTSel=2){
   // pTSel = 0 -> particle level
   // pTSel = 1 -> detector level
   // pTSel = 2 -> particle level for mass Part, detector level for mass detector
   
   if(pTSel == 0) suffixout+="pTP";
   if(pTSel == 1) suffixout+="pTD";
   if(pTSel == 2) suffixout+="pTPD";
   
   SetStyle(0);
   TList *lst=ReadFile(strIn, strLstIn);
   if(!lst) return;
   THnSparse *fhnMassResponse = (THnSparse*)lst->FindObject("fhnMassResponse");
   if(!fhnMassResponse){
      Printf("THnSparse not found");
      return;
   }
   TString titlemissing = "#it{p}_{T,part}"; //this is ad hoc, hopefully not needed for ever
   for(Int_t i=0; i<fhnMassResponse->GetNdimensions();i++){
      TString title = fhnMassResponse->GetAxis(i)->GetTitle();
      if(title.IsNull()) {
      	 Printf("Warning: axis %d has no title, I think it is %s", i, titlemissing.Data());
      	 title = titlemissing;
      	 fhnMassResponse->GetAxis(i)->SetTitle(titlemissing);
      }
      Printf("Axis %d is %s", i, title.Data());
      fhnMassResponse->GetAxis(i)->SetRange(0,-1);
   }
   TH2D *hpT   = (TH2D*)fhnMassResponse->Projection(3,2);  //ptpart (y) vs ptdet (x)
   TH2D *hm    = (TH2D*)fhnMassResponse->Projection(1,0);  //masspart (y) vs massdet (x)
   TH2D *hpTmP = (TH2D*)fhnMassResponse->Projection(3,1);  //ptpart (y) vs masspart (x)
   TH2D *hpTmD = (TH2D*)fhnMassResponse->Projection(2,0);  //ptdet (y) vs massdet (x)
   hpT   ->SetAxisRange(0.,1.e-5,"Z");;
   hm    ->SetAxisRange(0.,1.e-5,"Z");;
   hpTmP ->SetAxisRange(0.,1.e-5,"Z");;
   hpTmD ->SetAxisRange(0.,1.e-5,"Z");;   
   TCanvas *c = new TCanvas("c", "Saving...", 1000, 1000);
   c->Divide(2,2);
   c->cd(1);
   hpT   ->Draw("colz");
   c->cd(2);
   hm    ->Draw("colz");
   c->cd(3);
   hpTmP ->Draw("colz");
   c->cd(4);
   hpTmD ->Draw("colz");
   SaveCv(c);
   TString outname = Form("UnfoldingMatrix%s.root",suffixout.Data());
   TFile *fout = new TFile(outname, "recreate");
   fout->cd();
   hpT   ->Write();
   hm    ->Write();
   hpTmP ->Write();
   hpTmD ->Write();
   
   fhnMassResponse->Write();
   
   const Int_t npTjbins = 5;
   Double_t pTjbins[npTjbins+1] = {20., 40., 60., 80, 100., 120.};

   TCanvas *c1d[npTjbins];
   TCanvas *c1dlog[npTjbins];
   TCanvas *c1dRatio[npTjbins];
   
   Int_t nx, ny, dx, dy;
   CalculatePads(npTjbins, nx, ny, dx, dy);
   TCanvas *cRelUncRatio = new TCanvas("cRelUncRatio", "Rel Uncertainty on the Ratio", dy, dy);
   cRelUncRatio->Divide(nx,ny);
   TCanvas *cFinalRatio = new TCanvas("cFinalRatio", "Ratio used in the correction", dy, dy);
   cFinalRatio->Divide(nx,ny);
   TGraphErrors *gMeanMassDetc = new TGraphErrors(npTjbins);
   gMeanMassDetc->SetName("gMeanMassDetc");
   gMeanMassDetc->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
   gMeanMassDetc->SetMarkerStyle(25);
   TGraphErrors *gMeanMassPart = new TGraphErrors(npTjbins);
   gMeanMassPart->SetName("gMeanMassPart");
   gMeanMassPart->SetTitle(";#it{p}_{T,jet};#LT#it{M}_{jet}#GT");
   gMeanMassPart->SetMarkerStyle(25);
   
   for(Int_t iptj = 0; iptj < npTjbins; iptj++){
      c1d[iptj] = new TCanvas(Form("c1d%d",iptj), Form("Projection pT %.0f - %.0f", pTjbins[iptj], pTjbins[iptj+1]));
      c1dlog[iptj] = new TCanvas(Form("c1dlog%d",iptj), Form("Projection pT %.0f - %.0f, logy scale", pTjbins[iptj], pTjbins[iptj+1]));
      c1dRatio[iptj] = new TCanvas(Form("c1dRatio%d",iptj), Form("Ratio Part/Detc pT %.0f - %.0f", pTjbins[iptj], pTjbins[iptj+1]));
      
      TLegend *leg = new TLegend(0.6,0.5, 0.8, 0.8);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      
      //Particle level (pT = 3, Mass = 1)
      Int_t pTjbinselP[2] = {fhnMassResponse->GetAxis(3)->FindBin(pTjbins[iptj]), fhnMassResponse->GetAxis(3)->FindBin(pTjbins[iptj+1])-1};
      
      Int_t pTjbinselD[2] = {fhnMassResponse->GetAxis(2)->FindBin(pTjbins[iptj]), fhnMassResponse->GetAxis(2)->FindBin(pTjbins[iptj+1])-1};
      
      if(pTSel == 0){
      	 fhnMassResponse->GetAxis(3)->SetRange(pTjbinselP[0], pTjbinselP[1]); //pT Part?
      }
      if(pTSel == 1){
      	 fhnMassResponse->GetAxis(2)->SetRange(pTjbinselD[0], pTjbinselD[1]); //pT Det?
      }
      if(pTSel == 2){
      	 fhnMassResponse->GetAxis(3)->SetRange(pTjbinselP[0], pTjbinselP[1]); //pT Part?
      	 fhnMassResponse->GetAxis(2)->SetRange(pTjbinselD[0], pTjbinselD[1]); //pT Det?
      }
      
      TH1D* hMassPart = fhnMassResponse->Projection(1);
      hMassPart->SetName(Form("hMassPart_pTjPart%d", iptj));
      hMassPart->SetLineColor(colors[iptj]);
      hMassPart->SetMarkerColor(colors[iptj]);
      hMassPart->SetMarkerStyle(20);
      hMassPart->Scale(1./hMassPart->Integral());
      leg->AddEntry(hMassPart, "Part", "PL");
      
      fhnMassResponse->GetAxis(3)->SetRange(0, -1);
      c1d[iptj] ->cd();
      hMassPart->Draw("P");
      hMassPart->Write();
      c1dlog[iptj]->cd();
      gPad->SetLogy();
      hMassPart->Draw("P");
      
      gMeanMassPart->SetPoint(iptj, (pTjbins[iptj+1] - pTjbins[iptj])*0.5 + pTjbins[iptj], hMassPart->GetMean());
      gMeanMassPart->SetPointError(iptj, (pTjbins[iptj+1] - pTjbins[iptj])*0.5, hMassPart->GetMeanError());

      //Detector level (pT = 2, Mass = 0)
      
      TH1D* hMassDetc = fhnMassResponse->Projection(0);
      hMassDetc->SetName(Form("hMassDetc_pTjDetc%d", iptj));
      hMassDetc->SetLineColor(colors[iptj]);
      hMassDetc->SetMarkerColor(colors[iptj]);
      hMassDetc->SetMarkerStyle(24);
      hMassDetc->Scale(1./hMassDetc->Integral());
      leg->AddEntry(hMassDetc, "Detc", "PL");
      
      fhnMassResponse->GetAxis(2)->SetRange(0, -1);
      
      gMeanMassDetc->SetPoint(iptj, (pTjbins[iptj+1] - pTjbins[iptj])*0.5 + pTjbins[iptj], hMassDetc->GetMean());
      gMeanMassDetc->SetPointError(iptj, (pTjbins[iptj+1] - pTjbins[iptj])*0.5, hMassDetc->GetMeanError());

      c1d[iptj] ->cd();
      hMassDetc->Draw("Psames");
      leg->Draw();
      c1dlog[iptj]->cd();
      gPad->SetLogy();
      hMassDetc->Draw("Psames");
      leg->Draw();
      
      SaveCv(c1d[iptj]);
      SaveCv(c1dlog[iptj]);
      
      fout->cd();
      hMassDetc->Write();
      
      
      TH1D* hRatio = (TH1D*)hMassPart->Clone(Form("hMassPartToDetc_pTj%d", iptj));
      hRatio->Divide(hMassDetc);
      c1dRatio[iptj]->cd();
      hRatio->Draw("P");
      
      fout->cd();
      hRatio->Write();
      
      TH1D* hRelUncRatio = (TH1D*)hRatio->Clone(Form("hRelUncRatio_pTj%d", iptj));
      hRelUncRatio->SetMarkerStyle(24);
      Int_t binRelUncX=1; //bin corresponding at relative uncertainty X or grater
      Double_t X = 0.4; //to be decided
      for(Int_t k=0; k<hRelUncRatio->GetNbinsX(); k++){
      	 Double_t content = hRelUncRatio->GetBinContent(k+1), error = hRelUncRatio->GetBinError(k+1), newcontent = 0;
      	 if (content>0) newcontent = error/content;
      	 if ((binRelUncX == 1) && (newcontent > X) && (k>3)) {
      	    if((k+1) % 2 == 0) binRelUncX = k+1;
      	    else binRelUncX = k;
      	    Printf("Check k = %d, bin chosen = %d, content = %.1f", k, binRelUncX, newcontent);
      	 }
      	 hRelUncRatio->SetBinContent(k+1, newcontent );
      	 hRelUncRatio->SetBinError(k+1,0);
      }
      Printf("Stop at bin %d", binRelUncX);
      cRelUncRatio->cd(iptj+1);
      hRelUncRatio->Draw("P");
      SaveCv(cRelUncRatio);
      //remake fraction histograms
      TH1D *hRatioLessX = new TH1D(Form("hMassPartToDetcLess%.0f_pTj%d",X*100, iptj), Form("Fraction (with relative error < %.2f); #it{M}_{part}; Ratio", X), binRelUncX, hRatio->GetBinLowEdge(1), hRatio->GetBinLowEdge(binRelUncX+1));
      hRatioLessX->SetMarkerStyle(21);
      hRatioLessX->SetMarkerColor(colors[iptj]);
      hRatioLessX->Sumw2();
      for(Int_t k=0; k<binRelUncX; k++){
      	 hRatioLessX->SetBinContent(k+1, hRatio->GetBinContent(k+1));
      	 hRatioLessX->SetBinError(k+1, hRatio->GetBinError(k+1));
      }
      c1dRatio[iptj]->cd();
      hRatioLessX->Draw("sames");
      SaveCv(c1dRatio[iptj]);
      
      cFinalRatio->cd(iptj+1);
      hRatioLessX->Draw("sames");
      
      fout->cd();
      hRelUncRatio->Write();
      hRatioLessX->Write();
   }
   
   fout->cd();
   gMeanMassPart->Write();
   gMeanMassDetc->Write();
   SaveCv(cFinalRatio);
   Printf("Saved file %s",outname.Data());

}



TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Int_t* pTPBin, Int_t axisP, char* hprjname){
   fhnDeltaMass->    GetAxis(axisR)->SetRange(pTPBin[0],pTPBin[1]); //pT det
   TH1D* hDMpTP=(TH1D*)fhnDeltaMass->Projection(axisP);
   hDMpTP->SetName(hprjname);
   hDMpTP->SetLineWidth(2);
   return hDMpTP;
}

TH1D* DoProjections(THnSparse *fhnDeltaMass, Int_t axisR, Double_t* pTP, Int_t axisP, char* hprjname){
   Int_t pTPBin[2] = {
   fhnDeltaMass->GetAxis(axisR)->FindBin(pTP[0]), fhnDeltaMass->GetAxis(axisR)->FindBin(pTP[1])};
   fhnDeltaMass->    GetAxis(axisR)->SetRange(pTPBin[0],pTPBin[1]); //pT det
   TH1D* hDMpTP=(TH1D*)fhnDeltaMass->Projection(axisP);
   hDMpTP->SetName(hprjname);
   hDMpTP->SetLineWidth(2);
   return hDMpTP;
}

 void DefNewBins(TH1F* hTrig, Double_t* rebinpt, Int_t newnptbins){
   const Int_t nsteps=4;
   //Double_t steps[nsteps]={28,35,50,70}; //(in GeV/c)
   //Int_t deltapt[nsteps]={2,6,15,45};//rebinning at each step 
   Double_t steps[nsteps]={20,50,70,90}; //(in GeV/c)
   Int_t deltapt[nsteps]={2,6,15,45};//rebinning at each step
   
      Int_t ii=0,jj=0;
      for(Int_t k=0;k<newnptbins;k++){
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

Int_t GetJetPtBinTemp(Double_t ptj, Double_t delta, Double_t min){
   Int_t bin = TMath::FloorNint((ptj - min)/delta);
   return bin;
}

/*
TH1* CopyAndRebin(TH1* h, TString newname, Int_t nbinnew, Double_t binrangesnew[]){
   
   TH1* hrebin = h->Rebin(nbinnew, Form("%srb", h->GetName()), binrangesnew);


} */
