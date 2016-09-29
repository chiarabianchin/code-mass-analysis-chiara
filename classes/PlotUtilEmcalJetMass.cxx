#include "PlotUtilEmcalJetMass.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "Riostream.h"
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include <TRandom3.h>
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TList.h"

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

using std::cout;
using std::endl;

ClassImp(PlotUtilEmcalJetMass)

//----------------------------------------------------------------------------
PlotUtilEmcalJetMass::PlotUtilEmcalJetMass():
  fFileInName(0),
  fFileIn(0),
  fListArr(0),
  fRhoListArr(0),
  fRadius(0.4),
  fPtMin(60.),
  fPtMax(80.),
  fPtPMin(0.),
  fPtPMax(99990.),
  fTaskName(""),
  fJetName(""),
  fJetType("Charged"),
  fTag(""),
  fConstTag(""),
  fCentBin(0),
  fMinLeadTrkPt(-1.),
  fMaxLeadTrkPt(1e6),
  fScaleFactor(1),
  fListNameType(0),
  fStrTrk("PicoTracks_pT0150"),
  fStrClus("CaloClustersCorr"),
  fStrJetType("Charged")
{
  fListArr.SetOwner(kTRUE);
  fRhoListArr.SetOwner(kTRUE);
}

//----------------------------------------------------------------------------
Bool_t PlotUtilEmcalJetMass::LoadFile() {

  if(!fFileInName.IsNull())
    fFileIn = new TFile(fFileInName.Data());
  else {
    Printf("No file name provided: %s",fFileInName.Data());
    return kFALSE;
  }
  fFileIn->ls();
  // if(!LoadList()) return kFALSE;

  return kTRUE;
}

//----------------------------------------------------------------------------
Bool_t PlotUtilEmcalJetMass::LoadList() {

  if(!fFileIn) {
    Printf("Cannot load list without file.");
    return kFALSE;
  }

  //TString strTrk = "PicoTracks_pT0150";
  //if(fJetType.Contains("MCParticlesSelected")) strTrk = "MCParticlesSelected_pT0000";
  //if(fJetType.Contains("ThrmTracksEmb")) strTrk = "ThrmTracksEmb_pT0150";
  //TString strClus = "";
  //if(fJetType.Contains("Full")) strClus = "CaloClustersCorr";
  //TString strJetType = "Charged";
  //if(fJetType.Contains("Full")) strJetType = "Full";
  if(fJetType.Contains("Full")) fStrClus = "CaloClustersCorr";
  else fStrClus = "";
  // exemple: JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub_TC
  TString listName = Form("JetMass%s_Jet%s_AKT%sR%03d_%s_%sE_scheme%s%s",fTaskName.Data(), fJetName.Data(), fStrJetType.Data(), (Int_t)(fRadius*100),fStrTrk.Data(), fStrClus.Data(),fConstTag.Data(), fTag.Data());
  Printf("listName: %s",listName.Data());
  TList *lst = dynamic_cast<TList*>(fFileIn->Get(listName.Data()));
  if(!lst) {
    Printf("Couldn't find list %s",listName.Data());
    //fListArr.Add(0x0);
    return kFALSE;
  }
  fListArr.Add(lst);
  return kTRUE;
}
//----------------------------------------------------------------------------
Bool_t PlotUtilEmcalJetMass::LoadRhoLists(TString trkname) {

  if(!fFileIn) {
    Printf("Cannot load list without file.");
    return kFALSE;
  }
 Printf("Jet Type = %s", fJetType.Data());
  TString strTrk = Form("%s_pT0150", trkname.Data());
  if(fJetType.Contains("MCParticlesSelected")) strTrk = "MCParticlesSelected_pT0000";
  TString strClus = "";
  if(fJetType.Contains("Full")) strClus = "CaloClustersCorr";
  TString strJetType = "Charged";
  if(fJetType.Contains("Full")) strJetType = "Full";
  TString sparse = "Sparse", ending = "TPC_histos";
  //rho
  TString listName;
  if(!fListNameType) listName = Form("Rho%s%sR%03d_Jet%s_KT%sR%03d_%s_%sE_scheme_%s%s",sparse.Data() , fTaskName.Data(), (Int_t)(fRadius*100), fJetName.Data(), strJetType.Data(), (Int_t)(fRadius*100),strTrk.Data(), strClus.Data(),ending.Data(), fTag.Data());
  else listName = Form("Rho%s%s_Jet%s_KT%sR%03d_%s_%sE_scheme_%s%s",sparse.Data() , fTaskName.Data(), fJetName.Data(), strJetType.Data(), (Int_t)(fRadius*100),strTrk.Data(), strClus.Data(),ending.Data(), fTag.Data()); //figure out how to do handle the different names
  Printf("listName: %s",listName.Data());
  TList *lst = dynamic_cast<TList*>(fFileIn->Get(listName.Data()));
  if(!lst) {
    Printf("Couldn't find list %s",listName.Data());
    //fListArr.Add(0x0);
    return kFALSE;
  }
  fRhoListArr.Add(lst);
  
  //rho mass
  if(!fListNameType) listName = Form("RhoMass%s%sR%03d_Jet%s_KT%sR%03d_%s_%sE_scheme_%s%s", sparse.Data(), fTaskName.Data(), (Int_t)(fRadius*100), fJetName.Data(), strJetType.Data(), (Int_t)(fRadius*100),strTrk.Data(), strClus.Data(),ending.Data(), fTag.Data());
  else listName = Form("RhoMass%s%s_Jet%s_KT%sR%03d_%s_%sE_scheme_%s%s", sparse.Data(), fTaskName.Data(), fJetName.Data(), strJetType.Data(), (Int_t)(fRadius*100),strTrk.Data(), strClus.Data(),ending.Data(), fTag.Data());
  Printf("listName: %s",listName.Data());
  lst = dynamic_cast<TList*>(fFileIn->Get(listName.Data()));
  if(!lst) {
    Printf("Couldn't find list %s",listName.Data());
    //fListArr.Add(0x0);
    return kFALSE;
  }
  fRhoListArr.Add(lst);
  return kTRUE;
}

//----------------------------------------------------------------------------
TH3* PlotUtilEmcalJetMass::Scale(TH3* horig){
   // time consuming, but at the moment no better solution found
   
   TH3* hnew = (TH3D*)horig->Clone(Form("%sbis", horig->GetName()));
   Printf("Scale factor %f", fScaleFactor);
   hnew->Scale(1./fScaleFactor);
   Printf("Integral %.3f -> %.3f", horig->Integral(), hnew->Integral());
   return hnew;

}
//----------------------------------------------------------------------------
void PlotUtilEmcalJetMass::ClearScaledH(TH3* hnew){
   // time consuming, but at the moment no better solution found
   
   delete hnew;
   
   return;
}
//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetMassDistribution(JetSelection jSel, Int_t iList) {

  TString strJSel = "";
  if(jSel==kJetAll) strJSel = "AllSel";
  else if(jSel==kJetTagged) strJSel = "Tagged";
  else if(jSel==kJetTaggedMatch) strJSel = "TaggedMatch";
  else {
    Printf("Jet selection unknown. Returning null pointer.");
    return 0x0;
  }

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fh3PtJet1VsMassVsLeadPt%s_%d",strJSel.Data(),fCentBin);
  TH3 *fh3PtVsMassVsLeadPtJet1 = dynamic_cast<TH3*>(lst->FindObject(histName.Data()));
  if(!fh3PtVsMassVsLeadPtJet1) {
    Printf("Couldn't find fh3PtVsMassVsLeadPtJet1 %s. Returning null pointer.",histName.Data());
    return 0x0;
  }
  //if(!fh3PtVsMassVsLeadPtJet1->GetSumw2()) 
  fh3PtVsMassVsLeadPtJet1->Sumw2();
  TH3 *fh3PtVsMassVsLeadPtJet1bis = Scale(fh3PtVsMassVsLeadPtJet1);
  //if(!fh3PtVsMassVsLeadPtJet1bis) fh3PtVsMassVsLeadPtJet1bis=fh3PtVsMassVsLeadPtJet1;
  Int_t minz = fh3PtVsMassVsLeadPtJet1bis->GetZaxis()->FindBin(fMinLeadTrkPt+0.000001);
  Int_t maxz = fh3PtVsMassVsLeadPtJet1bis->GetZaxis()->FindBin(fMaxLeadTrkPt-0.000001);
  Int_t minx = fh3PtVsMassVsLeadPtJet1bis->GetXaxis()->FindBin(fPtMin+0.00001);
  Int_t maxx = fh3PtVsMassVsLeadPtJet1bis->GetXaxis()->FindBin(fPtMax-0.00001);
  Printf("Range for projection in pT jet %d - %d", minx, maxx);
  TH1D *hM = dynamic_cast<TH1D*>(fh3PtVsMassVsLeadPtJet1bis->ProjectionY(Form("fh1Mass%s_Lst%d_Cent%d_LeadTrk%.0f", strJSel.Data(),iList,fCentBin,fMinLeadTrkPt),minx,maxx,minz,maxz));
  ClearScaledH(fh3PtVsMassVsLeadPtJet1bis);
  return hM;
}

//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetPtDistribution(JetSelection jSel, Int_t iList) {

  TString strJSel = "";
  if(jSel==kJetAll) strJSel = "AllSel";
  else if(jSel==kJetTagged) strJSel = "Tagged";
  else if(jSel==kJetTaggedMatch) strJSel = "TaggedMatch";
  else {
     Printf("Jet selection %d unknown. Returning null pointer.", jSel);
     return 0x0;
  }

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fh3PtJet1VsMassVsLeadPt%s_%d",strJSel.Data(),fCentBin);
  TH3 *fh3PtVsMassVsLeadPtJet1 = dynamic_cast<TH3*>(lst->FindObject(histName.Data()));
  if(!fh3PtVsMassVsLeadPtJet1) {
    Printf("Couldn't find fh3PtVsMassVsLeadPtJet1 %s. Returning null pointer.",histName.Data());
    return 0x0;
  }

  Int_t minz = fh3PtVsMassVsLeadPtJet1->GetZaxis()->FindBin(fMinLeadTrkPt+0.000001);
  Int_t maxz = fh3PtVsMassVsLeadPtJet1->GetZaxis()->FindBin(fMaxLeadTrkPt-0.000001);
  
  TH1D *hPt = dynamic_cast<TH1D*>(fh3PtVsMassVsLeadPtJet1->ProjectionX(Form("fh1Pt%s_Lst%d_Cent%d_LeadTrk%.0f_%.0f",strJSel.Data(),iList,fCentBin,fMinLeadTrkPt,fMaxLeadTrkPt),1,-1,minz,maxz));
  hPt->Sumw2();
  return hPt;
}

//----------------------------------------------------------------------------
TH1* PlotUtilEmcalJetMass::GetCentralityHist(Int_t iList) const {

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fHistCentrality");
  TH1 *fHistCentrality = dynamic_cast<TH1*>(lst->FindObject(histName.Data()));
  return fHistCentrality;
}

//----------------------------------------------------------------------------
Double_t PlotUtilEmcalJetMass::GetNEvents(Int_t iList) const {

  // TList *lst = static_cast<TList*>(fListArr.At(iList));
  // if(!lst){
  //   Printf("Couldn't find list %d",iList);
  //   return 0x0;
  // }

  // TString histName = Form("fHistCentrality");
  // TH1 *fHistCentrality = dynamic_cast<TH1*>(lst->FindObject(histName.Data()));
  TH1 *fHistCentrality = GetCentralityHist(iList);
  Double_t centMin[4] = {0.,10.,30.,50.};
  Double_t centMax[4] = {10.,30.,50.,100.};
  Int_t min = fHistCentrality->GetXaxis()->FindBin(centMin[fCentBin]+0.000001);
  Int_t max = fHistCentrality->GetXaxis()->FindBin(centMax[fCentBin]-0.000001);
  Double_t nEvt = fHistCentrality->Integral(min,max);
  return nEvt;
}

//----------------------------------------------------------------------------
Double_t PlotUtilEmcalJetMass::GetNEventsAll(Int_t iList) const {

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fHistCentrality");
  TH1 *fHistCentrality = dynamic_cast<TH1*>(lst->FindObject(histName.Data()));
  Double_t nEvt = fHistCentrality->Integral();
  return nEvt;
}

//----------------------------------------------------------------------------
TH2D* PlotUtilEmcalJetMass::GetJetMassVsPt(JetSelection jSel, Int_t iList) {

  TString strJSel = "";
  if(jSel==kJetAll) strJSel = "AllSel";
  else if(jSel==kJetTagged) strJSel = "Tagged";
  else if(jSel==kJetTaggedMatch) strJSel = "TaggedMatch";
  else {
    Printf("Jet selection unknown. Returning null pointer.");
    return 0x0;
  }

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fh3PtJet1VsMassVsLeadPt%s_%d",strJSel.Data(),fCentBin);
  TH3 *fh3PtVsMassVsLeadPtJet1 = dynamic_cast<TH3*>(lst->FindObject(histName.Data()));
  if(!fh3PtVsMassVsLeadPtJet1) {
    Printf("Couldn't find fh3PtVsMassVsLeadPtJet1 %s. Returning null pointer.",histName.Data());
    return 0x0;
  }

  Int_t minz = fh3PtVsMassVsLeadPtJet1->GetZaxis()->FindBin(fMinLeadTrkPt+0.000001);
  Int_t maxz = fh3PtVsMassVsLeadPtJet1->GetZaxis()->FindBin(fMaxLeadTrkPt-0.000001);
  fh3PtVsMassVsLeadPtJet1->GetZaxis()->SetRange(minz,maxz);
  TH2D *hMPt = dynamic_cast<TH2D*>(fh3PtVsMassVsLeadPtJet1->Project3D("yx"));
  hMPt->SetName(Form("fh2MPt%s_Lst%d_Cent%d_LeadTrk%.0f_%.0f",strJSel.Data(),iList,fCentBin,fMinLeadTrkPt,fMaxLeadTrkPt));
  return hMPt;
}

//----------------------------------------------------------------------------
TGraphErrors* PlotUtilEmcalJetMass::GetMeanJetMassVsPt(JetSelection jSel, Int_t iList) {

  const Int_t nPtBins = 5;
  Double_t ptmin[nPtBins] = {20.,40.,60.,80.,100.};
  Double_t ptmax[nPtBins] = {40.,60.,80.,100.,120.};
  
  TGraphErrors *gr = new TGraphErrors();
  for(Int_t i = 0; i<nPtBins; i++) {
    SetJetPtRange(ptmin[i],ptmax[i]);
    TH1D *hM = GetJetMassDistribution(jSel,iList);
    gr->SetPoint(gr->GetN(),ptmin[i] + 0.5*(ptmax[i]-ptmin[i]), hM->GetMean());
    gr->SetPointError(gr->GetN()-1,0.5*(ptmax[i]-ptmin[i]), hM->GetMeanError());
  }
  return gr;
}

//----------------------------------------------------------------------------
TGraphErrors* PlotUtilEmcalJetMass::GetRMSJetMassVsPt(JetSelection jSel, Int_t iList) {

  const Int_t nPtBins = 5;
  Double_t ptmin[nPtBins] = {20.,40.,60.,80.,100.};
  Double_t ptmax[nPtBins] = {40.,60.,80.,100.,120.};
  
  TGraphErrors *gr = new TGraphErrors();
  for(Int_t i = 0; i<nPtBins; i++) {
    SetJetPtRange(ptmin[i],ptmax[i]);
    TH1D *hM = GetJetMassDistribution(jSel,iList);
    gr->SetPoint(gr->GetN(),ptmin[i] + 0.5*(ptmax[i]-ptmin[i]), hM->GetRMS());
    gr->SetPointError(gr->GetN()-1,0.5*(ptmax[i]-ptmin[i]), hM->GetRMSError());
  }
  return gr;
}

//----------------------------------------------------------------------------

TH1D* PlotUtilEmcalJetMass::GetMassInCentrality(Double_t c1, Double_t c2, JetSelection jSel, Int_t iList){
   
   TString strJSel = "";
   if(jSel==kJetAll) strJSel = "AllSel";
   else if(jSel==kJetTagged) strJSel = "Tagged";
   else if(jSel==kJetTaggedMatch) strJSel = "TaggedMatch";
   else {
      Printf("Jet selection %d unknown. Returning null pointer.", jSel);
      return 0x0;
   }
   TList *lst = static_cast<TList*>(fListArr.At(iList));
   if(!lst){
      Printf("Couldn't find list %d",iList);
      return 0x0;
   }
   TString histName = Form("fh3PtJet1VsMassVsCent%s",strJSel.Data());
   TH3 *fh3PtVsMassVsCent = dynamic_cast<TH3*>(lst->FindObject(histName.Data()));
   if(!fh3PtVsMassVsCent) {
    Printf("Couldn't find fh3PtVsMassVsCent %s. Returning null pointer.",histName.Data());
    return 0x0;
  }
   Int_t minx = fh3PtVsMassVsCent->GetXaxis()->FindBin(fPtMin);
   Int_t maxx = fh3PtVsMassVsCent->GetXaxis()->FindBin(fPtMax)-1;
   Int_t minz = fh3PtVsMassVsCent->GetZaxis()->FindBin(c1);
   Int_t maxz = fh3PtVsMassVsCent->GetZaxis()->FindBin(c2)-1;
   
   TH1D *hM = dynamic_cast<TH1D*>(fh3PtVsMassVsCent->ProjectionY(Form("fh1Mass%sLst%d_Cent%.0f-%.0f_pTj%.0f-%.0f", strJSel.Data(),iList,c1, c2,fPtMin,fPtMax),minx,maxx,minz,maxz));
   return hM;

}

TH2D* PlotUtilEmcalJetMass::GetRho(Int_t sel){
   //sel = 1 -> centrality (fHistRhovsCent)
   //sel = 2 -> Ntrack     (fHistRhovsNtrackvsV0Mult)
   //sel = 3 -> V0Mult     (fHistRhovsNtrackvsV0Mult)
   TList *lst = static_cast<TList*>(fRhoListArr.At(0)); //rho list
   if(!lst){
      Printf("Couldn't find list %d",0);
      return 0x0;
   }
   
   TString hname2D = "fHistRhovsCent", hname3D = "fHistRhovsNtrackvsV0Mult";
   TH2D *hRho = 0x0; 
   if(sel==1) {
      TH2 *fHistRhovsCent = dynamic_cast<TH2*>(lst->FindObject(hname2D.Data()));

      if(!fHistRhovsCent) {
      	 Printf("Couldn't find fHistRhovsCent %s. Returning null pointer.",hname2D.Data());
      	 lst->ls();
      	 return 0x0;
      }
      hRho = (TH2D*)fHistRhovsCent->Clone("hRho_Cent");
      
   }
   
   if(sel>1){
      TH3 *fHistRhovsNtrackvsV0Mult = dynamic_cast<TH3*>(lst->FindObject(hname3D.Data()));
      if(!fHistRhovsNtrackvsV0Mult) {
      	 Printf("Couldn't find fHistRhovsNtrackvsV0Mult %s. Returning null pointer.",hname3D.Data());
      	 return 0x0;
      }
      if(sel == 3){
      	 hRho = dynamic_cast<TH2D*>(fHistRhovsNtrackvsV0Mult->Project3D("zx"));
      	 hRho->SetName(Form("hRho_V0M"));
      	 
      }
      if(sel == 2){
      	 hRho = dynamic_cast<TH2D*>(fHistRhovsNtrackvsV0Mult->Project3D("yx"));
      	 hRho->SetName(Form("hRho_NTrack"));
      	 
      }
   }
  
   return hRho;
   
}

//----------------------------------------------------------------------------
TH2D* PlotUtilEmcalJetMass::GetRhoM(Int_t sel){
   //sel = 1 -> centrality (fHistRhoMassvsCent)
   //sel = 2 -> Ntrack     (fHistRhoMassvsNtrack)
   //
   TList *lst = static_cast<TList*>(fRhoListArr.At(1)); //rho list
   if(!lst){
    Printf("Couldn't find list %d",0);
    return 0x0;
  }
  
  TString hname2D = "fHistRhoMassvsCent", hname2D2 = "fHistRhoMassvsNtrack";
  TH2D *hRhoMass = 0x0; 
  if(sel<2) {
     TH2 *fHistRhoMassvsCent = dynamic_cast<TH2*>(lst->FindObject(hname2D.Data()));
     if(!fHistRhoMassvsCent) {
     	lst->ls();
     	Printf("Couldn't find fHistRhoMassvsCent %s. Returning null pointer.",hname2D.Data());
     	return 0x0;
     }
     hRhoMass = (TH2D*)fHistRhoMassvsCent->Clone("hRhom_Cent");
     
     
  }

  if(sel>1){
     TH2 *fHistRhoMassvsNtrack = dynamic_cast<TH2*>(lst->FindObject(hname2D2.Data()));
     if(!fHistRhoMassvsNtrack) {
     	Printf("Couldn't find fHistRhoMassvsNtrack %s. Returning null pointer.",hname2D.Data());
     	return 0x0;
     }

     hRhoMass = dynamic_cast<TH2D*>(fHistRhoMassvsNtrack->Clone("hRhoMass_Ntrack"));
  }
  return hRhoMass;
}

TH1D* PlotUtilEmcalJetMass::GetRhoProjection(Double_t c1, Double_t c2, Int_t sel){
   //sel = 0 -> no sel     (fHistRhovsCent)
   //sel = 1 -> centrality (fHistRhovsCent)
   //sel = 2 -> Ntrack     (fHistRhovsNtrackvsV0Mult)
   //sel = 3 -> V0Mult     (fHistRhovsNtrackvsV0Mult)
   TList *lst = static_cast<TList*>(fRhoListArr.At(0)); //rho list
   if(!lst){
    Printf("Couldn't find list %d",0);
    return 0x0;
  }
  
  TString hname2D = "fHistRhovsCent", hname3D = "fHistRhovsNtrackvsV0Mult";
  TH1D *hRho = 0x0; 
  if(sel<2) {
     TH2 *fHistRhovsCent = dynamic_cast<TH2*>(lst->FindObject(hname2D.Data()));
     if(!fHistRhovsCent) {
     	Printf("Couldn't find fHistRhovsCent %s. Returning null pointer.",hname2D.Data());
     	return 0x0;
     }
     
     //c1 and c2 are centrality percentiles
     Int_t minx = fHistRhovsCent->GetXaxis()->FindBin(c1);
     Int_t maxx = fHistRhovsCent->GetXaxis()->FindBin(c2)-1;
     if(sel == 0) {
     	minx = 0; maxx=-1;
     }
     hRho = dynamic_cast<TH1D*>(fHistRhovsCent->ProjectionY(Form("hRho_Cent%.0f-%.0f", c1, c2),minx,maxx));
     
     
  }

  if(sel>1){
     TH3 *fHistRhovsNtrackvsV0Mult = dynamic_cast<TH3*>(lst->FindObject(hname3D.Data()));
     if(!fHistRhovsNtrackvsV0Mult) {
     	Printf("Couldn't find fHistRhovsNtrackvsV0Mult %s. Returning null pointer.",hname3D.Data());
     	return 0x0;
     }
     if(sel == 3){
     	//c1 and c2 are V0 multiplicity
     	Int_t minz = fHistRhovsNtrackvsV0Mult->GetZaxis()->FindBin(c1);
     	Int_t maxz = fHistRhovsNtrackvsV0Mult->GetZaxis()->FindBin(c2)-1;
     	hRho = dynamic_cast<TH1D*>(fHistRhovsNtrackvsV0Mult->ProjectionZ(Form("hRho_V0M%.0f-%.0f", c1, c2),minz,maxz));
     	
     }
     if(sel == 2){
     	//c1 and c2 are number of tracks
     	Int_t minx = fHistRhovsNtrackvsV0Mult->GetXaxis()->FindBin(c1);
     	Int_t maxx = fHistRhovsNtrackvsV0Mult->GetXaxis()->FindBin(c2)-1;
     	hRho = dynamic_cast<TH1D*>(fHistRhovsNtrackvsV0Mult->ProjectionY(Form("hRho_NTrack%.0f-%.0f", c1, c2),minx,maxx));
     	
     }
  }
  hRho->Sumw2();
  return hRho;
}
//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetRhoMProjection(Double_t c1, Double_t c2, Int_t sel){
   //sel = 0 -> no sel     (fHistRhoMassvsCent)
   //sel = 1 -> centrality (fHistRhoMassvsCent)
   //sel = 2 -> Ntrack     (fHistRhoMassvsNtrack)
   //
   TList *lst = static_cast<TList*>(fRhoListArr.At(1)); //rho list
   if(!lst){
    Printf("Couldn't find list %d",0);
    return 0x0;
  }
  
  TString hname2D = "fHistRhoMassvsCent", hname2D2 = "fHistRhoMassvsNtrack";
  TH1D *hRhoMass = 0x0; 
  if(sel<2) {
     TH2 *fHistRhoMassvsCent = dynamic_cast<TH2*>(lst->FindObject(hname2D.Data()));
     if(!fHistRhoMassvsCent) {
     	lst->ls();
     	Printf("Couldn't find fHistRhoMassvsCent %s. Returning null pointer.",hname2D.Data());
     	return 0x0;
     }
     
     //c1 and c2 are centrality percentiles
     Int_t minx = fHistRhoMassvsCent->GetXaxis()->FindBin(c1);
     Int_t maxx = fHistRhoMassvsCent->GetXaxis()->FindBin(c2)-1;
     if(sel == 0) {
     	minx = 0; maxx=-1;
     }
     hRhoMass = dynamic_cast<TH1D*>(fHistRhoMassvsCent->ProjectionY(Form("hRhoMass_Cent%.0f-%.0f", c1, c2),minx,maxx));
     
     
  }

  if(sel>1){
     TH2 *fHistRhoMassvsNtrack = dynamic_cast<TH2*>(lst->FindObject(hname2D2.Data()));
     if(!fHistRhoMassvsNtrack) {
     	Printf("Couldn't find fHistRhoMassvsNtrack %s. Returning null pointer.",hname2D.Data());
     	return 0x0;
     }
     if(sel == 3){
     	//c1 and c2 are number of tracks
     	Int_t minx = fHistRhoMassvsNtrack->GetXaxis()->FindBin(c1);
     	Int_t maxx = fHistRhoMassvsNtrack->GetXaxis()->FindBin(c2)-1;
     	hRhoMass = dynamic_cast<TH1D*>(fHistRhoMassvsNtrack->ProjectionY(Form("hRhoMass_V0M%.0f-%.0f", c1, c2),minx,maxx));
     	
     }
  }
  hRhoMass->Sumw2();
  return hRhoMass;
}
//----------------------------------------------------------------------------

TH3* PlotUtilEmcalJetMass::GetEPCorr3D(JetSelection jSel, Int_t iList) {

  TString strJSel = "";
  if(jSel==kJetAll) strJSel = "AllSel";
  else if(jSel==kJetTagged) strJSel = "Tagged";
  else if(jSel==kJetTaggedMatch) strJSel = "TaggedMatch";
  else {
    Printf("Jet selection unknown. Returning null pointer.");
    return 0x0;
  }

  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fh3JetPtVsMassVsEPRel%s_%d",strJSel.Data(),fCentBin);
  TH3 *fh3JetPtVsMassVsEPRel = dynamic_cast<TH3*>(lst->FindObject(histName.Data()));
  if(!fh3JetPtVsMassVsEPRel) {
    Printf("Couldn't find fh3JetPtVsMassVsEPRel %s. Returning null pointer.",histName.Data());
    return 0x0;
  }
  return fh3JetPtVsMassVsEPRel;
}

//----------------------------------------------------------------------------
TH2D* PlotUtilEmcalJetMass::GetJetMassVsEP(JetSelection jSel, Int_t iList) {

  TH3 *fh3JetPtVsMassVsEPRel = GetEPCorr3D(jSel,iList);
  if(!fh3JetPtVsMassVsEPRel) return 0x0;

  Int_t minx = fh3JetPtVsMassVsEPRel->GetXaxis()->FindBin(fPtMin+0.00001);
  Int_t maxx = fh3JetPtVsMassVsEPRel->GetXaxis()->FindBin(fPtMax-0.00001);
  fh3JetPtVsMassVsEPRel->GetXaxis()->SetRange(minx,maxx);
  TH2D *h2 = dynamic_cast<TH2D*>(fh3JetPtVsMassVsEPRel->Project3D("yz"));
  return h2;
}

//----------------------------------------------------------------------------
TH2D* PlotUtilEmcalJetMass::GetJetPtVsEP(JetSelection jSel, Int_t iList) {

  TH3 *fh3JetPtVsMassVsEPRel = GetEPCorr3D(jSel,iList);
  if(!fh3JetPtVsMassVsEPRel) return 0x0;

  Int_t minx = fh3JetPtVsMassVsEPRel->GetXaxis()->FindBin(fPtMin+0.00001);
  Int_t maxx = fh3JetPtVsMassVsEPRel->GetXaxis()->FindBin(fPtMax-0.00001);
  fh3JetPtVsMassVsEPRel->GetXaxis()->SetRange(minx,maxx);
  TH2D *h2 = dynamic_cast<TH2D*>(fh3JetPtVsMassVsEPRel->Project3D("xz"));
  return h2;
}

//----------------------------------------------------------------------------
TProfile* PlotUtilEmcalJetMass::GetJetMassVsEPProf(JetSelection jSel, Int_t iList) {

  TH2D *h2 = GetJetMassVsEP(jSel,iList);
  if(!h2) return 0x0;

  TProfile *p = dynamic_cast<TProfile*>(h2->ProfileX(Form("hpM%d",iList),1,-1));//,"s"));
  return p;
}

//----------------------------------------------------------------------------
TProfile* PlotUtilEmcalJetMass::GetJetPtVsEPProf(JetSelection jSel, Int_t iList) {

  TH2D *h2 = GetJetPtVsEP(jSel,iList);
  if(!h2) return 0x0;

  TProfile *p = dynamic_cast<TProfile*>(h2->ProfileX(Form("hpPt%d",iList),1,-1));//,"s"));
  return p;
}

//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetMassDetFromTHnSparse(JetSelection jSel, Int_t iList){
  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fhnMassResponse");


  THnSparse* hns = (THnSparse*)lst->FindObject(histName);
  Int_t minx = hns->GetAxis(2)->FindBin(fPtMin+0.00001);
  Int_t maxx = hns->GetAxis(2)->FindBin(fPtMax-0.00001);
  hns->GetAxis(2)->SetRange(minx, maxx);
  TH1D* hMassDet = (TH1D*)hns->Projection(0);
  Printf("Mean value is %.3f", hMassDet->GetMean());
  return hMassDet;
}

//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetpTDetFromTHnSparse(JetSelection jSel, Int_t iList){
  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fhnMassResponse");


  THnSparse* hns = (THnSparse*)lst->FindObject(histName);
  Int_t minx = hns->GetAxis(3)->FindBin(fPtPMin+0.00001);
  Int_t maxx = hns->GetAxis(3)->FindBin(fPtPMax-0.00001);
  hns->GetAxis(3)->SetRange(minx, maxx);
  TH1D* hMassDet = (TH1D*)hns->Projection(2);
  
  return hMassDet;
}

//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetMassPartFromTHnSparse(JetSelection jSel, Int_t iList){
  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fhnMassResponse");


  THnSparseD* hns = (THnSparseD*)lst->FindObject(histName);
  if(!hns) {
     Printf("Not found");
     lst->ls();
  }
  Int_t minx = hns->GetAxis(2)->FindBin(fPtMin+0.00001);
  Int_t maxx = hns->GetAxis(2)->FindBin(fPtMax-0.00001);
  Printf("Select %f (%d) - %f (%d)", fPtMin, minx, fPtMax, maxx);
  hns->GetAxis(2)->SetRange(minx, maxx);
  TH1D* hMassPart = dynamic_cast<TH1D*>(hns->Projection(1));
  hMassPart->SetName(Form("hMassPart%d",iList));
  Printf("%p %s has %.0f entries",  hMassPart , hMassPart->GetName(), hMassPart->GetEntries());
  TCanvas *c=new TCanvas("c","c");
  c->cd(); hMassPart->Draw();
  Printf("Mean value is %.3f", hMassPart->GetMean());
  return hMassPart;
}
//----------------------------------------------------------------------------
TH1D* PlotUtilEmcalJetMass::GetJetpTPartFromTHnSparse(JetSelection jSel, Int_t iList){
  TList *lst = static_cast<TList*>(fListArr.At(iList));
  if(!lst){
    Printf("Couldn't find list %d",iList);
    return 0x0;
  }

  TString histName = Form("fhnMassResponse");


  THnSparse* hns = (THnSparse*)lst->FindObject(histName);
  Int_t minx = hns->GetAxis(2)->FindBin(fPtMin+0.00001);
  Int_t maxx = hns->GetAxis(2)->FindBin(fPtMax-0.00001);
  hns->GetAxis(2)->SetRange(minx, maxx);
  TH1D* hMassDet = (TH1D*)hns->Projection(3);
  return hMassDet;
}
