#ifndef __PlotUtilEmcalJetMass_h__
#define __PlotUtilEmcalJetMass_h__
#include "Rtypes.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TDirectoryFile.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TColor.h"

class TGraph;
class TGraphErrors;
class TArrayF;
class TArrayD;

//
// post processing AliAnalysisTaskEmcalJetMass
//

class PlotUtilEmcalJetMass {
 public:
  PlotUtilEmcalJetMass();
  //PlotUtilEmcalJetMass(const PlotUtilEmcalJetMass& obj); // copy constructor
  //PlotUtilEmcalJetMass& operator=(const PlotUtilEmcalJetMass& other); // assignment
  virtual ~PlotUtilEmcalJetMass() {;}

  enum JetSelection {
    kJetAll         = 0,
    kJetTagged      = 1,
    kJetTaggedMatch = 2
  };

  void SetInputFileName(TString str = "AnalysisResults.root") { fFileInName = str; }
  void SetJetRadius(Double_t r)                               { fRadius     = r;   }
  void SetTaskName(TString name= "Structure")                 { fTaskName    = name;}
  void SetJetName(TString name= "Rhosub")                     { fJetName    = name;}
  void SetJetType(TString type = "Charged")                   { fJetType    = type;}
  void SetTag(TString tag)                                    { fTag        = tag; }
  void SetConstTag(TString tag = "ConstSub")                  { fConstTag   = tag; }
  void SetCentBin(Int_t cb = 0)                               { fCentBin    = cb;  }

  void SetJetPtRange(Double_t min, Double_t max)              { fPtMin = min; fPtMax = max; }
  void SetJetPtPRange(Double_t min, Double_t max)             { fPtPMin = min; fPtPMax = max; }
  void SetMinLeadTrackPt(Double_t min)                        { fMinLeadTrkPt = min; }
  void SetMaxLeadTrackPt(Double_t max)                        { fMaxLeadTrkPt = max; }
  void SetScaleFactor(Double_t factor)                        { fScaleFactor = factor; }
  void SetListNameType(Int_t type)                            { fListNameType = type; }
  
  void SetTrackString(TString name)                           { fStrTrk = name; }
  void SetClusterString(TString name)                         { fStrClus = name; }
  void SetJTypeString(TString name)                           { fStrJetType = name; }
  
  Bool_t LoadFile();
  Bool_t LoadList();

  TFile* GetFileIn() const {return fFileIn;}
  TH1D *GetJetMassDistribution(JetSelection jSel, Int_t iList);
  TH1D *GetJetPtDistribution(JetSelection jSel, Int_t iList);
  
  TH1D* GetMassInCentrality(Double_t c1, Double_t c2, JetSelection jsel, Int_t iList);
  
  TH2D *GetJetMassVsEP(JetSelection jSel, Int_t iList);
  TH2D *GetJetPtVsEP(JetSelection jSel, Int_t iList);
  TH2D* GetJetMassVsPt(JetSelection jSel, Int_t iList);
  
  TGraphErrors* GetMeanJetMassVsPt(JetSelection jSel, Int_t iList);
  TGraphErrors* GetRMSJetMassVsPt(JetSelection jSel, Int_t iList);
  
  TProfile *GetJetMassVsEPProf(JetSelection jSel, Int_t iList);
  TProfile *GetJetPtVsEPProf(JetSelection jSel, Int_t iList);

  TH3 *GetEPCorr3D(JetSelection jSel, Int_t iList);

  TH1D* GetJetMassDetFromTHnSparse(JetSelection jSel, Int_t iList);
  TH1D* GetJetMassPartFromTHnSparse(JetSelection jSel, Int_t iList);
  TH1D* GetJetpTDetFromTHnSparse(JetSelection jSel, Int_t iList);
  TH1D* GetJetpTPartFromTHnSparse(JetSelection jSel, Int_t iList);
  Double_t GetScaleFactor() const    { return fScaleFactor; }
  
  Double_t GetNEvents(Int_t iList) const;
  Double_t GetNEventsAll(Int_t) const;
  TH1* GetCentralityHist(Int_t iList) const;
  
  //Rho analysis
  Bool_t LoadRhoLists(TString trkname = "PicoTracks");
  TH2D* GetRho(Int_t sel=1);
  TH2D* GetRhoM(Int_t sel=1);
  TH1D* GetRhoProjection(Double_t c1, Double_t c2, Int_t sel);
  TH1D* GetRhoMProjection(Double_t c1, Double_t c2, Int_t sel);
  
 private:
    TH3* Scale(TH3* horig);
    void  ClearScaledH(TH3* hnew);
 protected:
  TString          fFileInName;
  TFile           *fFileIn; //input file
  TObjArray        fListArr;
  TObjArray        fRhoListArr;
  Double_t         fRadius;
  Double_t         fPtMin;
  Double_t         fPtMax;
  Double_t         fPtPMin;
  Double_t         fPtPMax;
  TString          fTaskName;
  TString          fJetName;
  TString          fJetType;
  TString          fTag;
  TString          fConstTag;
  Int_t            fCentBin;
  Double_t         fMinLeadTrkPt;
  Double_t         fMaxLeadTrkPt;
  Double_t         fScaleFactor;
  Int_t            fListNameType;
  TString          fStrTrk;
  TString          fStrClus;
  TString          fStrJetType;
};
#endif
