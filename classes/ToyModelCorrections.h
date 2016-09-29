#ifndef __ToyModelCorrections_h__
#define __ToyModelCorrections_h__
#include "Rtypes.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TDirectoryFile.h"
#include "TRandom3.h"

class TGraph;
class TGraphErrors;
class TArrayF;
class TArrayD;

//
// toy model dijet response
//

class ToyModelCorrections {
 public:
  ToyModelCorrections();
  //ToyModelCorrections(const ToyModelCorrections& obj); // copy constructor
  //ToyModelCorrections& operator=(const ToyModelCorrections& other); // assignment
  virtual ~ToyModelCorrections() {;}

  enum DiJetVariable {
    kKt       = 0,
    kDiJetEta = 1
  };

  void     DoToy();

  void     WriteResponseToy();

  void     SetFileName(TString str)                {fFileName = str;}
  Bool_t   LoadFile();
  Bool_t   LoadList(TString strList);
  void     LoopLists();
  void     LoadDiJetResponse();

  Int_t    PickRandomSmearingDetector();

  void     SetFileNameDeltaPt(TString str)          {fFileNameDeltaPt = str;}
  void     LoadDeltaPtChris();
  Bool_t   LoadDeltaPt();
  Double_t PickRandomDeltaPt();
  void     PickRandomDeltaPtCorrelated(const Int_t nCentTot, Double_t &dptT, Double_t &dptA);

  void     PlotDiJetResponseKt();

  //Getters
  THnSparse *GetSmearingMatrixDetector();
  Double_t   GetRandomValueInBin(Int_t dim, Double_t *coord, THnSparse *hnsp, Double_t w);
  Double_t   GetRandomValueInBin(Int_t dim, Int_t bin, THnSparse *hnsp);

  //Setters
  void SetDiJetVariable(DiJetVariable v)   {fDiJetVariable = v; }
  void SetCentBin(Int_t i)                 {fCentBin    = i;}
  void SetPtMinTrack(Float_t p)            {fPtMinTrack = p;}
  void SetPtMinClus(Float_t p)             {fPtMinClus  = p;}
  void SetJetRadius(Float_t r)             {fJetRadius  = r;}
  void SetHadTrig(Int_t t)                 {fHadTrig    = t;}

  void SetNToyEvents(Int_t n)              {fNEvents    = n;}
  void SetDiJetType(Int_t t)               {fDiJetType  = t;}

  void SetToyResponseMatrix(THnSparse *h)  {fhnDiJetResponseToy = h;}

  void SetDeltaPtType(Int_t t)             {fDeltaPtType = t;}

 private:

 protected:
  DiJetVariable fDiJetVariable;               //variable to use
  Int_t         fCentBin;                     //centrality bin to use
  TString       fFileName;                    //name of root file with detector response
  TFile        *fFileIn;                      //root file
  TString       fListName;                    //name of list to use
  TList        *fListResponse;                //list with dijet response
  THnSparse    *fhnDiJetResponseCharged;      //6D dijet response matrix
  THnSparse    *fhnDiJetResponseFullCharged;  //6D dijet response matrix
  Float_t       fPtMinTrack;                  //min pT tracks
  Float_t       fPtMinClus;                   //min pT clusters
  Float_t       fJetRadius;                   //jet radius
  Int_t         fHadTrig;                     //leading track bias
  Int_t         fiPtTrigMC;                   //index pt trig MC
  Int_t         fiPtTrigDet;                  //index pt trig det
  Int_t         fiPtAssocMC;                  //index pt assoc MC
  Int_t         fiPtAssocDet;                 //index pt assoc det
  Int_t         fiDeltaPhiMC;                 //index delta phi dijet MC
  Int_t         fiDeltaPhiDet;                //index delta phi dijet Rec
  Int_t         fiVarMC;                      //index kt dijet MC
  Int_t         fiVarDet;                     //index kt dijet Rec
  TRandom3      fRandom;                      //random number generator
  TString       fFileNameDeltaPt;             //name of root file with delta pT
  TH1          *fh1DeltaPt;                   //DeltaPt Distribution
  TH1          *fh1DeltaPtCent[5];            //DeltaPt Distributions for 5 centrality bins (0-20,20-40,40-60,60-80,80-100)
  Int_t         fNCentBins;                   //#cent bins to consider for picking dpT
  TH1F         *fh1CentBinRndPick;            //bookkeep cent bins which are picked randomly for dpT
  Int_t         fNEvents;                     //number of toy event
  Int_t         fDiJetType;                   //type 0=charged-charged 1=full-charged
  THnSparse    *fhnDiJetResponseToy;          //6D dijet response matrix det+dpt smearing
  THnSparse    *fhnSparseReduced;             //!reduced matrix
  Int_t         fDeltaPtType;                 //deltapT type to use 0: RC  1: ExLJ 2: partial
};

#endif
