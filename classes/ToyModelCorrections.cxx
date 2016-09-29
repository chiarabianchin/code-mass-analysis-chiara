#include "ToyModelCorrections.h"
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
#include "TKey.h"

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

using std::cout;
using std::endl;

ClassImp(ToyModelCorrections)

//----------------------------------------------------------------------------
ToyModelCorrections::ToyModelCorrections():
  fDiJetVariable(kKt),
  fCentBin(0),
  fFileName(""),
  fFileIn(0x0),
  fListName(""),
  fListResponse(0x0),
  fhnDiJetResponseCharged(0x0),
  fhnDiJetResponseFullCharged(0x0),
  fPtMinTrack(0.15),
  fPtMinClus(0.3),
  fJetRadius(0.4),
  fHadTrig(0),
  fiPtTrigMC(0),
  fiPtTrigDet(1),
  fiPtAssocMC(2),
  fiPtAssocDet(3),
  fiDeltaPhiMC(4),
  fiDeltaPhiDet(5),
  fiVarMC(6),
  fiVarDet(7),
  fRandom(0),
  fFileNameDeltaPt(""),
  fh1DeltaPt(0x0),
  fNCentBins(1),
  fh1CentBinRndPick(0x0),
  fNEvents(1000),
  fDiJetType(0),
  fhnDiJetResponseToy(0x0),
  fhnSparseReduced(0x0),
  fDeltaPtType(2)
{
  for(Int_t i = 0; i<5; i++)
    fh1DeltaPtCent[i] = 0;

  fRandom.SetSeed(0);
}

//----------------------------------------------------------------------------
Bool_t ToyModelCorrections::LoadFile() {

  //Load input file

  if(!fFileIn)
    fFileIn = new TFile(fFileName.Data());

  fRandom.SetSeed(0);

  if(fFileIn) return kTRUE;
  else        return kFALSE;
}

//----------------------------------------------------------------------------
void ToyModelCorrections::LoopLists() {

  if(!LoadFile())
    return;

  fFileIn->ls();

  TIter next(fFileIn->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TString objname(key->GetClassName());
    TClass cls(objname);
    if (!cls.InheritsFrom("TList")) continue;
    TList *list = (TList*)key->ReadObj();
    Printf("list name: %s",list->GetName());
  }
}

//----------------------------------------------------------------------------
Bool_t ToyModelCorrections::LoadList(TString strList) {

  //Load requested list

  if(!LoadFile())
    return kFALSE;

  fListName = strList;

  fListResponse = (TList*)fFileIn->Get(fListName.Data());
  if(!fListResponse) {
    Printf("fListResponse %s not found",fListName.Data());
    return kFALSE;
  }
  else        
    return kTRUE;
}

//----------------------------------------------------------------------------
void ToyModelCorrections::LoadDiJetResponse() {
  //Load THnSparse with dijet response
  fListResponse->Print();
  fhnDiJetResponseCharged     = (THnSparse*)fListResponse->FindObject("fhnDiJetResponseCharged");
  fhnDiJetResponseFullCharged = (THnSparse*)fListResponse->FindObject("fhnDiJetResponseFullCharged");
}

//----------------------------------------------------------------------------
void ToyModelCorrections::LoadDeltaPtChris() {
  if(fh1DeltaPt) return;

  TFile *f = TFile::Open(fFileNameDeltaPt.Data());
  fh1DeltaPt = (TH1D*)f->Get("fh020DeltaPtSignal");
}

//----------------------------------------------------------------------------
Bool_t ToyModelCorrections::LoadDeltaPt() {
  /*
    fCentBin
    0:  0-20
    1:  20-40
    2:  40-60
    3:  60-80
    4:  80-100
    5:  0-100
    6:  40-100
    7:  0-2
    8:  2-5
    9:  5-10
    10: 10-20
    11: 0-40
   */

  if(fh1DeltaPt) return kTRUE;

  fh1CentBinRndPick = new TH1F("fh1CentBinRndPick","fh1CentBinRndPick",10,-0.5,9.5);

  TFile *f = TFile::Open(fFileNameDeltaPt.Data());

  TString strList = "";
  if(fDiJetType==0) 
    strList = Form("AliAnalysisTaskDeltaPtR%03dSparse_PicoTracks_RhoSparseR%03d_R%03d_TPC_histos",(Int_t)(100.*fJetRadius),(Int_t)(100.*fJetRadius),(Int_t)(100.*fJetRadius));
  else if(fDiJetType==1)
    strList = Form("AliAnalysisTaskDeltaPtR%03dSparse_PicoTracks_CaloClustersCorr_RhoSparseR%03d_Scaled_R%03d_EMCAL_histos",(Int_t)(100.*fJetRadius),(Int_t)(100.*fJetRadius),(Int_t)(100.*fJetRadius));
  TList *listDpt = (TList*)f->Get(strList.Data());
  if(!listDpt) {
    Printf("fDiJetType=%d, list for delta pT histos not found %s. Setting to 0x0",fDiJetType,strList.Data());
    f->ls();
    fh1DeltaPt = 0x0;
    return kFALSE;
  }

  TString strdpt = "fHistDeltaPtRC_";
  if(fDeltaPtType==0)
    strdpt = "fHistDeltaPtRC_";
  else if(fDeltaPtType==1)
    strdpt = "fHistDeltaPtRCExLJ_";
  else if(fDeltaPtType==2)
    strdpt = "fHistDeltaPtRCExPartialLJ_";

  fNCentBins = 1;
  fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),0));

  for(Int_t i = 0; i<5; i++)
    fh1DeltaPtCent[i] = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),i));

  if(fCentBin<=4)
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),fCentBin));
  else if(fCentBin==5) { //0-100%
    fh1DeltaPt = (TH1D*)(fh1DeltaPtCent[0]->Clone("fh1DeltaPt"));
    for(Int_t i = 1; i<5; i++)
      fh1DeltaPt->Add(fh1DeltaPtCent[i]);
    fNCentBins = 5;
  }
  else if(fCentBin==6) { //40-100%
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),2));
    TH1 *hAdd1 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),3));
    TH1 *hAdd2 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),4));
    fh1DeltaPt->Add(hAdd1);
    fh1DeltaPt->Add(hAdd2);
    fNCentBins = 3;
  }
  else if(fCentBin==7) {// 0-2%
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),0));
    TH1 *hAdd1 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),1));
    fh1DeltaPt->Add(hAdd1);
    fNCentBins = 1;
  }
  else if(fCentBin==8) {// 2-5%
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),2));
    TH1 *hAdd1 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),3));
    TH1 *hAdd2 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),4));
    fh1DeltaPt->Add(hAdd1);
    fh1DeltaPt->Add(hAdd2);
    fNCentBins = 1;
  }
  else if(fCentBin==9) {// 5-10%
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),0));
    fNCentBins = 1;
  }
  else if(fCentBin==10) {// 10-20%
    fh1DeltaPt = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),5));
    TH1 *hAdd1 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),6));
    TH1 *hAdd2 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),7));
    TH1 *hAdd3 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),8));
    TH1 *hAdd4 = (TH1*)listDpt->FindObject(Form("%s%d",strdpt.Data(),9));
    fh1DeltaPt->Add(hAdd1);
    fh1DeltaPt->Add(hAdd2);
    fh1DeltaPt->Add(hAdd3);
    fh1DeltaPt->Add(hAdd4);
    fNCentBins = 1;
  }
  else if(fCentBin==11)
    fNCentBins = 2;

  if(!fh1DeltaPt) {
    Printf("fh1DeltaPt could not be loaded %s... please check",Form("%s_%d",strdpt.Data(),fCentBin));
    listDpt->Print();
    return kFALSE;
  }
  else
    return kTRUE;
}


//----------------------------------------------------------------------------
THnSparse* ToyModelCorrections::GetSmearingMatrixDetector() {
  //
  if(!fhnSparseReduced) {
    THnSparse *hnsp = NULL;
    if(fDiJetType==0)
      hnsp = fhnDiJetResponseCharged;
    else if(fDiJetType==1)
      hnsp = fhnDiJetResponseFullCharged;

    if(!hnsp) {
      Printf("Cannot find detector response");
      return 0x0;
    }
    
    Int_t ndimVar[2] = {6,8}; //also store eta,dijet in case of dijet eta analysis
    const Int_t ndim = ndimVar[fDiJetVariable];
    Int_t dim[ndim];// = {0,1,2,3,4,5};
    for(Int_t i = 0; i<ndim; i++)
      dim[i] = i;
    
    fhnSparseReduced = static_cast<THnSparse*>(hnsp->ProjectionAny(ndim,dim,kTRUE,"E"));
  }

  return fhnSparseReduced;
}

//----------------------------------------------------------------------------
Int_t ToyModelCorrections::PickRandomSmearingDetector() {
  //
  // type 0=charged-charged 1=full-charged
  //
  THnSparse *hnsp = GetSmearingMatrixDetector();
  if(!hnsp) return -1;
  Int_t nBins = hnsp->GetNbins();
  Int_t bin = fRandom.Uniform(0,nBins);

  /*
  Int_t nDim = hnsp->GetNdimensions();
  Double_t* coord = new Double_t[nDim];

  hnsp->GetRandom(coord,kTRUE);
  Int_t bin = hnsp->GetBin(coord);
  delete [] coord;
  */

  /*
  hnsp->PrintBin(bin,"0");

  Int_t nDim = hnsp->GetNdimensions();
  Int_t* coord = new Int_t[nDim];
  Double_t v = hnsp->GetBinContent(bin,coord);
  */

  // for(Int_t i=0; i<nDim; i++)
  //   printf("coord[%d] = %d %f ",i,coord[i],hnsp->GetAxis(i)->GetBinCenter(coord[i]));
  // Printf("bin content %e",v);

  return bin;
}

//----------------------------------------------------------------------------
Double_t ToyModelCorrections::PickRandomDeltaPt() {
  //pick random delta pt value

  Bool_t b = LoadDeltaPt();
  if(b)
    return fh1DeltaPt->GetRandom();
  else
    return -999;
}

//----------------------------------------------------------------------------
void ToyModelCorrections::PickRandomDeltaPtCorrelated(const Int_t nCentTot, Double_t &dptT, Double_t &dptA) {

  //pick random delta pt value
  //delta pt of trigger and assoc drawn from same centrality bin

  Bool_t b = LoadDeltaPt();
  if(!b) {
    dptT = -999;
    dptA = -999;
  }

  if(nCentTot==1) {
    dptT = PickRandomDeltaPt();
    dptA = PickRandomDeltaPt();
  }
  else {
    Double_t rnd = fRandom.Rndm();
    Int_t centBin = TMath::FloorNint(nCentTot*rnd);
    if(fCentBin==6) centBin+=2;
    dptT = fh1DeltaPtCent[centBin]->GetRandom();
    dptA = fh1DeltaPtCent[centBin]->GetRandom();

    fh1CentBinRndPick->Fill(centBin);
  }

}

//----------------------------------------------------------------------------
Double_t ToyModelCorrections::GetRandomValueInBin(Int_t dim, Double_t *coord, THnSparse *hnsp, Double_t w) {
  // Get random interpolated value in bin

  Int_t nDim = hnsp->GetNdimensions();
  Int_t* coordLow = new Int_t[nDim];
  Int_t* coordUp = new Int_t[nDim];
  for(Int_t j=0; j<nDim; j++) {
    if(j==dim) {
      coordLow[j] = coord[j]-1;
      coordUp[j]  = coord[j]+1;
    }
    else {
      coordLow[j] = coord[j];
      coordUp[j]  = coord[j];
    }
  }
  
  Double_t dx = (hnsp->GetAxis(dim)->GetBinCenter(coordUp[dim])-hnsp->GetAxis(dim)->GetBinCenter(coordLow[dim]));
  Int_t binLow = hnsp->GetBin(coordLow);
  Int_t binUp = hnsp->GetBin(coordUp);
  Double_t dy = hnsp->GetBinContent(binUp,coordUp)-hnsp->GetBinContent(binLow,coordLow);

  Double_t b = 0.;
  if(dx<0. || dx>0.)
    b = dy/dx;
  Double_t a = w - b*hnsp->GetAxis(dim)->GetBinCenter(coord[dim]);
  TF1 func("func","[0]+[1]*x",0.,500.);
  func.SetParameters(a,b);
  Double_t val = func.GetRandom(hnsp->GetAxis(dim)->GetBinLowEdge(coord[dim]),hnsp->GetAxis(dim)->GetBinUpEdge(coord[dim]));

  delete [] coordLow;
  delete [] coordUp;

  return val;
}

//----------------------------------------------------------------------------
Double_t ToyModelCorrections::GetRandomValueInBin(Int_t dim, Int_t bin, THnSparse *hnsp) {
  // Get random interpolated value in bin

  Int_t nDim = hnsp->GetNdimensions();
  Int_t* coord = new Int_t[nDim];
  Int_t* coordLow = new Int_t[nDim];
  Int_t* coordUp = new Int_t[nDim];

  Double_t w = hnsp->GetBinContent(bin,coord);
 
  for(Int_t j=0; j<nDim; j++) {
    if(j==dim) {
      coordLow[j] = coord[j]-1;
      coordUp[j]  = coord[j]+1;
    }
    else {
      coordLow[j] = coord[j];
      coordUp[j]  = coord[j];
    }
  }
  
  /*
  Double_t dx = (hnsp->GetAxis(dim)->GetBinCenter(coordLow[dim]) - hnsp->GetAxis(dim)->GetBinCenter(coordUp[dim]));
  Int_t binLow = hnsp->GetBin(coordLow);
  Int_t binUp = hnsp->GetBin(coordUp);
  Double_t dy = hnsp->GetBinContent(binLow,coordLow) - hnsp->GetBinContent(binUp,coordUp);
  */

  TH1D *hproj = (TH1D*)hnsp->Projection(dim,"E");
  Double_t dx = (hproj->GetXaxis()->GetBinCenter(coordLow[dim]) - hproj->GetXaxis()->GetBinCenter(coordUp[dim]));
  Double_t dy = hproj->GetBinContent(coordLow[dim]) - hproj->GetBinContent(coordUp[dim]);

  Double_t b = 0.;
  if(dx<0. || dx>0.)
    b = dy/dx;
  Double_t a = w - b*hnsp->GetAxis(dim)->GetBinCenter(coord[dim]);
  //Printf("a=%f b=%f dx=%f",a,b,dx);
  TF1 func("func","[0]+[1]*x",hnsp->GetAxis(dim)->GetBinLowEdge(coord[dim]),hnsp->GetAxis(dim)->GetBinUpEdge(coord[dim]));
  func.SetParameters(a,b);
  Double_t val = func.GetRandom(hnsp->GetAxis(dim)->GetBinLowEdge(coord[dim]),hnsp->GetAxis(dim)->GetBinUpEdge(coord[dim]));

  delete    hproj;
  delete [] coord;
  delete [] coordLow;
  delete [] coordUp;

  return val;
}

//----------------------------------------------------------------------------
void ToyModelCorrections::DoToy() {
  //run toy

  THnSparse *hnsp;
  hnsp = GetSmearingMatrixDetector();
  if(!hnsp) return;
  if(fhnDiJetResponseToy) delete fhnDiJetResponseToy;

  if(fDiJetType==0)
    fhnDiJetResponseToy = static_cast<THnSparse*>(fhnDiJetResponseCharged->Clone("fhnDiJetResponseToy"));
  else if(fDiJetType==1)
    fhnDiJetResponseToy = static_cast<THnSparse*>(fhnDiJetResponseFullCharged->Clone("fhnDiJetResponseToy"));
  fhnDiJetResponseToy->Reset();
  fhnDiJetResponseToy->Sumw2();

  Int_t nDim = hnsp->GetNdimensions();
  Int_t* coord = new Int_t[nDim];
  for(Int_t i=0; i<fNEvents; i++) {
    if(i % 100000 == 0) cout << " " << i << flush;
    Int_t bin = PickRandomSmearingDetector();
    Double_t dptTrig  = 0.;//PickRandomDeltaPt();
    Double_t dptAssoc = 0.;//PickRandomDeltaPt();//should depend on cent bin of first
    PickRandomDeltaPtCorrelated(fNCentBins,dptTrig,dptAssoc);

    Double_t w = hnsp->GetBinContent(bin,coord);

    Double_t ptTrigGen = hnsp->GetAxis(0)->GetBinCenter(coord[0]) + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(0)->GetBinWidth(coord[0]) ;
    Double_t ptAssoGen = hnsp->GetAxis(2)->GetBinCenter(coord[2])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(2)->GetBinWidth(coord[2]);
    Double_t dphiGen = hnsp->GetAxis(4)->GetBinCenter(coord[4])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(4)->GetBinWidth(coord[4]);

    // Double_t ptTrigGen  = GetRandomValueInBin(fiPtTrigMC, bin, hnsp);//coord, hnsp, w);
    // Double_t ptAssoGen  = GetRandomValueInBin(fiPtAssocMC, bin, hnsp);//coord, hnsp, w);
    // Double_t dphiGen    = GetRandomValueInBin(fiDeltaPhiMC, bin, hnsp);//coord, hnsp, w);

    Double_t varGen = 0.;
    if(fDiJetVariable==kKt)
      varGen = ptTrigGen*TMath::Sin(dphiGen);
    else if(fDiJetVariable==kDiJetEta)
      varGen = hnsp->GetAxis(6)->GetBinCenter(coord[6])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(6)->GetBinWidth(coord[6]);      
    // varGen  = GetRandomValueInBin(fiVarMC, bin, hnsp);//coord, hnsp, w); 
    
    Double_t trueVar[4] = {
      ptTrigGen,
      ptAssoGen,
      dphiGen,
      varGen
    };

    Double_t ptTrigRec = hnsp->GetAxis(1)->GetBinCenter(coord[1]) + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(1)->GetBinWidth(coord[1]) + dptTrig;
    Double_t ptAssoRec = hnsp->GetAxis(3)->GetBinCenter(coord[3])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(3)->GetBinWidth(coord[3]) + dptAssoc;
    Double_t dphiRec = hnsp->GetAxis(5)->GetBinCenter(coord[5])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(5)->GetBinWidth(coord[5]);

    // Double_t ptTrigRec  = GetRandomValueInBin(fiPtTrigDet, bin, hnsp) + dptTrig;//coord, hnsp, w);
    // Double_t ptAssoRec  = GetRandomValueInBin(fiPtAssocDet, bin, hnsp) + dptAssoc;//coord, hnsp, w);
    // Double_t dphiRec    = GetRandomValueInBin(fiDeltaPhiDet, bin, hnsp);//coord, hnsp, w);

    Double_t varRec = 0.;
    if(fDiJetVariable==kKt)
      varRec = ptTrigRec*TMath::Sin(dphiRec);
    else if(fDiJetVariable==kDiJetEta)
      varRec = hnsp->GetAxis(7)->GetBinCenter(coord[7])  + (fRandom.Rndm() - 0.5) * hnsp->GetAxis(7)->GetBinWidth(coord[7]);
    //      varRec  = GetRandomValueInBin(fiVarDet, bin, hnsp);

    Double_t smearVar[4] = {
      ptTrigRec,
      ptAssoRec,
      dphiRec,
      varRec
    };
    
    Double_t diJetResponse[8] = {trueVar[0],smearVar[0],trueVar[1],smearVar[1],trueVar[2],smearVar[2],trueVar[3],smearVar[3]};
    fhnDiJetResponseToy->Fill(diJetResponse,w);
  }
  delete [] coord;
}

//----------------------------------------------------------------------------
void ToyModelCorrections::WriteResponseToy() {
  //Write response matrix from toy
  TString type = "";
  if(fDiJetType==0)
    type = "ChCh";
  else if(fDiJetType==1)
    type = "FuCh";

  TFile *respOut = new TFile(Form("DiJetResponseToyR%03dHadTrig%dRV%dCent%d%s.root",(Int_t)(100.*fJetRadius),fHadTrig,fDiJetVariable,fCentBin,type.Data()),"RECREATE");

  fhnDiJetResponseToy->Write(Form("fhnDiJetResponseToyR%03dHadTrig%dRV%dCent%d",(Int_t)(100.*fJetRadius),fHadTrig,fDiJetVariable,fCentBin));

  THnSparse *hnsp;
  hnsp = GetSmearingMatrixDetector();
  if(hnsp)
    hnsp->Write();

  fh1CentBinRndPick->Write();
  
  respOut->Write();
  respOut->Close();
}

//----------------------------------------------------------------------------
void ToyModelCorrections::PlotDiJetResponseKt() {

  //Some axis range definitions
  const Int_t nCorr = 4;
  Double_t ptTriggerMin[nCorr] = {20.,40.,60.,80.};
  Double_t ptTriggerMax[nCorr] = {40.,60.,80.,120};

  Double_t ptAssocChMin[nCorr] = {20.,20.,20.,20.};
  Double_t ptAssocChMax[nCorr] = {110.,110.,110.,110.};

  Double_t ptAssocFullChMin[nCorr] = {20.,20.,20.,20.};
  Double_t ptAssocFullChMax[nCorr] = {110.,110.,110.,110.};

  //Store dijet vars: pt,trig MC, pt,trig DET, pt,ass MC, pt,ass DET, dPhi MC, dPhi Det, kT MC, kT Det 
  const Int_t ptTrigMC       = 0;
  const Int_t ptTrigDet      = 1;
  const Int_t ptAssocMC      = 2;
  const Int_t ptAssocDet     = 3;
  const Int_t deltaPhiMC     = 4;
  const Int_t deltaPhiDet    = 5;
  const Int_t ktMC           = 6;
  const Int_t ktDet          = 7;

  TH1D *hPtTrigPart = static_cast<TH1D*>(fhnDiJetResponseToy->Projection(ptTrigMC,"E"));
  TH1D *hPtTrigDet = static_cast<TH1D*>(fhnDiJetResponseToy->Projection(ptTrigDet,"E"));

  Double_t nTriggersPart[nCorr] = {1};
  Double_t nTriggersDet[nCorr] = {1};

  Int_t binMin, binMax;

  for(Int_t i=0; i<nCorr; i++) {
    binMin = hPtTrigPart->FindBin(ptTriggerMin[i]);
    binMax = hPtTrigPart->FindBin(ptTriggerMax[i])-1;
    nTriggersPart[i] = hPtTrigPart->Integral(binMin,binMax);

    binMin = hPtTrigDet->FindBin(ptTriggerMin[i]);
    binMax = hPtTrigDet->FindBin(ptTriggerMax[i])-1;
    nTriggersDet[i] = hPtTrigDet->Integral(binMin,binMax);
  }

  TH1D *hKtPart[nCorr];
  TH1D *hKtDet[nCorr];
  TH1D *hKtDetFeedIn[nCorr];
  TH1D *hKtDetFeedOut[nCorr];
  TH2D *hKtDetPtTrigPart[nCorr];
  TLegend *legkt[nCorr];

  //Set DeltaPhi window for kT measurement
  Double_t dPhiWindow = 12.*fhnDiJetResponseToy->GetAxis(deltaPhiMC)->GetBinWidth(1);
  cout << "dPhiWindow: " << dPhiWindow << endl;
  fhnDiJetResponseToy->GetAxis(deltaPhiMC)->SetRangeUser(TMath::Pi()-dPhiWindow,TMath::Pi()+dPhiWindow);
 
  TCanvas *c2 = new TCanvas("c2","c2: Kt ",1000,800);
  c2->Divide(2,2);

  for(Int_t i=0; i<nCorr; i++) {
    fhnDiJetResponseToy->GetAxis(ptTrigDet)->SetRangeUser(0.,1e6);
    fhnDiJetResponseToy->GetAxis(ptAssocDet)->SetRangeUser(0.,1e6);
    fhnDiJetResponseToy->GetAxis(deltaPhiDet)->SetRangeUser(-0.5*TMath::Pi(),1.5*TMath::Pi());
    fhnDiJetResponseToy->GetAxis(ptTrigMC)->SetRangeUser(ptTriggerMin[i],ptTriggerMax[i]);
    fhnDiJetResponseToy->GetAxis(ptAssocMC)->SetRangeUser(ptAssocChMin[i],ptAssocChMax[i]);
    fhnDiJetResponseToy->GetAxis(deltaPhiMC)->SetRangeUser(TMath::Pi()-dPhiWindow,TMath::Pi()+dPhiWindow);
    
    hKtPart[i] = (TH1D*)fhnDiJetResponseToy->Projection(ktMC,"E");

    fhnDiJetResponseToy->GetAxis(ptTrigMC)->SetRangeUser(0.,1e6);
    fhnDiJetResponseToy->GetAxis(ptAssocMC)->SetRangeUser(0.,1e6);
    fhnDiJetResponseToy->GetAxis(ptTrigDet)->SetRangeUser(ptTriggerMin[i],ptTriggerMax[i]);
    fhnDiJetResponseToy->GetAxis(ptAssocDet)->SetRangeUser(ptAssocChMin[i],ptAssocChMax[i]);
    fhnDiJetResponseToy->GetAxis(deltaPhiMC)->SetRangeUser(-0.5*TMath::Pi(),1.5*TMath::Pi());
    fhnDiJetResponseToy->GetAxis(deltaPhiDet)->SetRangeUser(TMath::Pi()-dPhiWindow,TMath::Pi()+dPhiWindow);
    
    hKtDetPtTrigPart[i] = (TH2D*)fhnDiJetResponseToy->Projection(ktDet,ptTrigMC);
    
    binMin = hKtDetPtTrigPart[i]->GetXaxis()->FindBin(ptTriggerMin[i]);
    binMax = hKtDetPtTrigPart[i]->GetXaxis()->FindBin(ptTriggerMax[i])-1;
    hKtDet[i] = (TH1D*)hKtDetPtTrigPart[i]->ProjectionY(Form("hKtDet%d",i),binMin,binMax);
    hKtDetFeedIn[i] = (TH1D*)hKtDetPtTrigPart[i]->ProjectionY(Form("hKtDetFeedIn%d",i),1,binMin);
    hKtDetFeedOut[i] = (TH1D*)hKtDetPtTrigPart[i]->ProjectionY(Form("hKtDetFeedOut%d",i),binMax,hKtDetPtTrigPart[i]->GetXaxis()->GetNbins());


    if(nTriggersPart[i]>0) 
      hKtPart[i]->Scale(1./nTriggersPart[i]);
    if(nTriggersDet[i]>0) {
      hKtDet[i]->Scale(1./nTriggersDet[i]);
      hKtDetFeedIn[i]->Scale(1./nTriggersDet[i]);
      hKtDetFeedOut[i]->Scale(1./nTriggersDet[i]);
    }

    c2->cd(i+1);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    TH1F *hFrame = gPad->DrawFrame(-100,0.,100.,0.05);
    hFrame->SetXTitle("#it{k}_{T} = #it{p}_{T,jet} sin(#Delta#varphi)");
    hFrame->SetYTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}}{d#it{k}_{T}}");
    hFrame->GetXaxis()->SetLabelSize(0.05);
    hFrame->GetYaxis()->SetLabelSize(0.05);
    hFrame->GetXaxis()->SetTitleSize(0.06);
    hFrame->GetYaxis()->SetTitleSize(0.06);
    hFrame->GetXaxis()->SetTitleOffset(0.9);
    hFrame->GetYaxis()->SetTitleOffset(1.5);

    hKtDet[i]->SetLineColor(2);
    hKtDet[i]->SetMarkerColor(2);
    hKtDet[i]->SetLineWidth(2);
    hKtDet[i]->SetMarkerStyle(20);

    hKtDetFeedIn[i]->SetLineColor(4);
    hKtDetFeedIn[i]->SetMarkerColor(4);
    hKtDetFeedIn[i]->SetLineWidth(2);
    hKtDetFeedIn[i]->SetMarkerStyle(21);

    hKtDetFeedOut[i]->SetLineColor(8);
    hKtDetFeedOut[i]->SetMarkerColor(8);
    hKtDetFeedOut[i]->SetLineWidth(2);
    hKtDetFeedOut[i]->SetMarkerStyle(25);

    hKtPart[i]->SetLineColor(1);
    hKtPart[i]->SetMarkerColor(1);
    hKtPart[i]->SetLineWidth(2);
    hKtPart[i]->SetMarkerStyle(24);

    hKtPart[i]->DrawCopy("same histo");
    hKtDet[i]->DrawCopy("same histo");
    hKtDetFeedIn[i]->DrawCopy("same histo");
    hKtDetFeedOut[i]->DrawCopy("same histo");
       
    legkt[i] = new TLegend(0.22,0.65,0.6,0.95);
    legkt[i]->SetFillColor(10);
    legkt[i]->SetFillStyle(0);
    legkt[i]->SetBorderSize(0);
    legkt[i]->AddEntry(hKtPart[i],"Particle level","p");
    legkt[i]->AddEntry(hKtDet[i],"Detector level","p");
    legkt[i]->AddEntry(hKtDetFeedIn[i],"Detector level Feed-In low p_{T}","p");
    legkt[i]->AddEntry(hKtDetFeedOut[i],"Detector level Feed-In high p_{T}","p");
    legkt[i]->Draw();

    TLatex textTrig;
    textTrig.SetNDC();
    textTrig.DrawLatex(0.6,0.88,Form("%.0f<#it{p}_{T,jet,part}^{trigger}<%.0f GeV/#it{c}",ptTriggerMin[i],ptTriggerMax[i]));

    TLatex textAssoc;
    textAssoc.SetNDC();
    textAssoc.DrawLatex(0.6,0.77,Form("%.0f<#it{p}_{T,jet,part}^{assoc}<%.0f GeV/#it{c}",ptAssocFullChMin[i],ptAssocFullChMax[i]));
  }
}
